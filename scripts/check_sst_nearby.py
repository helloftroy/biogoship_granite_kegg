#!/usr/bin/env python3
"""Would moving the tile a little recover SST? A check; no tile is changed.

For every sample whose best tile (fixed or advected) has no SST, every scene of its
site is re-examined around the target position. A window qualifies if it is clear
by OLCI (at least --min-valid, as when tiles are cut) and has SST, which under the
nearest-neighbour resampling of tiles.SstSwath means at least one acceptable SLSTR
pixel inside the window's box widened by tiles.SST_SEARCH_PAD_DEG. For each scene
and SST source the script records how far the nearest qualifying window is from the
target, so the gain of allowing a move of 20 or 40 km can be read off.

SST sources tried per scene: the SST orbit paired with the scene, then the first
--fill-candidates daytime overpasses ranked by fill_missing_sst.py for that scene's
tile (from data/sst_fill_candidates.csv; scenes without a tile have none). Quality
level 4-5 is what tiles use; level 3 is counted separately to show what relaxing it
would add.

Outputs: data/sst_nearby_check.csv (one row per sample position, scene and SST
source) and a summary per cruise in the log.

Usage:
    python scripts/check_sst_nearby.py
    python scripts/check_sst_nearby.py --fill-candidates 0     # paired SST only, nothing downloaded
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship import tiles  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, latest_tile_log_rows, read_tile_log, wrap_lon_delta  # noqa: E402
from fill_missing_sst import Filler, extract, load_candidate_cache, local_archive  # noqa: E402

RADII_KM = (0.0, 10.0, 20.0, 40.0)
STRIDE = 3  # window centres tested every 3 pixels (0.9 km)


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest_dir / "biogoship_tile_manifest.csv")
    best = manifest[manifest["is_best"] == 1]
    lacking = best[pd.to_numeric(best["sst_valid_fraction"], errors="coerce").fillna(0) <= 0]
    # Samples with no clear tile at all are out of scope: there is no tile to move.
    wanted = lacking[["sample_id", "site_id", "position", "section"]].drop_duplicates()
    print(f"Sample positions whose best tile has no SST: {len(wanted):,} "
          f"({wanted['sample_id'].nunique():,} samples, {wanted[['site_id', 'position']].drop_duplicates().shape[0]:,} site positions)")  # fmt: skip

    log = pd.DataFrame(latest_tile_log_rows(read_tile_log(args.data_dir / "tile_log.csv")))
    log = log[log["status"].isin(["valid", "cloudy", "edge"]) & (log["row"] != "")]
    keys = set(map(tuple, wanted[["site_id", "position"]].to_numpy()))
    scenes = log[[(s, p) in keys for s, p in zip(log["site_id"], log["position"])]].copy()
    products = pd.read_csv(args.manifest_dir / "products.csv").set_index("product_id")
    candidates = load_candidate_cache(args.data_dir / "sst_fill_candidates.csv", args.candidate_hours)
    print(f"Scenes to re-examine: {len(scenes):,} over {scenes['product_id'].nunique():,} OLCI products")

    filler = Filler(argparse.Namespace(**vars(args), keep_work=False))
    filler.work = args.data_dir / "work" / "sst_nearby_check"
    filler.work.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, (product_id, group) in enumerate(scenes.groupby("product_id"), start=1):
        product = products.loc[product_id]
        try:
            lat, lon, clear = load_scene(args, filler.work, product["product_name"])
        except Exception as exc:  # noqa: BLE001 - one bad archive must not stop the check
            print(f"  {product['product_name'][:31]}: {exc}", flush=True)
            continue
        for _, scene in group.iterrows():
            sources = []
            if isinstance(product.get("sst_product_name"), str) and product["sst_product_name"]:
                sources.append(("paired", {"product_id": product["sst_product_id"], "product_name": product["sst_product_name"]}))
            daytime = [c for c in candidates.get(scene["tile_path"], []) if not c["night"]] if scene["tile_path"] else []
            sources += [(f"candidate{k + 1}", c) for k, c in enumerate(daytime[: args.fill_candidates])]
            for source, sst_product in sources:
                try:
                    swath = filler.swath(sst_product)
                except Exception as exc:  # noqa: BLE001
                    print(f"  SST {sst_product['product_name'][:31]}: {exc}", flush=True)
                    continue
                found = nearest_windows(lat, lon, clear, int(scene["row"]), int(scene["col"]), swath, args)
                rows.append({"site_id": scene["site_id"], "position": scene["position"], "product_id": product_id,
                             "time_gap_hours": scene["time_gap_hours"], "scene_status": scene["status"],
                             "sst_source": source, "sst_product": sst_product["product_name"],
                             "sst_gap_hours": sst_product.get("gap_hours", ""), **found})  # fmt: skip
        if i % 10 == 0 or i == scenes["product_id"].nunique():
            print(f"  {i:,}/{scenes['product_id'].nunique():,} products", flush=True)
    filler.cleanup()

    out = pd.DataFrame(rows)
    path = args.data_dir / "sst_nearby_check.csv"
    out.to_csv(path, index=False)
    summarise(wanted, out)
    print(f"\nWrote {path}")


def load_scene(args: argparse.Namespace, work: Path, product_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    archive = local_archive(args, product_name)
    if archive is None:
        raise RuntimeError("OLCI archive not on disk")
    folder = work / product_name
    folder.mkdir(parents=True, exist_ok=True)
    for name in tiles.OLCI_CHECK_FILES:
        if not (folder / name).exists():
            extract(archive, name, folder / name)
    lat, lon = tiles.load_geo(folder / "geo_coordinates.nc")
    clear = tiles.ocean_mask(tiles.load_wqsf(folder / "wqsf.nc"))
    for name in tiles.OLCI_CHECK_FILES:
        (folder / name).unlink(missing_ok=True)
    folder.rmdir()
    return lat, lon, clear


def good_sst_points(swath: tiles.SstSwath, lat_box: tuple[float, float], lon_centre: float, lon_half: float, min_quality: int):
    """Lat/lon of SLSTR pixels at or above min_quality inside a lat/lon box."""
    near = (swath.lat >= lat_box[0]) & (swath.lat <= lat_box[1]) & (np.abs(wrap_lon_delta(swath.lon - lon_centre)) <= lon_half)
    rows = np.flatnonzero(near.any(axis=1))
    if rows.size == 0:
        return np.empty(0, np.float32), np.empty(0, np.float32)
    r0, r1 = int(rows[0]), int(rows[-1]) + 1
    with h5py.File(swath.path, "r") as f:
        raw = f["sea_surface_temperature"][0, r0:r1, :]
        fill = int(np.ravel(f["sea_surface_temperature"].attrs.get("_FillValue", -32768))[0])
        quality = f["quality_level"][0, r0:r1, :]
    good = near[r0:r1] & (quality >= min_quality) & (raw != fill)
    return swath.lat[r0:r1][good], swath.lon[r0:r1][good]


def nearest_windows(lat, lon, clear, row: int, col: int, swath: tiles.SstSwath, args) -> dict:
    """Distance (km) from the target to the nearest OLCI-clear window with SST, for quality 4+ and 3+."""
    max_px = int(round(max(RADII_KM) / tiles.PIXEL_KM))
    r_lo, r_hi = max(row - tiles.HALF - max_px, 0), min(row - tiles.HALF + max_px, lat.shape[0] - tiles.TILE_SIZE)
    c_lo, c_hi = max(col - tiles.HALF - max_px, 0), min(col - tiles.HALF + max_px, lat.shape[1] - tiles.TILE_SIZE)
    out = {"olci_clear_km": np.nan, "sst_q4_km": np.nan, "sst_q3_km": np.nan}
    if r_lo > r_hi or c_lo > c_hi:
        return out
    region = clear[r_lo : r_hi + tiles.TILE_SIZE, c_lo : c_hi + tiles.TILE_SIZE]
    total = np.pad(region.astype(np.int32).cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    T = tiles.TILE_SIZE
    fraction = (total[T:, T:] - total[:-T, T:] - total[T:, :-T] + total[:-T, :-T]) / float(T * T)
    # Grid of window origins aligned so the window centred on the target is one of them.
    r0, c0 = row - tiles.HALF, col - tiles.HALF
    origins_r, origins_c = np.meshgrid(np.arange(r_lo + (r0 - r_lo) % STRIDE, r_hi + 1, STRIDE),
                                       np.arange(c_lo + (c0 - c_lo) % STRIDE, c_hi + 1, STRIDE), indexing="ij")  # fmt: skip
    frac = fraction[origins_r - r_lo, origins_c - c_lo]
    dist_km = np.hypot(origins_r - r0, origins_c - c0) * tiles.PIXEL_KM
    ok_olci = (frac >= args.min_valid) & (dist_km <= max(RADII_KM))
    if not ok_olci.any():
        return out
    out["olci_clear_km"] = round(float(dist_km[ok_olci].min()), 2)

    centre_lat = lat[origins_r + tiles.HALF, origins_c + tiles.HALF][ok_olci].astype(np.float64)
    centre_lon = lon[origins_r + tiles.HALF, origins_c + tiles.HALF][ok_olci].astype(np.float64)
    keep = np.isfinite(centre_lat) & np.isfinite(centre_lon)
    centre_lat, centre_lon, dists = centre_lat[keep], centre_lon[keep], dist_km[ok_olci][keep]
    if not len(dists):
        return out
    # Window box as in SstSwath.resample: half the tile (~6.3 km) plus the search pad.
    half_deg = tiles.HALF * tiles.PIXEL_KM / 111.195 + tiles.SST_SEARCH_PAD_DEG
    lon0 = float(lon[row, col])
    cos_lat = max(np.cos(np.radians(float(lat[row, col]))), 0.05)
    lat_box = (float(centre_lat.min()) - half_deg, float(centre_lat.max()) + half_deg)
    lon_half = float(np.abs(wrap_lon_delta(centre_lon - lon0)).max()) + half_deg / cos_lat
    for quality, key in ((tiles.SLSTR_MIN_QUALITY, "sst_q4_km"), (3, "sst_q3_km")):
        plat, plon = good_sst_points(swath, lat_box, lon0, lon_half, quality)
        if not len(plat):
            continue
        has = np.zeros(len(dists), bool)
        for start in range(0, len(dists), 512):
            dlat = np.abs(centre_lat[start : start + 512, None] - plat[None, :])
            dlon = np.abs(wrap_lon_delta(centre_lon[start : start + 512, None] - plon[None, :])) * cos_lat
            has[start : start + 512] = ((dlat <= half_deg) & (dlon <= half_deg)).any(axis=1)
        if has.any():
            out[key] = round(float(dists[has].min()), 2)
    return out


def summarise(wanted: pd.DataFrame, out: pd.DataFrame) -> None:
    if out.empty:
        print("Nothing examined.")
        return
    print("\nCheck: the current window (0 km) with paired SST should almost never have SST:")
    paired = out[out["sst_source"] == "paired"]
    print(f"  paired rows with SST at <= 1 km: {(paired['sst_q4_km'] <= 1).sum():,} of {len(paired):,}")
    for label, frame in (("paired SST only", paired), ("paired + fill candidates", out)):
        per = frame.groupby(["site_id", "position"])[["olci_clear_km", "sst_q4_km", "sst_q3_km"]].min()
        joined = wanted.merge(per.reset_index(), on=["site_id", "position"], how="left")
        print(f"\nSample positions recovered within a move of ... ({label}); n = {len(joined):,}")
        table = {}
        for radius in RADII_KM[1:]:
            table[f"q4 <= {radius:g} km"] = joined.groupby("section")["sst_q4_km"].apply(lambda d, r=radius: int((d <= r).sum()))
            table[f"q3 <= {radius:g} km"] = joined.groupby("section")["sst_q3_km"].apply(lambda d, r=radius: int((d <= r).sum()))
        table["total"] = joined.groupby("section").size()
        frame_out = pd.DataFrame(table)
        frame_out.loc["all"] = frame_out.sum()
        print(frame_out.to_string())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data")
    parser.add_argument("--products-dir", type=Path, default=REPO_ROOT / "data" / "products")
    parser.add_argument("--min-valid", type=float, default=0.8, help="OLCI clear fraction a window needs, as for tiles.")
    parser.add_argument("--fill-candidates", type=int, default=2, help="Daytime fill overpasses tried per scene besides the paired one.")
    parser.add_argument("--candidate-hours", type=float, default=120.0, help="Window of the fill candidate search to read.")
    parser.add_argument("--sst-cache", type=int, default=6, help="Downloaded SST files kept for reuse.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
