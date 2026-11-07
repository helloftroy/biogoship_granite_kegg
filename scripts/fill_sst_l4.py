#!/usr/bin/env python3
"""Copies of the tiles with every SST gap filled from OSTIA L4 SST, shifted toward daytime skin SST.

Why: SLSTR SST is missing on many tiles (cloud screening stricter than OLCI's;
I07N worst), so either those samples drop out or Granite sees a flat mean-filled
channel. OSTIA (Met Office, 0.05 deg daily, gap-free) is a foundation-SST analysis:
no daytime warming and no cool skin, so it runs a few tenths of a degree off the
daytime SLSTR skin SST that the rest of the tiles carry. It is shifted toward that:

  offset = median(SLSTR - OSTIA) over the tile's own SLSTR pixels, when the tile
           has at least --min-local-pixels of them (local warming that day) and it
           is within --max-offset-departure of the cruise's offset, else the median of those per-tile offsets over the cruise's tiles with paired
           SST (sst_source "paired": same daytime overpass as OLCI), else over all
           cruises in the manifest.

Every pixel without SLSTR SST then gets OSTIA + offset; SLSTR pixels are kept.
OSTIA is interpolated linearly in time (between daily fields) and space onto the
tile's pixels. Tiles store only their centre, so the pixel grid is placed north-up
at the nominal 300 m spacing around it; at OSTIA's 5 km grid (and a smoother
effective resolution) the orientation hardly matters.

Inputs:  a tile manifest (default data/manifests/biogoship_tile_manifest.csv) and
         its tiles; OSTIA subsets per cruise section, downloaded with the
         copernicusmarine toolbox into --l4-dir (same login as the currents).
Outputs: tiles under --out-root (same relative paths as under --tiles-root);
         --out-manifest, the manifest pointing at them with sst_valid_fraction
         updated and sst_slstr_fraction, sst_l4_fraction, sst_offset_k,
         sst_offset_kind added; --log-csv, one row per tile with the offsets and
         a check of OSTIA + offset against SLSTR on paired tiles (printed too).

Usage:
    python scripts/fill_sst_l4.py --dry-run
    python scripts/fill_sst_l4.py
    python scripts/fill_sst_l4.py --manifest-csv data/faire/manifests/biogoship_tile_manifest.csv \\
        --tiles-root data/faire/tiles --out-root data/faire/tiles_sst_l4 \\
        --out-manifest data/faire/manifests/biogoship_tile_manifest_sst_l4.csv \\
        --log-csv data/faire/sst_l4_log.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import tiles  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, from_repo  # noqa: E402

VARIABLE = "analysed_sst"
DATASETS = ("METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2", "METOFFICE-GLO-SST-L4-REP-OBS-SST")


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest_csv, dtype={"sst_source": str})
    tile_table = manifest.drop_duplicates("tile_path").set_index("tile_path")
    tile_table["time"] = pd.to_datetime(tile_table["tile_date"], utc=True).dt.tz_localize(None)
    tile_table["subset"] = subsets(tile_table, args.split_gap_days)
    print(f"{args.manifest_csv}: {len(tile_table):,} tiles in {tile_table['section'].nunique()} sections, "
          f"{tile_table['subset'].nunique()} OSTIA subsets")

    files = download(tile_table, args)
    if args.dry_run:
        return
    import xarray as xr

    opened = {path: xr.open_dataset(path)[VARIABLE] for paths in files.values() for path in paths}

    # Pass 1: OSTIA on each tile, SLSTR where present, the tile's own offset.
    l4_fields, slstr_fields, log = {}, {}, []
    for i, (path, tile) in enumerate(tile_table.iterrows(), start=1):
        with np.load(from_repo(path), allow_pickle=False) as npz:
            sst_norm = npz["data"][-1].astype(np.float32)
            lat_c, lon_c = float(npz["lat_center"]), float(npz["lon_center"])
        slstr = sst_norm * tiles.SST_STD_K + tiles.SST_MEAN_K
        grid = choose_file(files[tile["subset"]], lon_c, opened)
        l4 = ostia_on_tile(opened[grid], tile["time"], lat_c, lon_c)
        both = np.isfinite(slstr) & np.isfinite(l4)
        l4_fields[path], slstr_fields[path] = l4, slstr
        log.append({"tile_path": path, "section": tile["section"], "sst_source": tile.get("sst_source", ""),
                    "slstr_pixels": int(np.isfinite(slstr).sum()), "l4_centre_k": float(l4[tiles.HALF, tiles.HALF]),
                    "slstr_mean_k": float(np.nanmean(slstr)) if np.isfinite(slstr).any() else np.nan,
                    "local_offset_k": float(np.median(slstr[both] - l4[both])) if both.sum() >= args.min_local_pixels else np.nan})  # fmt: skip
        if i % 500 == 0:
            print(f"  OSTIA on {i:,}/{len(tile_table):,} tiles", flush=True)
    log = pd.DataFrame(log).set_index("tile_path")

    # Cruise offsets from paired tiles; the fallback is over all of them.
    paired = log[(log["sst_source"] == "paired") & log["local_offset_k"].notna()]
    cruise_offset = paired.groupby("section")["local_offset_k"].median()
    overall = float(paired["local_offset_k"].median())
    print(f"\nSLSTR (daytime skin) - OSTIA, per tile with paired SST; median offset applied where a tile has too little SLSTR:")
    print(paired.groupby("section")["local_offset_k"].describe(percentiles=[0.25, 0.5, 0.75])[["count", "25%", "50%", "75%"]]
          .round(2).to_string())  # fmt: skip
    print(f"All cruises: {overall:+.2f} K")
    check = paired.assign(pred=paired["l4_centre_k"] + paired["section"].map(cruise_offset), raw=paired["l4_centre_k"])
    for column, label in (("raw", "OSTIA"), ("pred", "OSTIA + cruise offset")):
        err = check[column] - check["slstr_mean_k"]
        print(f"Check on paired tiles, {label} vs tile SLSTR mean: bias {err.mean():+.2f} K, RMSE {np.sqrt((err ** 2).mean()):.2f} K")

    # Pass 2: fill and write.
    out_rows, written = {}, 0
    for path in tile_table.index:
        row = log.loc[path]
        fallback = float(cruise_offset.get(row["section"], overall))
        if np.isfinite(row["local_offset_k"]) and abs(row["local_offset_k"] - fallback) <= args.max_offset_departure:
            offset, kind = row["local_offset_k"], "tile"
        elif np.isfinite(row["local_offset_k"]):  # e.g. cloud-contaminated SLSTR at a cloud edge
            offset, kind = fallback, "cruise (tile offset rejected)"
        else:
            offset, kind = fallback, "cruise" if row["section"] in cruise_offset.index else "all"
        slstr, l4 = slstr_fields[path], l4_fields[path]
        gap = ~np.isfinite(slstr)
        filled = np.where(gap, l4 + offset, slstr).astype(np.float32)
        total = filled.size
        info = {"sst_slstr_fraction": round(float((~gap).sum()) / total, 4),
                "sst_l4_fraction": round(float((gap & np.isfinite(filled)).sum()) / total, 4),
                "sst_valid_fraction": round(float(np.isfinite(filled).sum()) / total, 4),
                "sst_offset_k": round(float(offset), 3), "sst_offset_kind": kind}  # fmt: skip
        new_path = variant_path(path, args)
        out_rows[path] = {**info, "tile_path": new_path}
        destination = from_repo(new_path)
        if destination.exists() and not args.overwrite:
            continue
        with np.load(from_repo(path), allow_pickle=False) as npz:
            contents = {key: npz[key] for key in npz.files}
        contents["data"] = contents["data"].copy()
        contents["data"][-1] = tiles.normalise_sst(filled)
        contents["sst_valid_fraction"] = np.array(info["sst_valid_fraction"], dtype=np.float32)
        contents["sst_l4_fraction"] = np.array(info["sst_l4_fraction"], dtype=np.float32)
        contents["sst_offset_k"] = np.array(offset, dtype=np.float32)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(".tmp.npz")
        np.savez_compressed(temporary, **contents)
        temporary.replace(destination)
        written += 1

    extra = pd.DataFrame.from_dict(out_rows, orient="index")
    variant = manifest.copy()
    for column in ("sst_slstr_fraction", "sst_l4_fraction", "sst_valid_fraction", "sst_offset_k", "sst_offset_kind"):
        variant[column] = variant["tile_path"].map(extra[column])
    source = variant["sst_source"].fillna("none").replace("", "none")
    variant["sst_source"] = np.where(variant["sst_l4_fraction"] <= 0, source,
                                     np.where(variant["sst_slstr_fraction"] > 0, source + "+ostia", "ostia"))  # fmt: skip
    variant["tile_path"] = variant["tile_path"].map(extra["tile_path"])
    args.out_manifest.parent.mkdir(parents=True, exist_ok=True)
    variant.to_csv(args.out_manifest, index=False)
    log.join(extra.drop(columns="tile_path")).to_csv(args.log_csv)
    print(f"\nTiles written: {written:,} (others already there); offset used: "
          + ", ".join(f"{k} {v:,}" for k, v in extra["sst_offset_kind"].value_counts().items()))
    print(f"Tiles still with SST gaps (OSTIA missing, e.g. at the coast): {int((extra['sst_valid_fraction'] < 1).sum()):,}")
    print(f"Wrote {args.out_manifest} and {args.log_csv}")


def subsets(tile_table: pd.DataFrame, split_gap_days: float) -> pd.Series:
    """OSTIA subset of each tile: its section, or with split_gap_days > 0 one per run of tiles
    in a section with no gap longer than that (a long expedition or a time series would
    otherwise need one subset spanning years and an ocean basin)."""
    if split_gap_days <= 0:
        return tile_table["section"].astype(str)
    out = pd.Series(index=tile_table.index, dtype=object)
    for section, group in tile_table.groupby("section"):
        times = group["time"].sort_values()
        run = (times.diff() > pd.Timedelta(days=split_gap_days)).cumsum()
        out[times.index] = [f"{section}_{k:03d}" for k in run]
    return out


def download(tile_table: pd.DataFrame, args: argparse.Namespace) -> dict[str, list[Path]]:
    """One OSTIA subset per cruise section, or per run of it (--split-gap-days), two if it crosses
    180 degrees; returns subset -> files."""
    args.l4_dir.mkdir(parents=True, exist_ok=True)
    files: dict[str, list[Path]] = {}
    for section, group in tile_table.groupby("subset"):
        lon = group["tile_lon"].to_numpy(dtype=float)
        start = (group["time"].min() - pd.Timedelta(days=1)).floor("D")
        end = (group["time"].max() + pd.Timedelta(days=2)).floor("D")
        lat_min = max(float(group["tile_lat"].min()) - args.pad_deg, -90.0)
        lat_max = min(float(group["tile_lat"].max()) + args.pad_deg, 90.0)
        stem = "ostia_" + "".join(ch if ch.isalnum() else "_" for ch in str(section))
        if lon.max() - lon.min() <= 180:
            pieces = [(stem, max(lon.min() - args.pad_deg, -180.0), min(lon.max() + args.pad_deg, 180.0))]
        else:  # crosses the date line: east of 180 and west of it separately
            east, west = lon[lon >= 0], lon[lon < 0]
            pieces = [(stem + "_e", max(east.min() - args.pad_deg, -180.0), 180.0),
                      (stem + "_w", -180.0, min(west.max() + args.pad_deg, 180.0))]  # fmt: skip
        files[section] = []
        for name, lon_min, lon_max in pieces:
            dest = args.l4_dir / f"{name}.nc"
            files[section].append(dest)
            print(f"{str(section):18s} {len(group):5,} tiles | lon {lon_min:7.2f} to {lon_max:7.2f}, lat {lat_min:6.2f} to "
                  f"{lat_max:6.2f} | {start:%Y-%m-%d} to {end:%Y-%m-%d} -> {dest.name}"
                  + (" (on disk)" if dest.exists() else ""), flush=True)  # fmt: skip
            if args.dry_run or dest.exists():
                continue
            fetch(name, lon_min, lon_max, lat_min, lat_max, start, end, args)
    return files


def fetch(name, lon_min, lon_max, lat_min, lat_max, start, end, args) -> None:
    import copernicusmarine

    errors = []
    for dataset in ([args.dataset_id] if args.dataset_id else DATASETS):
        try:
            copernicusmarine.subset(
                dataset_id=dataset, variables=[VARIABLE],
                minimum_longitude=lon_min, maximum_longitude=lon_max, minimum_latitude=lat_min, maximum_latitude=lat_max,
                start_datetime=start.strftime("%Y-%m-%dT%H:%M:%S"), end_datetime=end.strftime("%Y-%m-%dT%H:%M:%S"),
                output_directory=str(args.l4_dir), output_filename=f"{name}.nc", overwrite=True, disable_progress_bar=True,
            )  # fmt: skip
            if not (args.l4_dir / f"{name}.nc").exists():  # the toolbox can return without raising, e.g. no login
                raise RuntimeError("no file written (no Copernicus Marine login?)")
            print(f"{'':18s} from {dataset}", flush=True)
            return
        except Exception as exc:  # e.g. dates outside this dataset; try the next one
            errors.append(f"{dataset}: {exc}")
    raise SystemExit("OSTIA download failed:\n  " + "\n  ".join(errors))


def choose_file(paths: list[Path], lon: float, opened: dict) -> Path:
    for path in paths:
        lons = opened[path]["longitude"].values
        if lons.min() - 0.05 <= lon <= lons.max() + 0.05:
            return path
    return paths[0]


def ostia_on_tile(field, time: pd.Timestamp, lat_c: float, lon_c: float) -> np.ndarray:
    """OSTIA in kelvin on the tile's pixels (north-up grid around the centre), linear in time and space."""
    offsets_km = (np.arange(tiles.TILE_SIZE) - (tiles.TILE_SIZE - 1) / 2) * tiles.PIXEL_KM
    lats = lat_c - offsets_km / 111.195  # row 0 to the north
    lons = lon_c + offsets_km / (111.195 * max(np.cos(np.radians(lat_c)), 0.05))
    lons = (lons + 180.0) % 360.0 - 180.0
    pad = 0.2
    lat_axis = field["latitude"].values
    lat_slice = slice(lat_c - pad, lat_c + pad) if lat_axis[0] < lat_axis[-1] else slice(lat_c + pad, lat_c - pad)
    local = field.sel(time=slice(time - pd.Timedelta(days=1), time + pd.Timedelta(days=1)), latitude=lat_slice,
                      longitude=slice(lons.min() - pad, lons.max() + pad)).load()  # fmt: skip
    if local.sizes["time"] == 0 or local.sizes["latitude"] < 2 or local.sizes["longitude"] < 2:
        return np.full((tiles.TILE_SIZE, tiles.TILE_SIZE), np.nan, dtype=np.float32)
    if local.sizes["time"] == 1:
        at_time = local.isel(time=0)
    else:
        at_time = local.interp(time=np.datetime64(time, "ns"), method="linear")
        if not np.isfinite(at_time.values).any():  # outside the downloaded days: nearest day
            at_time = local.sel(time=np.datetime64(time, "ns"), method="nearest")
    values = at_time.interp(latitude=lats, longitude=lons, method="linear").values.astype(np.float32)
    if not np.isfinite(values).all():  # coast or ice edge: the mean of the surrounding OSTIA cells
        surrounding = float(np.nanmean(at_time.values)) if np.isfinite(at_time.values).any() else np.nan
        values = np.where(np.isfinite(values), values, surrounding)
    return values


def variant_path(path: str, args: argparse.Namespace) -> str:
    root = str(args.tiles_root).rstrip("/") + "/"
    if not path.startswith(root):
        raise SystemExit(f"{path} is not under --tiles-root {args.tiles_root}")
    return str(args.out_root).rstrip("/") + "/" + path[len(root):]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-csv", type=Path, default=MANIFEST_DIR / "biogoship_tile_manifest.csv")
    parser.add_argument("--tiles-root", type=Path, default=Path("data/tiles"), help="As written in the manifest's tile_path.")
    parser.add_argument("--out-root", type=Path, default=Path("data/tiles_sst_l4"))
    parser.add_argument("--out-manifest", type=Path, default=MANIFEST_DIR / "biogoship_tile_manifest_sst_l4.csv")
    parser.add_argument("--log-csv", type=Path, default=REPO_ROOT / "data" / "sst_l4_log.csv")
    parser.add_argument("--l4-dir", type=Path, default=REPO_ROOT / "data" / "sst_l4")
    parser.add_argument("--dataset-id", help=f"Copernicus Marine dataset. Default: {DATASETS[0]}, else {DATASETS[1]}.")
    parser.add_argument("--pad-deg", type=float, default=0.5)
    parser.add_argument("--split-gap-days", type=float, default=0.0,
                        help="Split a section's OSTIA subset wherever its tiles are more than this many days apart (0 = one per section).")
    parser.add_argument("--min-local-pixels", type=int, default=200,
                        help="SLSTR pixels a tile needs for its own offset (of 1,764). Default 200.")
    parser.add_argument("--max-offset-departure", type=float, default=2.0,
                        help="A tile's own offset further than this (K) from its cruise's is not trusted. Default 2.")
    parser.add_argument("--overwrite", action="store_true", help="Rewrite tiles that already exist under --out-root.")
    parser.add_argument("--dry-run", action="store_true", help="Print the OSTIA subsets; download and write nothing.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
