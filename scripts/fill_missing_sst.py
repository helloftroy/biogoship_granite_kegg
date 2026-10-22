#!/usr/bin/env python3
"""Fill tiles that have no SST with the closest-in-time SLSTR SST over the same spot.

A quarter of tiles have no SST: the SST orbit file paired with their OLCI scene
held no acceptable pixel there. Most are Sentinel-3B scenes from 2018, when no
Sentinel-3B SST product existed and the Sentinel-3A orbit file used instead was
over a different part of the ocean at that moment. This step looks for any SLSTR
SST overpass (either satellite) whose footprint contains the tile, within
--max-hours of the OLCI scene, and fills the tile's SST channel from the first
one that has acceptable SST there (quality level 4 or 5).

Order of preference: daytime overpasses (local solar time 06-18, like the OLCI
scene itself) closest in time first, then night-time ones; night SST lacks the
daytime warming of the skin layer, so it is only used when no daytime pass works.

The SST channel is resampled onto the tile's own OLCI pixels exactly as for
paired SST (nearest neighbour). The tile file is updated in place (its SST channel
was empty, so nothing is lost) and gains: sst_source = "filled", sst_fill_product,
sst_fill_gap_hours, sst_fill_night. Every tile's outcome goes to
data/sst_fill_log.csv. Re-running skips tiles already filled.

SST archives already in data/products/sst/ are used directly; others are fetched
as the single L2P file (credentials needed, see cdse_credentials.sh.example).
OLCI geolocation comes from data/products/olci/ archives.

Usage:
    python scripts/fill_missing_sst.py --dry-run          # search only; counts and download size
    python scripts/fill_missing_sst.py --max-hours 48
"""

from __future__ import annotations

import argparse
import shutil
import sys
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import cdse, tiles  # noqa: E402
from biogoship.tables import (  # noqa: E402
    REPO_ROOT,
    append_csv,
    from_repo,
    latest_tile_log_rows,
    parse_utc,
    read_csv,
    read_tile_log,
)

LOG_FIELDS = ["tile_path", "status", "sst_product", "satellite", "gap_hours", "night", "sst_valid_fraction", "candidates"]
MIN_FILLED_FRACTION = 0.5  # accept a candidate only if it covers at least half the tile with acceptable SST


def main() -> None:
    args = parse_args()
    log_path = args.data_dir / "sst_fill_log.csv"
    done = {row["tile_path"] for row in read_csv(log_path) if row["status"] == "filled"} if log_path.exists() else set()
    todo = tiles_without_sst(args, done)
    print(f"Tiles without SST: {len(todo):,} to fill ({len(done):,} filled earlier)")
    if not todo:
        return

    # Catalogue searches are saved, so a re-run (or the real run after a dry run) does not
    # repeat them: hundreds of searches in a row trip the catalogue's rate limit (HTTP 429).
    cache_path = args.data_dir / "sst_fill_candidates.csv"
    cached = load_candidate_cache(cache_path, args.max_hours)
    lock = threading.Lock()

    def search(tile: dict) -> list[dict] | None:
        if tile["tile_path"] in cached:
            return cached[tile["tile_path"]]
        window = timedelta(hours=args.max_hours)
        try:
            found = cdse.search_products(cdse.SLSTR_WST, tile["time"] - window, tile["time"] + window,
                                         point=(tile["lon"], tile["lat"]), retries=args.search_retries)  # fmt: skip
        except Exception as exc:  # noqa: BLE001 - leave this tile for the next run
            print(f"  search failed for {Path(tile['tile_path']).name}: {exc}", flush=True)
            return None
        finally:
            time.sleep(args.search_pause)
        ranked = rank_candidates(tile, found)
        with lock:
            save_candidates(cache_path, tile["tile_path"], args.max_hours, ranked)
        return ranked

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for tile, candidates in zip(todo, pool.map(search, todo)):
            tile["candidates"] = candidates
    unsearched = [tile for tile in todo if tile["candidates"] is None]
    todo = [tile for tile in todo if tile["candidates"] is not None]
    if unsearched:
        print(f"Catalogue search failed for {len(unsearched):,} tiles (rate limit or network); re-run later to retry them. "
              f"Searches done so far are saved in {cache_path.name}.")  # fmt: skip
    unique = {c["product_id"]: c for tile in todo for c in tile["candidates"][args.first_try : args.first_try + args.max_tries]}
    local = {pid for pid, c in unique.items() if local_archive(args, c["product_name"]) is not None}
    no_candidate = sum(1 for tile in todo if not tile["candidates"])
    print(f"Tiles with at least one SST overpass within +/-{args.max_hours:g} h: {len(todo) - no_candidate:,}; none: {no_candidate:,}")
    print(f"SST products among the first {args.max_tries} candidates per tile: {len(unique):,} ({len(local):,} already on disk; "
          f"up to {sum(float(c['size_bytes'] or 0) for p, c in unique.items() if p not in local) / 1e9:,.0f} GB to fetch, "
          "less since only the L2P file is fetched and tiles stop at the first usable one)")  # fmt: skip
    if args.dry_run:
        return

    filler = Filler(args)
    todo.sort(key=lambda tile: tile["time"])  # neighbouring tiles tend to share SST orbits; keeps the cache warm
    counts = {"filled": 0, "no_usable_sst": 0, "no_candidate": 0, "error": 0}
    for i, tile in enumerate(todo, start=1):
        row = {"tile_path": tile["tile_path"], "candidates": len(tile["candidates"])}
        try:
            row.update(filler.fill(tile))
        except Exception as exc:  # noqa: BLE001 - one bad tile must not stop the run
            row.update(status="error", sst_product=str(exc)[:200])
        counts[row["status"]] = counts.get(row["status"], 0) + 1
        append_csv(log_path, [row], LOG_FIELDS)
        if i % 25 == 0 or i == len(todo):
            print(f"  {i:,}/{len(todo):,}: " + ", ".join(f"{n:,} {k}" for k, n in counts.items() if n), flush=True)
    filler.cleanup()
    print(f"Log: {log_path}\nNext: build_tile_manifest.py, pool_tile_features.py and embed_tiles.py to refresh downstream tables.")


def tiles_without_sst(args: argparse.Namespace, done: set[str]) -> list[dict]:
    log = latest_tile_log_rows(read_tile_log(args.data_dir / "tile_log.csv"))
    out = []
    for row in log:
        if row["status"] != "valid" or row["tile_path"] in done:
            continue
        if row.get("sst_valid_fraction") and float(row["sst_valid_fraction"]) > 0:
            continue
        path = from_repo(row["tile_path"])
        if not path.exists():
            continue
        with np.load(path, allow_pickle=False) as tile:
            if np.isfinite(tile["data"][-1]).any():
                continue
            out.append({
                "tile_path": row["tile_path"], "time": parse_utc(str(tile["date"])), "granule": str(tile["granule_id"]),
                "lat": float(tile["lat_center"]), "lon": float(tile["lon_center"]),
            })  # fmt: skip
    return out


CANDIDATE_FIELDS = ["tile_path", "max_hours", "product_id", "product_name", "satellite", "start_utc", "end_utc",
                    "size_bytes", "gap_hours", "night"]  # fmt: skip


def load_candidate_cache(path: Path, max_hours: float) -> dict[str, list[dict]]:
    """tile_path -> ranked candidates from earlier searches with the same time window.

    A tile searched with no result is stored as a row with an empty product_id.
    """
    cache: dict[str, list[dict]] = {}
    if not path.exists():
        return cache
    for row in read_csv(path):
        if float(row["max_hours"]) != max_hours:
            continue
        found = cache.setdefault(row["tile_path"], [])
        if row["product_id"]:
            found.append({**row, "gap_hours": float(row["gap_hours"]), "night": row["night"] == "1"})
    return cache


def save_candidates(path: Path, tile_path: str, max_hours: float, ranked: list[dict]) -> None:
    rows = [{**c, "tile_path": tile_path, "max_hours": max_hours, "gap_hours": f"{c['gap_hours']:.3f}",
             "night": int(c["night"])} for c in ranked]  # fmt: skip
    append_csv(path, rows or [{"tile_path": tile_path, "max_hours": max_hours}], CANDIDATE_FIELDS)


def rank_candidates(tile: dict, products: list[dict]) -> list[dict]:
    """Daytime passes closest in time first, then night-time ones."""
    ranked = []
    for product in products:
        start, end = parse_utc(product["start_utc"]), parse_utc(product["end_utc"])
        middle = start + (end - start) / 2
        gap_hours = abs((middle - tile["time"]).total_seconds()) / 3600
        local_hour = (middle.hour + middle.minute / 60 + tile["lon"] / 15) % 24
        night = not 6 <= local_hour < 18
        ranked.append({**product, "gap_hours": gap_hours, "night": night})
    return sorted(ranked, key=lambda c: (c["night"], c["gap_hours"]))


def local_archive(args: argparse.Namespace, product_name: str) -> Path | None:
    for kind in ("sst", "olci"):
        path = args.products_dir / kind / f"{product_name}.zip"
        if path.exists():
            return path
    return None


class Filler:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.tokens = cdse.TokenProvider()
        self.work = args.data_dir / "work" / "sst_fill"
        self.work.mkdir(parents=True, exist_ok=True)
        self.swaths: dict[str, tiles.SstSwath] = {}
        self.geo: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def fill(self, tile: dict) -> dict:
        if not tile["candidates"]:
            return {"status": "no_candidate"}
        lat, lon = self.tile_grid(tile)
        for candidate in tile["candidates"][self.args.first_try : self.args.first_try + self.args.max_tries]:
            sst = self.swath(candidate).resample(lat, lon)
            fraction = float(np.isfinite(sst).mean())
            if fraction >= MIN_FILLED_FRACTION:
                self.write(tile, sst, candidate)
                return {
                    "status": "filled", "sst_product": candidate["product_name"], "satellite": candidate["satellite"],
                    "gap_hours": f"{candidate['gap_hours']:.2f}", "night": int(candidate["night"]),
                    "sst_valid_fraction": f"{fraction:.4f}",
                }  # fmt: skip
        return {"status": "no_usable_sst"}

    def tile_grid(self, tile: dict) -> tuple[np.ndarray, np.ndarray]:
        """Per-pixel lat/lon of the tile, from its OLCI scene's geolocation."""
        if tile["granule"] not in self.geo:
            if len(self.geo) >= 2:
                self.geo.pop(next(iter(self.geo)))
            self.geo[tile["granule"]] = tiles.load_geo(self.olci_geo_file(tile["granule"]))
        lat, lon = self.geo[tile["granule"]]
        row, col, km = tiles.locate_pixel(lat, lon, tile["lat"], tile["lon"])
        origin = tiles.window_origin(row, col, lat.shape)
        if origin is None or km > 1.0:
            raise RuntimeError(f"tile centre not found in its scene ({km:.2f} km)")
        window = (slice(origin[0], origin[0] + tiles.TILE_SIZE), slice(origin[1], origin[1] + tiles.TILE_SIZE))
        return lat[window], lon[window]

    def olci_geo_file(self, granule: str) -> Path:
        dest = self.work / granule / "geo_coordinates.nc"
        if not dest.exists():
            dest.parent.mkdir(parents=True, exist_ok=True)
            archive = local_archive(self.args, granule)
            if archive is None:
                raise RuntimeError(f"OLCI archive for {granule} not on disk; run download_products.py first")
            extract(archive, "geo_coordinates.nc", dest)
        return dest

    def swath(self, candidate: dict) -> tiles.SstSwath:
        pid = candidate["product_id"]
        if pid not in self.swaths:
            if len(self.swaths) >= 3:
                self.swaths.pop(next(iter(self.swaths)))
            self.swaths[pid] = tiles.SstSwath(self.l2p_file(candidate))
        return self.swaths[pid]

    def l2p_file(self, candidate: dict) -> Path:
        folder = self.work / candidate["product_name"]
        existing = [p for p in folder.glob("*.nc") if tiles.is_sst_l2p_file(p.name)] if folder.exists() else []
        if existing:
            return existing[0]
        folder.mkdir(parents=True, exist_ok=True)
        archive = local_archive(self.args, candidate["product_name"])
        if archive is not None:
            with zipfile.ZipFile(archive) as z:
                name = next(Path(m).name for m in z.namelist() if tiles.is_sst_l2p_file(Path(m).name))
            extract(archive, name, folder / name)
            return folder / name
        sizes = cdse.list_files(candidate["product_id"], candidate["product_name"])
        name = next(n for n in sizes if tiles.is_sst_l2p_file(n))
        print(f"  fetching {candidate['product_name'][:31]}... ({sizes[name] / 1e6:,.0f} MB)", flush=True)
        cdse.download(cdse.file_url(candidate["product_id"], candidate["product_name"], name), folder / name, self.tokens,
                      expected_bytes=sizes[name])  # fmt: skip
        self.trim_cache(keep=folder)
        return folder / name

    def trim_cache(self, keep: Path) -> None:
        folders = sorted((p for p in self.work.iterdir() if p.is_dir() and "_SL_2_WST_" in p.name), key=lambda p: p.stat().st_mtime)
        for old in folders[: max(0, len(folders) - self.args.sst_cache)]:
            if old != keep:
                shutil.rmtree(old, ignore_errors=True)

    def write(self, tile: dict, sst_kelvin: np.ndarray, candidate: dict) -> None:
        path = from_repo(tile["tile_path"])
        with np.load(path, allow_pickle=False) as npz:
            contents = {key: npz[key] for key in npz.files}
        contents["data"] = contents["data"].copy()
        contents["data"][-1] = tiles.normalise_sst(sst_kelvin)
        contents["sst_valid_fraction"] = np.array(np.isfinite(sst_kelvin).mean(), dtype=np.float32)
        contents.update(
            sst_source=np.array("filled"),
            sst_fill_product=np.array(candidate["product_name"]),
            sst_fill_gap_hours=np.array(candidate["gap_hours"], dtype=np.float32),
            sst_fill_night=np.array(bool(candidate["night"])),
        )
        temporary = path.with_name(path.stem + ".tmp.npz")
        np.savez_compressed(temporary, **contents)
        temporary.replace(path)

    def cleanup(self) -> None:
        if not self.args.keep_work:
            shutil.rmtree(self.work, ignore_errors=True)


def extract(archive: Path, filename: str, dest: Path) -> None:
    with zipfile.ZipFile(archive) as z:
        member = next(m for m in z.namelist() if Path(m).name == filename)
        with z.open(member) as src, dest.open("wb") as out:
            shutil.copyfileobj(src, out)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data")
    parser.add_argument("--products-dir", type=Path, default=REPO_ROOT / "data" / "products")
    parser.add_argument("--max-hours", type=float, default=48.0, help="Largest time gap to the OLCI scene. Default 48 h.")
    parser.add_argument("--max-tries", type=int, default=4, help="SST overpasses to try per tile before giving up.")
    parser.add_argument(
        "--first-try",
        type=int,
        default=0,
        help="Skip this many of each tile's best-ranked overpasses (e.g. 4 after a run with --max-tries 4).",
    )
    parser.add_argument("--sst-cache", type=int, default=6, help="Downloaded SST files kept for reuse.")
    parser.add_argument("--workers", type=int, default=2, help="Parallel catalogue searches.")
    parser.add_argument("--search-pause", type=float, default=1.0, help="Seconds each worker waits between searches.")
    parser.add_argument("--search-retries", type=int, default=10, help="Retries per search on rate limits (back-off up to 15 min).")
    parser.add_argument("--dry-run", action="store_true", help="Search and report; change nothing.")
    parser.add_argument("--keep-work", action="store_true", help="Keep extracted and downloaded files.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
