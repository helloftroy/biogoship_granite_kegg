#!/usr/bin/env python3
"""Download Sentinel-3 data and cut one Granite tile per (site, clear overpass).

Works through data/manifests/products.csv one OLCI scene at a time:

  1. Fetch only geo_coordinates.nc and wqsf.nc (~60 MB) and test every site the
     scene covers: is the 42x42 window around the site inside the swath, and is it
     at least --min-valid clear open ocean?
  2. Only if some site passes, fetch the 16 reflectance bands (~280 MB) and the
     paired SST orbit file, and write one .npz tile per passing site.
  3. Record the outcome for every (site, scene) in data/tile_log.csv, then delete
     the scene's files.

Re-running resumes: scenes already in data/product_log.csv are skipped. Scenes
that errored are retried on the next run.

Each tile .npz holds:
  data            float32 (17, 42, 42)  normalised; channels 0-15 OLCI, 16 SST
  valid_mask      bool (42, 42)         clear open-ocean pixels
  date            OLCI sensing start, ISO-8601 UTC
  lat_center, lon_center                 centre pixel position
  granule_id      source OL_2_WFR product name
  site_id, site_time, time_gap_hours, valid_fraction, sst_valid_fraction

Credentials: see cdse_credentials.sh.example.

Usage:
    python scripts/download_tiles.py --max-products 3     # smoke test
    python scripts/download_tiles.py
    python scripts/download_tiles.py --download-mode zip  # fallback: whole product archives
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import cdse, tiles  # noqa: E402
from biogoship.tables import (  # noqa: E402
    MANIFEST_DIR,
    REPO_ROOT,
    append_csv,
    format_utc,
    parse_utc,
    read_csv,
    repo_relative,
)

DATA_DIR = REPO_ROOT / "data"
MAX_PIXEL_KM = 1.0  # a site farther than this from its nearest pixel is outside the swath

TILE_LOG_FIELDS = [
    "site_id",
    "product_id",
    "product_name",
    "status",  # valid | cloudy | edge | out_of_swath
    "valid_fraction",
    "reject_flags",
    "time_gap_hours",
    "pixel_km",
    "row",
    "col",
    "tile_lat",
    "tile_lon",
    "sst_valid_fraction",
    "tile_path",
]
PRODUCT_LOG_FIELDS = ["product_id", "product_name", "status", "n_sites", "n_valid", "seconds", "message"]


def main() -> None:
    args = parse_args()
    sites = {row["site_id"]: row for row in read_csv(args.manifest_dir / "biogoship_sites.csv")}
    products = read_csv(args.manifest_dir / "products.csv")
    links: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_csv(args.manifest_dir / "site_products.csv"):
        links[row["product_id"]].append(row)

    product_log = args.data_dir / "product_log.csv"
    done = {row["product_id"] for row in read_csv(product_log) if row["status"] == "done"} if product_log.exists() else set()
    todo = [product for product in products if product["product_id"] not in done]
    if args.max_products:
        todo = todo[: args.max_products]
    print(f"OLCI scenes: {len(products):,} in catalogue, {len(done):,} already done, {len(todo):,} to process now")
    print(f"Download mode: {args.download_mode}; min valid fraction: {args.min_valid}; SST: {'off' if args.no_sst else 'on'}")

    worker = TileWorker(args)
    totals = {"valid": 0, "errors": 0}
    for i, product in enumerate(todo, start=1):
        started = time.time()
        print(f"[{i}/{len(todo)}] {product['product_name']}", flush=True)
        product_sites = [(sites[link["site_id"]], link) for link in links[product["product_id"]]]
        try:
            rows = worker.process(product, product_sites)
        except Exception as exc:  # noqa: BLE001 - one bad scene must not stop the run
            totals["errors"] += 1
            print(f"  ERROR: {exc}", flush=True)
            log_product(product_log, product, "error", len(product_sites), 0, started, str(exc)[:300])
            worker.cleanup(product)
            continue
        append_csv(args.data_dir / "tile_log.csv", rows, TILE_LOG_FIELDS)
        n_valid = sum(1 for row in rows if row["status"] == "valid")
        totals["valid"] += n_valid
        log_product(product_log, product, "done", len(product_sites), n_valid, started, "")
        worker.cleanup(product)
        counts = defaultdict(int)
        for row in rows:
            counts[row["status"]] += 1
        summary = ", ".join(f"{count} {status}" for status, count in sorted(counts.items()))
        print(f"  {summary} ({time.time() - started:.0f}s; {totals['valid']:,} tiles this run)", flush=True)

    print(f"\nFinished: {totals['valid']:,} tiles written, {totals['errors']:,} scenes failed (re-run to retry them).")
    print(f"Tiles: {args.data_dir / 'tiles'}\nLog:   {args.data_dir / 'tile_log.csv'}")
    if totals["errors"] and totals["errors"] == len(todo):
        raise SystemExit("Every scene failed; check credentials and network before re-running.")


def log_product(path: Path, product: dict, status: str, n_sites: int, n_valid: int, started: float, message: str) -> None:
    row = {
        "product_id": product["product_id"],
        "product_name": product["product_name"],
        "status": status,
        "n_sites": n_sites,
        "n_valid": n_valid,
        "seconds": f"{time.time() - started:.0f}",
        "message": message,
    }
    append_csv(path, [row], PRODUCT_LOG_FIELDS)


class TileWorker:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.tokens = cdse.TokenProvider()
        self.work_dir = args.data_dir / "work"
        self.sst_dir = self.work_dir / "sst_cache"
        self.tiles_dir = args.data_dir / "tiles"
        self.sst_swath: tiles.SstSwath | None = None
        self.sst_swath_id = ""

    # -- fetching -------------------------------------------------------------
    def product_dir(self, product_name: str) -> Path:
        return self.work_dir / product_name

    def fetch(self, product_id: str, product_name: str, filenames: list[str], dest_dir: Path) -> None:
        """Make the named files of a product available in dest_dir."""
        missing = [name for name in filenames if not (dest_dir / name).exists()]
        if not missing:
            return
        dest_dir.mkdir(parents=True, exist_ok=True)
        if self.args.download_mode == "zip" or self.local_zip(product_name) is not None:
            self._fetch_from_zip(product_id, product_name, missing, dest_dir)
            return
        sizes = cdse.list_files(product_id, product_name)
        absent = [name for name in missing if name not in sizes]
        if absent:
            raise RuntimeError(f"product has no file(s) {absent}")

        def fetch_one(name: str) -> None:
            cdse.download(
                cdse.file_url(product_id, product_name, name),
                dest_dir / name,
                self.tokens,
                expected_bytes=sizes[name],
                timeout=self.args.timeout,
                retries=self.args.retries,
                retry_sleep=self.args.retry_sleep,
            )

        megabytes = sum(sizes[name] for name in missing) / 1e6
        print(f"  downloading {len(missing)} file(s), {megabytes:,.0f} MB", flush=True)
        with ThreadPoolExecutor(max_workers=self.args.connections) as pool:
            list(pool.map(fetch_one, missing))

    def local_zip(self, product_name: str) -> Path | None:
        """Archive already fetched by download_products.py, if there is one."""
        if self.args.products_dir is None:
            return None
        matches = list(self.args.products_dir.glob(f"*/{product_name}.zip"))
        return matches[0] if matches else None

    def _fetch_from_zip(self, product_id: str, product_name: str, filenames: list[str], dest_dir: Path) -> None:
        zip_path = self.local_zip(product_name) or self.work_dir / f"{product_name}.zip"
        if not zip_path.exists():
            print("  downloading whole product archive", flush=True)
            cdse.download(
                cdse.product_zip_url(product_id),
                zip_path,
                self.tokens,
                timeout=self.args.timeout,
                retries=self.args.retries,
                retry_sleep=self.args.retry_sleep,
            )
        with zipfile.ZipFile(zip_path) as archive:
            members = {Path(member).name: member for member in archive.namelist()}
            for name in filenames:
                if name not in members:
                    raise RuntimeError(f"archive has no file {name}")
                with archive.open(members[name]) as src, (dest_dir / name).open("wb") as out:
                    shutil.copyfileobj(src, out)

    def cleanup(self, product: dict[str, str]) -> None:
        if self.args.keep_products:
            return
        shutil.rmtree(self.product_dir(product["product_name"]), ignore_errors=True)
        (self.work_dir / f"{product['product_name']}.zip").unlink(missing_ok=True)

    # -- sea surface temperature -----------------------------------------------
    def load_sst(self, product: dict[str, str]) -> tiles.SstSwath | None:
        """The SST orbit paired with this scene. Recent orbit files stay cached: neighbouring scenes share them."""
        sst_id, sst_name = product.get("sst_product_id", ""), product.get("sst_product_name", "")
        if self.args.no_sst or not sst_id:
            return None
        if self.sst_swath is not None and self.sst_swath_id == sst_id:
            return self.sst_swath
        sst_product_dir = self.sst_dir / sst_name
        try:
            l2p = next((p for p in sst_product_dir.glob("*.nc") if tiles.is_sst_l2p_file(p.name)), None)
            if l2p is None:
                filename = self._sst_filename(sst_id, sst_name)
                self.fetch(sst_id, sst_name, [filename], sst_product_dir)
                (self.work_dir / f"{sst_name}.zip").unlink(missing_ok=True)
                l2p = sst_product_dir / filename
            self.sst_swath, self.sst_swath_id = tiles.SstSwath(l2p), sst_id
            self._trim_sst_cache(keep=sst_product_dir)
            return self.sst_swath
        except Exception as exc:  # noqa: BLE001 - SST is optional; keep the OLCI tile
            print(f"  SST unavailable ({exc}); tiles get a NaN SST channel", flush=True)
            return None

    def _sst_filename(self, sst_id: str, sst_name: str) -> str:
        """Name of the single GHRSST L2P data file inside an SST product."""
        if self.args.download_mode == "zip" or self.local_zip(sst_name) is not None:
            zip_path = self.local_zip(sst_name) or self.work_dir / f"{sst_name}.zip"
            if not zip_path.exists():
                cdse.download(
                    cdse.product_zip_url(sst_id),
                    zip_path,
                    self.tokens,
                    timeout=self.args.timeout,
                    retries=self.args.retries,
                    retry_sleep=self.args.retry_sleep,
                )
            with zipfile.ZipFile(zip_path) as archive:
                names = [Path(member).name for member in archive.namelist()]
        else:
            names = list(cdse.list_files(sst_id, sst_name))
        l2p_names = [name for name in names if tiles.is_sst_l2p_file(name)]
        if not l2p_names:
            raise RuntimeError("no L2P file in SST product")
        return l2p_names[0]

    def _trim_sst_cache(self, keep: Path) -> None:
        cached = sorted((p for p in self.sst_dir.iterdir() if p.is_dir()), key=lambda p: p.stat().st_mtime)
        for old in cached[: max(0, len(cached) - self.args.sst_cache)]:
            if old != keep:
                shutil.rmtree(old, ignore_errors=True)

    # -- one scene --------------------------------------------------------------
    def process(self, product: dict[str, str], product_sites: list[tuple[dict, dict]]) -> list[dict]:
        product_id, product_name = product["product_id"], product["product_name"]
        product_dir = self.product_dir(product_name)
        self.fetch(product_id, product_name, tiles.OLCI_CHECK_FILES, product_dir)
        lat, lon = tiles.load_geo(product_dir / "geo_coordinates.nc")
        wqsf = tiles.load_wqsf(product_dir / "wqsf.nc")

        rows, passing = [], []
        for site, link in product_sites:
            row = {
                "site_id": site["site_id"],
                "product_id": product_id,
                "product_name": product_name,
                "time_gap_hours": link["time_gap_hours"],
            }
            rows.append(row)
            pixel_row, pixel_col, pixel_km = tiles.locate_pixel(lat, lon, float(site["lat"]), float(site["lon"]))
            row["pixel_km"] = f"{pixel_km:.3f}"
            if pixel_km > MAX_PIXEL_KM:
                row["status"] = "out_of_swath"
                continue
            row["row"], row["col"] = pixel_row, pixel_col
            origin = tiles.window_origin(pixel_row, pixel_col, lat.shape)
            if origin is None:
                row["status"] = "edge"
                continue
            r0, c0 = origin
            window = (slice(r0, r0 + tiles.TILE_SIZE), slice(c0, c0 + tiles.TILE_SIZE))
            valid = tiles.ocean_mask(wqsf[window])
            row["valid_fraction"] = f"{valid.mean():.4f}"
            if valid.mean() < self.args.min_valid:
                row["status"] = "cloudy"
                row["reject_flags"] = tiles.dominant_reject_flags(wqsf[window])
                continue
            passing.append((site, row, origin, valid, lat[window], lon[window]))

        if not passing:
            return rows

        self.fetch(product_id, product_name, tiles.OLCI_BAND_FILES, product_dir)
        sst_swath = self.load_sst(product)
        sensed = parse_utc(product["start_utc"])
        out_dir = self.tiles_dir / f"{sensed.year:04d}" / f"{sensed.month:02d}"
        out_dir.mkdir(parents=True, exist_ok=True)
        for site, row, (r0, c0), valid, tile_lat, tile_lon in passing:
            olci = tiles.read_olci_window(product_dir, r0, c0)
            sst = sst_swath.resample(tile_lat, tile_lon) if sst_swath is not None else None
            data = tiles.stack_tile(olci, sst)
            centre = tiles.HALF
            tile_path = out_dir / f"tile_{product['satellite']}_{sensed:%Y%m%dT%H%M%S}_{site['site_id']}.npz"
            np.savez_compressed(
                tile_path,
                data=data,
                valid_mask=valid,
                date=np.array(format_utc(sensed)),
                lat_center=np.array(tile_lat[centre, centre], dtype=np.float32),
                lon_center=np.array(tile_lon[centre, centre], dtype=np.float32),
                granule_id=np.array(product_name),
                site_id=np.array(site["site_id"]),
                site_time=np.array(site["time_utc"]),
                time_gap_hours=np.array(float(row["time_gap_hours"]), dtype=np.float32),
                valid_fraction=np.array(valid.mean(), dtype=np.float32),
                sst_valid_fraction=np.array(np.isfinite(data[-1]).mean(), dtype=np.float32),
            )
            row.update(
                status="valid",
                tile_lat=f"{float(tile_lat[centre, centre]):.6f}",
                tile_lon=f"{float(tile_lon[centre, centre]):.6f}",
                sst_valid_fraction=f"{np.isfinite(data[-1]).mean():.4f}",
                tile_path=repo_relative(tile_path),
            )
        return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR, help="Holds tiles/, work/ and the logs.")
    parser.add_argument("--min-valid", type=float, default=0.8, help="Minimum clear open-ocean fraction of a tile.")
    parser.add_argument("--max-products", type=int, default=0, help="Stop after this many scenes (0 = all).")
    parser.add_argument(
        "--download-mode",
        choices=["files", "zip"],
        default="files",
        help="'files' fetches only the files needed; 'zip' fetches whole product archives.",
    )
    parser.add_argument(
        "--products-dir",
        type=Path,
        default=DATA_DIR / "products",
        help="Archives from download_products.py; used instead of downloading when present.",
    )
    parser.add_argument("--no-sst", action="store_true", help="Skip SST; the 17th channel is NaN.")
    parser.add_argument("--connections", type=int, default=4, help="Parallel downloads (CDSE allows 4 per account).")
    parser.add_argument("--sst-cache", type=int, default=2, help="SST orbit files kept on disk for reuse.")
    parser.add_argument("--keep-products", action="store_true", help="Keep downloaded product files after tiling.")
    parser.add_argument("--timeout", type=float, default=1800.0, help="Socket timeout per download, seconds.")
    parser.add_argument("--retries", type=int, default=6)
    parser.add_argument("--retry-sleep", type=float, default=30.0)
    return parser.parse_args()


if __name__ == "__main__":
    main()
