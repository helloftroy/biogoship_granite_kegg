#!/usr/bin/env python3
"""Download surface currents for the back-tracking step.

Product: Copernicus Marine MULTIOBS_GLO_PHY_MYNRT_015_003 (GlobCurrent total
velocity from altimetry, wind-driven Ekman currents and drifters; 0.25 degree).
One NetCDF subset per cruise section, covering the sites plus a margin and the
full span of sampling and overpass times. The multi-year dataset is used where it
reaches, the near-real-time one after its end (2026 cruises). A section that
crosses the 180 degree meridian is fetched in two pieces and stored in 0 to 360
degrees longitude.

    data/currents/currents_<section>.nc     variables uo, vo at 0 m and 15 m

Needs a free Copernicus Marine account (a different account from the Copernicus
Data Space one used for Sentinel-3); see cdse_credentials.sh.example.

Usage:
    python scripts/download_currents.py --dry-run
    python scripts/download_currents.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import advection  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, parse_utc, read_csv  # noqa: E402


def main() -> None:
    args = parse_args()
    boxes = load_boxes(args.manifest_dir, args.margin_deg, args.pad_days)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for box in boxes:
        dest = args.out_dir / advection.section_filename(box["section"])
        state = "to download"
        if dest.exists():
            state = "on disk, covers the request" if covers(dest, box) else "on disk but too short: downloading again"
        print(
            f"{box['section']:12s} {len(box['sites']):4d} sites | lon {box['lon_min']:.1f} to {box['lon_max']:.1f}, "
            f"lat {box['lat_min']:.1f} to {box['lat_max']:.1f} | {box['start']:%Y-%m-%d} to {box['end']:%Y-%m-%d} "
            f"-> {dest.name} ({state})",
            flush=True,
        )
        dataset = args.dataset_id or advection.currents_dataset(box["start"], box["end"])
        print(f"{'':12s} dataset {dataset}" + (" | crosses 180 degrees: two pieces, stored in 0-360" if box["wraps"] else ""))
        if args.dry_run or (state == "on disk, covers the request" and not args.overwrite):
            continue
        if not box["wraps"]:
            fetch(dataset, box, box["lon_min"], box["lon_max"], args.out_dir, dest.name)
            continue
        west, east = dest.with_suffix(".west.nc"), dest.with_suffix(".east.nc")
        fetch(dataset, box, box["lon_min"], 180.0, args.out_dir, west.name)
        fetch(dataset, box, -180.0, box["lon_max"] - 360.0, args.out_dir, east.name)
        join_across_dateline(west, east, dest)
    print(f"Currents in {args.out_dir}")


def fetch(dataset: str, box: dict, lon_min: float, lon_max: float, out_dir: Path, name: str) -> None:
    import copernicusmarine

    copernicusmarine.subset(
        dataset_id=dataset,
        variables=advection.CURRENT_VARIABLES,
        minimum_longitude=lon_min,
        maximum_longitude=lon_max,
        minimum_latitude=box["lat_min"],
        maximum_latitude=box["lat_max"],
        start_datetime=box["start"].strftime("%Y-%m-%dT%H:%M:%S"),
        end_datetime=box["end"].strftime("%Y-%m-%dT%H:%M:%S"),
        output_directory=str(out_dir),
        output_filename=name,
        overwrite=True,
        disable_progress_bar=True,
    )


def join_across_dateline(west: Path, east: Path, dest: Path) -> None:
    """One file in 0-360 longitude from the pieces west (up to 180) and east (from -180) of the dateline."""
    import xarray as xr

    with xr.open_dataset(west) as w, xr.open_dataset(east) as e:
        e = e.assign_coords(longitude=e["longitude"] + 360.0)
        w = w.sel(longitude=w["longitude"] < e["longitude"].min())  # drop the shared 180 degree column
        joined = xr.concat([w, e], dim="longitude").sortby("longitude").load()
    joined.to_netcdf(dest)
    west.unlink()
    east.unlink()


def covers(path: Path, box: dict) -> bool:
    """Does an existing currents file span the box in time and space (to within one grid step)?"""
    import numpy as np
    import xarray as xr

    with xr.open_dataset(path) as ds:
        times = ds["time"].values
        lat, lon = ds["latitude"].values, ds["longitude"].values
    day = np.timedelta64(1, "D")
    start = np.datetime64(box["start"].replace(tzinfo=None), "ns")
    end = np.datetime64(box["end"].replace(tzinfo=None), "ns")
    return bool(
        times.min() <= start + day
        and times.max() >= end - day
        and lat.min() <= box["lat_min"] + 0.5
        and lat.max() >= box["lat_max"] - 0.5
        and lon.min() <= box["lon_min"] + 0.5
        and lon.max() >= box["lon_max"] - 0.5
    )


def load_boxes(manifest_dir: Path, margin_deg: float, pad_days: float) -> list[dict]:
    sites = read_csv(manifest_dir / "biogoship_sites.csv")
    site_products = read_csv(manifest_dir / "site_products.csv")
    scene_times = {row["product_id"]: parse_utc(row["start_utc"]) for row in read_csv(manifest_dir / "products.csv")}
    return advection.section_boxes(sites, site_products, scene_times, margin_deg, pad_days)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "currents")
    parser.add_argument(
        "--dataset-id", default=None, help="Force one currents dataset (default: multi-year, near-real-time past its end)."
    )
    parser.add_argument("--margin-deg", type=float, default=5.0, help="Room around the sites for water to drift into.")
    parser.add_argument(
        "--pad-days",
        type=float,
        default=8.0,
        help="Days of currents before the first and after the last sampling or overpass time. "
        "Must exceed backtrack_sites.py --track-days (7). Default 8.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Re-download subsets that already exist.")
    parser.add_argument("--dry-run", action="store_true", help="Print the subsets; download nothing.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
