#!/usr/bin/env python3
"""Download surface currents for the back-tracking step.

Product: Copernicus Marine MULTIOBS_GLO_PHY_MYNRT_015_003 (GlobCurrent total
velocity from altimetry, wind-driven Ekman currents and drifters; 0.25 degree).
One NetCDF subset per cruise section, covering the sites plus a margin and the
full span of sampling and overpass times:

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
        print(
            f"{box['section']:12s} {len(box['sites']):4d} sites | lon {box['lon_min']:.1f} to {box['lon_max']:.1f}, "
            f"lat {box['lat_min']:.1f} to {box['lat_max']:.1f} | {box['start']:%Y-%m-%d} to {box['end']:%Y-%m-%d} "
            f"-> {dest.name}{' (exists)' if dest.exists() else ''}",
            flush=True,
        )
        if args.dry_run or (dest.exists() and not args.overwrite):
            continue
        import copernicusmarine

        copernicusmarine.subset(
            dataset_id=args.dataset_id,
            variables=advection.CURRENT_VARIABLES,
            minimum_longitude=box["lon_min"],
            maximum_longitude=box["lon_max"],
            minimum_latitude=box["lat_min"],
            maximum_latitude=box["lat_max"],
            start_datetime=box["start"].strftime("%Y-%m-%dT%H:%M:%S"),
            end_datetime=box["end"].strftime("%Y-%m-%dT%H:%M:%S"),
            output_directory=str(args.out_dir),
            output_filename=dest.name,
            overwrite=True,
            disable_progress_bar=True,
        )
    print(f"Currents in {args.out_dir}")


def load_boxes(manifest_dir: Path, margin_deg: float, pad_days: float) -> list[dict]:
    sites = read_csv(manifest_dir / "biogoship_sites.csv")
    site_products = read_csv(manifest_dir / "site_products.csv")
    scene_times = {row["product_id"]: parse_utc(row["start_utc"]) for row in read_csv(manifest_dir / "products.csv")}
    return advection.section_boxes(sites, site_products, scene_times, margin_deg, pad_days)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "currents")
    parser.add_argument("--dataset-id", default=advection.CURRENTS_DATASET)
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
