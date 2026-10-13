#!/usr/bin/env python3
"""Work out where each site's water was at the time of every satellite overpass.

For each cruise section, one particle per site is released at the sampling
position and time and advected backward and forward with the surface currents
(OceanParcels, fourth-order Runge-Kutta, no diffusion). Each (site, scene) pair
then gets the particle's position at the scene's sensing time: before sampling
if the overpass came first, after it if the overpass came later.

Inputs:  data/currents/currents_<section>.nc   (scripts/download_currents.py)
Outputs: data/advection/trajectories.csv       hourly track per site
         data/manifests/advected_positions.csv one row per (site, scene)

`status` in advected_positions.csv:
  ok           position found
  stalled      the particle never moved: no current data at the site (land or coast)
  left_domain  the track left the currents subset before the overpass time
  no_currents  no currents file for the section

Usage:
    python scripts/backtrack_sites.py
    python scripts/backtrack_sites.py --depth-m 0
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
import warnings
from collections import defaultdict
from datetime import timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import advection  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, format_utc, haversine_km, parse_utc, read_csv, write_csv  # noqa: E402

POSITION_FIELDS = [
    "site_id",
    "product_id",
    "time_gap_hours",
    "site_lat",
    "site_lon",
    "advected_lat",
    "advected_lon",
    "displacement_km",
    "path_km",
    "status",
]
TRAJECTORY_FIELDS = ["site_id", "hours_from_sample", "time_utc", "lat", "lon"]


def main() -> None:
    args = parse_args()
    warnings.filterwarnings("ignore", module="parcels")
    logging.getLogger("parcels").setLevel(logging.ERROR)

    sites = read_csv(args.manifest_dir / "biogoship_sites.csv")
    site_products = read_csv(args.manifest_dir / "site_products.csv")
    scene_times = {row["product_id"]: parse_utc(row["start_utc"]) for row in read_csv(args.manifest_dir / "products.csv")}
    links_by_site: dict[str, list[dict[str, str]]] = defaultdict(list)
    for link in site_products:
        links_by_site[link["site_id"]].append(link)

    positions, trajectories = [], []
    for box in advection.section_boxes(sites, site_products, scene_times, margin_deg=0.0, pad_days=0.0):
        currents_path = args.currents_dir / advection.section_filename(box["section"])
        section_links = [(site, link) for site in box["sites"] for link in links_by_site[site["site_id"]]]
        if not currents_path.exists():
            print(f"{box['section']}: no currents file {currents_path.name}; run download_currents.py", flush=True)
            positions += [position_row(site, link, None, "no_currents") for site, link in section_links]
            continue

        started = time.time()
        max_hours = max(abs(float(link["time_gap_hours"])) for _site, link in section_links) + 1.0
        releases = [
            {"id": s["site_id"], "lat": float(s["lat"]), "lon": float(s["lon"]), "time": parse_utc(s["time_utc"])}
            for s in box["sites"]
        ]
        currents = advection.load_currents(currents_path, args.depth_m)
        tracks = advection.advect(currents, releases, max_hours, dt_minutes=args.dt_minutes)

        for site in box["sites"]:
            sampled = parse_utc(site["time_utc"])
            for hours, lat, lon in tracks[site["site_id"]]:
                trajectories.append(
                    {
                        "site_id": site["site_id"],
                        "hours_from_sample": f"{hours:.3f}",
                        "time_utc": format_utc(sampled + timedelta(hours=float(hours))),
                        "lat": f"{lat:.5f}",
                        "lon": f"{lon:.5f}",
                    }
                )
        section_rows = [position_row(site, link, tracks[site["site_id"]], "") for site, link in section_links]
        positions += section_rows
        moved = [float(row["displacement_km"]) for row in section_rows if row["status"] == "ok"]
        median = sorted(moved)[len(moved) // 2] if moved else float("nan")
        print(
            f"{box['section']:12s} {len(releases):4d} sites, {len(section_rows):5d} site x scene pairs, "
            f"{len(moved):5d} ok | displacement median {median:.1f} km, max {max(moved, default=float('nan')):.1f} km "
            f"({time.time() - started:.0f}s)",
            flush=True,
        )

    write_csv(args.out_dir / "trajectories.csv", trajectories, TRAJECTORY_FIELDS)
    write_csv(args.manifest_dir / "advected_positions.csv", positions, POSITION_FIELDS)
    counts = defaultdict(int)
    for row in positions:
        counts[row["status"]] += 1
    print("Pairs by status: " + ", ".join(f"{count:,} {status}" for status, count in sorted(counts.items())))
    print(f"Wrote {args.out_dir / 'trajectories.csv'}")
    print(f"Wrote {args.manifest_dir / 'advected_positions.csv'}")


def position_row(site: dict[str, str], link: dict[str, str], track, status: str) -> dict:
    lat, lon, hours = float(site["lat"]), float(site["lon"]), float(link["time_gap_hours"])
    row = {
        "site_id": site["site_id"],
        "product_id": link["product_id"],
        "time_gap_hours": link["time_gap_hours"],
        "site_lat": site["lat"],
        "site_lon": site["lon"],
        "status": status,
    }
    if track is None:
        return row
    where = advection.position_at(track, hours)
    if where is None:
        row["status"] = "left_domain"
        return row
    displacement = haversine_km(lat, lon, where[0], where[1])
    stalled = displacement == 0.0 and abs(hours) >= 1.0
    row.update(
        advected_lat=f"{where[0]:.6f}",
        advected_lon=f"{where[1]:.6f}",
        displacement_km=f"{displacement:.3f}",
        path_km=f"{advection.path_length_km(track, hours):.3f}",
        status="stalled" if stalled else "ok",
    )
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--currents-dir", type=Path, default=REPO_ROOT / "data" / "currents")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "advection")
    parser.add_argument("--depth-m", type=float, default=15.0, help="Current level to use: 0 or 15 m. Default 15.")
    parser.add_argument("--dt-minutes", type=float, default=20.0, help="Integration time step.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
