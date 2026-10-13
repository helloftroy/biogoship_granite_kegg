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

Tracks run at least --track-days each way (default 7), longer than the +/-3 day
overpass window, so that search_trajectory_products.py can look for overpasses
further out along the track.

`status` in advected_positions.csv:
  ok           position found
  stalled      the particle never moved: no current data at the site (land or coast)
  left_domain  the track ends before the overpass time (it left the currents subset)
  no_track     the site was not tracked (no currents file for its section)

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
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, format_utc, parse_utc, read_csv, write_csv  # noqa: E402

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

    tracks: dict = {}
    trajectories = []
    for box in advection.section_boxes(sites, site_products, scene_times, margin_deg=0.0, pad_days=0.0):
        currents_path = args.currents_dir / advection.section_filename(box["section"])
        section_links = [link for site in box["sites"] for link in links_by_site[site["site_id"]]]
        if not currents_path.exists():
            print(f"{box['section']}: no currents file {currents_path.name}; run download_currents.py", flush=True)
            continue

        started = time.time()
        # Long enough for every overpass already matched, and for the wider trajectory search.
        max_hours = max(args.track_days * 24.0, max(abs(float(link["time_gap_hours"])) for link in section_links) + 1.0)
        releases = [
            {"id": s["site_id"], "lat": float(s["lat"]), "lon": float(s["lon"]), "time": parse_utc(s["time_utc"])}
            for s in box["sites"]
        ]
        currents = advection.load_currents(currents_path, args.depth_m)
        section_tracks = advection.advect(currents, releases, max_hours, dt_minutes=args.dt_minutes)
        tracks.update(section_tracks)

        ends = []
        for site in box["sites"]:
            sampled = parse_utc(site["time_utc"])
            track = section_tracks[site["site_id"]]
            ends.append((track[0, 0], track[-1, 0]))
            for hours, lat, lon in track:
                trajectories.append(
                    {
                        "site_id": site["site_id"],
                        "hours_from_sample": f"{hours:.3f}",
                        "time_utc": format_utc(sampled + timedelta(hours=float(hours))),
                        "lat": f"{lat:.5f}",
                        "lon": f"{lon:.5f}",
                    }
                )
        full = sum(1 for start, end in ends if start <= -max_hours + 1 and end >= max_hours - 1)
        short_back = sorted(-start / 24 for start, _end in ends if start > -max_hours + 1)
        short_fwd = sorted(end / 24 for _start, end in ends if end < max_hours - 1)
        print(
            f"{box['section']:12s} {len(releases):4d} sites tracked +/-{max_hours / 24:.1f} days; "
            f"{full:4d} tracks complete in both directions ({time.time() - started:.0f}s)",
            flush=True,
        )
        for label, short in (("backward", short_back), ("forward", short_fwd)):
            if short:
                print(f"{'':12s} {len(short):4d} {label} tracks stop early, after {short[0]:.1f} to {short[-1]:.1f} days "
                      "(ran out of currents in time or space)")  # fmt: skip

    write_csv(args.out_dir / "trajectories.csv", trajectories, TRAJECTORY_FIELDS)
    positions = advection.advected_position_rows(sites, site_products, tracks)
    write_csv(args.manifest_dir / "advected_positions.csv", positions, advection.POSITION_FIELDS)

    section_of = {site["site_id"]: site["section"] for site in sites}
    print("\nDisplacement at overpass time, for the site x scene pairs resolved:")
    for section in sorted({section_of[row["site_id"]] for row in positions}):
        moved = sorted(
            float(row["displacement_km"])
            for row in positions
            if row["status"] == "ok" and section_of[row["site_id"]] == section
        )
        if moved:
            print(f"  {section:12s} {len(moved):5d} pairs | median {moved[len(moved) // 2]:6.1f} km, max {moved[-1]:6.1f} km")
    counts = defaultdict(int)
    for row in positions:
        counts[row["status"]] += 1
    print("Pairs by status: " + ", ".join(f"{count:,} {status}" for status, count in sorted(counts.items())))
    print(f"Wrote {args.out_dir / 'trajectories.csv'}")
    print(f"Wrote {args.manifest_dir / 'advected_positions.csv'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--currents-dir", type=Path, default=REPO_ROOT / "data" / "currents")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "advection")
    parser.add_argument(
        "--track-days",
        type=float,
        default=7.0,
        help="Track each site at least this many days backward and forward. Default 7.",
    )
    parser.add_argument("--depth-m", type=float, default=15.0, help="Current level to use: 0 or 15 m. Default 15.")
    parser.add_argument("--dt-minutes", type=float, default=20.0, help="Integration time step.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
