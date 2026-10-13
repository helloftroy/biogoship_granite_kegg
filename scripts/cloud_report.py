#!/usr/bin/env python3
"""Explain the tile yield: how many sites would get a tile under other settings.

Reads data/tile_log.csv only (no scenes), so it is instant. For each position it
reports, per section:

  * sites with a clear tile as the clear-fraction threshold is lowered
  * sites with a clear tile as the time window is narrowed
  * which quality flags spoil the cloudy windows

It cannot predict the effect of --max-offset-km (using a nearby clear window);
that needs a re-run of download_tiles.py with --recheck.

Usage:
    python scripts/cloud_report.py
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, latest_tile_log_rows, read_csv  # noqa: E402

THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.1]
WINDOWS_H = [6, 12, 24, 48, 72, 240]


def main() -> None:
    args = parse_args()
    section = {site["site_id"]: site["section"] for site in read_csv(args.manifest_dir / "biogoship_sites.csv")}
    log = latest_tile_log_rows(read_csv(args.tile_log))
    for position in sorted({row["position"] for row in log}):
        rows = [row for row in log if row["position"] == position and row["status"] in ("valid", "cloudy")]
        by_site: dict[str, list[tuple[float, float]]] = defaultdict(list)  # site -> [(clear fraction, |gap| hours)]
        for row in rows:
            by_site[row["site_id"]].append((float(row["valid_fraction"]), abs(float(row["time_gap_hours"]))))
        sections = sorted({section[site] for site in by_site})
        totals = {name: sum(1 for site in by_site if section[site] == name) for name in sections}
        print(f"\n=== position: {position} ===  ({len(rows):,} site x scene checks, {len(by_site):,} sites)")

        print("\nSites with a tile, by minimum clear fraction (any scene in the window):")
        print(f"{'section':12s} {'sites':>6s} " + " ".join(f"{'>=' + format(t, '.1f'):>7s}" for t in THRESHOLDS))
        for name in [*sections, "ALL"]:
            sites = [s for s in by_site if name == "ALL" or section[s] == name]
            best = [max(fraction for fraction, _gap in by_site[s]) for s in sites]
            print(f"{name:12s} {len(sites):6d} " + " ".join(f"{sum(1 for b in best if b >= t):7d}" for t in THRESHOLDS))

        print(f"\nSites with a tile (clear fraction >= {args.min_valid}), by largest allowed time gap:")
        print(f"{'section':12s} {'sites':>6s} " + " ".join(f"{'<=' + str(h) + 'h':>7s}" for h in WINDOWS_H))
        for name in [*sections, "ALL"]:
            sites = [s for s in by_site if name == "ALL" or section[s] == name]
            counts = [
                sum(1 for s in sites if any(f >= args.min_valid and gap <= hours for f, gap in by_site[s]))
                for hours in WINDOWS_H
            ]
            print(f"{name:12s} {len(sites):6d} " + " ".join(f"{count:7d}" for count in counts))

        flags: dict[str, float] = defaultdict(float)
        cloudy = [row for row in rows if row["status"] == "cloudy"]
        for row in cloudy:
            for item in filter(None, row.get("reject_flags", "").split(";")):
                name, _, fraction = item.partition("=")
                flags[name] += float(fraction)
        print(f"\nMean share of pixels carrying each flag, over the {len(cloudy):,} cloudy windows (top three flags per window):")
        for name, total in sorted(flags.items(), key=lambda item: -item[1]):
            print(f"  {name:16s} {total / max(len(cloudy), 1):.2f}")

        offsets = [float(row["offset_km"]) for row in log if row["position"] == position and row["status"] == "valid" and row.get("offset_km")]
        shifted = [o for o in offsets if o > 0.5]
        if shifted:
            shifted.sort()
            print(f"\nTiles taken from a nearby clear window: {len(shifted):,} of {len(offsets):,}; "
                  f"offset median {shifted[len(shifted) // 2]:.1f} km, max {shifted[-1]:.1f} km")  # fmt: skip


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--tile-log", type=Path, default=REPO_ROOT / "data" / "tile_log.csv")
    parser.add_argument("--min-valid", type=float, default=0.8)
    return parser.parse_args()


if __name__ == "__main__":
    main()
