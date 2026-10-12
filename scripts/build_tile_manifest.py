#!/usr/bin/env python3
"""Join samples to their satellite tiles and report coverage.

Reads data/tile_log.csv (written by download_tiles.py) and writes
data/manifests/biogoship_tile_manifest.csv: one row per (sample, clear tile).
`is_best` marks, for each sample, the tile closest in time to sampling; filter on
it for a one-tile-per-sample table. The column names follow the manifests of the
Florida project (tile_path, tile_date, tile_lat, tile_lon, match_day_gap,
match_distance_km), so the Granite embedding step can read it unchanged.

Usage:
    python scripts/build_tile_manifest.py
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, from_repo, haversine_km, read_csv, write_csv  # noqa: E402

SAMPLE_COLUMNS = [
    "biosample",
    "section",
    "depth_m",
    "sample_method",
    "metagenome_runs",
    "amplicon_runs",
    "temperature_c",
    "phosphate_umol_kg",
    "nitrate_nitrite_umol_l",
]
FIELDS = [
    "tile_path",
    "tile_date",
    "tile_lat",
    "tile_lon",
    "sample_id",
    "site_id",
    "sample_time",
    "sample_lat",
    "sample_lon",
    "match_day_gap",
    "match_hour_gap",
    "match_distance_km",
    "valid_fraction",
    "sst_valid_fraction",
    "is_best",
    "granule_id",
    *SAMPLE_COLUMNS,
]


def main() -> None:
    args = parse_args()
    samples = read_csv(args.manifest_dir / "biogoship_samples.csv")
    products = {row["product_id"]: row for row in read_csv(args.manifest_dir / "products.csv")}
    log = read_csv(args.tile_log)

    tiles_by_site: dict[str, list[dict[str, str]]] = defaultdict(list)
    missing = 0
    for row in log:
        if row["status"] != "valid":
            continue
        if not from_repo(row["tile_path"]).exists():
            missing += 1
            continue
        tiles_by_site[row["site_id"]].append(row)
    for site_tiles in tiles_by_site.values():
        site_tiles.sort(key=lambda row: abs(float(row["time_gap_hours"])))

    rows = []
    for sample in samples:
        for rank, tile in enumerate(tiles_by_site.get(sample["site_id"], [])):
            gap_hours = float(tile["time_gap_hours"])
            lat, lon = float(sample["lat"]), float(sample["lon"])
            rows.append(
                {
                    "tile_path": tile["tile_path"],
                    "tile_date": products[tile["product_id"]]["start_utc"],
                    "tile_lat": tile["tile_lat"],
                    "tile_lon": tile["tile_lon"],
                    "sample_id": sample["sample_id"],
                    "site_id": sample["site_id"],
                    "sample_time": sample["collection_time_utc"],
                    "sample_lat": sample["lat"],
                    "sample_lon": sample["lon"],
                    "match_day_gap": f"{abs(gap_hours) / 24:.3f}",
                    "match_hour_gap": f"{gap_hours:.2f}",
                    "match_distance_km": f"{haversine_km(lat, lon, float(tile['tile_lat']), float(tile['tile_lon'])):.3f}",
                    "valid_fraction": tile["valid_fraction"],
                    "sst_valid_fraction": tile["sst_valid_fraction"],
                    "is_best": "1" if rank == 0 else "0",
                    "granule_id": tile["product_name"],
                    **{column: sample.get(column, "") for column in SAMPLE_COLUMNS},
                }
            )
    write_csv(args.out_csv, rows, FIELDS)
    summarise(samples, log, tiles_by_site, rows, missing, args)


def summarise(samples, log, tiles_by_site, rows, missing, args) -> None:
    statuses = defaultdict(int)
    for row in log:
        statuses[row["status"]] += 1
    print("Site x scene outcomes: " + ", ".join(f"{count:,} {status}" for status, count in sorted(statuses.items())))
    if missing:
        print(f"WARNING: {missing:,} logged tiles are not on disk under {REPO_ROOT} and were left out.")
    print(f"\n{'section':12s} {'samples':>8s} {'with tile':>10s} {'metagenomes with tile':>22s}")
    for section in sorted({sample["section"] for sample in samples}):
        in_section = [sample for sample in samples if sample["section"] == section]
        with_tile = [sample for sample in in_section if sample["site_id"] in tiles_by_site]
        metagenomes = sum(1 for sample in with_tile if sample["metagenome_runs"])
        print(f"{section:12s} {len(in_section):8d} {len(with_tile):10d} {metagenomes:22d}")
    matched = {row["sample_id"] for row in rows}
    print(f"\nSamples with at least one clear tile: {len(matched):,} of {len(samples):,}")
    print(f"Manifest rows (sample x tile): {len(rows):,}")
    print(f"Wrote {args.out_csv}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--tile-log", type=Path, default=REPO_ROOT / "data" / "tile_log.csv")
    parser.add_argument("--out-csv", type=Path, default=MANIFEST_DIR / "biogoship_tile_manifest.csv")
    return parser.parse_args()


if __name__ == "__main__":
    main()
