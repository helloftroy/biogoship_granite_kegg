#!/usr/bin/env python3
"""Join samples to their satellite tiles and report coverage.

Reads data/tile_log.csv (written by download_tiles.py) and writes
data/manifests/biogoship_tile_manifest.csv: one row per (sample, clear tile).
`position` says whether the tile is centred on the sampling position (`fixed`) or
on where the water was at overpass time (`advected`). `is_best` marks, for each
sample and position, the tile closest in time to sampling; filter on it for a
one-tile-per-sample table. `same_scene_as_best_fixed` marks the advected tile cut
from the same overpass as the best fixed tile, for a like-for-like comparison. The column names follow the manifests of the
Florida project (tile_path, tile_date, tile_lat, tile_lon, match_day_gap,
match_distance_km), so the Granite embedding step can read it unchanged.

Also writes data/manifests/site_tile_status.csv: per site and position, how many
candidate scenes were clear, cloudy or did not cover the position, and a `flag`
naming the problem for sites left without a tile.

Usage:
    python scripts/build_tile_manifest.py
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import (  # noqa: E402
    MANIFEST_DIR,
    REPO_ROOT,
    from_repo,
    haversine_km,
    latest_tile_log_rows,
    read_csv,
    write_csv,
)

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
    "position",
    "sample_time",
    "sample_lat",
    "sample_lon",
    "match_day_gap",
    "match_hour_gap",
    "match_distance_km",
    "offset_km",
    "valid_fraction",
    "sst_valid_fraction",
    "is_best",
    "same_scene_as_best_fixed",
    "displacement_km",
    "granule_id",
    *SAMPLE_COLUMNS,
]


def main() -> None:
    args = parse_args()
    samples = read_csv(args.manifest_dir / "biogoship_samples.csv")
    products = {row["product_id"]: row for row in read_csv(args.manifest_dir / "products.csv")}
    log = latest_tile_log_rows(read_csv(args.tile_log))

    tiles_by_site: dict[str, list[dict[str, str]]] = defaultdict(list)
    missing = 0
    for row in log:
        if row["status"] != "valid":
            continue
        if not from_repo(row["tile_path"]).exists():
            missing += 1
            continue
        tiles_by_site[row["site_id"]].append(row)
    best_fixed_scene = {}
    for site_id, site_tiles in tiles_by_site.items():
        site_tiles.sort(key=lambda row: abs(float(row["time_gap_hours"])))
        seen = set()
        for tile in site_tiles:  # closest in time first, so the first of each position is its best
            tile["is_best"] = "0" if tile["position"] in seen else "1"
            seen.add(tile["position"])
            if tile["position"] == "fixed" and tile["is_best"] == "1":
                best_fixed_scene[site_id] = tile["product_id"]

    rows = []
    for sample in samples:
        for tile in tiles_by_site.get(sample["site_id"], []):
            gap_hours = float(tile["time_gap_hours"])
            lat, lon = float(sample["lat"]), float(sample["lon"])
            distance_km = haversine_km(lat, lon, float(tile["tile_lat"]), float(tile["tile_lon"]))
            same_scene = tile["position"] == "advected" and best_fixed_scene.get(sample["site_id"]) == tile["product_id"]
            rows.append(
                {
                    "tile_path": tile["tile_path"],
                    "tile_date": products[tile["product_id"]]["start_utc"],
                    "tile_lat": tile["tile_lat"],
                    "tile_lon": tile["tile_lon"],
                    "sample_id": sample["sample_id"],
                    "site_id": sample["site_id"],
                    "position": tile["position"],
                    "sample_time": sample["collection_time_utc"],
                    "sample_lat": sample["lat"],
                    "sample_lon": sample["lon"],
                    "match_day_gap": f"{abs(gap_hours) / 24:.3f}",
                    "match_hour_gap": f"{gap_hours:.2f}",
                    "match_distance_km": f"{distance_km:.3f}",
                    "offset_km": tile.get("offset_km") or "0.000",
                    "valid_fraction": tile["valid_fraction"],
                    "sst_valid_fraction": tile["sst_valid_fraction"],
                    "is_best": tile["is_best"],
                    "same_scene_as_best_fixed": "1" if same_scene else "0",
                    "displacement_km": f"{distance_km:.3f}" if tile["position"] == "advected" else "0.000",
                    "granule_id": tile["product_name"],
                    **{column: sample.get(column, "") for column in SAMPLE_COLUMNS},
                }
            )
    write_csv(args.out_csv, rows, FIELDS)
    write_csv(args.out_csv.with_name("site_tile_status.csv"), site_status(samples, log), STATUS_FIELDS)
    summarise(samples, log, tiles_by_site, rows, missing, args)


STATUS_FIELDS = [
    "site_id",
    "section",
    "position",
    "flag",
    "n_scenes",
    "n_valid",
    "n_cloudy",
    "n_edge",
    "n_out_of_swath",
    "n_no_position",
    "best_valid_gap_hours",
    "best_cloudy_valid_fraction",
]


def site_status(samples: list[dict[str, str]], log: list[dict[str, str]]) -> list[dict]:
    """Per site and position: how its candidate scenes turned out, and a flag saying what to do about it.

    flag: ok                  at least one clear tile
          all_cloudy          every scene covers the position but none is clear enough
          outside_all_scenes  the position is off the edge of, or outside, every scene
          cloudy_or_outside   a mix of the two above
          no_position         back-tracking gave no advected position for any scene
    """
    section = {sample["site_id"]: sample["section"] for sample in samples}
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in log:
        groups[(row["site_id"], row["position"])].append(row)
    out = []
    for (site_id, position), rows in sorted(groups.items()):
        n = {status: sum(1 for row in rows if row["status"] == status) for status in
             ("valid", "cloudy", "edge", "out_of_swath", "no_position")}  # fmt: skip
        outside = n["edge"] + n["out_of_swath"]
        if n["valid"]:
            flag = "ok"
        elif n["no_position"] == len(rows):
            flag = "no_position"
        elif n["cloudy"] and not outside:
            flag = "all_cloudy"
        elif outside and not n["cloudy"]:
            flag = "outside_all_scenes"
        else:
            flag = "cloudy_or_outside"
        valid_gaps = [abs(float(row["time_gap_hours"])) for row in rows if row["status"] == "valid"]
        cloudy = [float(row["valid_fraction"]) for row in rows if row["status"] == "cloudy"]
        out.append(
            {
                "site_id": site_id,
                "section": section.get(site_id, ""),
                "position": position,
                "flag": flag,
                "n_scenes": len(rows),
                "n_valid": n["valid"],
                "n_cloudy": n["cloudy"],
                "n_edge": n["edge"],
                "n_out_of_swath": n["out_of_swath"],
                "n_no_position": n["no_position"],
                "best_valid_gap_hours": f"{min(valid_gaps):.2f}" if valid_gaps else "",
                "best_cloudy_valid_fraction": f"{max(cloudy):.3f}" if cloudy else "",
            }
        )
    return out


def summarise(samples, log, tiles_by_site, rows, missing, args) -> None:
    statuses = defaultdict(int)
    for row in log:
        statuses[row["status"]] += 1
    print("Site x scene outcomes: " + ", ".join(f"{count:,} {status}" for status, count in sorted(statuses.items())))
    if missing:
        print(f"WARNING: {missing:,} logged tiles are not on disk under {REPO_ROOT} and were left out.")
    def sites_with(position: str) -> set[str]:
        return {site for site, site_tiles in tiles_by_site.items() if any(t["position"] == position for t in site_tiles)}

    fixed, advected = sites_with("fixed"), sites_with("advected")
    print(f"\n{'section':12s} {'samples':>8s} {'fixed tile':>11s} {'advected tile':>14s} {'both':>6s}")
    for section in sorted({sample["section"] for sample in samples}):
        in_section = [sample["site_id"] for sample in samples if sample["section"] == section]
        counts = [sum(1 for site in in_section if site in group) for group in (fixed, advected, fixed & advected)]
        print(f"{section:12s} {len(in_section):8d} {counts[0]:11d} {counts[1]:14d} {counts[2]:6d}")
    matched = {row["sample_id"] for row in rows}
    print(f"\nSamples with at least one clear tile: {len(matched):,} of {len(samples):,}")
    print(f"Manifest rows (sample x tile): {len(rows):,}")
    flags: dict[tuple[str, str], int] = defaultdict(int)
    for row in site_status(samples, log):
        flags[(row["position"], row["flag"])] += 1
    print("Sites by position and flag: " + ", ".join(f"{n:,} {p}/{flag}" for (p, flag), n in sorted(flags.items())))
    print(f"Wrote {args.out_csv}")
    print(f"Wrote {args.out_csv.with_name('site_tile_status.csv')}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--tile-log", type=Path, default=REPO_ROOT / "data" / "tile_log.csv")
    parser.add_argument("--out-csv", type=Path, default=MANIFEST_DIR / "biogoship_tile_manifest.csv")
    return parser.parse_args()


if __name__ == "__main__":
    main()
