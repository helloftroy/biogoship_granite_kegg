#!/usr/bin/env python3
"""Find more Sentinel-3 scenes for sites that still lack an advected tile.

The first search (search_products.py) only looked within a few days of sampling,
at the sampling position. This one follows each site's water along its track
(data/advection/trajectories.csv) further out in time: for every day between the
window already searched and --max-days, it asks the catalogue for OLCI scenes
covering the position the water occupied on that day.

By default only sites with no clear advected tile so far are searched, so the
window can be widened in stages (say 5 days, then 7) and each stage only pays for
the sites still missing.

It appends to data/manifests/site_products.csv (rows marked source=trajectory) and
products.csv, and rewrites advected_positions.csv to cover the new pairs. Then:

    download_products.py   fetches the new scenes
    download_tiles.py      tiles them, at the advected position only

No credentials are needed.

Usage:
    python scripts/search_trajectory_products.py --max-days 6
    python scripts/search_trajectory_products.py --max-days 6 --all-sites
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import search_products  # noqa: E402
from biogoship import advection, cdse  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, latest_tile_log_rows, parse_utc, read_csv, write_csv  # noqa: E402

SITE_PRODUCT_FIELDS = [*search_products.SITE_PRODUCT_FIELDS, "source"]


def main() -> None:
    args = parse_args()
    sites = read_csv(args.manifest_dir / "biogoship_sites.csv")
    site_products = read_csv(args.manifest_dir / "site_products.csv")
    products = read_csv(args.manifest_dir / "products.csv")
    if not args.trajectories_csv.exists():
        raise SystemExit(f"{args.trajectories_csv} not found: run scripts/backtrack_sites.py first.")
    tracks = advection.load_tracks(args.trajectories_csv)

    searched_days = searched_window_days(sites, site_products)
    todo = sites_to_search(sites, site_products, tracks, args)
    jobs = day_queries(todo, tracks, searched_days, args.max_days)
    print(f"Sites to search: {len(todo):,}; catalogue queries: {len(jobs):,} (one per site and day, out to +/-{args.max_days:g} days)")

    def run(job: dict) -> list[dict[str, str]]:
        return cdse.search_products(cdse.OLCI_WFR, job["start"], job["end"], point=(job["lon"], job["lat"]))

    known_pairs = {(row["site_id"], row["product_id"]) for row in site_products}
    known_products = {row["product_id"] for row in products}
    new_pairs: list[dict] = []
    new_olci: dict[str, dict] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for done, (job, found) in enumerate(zip(jobs, pool.map(run, jobs)), start=1):
            for product in found:
                pair = (job["site"]["site_id"], product["product_id"])
                if pair in known_pairs:
                    continue
                known_pairs.add(pair)
                gap_hours = (parse_utc(product["start_utc"]) - parse_utc(job["site"]["time_utc"])).total_seconds() / 3600
                new_pairs.append({**product, "site_id": pair[0], "time_gap_hours": f"{gap_hours:.2f}", "source": "trajectory"})
                if product["product_id"] not in known_products:
                    new_olci.setdefault(product["product_id"], product)
            if done % 200 == 0 or done == len(jobs):
                print(f"  {done:,}/{len(jobs):,} queries; {len(new_pairs):,} new site x scene pairs, {len(new_olci):,} new scenes", flush=True)

    if not new_pairs:
        print("No new scenes found; nothing changed.")
        return

    sst = search_products.search_sst(list(new_olci.values())) if new_olci else []
    for product in sorted(new_olci.values(), key=lambda row: row["start_utc"]):
        gaps = [abs(float(r["time_gap_hours"])) for r in new_pairs if r["product_id"] == product["product_id"]]
        products.append({**product, "n_sites": len(gaps), "min_abs_gap_hours": f"{min(gaps):.2f}", **search_products.pair_sst(product, sst)})
    products.sort(key=lambda row: row["start_utc"])
    for row in site_products:
        row.setdefault("source", "site")
    site_products += new_pairs
    site_products.sort(key=lambda row: (row["site_id"], abs(float(row["time_gap_hours"]))))

    write_csv(args.manifest_dir / "site_products.csv", site_products, SITE_PRODUCT_FIELDS)
    write_csv(args.manifest_dir / "products.csv", products, search_products.PRODUCT_FIELDS)
    positions = advection.advected_position_rows(sites, site_products, tracks)
    write_csv(args.manifest_dir / "advected_positions.csv", positions, advection.POSITION_FIELDS)
    summarise(todo, new_pairs, new_olci, products, sst)


def searched_window_days(sites: list[dict[str, str]], site_products: list[dict[str, str]]) -> dict[str, float]:
    """Per section, how many days around sampling the first search already covered."""
    section = {site["site_id"]: site["section"] for site in sites}
    days: dict[str, float] = defaultdict(float)
    for row in site_products:
        if row.get("source", "site") == "site":
            name = section[row["site_id"]]
            days[name] = max(days[name], math.ceil(abs(float(row["time_gap_hours"])) / 24.0))
    return days


def sites_to_search(sites, site_products, tracks, args) -> list[dict[str, str]]:
    tracked = [site for site in sites if site["site_id"] in tracks]
    if args.all_sites:
        return tracked
    if not args.tile_log.exists():
        raise SystemExit(
            f"{args.tile_log} not found. Tile the scenes already downloaded first (download_tiles.py), so that only "
            "sites still missing an advected tile are searched; or pass --all-sites."
        )
    log = latest_tile_log_rows(read_csv(args.tile_log))
    advected = [row for row in log if row["position"] == "advected"]
    if not advected:
        raise SystemExit(
            "The tile log has no advected-position results yet. Run download_tiles.py for the advected position "
            "first, or pass --all-sites to search for every tracked site."
        )
    have_tile = {row["site_id"] for row in advected if row["status"] == "valid"}
    return [site for site in tracked if site["site_id"] not in have_tile]


def day_queries(sites, tracks, searched_days: dict[str, float], max_days: float) -> list[dict]:
    """One catalogue query per site and day beyond the window already searched, at the water's position that day."""
    jobs = []
    for site in sites:
        sampled = parse_utc(site["time_utc"])
        first_day = int(searched_days.get(site["section"], 0))
        for day in range(first_day, int(math.ceil(max_days))):
            for sign in (-1, 1):
                lo, hi = sorted((sign * day * 24.0, sign * min(day + 1, max_days) * 24.0))
                where = advection.position_at(tracks[site["site_id"]], (lo + hi) / 2)
                if where is None:  # the track does not reach this day
                    continue
                jobs.append(
                    {
                        "site": site,
                        "lat": where[0],
                        "lon": where[1],
                        "start": sampled + timedelta(hours=lo),
                        "end": sampled + timedelta(hours=hi),
                    }
                )
    return jobs


def summarise(todo, new_pairs, new_olci, products, sst) -> None:
    gained = {row["site_id"] for row in new_pairs}
    olci_gb = sum(float(row["size_bytes"] or 0) for row in new_olci.values()) / 1e9
    new_ids = set(new_olci)
    new_sst = {row["sst_product_id"]: float(row["sst_size_bytes"] or 0) for row in products if row["product_id"] in new_ids and row["sst_product_id"]}
    print(f"\nSites searched: {len(todo):,}; sites that gained candidate scenes: {len(gained):,}")
    print(f"New site x scene pairs: {len(new_pairs):,}")
    print(
        f"New OLCI scenes: {len(new_olci):,} ({olci_gb:,.0f} GB) with up to {len(new_sst):,} SST files "
        f"({sum(new_sst.values()) / 1e9:,.0f} GB; fewer if some are already on disk)"
    )
    print("Updated site_products.csv, products.csv and advected_positions.csv")
    print("Next: download_products.py (new scenes only), then download_tiles.py")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--trajectories-csv", type=Path, default=REPO_ROOT / "data" / "advection" / "trajectories.csv")
    parser.add_argument("--tile-log", type=Path, default=REPO_ROOT / "data" / "tile_log.csv")
    parser.add_argument("--max-days", type=float, default=6.0, help="Search out to this many days from sampling. Default 6.")
    parser.add_argument("--all-sites", action="store_true", help="Search every tracked site, not only those without an advected tile.")
    parser.add_argument("--workers", type=int, default=4, help="Parallel catalogue requests.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
