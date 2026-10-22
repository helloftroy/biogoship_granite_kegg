#!/usr/bin/env python3
"""Find the Sentinel-3 products that cover each Bio-GO-SHIP site.

For every site this asks the Copernicus catalogue for OLCI level-2 water (OL_2_WFR)
scenes whose footprint contains the site and whose sensing time is within
--max-days of sampling. Each OLCI scene is then paired with the SLSTR sea surface
temperature (SL_2_WST) orbit file of the same overpass.

Re-running this rewrites both output files from scratch, which drops any pairs
added later by search_trajectory_products.py; re-run that afterwards if needed.

No credentials are needed: the catalogue is public. Outputs in data/manifests/:
  site_products.csv   one row per (site, candidate OLCI scene), with the time gap
  products.csv        one row per unique OLCI scene, with its paired SST product

Usage:
    python scripts/search_products.py
    python scripts/search_products.py --max-days 3 --section-max-days I09N=10
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import cdse  # noqa: E402
from biogoship.tables import MANIFEST_DIR, parse_utc, read_csv, write_csv  # noqa: E402

# No Sentinel-3 data exists before the Sentinel-3A launch, so earlier sites are not queried.
SENTINEL3_LAUNCH = datetime(2016, 2, 16, tzinfo=timezone.utc)
# The OLCI and SLSTR products of one overpass overlap in time; allow a little slack.
MAX_SST_GAP_MINUTES = 180.0

SITE_PRODUCT_FIELDS = ["site_id", "product_id", "product_name", "satellite", "start_utc", "time_gap_hours"]
PRODUCT_FIELDS = [
    "product_id",
    "product_name",
    "satellite",
    "start_utc",
    "end_utc",
    "size_bytes",
    "n_sites",
    "min_abs_gap_hours",
    "sst_product_id",
    "sst_product_name",
    "sst_size_bytes",
    "sst_gap_minutes",
    "sst_relation",
]


def main() -> None:
    args = parse_args()
    sites = read_csv(args.sites_csv)
    section_days = dict(args.section_max_days)

    def window(site: dict[str, str]) -> timedelta:
        return timedelta(days=section_days.get(site["section"], args.max_days))

    queryable = [site for site in sites if parse_utc(site["time_utc"]) + window(site) >= SENTINEL3_LAUNCH]
    print(f"Sites: {len(sites):,} total, {len(queryable):,} sampled close enough to the Sentinel-3 era to query")
    print(f"Window: +/-{args.max_days:g} days" + "".join(f"; {name} +/-{days:g}" for name, days in section_days.items()))

    def search_site(site: dict[str, str]) -> list[dict[str, str]]:
        when = parse_utc(site["time_utc"])
        time.sleep(args.pause)  # stay under the catalogue's rate limit on long runs
        return cdse.search_products(
            cdse.OLCI_WFR, when - window(site), when + window(site), point=(float(site["lon"]), float(site["lat"]))
        )

    site_products: list[dict] = []
    olci: dict[str, dict] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for done, (site, found) in enumerate(zip(queryable, pool.map(search_site, queryable)), start=1):
            when = parse_utc(site["time_utc"])
            for product in found:
                gap_hours = (parse_utc(product["start_utc"]) - when).total_seconds() / 3600
                site_products.append({**product, "site_id": site["site_id"], "time_gap_hours": f"{gap_hours:.2f}"})
                olci.setdefault(product["product_id"], product)
            if done % 50 == 0 or done == len(queryable):
                print(f"  searched {done:,}/{len(queryable):,} sites; {len(olci):,} unique OLCI scenes", flush=True)

    sst = search_sst(list(olci.values()))
    products = []
    for product in sorted(olci.values(), key=lambda row: row["start_utc"]):
        gaps = [abs(float(r["time_gap_hours"])) for r in site_products if r["product_id"] == product["product_id"]]
        products.append({**product, "n_sites": len(gaps), "min_abs_gap_hours": f"{min(gaps):.2f}", **pair_sst(product, sst)})

    site_products.sort(key=lambda row: (row["site_id"], abs(float(row["time_gap_hours"]))))
    write_csv(args.out_dir / "site_products.csv", site_products, SITE_PRODUCT_FIELDS)
    write_csv(args.out_dir / "products.csv", products, PRODUCT_FIELDS)
    summarise(sites, site_products, products, args)


def search_sst(olci: list[dict[str, str]]) -> list[dict[str, str]]:
    """All SST orbit files sensed around the OLCI scenes, fetched in a few date-range queries."""
    day = timedelta(days=1)
    days = sorted({parse_utc(row["start_utc"]).replace(hour=0, minute=0, second=0, microsecond=0) for row in olci})
    ranges: list[list[datetime]] = []
    for start in days:  # each scene needs SST from the day before through the day after
        if ranges and start - day <= ranges[-1][1]:
            ranges[-1][1] = start + 2 * day
        else:
            ranges.append([start - day, start + 2 * day])
    found: dict[str, dict[str, str]] = {}
    for start, end in ranges:
        for product in cdse.search_products(cdse.SLSTR_WST, start, end):
            found[product["product_id"]] = product
    print(f"SST orbit files in {len(ranges)} date range(s): {len(found):,}", flush=True)
    return list(found.values())


def pair_sst(olci: dict[str, str], sst: list[dict[str, str]]) -> dict[str, str]:
    """Closest SST orbit file in time: the same satellite if possible, else the other one."""
    olci_time = parse_utc(olci["start_utc"])

    def gap_minutes(row: dict[str, str]) -> tuple[float, float]:
        start, end = parse_utc(row["start_utc"]), parse_utc(row["end_utc"])
        to_midpoint = abs((olci_time - (start + (end - start) / 2)).total_seconds()) / 60
        if start <= olci_time <= end:
            return 0.0, to_midpoint
        return min(abs((olci_time - start).total_seconds()), abs((olci_time - end).total_seconds())) / 60, to_midpoint

    for relation, same in (("same_satellite", True), ("cross_satellite", False)):
        candidates = [row for row in sst if (row["satellite"] == olci["satellite"]) == same]
        if not candidates:
            continue
        best = min(candidates, key=gap_minutes)
        gap = gap_minutes(best)[0]
        if gap <= MAX_SST_GAP_MINUTES:
            return {
                "sst_product_id": best["product_id"],
                "sst_product_name": best["product_name"],
                "sst_size_bytes": best["size_bytes"],
                "sst_gap_minutes": f"{gap:.1f}",
                "sst_relation": relation,
            }
    return {"sst_product_id": "", "sst_product_name": "", "sst_size_bytes": "", "sst_gap_minutes": "", "sst_relation": ""}


def summarise(sites: list[dict], site_products: list[dict], products: list[dict], args: argparse.Namespace) -> None:
    covered = {row["site_id"] for row in site_products}
    print(f"\nSites with at least one candidate OLCI scene: {len(covered):,} of {len(sites):,}")
    for section in sorted({site["section"] for site in sites}):
        in_section = [site for site in sites if site["section"] == section]
        hit = sum(1 for site in in_section if site["site_id"] in covered)
        print(f"  {section:12s} {hit:4d} / {len(in_section):4d}")
    paired = sum(1 for row in products if row["sst_product_id"])
    olci_gb = sum(float(row["size_bytes"] or 0) for row in products) / 1e9
    print(f"Unique OLCI scenes: {len(products):,} ({olci_gb:,.0f} GB if downloaded whole); {paired:,} paired with SST")
    print(f"Wrote {args.out_dir / 'site_products.csv'}")
    print(f"Wrote {args.out_dir / 'products.csv'}")


def section_days_arg(text: str) -> tuple[str, float]:
    section, _, days = text.partition("=")
    return section, float(days)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites-csv", type=Path, default=MANIFEST_DIR / "biogoship_sites.csv")
    parser.add_argument("--out-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument(
        "--max-days", type=float, default=3.0, help="Largest gap between sampling and overpass. Default 3 days."
    )
    parser.add_argument(
        "--section-max-days",
        action="append",
        default=[],
        type=section_days_arg,
        metavar="SECTION=DAYS",
        help="Wider (or narrower) window for one section, e.g. I09N=10. May be repeated.",
    )
    parser.add_argument("--workers", type=int, default=4, help="Parallel catalogue requests.")
    parser.add_argument("--pause", type=float, default=0.25, help="Seconds each worker waits before a request.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
