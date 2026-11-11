#!/usr/bin/env python3
"""Collect the per-library taxonomy tables of the KO pipeline into one table per level.

Input: pipeline output folders named after the library (A392013_008_meta,
SRR12480057, ...), each with level_<n>_<rank>_abundance.tsv (taxon,
classification_method, read_count, relative_reads_mapped). The folders can be
given as directories or as the Google Drive zip downloads that contain them.

A taxon can appear twice in one library, once per classification method
(GTDB-Tk for prokaryotic MAGs, EUKulele otherwise); the read counts are summed.

Outputs in --out-dir (default data/targets/taxa/):
  taxa_long.csv.gz                    every row: library, level, taxon, method, read_count, relative_reads_mapped
  <rank>_counts.tsv                   wide read counts, libraries x taxa, in the layout of
                                      ko_counts_cov50_all.tsv (header one name short), without
                                      unclassified rows, ready for build_ko_targets.py:
                                        python scripts/build_ko_targets.py --counts data/targets/taxa/order_counts.tsv --name taxa_order
  libraries.csv                       per library: levels found, total reads, share unclassified (domain level)

Usage:
    python scripts/collect_taxa_tables.py drive-download-*.zip A392013_008_meta-*.zip
"""

from __future__ import annotations

import argparse
import io
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

PATTERN = re.compile(r"(?:^|/)([^/]+)/level_(\d)_([a-z_]+)_abundance\.tsv$")
RANKS = {"1": "domain", "2": "phylum", "3": "class", "4": "order", "5": "family", "6": "genus", "7": "species"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sources", type=Path, nargs="+", help="Zip files or folders holding <library>/level_*_abundance.tsv")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "targets" / "taxa")
    args = parser.parse_args()

    parts, seen = [], {}
    for source in args.sources:
        for name, read in files(source):
            m = PATTERN.search(name)
            if not m:
                continue
            library, level = m.group(1), m.group(2)
            if (library, level) in seen:
                continue  # the same library in two downloads: keep the first
            seen[(library, level)] = source.name
            t = pd.read_csv(io.BytesIO(read()), sep="\t")
            parts.append(t.assign(library=library, level=int(level)))
    long = pd.concat(parts, ignore_index=True)
    long = long[["library", "level", "taxon", "classification_method", "read_count", "relative_reads_mapped"]]
    print(f"{long['library'].nunique():,} libraries, {len(long):,} rows from {len(args.sources)} sources")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    long.to_csv(args.out_dir / "taxa_long.csv.gz", index=False)
    unclassified = long["taxon"].str.split(";").str[-1].str.contains(r"^unclassified|^Unclassified", regex=True)
    for level, rank in RANKS.items():
        t = long[(long["level"] == int(level)) & ~unclassified]
        wide = t.pivot_table(index="library", columns="taxon", values="read_count", aggfunc="sum", fill_value=0)
        wide.sort_index().to_csv(args.out_dir / f"{rank}_counts.tsv", sep="\t", index_label=False)
        print(f"  {rank:8s} {wide.shape[1]:6,} taxa")

    domain = long[long["level"] == 1]
    info = pd.DataFrame({
        "levels": long.groupby("library")["level"].nunique(),
        "total_reads": domain.groupby("library")["read_count"].sum(),
        "share_unclassified": domain[domain["taxon"].str.startswith("Unclassified")].groupby("library")["read_count"].sum()
                              / domain.groupby("library")["read_count"].sum(),
    }).fillna({"share_unclassified": 0})  # fmt: skip
    info.round(4).rename_axis("library").to_csv(args.out_dir / "libraries.csv")
    print(f"Wrote {args.out_dir}")


def files(source: Path):
    """(name, reader) for every file in a zip or under a folder."""
    if source.suffix == ".zip":
        z = zipfile.ZipFile(source)
        for name in z.namelist():
            yield name, (lambda n=name: z.read(n))
    else:
        for path in source.rglob("*.tsv"):
            yield str(path.relative_to(source.parent)), path.read_bytes


if __name__ == "__main__":
    main()
