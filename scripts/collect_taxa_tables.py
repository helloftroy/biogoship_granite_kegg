#!/usr/bin/env python3
"""Collect the KO pipeline's per-library taxonomy outputs and build summary tables.

The pipeline writes one folder per library (A392013_008_meta, SRR12480057, ...) with,
per taxonomic level (1 domain ... 7 species):
  level_<n>_<rank>_abundance.tsv   taxon, classification_method, read_count, relative_reads_mapped
  level_<n>_<rank>_ko_counts.tsv   taxon, ko, length, count (gene counts per taxon and KO)
  ..._og_counts.tsv, ..._pfam_counts.tsv
They arrive as Google Drive zips (one library's files can be split over several parts)
or as folders (rclone).

Step 1, stash (incremental; a file already stashed is not replaced):
  data/targets/taxa/raw/<library>/level_<n>_<rank>_abundance.tsv.gz        every level
  data/targets/taxa/ko_by_taxon/<library>/level_<n>_<rank>_ko_counts.tsv.gz  --ko-levels only
  (default 2 = phylum and 4 = order). OG and Pfam tables and other KO levels are not kept.
Step 2, build from everything stashed so far, in data/targets/taxa/:
  taxa_long.csv.gz      every abundance row: library, level, taxon, method, read_count, relative_reads_mapped
  <rank>_counts.tsv     wide read counts, libraries x taxa, unclassified rows dropped, in the layout of
                        ko_counts_cov50_all.tsv, for build_ko_targets.py:
                          python scripts/build_ko_targets.py --counts data/targets/taxa/order_counts.tsv \
                              --name taxa_order --out-dir data/targets/taxa
  libraries.csv         per library: abundance levels, KO-by-taxon levels, total reads, share unclassified
A taxon listed twice in one library (GTDB-Tk for MAGs, EUKulele otherwise) has its reads summed.

Usage:
    python scripts/collect_taxa_tables.py drive-download-*.zip      # stash, then build
    python scripts/collect_taxa_tables.py                            # rebuild from the stash only
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

PATTERN = re.compile(r"(?:^|/)([^/]+)/(level_(\d)_[a-z_]+?_(abundance|ko_counts|og_counts|pfam_counts))\.tsv(\.gz)?$")
RANKS = {1: "domain", 2: "phylum", 3: "class", 4: "order", 5: "family", 6: "genus", 7: "species"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sources", type=Path, nargs="*", help="Zip files or folders holding <library>/level_*.tsv")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "targets" / "taxa")
    parser.add_argument("--ko-levels", type=int, nargs="*", default=[2, 4], help="KO-by-taxon levels to keep")
    args = parser.parse_args()

    raw, kbt = args.out_dir / "raw", args.out_dir / "ko_by_taxon"
    for source in args.sources:
        added = skipped = 0
        for name, read in files(source):
            m = PATTERN.search(name)
            if not m:
                continue
            library, stem, level, kind, gz = m.group(1), m.group(2), int(m.group(3)), m.group(4), m.group(5)
            if kind == "abundance":
                target = raw / library / f"{stem}.tsv.gz"
            elif kind == "ko_counts" and level in args.ko_levels:
                target = kbt / library / f"{stem}.tsv.gz"
            else:
                continue
            if target.exists():
                skipped += 1
                continue
            data = read()
            target.parent.mkdir(parents=True, exist_ok=True)
            tmp = target.with_suffix(".tmp")
            tmp.write_bytes(data if gz else gzip.compress(data, 6))
            tmp.rename(target)
            added += 1
        print(f"{source.name}: {added} files stashed, {skipped} already there", flush=True)

    build(args.out_dir, raw, kbt)


def build(out: Path, raw: Path, kbt: Path) -> None:
    parts = []
    for path in sorted(raw.glob("*/level_*_abundance.tsv.gz")):
        level = int(path.name.split("_")[1])
        parts.append(pd.read_csv(path, sep="\t").assign(library=path.parent.name, level=level))
    if not parts:
        print("Nothing stashed yet.")
        return
    long = pd.concat(parts, ignore_index=True)
    long = long[["library", "level", "taxon", "classification_method", "read_count", "relative_reads_mapped"]]
    long.to_csv(out / "taxa_long.csv.gz", index=False)
    unclassified = long["taxon"].str.split(";").str[-1].str.contains(r"^unclassified|^Unclassified", regex=True)
    print(f"Abundance tables: {long['library'].nunique():,} libraries")
    for level, rank in RANKS.items():
        t = long[(long["level"] == level) & ~unclassified]
        wide = t.pivot_table(index="library", columns="taxon", values="read_count", aggfunc="sum", fill_value=0)
        wide.sort_index().to_csv(out / f"{rank}_counts.tsv", sep="\t", index_label=False)
        print(f"  {rank:8s} {wide.shape[1]:6,} taxa")

    domain = long[long["level"] == 1]
    total = domain.groupby("library")["read_count"].sum()
    ko_levels = pd.Series({d.name: ",".join(sorted(p.name.split("_")[1] for p in d.glob("*.tsv.gz")))
                           for d in kbt.glob("*") if d.is_dir()}, dtype=object)  # fmt: skip
    info = pd.DataFrame({
        "abundance_levels": long.groupby("library")["level"].nunique(),
        "ko_by_taxon_levels": ko_levels,
        "total_reads": total,
        "share_unclassified": domain[domain["taxon"].str.startswith("Unclassified")].groupby("library")["read_count"].sum() / total,
    })  # fmt: skip
    info["abundance_levels"] = info["abundance_levels"].fillna(0).astype(int)
    info.round(4).rename_axis("library").sort_index().to_csv(out / "libraries.csv")
    print(f"KO-by-taxon: {len(ko_levels):,} libraries; written to {out}")


def files(source: Path):
    """(name, reader) for every file in a zip or under a folder."""
    if source.suffix == ".zip":
        z = zipfile.ZipFile(source)
        for name in z.namelist():
            yield name, (lambda n=name: z.read(n))
    else:
        for path in source.rglob("level_*.tsv*"):
            yield str(path.relative_to(source.parent)), path.read_bytes


if __name__ == "__main__":
    main()
