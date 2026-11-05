#!/usr/bin/env python3
"""Per-sample KO (or Pfam) profiles as prediction targets: surface metagenomes, CLR.

Input: the gene-family count tables (rows = sequencing libraries, columns =
"count.ko.Kxxxxx" or "count.<pfam>"; coverage-weighted, so not whole numbers) and
the FAIRe workbook, whose experimentRunMetadata sheet links library names, run
accessions and file names to FAIRe samples.

Steps:
  1. Library -> sample. Original-cruise libraries named like C132020_021v2_meta
     are matched by name to biogoship_samples.csv; SRR accessions through its
     metagenome_runs; everything else through the FAIRe sample. (The FAIRe
     original_sample_id link is by time and place, so it can pick a neighbouring
     replicate; names are exact and win.)
  2. Drop controls (blanks, mocks, references), other projects, libraries with
     fewer than --min-total counts, and non-surface samples (depth > --max-depth;
     samples without a recorded depth are kept, as for the tiles).
  3. Proportions per library; libraries of one sample are averaged.
  4. Features present in at least --min-prevalence of the kept samples; CLR with
     half a count as pseudocount.

Sample key: our sample_id for the original cruises (as in data/embeddings), the
FAIRe samp_name otherwise (as in data/faire/embeddings).

Outputs in data/targets/: <name>_clr.csv (sample_key + one column per feature)
and <name>_samples.csv (cruise, libraries, total counts, sequencing instrument).

Usage:
    python scripts/build_ko_targets.py --counts ko_counts_cov50_all.tsv --name ko
    python scripts/build_ko_targets.py --counts pfam_counts_cov50_all.tsv --name pfam
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import MANIFEST_DIR, REPO_ROOT  # noqa: E402

CONTROL = re.compile(r"BLANK|MOCK|NEG|Reference|^NP_|^IC_")


def main() -> None:
    args = parse_args()
    counts = pd.read_csv(args.counts, sep="\t")  # header has one name fewer than the rows: first column is the index
    counts.columns = [re.sub(r"^count\.(ko\.)?", "", c) for c in counts.columns]
    print(f"{args.counts.name}: {counts.shape[0]:,} libraries x {counts.shape[1]:,} features")

    link = library_samples(counts.index, args)
    info = link.join(pd.Series(counts.sum(axis=1), name="total"))
    drop = info["samp_name"].isna() | info["samp_name"].fillna("").str.contains(CONTROL) | info["cruise"].isna()
    print(f"Libraries without a sample, controls or other projects: {int(drop.sum()):,}")
    info = info[~drop]
    low = info["total"] < args.min_total
    print(f"Libraries below {args.min_total:,.0f} counts: {int(low.sum()):,}")
    info = info[~low]
    deep = info["depth_m"] > args.max_depth
    print(f"Libraries from below {args.max_depth:g} m: {int(deep.sum()):,}")
    info = info[~deep]

    proportions = counts.loc[info.index].div(info["total"], axis=0)
    per_sample = proportions.groupby(info["sample_key"]).mean()
    totals = info.groupby("sample_key")["total"].sum()
    prevalence = (per_sample > 0).mean()
    keep = prevalence[prevalence >= args.min_prevalence].index
    print(f"Samples: {len(per_sample):,} ({int((info.groupby('sample_key').size() > 1).sum()):,} with more than one library, averaged); "
          f"features in >= {args.min_prevalence:.0%} of them: {len(keep):,} of {per_sample.shape[1]:,}")  # fmt: skip
    pseudo = (0.5 / totals).to_numpy()[:, None]
    logs = np.log(per_sample[keep].to_numpy() + pseudo)
    clr = pd.DataFrame(logs - logs.mean(axis=1, keepdims=True), index=per_sample.index, columns=keep)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    clr.round(5).rename_axis("sample_key").to_csv(args.out_dir / f"{args.name}_clr.csv")
    samples = info.groupby("sample_key").agg(cruise=("cruise", "first"), libraries=("samp_name", "size"), total=("total", "sum"),
                                             instrument=("instrument", "first"), depth_m=("depth_m", "first"))  # fmt: skip
    samples.to_csv(args.out_dir / f"{args.name}_samples.csv")
    print(samples.groupby("cruise").agg(samples=("libraries", "size"), instrument=("instrument", "first")).to_string())
    print(f"Wrote {args.out_dir / f'{args.name}_clr.csv'} and {args.name}_samples.csv")


def library_samples(libraries: pd.Index, args: argparse.Namespace) -> pd.DataFrame:
    """Per library: FAIRe samp_name, sample_key, cruise, depth and sequencing instrument."""
    runs = pd.read_excel(args.workbook, sheet_name="experimentRunMetadata", header=2)
    runs = runs[~runs["samp_name"].astype(str).str.startswith("#")]
    token_to_sample: dict[str, str] = {}
    for _, row in runs.iterrows():
        for column, value in row.items():
            if column == "samp_name" or pd.isna(value):
                continue
            for token in re.split(r"[\s|/,;]+", str(value)):
                token_to_sample.setdefault(re.sub(r"(_R[12])?(_001)?\.fastq\.gz$", "", token), row["samp_name"])
    out = pd.DataFrame(index=libraries)
    out["samp_name"] = [token_to_sample.get(lib) or token_to_sample.get(re.sub(r"_meta.*$", "", lib)) for lib in libraries]

    ours = pd.read_csv(args.manifest_dir / "biogoship_samples.csv")
    by_run = ours.assign(run=ours["metagenome_runs"].astype(str).str.split(r"[;, ]+")).explode("run").set_index("run")["sample_id"]
    by_name = pd.Series(libraries.str.replace(r"_meta.*$", "", regex=True), index=libraries)
    faire = pd.read_csv(args.faire_samples).set_index("samp_name")
    original = (by_name.where(by_name.isin(ours["sample_id"]))
                .fillna(pd.Series(libraries.map(by_run.to_dict()), index=libraries))
                .fillna(out["samp_name"].map(faire["original_sample_id"])))  # fmt: skip
    out["sample_key"] = original.fillna(out["samp_name"])
    out["cruise"] = out["samp_name"].map(faire["cruise"])
    depth = ours.set_index("sample_id")["depth_m"]
    out["depth_m"] = original.map(depth).fillna(out["samp_name"].map(faire["depth_min_m"]))
    instrument = runs.drop_duplicates("samp_name").set_index("samp_name")["instrument"].astype(str).str.extract(r"(OBI_\d+)")[0]
    out["instrument"] = out["samp_name"].map(instrument)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", type=Path, default=REPO_ROOT / "ko_counts_cov50_all.tsv")
    parser.add_argument("--name", default="ko")
    parser.add_argument("--workbook", type=Path, default=REPO_ROOT / "FAIRe-NOAA_uci-ess-biogoship.xlsx")
    parser.add_argument("--faire-samples", type=Path, default=REPO_ROOT / "data" / "faire" / "faire_samples.csv")
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "targets")
    parser.add_argument("--max-depth", type=float, default=10.0)
    parser.add_argument("--min-total", type=float, default=1e5, help="Smallest library total kept.")
    parser.add_argument("--min-prevalence", type=float, default=0.9, help="Share of samples a feature must be present in.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
