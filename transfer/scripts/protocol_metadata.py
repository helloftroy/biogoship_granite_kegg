#!/usr/bin/env python3
"""Sampling, extraction and sequencing protocol per modelled sample (RESULTS.md R34).

Every Bio-GO-SHIP library in ko_counts_cov50_all.tsv is linked to its FAIRe samp_name
(scripts/build_ko_targets.library_samples, the same link the KO targets use), then to the
FAIRe workbook:
  sampleMetadata         volume filtered (samp_size), size fraction, filter, storage solution,
                         collection method, extraction kit / protocol, extraction date,
                         extraction plate (extract_id 'PlateN_well'), DNA yield
  experimentRunMetadata  (assay MetaG_LibPrep_Illumina) sequencing run (seq_run_id: date,
                         flowcell, read layout), instrument, input read count
Values are per library; a sample_key with several libraries gets the first library's
categorical values and the sum of yields / read counts. Only the libraries that went
into data/targets/ko_samples.csv are kept.

Outputs (transfer/results/protocol/):
  sample_protocol.csv   one row per sample_key
  cruise_coverage.csv   per cruise: samples, share with each field, the distinct values

Usage:
    ~/miniforge3/bin/python transfer/scripts/protocol_metadata.py
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from biogoship.tables import MANIFEST_DIR  # noqa: E402
from build_ko_targets import library_samples  # noqa: E402

warnings.filterwarnings("ignore")
NA = ("not applicable", "not provided", "not collected", "missing")
SAMPLE_COLS = ["samp_size", "samp_size_unit", "size_frac", "size_frac_low", "filter_name", "filter_material",
               "samp_store_sol", "samp_store_temp", "samp_mat_process", "collection_method", "samp_collect_device",
               "samp_collect_method", "nucl_acid_ext_kit", "nucl_acid_ext", "nucl_acid_ext_lysis", "nucl_acid_ext_sep",
               "date_ext", "extract_id", "samp_vol_we_DNA_ext", "samp_vol_we_DNA_ext_unit", "dna_yield", "dna_yield_unit",
               "replicate_number"]  # fmt: skip


def clean(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip()
    return s.where(~s.str.lower().str.startswith(NA) & ~s.str.lower().isin(["nan", ""]))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", type=Path, default=REPO / "ko_counts_cov50_all.tsv")
    p.add_argument("--workbook", type=Path, default=REPO / "FAIRe-NOAA_uci-ess-biogoship.xlsx")
    p.add_argument("--faire-samples", type=Path, default=REPO / "data" / "faire" / "faire_samples.csv")
    p.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    p.add_argument("--ko-samples", type=Path, default=REPO / "data" / "targets" / "ko_samples.csv")
    p.add_argument("--out-dir", type=Path, default=REPO / "transfer" / "results" / "protocol")
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    libs = subprocess.run(["cut", "-f1", str(args.counts)], capture_output=True, text=True, check=True).stdout.split("\n")[1:]
    libs = pd.Index([x for x in libs if x])
    link = library_samples(libs, args)
    ks = pd.read_csv(args.ko_samples).set_index("sample_key")
    link = link[link["sample_key"].isin(ks.index) & (link["depth_m"].fillna(0) <= 10)]
    print(f"{len(libs):,} libraries; {len(link):,} link to the {len(ks):,} modelled samples")

    sm = pd.read_excel(args.workbook, sheet_name="sampleMetadata", header=2)
    sm = sm[~sm["samp_name"].astype(str).str.startswith("#")].set_index("samp_name")
    er = pd.read_excel(args.workbook, sheet_name="experimentRunMetadata", header=2)
    er = er[(er["assay_name"] == "MetaG_LibPrep_Illumina")]
    run = er.groupby("samp_name").agg(seq_run_id=("seq_run_id", "first"), seq_runs=("seq_run_id", "nunique"),
                                      instrument=("instrument", "first"),
                                      input_read_count=("input_read_count", lambda v: pd.to_numeric(v, errors="coerce").sum(min_count=1)))  # fmt: skip
    lib = link[["sample_key", "samp_name", "cruise"]].join(sm.reindex(columns=SAMPLE_COLS), on="samp_name").join(run, on="samp_name")
    for c in SAMPLE_COLS + ["seq_run_id", "instrument"]:
        lib[c] = clean(lib[c])

    # derived, numeric
    vol = pd.to_numeric(lib["samp_size"], errors="coerce")
    unit = lib["samp_size_unit"].str.lower()
    lib["volume_filtered_l"] = np.where(unit == "ml", vol / 1000, np.where(unit == "l", vol, np.nan))
    lib["dna_yield_ng"] = pd.to_numeric(lib["dna_yield"], errors="coerce")
    lib["extract_plate"] = lib["extract_id"].str.extract(r"^(Plate\d+)", expand=False)
    lib["extract_well_row"] = lib["extract_id"].str.extract(r"_(\d+)$", expand=False)
    lib["seq_run_date"] = lib["seq_run_id"].str.extract(r"^(\d{8})", expand=False)
    lib["seq_flowcell"] = lib["seq_run_id"].str.extract(r"^\d{8}_([A-Z0-9]+)_", expand=False)
    lib["seq_read_layout"] = lib["seq_run_id"].str.extract(r"_(PE\d+|SE\d+)$", expand=False)
    lib["instrument"] = lib["instrument"].str.extract(r"\[(OBI_\d+)\]", expand=False)
    lib["kit"] = lib["nucl_acid_ext_kit"]
    lib["lysis_buffer"] = lib["samp_store_sol"]

    first = lib.groupby("sample_key").first()
    out = first.drop(columns=["dna_yield_ng", "input_read_count", "volume_filtered_l"])
    out["dna_yield_ng"] = lib.groupby("sample_key")["dna_yield_ng"].sum(min_count=1)
    out["input_read_count"] = lib.groupby("sample_key")["input_read_count"].sum(min_count=1)
    out["volume_filtered_l"] = lib.groupby("sample_key")["volume_filtered_l"].mean()
    out["libraries"] = lib.groupby("sample_key").size()
    out["cruise"] = ks["cruise"].reindex(out.index)
    out["total_counts"] = ks["total"].reindex(out.index)
    out.to_csv(args.out_dir / "sample_protocol.csv")

    fields = ["volume_filtered_l", "size_frac", "filter_name", "lysis_buffer", "collection_method", "kit", "nucl_acid_ext",
              "date_ext", "extract_plate", "dna_yield_ng", "seq_run_id", "seq_flowcell", "seq_read_layout", "instrument",
              "input_read_count"]  # fmt: skip
    rows = []
    for cruise, g in out.groupby("cruise"):
        r = {"cruise": cruise, "samples": len(g)}
        for f in fields:
            v = g[f].dropna()
            r[f"{f} share"] = len(v) / len(g)
            if pd.api.types.is_numeric_dtype(g[f]):
                r[f"{f} values"] = f"{v.median():.3g} [{v.min():.3g}-{v.max():.3g}]" if len(v) else ""
            else:
                u = v.value_counts()
                r[f"{f} values"] = "; ".join(f"{k} ({n})" for k, n in u.head(6).items()) + (f"; +{len(u) - 6} more" if len(u) > 6 else "")
        rows.append(r)
    cov = pd.DataFrame(rows)
    cov.to_csv(args.out_dir / "cruise_coverage.csv", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 120)
    for f in fields:
        print(f"\n== {f}")
        print(cov[["cruise", "samples", f"{f} share", f"{f} values"]].to_string(index=False))


if __name__ == "__main__":
    main()
