#!/usr/bin/env python3
"""R70: what "depth" means. Raw sequencing effort vs KO-assigned reads, and their link to the eukaryote fraction.

The technical covariate everywhere so far (main model, hurdle models, R33/R34) is log10 KO-assigned reads (per sample
the summed KO counts of its libraries, ko_samples.csv "total"; per tile the mean of log10 over its samples). KO
assignment depends on the community: eukaryotic reads map and annotate poorly, so eukaryote-rich samples give fewer KO
reads from the same sequencing. This script builds the raw effort, read pairs sequenced (FAIRe input_read_count, from
transfer/results/library_metadata/library_metadata.csv), with the same library -> sample mapping and filters as the
targets (build_ko_targets.py), and splits log10 KO reads = log10 read pairs + log10 KO yield.

Outputs (data/results/sequencing_effort/): sample_effort.csv, tile_effort.csv, summary.txt
Usage: ~/miniforge3/bin/python3 scripts/sequencing_effort.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
OUT = REPO / "data" / "results" / "sequencing_effort"
BX = REPO / "data" / "results" / "benchmark_export"
LOG: list = []


def say(s=""):
    print(s, flush=True)
    LOG.append(str(s))


def build():
    from biogoship.tables import MANIFEST_DIR, REPO_ROOT
    from build_ko_targets import CONTROL, library_samples

    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv", index_col="sample_key")
    tot = pd.read_csv(REPO / "ko_counts_cov50_all.tsv", sep="\t", index_col=0).sum(axis=1).rename("total")
    ns = argparse.Namespace(workbook=REPO_ROOT / "FAIRe-NOAA_uci-ess-biogoship.xlsx",
                            faire_samples=REPO / "data" / "faire" / "faire_samples.csv", manifest_dir=MANIFEST_DIR)  # fmt: skip
    info = library_samples(tot.index, ns).join(tot)
    keep = (info["samp_name"].notna() & ~info["samp_name"].fillna("").str.contains(CONTROL) & info["cruise"].notna()
            & info["sample_key"].isin(ks.index) & (info["total"] >= 1e5) & ~(info["depth_m"] > 10))  # fmt: skip
    info = info[keep]
    meta = pd.read_csv(REPO / "transfer/results/library_metadata/library_metadata.csv").drop_duplicates("lib_id").set_index("lib_id")
    info["read_pairs"] = meta["total_reads"].reindex(info.index)
    say(f"libraries {len(info):,} for {info['sample_key'].nunique():,} of {len(ks):,} target samples; "
        f"with read pairs {info['read_pairs'].notna().sum():,}")  # fmt: skip
    s = info.groupby("sample_key").agg(ko_reads=("total", "sum"), read_pairs=("read_pairs", lambda x: x.sum(min_count=len(x))),
                                       libraries=("total", "size"))  # fmt: skip
    chk = np.log10(s["ko_reads"]) - np.log10(ks["total"].reindex(s.index))
    say(f"check: per-sample KO reads vs ko_samples.csv total, max |log10 diff| {chk.abs().max():.2e} "
        f"(samples missing in ko_samples: {chk.isna().sum()})")  # fmt: skip
    s["log10_read_pairs"] = np.log10(s["read_pairs"])
    s["log10_ko_reads"] = np.log10(s["ko_reads"])
    s["log10_ko_yield"] = s["log10_ko_reads"] - s["log10_read_pairs"]
    say(f"samples without read pairs for every library: {s['read_pairs'].isna().sum():,}")
    OUT.mkdir(parents=True, exist_ok=True)
    s.round(5).to_csv(OUT / "sample_effort.csv")

    tiles = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    mem = tiles["sample_keys"].str.split(";").explode()
    cols = ["log10_read_pairs", "log10_ko_reads", "log10_ko_yield"]
    t = s.reindex(mem.values)[cols].set_axis(mem.index).groupby(level=0).mean().reindex(tiles.index)  # mean of log10, as main_model
    t["samples_missing_pairs"] = s.reindex(mem.values)["read_pairs"].isna().set_axis(mem.index).groupby(level=0).sum()
    t = t.join(tiles[["cruise", "platform_group", "log10_ko_reads"]].rename(columns={"log10_ko_reads": "export_log10_ko_reads"}))
    C = pd.read_csv(BX / "targets" / "C.csv.gz", index_col="tile_key")["euk_prok_log_ratio"]
    t["C"] = C.reindex(t.index)
    say(f"tiles {len(t):,}; with read pairs for all samples {(t['samples_missing_pairs'] == 0).sum():,}; "
        f"r(log10 KO reads rebuilt, export) {t[['log10_ko_reads', 'export_log10_ko_reads']].corr().iloc[0, 1]:.4f}")  # fmt: skip
    t.round(5).to_csv(OUT / "tile_effort.csv")
    return t


def summarise(t):
    t = t[t["samples_missing_pairs"] == 0]
    v = ["log10_read_pairs", "log10_ko_reads", "log10_ko_yield", "C"]
    w = t[v] - t.groupby("cruise")[v].transform("mean")
    say(f"\ntiles used {len(t):,}")
    say("within-cruise Pearson r (tile values minus cruise mean):")
    say(w.corr().round(3).to_string())
    say("\npooled Pearson r:")
    say(t[v].corr().round(3).to_string())
    say(f"\nwithin-cruise SD: " + ", ".join(f"{c} {w[c].std():.3f}" for c in v))
    say(f"share of within-cruise variance of log10 KO reads shared with log10 read pairs (r^2): "
        f"{w.corr().loc['log10_ko_reads', 'log10_read_pairs'] ** 2:.3f}")  # fmt: skip
    pc = t.groupby("cruise").apply(lambda g: pd.Series({"tiles": len(g), "r_pairs_C": g["log10_read_pairs"].corr(g["C"]),
                                                       "r_koreads_C": g["log10_ko_reads"].corr(g["C"]),
                                                       "r_yield_C": g["log10_ko_yield"].corr(g["C"])}))  # fmt: skip
    say("\nper cruise r with C:")
    say(pc.round(2).to_string())
    say(f"cruises with r(yield, C) < 0: {(pc['r_yield_C'] < 0).sum()} of {len(pc)}")
    (OUT / "summary.txt").write_text("\n".join(LOG) + "\n")


if __name__ == "__main__":
    summarise(build())
