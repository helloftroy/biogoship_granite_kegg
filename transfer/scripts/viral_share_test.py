#!/usr/bin/env python3
"""Viral read share against the residual axes (RESULTS.md R42). Descriptive, not causal: viral reads
come from the same libraries (a lysis-pressure proxy at best; most viral reads in a 0.22 um
Sterivex metagenome are prophages, large viruses on cells, or virions retained on the filter).

Per sample (data/targets/taxa/domain_counts.tsv, libraries summed per sample via ko_samples):
  Viruses CLR (of the 4 domains), log10 Viruses/Bacteria, log10 Viruses/Cyanobacteriota reads
  (phylum counts), unclassified-free viral share.
Per tile: mean over its samples (real model tiles, transfer/results/residual_modules_noA/scores.csv).
Spearman all / within cruise (R33 one_driver) for F1, F2, W1, W2, W3; and the R40 platform eta2 of
the per-cruise means (does the viral share itself carry the platform batch?).

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/viral_share_test.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from residual_modules import one_driver  # noqa: E402
from tune_figure_kos import load  # noqa: E402
import argparse  # noqa: E402

OUT = REPO / "transfer" / "results" / "viral_share"
OUT.mkdir(parents=True, exist_ok=True)
args = argparse.Namespace(targets=REPO / "data" / "targets" / "ko_clr.csv",
                          embeddings_dirs=[REPO / "data" / "embeddings_sst_l4", REPO / "data" / "faire" / "embeddings_sst_l4"])
kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()[:3]
_, _, _, _, cruises, unit, rows = load(args, kos, return_rows=True)
sc = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "scores.csv", index_col="tile_key").loc[unit.index]

lib = pd.read_csv(REPO / "transfer" / "results" / "library_metadata" / "library_metadata.csv")
lib = lib[lib["in_real_targets"]].drop_duplicates("lib_id").set_index("lib_id")
dom = pd.read_csv(REPO / "data" / "targets" / "taxa" / "domain_counts.tsv", sep="\t", index_col=0)
phy = pd.read_csv(REPO / "data" / "targets" / "taxa" / "phylum_counts.tsv", sep="\t", index_col=0)
cyan = [c for c in phy.columns if c.endswith("Cyanobacteriota")]
d = dom.reindex(lib.index).groupby(lib["sample_id"]).sum(min_count=1)
cy = phy[cyan].sum(axis=1).reindex(lib.index).groupby(lib["sample_id"]).sum(min_count=1)
lg = np.log(d + 0.5)
feat = pd.DataFrame({
    "Viruses CLR (4 domains)": lg["Viruses"] - lg.mean(axis=1),
    "log10 Viruses / Bacteria": np.log10((d["Viruses"] + 0.5) / (d["Bacteria"] + 0.5)),
    "log10 Viruses / Cyanobacteriota": np.log10((d["Viruses"] + 0.5) / (cy + 0.5)),
    "viral share of classified reads": d["Viruses"] / d.sum(axis=1),
})
print(f"samples with domain counts: {feat.dropna().shape[0]:,}")
tile = {c: rows["sample_id"].map(v).astype(float).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy() for c, v in feat.items()}
out = []
log = []
for c, x in tile.items():
    parts = []
    for comp in ["F1", "F2", "W1", "W2", "W3"]:
        r = one_driver(c, sc[comp].to_numpy(), x, cruises)
        if r:
            out.append({**r, "component": comp})
            parts.append(f"{comp} {r['all_r']:+.2f}/{r['within_r']:+.2f} ({r['cruises_same_sign']}/{r['cruises_tested']})")
    line = f"{c:34s} n {int(np.isfinite(x).sum())}  " + "  ".join(parts)
    print(line); log.append(line)
pd.DataFrame(out).round(4).to_csv(OUT / "drivers.csv", index=False)
# platform batch on the viral share itself
plat = {"I09N": "6000", "P18": "6000", "I07N": "6000", "AMT28": "6000", "C13.5/A13.5": "6000", "HLY2202-2022": "6000"}
cm = pd.DataFrame(tile, index=unit.index).groupby(cruises).mean()
cm["platform"] = [plat.get(c, "X") for c in cm.index]
rng = np.random.default_rng(0)
for c in tile:
    x = cm[c]; g = cm["platform"]
    def eta(gg):
        m = x.groupby(gg.values).transform("mean"); return ((m - x.mean()) ** 2).sum() / ((x - x.mean()) ** 2).sum()
    e = eta(g); null = [eta(pd.Series(rng.permutation(g.values), index=g.index)) for _ in range(9999)]
    line = f"platform eta2 of cruise means, {c}: {e:.2f} (p {(1 + (np.array(null) >= e).sum()) / 10000:.3f}); NovaSeq 6000 mean {x[g == '6000'].mean():.3f}, X mean {x[g == 'X'].mean():.3f}"
    print(line); log.append(line)
cm.round(4).to_csv(OUT / "cruise_means.csv")
(OUT / "summary.txt").write_text("\n".join(log) + "\n")
