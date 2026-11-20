#!/usr/bin/env python3
"""Is the platform / axis-1 offset a sequencing GC bias? Per-KO GC vs per-KO effects (RESULTS.md R44).

Per-KO GC: transfer/results/ko_gc/ko_gc.csv (transfer/scripts/ko_gc_content.py on Morrill: 1-in-10 crc32
sample of the reference genes with a KO; gene-weighted mean GC, overall and per domain; eukaryote share
= genes_Eukaryota / genes).

Per-KO effects, over the real-model KO set (data/targets/ko_clr.csv; 3,933 KOs), in units of the KO's
pooled within-cruise SD (as R39):
  loadings     F1, F2, F3, W1, W2, W3 (transfer/results/residual_modules_noA)
  hiseq        HiSeq 4000 - NovaSeq 6000, different samples within a cruise (I09N-2016, AE1319-2013,
               BVAL46-2011, NH1418-2014; surface libraries); per-cruise mean difference, pooled weighted by
               the smaller group. Library CLR after rarefying to 250,000 KO reads (as sequencer_contrasts.py).
  era          NovaSeq 6000 cruises - NovaSeq X ("...LT4") cruises: mean of cruise-mean library CLR
               (modelled surface libraries; platform + years + regions together).
  rerun        HLY2202 same DNA on two NovaSeq 6000 runs (HMY7WDSX7 - HFJLFDSX7), mean paired difference
               (control: same chemistry, expect no GC pattern).
  protocolA    A16N Protocol A - B, mean paired difference (transfer/results/axis2_protocol/pairs_list.csv).
  c13_v2       C13.5 v2 - v1 (duplicates), mean paired difference.
Statistics: Spearman rho of each effect with gc_mean; partial Spearman (ranks residualised on the ranks of
eukaryote share and mean gene length); the same with the bacterial-gene GC (gc_mean_Bacteria) among KOs
with >= 20 sampled bacterial genes. KOs are not independent (shared pathways, shared taxa), so the p values
are for orientation; the comparison between contrasts (and with the rerun control) carries the argument.
Also: mean effect per GC decile (figure).

Outputs: transfer/results/ko_gc/gc_effects.csv (per KO), gc_correlations.csv, summary.txt;
transfer/figures/ko_gc_effects.png.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/ko_gc_test.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "transfer" / "results" / "ko_gc"
FIG = REPO / "transfer" / "figures"
LIBDIR = REPO / "transfer" / "results" / "library_metadata"
N = 250_000


def rarefied_clr(lib: pd.DataFrame, kos: list[str]) -> pd.DataFrame:
    cache = OUT / "rarefied_clr.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        return pd.DataFrame(z["x"], index=z["index"], columns=z["columns"])
    rng = np.random.default_rng(0)
    header = open(REPO / "ko_counts_cov50_all.tsv").readline().rstrip("\n").split("\t")
    col_idx = {c.replace("count.ko.", ""): i for i, c in enumerate(header)}
    sel = [col_idx[k] for k in kos]
    want, rows, keys = set(lib.index), [], []
    with open(REPO / "ko_counts_cov50_all.tsv") as fh:
        fh.readline()
        for line in fh:
            name, rest = line.split("\t", 1)
            if name not in want:
                continue
            v = np.fromstring(rest, sep="\t")
            if v.sum() < N:
                continue
            p = rng.multinomial(N, v / v.sum())[sel] / N
            lg = np.log(p + 0.5 / N)
            rows.append(lg - lg.mean())
            keys.append(name)
    clr = pd.DataFrame(rows, index=keys, columns=kos)
    np.savez_compressed(cache, x=clr.to_numpy(np.float32), index=np.array(keys, dtype=object), columns=np.array(kos, dtype=object))
    return clr


def partial_spearman(x, y, covs):
    rx, ry = rankdata(x), rankdata(y)
    c = np.column_stack([np.ones(len(x))] + [rankdata(v) for v in covs])
    ex = rx - c @ np.linalg.lstsq(c, rx, rcond=None)[0]
    ey = ry - c @ np.linalg.lstsq(c, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def main() -> None:
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    OUT.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", index_col="sample_key")
    kos = list(targets.columns)
    ids = [k.replace("ko:", "") for k in kos]
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    within = targets.groupby(ks["cruise"].reindex(targets.index).to_numpy()).transform(lambda x: x - x.mean())
    sd = within.std().to_numpy()
    sd = np.where(sd > 0, sd, 1)

    lib = pd.read_csv(LIBDIR / "library_metadata.csv")
    lib = lib[lib["in_ko_table"] & lib["sample_id"].notna()].drop_duplicates("lib_id").set_index("lib_id")
    clr = rarefied_clr(lib, kos)
    lib = lib.loc[clr.index]
    depth = pd.read_csv(REPO / "data" / "manifests" / "biogoship_samples.csv").drop_duplicates("sample_id").set_index("sample_id")["depth_m"]
    faire = pd.read_csv(REPO / "data" / "faire" / "faire_samples.csv").drop_duplicates("samp_name").set_index("samp_name")["depth_min_m"]
    surface = (lib["sample_id"].map(depth).fillna(lib["samp_name"].map(faire)).fillna(0) <= 10).to_numpy()
    fc = lib["flowcell_id"].astype(str)
    eff = {}

    # HiSeq - NovaSeq within cruise
    num, den = np.zeros(len(kos)), 0.0
    for cruise in ["I09N-2016", "AE1319-2013", "BVAL46-2011", "NH1418-2014"]:
        m = (lib["cruise"] == cruise).to_numpy() & surface
        hi, no = m & fc.str.endswith("BBXX").to_numpy(), m & fc.str.contains("DSX").to_numpy()
        if hi.sum() < 2 or no.sum() < 2:
            continue
        wgt = min(hi.sum(), no.sum())
        num += wgt * (clr[hi].mean().to_numpy() - clr[no].mean().to_numpy())
        den += wgt
    eff["hiseq"] = num / den / sd

    # era: NovaSeq 6000 cruises - LT4 cruises
    mod = lib["in_real_targets"].to_numpy() & surface
    cm = clr[mod].groupby(lib.loc[mod, "cruise"].to_numpy()).mean()
    seq = lib.loc[mod].groupby("cruise")["flowcell_id"].agg(
        lambda v: "LT4" if v.astype(str).str.endswith("LT4").all() else ("6000" if v.astype(str).str.contains("DS|DR").all() else "mixed"))
    eff["era"] = (cm[seq.reindex(cm.index) == "6000"].mean() - cm[seq.reindex(cm.index) == "LT4"].mean()).to_numpy() / sd

    # HLY2202 rerun: same sample, two NovaSeq 6000 flowcells
    h = lib[(lib["cruise"] == "HLY2202-2022")]
    a = h[h["flowcell_id"] == "HMY7WDSX7"].groupby("sample_id").apply(lambda g: g.index[0])
    b = h[h["flowcell_id"] == "HFJLFDSX7"].groupby("sample_id").apply(lambda g: g.index[0])
    both = a.index.intersection(b.index)
    eff["rerun"] = (clr.loc[a[both]].to_numpy() - clr.loc[b[both]].to_numpy()).mean(0) / sd
    say(f"HLY2202 rerun pairs: {len(both)}")

    # paired protocol contrasts on the sample-level targets
    allc = pd.read_csv(REPO / "data" / "targets" / "ko_clr_incl_protocolA.csv", index_col="sample_key")
    pairs = pd.read_csv(REPO / "transfer" / "results" / "axis2_protocol" / "pairs_list.csv")
    for contrast, name in (("A16N Protocol A - B", "protocolA"), ("C13.5 v2 - v1", "c13_v2")):
        g = pairs[pairs["contrast"] == contrast]
        d = allc.reindex(columns=kos).loc[g["treated"]].to_numpy() - allc.reindex(columns=kos).loc[g["partner"]].to_numpy()
        eff[name] = np.nanmean(d, 0) / sd
        say(f"{contrast}: {len(g)} pairs")

    load = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "loadings.csv.gz", index_col=0).loc[ids]
    table = pd.DataFrame({c: load[c].to_numpy() for c in ["F1", "F2", "F3", "W1", "W2", "W3"]}, index=ids)
    for k, v in eff.items():
        table[k] = v
    gc = pd.read_csv(OUT / "ko_gc.csv").set_index("ko")
    gc["euk_share"] = gc["genes_Eukaryota"] / gc["genes"]
    table = table.join(gc[["genes", "gc_mean", "gc_mean_Bacteria", "genes_Bacteria", "gc_mean_Eukaryota", "euk_share", "length_mean"]])
    say(f"KOs in the real model: {len(ids):,}; with GC: {table['gc_mean'].notna().sum():,} "
        f"(median sampled genes {table['genes'].median():.0f}); GC median {table['gc_mean'].median():.3f}, "
        f"IQR {table['gc_mean'].quantile(.25):.3f}-{table['gc_mean'].quantile(.75):.3f}")
    table.round(5).to_csv(OUT / "gc_effects.csv", index_label="ko")

    # correlations between the effects themselves (for context)
    effs = ["F1", "W1", "F2", "W2", "F3", "W3", "hiseq", "era", "rerun", "protocolA", "c13_v2"]
    say("\nSpearman between effect vectors (KOs): hiseq/era/rerun vs F1/W1")
    for e in ["hiseq", "era", "rerun", "protocolA"]:
        say(f"   {e:10s} " + "  ".join(f"{c} {spearmanr(table[e], table[c], nan_policy='omit')[0]:+.2f}" for c in ["F1", "W1", "F3", "hiseq", "era"]))

    rows = []
    ok = table.dropna(subset=["gc_mean", "euk_share", "length_mean"])
    okb = table.dropna(subset=["gc_mean_Bacteria"])
    okb = okb[okb["genes_Bacteria"] >= 20]
    dec0 = pd.qcut(ok["gc_mean"], 10, labels=False)
    say("\nGC vs per-KO effect (Spearman rho; partial = controlling eukaryote share and gene length; bacterial-gene GC)")
    say(f"   {'effect':10s} {'rho GC':>8s} {'partial':>8s} {'rho bactGC':>11s}   n {len(ok):,} / {len(okb):,}")
    for e in effs:
        o, ob = ok.dropna(subset=[e]), okb.dropna(subset=[e])
        r, p = spearmanr(o["gc_mean"], o[e])
        pr = partial_spearman(o["gc_mean"].to_numpy(), o[e].to_numpy(), [o["euk_share"].to_numpy(), o["length_mean"].to_numpy()])
        rb, pb = spearmanr(ob["gc_mean_Bacteria"], ob[e])
        re, _ = spearmanr(o["euk_share"], o[e])
        rows.append({"effect": e, "rho_gc": r, "p_gc": p, "partial_rho_gc": pr, "rho_bacterial_gc": rb, "p_bacterial_gc": pb,
                     "rho_euk_share": re, "n": len(o), "n_bacterial": len(ob)})
        say(f"   {e:10s} {r:+8.3f} {pr:+8.3f} {rb:+11.3f}   (euk share rho {re:+.2f})")
    say(f"   GC vs eukaryote share: rho {spearmanr(ok['gc_mean'], ok['euk_share'])[0]:+.3f}")

    # does GC mediate the resemblance of the technical contrasts to axis 1?
    say("\nResemblance to the axis-1 loading with and without GC (Spearman; partial on GC, eukaryote share, length)")
    cov = [ok["gc_mean"].to_numpy(), ok["euk_share"].to_numpy(), ok["length_mean"].to_numpy()]
    from sklearn.preprocessing import SplineTransformer

    basis = SplineTransformer(n_knots=6, degree=3).fit_transform(ok[["gc_mean"]])
    for e in ["hiseq", "era", "rerun"]:
        for c in ["F1", "W1"]:
            r0 = spearmanr(ok[e], ok[c])[0]
            r1 = partial_spearman(ok[e].to_numpy(), ok[c].to_numpy(), cov)
            say(f"   {e:6s} vs {c}: rho {r0:+.3f}; partial {r1:+.3f}")
        y = ok[e].to_numpy()
        fit = basis @ np.linalg.lstsq(basis, y, rcond=None)[0]
        say(f"   {e:6s}: share of its per-KO variance explained by a GC spline: {1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum():.3f}; "
            f"RMS {np.sqrt((y ** 2).mean()):.3f} SD; decile range {np.ptp(ok.groupby(dec0)[e].mean().to_numpy()):.3f}")
    pd.DataFrame(rows).round(4).to_csv(OUT / "gc_correlations.csv", index=False)

    dec = pd.qcut(ok["gc_mean"], 10, labels=False)
    say("\nMean effect per GC decile (KO SD units; loadings as is)")
    md = ok.groupby(dec)[["gc_mean"] + effs].mean()
    say(md.round(3).to_string())
    md.round(4).to_csv(OUT / "gc_deciles.csv")
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")
    figure(ok, effs)


def figure(ok, effs):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6})
    panels = [("hiseq", "HiSeq 4000 − NovaSeq 6000\n(within cruise)"), ("rerun", "Same DNA, two NovaSeq 6000\nruns (control)"),
              ("era", "NovaSeq 6000 − NovaSeq X\ncruises"), ("F1", "Axis 1 loading\n(F1)")]
    fig, axes = plt.subplots(1, 4, figsize=(7.4, 2.3), sharex=True)
    bins = np.quantile(ok["gc_mean"], np.linspace(0, 1, 21))
    mid = ok.groupby(pd.cut(ok["gc_mean"], bins, include_lowest=True), observed=True)
    for ax, (e, lab) in zip(axes, panels):
        ax.scatter(ok["gc_mean"] * 100, ok[e], s=2, color="#999999", alpha=0.25, linewidths=0, rasterized=True)
        m = mid[["gc_mean", e]].mean()
        ax.plot(m["gc_mean"] * 100, m[e], color="#0072B2", lw=1.4)
        ax.axhline(0, color="#333333", lw=0.5)
        lo, hi = np.nanpercentile(ok[e], [1, 99])
        ax.set_ylim(lo, hi)
        ax.set_xlabel("Gene GC (%)")
        ax.set_ylabel(lab + ("" if e == "F1" else " (SD)"), fontsize=7)
        ax.text(0.03, 0.97, f"ρ = {spearmanr(ok['gc_mean'], ok[e])[0]:+.2f}", transform=ax.transAxes, va="top", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "ko_gc_effects.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
