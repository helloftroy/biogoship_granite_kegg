#!/usr/bin/env python3
"""Do plastid (and mitochondrial) copies inflate the prokaryote genome-equivalent denominator? (RESULTS.md R49)

The R46 denominator pools 54 prokaryote-type ribosomal-protein KOs. Plastid genomes encode many of these
(red-lineage plastids, e.g. diatoms and haptophytes, nearly all; green plastids fewer), and nucleus-encoded
organellar ribosomal proteins share the bacterial-type KOs. Plastids are high-copy, so in eukaryote-rich
samples "prokaryote genome equivalents" could be inflated. Cyanobacteria are affected the same way by
nature (plastids descend from them) but are genuine prokaryotes.

1. Literature split of the 54 markers (gene symbol):
   plastid-encoded in at least one major algal lineage (green or red plastid genomes; union):
     L1-L6, L11, L7/L12, L13, L14, L16, L18-L24, L27, L29, L31-L36, S2-S20
   rarely or never plastid-encoded: L9, L10, L15, L17, L25, L28, L30, S1, S21
   (a union of lineages, so a KO in the first group need not be plastid-encoded in a given alga; the
   empirical split in 2 is the data's own answer).
2. Empirical, from the per-taxon KO counts (data/targets/taxa/ko_by_taxon, phylum level, 1,379 libraries):
   per marker KO, the share of its reads assigned to eukaryotic taxa (plastid, mitochondrial and nuclear
   organellar copies in eukaryotic reference bins) and to Cyanobacteriota; overall and in eukaryote-rich
   libraries (top quartile of the domain eukaryote fraction). Markers whose eukaryote-assigned share stays
   below --max-euk (default 2 %) in the eukaryote-rich quartile form the "clean" subset.
3. Divergence test on all libraries (ko_counts_cov50_all.tsv, non-blank, > 100,000 KO reads):
   per subset S, log denominator D_S = mean over KOs in S of log(count + 0.5). The difference D_A - D_B
   between two subsets is the log ratio of the genome equivalents they imply (up to a constant).
   Spearman of D_A - D_B with the eukaryote fraction (domain reads) and with the Cyanobacteriota share of
   reads (phylum_counts.tsv), overall and within cruise; mean difference by eukaryote-fraction decile
   (how much the plastid-prone set inflates genome equivalents, in %).
4. Truth check where per-taxon counts exist: D_prok-assigned = log of the 54 markers' reads assigned to
   Bacteria + Archaea only; bias of D_all and D_clean against it as the eukaryote fraction rises.

Outputs: data/results/plastid_markers/markers.csv, library_denominators.csv, deciles.csv, summary.txt;
figures/output/plastid_marker_check.png.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/plastid_marker_check.py
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "results" / "plastid_markers"
FIG = REPO / "figures" / "output"
NON_PLASTID = {"L9", "L10", "L15", "L17", "L25", "L28", "L30", "S1", "S21"}


def protein(name: str) -> str:
    m = re.search(r"(large|small) subunit ribosomal protein (L7/L12|[LS]\d+)", name)
    return m.group(2) if m else ""


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--max-euk", type=float, default=0.02)
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    names = pd.read_csv(REPO / "data" / "reference" / "kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "name"])
    names = names.assign(ko=names["ko"].str.replace("ko:", "")).set_index("ko")["name"]
    targets = [k.replace("ko:", "") for k in pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", nrows=0).columns.drop("sample_key")]
    sym = names.reindex(targets).fillna("").str.split(";").str[0]
    mk = [k for k in targets if re.search(r"\b(rpl[A-Z]|rps[A-Z]|rpm[A-Z])\b", sym[k]) and "ribosomal protein" in names.get(k, "")]
    m = pd.DataFrame({"ko": mk, "name": [names[k] for k in mk]})
    m["protein"] = m["name"].map(protein)
    m["literature"] = np.where(m["protein"].isin(NON_PLASTID), "non-plastid", "plastid-encoded (some lineage)")
    say(f"markers: {len(m)}; literature split: " + ", ".join(f"{k} {v}" for k, v in m["literature"].value_counts().items()))

    # domain fraction, cyanobacterial share
    dom = pd.read_csv(REPO / "data" / "targets" / "taxa" / "domain_counts.tsv", sep="\t", index_col=0)
    dsum = dom.sum(1)
    euk_frac = dom["Eukaryota"] / dsum
    phy = pd.read_csv(REPO / "data" / "targets" / "taxa" / "phylum_counts.tsv", sep="\t", index_col=0)
    cyano = phy.filter(like="Cyanobacteri").sum(1) / phy.sum(1)

    # 2. empirical per-KO eukaryote-assigned share
    stash = REPO / "data" / "targets" / "taxa" / "ko_by_taxon"
    want = {"ko:" + k for k in mk}
    rows = []
    for d in sorted(stash.iterdir()):
        f = d / "level_2_phylum_division_ko_counts.tsv.gz"
        if not f.exists():
            continue
        k = pd.read_csv(f, sep="\t", usecols=["taxon", "ko", "count"])
        k = k[k["ko"].isin(want)]
        k["dom"] = k["taxon"].str.split(";").str[0]
        k["cy"] = k["taxon"].str.contains("Cyanobacteri")
        k["total"] = k["count"]
        k["euk"] = k["count"].where(k["dom"] == "Eukaryota", 0)
        k["prok"] = k["count"].where(k["dom"].isin(["Bacteria", "Archaea"]), 0)
        k["cyano"] = k["count"].where(k["cy"], 0)
        g = k.groupby("ko")[["total", "euk", "prok", "cyano"]].sum()
        g["lib"] = d.name
        rows.append(g.reset_index())
    kt = pd.concat(rows)
    kt["ko"] = kt["ko"].str.replace("ko:", "")
    kt["euk_frac"] = kt["lib"].map(euk_frac)
    libs_t = kt["lib"].unique()
    q75 = euk_frac.reindex(libs_t).quantile(0.75)
    say(f"per-taxon libraries: {len(libs_t):,}; eukaryote-fraction top quartile cut {q75:.3f}")
    agg = []
    for ko, g in kt.groupby("ko"):
        rich = g[g["euk_frac"] >= q75]
        agg.append({"ko": ko, "euk_share_all": g["euk"].sum() / g["total"].sum(), "euk_share_rich": rich["euk"].sum() / rich["total"].sum(),
                    "cyano_share_all": g["cyano"].sum() / g["total"].sum(),
                    "euk_share_rich_median_lib": (rich["euk"] / rich["total"]).median()})
    m = m.merge(pd.DataFrame(agg), on="ko", how="left")
    m["empirical"] = np.where(m["euk_share_rich"] < args.max_euk, "clean", "eukaryote-carried")
    m = m.sort_values("euk_share_rich")
    m.round(4).to_csv(OUT / "markers.csv", index=False)
    say(m[["ko", "protein", "literature", "euk_share_all", "euk_share_rich", "cyano_share_all", "empirical"]].round(3).to_string(index=False))
    say("\neukaryote-assigned share in eukaryote-rich libraries, by literature class (median over KOs): "
        + ", ".join(f"{c} {g['euk_share_rich'].median():.3f}" for c, g in m.groupby("literature")))
    say(f"empirical split at {args.max_euk:.0%}: " + ", ".join(f"{k} {v}" for k, v in m["empirical"].value_counts().items()))

    # 3. divergence on all libraries
    h = [c.replace("count.ko.", "") for c in open(REPO / "ko_counts_cov50_all.tsv").readline().rstrip("\n").split("\t")]
    col = {k: i for i, k in enumerate(h)}
    sets = {"all": list(m["ko"]), "lit_plastid": list(m.loc[m["literature"] != "non-plastid", "ko"]),
            "lit_nonplastid": list(m.loc[m["literature"] == "non-plastid", "ko"]),
            "emp_clean": list(m.loc[m["empirical"] == "clean", "ko"]), "emp_euk": list(m.loc[m["empirical"] != "clean", "ko"])}
    idx = {s: [col[k] for k in v] for s, v in sets.items()}
    D, L, T = [], [], []
    with open(REPO / "ko_counts_cov50_all.tsv") as fh:
        fh.readline()
        for line in fh:
            n, r = line.split("\t", 1)
            v = np.fromstring(r, sep="\t")
            D.append([np.log(v[i] + 0.5).mean() for i in idx.values()])
            L.append(n)
            T.append(v.sum())
    den = pd.DataFrame(D, index=L, columns=list(idx))
    den["ko_reads"] = T
    den = den[~den.index.str.contains("BLANK|NTC|BLNK|CTRL|MOCK", case=False) & (den["ko_reads"] > 1e5)]
    den["euk_frac"] = euk_frac.reindex(den.index)
    den["cyano"] = cyano.reindex(den.index)
    lib = pd.read_csv(REPO / "transfer" / "results" / "library_metadata" / "library_metadata.csv").drop_duplicates("lib_id").set_index("lib_id")
    den["cruise"] = lib["cruise"].reindex(den.index)
    # truth (prokaryote-assigned reads of all 54 markers), where per-taxon counts exist
    pk = kt.pivot_table(index="lib", columns="ko", values="prok", aggfunc="sum").reindex(columns=m["ko"]).fillna(0)
    den["prok_assigned"] = np.log(pk + 0.5).mean(1).reindex(den.index)
    den = den.dropna(subset=["euk_frac"])
    den.round(5).to_csv(OUT / "library_denominators.csv", index_label="lib_id")

    def rel(a, b):
        x = den[a] - den[b]
        xc = x - x.groupby(den["cruise"]).transform("mean")
        ec = den["euk_frac"] - den["euk_frac"].groupby(den["cruise"]).transform("mean")
        cc = den["cyano"] - den["cyano"].groupby(den["cruise"]).transform("mean")
        return x, spearmanr(x, den["euk_frac"], nan_policy="omit")[0], spearmanr(xc, ec, nan_policy="omit")[0], spearmanr(x, den["cyano"], nan_policy="omit")[0], spearmanr(xc, cc, nan_policy="omit")[0]

    say(f"\n== Divergence of denominators with the eukaryote fraction ({len(den):,} libraries)")
    dec = pd.qcut(den["euk_frac"], 10, labels=False, duplicates="drop")
    dtab = {"euk_frac": den.groupby(dec)["euk_frac"].median()}
    for a, b in (("lit_plastid", "lit_nonplastid"), ("emp_euk", "emp_clean"), ("all", "emp_clean"), ("all", "lit_nonplastid"),
                 ("all", "prok_assigned"), ("emp_clean", "prok_assigned"), ("lit_nonplastid", "prok_assigned")):
        x, r, rw, rc, rcw = rel(a, b)
        dm = x.groupby(dec).mean()
        dtab[f"{a} - {b}"] = dm - dm.iloc[0]
        say(f"   {a:15s} - {b:15s}: rho with euk fraction {r:+.3f} (within cruise {rw:+.3f}); with cyano share {rc:+.3f} "
            f"(within cruise {rcw:+.3f}); top - bottom decile {dm.iloc[-1] - dm.iloc[0]:+.3f} log = {np.expm1(dm.iloc[-1] - dm.iloc[0]):+.0%}"
            + (f"  (n {x.notna().sum():,})" if b == "prok_assigned" or a == "prok_assigned" else ""))
    dt = pd.DataFrame(dtab)
    dt.round(4).to_csv(OUT / "deciles.csv")
    say("\n   mean log difference by eukaryote-fraction decile, relative to the lowest decile")
    say(dt.round(3).to_string())
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")
    figure(den)


def figure(den):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6})
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.6))
    x = np.log10(den["euk_frac"].clip(lower=1e-3))
    for ax, (a, b, lab) in zip(axes, (("lit_plastid", "lit_nonplastid", "Plastid-encodable − non-plastid\nribosomal proteins (log)"),
                                      ("all", "prok_assigned", "All 54 markers − prokaryote-\nassigned reads only (log)"))):
        y = den[a] - den[b]
        ok = y.notna()
        ax.scatter(x[ok], y[ok] - y[ok].median(), s=3, color="#0072B2", alpha=0.3, linewidths=0, rasterized=True)
        bins = pd.qcut(x[ok], 15, duplicates="drop")
        mm = pd.DataFrame({"x": x[ok], "y": y[ok] - y[ok].median()}).groupby(bins, observed=True).mean()
        ax.plot(mm["x"], mm["y"], color="#D55E00", lw=1.4)
        ax.axhline(0, color="#999999", lw=0.5)
        ax.set_xlabel("Eukaryote fraction of mapped reads (log10)")
        ax.set_ylabel(lab, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "plastid_marker_check.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
