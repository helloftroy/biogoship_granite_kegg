#!/usr/bin/env python3
"""R67: which domain carries each KO in our samples? Read-based domain share vs the reference-gene share.

Until now a KO's "eukaryote share" (R44 transfer/results/ko_gc/ko_gc.csv; used in R48 / R51 classes and R66) was
the eukaryotic fraction of a 1-in-10 sample of the KO's genes in the Bio-GO-SHIP 2026 reference catalogue
(non-redundant cluster representatives from several sources). It counts distinct catalogue genes, not what is in
the water.

Here each KO's reads are split by the domain of the reference sequence they map to, from the per-library
phylum-level per-taxon KO counts (R50, data/targets/taxa/ko_by_taxon/<lib>/level_2_phylum_division_ko_counts.tsv.gz;
their sum over taxa is ko_counts_cov50_all.tsv). Libraries and samples: as build_ko_targets.py /
domain_split_targets.py (cov50, surface, non-control, >= 1e5 KO reads), restricted to the hurdle samples.

Per KO:
  reads, reads_prok (Bacteria + Archaea), reads_euk, reads_vir      pooled over all libraries
  euk_read_share       reads_euk / reads (pooled)
  euk_read_share_med   median over samples with >= MIN reads of that sample's eukaryotic share
  euk_read_share_iqr   IQR of that per-sample share (is the carrier domain stable across samples?)
  n_samples_min        samples with >= MIN reads
  euk_gene_share       R44 reference share, Eukaryota / (Eukaryota + Bacteria + Archaea) (as R48)

Outputs: data/results/ko_domain/ko_domain_reads.csv, library_domain.npz (cache), summary.txt,
figures/output/ko_domain_reads.{png,pdf}
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO))
TAXDIR = REPO / "data" / "targets" / "taxa" / "ko_by_taxon"
OUT = REPO / "data" / "results" / "ko_domain"
FIG = REPO / "figures" / "output"
MIN = 10  # reads for a per-sample share
LOG: list = []


def say(s=""):
    print(s, flush=True)
    LOG.append(str(s))


def split_one(args):
    lib, kos = args
    k = pd.read_csv(TAXDIR / lib / "level_2_phylum_division_ko_counts.tsv.gz", sep="\t", usecols=["taxon", "ko", "count"])
    k["ko"] = k["ko"].str.replace("ko:", "", regex=False)
    dom = k["taxon"].str.split(";").str[0]
    k["dom"] = np.where(dom.isin(["Bacteria", "Archaea"]), "prok", np.where(dom == "Eukaryota", "euk",
                        np.where(dom == "Viruses", "vir", "other")))  # fmt: skip
    g = k.groupby(["dom", "ko"])["count"].sum()
    return lib, {d: g.get(d, pd.Series(dtype=float)).reindex(kos).fillna(0).to_numpy(np.float32) for d in ("prok", "euk", "vir", "other")}


def build(jobs):
    from biogoship.tables import MANIFEST_DIR, REPO_ROOT
    from build_ko_targets import CONTROL, library_samples

    head = open(REPO / "ko_counts_cov50_all.tsv").readline().rstrip("\n").split("\t")
    kos = [c.replace("count.ko.", "") for c in head]
    counts = pd.read_csv(REPO / "ko_counts_cov50_all.tsv", sep="\t", header=None, skiprows=1, index_col=0)
    counts.columns = kos
    samples = pd.read_pickle(REPO / "data/results/hurdle/sample_counts.pkl").index
    ns = argparse.Namespace(workbook=REPO_ROOT / "FAIRe-NOAA_uci-ess-biogoship.xlsx",
                            faire_samples=REPO / "data" / "faire" / "faire_samples.csv", manifest_dir=MANIFEST_DIR)  # fmt: skip
    info = library_samples(counts.index, ns).join(pd.Series(counts.sum(axis=1), name="total"))
    keep = (info["samp_name"].notna() & ~info["samp_name"].fillna("").str.contains(CONTROL) & info["cruise"].notna()
            & info["sample_key"].isin(samples) & (info["total"] >= 1e5) & ~(info["depth_m"] > 10))  # fmt: skip
    info = info[keep]
    libs = list(info.index)
    say(f"libraries {len(libs):,} for {info['sample_key'].nunique():,} of {len(samples):,} hurdle samples; KOs {len(kos):,}")
    cache = OUT / "library_domain.npz"
    if cache.exists() and list(np.load(cache, allow_pickle=True)["libs"]) == libs:
        z = np.load(cache, allow_pickle=True)
        M = {d: z[d] for d in ("prok", "euk", "vir", "other")}
    else:
        M = {d: np.zeros((len(libs), len(kos)), np.float32) for d in ("prok", "euk", "vir", "other")}
        pos = {lib: i for i, lib in enumerate(libs)}
        with ProcessPoolExecutor(jobs) as ex:
            for j, (lib, r) in enumerate(ex.map(split_one, [(lib, kos) for lib in libs], chunksize=8)):
                for d in M:
                    M[d][pos[lib]] = r[d]
                if (j + 1) % 250 == 0:
                    print(f"  {j + 1:,}/{len(libs):,}", flush=True)
        np.savez_compressed(cache, libs=np.array(libs, dtype=object), kos=np.array(kos, dtype=object), **M)
    tot = counts.loc[libs, kos].to_numpy(np.float64)
    split = M["prok"].astype(np.float64) + M["euk"] + M["vir"] + M["other"]
    say(f"check: level-2 domain split vs official table, max |diff| {np.abs(split - tot).max():.4f} "
        f"(total reads {tot.sum():,.0f}); domain shares of all KO reads: prok {M['prok'].sum() / tot.sum():.4f}, "
        f"euk {M['euk'].sum() / tot.sum():.4f}, viruses {M['vir'].sum() / tot.sum():.4f}, other {M['other'].sum() / tot.sum():.4f}")  # fmt: skip
    # per KO
    r = pd.DataFrame(index=pd.Index(kos, name="ko"))
    r["reads"] = tot.sum(0)
    for d in ("prok", "euk", "vir"):
        r[f"reads_{d}"] = M[d].sum(0)
    r["euk_read_share"] = (r["reads_euk"] / r["reads"]).clip(0, 1)  # float32 sums can put all-eukaryotic KOs a hair above 1
    # per sample (libraries summed)
    sk = info["sample_key"].to_numpy()
    E = pd.DataFrame(M["euk"], index=sk).groupby(level=0).sum()
    T = pd.DataFrame(tot, index=sk).groupby(level=0).sum()
    sh = (E / T).clip(0, 1).where(T >= MIN)
    r["n_samples_min"] = (T >= MIN).sum().to_numpy()
    r["euk_read_share_med"] = sh.median().to_numpy()
    r["euk_read_share_iqr"] = (sh.quantile(0.75) - sh.quantile(0.25)).to_numpy()
    gc = pd.read_csv(REPO / "transfer/results/ko_gc/ko_gc.csv", index_col="ko")
    r["euk_gene_share"] = (gc["genes_Eukaryota"] / gc[["genes_Eukaryota", "genes_Bacteria", "genes_Archaea"]].sum(1)).reindex(r.index)
    r["ref_genes"] = gc["genes"].reindex(r.index)
    names = pd.read_csv(REPO / "data/reference/kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "name"])
    r["name"] = names.assign(ko=names["ko"].str.replace("ko:", "")).set_index("ko")["name"].reindex(r.index)
    r.to_csv(OUT / "ko_domain_reads.csv", float_format="%.6g")
    return r


def bins(x):
    return pd.cut(x, [-1e-9, 0.05, 0.1, 0.5, 0.9, 0.95, 1.0], labels=["≤0.05", "0.05–0.1", "0.1–0.5", "0.5–0.9", "0.9–0.95", ">0.95"])


def summarise(r):
    hur = pd.read_csv(REPO / "data/results/hurdle/hurdle_all.csv", index_col="ko").index
    core = pd.read_csv(REPO / "data/results/benchmark_export/targets/clr_all.csv.gz", nrows=0).columns.drop("tile_key")
    cls = pd.read_csv(REPO / "data/targets/domain_split/ko_classes.csv", index_col="ko")["class"]
    for nm, idx in (("all KOs with reads", r.index[r["reads"] > 0]), ("core 3,933", core), ("hurdle 5,808", hur)):
        g = r.loc[idx].dropna(subset=["euk_gene_share"])
        rho = spearmanr(g["euk_gene_share"], g["euk_read_share"])[0]
        say(f"\n## {nm}: n {len(g):,} (with a reference share); Spearman(reference share, read share) {rho:+.3f}")
        say(f"  median reference share {g['euk_gene_share'].median():.3f}; median read share {g['euk_read_share'].median():.3f}; "
            f"read share pooled over these KOs' reads {g['reads_euk'].sum() / g['reads'].sum():.3f}")  # fmt: skip
        say("  rows = reference-gene share, columns = read share (pooled):\n"
            + pd.crosstab(bins(g["euk_gene_share"]), bins(g["euk_read_share"]), margins=True).to_string())  # fmt: skip
        hi, lo = g["euk_gene_share"] > 0.9, g["euk_gene_share"] <= 0.1
        say(f"  reference > 0.9 ({hi.sum():,}): read share median {g.loc[hi, 'euk_read_share'].median():.3f}; "
            f"> 0.9 by reads {(g.loc[hi, 'euk_read_share'] > 0.9).mean():.1%}; <= 0.1 by reads {(g.loc[hi, 'euk_read_share'] <= 0.1).mean():.1%}")  # fmt: skip
        say(f"  reference <= 0.1 ({lo.sum():,}): read share median {g.loc[lo, 'euk_read_share'].median():.3f}; "
            f"<= 0.1 by reads {(g.loc[lo, 'euk_read_share'] <= 0.1).mean():.1%}; > 0.9 by reads {(g.loc[lo, 'euk_read_share'] > 0.9).mean():.1%}")  # fmt: skip
    say("\n## R48 / R51 classes (core KOs) by read share")
    g = r.loc[cls.index].assign(cls=cls)
    say(g.groupby("cls")["euk_read_share"].describe(percentiles=[.1, .5, .9]).round(3).to_string())
    say("  read share of pooled reads per class: " + ", ".join(
        f"{c} {gg['reads_euk'].sum() / gg['reads'].sum():.3f}" for c, gg in g.groupby("cls")))  # fmt: skip
    say("\n## Stability across samples (KOs with >= 100 samples of >= 10 reads)")
    s = r[r["n_samples_min"] >= 100]
    say(f"  n {len(s):,}; IQR of per-sample eukaryotic read share: median {s['euk_read_share_iqr'].median():.3f}, 90% {s['euk_read_share_iqr'].quantile(.9):.3f}")
    mid = s[(s["euk_read_share_med"] > 0.1) & (s["euk_read_share_med"] < 0.9)]
    say(f"  KOs with a mixed median share (0.1–0.9): {len(mid):,}; their IQR median {mid['euk_read_share_iqr'].median():.3f}")
    say("\n## Sanity: known KOs")
    known = {"K02588": "nifH", "K02586": "nifD", "K02591": "nifK", "K10944": "amoA", "K10945": "amoB", "K00399": "mcrA",
             "K06163": "phnJ", "K00372": "nasA", "K02703": "psbA", "K02689": "psaA", "K01601": "rbcL",
             "K02950": "rpsL (prok RP)", "K02886": "rplB (prok RP)", "K07374": "TUBA", "K07375": "TUBB",
             "K05692": "ACTB_G", "K11253": "H3", "K11254": "H4", "K04875": "KCNA2", "K04877": "KCNA4",
             "K15670": "rifA", "K18925": "ykkD", "K09431": "ETV1", "K13871": "SLC7A14", "K08671": "PCSK4",
             "K06425": "sspH", "K07785": "nrsD"}  # fmt: skip
    k = r.reindex(list(known)).assign(gene=list(known.values()))
    k["in_hurdle"] = k.index.isin(hur)
    say(k[["gene", "reads", "euk_gene_share", "euk_read_share", "euk_read_share_med", "euk_read_share_iqr", "n_samples_min", "in_hurdle"]]
        .round(3).to_string())  # fmt: skip
    (OUT / "summary.txt").write_text("\n".join(LOG) + "\n")


def figure(r):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})  # fmt: skip
    hur = set(pd.read_csv(REPO / "data/results/hurdle/hurdle_all.csv")["ko"])
    g = r[(r["reads"] > 0) & r["euk_gene_share"].notna()]
    fig, ax = plt.subplots(figsize=(3.5, 3.3))
    plt.subplots_adjust(left=0.15, right=0.97, top=0.97, bottom=0.14)
    o = ~g.index.isin(hur)
    ax.scatter(g.loc[o, "euk_gene_share"], g.loc[o, "euk_read_share"], s=1.5, lw=0, color="#b0b0b0", alpha=0.5, rasterized=True,
               label=f"Other KOs ({o.sum():,})")  # fmt: skip
    ax.scatter(g.loc[~o, "euk_gene_share"], g.loc[~o, "euk_read_share"], s=1.5, lw=0, color="#0072B2", alpha=0.5, rasterized=True,
               label=f"Patchy (hurdle) KOs ({(~o).sum():,})")  # fmt: skip
    ax.plot([0, 1], [0, 1], color="#555555", lw=0.6, ls=(0, (3, 2)))
    lab = {"K02703": "psbA", "K01601": "rbcL", "K02950": "rpsL", "K04875": "KCNA2", "K13871": "SLC7A14", "K15670": "rifA",
           "K07374": "TUBA", "K02588": "nifH"}  # fmt: skip
    for k, nm in lab.items():
        if k in g.index:
            x, y = g.loc[k, "euk_gene_share"], g.loc[k, "euk_read_share"]
            ax.scatter(x, y, s=14, facecolor="none", edgecolor="black", lw=0.7, zorder=4)
            ax.annotate(nm, (x, y), xytext=(-4, 4) if x > 0.5 else (4, 4), ha="right" if x > 0.5 else "left",
                        textcoords="offset points", fontsize=6, style="italic")  # fmt: skip
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Eukaryotic share of the KO's reference-catalogue genes")
    ax.set_ylabel("Eukaryotic share of the KO's reads in our samples")
    ax.legend(frameon=False, fontsize=6, loc="upper left", markerscale=4)
    for suf in ("png", "pdf"):
        fig.savefig(FIG / f"ko_domain_reads.{suf}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--figure-only", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if "--figure-only" in sys.argv:
        figure(pd.read_csv(OUT / "ko_domain_reads.csv", index_col="ko"))
        sys.exit()
    r = build(a.jobs)
    summarise(r)
    figure(r)
