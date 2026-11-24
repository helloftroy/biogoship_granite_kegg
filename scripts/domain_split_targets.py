#!/usr/bin/env python3
"""Domain-split KO targets A / B / C / D (RESULTS.md R51; design R48, denominator R49).

Per library (cov50, the libraries behind data/targets/ko_clr.csv), from the phylum-level per-taxon KO counts
(data/targets/taxa/ko_by_taxon/<lib>/level_2_phylum_division_ko_counts.tsv.gz; their sum over taxa is
ko_counts_cov50_all.tsv), each of the 3,933 model KOs is split into reads assigned to Bacteria + Archaea
("prok") and to Eukaryota ("euk"). Proportions are per library (of its total KO reads), libraries of one
sample averaged (as scripts/build_ko_targets.py); pseudocount 0.5 / sample total.

KO classes (5 % / 95 % eukaryote share of reference genes, transfer/results/ko_gc/ko_gc.csv, R48):
prokaryote-specific, eukaryote-specific, shared. The 54 prokaryote-type ribosomal-protein markers are the
denominator and not targets.

Targets (sample x KO, log scale):
  A   prokaryote-specific KOs: log copies per prokaryote genome equivalent
        = log(p_k) - mean over the 54 markers of log(prok-assigned marker proportion)        (R49 fix)
  B   eukaryote-specific KOs: within-eukaryote composition
        = log(p_k) - log(pooled eukaryote-specific proportion over all 8,818 eukaryote-specific KOs of the
          full table); samples with < --min-euk (5,000) pooled eukaryote-specific reads dropped
  C   eukaryotic amount (one column): log(eukaryote-assigned / prokaryote-assigned KO reads), all KOs (R68).
        The R51 version, log(pooled eukaryote-specific / pooled prokaryote-specific reads) with classes from the
        reference-gene share, is written as C_refclass.csv.
  D   shared KOs, whole community: the existing CLR values (ko_clr.csv columns)
  Dp  shared KOs, prokaryote-assigned reads per prokaryote genome equivalent (as A)
  De  shared KOs, eukaryote-assigned reads within eukaryotes (as B)
  Dp / De keep KOs whose assigned part is non-zero in >= 90 % of samples.

Outputs in data/targets/domain_split/: A.csv B.csv C.csv D.csv Dp.csv De.csv (sample_key + columns),
ko_classes.csv, library_split.npz (cache of the per-library prok / euk counts), summary.txt.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/domain_split_targets.py
"""

from __future__ import annotations

import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO))
OUT = REPO / "data" / "targets" / "domain_split"
TAXDIR = REPO / "data" / "targets" / "taxa" / "ko_by_taxon"


def split_one(args):
    lib, kos = args
    f = TAXDIR / lib / "level_2_phylum_division_ko_counts.tsv.gz"
    k = pd.read_csv(f, sep="\t", usecols=["taxon", "ko", "count"])
    k["ko"] = k["ko"].str.replace("ko:", "", regex=False)
    dom = k["taxon"].str.split(";").str[0]
    k = k.assign(dom=np.where(dom.isin(["Bacteria", "Archaea"]), "prok", np.where(dom == "Eukaryota", "euk", "other")))
    g = k[k["ko"].isin(kos)].groupby(["dom", "ko"])["count"].sum()
    out = {}
    for d in ("prok", "euk"):
        out[d] = g.get(d, pd.Series(dtype=float)).reindex(kos).fillna(0).to_numpy()
    return lib, out["prok"], out["euk"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--min-euk", type=float, default=5000)
    p.add_argument("--jobs", type=int, default=6)
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    # libraries -> samples, as build_ko_targets.py
    from build_ko_targets import CONTROL, library_samples
    from biogoship.tables import MANIFEST_DIR, REPO_ROOT

    head = open(REPO / "ko_counts_cov50_all.tsv").readline().rstrip("\n").split("\t")
    allk = [c.replace("count.ko.", "") for c in head]
    counts = pd.read_csv(REPO / "ko_counts_cov50_all.tsv", sep="\t", header=None, skiprows=1, index_col=0)
    counts.columns = allk
    clr = pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", index_col="sample_key")
    kos = list(clr.columns)
    ns = argparse.Namespace(workbook=REPO_ROOT / "FAIRe-NOAA_uci-ess-biogoship.xlsx", faire_samples=REPO / "data" / "faire" / "faire_samples.csv",
                            manifest_dir=MANIFEST_DIR)
    link = library_samples(counts.index, ns)
    info = link.join(pd.Series(counts.sum(axis=1), name="total"))
    keep = (info["samp_name"].notna() & ~info["samp_name"].fillna("").str.contains(CONTROL) & info["cruise"].notna()
            & info["sample_key"].isin(clr.index) & (info["total"] >= 1e5) & ~(info["depth_m"] > 10))
    info = info[keep]
    say(f"libraries: {len(info):,} for {info['sample_key'].nunique():,} of {len(clr):,} target samples")

    # check: rebuilding the CLR from these libraries reproduces ko_clr.csv
    prop_all = counts.loc[info.index].div(info["total"], axis=0)
    per_sample = prop_all.groupby(info["sample_key"]).mean()
    totals = info.groupby("sample_key")["total"].sum()
    pseudo = (0.5 / totals).reindex(per_sample.index).to_numpy()[:, None]
    lg = np.log(per_sample[kos].to_numpy() + pseudo)
    re_clr = lg - lg.mean(1, keepdims=True)
    diff = np.abs(re_clr - clr.loc[per_sample.index, kos].to_numpy()).max()
    say(f"check: rebuilt CLR vs ko_clr.csv max |diff| {diff:.2e}")

    # classes
    gc = pd.read_csv(REPO / "transfer" / "results" / "ko_gc" / "ko_gc.csv", index_col="ko")
    share = gc["genes_Eukaryota"] / gc[["genes_Eukaryota", "genes_Bacteria", "genes_Archaea"]].sum(1)
    cls_all = pd.Series(np.where(share <= 0.05, "prok", np.where(share >= 0.95, "euk", "shared")), index=share.index)
    names = pd.read_csv(REPO / "data" / "reference" / "kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "name"])
    names = names.assign(ko=names["ko"].str.replace("ko:", "")).set_index("ko")["name"]
    sym = names.reindex(kos).fillna("").str.split(";").str[0]
    markers = [k for k in kos if re.search(r"\b(rpl[A-Z]|rps[A-Z]|rpm[A-Z])\b", sym[k]) and "ribosomal protein" in names.get(k, "")]
    cls = cls_all.reindex(kos).fillna("unknown")
    cls[markers] = "marker"
    pd.DataFrame({"class": cls, "euk_gene_share": share.reindex(kos)}).rename_axis("ko").to_csv(OUT / "ko_classes.csv")
    say("classes over model KOs: " + ", ".join(f"{k} {v}" for k, v in cls.value_counts().items()))
    full_euk = [k for k in allk if cls_all.get(k) == "euk"]
    full_prok = [k for k in allk if cls_all.get(k) == "prok"]
    say(f"full table: eukaryote-specific {len(full_euk):,}, prokaryote-specific {len(full_prok):,}")

    # per-library domain split
    cache = OUT / "library_split.npz"
    libs = list(info.index)
    if cache.exists() and list(np.load(cache, allow_pickle=True)["libs"]) == libs:
        z = np.load(cache, allow_pickle=True)
        prok, euk = z["prok"], z["euk"]
    else:
        miss = [lib for lib in libs if not (TAXDIR / lib / "level_2_phylum_division_ko_counts.tsv.gz").exists()]
        if miss:
            raise SystemExit(f"{len(miss)} libraries lack level-2 files: {miss[:5]}")
        prok, euk = np.zeros((len(libs), len(kos))), np.zeros((len(libs), len(kos)))
        pos = {lib: i for i, lib in enumerate(libs)}
        with ProcessPoolExecutor(args.jobs) as ex:
            for j, (lib, pr, eu) in enumerate(ex.map(split_one, [(lib, kos) for lib in libs], chunksize=8)):
                prok[pos[lib]], euk[pos[lib]] = pr, eu
                if (j + 1) % 500 == 0:
                    print(f"  {j + 1:,}/{len(libs):,}", flush=True)
        np.savez_compressed(cache, libs=np.array(libs, dtype=object), kos=np.array(kos, dtype=object), prok=prok, euk=euk)
    tot_k = counts.loc[libs, kos].to_numpy()
    say(f"assigned share of model-KO reads: prok {prok.sum() / tot_k.sum():.3f}, euk {euk.sum() / tot_k.sum():.3f}, "
        f"other/unassigned {1 - (prok.sum() + euk.sum()) / tot_k.sum():.3f}")

    def sample_mean(m):
        return pd.DataFrame(m / info["total"].to_numpy()[:, None], index=libs).groupby(info["sample_key"].to_numpy()).mean()

    P = per_sample[kos]
    Pp, Pe = sample_mean(prok).loc[P.index].set_axis(kos, axis=1), sample_mean(euk).loc[P.index].set_axis(kos, axis=1)
    ps = pseudo  # 0.5 / sample total
    mk = [kos.index(k) for k in markers]
    logG = np.log(Pp.to_numpy()[:, mk] + ps).mean(1, keepdims=True)
    logG_all = np.log(P.to_numpy()[:, mk] + ps).mean(1, keepdims=True)
    say(f"prokaryote-assigned share of marker reads: median over samples {np.median(Pp.to_numpy()[:, mk].sum(1) / P.to_numpy()[:, mk].sum(1)):.3f}")
    Esum = prop_all[full_euk].sum(1).groupby(info["sample_key"]).mean().loc[P.index]
    Psum = prop_all[full_prok].sum(1).groupby(info["sample_key"]).mean().loc[P.index]
    Ereads = counts.loc[libs, full_euk].sum(1).groupby(info["sample_key"]).sum().loc[P.index]

    def write(name, df, note):
        df.round(5).rename_axis("sample_key").to_csv(OUT / f"{name}.csv")
        say(f"{name}: {df.shape[0]:,} samples x {df.shape[1]:,}  ({note})")

    A_k = [k for k in kos if cls[k] == "prok"]
    B_k = [k for k in kos if cls[k] == "euk"]
    S_k = [k for k in kos if cls[k] == "shared"]
    write("A", pd.DataFrame(np.log(P[A_k].to_numpy() + ps) - logG, index=P.index, columns=A_k), "log per prokaryote GE, prok-assigned markers")
    okB = (Ereads >= args.min_euk).to_numpy()
    say(f"B: samples with < {args.min_euk:,.0f} pooled eukaryote-specific reads dropped: {int((~okB).sum())}")
    write("B", pd.DataFrame((np.log(P[B_k].to_numpy() + ps) - np.log(Esum.to_numpy()[:, None]))[okB], index=P.index[okB], columns=B_k),
          "log share of pooled eukaryote-specific reads")
    # C (R68): direct read ratio, log(eukaryote-assigned / prokaryote-assigned KO reads) over all KOs, per library
    # (as shares of its KO reads), libraries of one sample averaged. From the per-library domain split of every KO
    # (scripts/ko_domain_reads.py cache). The R51 version (pooled reference-class KOs) is kept as C_refclass.csv.
    write("C_refclass", pd.DataFrame({"euk_prok_log_ratio": np.log(Esum / Psum)}, index=P.index),
          "R51 definition: log pooled euk-specific / prok-specific (classes from the reference-gene share)")
    zc = np.load(REPO / "data" / "results" / "ko_domain" / "library_domain.npz", allow_pickle=True)
    pos = {lib: i for i, lib in enumerate(zc["libs"])}
    if not set(libs) <= set(pos):
        raise SystemExit("library_domain.npz lacks some libraries; run scripts/ko_domain_reads.py first")
    ix = [pos[lib] for lib in libs]
    e_all = zc["euk"][ix].astype(np.float64).sum(1) / info["total"].to_numpy()
    p_all = zc["prok"][ix].astype(np.float64).sum(1) / info["total"].to_numpy()
    E_all = pd.Series(e_all, index=libs).groupby(info["sample_key"]).mean().loc[P.index]
    P_all = pd.Series(p_all, index=libs).groupby(info["sample_key"]).mean().loc[P.index]
    say(f"C: eukaryote-assigned share of KO reads, median over samples {E_all.median():.4f} (min {E_all.min():.5f}); "
        f"r with the R51 version {np.corrcoef(np.log(E_all / P_all), np.log(Esum / Psum))[0, 1]:.4f}")
    write("C", pd.DataFrame({"euk_prok_log_ratio": np.log(E_all / P_all)}, index=P.index), "log eukaryote / prokaryote KO reads (all KOs)")
    write("D", clr.loc[P.index, S_k], "whole-community CLR, unchanged")
    prev_p = (Pp[S_k] > 0).mean()
    prev_e = (Pe[S_k] > 0).mean()
    Dp_k = [k for k in S_k if prev_p[k] >= 0.9]
    De_k = [k for k in S_k if prev_e[k] >= 0.9]
    write("Dp", pd.DataFrame(np.log(Pp[Dp_k].to_numpy() + ps) - logG, index=P.index, columns=Dp_k), "prok-assigned part per prokaryote GE")
    write("De", pd.DataFrame((np.log(Pe[De_k].to_numpy() + ps) - np.log(Esum.to_numpy()[:, None]))[okB], index=P.index[okB], columns=De_k),
          "euk-assigned part within eukaryotes")
    # how the shared KOs split
    sp = Pp[S_k].sum() / P[S_k].sum()
    say(f"shared KOs: prokaryote-assigned share of reads, median over KOs {sp.median():.3f} (IQR {sp.quantile(.25):.3f}-{sp.quantile(.75):.3f}); "
        f"> 90 % prok {int((sp > 0.9).sum())}, < 10 % prok {int((sp < 0.1).sum())}")
    pd.DataFrame({"prok_share_reads": sp, "prev_prok": prev_p, "prev_euk": prev_e}).rename_axis("ko").round(4).to_csv(OUT / "shared_split.csv")
    # denominator: R49 correction size here
    eukfrac = Esum / (Esum + Psum)
    from scipy.stats import spearmanr
    say(f"denominator: all-54 minus prok-assigned log GE, median {np.median(logG_all - logG):.3f}; rho with euk fraction {spearmanr(logG_all[:, 0] - logG[:, 0], eukfrac)[0]:.2f}")
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
