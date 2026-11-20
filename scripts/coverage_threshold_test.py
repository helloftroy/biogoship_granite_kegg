#!/usr/bin/env python3
"""Does the 50 % breadth-of-coverage filter turn library complexity into missing KOs? cov0 vs cov50 (RESULTS.md R47)

Hypothesis (R45 / R46): a library with fewer unique fragments covers fewer rare genes to 50 % breadth; under
the cov50 filter their reads are dropped, the KOs go to zero, the KO profile becomes less even, and (through
the CLR) axis 1 rises. HiSeq 4000 libraries of the four mixed cruises show exactly this. If the HiSeq -
NovaSeq gap in unevenness and axis 1 mostly disappears with no coverage filter (cov0), the filter is the
converter; duplicate rates would then only explain why some libraries have lower complexity.

Per coverage table (--tables label=path ...; default cov50 = ko_counts_cov50_all.tsv and
cov0 = ko_counts_cov0_all.tsv.gz if present), libraries of the KO table (blanks dropped), rarefied to
250,000 KO reads (seed 0; scripts/marker_normalisation.rarefied_logp, cached in
data/results/marker_norm/<label>/), over the 3,933 real-model KOs:
  clr_sd       SD of the CLR over KOs (KO-profile spread, R45)
  undetected   share of KOs with zero rarefied count
  shannon      Shannon diversity of the KO proportions
  f1_clr       projection of the CLR on the unit F1 loading (transfer/results/residual_modules_noA), per-KO
               SD units (pooled within-cruise SD of the real-model targets; the same scaling for both tables)
  f1_ge        the same for log copies per prokaryote genome equivalent (54 ribosomal-protein markers, R46)
  map_rate     domain-assigned reads / read pairs, if a domain table exists (the abundance counts also pass
               the coverage filter; cov0 gives all mapped reads)
Contrasts on the libraries present in both tables:
  HiSeq 4000 - NovaSeq 6000 within the four mixed cruises (surface): pooled, weighted by the smaller group;
  permutation within cruise (1,999); also in units of the metric's within-cruise SD (all libraries), so
  cov0 and cov50 compare on their own scales; per-cruise differences.
  HLY2202 same DNA on two NovaSeq 6000 flowcells: mean and SD of paired differences (control).
  Per library: cov0 - cov50 change of each metric, by platform (mixed cruises).

Outputs: data/results/coverage_threshold/library_metrics.csv, contrasts.csv, summary.txt;
figures/output/coverage_threshold_gap.png.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/coverage_threshold_test.py
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/coverage_threshold_test.py --tables cov50=ko_counts_cov50_all.tsv cov0=ko_counts_cov0_all.tsv.gz
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from marker_normalisation import LIBMETA, MIXED, N, platform, rarefied_logp  # noqa: E402

OUT = REPO / "data" / "results" / "coverage_threshold"
FIG = REPO / "figures" / "output"
METRICS = ["clr_sd", "undetected", "shannon", "f1_clr", "f1_ge", "map_rate"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tables", nargs="+", default=None, help="label=path pairs")
    args = p.parse_args()
    tables = dict(t.split("=", 1) for t in args.tables) if args.tables else {"cov50": "ko_counts_cov50_all.tsv", "cov0": "ko_counts_cov0_all.tsv.gz"}
    tables = {k: (Path(v) if Path(v).is_absolute() else REPO / v) for k, v in tables.items()}
    missing = [k for k, v in tables.items() if not v.exists()]
    for k in missing:
        print(f"{k}: {tables[k]} not found; skipped")
        tables.pop(k)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    kos = [k.replace("ko:", "") for k in pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", nrows=0).columns.drop("sample_key")]
    targets = pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", index_col="sample_key")
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    tw = targets.groupby(ks["cruise"].reindex(targets.index).to_numpy()).transform(lambda x: x - x.mean())
    sd = tw.std().to_numpy()
    sd = np.where(sd > 0, sd, 1)
    w = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "loadings.csv.gz", index_col=0)["F1"].reindex(kos).to_numpy()
    w = w / np.linalg.norm(w)
    names = pd.read_csv(REPO / "data" / "reference" / "kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "name"])
    names = names.assign(ko=names["ko"].str.replace("ko:", "")).set_index("ko")["name"]
    sym = names.reindex(kos).fillna("").str.split(";").str[0]
    prok = [i for i, k in enumerate(kos) if re.search(r"\b(rpl[A-Z]|rps[A-Z]|rpm[A-Z])\b", sym[k]) and "ribosomal protein" in names.get(k, "")]

    lib = pd.read_csv(LIBMETA)
    lib = lib[lib["in_ko_table"] & lib["sample_id"].notna()].drop_duplicates("lib_id").set_index("lib_id")
    lib = lib[~lib.index.str.contains("BLANK|NTC|BLNK|CTRL|MOCK", case=False) & lib["cruise"].notna()]
    lib["platform"] = lib["sequencer_model"].map(platform)
    depth = pd.read_csv(REPO / "data" / "manifests" / "biogoship_samples.csv").drop_duplicates("sample_id").set_index("sample_id")["depth_m"]
    faire = pd.read_csv(REPO / "data" / "faire" / "faire_samples.csv").drop_duplicates("samp_name").set_index("samp_name")["depth_min_m"]
    lib["surface"] = (lib["sample_id"].map(depth).fillna(lib["samp_name"].map(faire)).fillna(0) <= 10).to_numpy()

    frames = []
    for label, path in tables.items():
        lp = rarefied_logp(path, set(lib.index), kos, REPO / "data" / "results" / "marker_norm" / label / "rarefied_logp.npz")
        x = lp.to_numpy().astype(float)
        floor = np.log(0.5 / N)
        clr = x - x.mean(1, keepdims=True)
        ge = x - x[:, prok].mean(1, keepdims=True)
        pr = np.exp(x)
        pr /= pr.sum(1, keepdims=True)
        m = pd.DataFrame({"clr_sd": clr.std(1), "undetected": (x <= floor + 1e-6).mean(1), "shannon": -(pr * np.log(pr)).sum(1),
                          "f1_clr": (clr / sd) @ w, "f1_ge": (ge / sd) @ w}, index=lp.index)
        dpath = path.with_name(path.name.replace(".tsv.gz", "").replace(".tsv", "") + "_domain_counts.tsv")
        if label == "cov50":
            dpath = REPO / "data" / "targets" / "taxa" / "domain_counts.tsv"
        if dpath.exists():
            dom = pd.read_csv(dpath, sep="\t", index_col=0)
            m["map_rate"] = dom.sum(1).reindex(m.index) / lib["total_reads"].reindex(m.index)
        else:
            m["map_rate"] = np.nan
        m["coverage"] = label
        frames.append(m)
        say(f"{label}: {len(m):,} libraries rarefied ({path.name})")
    met = pd.concat(frames)
    common = set.intersection(*[set(f.index) for f in frames])
    met = met[met.index.isin(common)]
    met = met.join(lib[["cruise", "platform", "flowcell_id", "sample_id", "surface", "total_reads"]])
    met.round(5).to_csv(OUT / "library_metrics.csv", index_label="lib_id")
    say(f"libraries in all tables: {len(common):,}")

    def gap(df, col):
        num = den = 0.0
        per = {}
        for c in MIXED:
            g = df[(df["cruise"] == c) & df["surface"]]
            hi, no = g[g["platform"] == "HiSeq 4000"][col].dropna(), g[g["platform"] == "NovaSeq 6000"][col].dropna()
            if len(hi) < 2 or len(no) < 2:
                continue
            per[c] = hi.mean() - no.mean()
            wt = min(len(hi), len(no))
            num += wt * per[c]
            den += wt
        return (num / den if den else np.nan), per

    rows = []
    for label in tables:
        d = met[met["coverage"] == label]
        say(f"\n== {label}: HiSeq 4000 - NovaSeq 6000 within the mixed cruises (surface)")
        for col in METRICS:
            if d[col].isna().all():
                continue
            obs, per = gap(d, col)
            wsd = (d[col] - d.groupby("cruise")[col].transform("mean")).std()
            g = d[d["cruise"].isin(MIXED) & d["surface"] & d["platform"].isin(["HiSeq 4000", "NovaSeq 6000"])].copy()
            cnt = 0
            for _ in range(1999):
                g["platform"] = g.groupby("cruise")["platform"].transform(lambda v: rng.permutation(v.to_numpy()))
                cnt += abs(gap(g, col)[0]) >= abs(obs)
            pv = (cnt + 1) / 2000
            h = d[d["cruise"] == "HLY2202-2022"]
            a = h[h["flowcell_id"] == "HMY7WDSX7"].groupby("sample_id")[col].mean()
            b = h[h["flowcell_id"] == "HFJLFDSX7"].groupby("sample_id")[col].mean()
            rr = (a - b).dropna()
            rows.append({"coverage": label, "metric": col, "hiseq_gap": obs, "gap_in_within_cruise_sd": obs / wsd, "p": pv,
                         "within_cruise_sd": wsd, "rerun_mean": rr.mean(), "rerun_sd": rr.std(), "rerun_n": len(rr),
                         **{f"gap_{c}": v for c, v in per.items()}})
            say(f"   {col:11s} gap {obs:+9.4f} = {obs / wsd:+.2f} within-cruise SD (p {pv:.4f}); per cruise "
                + ", ".join(f"{c.split('-')[0]} {v:+.3f}" for c, v in per.items())
                + f"; rerun {rr.mean():+.4f} (SD {rr.std():.4f}, n {len(rr)})")
    con = pd.DataFrame(rows)
    con.round(5).to_csv(OUT / "contrasts.csv", index=False)
    if len(tables) >= 2 and "cov0" in tables and "cov50" in tables:
        say("\n== cov0 vs cov50: gap (in within-cruise SD units) and share of the cov50 gap left at cov0")
        for col in METRICS:
            r50 = con[(con["coverage"] == "cov50") & (con["metric"] == col)]
            r0 = con[(con["coverage"] == "cov0") & (con["metric"] == col)]
            if r50.empty or r0.empty:
                continue
            g50, g0 = r50["gap_in_within_cruise_sd"].iloc[0], r0["gap_in_within_cruise_sd"].iloc[0]
            say(f"   {col:11s} cov50 {g50:+.2f} SD -> cov0 {g0:+.2f} SD ({g0 / g50:.0%} left); raw {r50['hiseq_gap'].iloc[0]:+.4f} -> {r0['hiseq_gap'].iloc[0]:+.4f}")
        wide = met.rename_axis("lib_id").reset_index().pivot_table(index="lib_id", columns="coverage", values=METRICS, dropna=False)
        ch = pd.DataFrame({c: wide[(c, "cov0")] - wide[(c, "cov50")] for c in METRICS if (c, "cov0") in wide and (c, "cov50") in wide}).join(lib[["cruise", "platform", "surface"]])
        mx = ch[ch["cruise"].isin(MIXED) & ch["surface"]]
        say("\n   per-library change cov0 - cov50 (mixed cruises, surface; mean by platform)")
        say(mx.groupby("platform")[[c for c in METRICS if c in mx]].mean().round(4).to_string())
        figure(met)
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")


def figure(met):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6})
    col = {"HiSeq 4000": "#E69F00", "NovaSeq 6000": "#0072B2"}
    panels = [("undetected", "KOs undetected (share)"), ("clr_sd", "KO-profile spread (SD of CLR)"), ("f1_clr", "Axis-1 score (CLR)")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    g = met[met["cruise"].isin(MIXED) & met["surface"] & met["platform"].isin(col)]
    for ax, (m, lab) in zip(axes, panels):
        for i, cov in enumerate(["cov50", "cov0"]):
            for j, pf in enumerate(col):
                v = g[(g["coverage"] == cov) & (g["platform"] == pf)][m]
                # centre within cruise so cruises pool
                cen = g[g["coverage"] == cov].groupby("cruise")[m].transform("mean")
                v = (v - cen.reindex(v.index)).dropna()
                x = i + (j - 0.5) * 0.36 + np.random.default_rng(i * 2 + j).uniform(-0.08, 0.08, len(v))
                ax.scatter(x, v, s=6, color=col[pf], alpha=0.5, linewidths=0, label=pf if (i == 0 and m == "undetected") else None)
                ax.plot([i + (j - 0.5) * 0.36 - 0.13, i + (j - 0.5) * 0.36 + 0.13], [v.mean()] * 2, color="#222222", lw=1.2)
        ax.axhline(0, color="#999999", lw=0.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Coverage ≥ 50%", "No coverage filter"])
        ax.set_ylabel(lab + ",\ncentred within cruise", fontsize=7)
    axes[0].legend(frameon=False, fontsize=6.5, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, handletextpad=0.2)
    fig.tight_layout()
    fig.savefig(FIG / "coverage_threshold_gap.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
