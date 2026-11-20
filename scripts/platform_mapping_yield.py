#!/usr/bin/env python3
"""Is the platform offset on axis 1 a mapping / annotation effect? (RESULTS.md R45)

Mechanism to test: reads that map to the reference (and to KO-annotated genes) at different rates on
different platforms would shift the KO composition, especially of eukaryotic and divergent genes,
without any sequencing-chemistry bias.

Pipeline facts (biogoship-metagenomics/mapping_pipeline, read-only): minimap2 -ax sr against the
reference, no identity / MAPQ filter; reads per reference gene from samtools coverage (numreads: mapped
reads, both mates, secondary / QC-fail / duplicate excluded); genes kept if their breadth of coverage
passes --coverage_threshold (50 % for ko_counts_cov50_all.tsv); KO counts from the eggNOG table. No
per-read identity is kept, so the KOs' hit identity is not available.

Per library (transfer/results/library_metadata/library_metadata.csv; blanks and controls dropped):
  read length         nominal (FAIRe seq_run_id) and ENA base_count / read_count where deposited
  total_reads         FAIRe input_read_count (read pairs; equals ENA spots where both exist)
  map_rate            reads with a domain assignment (sum of data/targets/taxa/domain_counts.tsv) / total_reads;
                      domain counts include both mates, total_reads are pairs, so this is mapped reads per pair (0-2)
  ko_yield            KO-assigned reads (row sum of ko_counts_cov50_all.tsv) / total_reads
  ko_per_mapped       KO-assigned / domain-assigned reads
  euk_frac            Eukaryota / (Archaea + Bacteria + Eukaryota) domain reads
  f1                  axis-1 projection of the library's KO CLR, rarefied to 250,000 KO reads (as R39 / R44;
                      cache transfer/results/ko_gc/rarefied_clr.npz), per-KO SD units, unit F1 loading
Tests:
  1  platform medians (all libraries) and cruise-level medians
  2  HiSeq 4000 - NovaSeq 6000 within the four mixed cruises (surface): pooled difference (weighted by the
     smaller group), permutation of platform within cruise (9,999)
  3  HLY2202, same DNA on two NovaSeq 6000 runs: paired differences (control)
  4  within cruise (cruise-centred, all libraries): Spearman of f1 and euk_frac with map_rate, ko_yield,
     ko_per_mapped, log10 total_reads
  5  mediation: the HiSeq - NovaSeq f1 difference before and after removing a within-cruise linear fit of f1
     on map_rate, ko_per_mapped and log10 total_reads (fit over all libraries, cruise-centred; then the
     contrast of 2 on the residual)
  7  KO-profile spread per library (after rarefying): clr_sd (SD of the CLR over KOs), undetected (share of
     KOs with zero rarefied count), Shannon diversity of the KO proportions; within-cruise association with f1,
     the HiSeq contrast, and as mediators in 5; cruise-level Spearman with f1
  6  per KO: platform effects (transfer/results/ko_gc/gc_effects.csv) vs the KO's eukaryote share of genes,
     its mean abundance (mean CLR over targets) and gene length

Outputs: data/results/platform_mapping/library_yield.csv, tests.csv, summary.txt;
figures/output/platform_mapping_yield.png.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/platform_mapping_yield.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "results" / "platform_mapping"
FIG = REPO / "figures" / "output"
LIBMETA = REPO / "transfer" / "results" / "library_metadata" / "library_metadata.csv"
MIXED = ["I09N-2016", "AE1319-2013", "BVAL46-2011", "NH1418-2014"]
METRICS = ["map_rate", "ko_yield", "ko_per_mapped", "euk_frac", "log_reads", "clr_sd", "undetected", "shannon", "f1"]


def platform(m: str) -> str:
    m = str(m)
    return "HiSeq 4000" if "HiSeq 4000" in m else ("NovaSeq X" if "NovaSeq X" in m else ("NovaSeq 6000" if "6000" in m else "other"))


def main() -> None:
    rng = np.random.default_rng(0)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    OUT.mkdir(parents=True, exist_ok=True)
    lib = pd.read_csv(LIBMETA)
    lib = lib[lib["in_ko_table"] & lib["sample_id"].notna()].drop_duplicates("lib_id").set_index("lib_id")
    lib = lib[~lib.index.str.contains("BLANK|NTC|BLNK|CTRL|MOCK", case=False)]
    lib["platform"] = lib["sequencer_model"].map(platform)
    dom = pd.read_csv(REPO / "data" / "targets" / "taxa" / "domain_counts.tsv", sep="\t", index_col=0)
    lib = lib.join(dom, how="inner")
    lib["domain_reads"] = dom.reindex(lib.index).sum(1)
    lib["map_rate"] = lib["domain_reads"] / lib["total_reads"]
    lib["ko_yield"] = lib["ko_assigned_reads"] / lib["total_reads"]
    lib["ko_per_mapped"] = lib["ko_assigned_reads"] / lib["domain_reads"]
    lib["euk_frac"] = lib["Eukaryota"] / lib[["Archaea", "Bacteria", "Eukaryota"]].sum(1)
    lib["log_reads"] = np.log10(lib["total_reads"])
    ena = pd.read_csv(REPO / "transfer" / "data" / "ena_PRJNA656268_runs.tsv", sep="\t").set_index("run_accession")
    lib["ena_bases_per_spot"] = lib["sra_run"].map(ena["base_count"] / ena["read_count"])
    lib = lib[np.isfinite(lib["map_rate"]) & (lib["total_reads"] > 1e5)]

    # f1 projection of the rarefied library CLR
    z = np.load(REPO / "transfer" / "results" / "ko_gc" / "rarefied_clr.npz", allow_pickle=True)
    clr = pd.DataFrame(z["x"], index=z["index"], columns=z["columns"])
    targets = pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", index_col="sample_key")
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    within = targets.groupby(ks["cruise"].reindex(targets.index).to_numpy()).transform(lambda x: x - x.mean())
    sd = within.std().to_numpy()
    sd = np.where(sd > 0, sd, 1)
    load = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "loadings.csv.gz", index_col=0)
    w = load.loc[[k.replace("ko:", "") for k in targets.columns], "F1"].to_numpy()
    w = w / np.linalg.norm(w)
    clr = clr[list(targets.columns)]
    f1 = pd.Series(((clr.to_numpy() - clr.to_numpy().mean(0)) / sd) @ w, index=clr.index)
    lib["f1"] = f1.reindex(lib.index)
    x = clr.to_numpy()
    pr = np.exp(x)
    pr /= pr.sum(1, keepdims=True)
    lib["clr_sd"] = pd.Series(x.std(1), index=clr.index).reindex(lib.index)
    lib["undetected"] = pd.Series(np.isclose(x, x.min(1, keepdims=True), atol=1e-4).mean(1), index=clr.index).reindex(lib.index)
    lib["shannon"] = pd.Series(-(pr * np.log(pr)).sum(1), index=clr.index).reindex(lib.index)
    depth = pd.read_csv(REPO / "data" / "manifests" / "biogoship_samples.csv").drop_duplicates("sample_id").set_index("sample_id")["depth_m"]
    faire = pd.read_csv(REPO / "data" / "faire" / "faire_samples.csv").drop_duplicates("samp_name").set_index("samp_name")["depth_min_m"]
    lib["surface"] = lib["sample_id"].map(depth).fillna(lib["samp_name"].map(faire)).fillna(0) <= 10
    keep = ["sample_id", "cruise", "platform", "flowcell_id", "read_length", "ena_bases_per_spot", "surface", "total_reads",
            "domain_reads", "ko_assigned_reads", "map_rate", "ko_yield", "ko_per_mapped", "euk_frac", "log_reads", "clr_sd", "undetected", "shannon", "f1"]
    lib[keep].round(5).to_csv(OUT / "library_yield.csv", index_label="lib_id")
    rows = []

    # 1. read length and platform medians
    say("== 1. Read length and yield by platform (all libraries in the KO table, blanks dropped)")
    say("   nominal read length: " + "; ".join(f"{p}: {', '.join(sorted(set(g['read_length'].astype(str).str.split(' ').str[0])))}"
                                              for p, g in lib.groupby("platform")))
    say("   ENA bases per spot: " + "; ".join(f"{p}: median {g['ena_bases_per_spot'].median():.1f} (n {g['ena_bases_per_spot'].notna().sum()})"
                                             for p, g in lib.groupby("platform") if g["ena_bases_per_spot"].notna().any()))
    med = lib.groupby("platform")[["total_reads", "map_rate", "ko_yield", "ko_per_mapped", "euk_frac"]].median()
    med.insert(0, "libraries", lib.groupby("platform").size())
    say(med.round(3).to_string())
    cmed = lib.groupby("cruise").agg(platform=("platform", lambda v: "/".join(sorted(set(v)))), n=("f1", "size"),
                                     map_rate=("map_rate", "median"), ko_yield=("ko_yield", "median"),
                                     ko_per_mapped=("ko_per_mapped", "median"), euk_frac=("euk_frac", "median"),
                                     reads_M=("total_reads", lambda v: v.median() / 1e6), f1=("f1", "mean"))
    say("\n   cruise medians (f1: mean)")
    say(cmed.sort_values("f1").round(3).to_string())
    cmed.round(4).to_csv(OUT / "cruise_yield.csv")
    one = cmed[~cmed["platform"].str.contains("/")]
    for m in ["map_rate", "ko_yield", "ko_per_mapped", "euk_frac", "reads_M"]:
        r = spearmanr(one["f1"], one[m])[0]
        rows.append({"test": "cruise level (single-platform cruises)", "metric": m, "value": r, "n": len(one)})
    say("   cruise-level Spearman with mean f1 (single-platform cruises, n %d): " % len(one)
        + ", ".join(f"{m} {r['value']:+.2f}" for m, r in zip(["map_rate", "ko_yield", "ko_per_mapped", "euk_frac", "reads_M"], rows[-5:])))

    # 2. HiSeq - NovaSeq within mixed cruises
    def contrast(df, col):
        num = den = 0.0
        for c in MIXED:
            g = df[(df["cruise"] == c) & df["surface"]]
            hi, no = g[g["platform"] == "HiSeq 4000"][col].dropna(), g[g["platform"] == "NovaSeq 6000"][col].dropna()
            if len(hi) < 2 or len(no) < 2:
                continue
            wt = min(len(hi), len(no))
            num += wt * (hi.mean() - no.mean())
            den += wt
        return num / den

    def perm_p(df, col, obs, n=9999):
        g = df[df["cruise"].isin(MIXED) & df["surface"] & df["platform"].isin(["HiSeq 4000", "NovaSeq 6000"])].copy()
        cnt = 0
        for _ in range(n):
            g["platform"] = g.groupby("cruise")["platform"].transform(lambda v: rng.permutation(v.to_numpy()))
            cnt += abs(contrast(g, col)) >= abs(obs)
        return (cnt + 1) / (n + 1)

    say("\n== 2. HiSeq 4000 - NovaSeq 6000 within the mixed cruises (surface libraries)")
    for c in MIXED:
        g = lib[(lib["cruise"] == c) & lib["surface"]]
        say(f"   {c}: " + "; ".join(f"{p} n {len(h)}: map {h['map_rate'].median():.3f}, KO yield {h['ko_yield'].median():.3f}, "
                                    f"KO/mapped {h['ko_per_mapped'].median():.3f}, euk {h['euk_frac'].median():.3f}, "
                                    f"reads {h['total_reads'].median() / 1e6:.1f} M"
                                    for p, h in g.groupby("platform")))
    for m in METRICS:
        obs = contrast(lib, m)
        p = perm_p(lib, m, obs, 1999)
        rows.append({"test": "HiSeq - NovaSeq 6000 within cruise", "metric": m, "value": obs, "p": p})
        say(f"   pooled {m:14s} {obs:+.4f}  (permutation p {p:.4f})")

    # 3. HLY2202 rerun
    say("\n== 3. HLY2202 same DNA, HMY7WDSX7 - HFJLFDSX7 (control)")
    h = lib[lib["cruise"] == "HLY2202-2022"]
    a = h[h["flowcell_id"] == "HMY7WDSX7"].groupby("sample_id")[METRICS].mean()
    b = h[h["flowcell_id"] == "HFJLFDSX7"].groupby("sample_id")[METRICS].mean()
    both = a.index.intersection(b.index)
    d = a.loc[both] - b.loc[both]
    for m in METRICS:
        rows.append({"test": "HLY2202 rerun (paired)", "metric": m, "value": d[m].mean(), "sd": d[m].std(), "n": len(both)})
        say(f"   {m:14s} mean {d[m].mean():+.4f} (SD {d[m].std():.4f}, n {len(both)})")

    # 4. within-cruise associations
    say("\n== 4. Within cruise (cruise-centred; all libraries with f1): Spearman")
    cc = lib.dropna(subset=["f1", "cruise"]).copy()
    cc = cc[np.isfinite(cc[METRICS]).all(1)]
    for m in METRICS:
        cc[m + "_c"] = cc[m] - cc.groupby("cruise")[m].transform("mean")
    for y in ["f1", "euk_frac", "clr_sd", "undetected"]:
        parts = []
        for x in ["map_rate", "ko_yield", "ko_per_mapped", "log_reads", "clr_sd", "undetected", "shannon"] + (["euk_frac"] if y == "f1" else []):
            r = spearmanr(cc[y + "_c"], cc[x + "_c"])[0]
            rows.append({"test": f"within cruise: {y} vs", "metric": x, "value": r, "n": len(cc)})
            parts.append(f"{x} {r:+.3f}")
        say(f"   {y}: " + ", ".join(parts) + f"  (n {len(cc):,})")

    # 5. mediation
    say("\n== 5. HiSeq - NovaSeq f1 difference with mapping metrics removed (within-cruise linear fit)")
    X = cc[["map_rate_c", "ko_per_mapped_c", "log_reads_c"]].to_numpy()
    beta = np.linalg.lstsq(X, cc["f1_c"].to_numpy(), rcond=None)[0]
    cc["f1_resid"] = cc["f1"] - X @ beta
    r2 = 1 - ((cc["f1_c"] - X @ beta) ** 2).sum() / (cc["f1_c"] ** 2).sum()
    say(f"   fit: f1 ~ map_rate + ko_per_mapped + log10 reads (cruise-centred), R2 {r2:.3f}; "
        f"coefficients {beta[0]:+.2f}, {beta[1]:+.2f}, {beta[2]:+.2f}")
    for sub, cols in (("map_rate only", ["map_rate_c"]), ("ko_per_mapped only", ["ko_per_mapped_c"]), ("log reads only", ["log_reads_c"]),
                      ("all three", ["map_rate_c", "ko_per_mapped_c", "log_reads_c"]),
                      ("KO-profile spread (clr_sd)", ["clr_sd_c"]), ("Shannon", ["shannon_c"]), ("undetected KOs", ["undetected_c"]),
                      ("spread + euk_frac", ["clr_sd_c", "euk_frac_c"]), ("euk_frac only", ["euk_frac_c"])):
        Xs = cc[cols].to_numpy()
        bs = np.linalg.lstsq(Xs, cc["f1_c"].to_numpy(), rcond=None)[0]
        cc["tmp"] = cc["f1"] - Xs @ bs
        before, after = contrast(cc, "f1"), contrast(cc, "tmp")
        r2s = 1 - ((cc["f1_c"] - Xs @ bs) ** 2).sum() / (cc["f1_c"] ** 2).sum()
        rows.append({"test": "mediation: HiSeq - NovaSeq f1", "metric": sub, "value": after, "before": before, "within_cruise_r2": r2s})
        say(f"   {sub:28s} within-cruise R2 {r2s:.3f}; f1 difference {before:+.3f} -> {after:+.3f} ({1 - after / before:+.0%} explained)")

    # 6. per KO
    say("\n== 6. Per KO: platform effects vs KO properties (Spearman)")
    ge = pd.read_csv(REPO / "transfer" / "results" / "ko_gc" / "gc_effects.csv").set_index("ko")
    gc = pd.read_csv(REPO / "transfer" / "results" / "ko_gc" / "ko_gc.csv").set_index("ko")
    ge["euk_share"] = gc["genes_Eukaryota"] / gc["genes"]
    raw = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv", nrows=0)  # noqa: F841
    ge["mean_clr"] = targets.mean().rename(lambda k: k.replace("ko:", ""))
    for e in ["hiseq", "era", "rerun", "F1"]:
        parts = []
        for x in ["euk_share", "mean_clr", "length_mean"]:
            r = spearmanr(ge[e], ge[x], nan_policy="omit")[0]
            rows.append({"test": f"per KO: {e} vs", "metric": x, "value": r, "n": int(ge[[e, x]].notna().all(1).sum())})
            parts.append(f"{x} {r:+.3f}")
        say(f"   {e:6s}: " + ", ".join(parts))

    say("\n== 7. Cruise level (all cruises): Spearman of mean f1 with mean spread / Shannon / undetected")
    cm = cc.groupby("cruise")[["f1", "clr_sd", "shannon", "undetected", "euk_frac"]].mean()
    for m in ["clr_sd", "shannon", "undetected", "euk_frac"]:
        r = spearmanr(cm["f1"], cm[m])[0]
        rows.append({"test": "cruise level (all cruises): f1 vs", "metric": m, "value": r, "n": len(cm)})
        say(f"   {m:12s} {r:+.2f} (n {len(cm)})")
    pd.DataFrame(rows).round(5).to_csv(OUT / "tests.csv", index=False)
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")
    figure(lib, cc)


def figure(lib, cc):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6})
    col = {"HiSeq 4000": "#E69F00", "NovaSeq 6000": "#0072B2", "NovaSeq X": "#009E73"}
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5))
    g = lib[lib["cruise"].isin(MIXED) & lib["surface"]]
    for ax, m, lab in ((axes[0], "map_rate", "Mapped reads per read pair"),
                       (axes[1], "f1", "Axis-1 score of the library (F1)")):
        for i, c in enumerate(MIXED):
            for j, p in enumerate(["HiSeq 4000", "NovaSeq 6000"]):
                v = g[(g["cruise"] == c) & (g["platform"] == p)][m].dropna()
                x = i + (j - 0.5) * 0.36 + np.random.default_rng(i * 2 + j).uniform(-0.08, 0.08, len(v))
                ax.scatter(x, v, s=7, color=col[p], alpha=0.6, linewidths=0, label=p if i == 0 else None)
                if len(v):
                    ax.plot([i + (j - 0.5) * 0.36 - 0.13, i + (j - 0.5) * 0.36 + 0.13], [v.median()] * 2, color="#222222", lw=1.2)
        ax.set_xticks(range(len(MIXED)))
        ax.set_xticklabels([c.replace("-", "\n") for c in MIXED], fontsize=7)
        ax.set_ylabel(lab)
    axes[0].legend(frameon=False, fontsize=6.5, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, handletextpad=0.2)
    ax = axes[2]
    for p, h in cc.groupby("platform"):
        if p in col:
            ax.scatter(h["clr_sd_c"], h["f1_c"], s=3, color=col[p], alpha=0.35, linewidths=0, label=p, rasterized=True)
    ax.axhline(0, color="#999999", lw=0.5)
    ax.axvline(0, color="#999999", lw=0.5)
    ax.set_xlabel("KO-profile spread (SD of CLR),\ncentred within cruise")
    ax.set_ylabel("Axis-1 score, centred within cruise")
    ax.legend(frameon=False, fontsize=6.5, markerscale=2.5, loc="lower right", handletextpad=0.2)
    fig.tight_layout()
    fig.savefig(FIG / "platform_mapping_yield.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
