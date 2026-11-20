#!/usr/bin/env python3
"""Library-level technical contrasts along axis 1 (RESULTS.md R39).

Uses transfer/results/library_metadata/library_metadata.csv (transfer/scripts/library_metadata.py).
Per library: KO CLR over the real-model KO set (data/targets/ko_clr.csv columns) from
ko_counts_cov50_all.tsv, after RAREFYING every library to N = 250,000 KO-assigned reads (one
multinomial draw over all 18,836 KO features, seed 0; libraries below N are left out), then
proportions, pseudocount 0.5 / N, CLR. Rarefying is needed: library depth differs several-fold
between runs (e.g. HLY2202's second run is a shallow top-up), and through the CLR pseudocount depth
moves the projection (axis 3, R33). Domain CLR from data/targets/taxa/domain_counts.tsv (not rarefied;
domain counts are large).

Axis-1 projection: library CLR, centred on the reference group of the contrast, divided by the KO's
SD (pooled within-cruise SD of the tile-level observed CLR of the real model), projected on the unit
F1 loading (transfer/results/residual_modules_noA). The same projection of the cruise mean CLRs gives
the scale of the era offset: mean of the NovaSeq 6000 cruises minus mean of the "...LT4" cruises
(modelled cruises only).

Contrasts (surface libraries, depth <= 10 m):
  a. HiSeq 4000 vs NovaSeq 6000, different samples of the same cruise (I09N-2016, AE1319-2013,
     BVAL46-2011, NH1418-2014): difference of means per cruise, pooled (weighted by the smaller n);
     permutation of the platform labels inside each cruise (9,999).
  b. The same DNA on two NovaSeq 6000 runs (HLY2202: HFJLFDSX7 vs HMY7WDSX7): paired differences.
  c. Two libraries of one sample on the same flowcell (I09N-2025): paired differences (floor).
  d. Two "...LT4" flowcells within A16N-2023 (22NCJYLT4 vs 23LT2GLT4), different samples.
Reported per contrast: axis-1 difference (and in units of the era offset), Bacteria and Eukaryota
CLR differences.

Output: transfer/results/library_metadata/sequencer_contrasts.csv, summary.txt.

Usage:
    ~/miniforge3/bin/python transfer/scripts/sequencer_contrasts.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "transfer" / "results" / "library_metadata"


def main() -> None:
    rng = np.random.default_rng(0)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    lib = pd.read_csv(OUT / "library_metadata.csv")
    lib = lib[lib["in_ko_table"] & lib["sample_id"].notna()].drop_duplicates("lib_id").set_index("lib_id")
    targets = pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", index_col="sample_key")
    kos = list(targets.columns)
    want = set(lib.index)
    # library rows of the KO table, KO-set columns only
    header = open(REPO / "ko_counts_cov50_all.tsv").readline().rstrip("\n").split("\t")
    names = [c.replace("count.ko.", "") for c in header]
    col_idx = {n: i for i, n in enumerate(names)}
    sel = [col_idx[k] for k in kos]
    N = 250_000
    rows, keys, small = [], [], 0
    with open(REPO / "ko_counts_cov50_all.tsv") as fh:
        fh.readline()
        for line in fh:
            name, rest = line.split("\t", 1)
            if name not in want:
                continue
            v = np.fromstring(rest, sep="\t")
            total = v.sum()
            if total < N:
                small += 1
                continue
            r = rng.multinomial(N, v / total)
            p = r[sel] / N
            lg = np.log(p + 0.5 / N)
            rows.append(lg - lg.mean())
            keys.append(name)
    clr = pd.DataFrame(rows, index=keys, columns=kos)
    say(f"Libraries rarefied to {N:,} KO-assigned reads: {len(clr):,}; below that, left out: {small:,}")
    lib = lib.loc[clr.index]
    depth = pd.read_csv(REPO / "data" / "manifests" / "biogoship_samples.csv").drop_duplicates("sample_id").set_index("sample_id")["depth_m"]
    faire = pd.read_csv(REPO / "data" / "faire" / "faire_samples.csv").drop_duplicates("samp_name").set_index("samp_name")["depth_min_m"]
    d = lib["sample_id"].map(depth).fillna(lib["samp_name"].map(faire))
    surface = (d.fillna(0) <= 10).to_numpy()
    dom = pd.read_csv(REPO / "data" / "targets" / "taxa" / "domain_counts.tsv", sep="\t", index_col=0)
    ld = np.log(dom.add(0.5))
    dclr = ld.sub(ld.mean(axis=1), axis=0).reindex(clr.index)

    load = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "loadings.csv.gz", index_col=0)
    w = load.loc[[k.replace("ko:", "") for k in kos], "F1"].to_numpy()
    w = w / np.linalg.norm(w)
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    within = targets.groupby(ks["cruise"].reindex(targets.index).to_numpy()).transform(lambda x: x - x.mean())
    sd = within.std().to_numpy()
    sd = np.where(sd > 0, sd, 1)

    def proj(m, ref):
        return ((m - ref) / sd) @ w

    # era scale: cruise means of the modelled libraries
    mod = lib["in_real_targets"].to_numpy() & surface
    cm = clr[mod].groupby(lib.loc[mod, "cruise"].to_numpy()).mean()
    seq = lib.loc[mod].groupby("cruise")["flowcell_id"].agg(lambda v: "LT4" if v.astype(str).str.endswith("LT4").all() else
                                                           ("6000" if v.astype(str).str.contains("DS|DR").all() else "mixed"))
    ref = cm.mean().to_numpy()
    cp = pd.Series(proj(cm.to_numpy(), ref), index=cm.index)
    era = cp[seq.reindex(cp.index) == "6000"].mean() - cp[seq.reindex(cp.index) == "LT4"].mean()
    say("Cruise means of library CLR projected on F1 (modelled surface libraries):")
    say("   " + ", ".join(f"{c} {v:+.2f} ({seq.get(c)})" for c, v in cp.sort_values().items()))
    say(f"Era offset (NovaSeq 6000 cruises - LT4 cruises) along F1: {era:+.3f}")

    out = []

    def report(name, diffs_axis, diffs_b, diffs_e, n, p=np.nan, kind=""):
        r = {"contrast": name, "kind": kind, "n": n, "axis1_diff": diffs_axis, "axis1_diff_in_era_units": diffs_axis / era,
             "bacteria_clr_diff": diffs_b, "eukaryota_clr_diff": diffs_e, "p": p}
        out.append(r)
        say(f"   {name:58s} n {n:4}  F1 {diffs_axis:+.3f} ({diffs_axis / era:+.2f} x era)  Bacteria {diffs_b:+.3f}  "
            f"Eukaryota {diffs_e:+.3f}" + (f"  p {p:.3f}" if np.isfinite(p) else ""))

    say("\n== a. HiSeq 4000 - NovaSeq 6000, different samples within a cruise (surface libraries)")
    pooled_num = pooled_den = 0.0
    strata = []
    for cruise in ["I09N-2016", "AE1319-2013", "BVAL46-2011", "NH1418-2014"]:
        m = (lib["cruise"] == cruise).to_numpy() & surface
        hi = m & lib["flowcell_id"].astype(str).str.endswith("BBXX").to_numpy()
        no = m & lib["flowcell_id"].astype(str).str.contains("DSX").to_numpy()
        if hi.sum() < 2 or no.sum() < 2:
            say(f"   {cruise}: HiSeq {hi.sum()}, NovaSeq {no.sum()} surface libraries: skipped")
            continue
        r0 = clr[no].mean().to_numpy()
        a = proj(clr[hi].to_numpy(), r0).mean() - proj(clr[no].to_numpy(), r0).mean()
        b = dclr.loc[hi, "Bacteria"].mean() - dclr.loc[no, "Bacteria"].mean()
        e = dclr.loc[hi, "Eukaryota"].mean() - dclr.loc[no, "Eukaryota"].mean()
        report(f"{cruise}: HiSeq ({hi.sum()}) - NovaSeq ({no.sum()})", a, b, e, int(hi.sum() + no.sum()), kind="different samples")
        wgt = min(hi.sum(), no.sum())
        pooled_num += wgt * a
        pooled_den += wgt
        strata.append((proj(clr[m].to_numpy(), r0), hi[m]))
    if strata:
        obs = pooled_num / pooled_den
        null = []
        for _ in range(9999):
            num = den = 0.0
            for s, h in strata:
                hp = rng.permutation(h)
                wgt = min(hp.sum(), (~hp).sum())
                num += wgt * (s[hp].mean() - s[~hp].mean())
                den += wgt
            null.append(num / den)
        pv = (1 + (np.abs(null) >= abs(obs)).sum()) / 10000
        report("pooled HiSeq - NovaSeq (weighted by the smaller group)", obs, np.nan, np.nan,
               int(sum(len(s) for s, _ in strata)), pv, "different samples")

    def paired(name, pairs, kind):
        if not pairs:
            say(f"   {name}: no pairs")
            return
        a_ = np.array([proj(clr.loc[[x]].to_numpy(), clr.loc[[y]].to_numpy()[0])[0] for x, y in pairs])
        b_ = np.array([dclr.loc[x, "Bacteria"] - dclr.loc[y, "Bacteria"] for x, y in pairs])
        e_ = np.array([dclr.loc[x, "Eukaryota"] - dclr.loc[y, "Eukaryota"] for x, y in pairs])
        report(name + " (mean)", a_.mean(), np.nanmean(b_), np.nanmean(e_), len(pairs), kind=kind)
        report(name + " (SD of the pair differences)", a_.std(ddof=1), np.nanstd(b_, ddof=1), np.nanstd(e_, ddof=1), len(pairs), kind=kind)

    say("\n== b. Same DNA, two NovaSeq 6000 runs (HLY2202: HMY7WDSX7 - HFJLFDSX7)")
    h = lib[(lib["cruise"] == "HLY2202-2022") & surface]
    pairs = []
    for s, g in h.groupby("sample_id"):
        x = g[g["flowcell_id"] == "HMY7WDSX7"].index
        y = g[g["flowcell_id"] == "HFJLFDSX7"].index
        if len(x) and len(y):
            pairs.append((x[0], y[0]))
    paired("HLY2202 HMY7WDSX7 - HFJLFDSX7, same samples", pairs, "same DNA, two runs")

    say("\n== c. Two libraries of one sample on the same flowcell (I09N-2025)")
    g9 = lib[(lib["cruise"] == "I09N-2025") & surface]
    pairs = [tuple(g.index[:2]) for _, g in g9.groupby("sample_id") if len(g) >= 2]
    paired("I09N-2025 library 2 - library 1, same flowcell", pairs, "same sample, same run")

    say("\n== d. Two LT4 flowcells within A16N-2023 (different samples; Protocol A excluded)")
    a16 = lib[(lib["cruise"] == "A16N-2023") & surface & lib["in_real_targets"]]
    x, y = a16[a16["flowcell_id"] == "23LT2GLT4"].index, a16[a16["flowcell_id"] == "22NCJYLT4"].index
    r0 = clr.loc[y].mean().to_numpy()
    report(f"A16N 23LT2GLT4 ({len(x)}) - 22NCJYLT4 ({len(y)})", proj(clr.loc[x].to_numpy(), r0).mean() - 0.0,
           dclr.loc[x, "Bacteria"].mean() - dclr.loc[y, "Bacteria"].mean(),
           dclr.loc[x, "Eukaryota"].mean() - dclr.loc[y, "Eukaryota"].mean(), len(x) + len(y), kind="different samples")
    pd.DataFrame(out).round(4).to_csv(OUT / "sequencer_contrasts.csv", index=False)
    (OUT / "sequencer_contrasts_summary.txt").write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
