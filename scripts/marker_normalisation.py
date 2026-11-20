#!/usr/bin/env python3
"""CLR vs a single-copy-marker denominator ("per genome equivalent"): does the CLR amplify KO dropouts? (RESULTS.md R46)

Concern: in the CLR, log(p_k) - mean over ALL KOs of log(p_j). When rare KOs drop to zero (pseudocount),
the mean falls and every detected KO's CLR rises, so a library that detects fewer KOs looks "common KOs up"
(the axis-1 / platform pattern of R45). With a fixed reference set M of universal single-copy genes,
  GE_k = log(p_k + e) - mean over m in M of log(p_m + e)          (log copies per genome equivalent)
rare-KO dropouts cannot move common KOs. Computed from the same rarefied proportions, GE_k = CLR_k - mean_M CLR_m.

Markers M: prokaryote-type ribosomal proteins of the KO set (KEGG names with an rpl* / rps* / rpm* gene
symbol; modifying enzymes rim* / prm* excluded) that are detected in every rarefied library. They are
universal single-copy in Bacteria and Archaea (the mitochondrial / plastid MRPL / MRPS copies share the KO;
eukaryotic cytosolic ribosomal proteins are separate "...e" KOs). GE is therefore per PROKARYOTE genome
equivalent; eukaryotic KOs rise with the eukaryote : prokaryote ratio rather than being normalised by it.
The eukaryote-type ribosomal proteins (RP-...e) are missing from 6-30 % of libraries and are reported, not
used. KOs that are not 100 % prevalent: a zero stays at the rarefied detection limit log(0.5 / N) minus the
marker mean; it is the same limit in every library (all rarefied to N), and it no longer shifts other KOs.

Libraries: KO table rows (--counts; ko_counts_cov50_all.tsv by default, or the cov0 table), blanks dropped,
rarefied to N = 250,000 KO reads (seed 0), proportions + 0.5 / N, over the 3,933 real-model KOs.

Per transform (CLR, GE), in units of the KO's within-cruise library SD under that transform:
  hiseq    HiSeq 4000 - NovaSeq 6000 within the mixed cruises (surface; as R45), per KO
  era      NovaSeq 6000 cruises - NovaSeq X cruises (cruise means of modelled surface libraries), per KO
  rerun    HLY2202 same DNA on two NovaSeq 6000 flowcells, mean paired difference per KO (control)
  summary  RMS over KOs; share |effect| > 0.5 SD; Spearman with KO mean abundance (dropout signature);
           cosine with the R45 F1 loading
  spread   per KO, within-cruise Spearman with the library KO-profile spread (SD of CLR; R45): median over
           KOs, and the median over the 500 most abundant KOs (what a CLR artefact would push up)
  markers  within-library SD of the markers' log proportions (how tight the reference is), per platform

Outputs: data/results/marker_norm/<label>/ ko_effects.csv, summary.csv, summary.txt, markers.csv.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/marker_normalisation.py                     # cov50
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/marker_normalisation.py --counts <cov0.tsv> --label cov0
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
LIBMETA = REPO / "transfer" / "results" / "library_metadata" / "library_metadata.csv"
MIXED = ["I09N-2016", "AE1319-2013", "BVAL46-2011", "NH1418-2014"]
N = 250_000


def platform(m) -> str:
    m = str(m)
    return "HiSeq 4000" if "HiSeq 4000" in m else ("NovaSeq X" if "NovaSeq X" in m else ("NovaSeq 6000" if "6000" in m else "other"))


def rarefied_logp(counts: Path, libs: set, kos: list[str], cache: Path) -> pd.DataFrame:
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        return pd.DataFrame(z["x"], index=z["index"], columns=z["columns"])
    rng = np.random.default_rng(0)
    header = open(counts).readline().rstrip("\n").split("\t")
    col = {re.sub(r"^count\.ko\.", "", c).replace("ko:", ""): i for i, c in enumerate(header)}
    missing = [k for k in kos if k not in col]
    if missing:
        print(f"{len(missing)} KOs not in {counts.name}; set to zero")
    rows, keys = [], []
    with open(counts) as fh:
        fh.readline()
        for line in fh:
            name, rest = line.split("\t", 1)
            if name not in libs:
                continue
            v = np.fromstring(rest, sep="\t")
            if v.sum() < N:
                continue
            r = rng.multinomial(N, v / v.sum())
            # the header has no row-name column: header position = value position
            p = np.array([r[col[k]] if k in col else 0 for k in kos]) / N
            rows.append(np.log(p + 0.5 / N))
            keys.append(name)
    x = np.array(rows, dtype=np.float32)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, x=x, index=np.array(keys, dtype=object), columns=np.array(kos, dtype=object))
    return pd.DataFrame(x, index=keys, columns=kos)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", type=Path, default=REPO / "ko_counts_cov50_all.tsv")
    p.add_argument("--label", default="cov50")
    args = p.parse_args()
    out = REPO / "data" / "results" / "marker_norm" / args.label
    out.mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    kos = [k.replace("ko:", "") for k in pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", nrows=0).columns.drop("sample_key")]
    lib = pd.read_csv(LIBMETA)
    lib = lib[lib["in_ko_table"] & lib["sample_id"].notna()].drop_duplicates("lib_id").set_index("lib_id")
    lib = lib[~lib.index.str.contains("BLANK|NTC|BLNK|CTRL|MOCK", case=False) & lib["cruise"].notna()]
    lib["platform"] = lib["sequencer_model"].map(platform)
    lp = rarefied_logp(args.counts, set(lib.index), kos, out / "rarefied_logp.npz")
    lib = lib.loc[lp.index]
    depth = pd.read_csv(REPO / "data" / "manifests" / "biogoship_samples.csv").drop_duplicates("sample_id").set_index("sample_id")["depth_m"]
    faire = pd.read_csv(REPO / "data" / "faire" / "faire_samples.csv").drop_duplicates("samp_name").set_index("samp_name")["depth_min_m"]
    lib["surface"] = (lib["sample_id"].map(depth).fillna(lib["samp_name"].map(faire)).fillna(0) <= 10).to_numpy()
    say(f"{args.label}: {len(lp):,} libraries rarefied to {N:,} KO reads; {len(kos):,} KOs")

    # markers
    names = pd.read_csv(REPO / "data" / "reference" / "kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "name"])
    names["ko"] = names["ko"].str.replace("ko:", "")
    names = names.set_index("ko")["name"]
    x = lp.to_numpy()
    floor = np.log(0.5 / N)
    prev = pd.Series((x > floor + 1e-6).mean(0), index=kos)
    sym = names.reindex(kos).fillna("").str.split(";").str[0]
    prok = [k for k in kos if re.search(r"\b(rpl[A-Z]|rps[A-Z]|rpm[A-Z])\b", sym[k]) and "ribosomal protein" in names.get(k, "")]
    euk = [k for k in kos if re.match(r"RP-[LS]\d+[A-Z]*e\b", sym[k])]
    markers = [k for k in prok if prev[k] == 1.0]
    mk = pd.DataFrame({"ko": prok + euk, "set": ["prokaryote"] * len(prok) + ["eukaryote (not used)"] * len(euk)})
    mk["name"] = mk["ko"].map(names)
    mk["prevalence"] = mk["ko"].map(prev)
    mk["used"] = mk["ko"].isin(markers)
    mk.round(4).to_csv(out / "markers.csv", index=False)
    say(f"markers: {len(markers)} prokaryote ribosomal-protein KOs at 100 % prevalence (of {len(prok)}); eukaryote-type "
        f"RP-...e KOs: {len(euk)}, prevalence {prev[euk].min():.2f}-{prev[euk].max():.2f} (not used)")
    mi = [kos.index(k) for k in markers]
    msd = pd.Series(x[:, mi].std(1), index=lp.index)
    say("within-library SD of marker log proportions (median by platform): "
        + ", ".join(f"{pf} {msd[lib['platform'] == pf].median():.3f}" for pf in ["HiSeq 4000", "NovaSeq 6000", "NovaSeq X"]))

    clr = x - x.mean(1, keepdims=True)
    ge = x - x[:, mi].mean(1, keepdims=True)
    spread = clr.std(1)
    cr = lib["cruise"].to_numpy()

    load = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "loadings.csv.gz", index_col=0)["F1"].reindex(kos).to_numpy()
    mean_ab = clr.mean(0)
    top = np.argsort(-mean_ab)[:500]
    rows, eff_cols = [], {}
    for tname, m in (("CLR", clr), ("GE (per prokaryote genome equivalent)", ge)):
        df = pd.DataFrame(m, index=lp.index, columns=kos)
        within = df - df.groupby(cr).transform("mean")
        sd = within.std().to_numpy()
        sd = np.where(sd > 0, sd, 1)
        # hiseq
        num, den = np.zeros(len(kos)), 0.0
        for c in MIXED:
            g = (cr == c) & lib["surface"].to_numpy()
            hi = g & (lib["platform"] == "HiSeq 4000").to_numpy()
            no = g & (lib["platform"] == "NovaSeq 6000").to_numpy()
            if hi.sum() < 2 or no.sum() < 2:
                continue
            wt = min(hi.sum(), no.sum())
            num += wt * (m[hi].mean(0) - m[no].mean(0))
            den += wt
        effs = {"hiseq": num / den / sd}
        # era
        mod = lib["in_real_targets"].to_numpy() & lib["surface"].to_numpy()
        cm = df[mod].groupby(cr[mod]).mean()
        seq = lib.loc[mod].groupby("cruise")["platform"].agg(lambda v: v.iloc[0] if v.nunique() == 1 else "mixed")
        effs["era"] = (cm[seq.reindex(cm.index) == "NovaSeq 6000"].mean() - cm[seq.reindex(cm.index) == "NovaSeq X"].mean()).to_numpy() / sd
        # rerun
        h = lib[lib["cruise"] == "HLY2202-2022"]
        a = h[h["flowcell_id"] == "HMY7WDSX7"].groupby("sample_id").apply(lambda q: q.index[0])
        b = h[h["flowcell_id"] == "HFJLFDSX7"].groupby("sample_id").apply(lambda q: q.index[0])
        both = a.index.intersection(b.index)
        effs["rerun"] = (df.loc[a[both]].to_numpy() - df.loc[b[both]].to_numpy()).mean(0) / sd
        # spread coupling
        sp_c = spread - pd.Series(spread).groupby(cr).transform("mean").to_numpy()
        wz = within.to_numpy()
        r_sp = np.array([spearmanr(wz[:, j], sp_c)[0] if np.ptp(wz[:, j]) > 0 else np.nan for j in range(len(kos))])
        for e, v in effs.items():
            eff_cols[f"{tname.split(' ')[0]}_{e}"] = v
            rows.append({"transform": tname, "contrast": e, "rms_sd": np.sqrt(np.nanmean(v ** 2)),
                         "share_abs_gt_0.5sd": np.nanmean(np.abs(v) > 0.5),
                         "rho_mean_abundance": spearmanr(v, mean_ab, nan_policy="omit")[0],
                         "cos_F1": float(np.nansum(v * load) / np.sqrt(np.nansum(v ** 2) * np.nansum(load ** 2)))})
        eff_cols[f"{tname.split(' ')[0]}_rho_spread"] = r_sp
        rows.append({"transform": tname, "contrast": "spread coupling (median rho, all KOs / top 500)",
                     "rms_sd": np.nanmedian(r_sp), "share_abs_gt_0.5sd": np.nanmedian(r_sp[top]),
                     "rho_mean_abundance": spearmanr(r_sp, mean_ab, nan_policy="omit")[0], "cos_F1": np.nan})
    summ = pd.DataFrame(rows)
    summ.round(4).to_csv(out / "summary.csv", index=False)
    pd.set_option("display.width", 200)
    say(summ.round(3).to_string(index=False))
    ke = pd.DataFrame(eff_cols, index=kos)
    ke["mean_clr"] = mean_ab
    ke["prevalence"] = prev
    ke["marker"] = ke.index.isin(markers)
    ke.round(5).to_csv(out / "ko_effects.csv", index_label="ko")
    (out / "summary.txt").write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
