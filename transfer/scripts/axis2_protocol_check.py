#!/usr/bin/env python3
"""Is residual axis 2 biology or protocol? And are the C13.5 v1/v2 pairs replicates? (RESULTS.md R36)

Uses the real-model targets, without the A16N Protocol A samples (data/targets/ko_clr.csv
built with data/manifests/excluded_samples.csv), and the residual PCA refitted on them
(scripts/residual_modules.py --out-dir transfer/results/residual_modules_noA).

1. Refit vs R33: cosine of each refitted component's loadings with the R33 component of the
   same name, on the 3,921 KOs the two target sets share.
2. Protocol-group PCA: the within-cruise residuals PCA'd separately in
     old workflow  enzymatic extraction, lysis buffer, Nextera tagmentation, HiSeq/NovaSeq 6000:
                   I09N, P18, I07N, AMT28, C13.5/A13.5, HLY2202-2022 (HLY: DNA Flex, NovaSeq 6000)
     new workflow  bead beating, DNA/RNA Shield (GOMECC4: beads, no preservative), NovaSeq X:
                   A16N-2023, I05-2023, A13.5-2024, I08S-2024, GOMECC4
   (P02-2022 and I09N-2025, workflow not verified, in neither). For W1-W3 of the full fit:
   the best-matching component in each group (rank, share, |cosine| with the full-fit loading)
   and the congruence between the two groups' matches. Reference distribution: the same
   statistic for 200 random splits of the same 11 cruises into 6 + 5. If the protocol split
   gives the same congruence as random splits, the axis does not depend on the workflow.
3. Paired samples (same time and place, both in one tile): treated minus partner, from the
   targets that keep every sample (data/targets/ko_clr_incl_protocolA.csv; common KOs):
     A16N Protocol A - Protocol B                (technical contrast, cruise report)
     C13.5/A13.5 v2 - v1 (or unsuffixed)         (duplicates, no protocol difference recorded)
   Per pair the difference of standardised CLRs (KO SD of the within-cruise residuals),
   projected on the axis loadings (in SD of the tile scores); domain CLR differences. KO-level
   systematic effect: RMS over KOs of the mean paired difference, against its sign-flip null
   (whole pairs flipped, 999 draws); share of KOs with |t| > 3; cosine with the loadings.

Outputs: transfer/results/axis2_protocol/ (refit_vs_r33.csv, group_pca.csv, pairs_tests.csv,
pairs_ko_effect.csv, summary.txt).

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/axis2_protocol_check.py
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from added_fields_all_kos import cv  # noqa: E402
from residual_modules import standardise, top_pca  # noqa: E402
from tune_figure_kos import load  # noqa: E402

OLD = ["I09N", "P18", "I07N", "AMT28", "C13.5/A13.5", "HLY2202-2022"]
NEW = ["A16N-2023", "I05-2023", "A13.5-2024", "I08S-2024", "GOMECC4"]
OUT = REPO / "transfer" / "results" / "axis2_protocol"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--targets", type=Path, default=REPO / "data" / "targets" / "ko_clr.csv")
    p.add_argument("--targets-all", type=Path, default=REPO / "data" / "targets" / "ko_clr_incl_protocolA.csv")
    p.add_argument("--embeddings-dirs", type=Path, nargs="+",
                   default=[REPO / "data" / "embeddings_sst_l4", REPO / "data" / "faire" / "embeddings_sst_l4"])
    p.add_argument("--refit", type=Path, default=REPO / "transfer" / "results" / "residual_modules_noA")
    p.add_argument("--r33", type=Path, default=REPO / "data" / "results" / "residual_modules")
    p.add_argument("--random-splits", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    new_l = pd.read_csv(args.refit / "loadings.csv.gz", index_col=0)
    old_l = pd.read_csv(args.r33 / "loadings.csv.gz", index_col=0)
    common = new_l.index.intersection(old_l.index)

    # 1. refit vs R33
    rows = []
    for c in [f"F{i}" for i in range(1, 6)] + [f"W{i}" for i in range(1, 6)]:
        a, b = new_l.loc[common, c], old_l.loc[common, c]
        rows.append({"component": c, "cosine_with_r33": float(a @ b / np.linalg.norm(a) / np.linalg.norm(b))})
    refit = pd.DataFrame(rows)
    comp = pd.read_csv(args.refit / "components.csv").set_index("component")
    comp33 = pd.read_csv(args.r33 / "components.csv").set_index("component")
    refit = refit.join(comp[["share", "split_half_congruence"]], on="component").join(
        comp33[["share", "split_half_congruence"]].add_suffix("_r33"), on="component")
    refit.round(4).to_csv(OUT / "refit_vs_r33.csv", index=False)
    say(f"== 1. Refit without Protocol A vs R33 ({len(common):,} common KOs)")
    say(refit.round(3).to_string(index=False))

    # 2. protocol-group PCA
    kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, _, y, cruises, unit, rows_ = load(args, kos, return_rows=True)
    _, pred = cv(emb, temp[:, 0], None, y, cruises, np.arange(len(unit)))
    res = y - pred
    within = res - pd.DataFrame(res).groupby(cruises).transform("mean").to_numpy()
    ko_ids = [k.replace("count.ko.", "").replace("ko:", "") for k in kos]
    say(f"\n== 2. Protocol-group PCA of the within-cruise residuals ({len(unit):,} tiles, {len(kos):,} KOs)")

    def group_vecs(names):
        m = np.isin(cruises, names)
        ev, vec = top_pca(standardise(within[m]))
        return ev[:30] / ev.sum(), vec[:, :30]

    def match(vecs, target):
        cos = np.abs(vecs.T @ target) / np.linalg.norm(target)
        j = int(cos.argmax())
        return j, cos[j]

    targets = {c: new_l.loc[ko_ids, c].to_numpy() for c in ["W1", "W2", "W3"]}
    so, vo = group_vecs(OLD)
    sn, vn = group_vecs(NEW)
    grp_rows = []
    pool = OLD + NEW
    rand = {c: [] for c in targets}
    for _ in range(args.random_splits):
        perm = list(rng.permutation(pool))
        _, va = group_vecs(perm[:6])
        _, vb = group_vecs(perm[6:])
        for c, t in targets.items():
            ja, _ = match(va, t)
            jb, _ = match(vb, t)
            rand[c].append(abs(va[:, ja] @ vb[:, jb]))
    for c, t in targets.items():
        jo, co = match(vo, t)
        jn, cn = match(vn, t)
        cong = abs(vo[:, jo] @ vn[:, jn])
        r = np.array(rand[c])
        grp_rows.append({"component": c, "old_rank": jo + 1, "old_share": so[jo], "old_cos_full": co, "new_rank": jn + 1,
                         "new_share": sn[jn], "new_cos_full": cn, "old_vs_new_congruence": cong,
                         "random_split_median": np.median(r), "random_split_p05": np.percentile(r, 5),
                         "share_random_below": (r < cong).mean()})
        say(f"   {c}: old workflow rank {jo + 1} ({so[jo]:.1%}, cos with full {co:.2f}); new workflow rank {jn + 1} "
            f"({sn[jn]:.1%}, cos {cn:.2f}); old vs new congruence {cong:.2f} | random 6+5 splits median {np.median(r):.2f} "
            f"(5th pct {np.percentile(r, 5):.2f}; {(r < cong).mean():.0%} of splits lower)")
    pd.DataFrame(grp_rows).round(4).to_csv(OUT / "group_pca.csv", index=False)

    # 3. paired samples
    allc = pd.read_csv(args.targets_all, index_col="sample_key")
    kc = [k for k in kos if k in allc.columns]
    idx = [kos.index(k) for k in kc]
    sd = within.std(0)[idx]
    dom = pd.read_csv(REPO / "data" / "targets" / "taxa" / "taxa_domain_clr.csv").set_index("sample_key")
    wl = {c: new_l.loc[[ko_ids[i] for i in idx], c].to_numpy() for c in ["F1", "F2", "F3", "W1", "W2", "W3"]}
    wl = {c: v / np.linalg.norm(v) for c, v in wl.items()}
    tsd = {c: ((standardise(within) if c[0] == "W" else standardise(res))[:, idx] @ wl[c]).std() for c in wl}
    sdf = res.std(0)[idx]

    pairs = []
    excl = pd.read_csv(REPO / "data" / "manifests" / "excluded_samples.csv")
    import warnings

    warnings.filterwarnings("ignore")
    sm = pd.read_excel(REPO / "FAIRe-NOAA_uci-ess-biogoship.xlsx", sheet_name="sampleMetadata", header=2)
    sm = sm[sm["samp_name"].astype(str).str.endswith("_MG")]
    for a in excl["sample_key"]:
        ev = sm.loc[sm["samp_name"] == a, "eventDate"]
        if ev.empty:
            continue
        part = sm[(sm["eventDate"] == ev.iloc[0]) & (sm["samp_name"] != a) & (sm["short_name"] == "A16N-2023")]["samp_name"]
        for b in part:
            pairs.append(("A16N Protocol A - B", a, b))
    man = pd.read_csv(REPO / "data" / "manifests" / "biogoship_samples.csv")
    c13 = man[man["section"] == "C13.5/A13.5"]
    for _, g in c13.groupby(["collection_time_utc", "lat", "lon"]):
        ids = sorted(g["sample_id"])
        v2 = [s for s in ids if s.endswith("v2")]
        v1 = [s for s in ids if not s.endswith("v2")]
        if v2 and v1:
            pairs.append(("C13.5 v2 - v1", v2[0], v1[0]))
    pairs = pd.DataFrame(pairs, columns=["contrast", "treated", "partner"])
    pairs = pairs[pairs["treated"].isin(allc.index) & pairs["partner"].isin(allc.index)]
    say(f"\n== 3. Paired samples with both libraries in the targets: {pairs['contrast'].value_counts().to_dict()}")

    from scipy.stats import wilcoxon

    test_rows, eff_rows = [], []
    for contrast, g in pairs.groupby("contrast"):
        d = allc.loc[g["treated"], kc].to_numpy() - allc.loc[g["partner"], kc].to_numpy()
        feats = {}
        for c in wl:
            feats[c] = (d / (sd if c[0] == "W" else sdf)) @ wl[c] / tsd[c]
        for c in ["Archaea", "Bacteria", "Eukaryota"]:
            feats[f"{c} CLR"] = dom.loc[g["treated"], c].to_numpy() - dom.loc[g["partner"], c].to_numpy()
        say(f"   {contrast} ({len(g)} pairs):")
        for c, v in feats.items():
            v = v[np.isfinite(v)]
            pv = wilcoxon(v).pvalue if len(v) >= 5 else np.nan
            test_rows.append({"contrast": contrast, "variable": c, "n": len(v), "mean_diff": v.mean(), "sd_diff": v.std(ddof=1),
                              "same_sign": int((np.sign(v) == np.sign(v.mean())).sum()), "wilcoxon_p": pv})
            say(f"      {c:14s} mean {v.mean():+.3f} (SD {v.std(ddof=1):.3f}; {test_rows[-1]['same_sign']}/{len(v)} same sign; p {pv:.3f})")
        z = d / sd
        eff = z.mean(0)
        t = eff / (z.std(0, ddof=1) / np.sqrt(len(z)))
        rms = np.sqrt((eff ** 2).mean())
        null = np.array([np.sqrt(((z * rng.choice([-1, 1], len(z))[:, None]).mean(0) ** 2).mean()) for _ in range(999)])
        cos = {c: float(eff @ wl[c] / np.linalg.norm(eff)) for c in wl}
        eff_rows.append({"contrast": contrast, "pairs": len(g), "rms_effect": rms, "rms_null_median": np.median(null),
                         "p_rms": (1 + (null >= rms).sum()) / 1000, "share_abs_t_gt_3": (np.abs(t) > 3).mean(),
                         **{f"cos {c}": v for c, v in cos.items()}})
        say(f"      KO-level: RMS mean effect {rms:.3f} SD vs sign-flip null median {np.median(null):.3f} (p {eff_rows[-1]['p_rms']:.3f}); "
            f"|t| > 3: {(np.abs(t) > 3).mean():.1%}; cosine " + ", ".join(f"{c} {v:+.2f}" for c, v in cos.items()))
    pd.DataFrame(test_rows).round(4).to_csv(OUT / "pairs_tests.csv", index=False)
    pd.DataFrame(eff_rows).round(4).to_csv(OUT / "pairs_ko_effect.csv", index=False)
    pairs.to_csv(OUT / "pairs_list.csv", index=False)
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
