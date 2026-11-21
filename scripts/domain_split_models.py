#!/usr/bin/env python3
"""Satellite models and variance budgets for the domain-split targets (RESULTS.md R51).

Targets (scripts/domain_split_targets.py, data/targets/domain_split/): A prokaryote-specific KOs per
prokaryote genome equivalent, B eukaryote-specific KOs within eukaryotes, C log eukaryote / prokaryote,
D shared KOs whole-community CLR, Dp / De the prokaryote- / eukaryote-assigned parts of the shared KOs.
Reference: the current CLR (data/targets/ko_clr.csv) on all 3,933 KOs (the R43 16 %) and on the same KOs and
samples as each new target, so a change is the transform, not the KO set.

Same samples / tiles as R43 (tune_figure_kos.load: tiles merged by content, targets averaged per tile), one
cruise held out, ridge with per-KO alpha (GCV on the training rows). Feature sets (all fitted on training
cruises only):
  temp      in-situ temperature, spline (not a satellite input; reference)
  sst       satellite SST of the tile, spline
  olci      16 OLCI band means of the tile, splines
  raw_sat   OLCI + SST, splines (R17 "raw satellite")
  granite   32 Granite PCs
  baseline  32 Granite PCs + in-situ temperature spline (the R43 model)
For clr_all, A, D and Dp also: baseline + C observed (spline; an oracle, C is not observed from space) and baseline + C
predicted (the out-of-fold baseline prediction of C, a legitimate two-stage model).

Scores per target set and feature set: explained share = 1 - SS(residual) / SS(total), pooled over tiles and
KOs, raw units and with every KO standardised; median per-KO R2; per KO, which feature set scores best.
Variance budget (baseline), sequential as R43: explained; platform (NovaSeq 6000 vs X cruise groups);
other cruise offsets; sequencing depth (spline of log10 KO reads, per-cruise centred); first two
within-cruise residual PCs (standardised PCA; their raw share); replicate noise (same-tile sample pairs,
capped); remaining.

Outputs: data/results/domain_split/{scores,per_ko,budget,best_feature}.csv, summary.txt;
figures/output/domain_split_skill.png, domain_split_budget.png.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/domain_split_models.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from tune_figure_kos import build, load, ridge_gcv  # noqa: E402

OUT = REPO / "data" / "results" / "domain_split"
FIG = REPO / "figures" / "output"
TD = REPO / "data" / "targets" / "domain_split"
EMB = [REPO / "data" / "embeddings_sst_l4", REPO / "data" / "faire" / "embeddings_sst_l4"]
NOVASEQ_6000 = ["I09N", "P18", "I07N", "AMT28", "C13.5/A13.5", "HLY2202-2022"]
FSETS = {"temp": ["T~"], "sst": ["SST~"], "olci": ["B~"], "raw_sat": ["S~"], "granite": ["P32"], "baseline": ["P32", "T~"]}
LABEL = {"clr_all": "Current CLR, all 3,933 KOs", "A": "A  prokaryote-specific, per genome equivalent",
         "clr_A": "   same KOs, current CLR", "B": "B  eukaryote-specific, within eukaryotes", "clr_B": "   same KOs, current CLR",
         "C": "C  eukaryote / prokaryote amount", "D": "D  shared KOs, whole-community CLR",
         "Dp": "Dp shared, prokaryote part per GE", "clr_Dp": "   same KOs, current CLR",
         "De": "De shared, eukaryote part within eukaryotes", "clr_De": "   same KOs, current CLR"}  # fmt: skip


def cv_predict(blocks, emb, temp, sat, y, cr, extra=None):
    pred = np.full_like(y, np.nan)
    for c in np.unique(cr):
        te, tr = np.where(cr == c)[0], np.where(cr != c)[0]
        cache = {}
        xtr, xte = build(blocks, tr, te, emb, temp, sat, cache)
        if extra is not None:
            from sklearn.preprocessing import SplineTransformer, StandardScaler

            s = StandardScaler().fit(extra[tr])
            sp = SplineTransformer(n_knots=8, degree=3, extrapolation="linear").fit(s.transform(extra[tr]))
            xtr, xte = np.hstack([xtr, sp.transform(s.transform(extra[tr]))]), np.hstack([xte, sp.transform(s.transform(extra[te]))])
        pred[te] = np.asarray(ridge_gcv(xtr, y[tr], xte)).reshape(len(te), -1)
    return pred


def scores(y, pred):
    yc = y - y.mean(0)
    r = y - pred
    sd = yc.std(0)
    sd = np.where(sd > 0, sd, 1)
    per_ko = 1 - (r**2).sum(0) / np.maximum((yc**2).sum(0), 1e-12)
    return {"explained_raw": 1 - (r**2).sum() / (yc**2).sum(), "explained_std": 1 - ((r / sd) ** 2).sum() / ((yc / sd) ** 2).sum(),
            "median_ko_r2": float(np.median(per_ko)), "share_ko_r2_gt_0.2": float((per_ko > 0.2).mean())}, per_ko  # fmt: skip


def centre(m, g):
    return m - pd.DataFrame(m).groupby(g).transform("mean").to_numpy()


def budget(y, pred, cr, unit, rows, sample_targets, scale):
    from sklearn.preprocessing import SplineTransformer

    yc = y - y.mean(0)
    res = y - pred
    res = res - res.mean(0)
    ss = lambda m: ((m * scale) ** 2).sum()  # noqa: E731
    tot = ss(yc)
    b = {"explained": tot - ss(res)}
    plat = np.isin(cr, NOVASEQ_6000)
    pm = pd.DataFrame(res).groupby(plat).transform("mean").to_numpy()
    cm = pd.DataFrame(res).groupby(cr).transform("mean").to_numpy()
    b["platform"] = ss(pm)
    b["cruise offsets"] = ss(cm - pm)
    w = res - cm
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    logtot = rows["sample_id"].map(np.log10(ks["total"])).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()
    basis = centre(SplineTransformer(n_knots=8, degree=3).fit_transform(logtot[:, None]), cr)
    beta, *_ = np.linalg.lstsq(basis, w, rcond=None)
    d = basis @ beta
    b["sequencing depth"] = ss(d)
    w2 = w - d
    if y.shape[1] > 2:
        z = w2 * scale
        u, s, vt = np.linalg.svd(z / np.where(z.std(0) > 0, z.std(0), 1), full_matrices=False)
        for i in range(2):
            x = u[:, [i]]
            f = x @ (x.T @ w2)
            b[f"within-cruise PC{i + 1}"] = ss(f)
            w2 = w2 - f
    multi = rows.groupby("tile_key")["sample_id"].agg(list)
    diffs = []
    for samples in multi[multi.str.len() >= 2]:
        v = sample_targets.reindex([x for x in samples if x in sample_targets.index]).to_numpy()
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                diffs.append(v[i] - v[j])
    sigma2 = (np.array(diffs) ** 2).mean(0) / 2
    k = multi.reindex(unit.index).str.len().fillna(1).to_numpy()
    noise = (sigma2 * scale**2).sum() * (1 / k).sum()
    left = ss(w2)
    b["replicate noise"] = min(noise, left)
    b["remaining"] = left - min(noise, left)
    return pd.Series(b) / tot


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    clr = pd.read_csv(REPO / "data" / "targets" / "ko_clr.csv", index_col="sample_key")
    sets = {k: pd.read_csv(TD / f"{k}.csv", index_col="sample_key") for k in ["A", "B", "C", "D", "Dp", "De"]}
    ref = {"clr_all": clr, "clr_A": clr[sets["A"].columns], "clr_B": clr.loc[sets["B"].index, sets["B"].columns],
           "clr_Dp": clr[sets["Dp"].columns], "clr_De": clr.loc[sets["De"].index, sets["De"].columns]}  # fmt: skip
    order = ["C", "clr_all", "A", "clr_A", "B", "clr_B", "D", "Dp", "clr_Dp", "De", "clr_De"]
    allsets = {**sets, **ref}
    tmp = OUT / "_tmp_targets.csv"
    rows_sc, rows_ko, rows_b, best = [], [], [], []
    c_oof = None
    for name in order:
        t = allsets[name]
        t.rename_axis("sample_key").to_csv(tmp)
        args = argparse.Namespace(targets=tmp, embeddings_dirs=EMB)
        emb, temp, sat, y, cr, unit, rows = load(args, list(t.columns), return_rows=True)
        say(f"\n{name}: {LABEL[name].strip()} - {y.shape[0]:,} tiles x {y.shape[1]:,}")
        per = {}
        preds = {}
        for fs, blocks in FSETS.items():
            preds[fs] = cv_predict(blocks, emb, temp, sat, y, cr)
            sc, per[fs] = scores(y, preds[fs])
            rows_sc.append({"target": name, "features": fs, "tiles": y.shape[0], "kos": y.shape[1], **sc})
            say(f"   {fs:9s} explained raw {sc['explained_raw']:6.1%}  std {sc['explained_std']:6.1%}  median KO R2 {sc['median_ko_r2']:+.3f}  "
                f"KOs R2>0.2 {sc['share_ko_r2_gt_0.2']:.0%}")
        if name == "C":
            c_oof = pd.Series(preds["baseline"][:, 0], index=unit.index)
            c_tile = pd.Series(y[:, 0], index=unit.index)
            c_rows = rows[["sample_id", "tile_key"]]
        if name in ("clr_all", "A", "D", "Dp"):
            # C as covariate: observed (oracle) and out-of-fold predicted, matched by sample through the tiles
            csamp = sets["C"]["euk_prok_log_ratio"]
            c_obs = rows["sample_id"].map(csamp).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()[:, None]
            c_by_sample = c_rows.assign(p=c_rows["tile_key"].map(c_oof)).drop_duplicates("sample_id").set_index("sample_id")["p"]
            c_pred = rows["sample_id"].map(c_by_sample).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()[:, None]
            for fs, extra in (("baseline+C_observed", c_obs), ("baseline+C_predicted", c_pred)):
                ok = ~np.isnan(extra[:, 0])
                p = cv_predict(FSETS["baseline"], emb[ok], temp[ok], sat[ok], y[ok], cr[ok], extra=extra[ok])
                sc, per[fs] = scores(y[ok], p)
                rows_sc.append({"target": name, "features": fs, "tiles": int(ok.sum()), "kos": y.shape[1], **sc})
                say(f"   {fs:21s} explained raw {sc['explained_raw']:6.1%}  std {sc['explained_std']:6.1%}  median KO R2 {sc['median_ko_r2']:+.3f}  "
                    f"(tiles {ok.sum():,})")
        pk = pd.DataFrame(per, index=t.columns)
        pk.insert(0, "target", name)
        rows_ko.append(pk)
        if y.shape[1] > 1:
            core = pk[list(FSETS)]
            wins = core.idxmax(1).where(core.max(1) > 0.05, "none > 0.05")
            vc = wins.value_counts(normalize=True)
            best.append(pd.DataFrame({"target": name, "best": vc.index, "share": vc.to_numpy()}))
            say("   best feature set per KO (R2 > 0.05): " + ", ".join(f"{k} {v:.0%}" for k, v in vc.items()))
        sd = (y - y.mean(0)).std(0)
        for scale_name, s in (("raw", np.ones(y.shape[1])), ("standardised", 1 / np.where(sd > 0, sd, 1))):
            b = budget(y, preds["baseline"], cr, unit, rows, t, s)
            rows_b.append(pd.DataFrame({"target": name, "scale": scale_name, "component": b.index, "share": b.to_numpy()}))
            say(f"   budget ({scale_name}): " + "; ".join(f"{k} {v:.1%}" for k, v in b.items()))
    tmp.unlink()
    sc = pd.DataFrame(rows_sc)
    sc.round(4).to_csv(OUT / "scores.csv", index=False)
    pd.concat(rows_ko).round(4).rename_axis("ko").to_csv(OUT / "per_ko.csv")
    bd = pd.concat(rows_b)
    bd.round(4).to_csv(OUT / "budget.csv", index=False)
    pd.concat(best).round(4).to_csv(OUT / "best_feature.csv", index=False)
    say("\nExplained share (raw units), target x features:")
    say((sc.pivot(index="target", columns="features", values="explained_raw").reindex(order) * 100).round(1).to_string())
    say("\nExplained share (KOs standardised):")
    say((sc.pivot(index="target", columns="features", values="explained_std").reindex(order) * 100).round(1).to_string())
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")
    figures(sc, bd, order)


def figures(sc, bd, order):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    main_sets = ["clr_all", "A", "B", "C", "D", "Dp", "De"]
    names = {"clr_all": "Current CLR\n(all KOs)", "A": "A\nProk-specific\nper GE", "B": "B\nEuk-specific\nwithin euks",
             "C": "C\nEuk / prok\namount", "D": "D\nShared\nCLR", "Dp": "Dp\nShared, prok\nper GE", "De": "De\nShared, euk\nwithin euks"}  # fmt: skip
    fs = ["temp", "sst", "olci", "raw_sat", "granite", "baseline"]
    flab = {"temp": "In-situ temperature", "sst": "Satellite SST", "olci": "OLCI bands", "raw_sat": "OLCI + SST",
            "granite": "Granite PCs", "baseline": "Granite PCs + temperature"}
    cols = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#000000"]
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    piv = sc.pivot(index="target", columns="features", values="explained_raw")
    w = 0.13
    for j, f in enumerate(fs):
        v = piv.reindex(main_sets)[f].to_numpy() * 100
        ax.bar(np.arange(len(main_sets)) + (j - 2.5) * w, v, w * 0.92, color=cols[j], label=flab[f])
    ax.axhline(0, color="#666666", lw=0.5)
    ax.set_xticks(range(len(main_sets)))
    ax.set_xticklabels([names[s] for s in main_sets], fontsize=7)
    ax.set_ylabel("Variance explained, cruise held out (%)")
    ax.legend(frameon=False, fontsize=6.5, ncol=3, loc="upper left", bbox_to_anchor=(0, 1.18))
    lo = min(-5, np.nanmin(piv.reindex(main_sets)[fs].to_numpy()) * 100 - 2)
    ax.set_ylim(lo, None)
    fig.tight_layout()
    fig.savefig(FIG / "domain_split_skill.png", dpi=300)
    plt.close(fig)

    comp = ["explained", "platform", "cruise offsets", "sequencing depth", "within-cruise PC1", "within-cruise PC2", "replicate noise", "remaining"]
    ccol = {"explained": "#0072B2", "platform": "#E69F00", "cruise offsets": "#BBBBBB", "sequencing depth": "#F0C674",
            "within-cruise PC1": "#009E73", "within-cruise PC2": "#8FD1B5", "replicate noise": "#555555", "remaining": "#DDDDDD"}
    b = bd[bd["scale"] == "raw"].pivot(index="target", columns="component", values="share").reindex(main_sets).reindex(columns=comp).fillna(0) * 100
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    left = np.zeros(len(main_sets))
    for c in comp:
        ax.barh(range(len(main_sets)), b[c], left=left, height=0.62, color=ccol[c], label=c[0].upper() + c[1:], edgecolor="white", linewidth=0.5)
        for i, v in enumerate(b[c]):
            if v >= 4:
                ax.text(left[i] + v / 2, i, f"{v:.0f}", ha="center", va="center", fontsize=6,
                        color="white" if c in ("explained", "replicate noise") else "black")
        left += b[c].to_numpy()
    ax.set_yticks(range(len(main_sets)))
    ax.set_yticklabels([names[s].replace("\n", " ") for s in main_sets], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of target variance (%), baseline model, cumulative")
    ax.legend(frameon=False, fontsize=6, ncol=4, loc="upper left", bbox_to_anchor=(0, 1.2))
    fig.tight_layout()
    fig.savefig(FIG / "domain_split_budget.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
