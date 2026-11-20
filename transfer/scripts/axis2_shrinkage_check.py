#!/usr/bin/env python3
"""Is residual axis 2 ridge shrinkage (attenuated predictions) rather than a missing driver? (RESULTS.md R38)

Real model (data/targets/ko_clr.csv, A16N Protocol A excluded; baseline = 32 Granite PCs +
temperature spline, ridge, one cruise held out). W2 = second component of the within-cruise
residual PCA (as R33/R36; sign and identity checked against transfer/results/residual_modules_noA).

If the predictions are compressed (calibration slope b > 1), residual ~ (b - 1) * prediction + noise,
so a residual axis would correlate POSITIVELY with the predicted signal, and rescaling the
predictions would remove it. Tests:

1. Calibration slopes per KO: observed on predicted, pooled over all tiles and within cruise
   (both centred per cruise). Summaries for the 50 KOs loading most positively and most
   negatively on W2, the KEGG photosynthesis KOs (BRITE C 'Photosynthesis', 'Photosynthesis -
   antenna proteins'), and all KOs; and against KOs of the same R2 decile.
2. Axis-2 composite: mean z-score (within cruise) of the 50 + KOs minus the 50 - KOs, from the
   observed CLRs and from the predictions (z with the observed SDs). Pearson r of the W2 score
   with each, within cruise; R2 of the W2 score on the predicted composite (= the share that
   compression could account for). The same for Cyanobacteriota and Cryptophyta phylum CLRs,
   predicted by the same baseline.
3. Recalibration: per KO, the within-cruise slope estimated on the 12 other cruises (their
   out-of-fold predictions) is applied to the held-out cruise's within-cruise prediction; the
   within-cruise residual PCA is redone. Note: those other cruises' predictions came from models
   that saw the held-out cruise, a mild leak that can only favour removing the axis.
   Reported: W2 share, cosine with the original W2, split-half, photosynthesis z.

Outputs: transfer/results/axis2_shrinkage/ (ko_slopes.csv, sets.csv, composite.csv, recalibrated.csv,
summary.txt).

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/axis2_shrinkage_check.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from added_fields_all_kos import cv  # noqa: E402
from residual_modules import categories, split_half, standardise, top_pca  # noqa: E402
from tune_figure_kos import load  # noqa: E402

OUT = REPO / "transfer" / "results" / "axis2_shrinkage"


def centre(m, groups):
    return m - pd.DataFrame(m).groupby(groups).transform("mean").to_numpy()


def slopes(obs, pred):
    """OLS slope of obs on pred per column (both already centred)."""
    den = (pred ** 2).sum(0)
    return np.where(den > 0, (obs * pred).sum(0) / np.where(den > 0, den, 1), np.nan)


def corr(a, b):
    a, b = a - a.mean(), b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def w2_of(within_res, ref, rng, cruises, splits=20):
    z = standardise(within_res)
    ev, vec = top_pca(z)
    # the component that matches the reference W2 best among the first 10
    cos = np.abs(vec[:, :10].T @ ref) / np.linalg.norm(ref)
    j = int(cos.argmax())
    v = vec[:, j] * np.sign(vec[:, j] @ ref)
    score = z @ v / np.sqrt(ev[j])
    loading = v * np.sqrt(ev[j] / (len(z) - 1))
    rep, _ = split_half(within_res, cruises, v[:, None], rng, splits)
    return {"rank": j + 1, "share": ev[j] / ev.sum(), "cos_ref": cos[j], "split_half": rep[0]}, score, loading


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--targets", type=Path, default=REPO / "data" / "targets" / "ko_clr.csv")
    p.add_argument("--embeddings-dirs", type=Path, nargs="+",
                   default=[REPO / "data" / "embeddings_sst_l4", REPO / "data" / "faire" / "embeddings_sst_l4"])
    p.add_argument("--ref", type=Path, default=REPO / "transfer" / "results" / "residual_modules_noA")
    p.add_argument("--top", type=int, default=50)
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, _, y, cruises, unit, rows = load(args, kos, return_rows=True)
    ko_ids = [k.replace("count.ko.", "").replace("ko:", "") for k in kos]
    r2, pred = cv(emb, temp[:, 0], None, y, cruises, np.arange(len(unit)))
    yo, yp = centre(y, cruises), centre(pred, cruises)
    within = yo - yp
    ref = pd.read_csv(args.ref / "loadings.csv.gz", index_col=0).loc[ko_ids, "W2"].to_numpy()
    info, score, loading = w2_of(within, ref, rng, cruises)
    say(f"{len(unit):,} tiles, {len(kos):,} KOs; W2 = within-cruise component {info['rank']} ({info['share']:.1%}; "
        f"cos with R36 W2 {info['cos_ref']:.3f}; split-half {info['split_half']:.2f})")

    # 1. calibration slopes
    b_within = slopes(yo, yp)
    y0, p0 = y - y.mean(0), pred - pred.mean(0)
    b_pooled = slopes(y0, p0)
    brite = pd.read_csv(REPO / "data" / "reference" / "kegg_brite_ko00001.csv")
    photo = set(brite.loc[brite["C"].astype(str).str.startswith("Photosynthesis"), "ko"])
    order = np.argsort(loading)
    sets = {f"W2 + pole (top {args.top})": order[::-1][: args.top], f"W2 - pole (top {args.top})": order[: args.top],
            "KEGG photosynthesis KOs": np.array([i for i, k in enumerate(ko_ids) if k in photo]), "all KOs": np.arange(len(kos))}
    dec = pd.qcut(r2, 10, labels=False, duplicates="drop")
    tab = pd.DataFrame({"ko": ko_ids, "r2": r2, "slope_within": b_within, "slope_pooled": b_pooled, "w2_loading": loading,
                        "r2_decile": dec, "photosynthesis": [k in photo for k in ko_ids]})
    tab.round(4).to_csv(OUT / "ko_slopes.csv", index=False)
    # expected slope at the same R2: median slope of the KO's R2 decile
    dec_med = tab.groupby("r2_decile")["slope_within"].median()
    say("\n== 1. Calibration slopes (observed on cross-validated predicted; 1 = calibrated, > 1 = compressed)")
    set_rows = []
    for name, ii in sets.items():
        t = tab.iloc[ii]
        exp = dec_med.loc[t["r2_decile"]].to_numpy()
        set_rows.append({"set": name, "kos": len(t), "median_r2": t["r2"].median(), "median_slope_within": t["slope_within"].median(),
                         "median_slope_pooled": t["slope_pooled"].median(), "median_slope_same_r2_decile": np.median(exp),
                         "share_slope_within_gt_1.5": (t["slope_within"] > 1.5).mean()})
        r = set_rows[-1]
        say(f"   {name:28s} {r['kos']:5d} KOs  median R2 {r['median_r2']:.3f}  slope within {r['median_slope_within']:.2f}  "
            f"pooled {r['median_slope_pooled']:.2f}  (same-R2 KOs {r['median_slope_same_r2_decile']:.2f}); slope > 1.5: "
            f"{r['share_slope_within_gt_1.5']:.0%}")
    pd.DataFrame(set_rows).round(4).to_csv(OUT / "sets.csv", index=False)

    # 2. composites
    sd = yo.std(0)
    plus, minus = sets[f"W2 + pole (top {args.top})"], sets[f"W2 - pole (top {args.top})"]
    comp_obs = (yo[:, plus] / sd[plus]).mean(1) - (yo[:, minus] / sd[minus]).mean(1)
    comp_pred = (yp[:, plus] / sd[plus]).mean(1) - (yp[:, minus] / sd[minus]).mean(1)
    s = centre(score[:, None], cruises)[:, 0]
    b_comp = (comp_obs * comp_pred).sum() / (comp_pred ** 2).sum()
    say("\n== 2. W2 score vs the axis-2 signal, within cruise (Pearson r)")
    crow = [{"series": "KO composite (+50 minus -50)", "r_observed": corr(s, comp_obs), "r_predicted": corr(s, comp_pred),
             "r2_score_on_predicted": corr(s, comp_pred) ** 2, "calibration_slope": b_comp,
             "sd_predicted_over_observed": comp_pred.std() / comp_obs.std(), "r_obs_pred": corr(comp_obs, comp_pred)}]
    phy = pd.read_csv(REPO / "data" / "targets" / "taxa" / "taxa_phylum_clr.csv").set_index("sample_key")
    for col in [c for c in phy.columns if c.endswith("Cyanobacteriota") or c.endswith("Cryptophyta")]:
        v = rows["sample_id"].map(phy[col]).astype(float).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()
        ok = np.isfinite(v)
        _, pv = cv(emb, temp[:, 0], None, np.column_stack([v, v]), cruises, np.flatnonzero(ok))  # cv needs >= 2 targets
        pv = pv[:, :1]
        vo, vp = centre(v[ok, None], cruises[ok])[:, 0], centre(pv[ok], cruises[ok])[:, 0]
        so = s[ok]
        crow.append({"series": f"{col} CLR", "r_observed": corr(so, vo), "r_predicted": corr(so, vp),
                     "r2_score_on_predicted": corr(so, vp) ** 2, "calibration_slope": (vo * vp).sum() / (vp ** 2).sum(),
                     "sd_predicted_over_observed": vp.std() / vo.std(), "r_obs_pred": corr(vo, vp)})
    for r in crow:
        say(f"   {r['series']:40s} r(W2, observed) {r['r_observed']:+.2f}  r(W2, predicted) {r['r_predicted']:+.2f}  "
            f"[obs-pred r {r['r_obs_pred']:+.2f}; SD pred/obs {r['sd_predicted_over_observed']:.2f}; calibration slope "
            f"{r['calibration_slope']:.2f}]")
    pd.DataFrame(crow).round(4).to_csv(OUT / "composite.csv", index=False)

    # 3. recalibration
    say("\n== 3. Recalibrated predictions (per-KO within-cruise slope from the other 12 cruises), residual PCA redone")
    cal = np.empty_like(yp)
    for c in np.unique(cruises):
        te, tr = cruises == c, cruises != c
        b = slopes(yo[tr], yp[tr])
        cal[te] = yp[te] * np.nan_to_num(b, nan=1.0)
    within_cal = yo - cal
    r2_w = 1 - (within ** 2).sum(0) / (yo ** 2).sum(0)
    r2_wc = 1 - (within_cal ** 2).sum(0) / (yo ** 2).sum(0)
    info_c, score_c, loading_c = w2_of(within_cal, ref, rng, cruises)
    cat0 = categories(loading, ko_ids, brite).set_index(["level", "category"])["z"]
    cat1 = categories(loading_c, ko_ids, brite).set_index(["level", "category"])["z"]
    ph0, ph1 = cat0.get(("C", "Photosynthesis"), np.nan), cat1.get(("C", "Photosynthesis"), np.nan)
    rec = [{"residuals": "baseline", **info, "photosynthesis_z": ph0, "median_within_r2": np.median(r2_w)},
           {"residuals": "recalibrated", **info_c, "photosynthesis_z": ph1, "median_within_r2": np.median(r2_wc),
            "cos_loadings_with_baseline_w2": float(loading_c @ loading / np.linalg.norm(loading_c) / np.linalg.norm(loading)),
            "r_scores_with_baseline_w2": corr(score_c, score)}]
    for r in rec:
        say(f"   {r['residuals']:13s} W2 = component {r['rank']} ({r['share']:.1%}); cos with R36 W2 {r['cos_ref']:.3f}; split-half "
            f"{r['split_half']:.2f}; Photosynthesis z {r['photosynthesis_z']:+.1f}; median within-cruise KO R2 {r['median_within_r2']:.3f}")
    say(f"   recalibrated vs baseline W2: loading cosine {rec[1]['cos_loadings_with_baseline_w2']:.3f}, score r {rec[1]['r_scores_with_baseline_w2']:.3f}")
    pd.DataFrame(rec).round(4).to_csv(OUT / "recalibrated.csv", index=False)
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
