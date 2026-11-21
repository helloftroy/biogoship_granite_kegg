#!/usr/bin/env python3
"""Two-stage model: satellite (+ mixing history) -> domain balance C -> gene content; scores with technical
variance removed from the denominator (RESULTS.md R54).

Story being tested: environment selects who is there, and who is there determines which genes are present.
  stage 1  predict C (log eukaryote-specific / prokaryote-specific reads, R51) from
             base      32 Granite PCs + in-situ temperature spline (the R51 two-stage)
             sat       16 OLCI band means + satellite SST, splines (best R51 single predictor of C)
             sat_mld   sat + 30-day mean mixed-layer depth along the back-trajectory (spline; R53)
  stage 2  each KO target from Granite features + a spline of the stage-1 C prediction.
           Stage-2 features: baseline (32 Granite PCs + temperature spline) or granite (32 Granite PCs).
Nested, no leakage: for held-out cruise c, stage 1 is fitted on the other cruises; the stage-1 values used
for the TRAINING rows are themselves leave-one-cruise-out predictions within those training cruises (R51 used
the outer out-of-fold C for training rows, so the test cruise had informed them; both are reported for the
R51 configuration).

Technical variance out of the denominator ("tech" mode): within each fold, on the training tiles only, an OLS
per target on platform (NovaSeq 6000 vs NovaSeq X cruise groups, as R43) + a cubic spline (8 knots) of log10
KO reads per tile; the fitted technical part is subtracted from training AND test targets (test with the
training coefficients); the model is fitted and scored on the adjusted targets. Explained = 1 - SS(resid) /
SS(adjusted targets - mean), pooled over the held-out predictions. The denominator shrinks by the technical
variance; the share is "of the biological (non-technical) variance". Platform is cruise-level and confounded
with year and region, so this removes an upper bound of the platform effect.

Targets (R51): clr_all (current CLR, 3,933 KOs), A, B, C, D, Dp, De.

Step "spline" (default part 1): the nested stage 2 above, C prediction entered as a spline. It turned out
fragile (R54): the stage-2 spline is fitted on inner-fold stage-1 values whose distribution differs from the
test values, and a placebo (held-out cruise's C scrambled inside the R51 non-nested scheme) drops skill to
0.4 %. So part 2 compares stage-2 forms:
  nested_linear   stage-2 term = standardised (not splined) nested stage-1 prediction
  plugin_spline   stage 2 fitted on the OBSERVED C of the training tiles (known there), spline; the test
                  cruise gets its stage-1 prediction (fitted on the training cruises)
  plugin_linear   the same, linear term  <- the form used for the headline
and part 3 gives the headline (Granite PCs vs plug-in linear two-stage with C from OLCI + SST + 30-day MLD),
plain and with technical variance removed.

Outputs: data/results/two_stage/{scores.csv,stage1_C.csv,summary.txt} (part 1), stage2_forms.csv (part 2),
plugin_headline.csv (part 3); figures/output/two_stage.png (from part 3 and part 1 baselines).

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/two_stage_models.py --parts spline forms headline
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import SplineTransformer, StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from domain_split_models import EMB, FSETS, NOVASEQ_6000  # noqa: E402
from tune_figure_kos import build, load, ridge_gcv  # noqa: E402
from water_mass_test import features  # noqa: E402

OUT = REPO / "data" / "results" / "two_stage"
FIG = REPO / "figures" / "output"
TD = REPO / "data" / "targets" / "domain_split"
TARGETS = {"clr_all": REPO / "data" / "targets" / "ko_clr.csv", **{k: TD / f"{k}.csv" for k in ["A", "B", "C", "D", "Dp", "De"]}}
S1_CACHE: dict = {}
STAGE1 = {"base": (FSETS["baseline"], False), "sat": (FSETS["raw_sat"], False), "sat_mld": (FSETS["raw_sat"], True)}


def spl(xtr, xte):
    s = StandardScaler().fit(xtr)
    sp = SplineTransformer(n_knots=8, degree=3, extrapolation="linear").fit(s.transform(xtr))
    return sp.transform(s.transform(xtr)), sp.transform(s.transform(xte))


def fit_predict(blocks, tr, te, D, y, extra=None):
    """Ridge (per-target GCV alpha) on feature blocks (+ splines of extra columns), trained on rows tr."""
    emb, temp, sat = D
    xtr, xte = build(blocks, tr, te, emb, temp, sat, {})
    if extra is not None:
        etr, ete = extra
        a, b = spl(etr, ete)
        xtr, xte = np.hstack([xtr, a]), np.hstack([xte, b])
    return np.asarray(ridge_gcv(xtr, y[tr], xte)).reshape(len(te), -1)


def stage1(kind, tr, te, D, c, mld, cr):
    """Stage-1 C for training rows (inner leave-one-cruise-out) and test rows (fit on all training rows)."""
    blocks, use_mld = STAGE1[kind]
    ex = (lambda a, b: (mld[a], mld[b])) if use_mld else (lambda a, b: None)
    pte = fit_predict(blocks, tr, te, D, c, ex(tr, te))[:, 0]
    ptr = np.empty(len(tr))
    for k in np.unique(cr[tr]):
        inner_te = tr[cr[tr] == k]
        inner_tr = tr[cr[tr] != k]
        ptr[cr[tr] == k] = fit_predict(blocks, inner_tr, inner_te, D, c, ex(inner_tr, inner_te))[:, 0]
    return ptr, pte


def tech_adjust(y, tr, te, plat, logtot):
    """Subtract platform + depth-spline fit (training rows only) from y; returns adjusted (y_tr, y_te)."""
    sp = SplineTransformer(n_knots=8, degree=3).fit(logtot[tr, None])
    X = lambda i: np.hstack([np.ones((len(i), 1)), plat[i, None].astype(float), sp.transform(logtot[i, None])])  # noqa: E731
    beta, *_ = np.linalg.lstsq(X(tr), y[tr], rcond=None)
    mu = y[tr].mean(0)
    return y[tr] - (X(tr) @ beta - mu), y[te] - (X(te) @ beta - mu)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    fsamp = features()
    csamp = pd.read_csv(TD / "C.csv", index_col="sample_key")["euk_prok_log_ratio"]
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    rows_out, c1_rows = [], []
    for name, path in TARGETS.items():
        kos = pd.read_csv(path, nrows=0).columns.drop("sample_key").tolist()
        emb, temp, sat, y, cr, unit, rows = load(argparse.Namespace(targets=path, embeddings_dirs=EMB), kos, return_rows=True)
        D = (emb, temp, sat)
        tile = lambda s: s.reindex(rows["sample_id"]).set_axis(rows.index).groupby(rows["tile_key"]).mean().reindex(unit.index)  # noqa: E731
        mld = tile(fsamp["mld_mean30"])
        mld = mld.fillna(mld.median()).to_numpy()[:, None]
        c = tile(csamp).to_numpy()[:, None]
        logtot = tile(np.log10(ks["total"])).to_numpy()
        plat = np.isin(cr, NOVASEQ_6000)
        say(f"\n{name}: {y.shape[0]:,} tiles x {y.shape[1]:,}")
        configs = [("baseline", None), ("granite", None)]
        if name == "C":
            configs = [("baseline", None), ("granite", None), ("raw_sat", None), ("raw_sat+mld30", None)]
        else:
            configs += [(f"{s2} + C[{k}]", k) for s2 in ("baseline", "granite") for k in STAGE1]
            if name == "clr_all":
                configs += [("baseline + C[base] (R51, not nested)", "base_leaky")]
        for mode in ("plain", "tech"):
            preds = {cfg: np.full_like(y, np.nan) for cfg, _ in configs}
            yadj = np.full_like(y, np.nan)
            c1 = {k: np.full(len(y), np.nan) for k in STAGE1}
            leaky = None
            if any(k == "base_leaky" for _, k in configs):
                leaky = np.full(len(y), np.nan)
                for k in np.unique(cr):
                    te, tr = np.where(cr == k)[0], np.where(cr != k)[0]
                    leaky[te] = fit_predict(FSETS["baseline"], tr, te, D, c)[:, 0]
            for k in np.unique(cr):
                te, tr = np.where(cr == k)[0], np.where(cr != k)[0]
                if mode == "tech":
                    ytr, yte = tech_adjust(y, tr, te, plat, logtot)
                else:
                    ytr, yte = y[tr], y[te]
                yadj[te] = yte
                Y = np.zeros_like(y)
                Y[tr] = ytr
                s1 = {}
                for cfg, kind in configs:
                    if kind is None:
                        blocks = {"baseline": FSETS["baseline"], "granite": FSETS["granite"], "raw_sat": FSETS["raw_sat"],
                                  "raw_sat+mld30": FSETS["raw_sat"]}[cfg]  # fmt: skip
                        ex = (mld[tr], mld[te]) if cfg == "raw_sat+mld30" else None
                        preds[cfg][te] = fit_predict(blocks, tr, te, D, Y, ex)
                        continue
                    if kind == "base_leaky":
                        ex = (leaky[tr, None], leaky[te, None])
                    else:
                        if kind not in s1:
                            key = (len(y), kind, k)
                            if key not in S1_CACHE:
                                S1_CACHE[key] = stage1(kind, tr, te, D, c, mld, cr)
                            s1[kind] = S1_CACHE[key]
                            c1[kind][te] = s1[kind][1]
                        ex = (s1[kind][0][:, None], s1[kind][1][:, None])
                    s2 = "baseline" if cfg.startswith("baseline") else "granite"
                    preds[cfg][te] = fit_predict(FSETS[s2], tr, te, D, Y, ex)
            yc = yadj - yadj.mean(0)
            sd = np.where(yc.std(0) > 0, yc.std(0), 1)
            tot_raw = ((y - y.mean(0)) ** 2).sum()
            for cfg, _ in configs:
                r = yadj - preds[cfg]
                e_raw = 1 - (r**2).sum() / (yc**2).sum()
                e_std = 1 - ((r / sd) ** 2).sum() / ((yc / sd) ** 2).sum()
                rows_out.append({"target": name, "mode": mode, "model": cfg, "tiles": len(y), "kos": y.shape[1],
                                 "explained_raw": e_raw, "explained_std": e_std, "denominator_share": (yc**2).sum() / tot_raw})  # fmt: skip
                say(f"   [{mode}] {cfg:38s} explained raw {e_raw:6.2%}  std {e_std:6.2%}   (denominator {(yc ** 2).sum() / tot_raw:.1%} of total)")
            if mode == "plain" and name == "clr_all":
                for kind, v in c1.items():
                    ok = ~np.isnan(v)
                    e = 1 - ((c[ok, 0] - v[ok]) ** 2).sum() / ((c[ok, 0] - c[ok, 0].mean()) ** 2).sum()
                    c1_rows.append({"stage1": kind, "C_explained_outer": e})
                    say(f"   stage-1 C [{kind}] explained (outer, held-out cruise) {e:.2%}")
    sc = pd.DataFrame(rows_out)
    sc.round(5).to_csv(OUT / "scores.csv", index=False)
    pd.DataFrame(c1_rows).round(5).to_csv(OUT / "stage1_C.csv", index=False)
    for mode in ("plain", "tech"):
        say(f"\nExplained share, raw units, {mode}:")
        say((sc[sc["mode"] == mode].pivot(index="model", columns="target", values="explained_raw")[list(TARGETS)] * 100).round(1).to_string())
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")


def figure(sc):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    tg = ["clr_all", "D", "A", "Dp", "B", "De"]
    tl = {"clr_all": "Current CLR\n(all KOs)", "D": "D\nShared CLR", "A": "A\nProk-specific\nper GE", "Dp": "Dp\nShared, prok\nper GE",
          "B": "B\nEuk-specific\nwithin euks", "De": "De\nShared, euk\nwithin euks"}  # fmt: skip
    models = ["granite", "granite + C[base]", "granite + C[sat]", "granite + C[sat_mld]"]
    ml = {"granite": "Granite PCs", "granite + C[base]": "+ C from Granite + temperature", "granite + C[sat]": "+ C from OLCI + SST",
          "granite + C[sat_mld]": "+ C from OLCI + SST + 30-d MLD"}  # fmt: skip
    cols = ["#0072B2", "#999999", "#56B4E9", "#E69F00"]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for ax, mode, lab in zip(axs, ("plain", "tech"), ("of all variance", "of non-technical variance")):
        s = sc[sc["mode"] == mode]
        for j, m in enumerate(models):
            v = s[s["model"] == m].set_index("target")["explained_raw"].reindex(tg).to_numpy() * 100
            ax.bar(np.arange(len(tg)) + (j - 1.5) * 0.2, v, 0.19, color=cols[j], label=ml[m])
        ax.axhline(0, color="#666666", lw=0.5)
        ax.set_xticks(range(len(tg)))
        ax.set_xticklabels([tl[t] for t in tg], fontsize=6.5)
        ax.set_ylabel(f"Variance explained, cruise held out\n(% {lab})")
    axs[0].legend(frameon=False, fontsize=6.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "two_stage.png", dpi=300)
    plt.close(fig)


def _setup(tname, fs, ks):
    path = TARGETS[tname]
    kos = pd.read_csv(path, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, sat, y, cr, unit, rows = load(argparse.Namespace(targets=path, embeddings_dirs=EMB), kos, return_rows=True)
    tile = lambda s: s.reindex(rows["sample_id"]).set_axis(rows.index).groupby(rows["tile_key"]).mean().reindex(unit.index)  # noqa: E731
    c = tile(pd.read_csv(TD / "C.csv", index_col="sample_key")["euk_prok_log_ratio"]).to_numpy()[:, None]
    mld = tile(fs["mld_mean30"])
    mld = mld.fillna(mld.median()).to_numpy()[:, None]
    return emb, temp, sat, y, cr, c, mld, tile(np.log10(ks["total"])).to_numpy()


def forms() -> None:
    fs = features()
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    out = []
    for tname in ["clr_all", "D", "A"]:
        emb, temp, sat, y, cr, c, mld, _ = _setup(tname, fs, ks)
        D = (emb, temp, sat)
        preds = {}
        for kind in ["base", "sat_mld"]:
            for form in ["nested_linear", "plugin_spline", "plugin_linear"]:
                for s2 in ["baseline", "granite"]:
                    preds[(kind, form, s2)] = np.full_like(y, np.nan)
            for k in np.unique(cr):
                te, tr = np.where(cr == k)[0], np.where(cr != k)[0]
                ptr, pte = stage1(kind, tr, te, D, c, mld, cr)
                for s2 in ["baseline", "granite"]:
                    xtr, xte = build(FSETS[s2], tr, te, emb, temp, sat, {})
                    for form in ["nested_linear", "plugin_spline", "plugin_linear"]:
                        ctr = ptr[:, None] if form == "nested_linear" else c[tr]
                        cte = pte[:, None]
                        if form.endswith("linear"):
                            sc = StandardScaler().fit(ctr)
                            a, b = sc.transform(ctr), sc.transform(cte)
                        else:
                            a, b = spl(ctr, cte)
                        preds[(kind, form, s2)][te] = np.asarray(ridge_gcv(np.hstack([xtr, a]), y[tr], np.hstack([xte, b]))).reshape(len(te), -1)
        yc = y - y.mean(0)
        for key, pr in preds.items():
            e = 1 - ((y - pr) ** 2).sum() / (yc**2).sum()
            out.append((tname, *key, e))
            print(tname, key, f"{e:.2%}", flush=True)
    pd.DataFrame(out, columns=["target", "stage1", "form", "stage2", "explained_raw"]).round(5).to_csv(OUT / "stage2_forms.csv", index=False)


def headline() -> None:
    fs = features()
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    out = []
    for tname in ["clr_all", "D", "A", "Dp", "B", "De", "C"]:
        emb, temp, sat, y, cr, c, mld, logtot = _setup(tname, fs, ks)
        D = (emb, temp, sat)
        plat = np.isin(cr, NOVASEQ_6000)
        for mode in ["plain", "tech"]:
            P = {m: np.full_like(y, np.nan) for m in ["granite", "two_stage"]}
            ya = np.full_like(y, np.nan)
            for k in np.unique(cr):
                te, tr = np.where(cr == k)[0], np.where(cr != k)[0]
                ytr, yte = tech_adjust(y, tr, te, plat, logtot) if mode == "tech" else (y[tr], y[te])
                ya[te] = yte
                xtr, xte = build(FSETS["granite"], tr, te, emb, temp, sat, {})
                P["granite"][te] = np.asarray(ridge_gcv(xtr, ytr, xte)).reshape(len(te), -1)
                if tname != "C":
                    key = (len(y), "sat_mld", k)
                    if key not in S1_CACHE:
                        S1_CACHE[key] = stage1("sat_mld", tr, te, D, c, mld, cr)
                    pte = S1_CACHE[key][1]
                    sc = StandardScaler().fit(c[tr])
                    P["two_stage"][te] = np.asarray(ridge_gcv(np.hstack([xtr, sc.transform(c[tr])]), ytr,
                                                              np.hstack([xte, sc.transform(pte[:, None])]))).reshape(len(te), -1)  # fmt: skip
            yc = ya - ya.mean(0)
            sd = np.where(yc.std(0) > 0, yc.std(0), 1)
            for m, pr in P.items():
                if np.isnan(pr).all():
                    continue
                r = ya - pr
                out.append((tname, mode, m, 1 - (r**2).sum() / (yc**2).sum(), 1 - ((r / sd) ** 2).sum() / ((yc / sd) ** 2).sum(),
                            (yc**2).sum() / ((y - y.mean(0)) ** 2).sum()))  # fmt: skip
                print(tname, mode, m, f"{out[-1][3]:.2%} std {out[-1][4]:.2%}", flush=True)
    pd.DataFrame(out, columns=["target", "mode", "model", "explained_raw", "explained_std", "denominator_share"]).round(5).to_csv(
        OUT / "plugin_headline.csv", index=False)
    headline_figure()


def headline_figure() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    h = pd.read_csv(OUT / "plugin_headline.csv")
    b = pd.read_csv(OUT / "scores.csv")
    b = b[b["model"] == "baseline"].assign(model="baseline")[["target", "mode", "model", "explained_raw"]]
    sc = pd.read_csv(OUT / "scores.csv")
    c1 = sc[(sc["target"] == "C") & (sc["model"] == "raw_sat+mld30")].assign(model="two_stage")[["target", "mode", "model", "explained_raw"]]
    h = pd.concat([b, h[["target", "mode", "model", "explained_raw"]], c1])
    tg = ["clr_all", "C", "D", "A", "Dp", "B", "De"]
    tl = {"clr_all": "All KOs\n(CLR)", "C": "C", "D": "D", "A": "A", "Dp": "Dp", "B": "B", "De": "De"}
    models = [("baseline", "Granite PCs + temperature", "#999999"), ("granite", "Granite PCs", "#0072B2"),
              ("two_stage", "Two-stage, C from OLCI + SST + 30-d MLD\n(for C: that model itself)", "#E69F00")]  # fmt: skip
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1))
    for ax, mode, lab in zip(axs, ("plain", "tech"), ("of all variance", "of non-technical variance")):
        s = h[h["mode"] == mode]
        for j, (m, ml, col) in enumerate(models):
            v = s[s["model"] == m].set_index("target")["explained_raw"].reindex(tg).to_numpy() * 100
            ax.bar(np.arange(len(tg)) + (j - 1) * 0.26, v, 0.24, color=col, label=ml)
        ax.axhline(0, color="#666666", lw=0.5)
        ax.set_xticks(range(len(tg)))
        ax.set_xticklabels([tl[t] for t in tg], fontsize=7)
        ax.set_ylim(0, 65)
        ax.set_xlabel("Target")
        ax.set_ylabel(f"Variance explained, cruise held out\n(% {lab})")
    axs[0].legend(frameon=False, fontsize=6.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "two_stage.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", default=["spline", "forms", "headline"], choices=["spline", "forms", "headline", "figure"])
    a = ap.parse_args()
    if "spline" in a.parts:
        main()
    if "forms" in a.parts:
        forms()
    if "headline" in a.parts:
        headline()
    if "figure" in a.parts:
        headline_figure()
