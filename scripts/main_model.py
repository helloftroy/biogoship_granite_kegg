#!/usr/bin/env python3
"""Main model for the paper: Granite PCs alone, ridge, cruise held out; C as its own result (RESULTS.md R59).

Replaces the chained (two-stage) models of R51 / R54 / R56 as the headline.

Tiles: the R43 set (1,215 tiles; 1,169 for B and De), samples sharing a Granite tile averaged into one row.
Targets: clr_all (all 3,933 KOs, CLR; comparison row), A, B, C, D, Dp, De (R51 domain split).
Technical removal: in-sample, as R56: per target column, OLS on platform (NovaSeq 6000 cruise group) + cubic
spline (8 knots) of log10 KO reads per tile, fitted on all tiles; the fitted part is removed and every model is
fitted and scored on the adjusted targets.

Models, each cruise held out in turn (13 folds); features standardised inside each fold:
  granite        32 Granite PCs (PCA fitted on the training tiles)                      <- KO targets
  baseline       32 Granite PCs + in-situ temperature spline (the old R17 baseline, uses ship data)
  temp           in-situ temperature spline alone (for the per-KO scatter and heatmap)
  sat            16 OLCI band means + satellite SST, quantile-normalised, splines (R51 raw_sat)
  sat_mld        sat + spline of the 30-day mean mixed-layer depth along the back-trajectory  <- C
Ridge penalty, three ways:
  gcv            per target column, efficient leave-one-tile-out GCV (RidgeCV; R51-R56)
  cv             per target column, leave-one-cruise-out inside the training cruises
  cv_shared      one penalty for the whole target set, leave-one-cruise-out inside the training cruises
R38 found the held-out predictions over-dispersed (observed-on-predicted slope < 1). The calibration slope is
reported for each.

Parts:
  fit        scores.csv (pooled explained, calibration slope, penalty), per-KO and per-cruise R2 of clr_all in
             the ko_temperature.csv layout (for figures/make_temperature_gap.py and make_ko_heatmap.py),
             predictions_*.npz (clr_all and C)
  partition  partition.csv: technical / explained / cruise offsets / within-tile variability / within-cruise,
             summing to 100 %, for the main model of each target; within_tile_k.csv: within-tile variability and
             held-out error by samples per tile
  figure     figures/output/variance_partition.png (main figure), figures/output/c_prediction.png (C)

Outputs: data/results/main_model/
Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/main_model.py --parts fit partition figure
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import SplineTransformer

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from domain_split_models import EMB, NOVASEQ_6000  # noqa: E402
from mld_platform_partition import replicate_noise  # noqa: E402
from screen_ko_predictability import ALPHAS, metadata_features  # noqa: E402
from tune_figure_kos import build, load, ridge_path  # noqa: E402
from two_stage_models import TARGETS, spl  # noqa: E402
from water_mass_test import features  # noqa: E402

OUT = REPO / "data" / "results" / "main_model"
FIG = REPO / "figures" / "output"
TD = REPO / "data" / "targets" / "domain_split"
ORDER = ["clr_all", "C", "D", "A", "Dp", "B", "De"]
BLOCKS = {"granite": ["P32"], "baseline": ["P32", "T~"], "temp": ["T~"], "sat": ["S~"], "sat_mld": ["S~"]}
# "place": time & place (lat / 90, |lat| / 90, sin and cos of longitude, sin and cos of day of year; standardised,
# linear). R5's version without |lat| cannot follow poleward increases and scores below 0 (R60)
NAMES = {"clr_all": "All KOs", "C": "Eukaryote :\nprokaryote\nreads", "D": "Shared KOs", "A": "Prokaryote-\nspecific KOs",
         "Dp": "Shared KOs,\nprokaryote\npart", "B": "Eukaryote-\nspecific KOs", "De": "Shared KOs,\neukaryote\npart"}  # fmt: skip
PENALTIES = ["gcv", "cv", "cv_shared"]
LOG: list[str] = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def setup(tname, fs, ks):
    path = TARGETS[tname]
    kos = pd.read_csv(path, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, sat, y, cr, unit, rows = load(argparse.Namespace(targets=path, embeddings_dirs=EMB), kos, return_rows=True)
    tile = lambda s: s.reindex(rows["sample_id"]).set_axis(rows.index).groupby(rows["tile_key"]).mean().reindex(unit.index)  # noqa: E731
    mld = tile(fs["mld_mean30"])
    mld = mld.fillna(mld.median()).to_numpy()[:, None]
    logtot = tile(np.log10(ks["total"])).to_numpy()
    plat = np.isin(cr, NOVASEQ_6000)
    X = np.hstack([np.ones((len(y), 1)), plat[:, None].astype(float), SplineTransformer(n_knots=8, degree=3).fit_transform(logtot[:, None])])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ beta
    ya = y - (fit - y.mean(0))
    k = rows.groupby("tile_key").size().reindex(unit.index).to_numpy()
    return dict(emb=emb, temp=temp, sat=sat, y=y, ya=ya, tech=fit - fit.mean(0), cr=cr, unit=unit, rows=rows, mld=mld,
                logtot=logtot, plat=plat, kos=kos, k=k)  # fmt: skip


def features_for(model, d, tr, te, cache):
    if model == "place":
        x = np.hstack([metadata_features(d["unit"]), np.abs(d["unit"]["sample_lat"].to_numpy())[:, None] / 90])
        mu, sd = x[tr].mean(0), x[tr].std(0)
        return (x[tr] - mu) / sd, (x[te] - mu) / sd
    xtr, xte = build(BLOCKS[model], tr, te, d["emb"], d["temp"], d["sat"], cache)
    if model == "sat_mld":
        a, b = spl(d["mld"][tr], d["mld"][te])
        xtr, xte = np.hstack([xtr, a]), np.hstack([xte, b])
    return xtr, xte


def ridge(xtr, ytr, xte, penalty, groups):
    """Predictions for the test rows and the penalty used per target column."""
    if penalty == "gcv":
        m = RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(xtr, ytr)
        return m.predict(xte).reshape(len(xte), -1), np.broadcast_to(m.alpha_, (ytr.shape[1],))
    sse = np.zeros((len(ALPHAS), ytr.shape[1]))
    for g in np.unique(groups):
        a, b = groups != g, groups == g
        sse += ((ytr[b][None] - ridge_path(xtr[a], ytr[a], xtr[b])) ** 2).sum(axis=1)
    idx = sse.argmin(0) if penalty == "cv" else np.full(ytr.shape[1], sse.sum(1).argmin())
    path = ridge_path(xtr, ytr, xte)
    return path[idx, :, np.arange(ytr.shape[1])].T, ALPHAS[idx]


def loco(d, y, models):
    """Leave-one-cruise-out predictions for every (model, penalty); also the chosen penalties."""
    cr = d["cr"]
    pred = {mp: np.full_like(y, np.nan) for mp in models}
    alpha = {mp: [] for mp in models}
    for c in np.unique(cr):
        te, tr = np.where(cr == c)[0], np.where(cr != c)[0]
        cache: dict = {}
        for model, pen in models:
            xtr, xte = features_for(model, d, tr, te, cache)
            p, a = ridge(xtr, y[tr], xte, pen, cr[tr])
            pred[(model, pen)][te] = p
            alpha[(model, pen)].append(a)
    return pred, {mp: np.log10(np.vstack(v)) for mp, v in alpha.items()}


def explained(y, p):
    return 1 - ((y - p) ** 2).sum() / ((y - y.mean(0)) ** 2).sum()


def calib_slope(y, p):
    """Slope of observed on predicted, per column (centred); < 1 means the predictions are over-dispersed."""
    pc, yc = p - p.mean(0), y - y.mean(0)
    return (pc * yc).sum(0) / np.maximum((pc**2).sum(0), 1e-12)


def per_ko_tables(y, preds, kos, cr):
    pooled, cruise = [], []
    for (fs, learner), p in preds.items():
        r2 = 1 - ((y - p) ** 2).sum(0) / ((y - y.mean(0)) ** 2).sum(0)
        rows = []
        for c in np.unique(cr):
            te = cr == c
            r = 1 - ((y[te] - p[te]) ** 2).sum(0) / ((y[te] - y[~te].mean(0)) ** 2).sum(0)
            rows.append(pd.DataFrame({"ko": kos, "feature_set": fs, "learner": learner, "held_out": c, "n": int(te.sum()), "r2": r}))
        rows = pd.concat(rows)
        cruise.append(rows)
        med = rows.groupby("ko", sort=False)["r2"].median().reindex(kos).to_numpy()
        pooled.append(pd.DataFrame({"ko": kos, "feature_set": fs, "learner": learner, "r2_pooled": r2, "r2_median_cruise": med}))
    return pd.concat(pooled), pd.concat(cruise)


def fit_part(fs, ks):
    out = []
    for tname in ORDER:
        d = setup(tname, fs, ks)
        if tname == "C":
            models = [(m, "cv") for m in ("granite", "sat", "sat_mld", "place", "temp", "baseline")] + [("granite", "gcv")]
        else:
            models = [("granite", p) for p in PENALTIES] + [(m, "cv") for m in ("baseline", "temp", "sat", "place")]
        pred, alpha = loco(d, d["ya"], models)
        for (m, p), v in pred.items():
            s = calib_slope(d["ya"], v)
            out.append({"target": tname, "model": m, "penalty": p, "tiles": len(v), "columns": v.shape[1],
                        "explained": explained(d["ya"], v), "explained_unadjusted_targets": np.nan,
                        "calib_slope_pooled": float(calib_slope(d["ya"].ravel()[:, None], v.ravel()[:, None])[0]),
                        "calib_slope_median": float(np.median(s)), "log10_alpha_median": float(np.median(alpha[(m, p)]))})  # fmt: skip
            say(f"{tname:8s} {m:9s} {p:9s} explained {out[-1]['explained']:7.2%}  slope pooled "
                f"{out[-1]['calib_slope_pooled']:.3f} (median {out[-1]['calib_slope_median']:.3f})  "
                f"log10 alpha median {out[-1]['log10_alpha_median']:.2f}")  # fmt: skip
        if tname == "clr_all":
            keep = {("embedding", "ridge"): pred[("granite", "cv")], ("embedding", "ridge_gcv"): pred[("granite", "gcv")],
                    ("embedding", "ridge_cv_shared"): pred[("granite", "cv_shared")], ("temp", "spline"): pred[("temp", "cv")],
                    ("metadata", "ridge"): pred[("place", "cv")], ("satellite", "spline"): pred[("sat", "cv")],
                    ("pca32_temp", "ridge"): pred[("baseline", "cv")]}  # fmt: skip
            pooled, cruise = per_ko_tables(d["ya"], keep, d["kos"], d["cr"])
            pooled.round(4).to_csv(OUT / "ko_temperature.csv", index=False)
            cruise.round(4).to_csv(OUT / "ko_temperature_per_cruise.csv", index=False)
        for (m, p), v in pred.items():
            np.savez_compressed(OUT / f"predictions_{tname}_{m}_{p}.npz", predicted=v.astype(np.float32),
                                observed=d["ya"].astype(np.float32), ko=np.array(d["kos"]), cruise=d["cr"],
                                tile=d["unit"].index.to_numpy(), k=d["k"])  # fmt: skip
    t = pd.DataFrame(out).drop(columns="explained_unadjusted_targets")
    t.round(5).to_csv(OUT / "scores.csv", index=False)
    say("\n" + (t.pivot_table(index=["target", "model"], columns="penalty", values="explained", sort=False) * 100).round(2).to_string())


MAIN = {"C": ("granite", "cv")}  # set after the fit part; KO targets: see main_choice()


def main_choice(tname):
    return ("granite", MAIN.get("ko_penalty", "cv"))


def partition_part(fs, ks):
    out, kk = [], []
    for label, tname, (m, p) in [(t, t, main_choice(t)) for t in ORDER]:
        d = setup(tname, fs, ks)
        z = np.load(OUT / f"predictions_{tname}_{m}_{p}.npz")
        pred, ya, y, cr = z["predicted"].astype(float), d["ya"], d["y"], d["cr"]
        tot = ((y - y.mean(0)) ** 2).sum()
        tech = (d["tech"] ** 2).sum()
        yac = ya - ya.mean(0)
        res = ya - pred
        expl = (yac**2).sum() - (res**2).sum()
        res = res - res.mean(0)
        cm = pd.DataFrame(res).groupby(cr).transform("mean").to_numpy()
        offs = (cm**2).sum()
        within = ((res - cm) ** 2).sum()
        st = pd.read_csv(TARGETS[tname], index_col="sample_key")
        st = st["euk_prok_log_ratio"].to_frame() if tname == "C" else st[d["kos"]]
        noise, npairs = replicate_noise(y, d, st)
        noise = min(noise, within)
        r = {"target": label, "model": m, "penalty": p, "tiles": len(y), "columns": y.shape[1], "technical": tech / tot,
             "explained": expl / tot, "cruise_offsets": offs / tot, "within_tile": noise / tot, "within_cruise": (within - noise) / tot,
             "within_tile_pairs": npairs}  # fmt: skip
        r["sum"] = r["technical"] + r["explained"] + r["cruise_offsets"] + r["within_tile"] + r["within_cruise"]
        r["of_nontechnical"] = r["explained"] / (1 - r["technical"])
        r["of_reproducible_nontechnical"] = r["explained"] / (1 - r["technical"] - r["within_tile"])
        out.append(r)
        say(f"{label:9s} {m}/{p}: tech {r['technical']:6.2%}  explained {r['explained']:6.2%}  cruise offsets "
            f"{r['cruise_offsets']:6.2%}  within-tile {r['within_tile']:6.2%}  within-cruise {r['within_cruise']:6.2%}  "
            f"(sum {r['sum']:.4f}) -> {r['of_nontechnical']:6.2%} of non-technical, "
            f"{r['of_reproducible_nontechnical']:6.2%} of reproducible non-technical")  # fmt: skip
        if label == tname:
            kk.append(by_k(tname, d, st, res + 0, pred))
    pd.DataFrame(out).round(5).to_csv(OUT / "partition.csv", index=False)
    kk = pd.concat(kk)
    kk.round(6).to_csv(OUT / "within_tile_k.csv", index=False)
    say("\nWithin-tile variability and held-out error by samples per tile (sums over target columns; "
        "sigma2 = single-sample within-tile variance from the pairs of that group; 'expected' = sigma2 / k):")
    say(kk.round(4).to_string(index=False))


def by_k(tname, d, st, _res, pred):
    """Per k group: single-sample sigma2 from that group's pairs, sigma2/k, and the held-out MSE of the tile mean
    (main model, adjusted target), all tiles and GOMECC4 only (most multi-sample tiles are GOMECC4)."""
    rows, unit, k, cr = d["rows"], d["unit"], d["k"], d["cr"]
    members = rows.groupby("tile_key")["sample_id"].agg(list).reindex(unit.index)
    err = ((d["ya"] - pred) ** 2).sum(1)
    grp = np.select([k == 1, k == 2, k == 3, k >= 4], ["1", "2", "3", "4-6"])
    out = []
    for g in ["1", "2", "3", "4-6"]:
        sel = grp == g
        diffs = []
        for samples in members[sel]:
            v = st.reindex([x for x in samples if x in st.index]).to_numpy()
            for i in range(len(v)):
                for j in range(i + 1, len(v)):
                    diffs.append(v[i] - v[j])
        s2 = float((np.array(diffs) ** 2).mean(0).sum() / 2) if diffs else np.nan
        kmean = float(k[sel].mean())
        go = sel & (cr == "GOMECC4")
        out.append({"target": tname, "samples_per_tile": g, "tiles": int(sel.sum()), "pairs": len(diffs), "sigma2_single": s2,
                    "sigma2_over_k": s2 / kmean if diffs else np.nan, "heldout_mse": float(err[sel].mean()),
                    "tiles_gomecc4": int(go.sum()), "heldout_mse_gomecc4": float(err[go].mean()) if go.any() else np.nan})  # fmt: skip
    return pd.DataFrame(out)


def figure_part():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})  # fmt: skip
    p = pd.read_csv(OUT / "partition.csv").set_index("target")
    tl = NAMES
    parts = [("explained", "Explained by satellite (cruise held out)", "#E69F00"),
             ("cruise_offsets", "Cruise offsets, not explained", "#56B4E9"),
             ("within_cruise", "Within-cruise, not explained", "#0072B2"),
             ("within_tile", "Within-tile variability", "#BBBBBB"),
             ("technical", "Technical (platform, depth)", "#666666")]  # fmt: skip
    fig, ax = plt.subplots(figsize=(7.2, 3.3))
    x = np.arange(len(ORDER))
    bottom = np.zeros(len(ORDER))
    for col, lab, colr in parts:
        v = p.loc[ORDER, col].clip(lower=0).to_numpy() * 100
        ax.bar(x, v, 0.62, bottom=bottom, color=colr, label=lab, edgecolor="white", linewidth=0.8)
        if col == "explained":
            for xi, vi in zip(x, v):
                ax.text(xi, vi / 2, f"{vi:.0f}", ha="center", va="center", fontsize=7, color="black")
        bottom += v
    for xi, t in zip(x, ORDER):
        ax.text(xi, 101.5, f"{p.loc[t, 'of_nontechnical'] * 100:.0f}% | {p.loc[t, 'of_reproducible_nontechnical'] * 100:.0f}%",
                ha="center", va="bottom", fontsize=6.5)  # fmt: skip
    ax.axvline(1.5, color="#999999", lw=0.6, ls=(0, (2, 2)))
    ax.set_xticks(x)
    ax.set_xticklabels([tl[t] for t in ORDER], fontsize=6.5, linespacing=1.05)
    ax.set_ylim(0, 108)
    ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel("Share of total variance (%)")
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[::-1], l[::-1], frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    for suf in ("png", "pdf"):
        fig.savefig(FIG / f"variance_partition.{suf}", dpi=300)
    plt.close(fig)
    c_figure(plt)


def c_figure(plt):
    """Eukaryote : prokaryote reads (C): a observed vs held-out predicted (Granite PCs), b per-cruise skill of Granite PCs,
    OLCI + SST and time & place, c observed and predicted against latitude."""
    m, p = MAIN["C"]
    z = {k: np.load(OUT / f"predictions_C_{k}_{p}.npz", allow_pickle=True) for k in ("granite", "sat", "place")}
    obs = z[m]["observed"][:, 0].astype(float)
    cr = z[m]["cruise"]
    pr = {k: v["predicted"][:, 0].astype(float) for k, v in z.items()}
    unit = pd.read_csv(OUT / "c_tiles.csv")
    lat = unit["lat"].to_numpy()
    plat = np.isin(cr, NOVASEQ_6000)
    BLUE, ORANGE, GREY = "#0072B2", "#E69F00", "#999999"
    fig = plt.figure(figsize=(7.2, 2.6))
    a = fig.add_axes([0.065, 0.17, 0.25, 0.75])
    b = fig.add_axes([0.40, 0.17, 0.25, 0.75])
    c = fig.add_axes([0.73, 0.17, 0.25, 0.75])
    lo, hi = np.r_[obs, pr[m]].min() - 0.1, np.r_[obs, pr[m]].max() + 0.1
    for msk, col, lab in ((plat, ORANGE, "NovaSeq 6000 cruises"), (~plat, BLUE, "NovaSeq X cruises")):
        a.scatter(pr[m][msk], obs[msk], s=4, lw=0, alpha=0.6, color=col, label=lab)
    a.plot([lo, hi], [lo, hi], color="#555555", lw=0.6, ls=(0, (3, 2)))
    a.set_xlim(lo, hi)
    a.set_ylim(lo, hi)
    a.set_xlabel("Predicted, cruise held out")
    a.set_ylabel("Observed log(eukaryote / prokaryote reads)")
    a.legend(frameon=False, fontsize=6, loc="upper left", markerscale=2, handletextpad=0.2)
    e = 1 - ((obs - pr[m]) ** 2).sum() / ((obs - obs.mean()) ** 2).sum()
    a.text(0.97, 0.04, f"R² = {e:.2f}\n{len(obs):,} tiles", transform=a.transAxes, ha="right", va="bottom", fontsize=6.5)
    cruises = pd.Series(obs).groupby(cr).size().index
    rows = []
    for k in cruises:
        te = cr == k
        sst = ((obs[te] - obs[~te].mean()) ** 2).sum()
        rows.append({kk: 1 - ((obs[te] - pr[kk][te]) ** 2).sum() / sst for kk in pr} | {"cruise": k, "plat": bool(plat[te][0])})
    pc = pd.DataFrame(rows).sort_values("granite")
    yy = np.arange(len(pc))
    for kk, col, mk, lab in (("place", GREY, "s", "Time & place"), ("sat", "#56B4E9", "^", "OLCI + SST"),
                             ("granite", BLUE, "o", "Granite PCs")):  # fmt: skip
        b.scatter(pc[kk].clip(lower=-1), yy, s=16 if kk == "granite" else 10, marker=mk, color=col, lw=0, label=lab, zorder=3)
    for i, (_, r) in enumerate(pc.iterrows()):
        v = [max(r[kk], -1) for kk in ("place", "sat", "granite")]
        b.plot([min(v), max(v)], [i, i], color="#dddddd", lw=0.6, zorder=1)
    b.axvline(0, color="#cccccc", lw=0.6)
    b.set_yticks(yy)
    b.set_yticklabels([("* " if r else "") + k for k, r in zip(pc["cruise"], pc["plat"])], fontsize=5.5)
    b.set_xlabel("R² within the held-out cruise")
    b.set_xlim(-1.05, 1)
    b.legend(frameon=False, fontsize=5.5, loc="center left", handletextpad=0.1)
    pc.round(4).to_csv(OUT / "c_per_cruise.csv", index=False)
    o = np.argsort(lat)
    c.scatter(lat, obs, s=3, lw=0, color="#bbbbbb", label="Observed")
    sm = pd.Series(pr[m][o]).rolling(41, center=True, min_periods=10).mean()
    so = pd.Series(obs[o]).rolling(41, center=True, min_periods=10).mean()
    c.plot(lat[o], so, color="#333333", lw=1.0, label="Observed, running mean")
    c.plot(lat[o], sm, color=BLUE, lw=1.0, label="Predicted, running mean")
    c.set_xlabel("Latitude (°)")
    c.set_ylabel("log(eukaryote / prokaryote reads)")
    c.legend(frameon=False, fontsize=5.5, loc="upper left", markerscale=2, handletextpad=0.2)
    for ax, letter in zip((a, b, c), "abc"):
        ax.text(-0.3 if letter != "b" else -0.2, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for suf in ("png", "pdf"):
        fig.savefig(FIG / f"c_prediction.{suf}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", default=["fit", "partition", "figure"], choices=["fit", "partition", "figure"])
    ap.add_argument("--ko-penalty", default="cv", choices=PENALTIES, help="penalty of the main KO model")
    ap.add_argument("--c-penalty", default="cv", choices=["gcv", "cv"], help="penalty of the main C model")
    a = ap.parse_args()
    MAIN["ko_penalty"], MAIN["C"] = a.ko_penalty, ("granite", a.c_penalty)
    OUT.mkdir(parents=True, exist_ok=True)
    fs = features() if {"fit", "partition"} & set(a.parts) else None
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    if "fit" in a.parts:
        fit_part(fs, ks)
        d = setup("C", fs, ks)
        d["unit"][["section", "sample_lat", "sample_lon"]].rename(columns={"sample_lat": "lat", "sample_lon": "lon"}).assign(
            k=d["k"]).to_csv(OUT / "c_tiles.csv")  # fmt: skip
    if "partition" in a.parts:
        say(f"\nVariance partition (main models: KO targets granite/{MAIN['ko_penalty']}, C granite/{MAIN['C'][1]}):")
        partition_part(fs, ks)
    if "figure" in a.parts:
        figure_part()
    if LOG:
        (OUT / f"summary_{'_'.join(p for p in a.parts if p != 'figure')}.txt").write_text("\n".join(LOG) + "\n")
