#!/usr/bin/env python3
"""Model development on the KOs in the two figures (heatmap rows; the cruise-line genes are among them).

Same samples, targets and cruise-held-out CV as scripts/test_ko_temperature.py (R17):
1,429 samples / 1,216 distinct tiles / 13 cruises, OSTIA-filled embeddings.

Per fold, on the training cruises only:
  E      the Granite embedding (512), standardised
  P<k>   the first k principal components of E, standardised
  T      in-situ temperature, standardised;  T~ its cubic B-spline basis (8 knots)
  S      raw satellite: 16 OLCI band means + tile SST mean, rank-normal, median-imputed
         (as R17);  S~ a spline basis of each;  SST~ / B~ the SST mean or the 16 bands alone
  ~P<k>  a spline basis of each of the k PCs (an additive curve per PC)

Learners: ridge with the penalty chosen per KO by leave-one-out (GCV, as before), or
"grp": chosen per KO by leave-one-cruise-out within the training cruises. GCV holds out
single samples, whose neighbours on the same cruise stay in training, so it tends to pick
a penalty that is too weak for a new cruise. "avg(a, b)" averages two models' predictions.

Per-KO choice ("best per KO") is nested: for each held-out cruise, every model is scored
by cruise-held-out CV on the other 12 cruises only, the best is refitted on those 12 and
predicts the held-out cruise. So the reported R2 of the choice is not chosen on itself.

Outputs in --out-dir (default data/results/figure_ko_models/):
  model_r2.csv        KO x model: pooled R2 and median R2 within a held-out cruise
  model_r2_cruise.csv KO x model x held-out cruise: R2 (against the training-cruise mean)
  best_per_ko.csv     per KO and held-out cruise: the model chosen
  predictions.npz     held-out predictions of every model and of the per-KO choice

Usage:
    python scripts/tune_figure_kos.py
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import QuantileTransformer, SplineTransformer, StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship.tables import REPO_ROOT  # noqa: E402
from screen_ko_predictability import ALPHAS, load_embeddings  # noqa: E402
from test_ko_temperature import satellite_features  # noqa: E402

RESULTS = REPO_ROOT / "data" / "results"

# name: (blocks, learner). Blocks are built per fold by build().
MODELS = {
    "emb": (["E"], "gcv"),
    "emb grp": (["E"], "grp"),
    "temp~": (["T~"], "gcv"),
    "sat~": (["S~"], "gcv"),
    "sst~": (["SST~"], "gcv"),
    "bands~": (["B~"], "gcv"),
    "temp~ + sat~": (["T~", "S~"], "gcv"),
    "pca8 + temp": (["P8", "T"], "gcv"),
    "pca16 + temp": (["P16", "T"], "gcv"),
    "pca32 + temp": (["P32", "T"], "gcv"),
    "pca64 + temp": (["P64", "T"], "gcv"),
    "pca8 + temp~": (["P8", "T~"], "gcv"),
    "pca16 + temp~": (["P16", "T~"], "gcv"),
    "pca32 + temp~": (["P32", "T~"], "gcv"),
    "pca32 + temp~ grp": (["P32", "T~"], "grp"),
    "emb + temp~ grp": (["E", "T~"], "grp"),
    "~pca8 + temp~": (["~P8", "T~"], "gcv"),
    "~pca16 + temp~": (["~P16", "T~"], "gcv"),
    "pca16 + temp~ + sat~": (["P16", "T~", "S~"], "gcv"),
    "~pca8 + temp~ + sat~": (["~P8", "T~", "S~"], "gcv"),
}
AVERAGES = {
    "avg(emb, temp~)": ("emb", "temp~"),
    "avg(pca32 + temp~, sat~)": ("pca32 + temp~", "sat~"),
    "avg(emb, temp~ + sat~)": ("emb", "temp~ + sat~"),
}


def main() -> None:
    args = parse_args()
    kos = pd.read_csv(args.kos)["ko"].tolist()
    emb, temp, sat, y, cruises, unit = load(args, kos)
    global CRUISES
    CRUISES = cruises
    names = list(MODELS) + list(AVERAGES)
    print(f"{len(unit):,} tiles, {len(set(cruises))} cruises, {len(kos)} KOs, {len(names)} models", flush=True)

    t0 = time.time()
    outer = cv_predict(np.arange(len(y)), emb, temp, sat, y, cruises)  # model -> held-out predictions
    print(f"Outer CV done ({time.time() - t0:.0f} s); nested choice per KO ...", flush=True)

    chosen_pred = np.full_like(y, np.nan)
    chosen = []
    for cruise in np.unique(cruises):
        train = np.flatnonzero(cruises != cruise)
        test = np.flatnonzero(cruises == cruise)
        inner = cv_predict(train, emb, temp, sat, y, cruises)
        inner_r2 = np.array([pooled_r2(y[train], inner[m][train]) for m in names])  # models x KOs
        pick = inner_r2.argmax(axis=0)
        fitted = fit_all(train, test, emb, temp, sat, y, needed={names[i] for i in set(pick)})
        for j, i in enumerate(pick):
            chosen_pred[test, j] = fitted[names[i]][:, j]
        chosen.append(pd.DataFrame({"ko": kos, "held_out": cruise, "model": [names[i] for i in pick],
                                    "inner_r2": inner_r2[pick, np.arange(len(kos))]}))  # fmt: skip
        print(f"  {cruise:14s} ({time.time() - t0:.0f} s)", flush=True)
    outer["best per KO"] = chosen_pred
    names.append("best per KO")

    rows, cruise_rows = [], []
    for m in names:
        rows.append(pd.DataFrame({"ko": kos, "model": m, "r2_pooled": pooled_r2(y, outer[m])}))
        for cruise in np.unique(cruises):
            test = cruises == cruise
            ss_res = ((y[test] - outer[m][test]) ** 2).sum(axis=0)
            ss_tot = ((y[test] - y[~test].mean(axis=0)) ** 2).sum(axis=0)
            cruise_rows.append(pd.DataFrame({"ko": kos, "model": m, "held_out": cruise, "r2": 1 - ss_res / ss_tot}))
    per_cruise = pd.concat(cruise_rows)
    table = pd.concat(rows).merge(per_cruise.groupby(["ko", "model"])["r2"].median().rename("r2_median_cruise").reset_index())

    args.out_dir.mkdir(parents=True, exist_ok=True)
    table.round(4).to_csv(args.out_dir / "model_r2.csv", index=False)
    per_cruise.round(4).to_csv(args.out_dir / "model_r2_cruise.csv", index=False)
    pd.concat(chosen).round(4).to_csv(args.out_dir / "best_per_ko.csv", index=False)
    np.savez_compressed(args.out_dir / "predictions.npz", observed=y.astype(np.float32), ko=np.array(kos),
                        models=np.array(names), predicted=np.stack([outer[m] for m in names]).astype(np.float32),
                        tile_path=unit["tile_path"].to_numpy().astype(str), section=cruises)  # fmt: skip

    groups = pd.read_csv(args.kos).set_index("ko")["category"]
    wide = table.pivot(index="ko", columns="model", values="r2_pooled")[names]
    within = table.pivot(index="ko", columns="model", values="r2_median_cruise")[names]
    summary = pd.concat({"all": wide.median(), **{g: wide[groups.reindex(wide.index) == g].median()
                                                   for g in ["granite", "same", "temp"]},
                         "within": within.median()}, axis=1)  # fmt: skip
    pd.set_option("display.width", 200)
    print("\nMedian pooled R2 over the figure KOs (category from the heatmap), and median R2 within a held-out cruise:")
    print(summary.round(3).sort_values("all", ascending=False).to_string())
    print("\nModel chosen per KO (count over KO x held-out cruise):")
    print(pd.concat(chosen)["model"].value_counts().to_string())
    print(f"\nWrote {args.out_dir}")


def load(args: argparse.Namespace, kos: list[str]):
    """The R17 samples: tiles merged by content, targets averaged over merged samples."""
    targets = pd.read_csv(args.targets, index_col="sample_key", usecols=["sample_key", *kos])[kos]
    ns = argparse.Namespace(embeddings_dirs=args.embeddings_dirs, position="fixed", best_tile="prefer-sst")
    rows, X_embed = load_embeddings(ns)
    rows = rows[rows["sample_id"].isin(targets.index) & rows["temperature_c"].notna()]
    content = pd.Series([v.tobytes() for v in X_embed[rows["embedding_row"].to_numpy()].round(5)], index=rows.index)
    rows = rows.assign(tile_key=content.map({k: i for i, k in enumerate(content.unique())}))
    Y = targets.loc[rows["sample_id"]].set_axis(rows.index).groupby(rows["tile_key"]).mean()
    temp = rows.groupby("tile_key")["temperature_c"].mean().loc[Y.index].to_numpy()[:, None]
    unit = rows.drop_duplicates("tile_key").set_index("tile_key").loc[Y.index]
    emb = X_embed[unit["embedding_row"].to_numpy()]
    sat = satellite_features(args.embeddings_dirs, unit["tile_path"])
    return emb, temp, sat, Y.to_numpy(), unit["section"].astype(str).to_numpy(), unit


def build(blocks: list[str], train: np.ndarray, test: np.ndarray, emb, temp, sat, cache: dict):
    """Feature blocks fitted on the training rows; cached per fold."""
    out_tr, out_te = [], []
    for b in blocks:
        if b not in cache:
            cache[b] = make_block(b, train, test, emb, temp, sat, cache)
        out_tr.append(cache[b][0])
        out_te.append(cache[b][1])
    return np.hstack(out_tr), np.hstack(out_te)


def make_block(b: str, train, test, emb, temp, sat, cache):
    def spline(tr, te):
        s = SplineTransformer(n_knots=8, degree=3, extrapolation="linear").fit(tr)
        return s.transform(tr), s.transform(te)

    def scale(tr, te):
        s = StandardScaler().fit(tr)
        return s.transform(tr), s.transform(te)

    if b == "E":
        return scale(emb[train], emb[test])
    if b == "T":
        return scale(temp[train], temp[test])
    if b == "T~":
        return spline(*cache.setdefault("T", scale(temp[train], temp[test])))
    if b.startswith("P"):
        e_tr, e_te = cache.setdefault("E", scale(emb[train], emb[test]))
        if "pca" not in cache:
            cache["pca"] = PCA(n_components=64, random_state=0).fit(e_tr)
        k = int(b[1:])
        return scale(cache["pca"].transform(e_tr)[:, :k], cache["pca"].transform(e_te)[:, :k])
    if b.startswith("~P"):
        return spline(*cache.setdefault(b[1:], make_block(b[1:], train, test, emb, temp, sat, cache)))
    if b in ("S~", "SST~", "B~"):
        if "S" not in cache:
            tr, te = sat[train], sat[test]
            med = np.nanmedian(tr, axis=0)
            tr, te = np.where(np.isnan(tr), med, tr), np.where(np.isnan(te), med, te)
            q = QuantileTransformer(n_quantiles=min(1000, len(tr)), output_distribution="normal", random_state=0).fit(tr)
            cache["S"] = scale(q.transform(tr), q.transform(te))
        tr, te = cache["S"]
        cols = {"S~": slice(None), "SST~": slice(-1, None), "B~": slice(0, -1)}[b]
        return spline(tr[:, cols], te[:, cols])
    raise ValueError(b)


def ridge_gcv(xtr, ytr, xte):
    return RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(xtr, ytr).predict(xte)


def ridge_grouped(xtr, ytr, xte, groups):
    """Ridge with the penalty per KO chosen by leave-one-cruise-out within the training rows."""
    sse = np.zeros((len(ALPHAS), ytr.shape[1]))
    for g in np.unique(groups):
        a, b = groups != g, groups == g
        sse += ((ytr[b][None] - ridge_path(xtr[a], ytr[a], xtr[b])) ** 2).sum(axis=1)
    best = ALPHAS[sse.argmin(axis=0)]
    path = ridge_path(xtr, ytr, xte)
    return path[np.searchsorted(ALPHAS, best), :, np.arange(ytr.shape[1])].T


def ridge_path(xtr, ytr, xte):
    """Predictions for every alpha in ALPHAS: (alphas, test rows, KOs)."""
    xm, ym = xtr.mean(axis=0), ytr.mean(axis=0)
    u, s, vt = np.linalg.svd(xtr - xm, full_matrices=False)
    uty = u.T @ (ytr - ym)
    proj = (xte - xm) @ vt.T
    shrink = s[None, :] / (s[None, :] ** 2 + ALPHAS[:, None])  # alphas x components
    return ym + np.einsum("tc,ac,ck->atk", proj, shrink, uty)


def fit_all(train, test, emb, temp, sat, y, needed=None):
    """Fit every model (or those needed, with their averaged parts) on train, predict test."""
    needed = set(MODELS) | set(AVERAGES) if needed is None else needed
    base = {m for m in needed if m in MODELS} | {p for m in needed if m in AVERAGES for p in AVERAGES[m]}
    cache, out = {}, {}
    groups = None
    for m in MODELS:
        if m not in base:
            continue
        blocks, learner = MODELS[m]
        xtr, xte = build(blocks, train, test, emb, temp, sat, cache)
        if learner == "gcv":
            out[m] = ridge_gcv(xtr, y[train], xte)
        else:
            groups = CRUISES[train] if groups is None else groups
            out[m] = ridge_grouped(xtr, y[train], xte, groups)
    for m, (a, b) in AVERAGES.items():
        if m in needed:
            out[m] = (out[a] + out[b]) / 2
    return out


def cv_predict(rows: np.ndarray, emb, temp, sat, y, cruises):
    """Leave-one-cruise-out over the given rows; full-length arrays, NaN outside rows."""
    out = {m: np.full_like(y, np.nan) for m in list(MODELS) + list(AVERAGES)}
    for cruise in np.unique(cruises[rows]):
        test = rows[cruises[rows] == cruise]
        train = rows[cruises[rows] != cruise]
        for m, pred in fit_all(train, test, emb, temp, sat, y).items():
            out[m][test] = pred
    return out


def pooled_r2(y, pred):
    return 1 - ((y - pred) ** 2).sum(axis=0) / ((y - y.mean(axis=0)) ** 2).sum(axis=0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kos", type=Path, default=RESULTS / "ko_temperature_all_ostia_v2" / "ko_heatmap_rows.csv",
                        help="CSV with columns ko and category (the heatmap rows).")
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings_sst_l4", REPO_ROOT / "data" / "faire" / "embeddings_sst_l4"])
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "figure_ko_models")
    return parser.parse_args()


CRUISES: np.ndarray = np.array([])

if __name__ == "__main__":
    main()
