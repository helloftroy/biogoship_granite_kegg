#!/usr/bin/env python3
"""How good can the KO models get? Reductions of the embedding, other learners, and in-situ metadata.

KOs: 5 from each category of the all-KO comparison with in-situ temperature (R6,
data/results/ko_temperature_all/ko_temperature_gap.csv), evenly spread in rank:
  clear     embedding clearly better than temperature (>= 0.1, >= 9 of 13 cruises),
            embedding R2 >= 0.2
  ahead     embedding ahead by >= 0.1 pooled, but not in 9 of 13 cruises;
            embedding R2 >= 0.2
  similar   |difference| < 0.03, both R2 >= 0.3 (mostly the cyanobacterial KOs)
  worse     temperature clearly better, temperature R2 >= 0.2
(Chosen from the same cross-validation, so the absolute scores of the embedding
ridge baseline are flattering within "clear" and "ahead"; comparisons between
methods on the same KOs are fair.)

Feature blocks (all standardised on the training cruises; PCA and PLS are fitted
on the training cruises only, so the held-out cruise never shapes the features):
  emb        Granite embedding, 512 numbers
  pcaK       first K principal components of the embedding (K = 4, 8, 16, 32, 64)
  plsK       K partial-least-squares components of the embedding against the KOs
             being predicted (fitted to the training cruises' KOs)
  meta       position and season (lat/90, sine and cosine of longitude and day of
             year), in-situ temperature and salinity
  temp       in-situ temperature
  ts         in-situ temperature and salinity
  season     sine and cosine of day of year
  chl        band-ratio chlorophyll index log10(max(Oa03-05) / Oa06) (sample set col)
  bands      the 16 OLCI band means over clear pixels (sample set col)
  np         in-situ nitrate + nitrite and phosphate (only some samples: set "np")

Learners:
  ridge      RidgeCV, alpha per KO
  spline     cubic B-spline of each input (8 knots), then ridge (temperature alone
             only: on 6+ inputs it extrapolates wildly to new cruises)
  krr        kernel ridge, RBF kernel; alpha and gamma chosen by an inner
             leave-cruises-out split (3 groups) of the training cruises
  mlp        one hidden layer of 64, L2 penalty chosen like krr, early stopping
  catboost   500 trees, depth 6, learning rate 0.05, MultiRMSE

Sample sets:
  all   samples with a KO profile, an in-situ temperature and salinity, and a
        fixed tile with SST (option A); identical tiles merged
  col   the subset whose tile has band means (tile_raw_features.csv in both
        embeddings folders)
  np    the subset that also has nitrate + nitrite and phosphate (FAIRe sheet
        values where the manifest has none)

Outputs in data/results/model_experiments/: scores.csv (per KO, sample set,
experiment: pooled R2 and median R2 within a held-out cruise) and per_cruise.csv.

Usage:
    python scripts/model_experiments.py
    python scripts/model_experiments.py --only pca32+meta:krr emb:ridge
"""

from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import RidgeCV
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship.tables import REPO_ROOT  # noqa: E402
from screen_ko_predictability import ALPHAS, load_embeddings, metadata_features  # noqa: E402
from test_ko_temperature import BAND_MEANS, chlorophyll_index, load_raw_features  # noqa: E402

# (sample set, feature blocks, learner)
EXPERIMENTS = [
    ("all", "emb", "ridge"),
    *[("all", f"pca{k}", learner) for k in (16, 32, 64) for learner in ("ridge", "krr", "mlp")],
    ("all", "pca32", "catboost"),
    ("all", "pls32", "ridge"), ("all", "pls32", "krr"), ("all", "pls32", "catboost"),
    ("all", "temp", "spline"),
    ("all", "meta", "ridge"), ("all", "meta", "krr"), ("all", "meta", "catboost"),
    ("all", "temp+pca4", "krr"), ("all", "temp+pca8", "krr"), ("all", "temp+pca8", "catboost"),
    ("all", "emb+meta", "ridge"),
    ("all", "pca32+meta", "ridge"), ("all", "pca32+meta", "krr"), ("all", "pca32+meta", "mlp"), ("all", "pca32+meta", "catboost"),
    ("all", "pls32+meta", "krr"), ("all", "pls32+meta", "catboost"),
    # metadata without position (ts+season) and without position and season (ts), ridge only (R9)
    *[("all", f"{base}{extra}", "ridge") for base in ("emb", "pca16", "pca32")
      for extra in ("+ts", "+ts+season", "+temp", "+temp+season")],
    ("all", "ts", "ridge"), ("all", "ts+season", "ridge"), ("all", "pca16+meta", "ridge"),
    # temperature + chlorophyll index (+ PCs) on the tiles with band means (R10)
    *[("col", f, "ridge") for f in ("emb", "pca16", "pca32", "pca16+temp", "pca32+temp", "chl", "temp+chl", "bands", "temp+bands",
                                    "emb+chl", "pca16+temp+chl", "pca32+temp+chl", "pca32+chl", "temp+pca4+chl", "temp+pca8+chl")],
    ("col", "temp", "spline"), ("col", "temp+chl", "spline"), ("col", "temp+pca8", "catboost"), ("col", "temp+pca8+chl", "catboost"),
    ("col", "temp+pca4+chl", "catboost"),
    ("np", "emb", "ridge"), ("np", "temp", "spline"), ("np", "meta", "catboost"),
    ("np", "meta+np", "ridge"), ("np", "meta+np", "krr"), ("np", "meta+np", "catboost"),
    ("np", "emb+meta+np", "ridge"), ("np", "pca32+meta+np", "krr"), ("np", "pca32+meta+np", "catboost"),
]  # fmt: skip
KRR_ALPHAS = (0.03, 0.1, 0.3, 1.0, 3.0)
KRR_GAMMA_SCALES = (0.3, 1.0, 3.0)  # times 1 / number of features
MLP_ALPHAS = (0.1, 1.0, 10.0)


def main() -> None:
    args = parse_args()
    experiments = [e for e in EXPERIMENTS if not args.only or f"{e[0]}:{e[1]}:{e[2]}" in args.only
                   or f"{e[1]}:{e[2]}" in args.only]  # fmt: skip
    kos = pick_kos(args)
    data = load(args, kos)
    rng = np.random.default_rng(args.seed)
    out, per_cruise = [], []
    for sample_set, blocks, learner in experiments:
        d = data[sample_set]
        y, cruises = d["Y"].to_numpy(), d["cruises"]
        pred = np.full_like(y, np.nan)
        started = time.time()
        for cruise in np.unique(cruises):
            test = cruises == cruise
            train_x, test_x = features(blocks, d, ~test, test, y[~test])
            pred[test] = fit_predict(learner, train_x, y[~test], test_x, cruises[~test], args.seed + int(rng.integers(1000)))
            ss_res = ((y[test] - pred[test]) ** 2).sum(axis=0)
            ss_tot = ((y[test] - y[~test].mean(axis=0)) ** 2).sum(axis=0)
            per_cruise.append(pd.DataFrame({"ko": d["Y"].columns, "sample_set": sample_set, "features": blocks, "learner": learner,
                                            "held_out": cruise, "n": int(test.sum()), "r2": 1 - ss_res / ss_tot}))  # fmt: skip
        pooled = 1 - ((y - pred) ** 2).sum(axis=0) / ((y - y.mean(axis=0)) ** 2).sum(axis=0)
        out.append(pd.DataFrame({"ko": d["Y"].columns, "sample_set": sample_set, "features": blocks, "learner": learner,
                                 "r2_pooled": pooled}))  # fmt: skip
        print(f"  {sample_set:3s} {blocks:15s} {learner:8s} median pooled R2 {np.median(pooled):6.3f}  ({time.time() - started:4.0f} s)",
              flush=True)  # fmt: skip

    per_cruise = pd.concat(per_cruise)
    scores = pd.concat(out).merge(per_cruise.groupby(["ko", "sample_set", "features", "learner"])["r2"].median()
                                  .rename("r2_median_cruise").reset_index())  # fmt: skip
    scores["category"] = scores["ko"].map(kos.set_index("ko")["category"])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scores.round(4).to_csv(args.out_dir / "scores.csv", index=False)
    per_cruise.round(4).to_csv(args.out_dir / "per_cruise.csv", index=False)
    pd.set_option("display.width", 250)
    for value, label in (("r2_pooled", "pooled R2"), ("r2_median_cruise", "median R2 within a held-out cruise")):
        table = scores.pivot_table(index=["sample_set", "features", "learner"], columns="category", values=value, aggfunc="median")
        table["all 20"] = scores.groupby(["sample_set", "features", "learner"])[value].median()
        print(f"\nMedian over KOs, {label}:")
        print(table[["clear", "ahead", "similar", "worse", "all 20"]].round(2).to_string())
    print(f"\nWrote {args.out_dir / 'scores.csv'}")


def load(args: argparse.Namespace, kos: pd.DataFrame) -> dict:
    targets = pd.read_csv(args.targets, index_col="sample_key", usecols=["sample_key", *kos["ko"]])[list(kos["ko"])]
    rows, X_embed = load_embeddings(args)
    rows = rows[rows["sample_id"].isin(targets.index)].copy()
    extra = in_situ(args)
    rows["salinity"] = rows["sample_id"].map(extra["salinity"])
    rows["nitrate_nitrite"] = pd.to_numeric(rows["nitrate_nitrite_umol_l"], errors="coerce").fillna(rows["sample_id"].map(extra["nitrate_nitrite"]))
    rows["phosphate"] = pd.to_numeric(rows["phosphate_umol_kg"], errors="coerce").fillna(rows["sample_id"].map(extra["phosphate"]))
    rows = rows[rows["temperature_c"].notna() & rows["salinity"].notna()]
    raw = load_raw_features(args.embeddings_dirs)
    has_colour = rows["tile_path"].isin(raw.index)
    data = {}
    for name, keep in (("all", rows.index), ("np", rows.index[rows[["nitrate_nitrite", "phosphate"]].notna().all(axis=1)]),
                       ("col", rows.index[has_colour])):  # fmt: skip
        part = rows.loc[keep]
        content = pd.Series([v.tobytes() for v in X_embed[part["embedding_row"].to_numpy()].round(5)], index=part.index)
        part = part.assign(tile_key=content.map({k: i for i, k in enumerate(content.unique())}))
        Y = targets.loc[part["sample_id"]].set_axis(part.index).groupby(part["tile_key"]).mean()
        means = part.groupby("tile_key")[["temperature_c", "salinity", "nitrate_nitrite", "phosphate"]].mean().loc[Y.index]
        unit = part.drop_duplicates("tile_key").set_index("tile_key").loc[Y.index]
        data[name] = {"Y": Y, "cruises": unit["section"].astype(str).to_numpy(),
                      "emb": X_embed[unit["embedding_row"].to_numpy()],
                      "temp": means[["temperature_c"]].to_numpy(),
                      "meta": np.hstack([metadata_features(unit), means[["temperature_c", "salinity"]].to_numpy()]),
                      "ts": means[["temperature_c", "salinity"]].to_numpy(),
                      "season": metadata_features(unit)[:, 3:5],
                      "np": means[["nitrate_nitrite", "phosphate"]].to_numpy()}  # fmt: skip
        if name == "col":
            bands = raw.loc[unit["tile_path"], BAND_MEANS].to_numpy()
            data[name].update({"bands": bands, "chl": chlorophyll_index(bands)[:, None]})
        print(f"Sample set {name}: {len(part):,} samples, {len(unit):,} tiles, {unit['section'].nunique()} cruises")
    return data


def in_situ(args: argparse.Namespace) -> pd.DataFrame:
    """Salinity, nitrate + nitrite and phosphate from the FAIRe sheet, keyed by FAIRe samp_name and original sample id."""
    faire = pd.read_csv(args.faire_samples)
    nn = pd.to_numeric(faire["nitrate_plus_nitrite"], errors="coerce")
    nn = nn.fillna(pd.to_numeric(faire["nitrate"], errors="coerce") + pd.to_numeric(faire["nitrite"], errors="coerce").fillna(0))
    table = pd.DataFrame({"salinity": pd.to_numeric(faire["salinity"], errors="coerce"), "nitrate_nitrite": nn,
                          "phosphate": pd.to_numeric(faire["phosphate"], errors="coerce")})  # fmt: skip
    by_name = table.set_index(faire["samp_name"])
    by_original = table[faire["original_sample_id"].notna().to_numpy()].set_index(faire["original_sample_id"].dropna())
    both = pd.concat([by_name, by_original])
    return both[~both.index.duplicated()]


def features(blocks: str, d: dict, train: np.ndarray, test: np.ndarray, y_train: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    parts_train, parts_test = [], []
    for block in blocks.split("+"):
        if block.startswith(("pca", "pls")):
            k = int(block[3:])
            scaler = StandardScaler().fit(d["emb"][train])
            e_train, e_test = scaler.transform(d["emb"][train]), scaler.transform(d["emb"][test])
            if block.startswith("pca"):
                model = PCA(n_components=k, random_state=0).fit(e_train)
                a, b = model.transform(e_train), model.transform(e_test)
            else:
                model = PLSRegression(n_components=k, scale=True).fit(e_train, y_train)
                a, b = model.transform(e_train), model.transform(e_test)
        else:
            a, b = d[block][train], d[block][test]
        scaler = StandardScaler().fit(a)
        parts_train.append(scaler.transform(a))
        parts_test.append(scaler.transform(b))
    return np.hstack(parts_train), np.hstack(parts_test)


def fit_predict(learner: str, x: np.ndarray, y: np.ndarray, x_test: np.ndarray, groups: np.ndarray, seed: int) -> np.ndarray:
    if learner == "ridge":
        return RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(x, y).predict(x_test)
    if learner == "spline":
        model = make_pipeline(SplineTransformer(n_knots=8, degree=3, extrapolation="linear"), RidgeCV(alphas=ALPHAS, alpha_per_target=True))
        return model.fit(x, y).predict(x_test)
    if learner == "krr":
        params = [{"alpha": a, "gamma": g / x.shape[1]} for a, g in itertools.product(KRR_ALPHAS, KRR_GAMMA_SCALES)]
        best = inner_select(lambda p: KernelRidge(kernel="rbf", **p), params, x, y, groups)
        return KernelRidge(kernel="rbf", **best).fit(x, y).predict(x_test)
    if learner == "mlp":
        make = lambda p: MLPRegressor(hidden_layer_sizes=(64,), alpha=p["alpha"], early_stopping=True, max_iter=1000,  # noqa: E731
                                      random_state=seed)  # fmt: skip
        best = inner_select(make, [{"alpha": a} for a in MLP_ALPHAS], x, y, groups)
        return make(best).fit(x, y).predict(x_test)
    from catboost import CatBoostRegressor

    model = CatBoostRegressor(iterations=500, depth=6, learning_rate=0.05, loss_function="MultiRMSE",
                              random_seed=seed, verbose=False, thread_count=-1)  # fmt: skip
    return model.fit(x, y).predict(x_test)


def inner_select(make, params: list[dict], x: np.ndarray, y: np.ndarray, groups: np.ndarray) -> dict:
    """Parameters with the lowest squared error over 3 inner splits of the training cruises (whole cruises held out)."""
    names = np.unique(groups)
    split = np.random.default_rng(0).permutation(len(names)) % 3
    folds = [np.isin(groups, names[split == k]) for k in range(3)]
    errors = []
    for p in params:
        err = 0.0
        for held in folds:
            if held.all() or not held.any():
                continue
            err += ((make(p).fit(x[~held], y[~held]).predict(x[held]) - y[held]) ** 2).sum()
        errors.append(err)
    return params[int(np.argmin(errors))]


def pick_kos(args: argparse.Namespace) -> pd.DataFrame:
    gap = pd.read_csv(args.gap, index_col="ko")
    clear_emb = gap["verdict"] == "embedding clearly better"
    categories = {
        "clear": gap[clear_emb & (gap["r2_embedding"] >= 0.2)].sort_values("r2_embedding", ascending=False),
        "ahead": gap[~clear_emb & (gap["difference"] >= 0.1) & (gap["r2_embedding"] >= 0.2)].sort_values("r2_embedding", ascending=False),
        "similar": gap[(gap["difference"].abs() < 0.03) & (gap[["r2_temp", "r2_embedding"]].min(axis=1) >= 0.3)]
        .sort_values("r2_embedding", ascending=False),
        "worse": gap[(gap["verdict"] == "temperature clearly better") & (gap["r2_temp"] >= 0.2)].sort_values("r2_temp", ascending=False),
    }
    frames = []
    for name, table in categories.items():
        pick = table.iloc[np.linspace(0, len(table) - 1, args.per_category).round().astype(int)]
        frames.append(pick.assign(category=name).reset_index()[["ko", "category", "r2_temp", "r2_embedding"]])
    kos = pd.concat(frames, ignore_index=True)
    print(kos.round(2).to_string(index=False))
    return kos


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--gap", type=Path, default=REPO_ROOT / "data" / "results" / "ko_temperature_all" / "ko_temperature_gap.csv")
    parser.add_argument("--faire-samples", type=Path, default=REPO_ROOT / "data" / "faire" / "faire_samples.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings", REPO_ROOT / "data" / "faire" / "embeddings"])  # fmt: skip
    parser.add_argument("--position", choices=["fixed", "advected"], default="fixed")
    parser.add_argument("--best-tile", choices=["prefer-sst", "closest"], default="prefer-sst")
    parser.add_argument("--per-category", type=int, default=5)
    parser.add_argument("--only", nargs="+", help="Run only these experiments: 'features:learner' or 'set:features:learner'.")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "results" / "model_experiments")
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    main()
