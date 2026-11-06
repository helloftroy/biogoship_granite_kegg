#!/usr/bin/env python3
"""Does in-situ temperature, added to position + season, predict KOs as well as Granite?

Takes 20 KOs from the screen (data/results/ko_screen/ko_screen.csv): the 10 best
predicted by the embedding, and 10 middling ones (embedding pooled R2 0.25-0.40
and at least 0.15 above position + season), evenly spaced in rank. Then predicts
them, one cruise held out at a time, from:

  metadata          position and season (as in the screen: lat/90, sine and
                    cosine of longitude and of day of year)
  metadata_temp     the same + in-situ temperature (manifest temperature_c)
  temp              in-situ temperature alone
  embedding         Granite embedding of the sample's tile (512 numbers)
  embedding_temp    embedding + in-situ temperature

Learners: ridge, k nearest neighbours (k = 15, distance-weighted) and CatBoost
(500 trees, depth 6, learning rate 0.05; bends freely on few features, which ridge
cannot). The score per KO and feature set is the best of the three learners.

Samples: as in the screen (fixed position by default, prefer-SST tile, SST-only,
identical tiles merged), restricted to samples with an in-situ temperature, so all
feature sets see the same rows. Temperature is averaged over merged samples.

Outputs in data/results/ko_temperature/: ko_temperature.csv (one row per KO and
feature set, each learner) and ko_temperature_per_cruise.csv.

Usage:
    python scripts/test_ko_temperature.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from catboost import CatBoostRegressor
from sklearn.linear_model import RidgeCV
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship.tables import REPO_ROOT  # noqa: E402
from screen_ko_predictability import ALPHAS, load_embeddings, metadata_features  # noqa: E402

FEATURE_SETS = ("metadata", "metadata_temp", "temp", "embedding", "embedding_temp")
LEARNERS = ("ridge", "knn", "catboost")


def main() -> None:
    args = parse_args()
    kos = pick_kos(args)
    targets = pd.read_csv(args.targets, index_col="sample_key", usecols=["sample_key", *kos["ko"]])[list(kos["ko"])]
    rows, X_embed = load_embeddings(args)
    rows = rows[rows["sample_id"].isin(targets.index) & rows["temperature_c"].notna()]
    content = pd.Series([v.tobytes() for v in X_embed[rows["embedding_row"].to_numpy()].round(5)], index=rows.index)
    rows = rows.assign(tile_key=content.map({k: i for i, k in enumerate(content.unique())}))
    Y = targets.loc[rows["sample_id"]].set_axis(rows.index).groupby(rows["tile_key"]).mean()
    temp = rows.groupby("tile_key")["temperature_c"].mean().loc[Y.index].to_numpy()[:, None]
    unit = rows.drop_duplicates("tile_key").set_index("tile_key").loc[Y.index]
    meta, emb = metadata_features(unit), X_embed[unit["embedding_row"].to_numpy()]
    X = {"metadata": meta, "metadata_temp": np.hstack([meta, temp]), "temp": temp,
         "embedding": emb, "embedding_temp": np.hstack([emb, temp])}  # fmt: skip
    cruises = unit["section"].astype(str).to_numpy()
    print(f"Samples with a KO profile, an in-situ temperature and a {args.position} tile with SST: {len(rows):,}; "
          f"distinct tiles: {len(unit):,}; cruises: {len(set(cruises))}")  # fmt: skip

    y = Y.to_numpy()
    results, per_cruise = [], []
    for name in FEATURE_SETS:
        for learner in LEARNERS:
            pred = np.full_like(y, np.nan)
            for cruise in np.unique(cruises):
                test = cruises == cruise
                scaler = StandardScaler().fit(X[name][~test])
                train_x, test_x = scaler.transform(X[name][~test]), scaler.transform(X[name][test])
                pred[test] = fit_predict(learner, train_x, y[~test], test_x, args.seed)
                ss_res = ((y[test] - pred[test]) ** 2).sum(axis=0)
                ss_tot = ((y[test] - y[~test].mean(axis=0)) ** 2).sum(axis=0)
                per_cruise.append(pd.DataFrame({"ko": Y.columns, "feature_set": name, "learner": learner, "held_out": cruise,
                                                "n": int(test.sum()), "r2": 1 - ss_res / ss_tot}))  # fmt: skip
            pooled = 1 - ((y - pred) ** 2).sum(axis=0) / ((y - y.mean(axis=0)) ** 2).sum(axis=0)
            results.append(pd.DataFrame({"ko": Y.columns, "feature_set": name, "learner": learner, "r2_pooled": pooled}))
            print(f"  {name:15s} {learner:9s} median pooled R2 over the 20 KOs {np.median(pooled):6.3f}", flush=True)

    per_cruise = pd.concat(per_cruise)
    out = pd.concat(results).merge(
        per_cruise.groupby(["ko", "feature_set", "learner"])["r2"].median().rename("r2_median_cruise").reset_index())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out.round(4).to_csv(args.out_dir / "ko_temperature.csv", index=False)
    per_cruise.round(4).to_csv(args.out_dir / "ko_temperature_per_cruise.csv", index=False)

    best = out.groupby(["ko", "feature_set"])["r2_pooled"].max().unstack()[list(FEATURE_SETS)]
    best = kos.set_index("ko")[["group", "name"]].join(best)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 38)
    print("\nPooled R2 (best learner per KO and feature set):")
    print(best.round(2).to_string())
    print("\nMedian over KOs:")
    print(best.groupby("group")[list(FEATURE_SETS)].median().round(2).to_string())
    print(f"\nWrote {args.out_dir / 'ko_temperature.csv'}")


def fit_predict(learner: str, train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, seed: int) -> np.ndarray:
    if learner == "ridge":
        return RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(train_x, train_y).predict(test_x)
    if learner == "knn":
        return KNeighborsRegressor(n_neighbors=15, weights="distance").fit(train_x, train_y).predict(test_x)
    model = CatBoostRegressor(iterations=500, depth=6, learning_rate=0.05, loss_function="MultiRMSE",
                              random_seed=seed, verbose=False, thread_count=-1)  # fmt: skip
    return model.fit(train_x, train_y).predict(test_x)


def pick_kos(args: argparse.Namespace) -> pd.DataFrame:
    screen = pd.read_csv(args.screen).sort_values("r2_pooled_embedding", ascending=False)
    top = screen.head(10).assign(group="top")
    margin = screen["r2_pooled_embedding"] - screen["r2_pooled_metadata"]
    middle = screen[screen["r2_pooled_embedding"].between(0.25, 0.40) & (margin >= 0.15)]
    middle = middle.iloc[np.linspace(0, len(middle) - 1, 10).round().astype(int)].assign(group="middle")
    kos = pd.concat([top, middle])
    names = pd.read_csv(args.ko_names, sep="\t", header=None, names=["ko", "name"]).set_index("ko")["name"]
    kos["name"] = kos["ko"].map(names).fillna("").str.replace(r"\s*\[EC:[^\]]*\]", "", regex=True)
    return kos[["ko", "group", "name", "r2_pooled_embedding", "r2_pooled_metadata"]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--screen", type=Path, default=REPO_ROOT / "data" / "results" / "ko_screen" / "ko_screen.csv")
    parser.add_argument("--ko-names", type=Path, default=REPO_ROOT / "data" / "reference" / "kegg_ko_list.tsv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings", REPO_ROOT / "data" / "faire" / "embeddings"])  # fmt: skip
    parser.add_argument("--position", choices=["fixed", "advected"], default="fixed")
    parser.add_argument("--best-tile", choices=["prefer-sst", "closest"], default="prefer-sst")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "results" / "ko_temperature")
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    main()
