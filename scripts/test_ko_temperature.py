#!/usr/bin/env python3
"""Does in-situ temperature (with position, season or ocean colour) predict KOs as well as Granite?

KOs (--kos):
  selected  20 KOs from the screen (data/results/ko_screen/ko_screen.csv): the 10
            best predicted by the embedding, and 10 middling ones (embedding pooled
            R2 0.25-0.40 and at least 0.15 above position + season), evenly spaced
            in rank.
  all       every KO in the targets table.

Feature sets (--feature-sets; any of):
  metadata          position and season (as in the screen: lat/90, sine and cosine
                    of longitude and of day of year)
  metadata_temp     the same + in-situ temperature (manifest temperature_c)
  temp              in-situ temperature alone
  colour            mean of the 16 OLCI bands over the tile's clear pixels
                    (Granite-normalised; from tile_raw_features.csv)
  temp_colour       temperature + the 16 band means
  temp_chl          temperature + a band-ratio chlorophyll index,
                    log10(max(Rrs443, Rrs490, Rrs510) / Rrs560) from the band means
                    (the ratio behind OC4ME; higher = less chlorophyll)
  embedding         Granite embedding of the sample's tile (512 numbers)
  embedding_temp    embedding + in-situ temperature
All inputs are standardised on the training cruises.

Learners (--learners): ridge, spline (a cubic B-spline of each input, 8 knots,
then ridge: a smooth curve, cheap for thousands of KOs), k nearest neighbours (k = 15, distance-weighted),
CatBoost (500 trees, depth 6, learning rate 0.05, MultiRMSE; skipped for feature
sets wider than --catboost-max-features, i.e. the embedding; so is spline). Per KO and feature
set the summary takes the best learner.

Samples: as in the screen (fixed position by default, prefer-SST tile, SST-only,
identical tiles merged), restricted to samples that have every requested feature
(in-situ temperature; band means when a colour set is requested), so all feature
sets see the same rows. Temperature is averaged over merged samples. One cruise is
held out at a time.

Outputs in --out-dir (default data/results/ko_temperature/): ko_temperature.csv
(KO x feature set x learner: pooled R2 and median R2 within a held-out cruise) and
ko_temperature_per_cruise.csv.

Usage:
    python scripts/test_ko_temperature.py
    python scripts/test_ko_temperature.py --kos all --learners ridge spline knn --feature-sets temp embedding \\
        --out-dir data/results/ko_temperature_all
    python scripts/test_ko_temperature.py --learners ridge knn catboost \\
        --feature-sets temp colour temp_colour temp_chl embedding --out-dir data/results/ko_colour
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship import tiles  # noqa: E402
from biogoship.tables import REPO_ROOT  # noqa: E402
from screen_ko_predictability import ALPHAS, load_embeddings, metadata_features  # noqa: E402

FEATURE_SETS = ("metadata", "metadata_temp", "temp", "colour", "temp_colour", "temp_chl", "embedding", "embedding_temp")
DEFAULT_SETS = ("metadata", "metadata_temp", "temp", "embedding", "embedding_temp")
LEARNERS = ("ridge", "spline", "knn", "catboost")
BAND_MEANS = [f"{band}_mean" for band in tiles.OLCI_BANDS]


def main() -> None:
    args = parse_args()
    kos = pick_kos(args)
    targets = pd.read_csv(args.targets, index_col="sample_key", usecols=["sample_key", *kos["ko"]])[list(kos["ko"])]
    rows, X_embed = load_embeddings(args)
    rows = rows[rows["sample_id"].isin(targets.index) & rows["temperature_c"].notna()]
    needs_colour = any("colour" in name or "chl" in name for name in args.feature_sets)
    if needs_colour:
        raw = load_raw_features(args.embeddings_dirs)
        rows = rows[rows["tile_path"].isin(raw.index)]
    content = pd.Series([v.tobytes() for v in X_embed[rows["embedding_row"].to_numpy()].round(5)], index=rows.index)
    rows = rows.assign(tile_key=content.map({k: i for i, k in enumerate(content.unique())}))
    Y = targets.loc[rows["sample_id"]].set_axis(rows.index).groupby(rows["tile_key"]).mean()
    temp = rows.groupby("tile_key")["temperature_c"].mean().loc[Y.index].to_numpy()[:, None]
    unit = rows.drop_duplicates("tile_key").set_index("tile_key").loc[Y.index]
    meta, emb = metadata_features(unit), X_embed[unit["embedding_row"].to_numpy()]
    X = {"metadata": meta, "metadata_temp": np.hstack([meta, temp]), "temp": temp,
         "embedding": emb, "embedding_temp": np.hstack([emb, temp])}  # fmt: skip
    if needs_colour:
        colour = raw.loc[unit["tile_path"], BAND_MEANS].to_numpy()
        chl = chlorophyll_index(colour)[:, None]
        X.update({"colour": colour, "temp_colour": np.hstack([temp, colour]), "temp_chl": np.hstack([temp, chl])})
    cruises = unit["section"].astype(str).to_numpy()
    print(f"Samples with a KO profile, an in-situ temperature{', band means' if needs_colour else ''} and a "
          f"{args.position} tile with SST: {len(rows):,}; distinct tiles: {len(unit):,}; cruises: {len(set(cruises))}; "
          f"KOs: {Y.shape[1]:,}")  # fmt: skip

    y = Y.to_numpy()
    results, per_cruise = [], []
    for name in args.feature_sets:
        for learner in args.learners:
            if learner in ("catboost", "spline") and X[name].shape[1] > args.catboost_max_features:
                continue  # slow or pointless on the 512-number embedding; ridge wins there (R5)
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
            if args.save_predictions:
                args.out_dir.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(args.out_dir / f"predictions_{name}_{learner}.npz", predicted=pred.astype(np.float32),
                                    observed=y.astype(np.float32), ko=np.array(Y.columns))  # fmt: skip
            results.append(pd.DataFrame({"ko": Y.columns, "feature_set": name, "learner": learner, "r2_pooled": pooled}))
            print(f"  {name:15s} {learner:9s} median pooled R2 over KOs {np.median(pooled):6.3f}", flush=True)

    per_cruise = pd.concat(per_cruise)
    out = pd.concat(results).merge(
        per_cruise.groupby(["ko", "feature_set", "learner"])["r2"].median().rename("r2_median_cruise").reset_index())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out.round(4).to_csv(args.out_dir / "ko_temperature.csv", index=False)
    if args.save_predictions:  # the rows of the prediction files: one per distinct tile
        cols = ["sample_id", "section", "sample_time", "sample_lat", "sample_lon", "tile_path"]
        table = unit[[c for c in cols if c in unit]].assign(temperature_c=temp[:, 0])
        table.assign(n_samples=rows.groupby("tile_key").size().loc[Y.index].to_numpy()).to_csv(
            args.out_dir / "prediction_rows.csv", index=False)
    per_cruise.round(4).to_csv(args.out_dir / "ko_temperature_per_cruise.csv", index=False)

    best = out.loc[out.groupby(["ko", "feature_set"])["r2_pooled"].idxmax()]
    pooled = best.pivot(index="ko", columns="feature_set", values="r2_pooled")[list(args.feature_sets)]
    within = best.pivot(index="ko", columns="feature_set", values="r2_median_cruise")[list(args.feature_sets)]
    groups = kos.set_index("ko")["group"]
    pd.set_option("display.width", 250)
    print("\nMedian over KOs, pooled R2 (best learner per KO and feature set):")
    print(pooled.groupby(groups).median().round(2).to_string())
    print("\nMedian over KOs, median R2 within a held-out cruise (same learner):")
    print(within.groupby(groups).median().round(2).to_string())
    print("\nWinning learner (KOs):")
    print(best.groupby(["feature_set", "learner"]).size().unstack(fill_value=0).loc[list(args.feature_sets)].to_string())
    if len(kos) <= 40:
        print("\nPooled R2 per KO:")
        print(kos.set_index("ko")[["group", "name"]].join(pooled.round(2)).to_string())
    print(f"\nWrote {args.out_dir / 'ko_temperature.csv'}")


def load_raw_features(folders: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(folder / "tile_raw_features.csv") for folder in folders if (folder / "tile_raw_features.csv").exists()]
    if not frames:
        raise SystemExit("No tile_raw_features.csv in the embeddings folders (scripts/pool_tile_features.py makes it).")
    raw = pd.concat(frames).drop_duplicates("tile_path").set_index("tile_path")
    return raw[raw[BAND_MEANS].notna().all(axis=1)]


def chlorophyll_index(band_means: np.ndarray) -> np.ndarray:
    """log10(max blue / green) from Granite-normalised band means (Oa03, Oa04, Oa05 against Oa06)."""
    reflectance = band_means * tiles.OLCI_STDS + tiles.OLCI_MEANS
    blue = reflectance[:, [2, 3, 4]].max(axis=1)
    green = reflectance[:, 5]
    return np.log10(np.clip(blue, 1e-5, None) / np.clip(green, 1e-5, None))


def fit_predict(learner: str, train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, seed: int) -> np.ndarray:
    if learner == "ridge":
        return RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(train_x, train_y).predict(test_x)
    if learner == "spline":  # a smooth curve per input (cubic B-splines, 8 knots), then ridge: bends, cheap for many KOs
        model = make_pipeline(SplineTransformer(n_knots=8, degree=3, extrapolation="linear"),
                              RidgeCV(alphas=ALPHAS, alpha_per_target=True))  # fmt: skip
        return model.fit(train_x, train_y).predict(test_x)
    if learner == "knn":
        return KNeighborsRegressor(n_neighbors=15, weights="distance").fit(train_x, train_y).predict(test_x)
    from catboost import CatBoostRegressor

    model = CatBoostRegressor(iterations=500, depth=6, learning_rate=0.05, loss_function="MultiRMSE",
                              random_seed=seed, verbose=False, thread_count=-1)  # fmt: skip
    return model.fit(train_x, train_y).predict(test_x)


def pick_kos(args: argparse.Namespace) -> pd.DataFrame:
    screen = pd.read_csv(args.screen).sort_values("r2_pooled_embedding", ascending=False)
    if args.kos == "all":
        header = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key")
        kos = pd.DataFrame({"ko": header, "group": "all"})
    else:
        top = screen.head(10).assign(group="top")
        margin = screen["r2_pooled_embedding"] - screen["r2_pooled_metadata"]
        middle = screen[screen["r2_pooled_embedding"].between(0.25, 0.40) & (margin >= 0.15)]
        middle = middle.iloc[np.linspace(0, len(middle) - 1, 10).round().astype(int)].assign(group="middle")
        kos = pd.concat([top, middle])[["ko", "group"]]
    names = pd.read_csv(args.ko_names, sep="\t", header=None, names=["ko", "name"]).set_index("ko")["name"]
    kos["name"] = kos["ko"].map(names).fillna("").str.replace(r"\s*\[EC:[^\]]*\]", "", regex=True)
    return kos.reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--screen", type=Path, default=REPO_ROOT / "data" / "results" / "ko_screen" / "ko_screen.csv")
    parser.add_argument("--ko-names", type=Path, default=REPO_ROOT / "data" / "reference" / "kegg_ko_list.tsv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings", REPO_ROOT / "data" / "faire" / "embeddings"])  # fmt: skip
    parser.add_argument("--position", choices=["fixed", "advected"], default="fixed")
    parser.add_argument("--best-tile", choices=["prefer-sst", "closest"], default="prefer-sst")
    parser.add_argument("--kos", choices=["selected", "all"], default="selected")
    parser.add_argument("--feature-sets", nargs="+", choices=FEATURE_SETS, default=list(DEFAULT_SETS))
    parser.add_argument("--learners", nargs="+", choices=LEARNERS, default=list(LEARNERS))
    parser.add_argument("--catboost-max-features", type=int, default=64,
                        help="Skip CatBoost for feature sets wider than this (the embedding). R5 ran it on all (1e9).")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "results" / "ko_temperature")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--save-predictions", action="store_true",
                        help="Write the held-out predictions (predictions_<features>_<learner>.npz, prediction_rows.csv).")
    return parser.parse_args()


if __name__ == "__main__":
    main()
