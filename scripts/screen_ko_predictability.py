#!/usr/bin/env python3
"""Which KOs can the Granite embedding predict? A fast screen of every KO at once.

Joins the CLR targets from build_ko_targets.py to the Granite embeddings of both
sample sets (original cruises in data/embeddings, FAIRe cruises in
data/faire/embeddings) and fits ridge regressions, one cruise held out at a time.
Ridge solves all targets in one go (alpha chosen per target by efficient
leave-one-out on the training cruises), so thousands of KOs take minutes; the
best ones are then worth a CatBoost run with compare_feature_sets.py.

Feature sets, each with two learners (ridge = linear; knn = average of the 15 most
similar training samples, which bends freely but cannot extrapolate):
  embedding  Granite embedding of the sample's tile (512 numbers)
  metadata   position and season: latitude, longitude and day of year (longitude
             and day of year as sine and cosine)

Tile per sample (as in compare_feature_sets.py): --best-tile prefer-sst takes the
closest-in-time tile that has SST, closest the closest tile; either way tiles
without SST are not used (the sample drops out). Samples sharing an identical tile
(replicates, or neighbours snapped to one clear window) become one row, targets
averaged.

Scores per KO: R2 on each held-out cruise against the training mean, their median,
and the pooled R2 over all held-out predictions. A null run with the embeddings
shuffled between samples gives the R2 to expect by chance.

Outputs in data/results/ko_screen/: ko_screen.csv (one row per KO) and
ko_screen_per_cruise.csv.

Usage:
    python scripts/screen_ko_predictability.py
    python scripts/screen_ko_predictability.py --position advected --best-tile closest
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

ALPHAS = np.logspace(-1, 5, 25)


def main() -> None:
    args = parse_args()
    targets = pd.read_csv(args.targets, index_col="sample_key")
    if args.instrument:  # one sequencing batch only: is a KO predictable without batch differences between cruises?
        samples = pd.read_csv(args.targets.with_name(args.targets.name.replace("_clr", "_samples")), index_col="sample_key")
        targets = targets[targets.index.isin(samples.index[samples["instrument"] == args.instrument])]
    rows, X_embed = load_embeddings(args)
    rows = rows[rows["sample_id"].isin(targets.index)]
    # Identical tiles -> one row, targets averaged.
    content = pd.Series([v.tobytes() for v in X_embed[rows["embedding_row"].to_numpy()].round(5)], index=rows.index)
    rows = rows.assign(tile_key=content.map({k: i for i, k in enumerate(content.unique())}))
    Y = targets.loc[rows["sample_id"]].set_axis(rows.index).groupby(rows["tile_key"]).mean()
    unit = rows.drop_duplicates("tile_key").set_index("tile_key").loc[Y.index]
    X = {"embedding": X_embed[unit["embedding_row"].to_numpy()], "metadata": metadata_features(unit)}
    cruises = unit["section"].astype(str).to_numpy()
    print(f"Samples with a target and a {args.position} tile with SST ({args.best_tile}): {len(rows):,}; "
          f"distinct tiles: {len(unit):,}; KOs: {Y.shape[1]:,}")  # fmt: skip
    print(pd.Series(cruises).value_counts().to_string())

    rng = np.random.default_rng(args.seed)
    X["embedding_shuffled"] = X["embedding"][rng.permutation(len(unit))]
    runs = {f"{name}_{learner}": (features, learner) for name, features in X.items() for learner in ("ridge", "knn")}
    y = Y.to_numpy()
    results, per_cruise = {}, []
    for name, (features, learner) in runs.items():
        pred = np.full_like(y, np.nan)
        for cruise in np.unique(cruises):
            test = cruises == cruise
            scaler = StandardScaler().fit(features[~test])
            if learner == "ridge":
                model = RidgeCV(alphas=ALPHAS, alpha_per_target=True)
            else:
                model = KNeighborsRegressor(n_neighbors=15, weights="distance")
            model.fit(scaler.transform(features[~test]), y[~test])
            pred[test] = model.predict(scaler.transform(features[test]))
            ss_res = ((y[test] - pred[test]) ** 2).sum(axis=0)
            ss_tot = ((y[test] - y[~test].mean(axis=0)) ** 2).sum(axis=0)
            per_cruise.append(pd.DataFrame({"ko": Y.columns, "feature_set": name, "held_out": cruise, "n": int(test.sum()),
                                            "r2": 1 - ss_res / ss_tot}))  # fmt: skip
        pooled = 1 - ((y - pred) ** 2).sum(axis=0) / ((y - y.mean(axis=0)) ** 2).sum(axis=0)
        results[name] = pooled
        print(f"  {name:19s} pooled R2 over KOs: median {np.median(pooled):6.3f}, 90th pct {np.quantile(pooled, 0.9):6.3f}, "
              f"max {pooled.max():6.3f}", flush=True)  # fmt: skip

    per_cruise = pd.concat(per_cruise)
    median = per_cruise.groupby(["ko", "feature_set"])["r2"].median().unstack()
    out = pd.DataFrame({"ko": Y.columns, "mean_clr": y.mean(axis=0), "sd_clr": y.std(axis=0)}).set_index("ko")
    for name in runs:
        out[f"r2_pooled_{name}"] = results[name]
        out[f"r2_median_cruise_{name}"] = median[name].reindex(out.index)
    for name in ("embedding", "metadata", "embedding_shuffled"):  # the better learner per KO
        out[f"r2_pooled_{name}"] = out[[f"r2_pooled_{name}_ridge", f"r2_pooled_{name}_knn"]].max(axis=1)
        out[f"r2_median_cruise_{name}"] = out[[f"r2_median_cruise_{name}_ridge", f"r2_median_cruise_{name}_knn"]].max(axis=1)
    null_cut = np.quantile(out["r2_pooled_embedding_shuffled"], 0.99)
    out = out.sort_values("r2_pooled_embedding", ascending=False)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out.round(4).to_csv(args.out_dir / "ko_screen.csv")
    per_cruise.round(4).to_csv(args.out_dir / "ko_screen_per_cruise.csv", index=False)

    print("Scores below take the better of the two learners per KO for each feature set.")
    print(f"\nShuffled-embedding null: 99th percentile of pooled R2 = {null_cut:.3f}")
    for cut in (0.1, 0.25, 0.5):
        e, m = out["r2_pooled_embedding"] >= cut, out["r2_pooled_metadata"] >= cut
        print(f"KOs with pooled R2 >= {cut}: embedding {int(e.sum()):,}, metadata {int(m.sum()):,}, "
              f"embedding but not metadata {int((e & ~m).sum()):,}")  # fmt: skip
    better = (out["r2_pooled_embedding"] - out["r2_pooled_metadata"])
    print(f"Embedding beats metadata for {int((better > 0).sum()):,} of {len(out):,} KOs (median difference {better.median():+.3f})")
    print(f"\nTop 15 by embedding:\n{out.head(15)[['r2_pooled_embedding', 'r2_median_cruise_embedding', 'r2_pooled_metadata', 'r2_median_cruise_metadata']].to_string()}")
    print(f"\nWrote {args.out_dir / 'ko_screen.csv'}")


def load_embeddings(args: argparse.Namespace) -> tuple[pd.DataFrame, np.ndarray]:
    """Chosen tile per sample from both sample sets, with one shared embedding matrix."""
    frames, matrices, offset = [], [], 0
    for folder in args.embeddings_dirs:
        rows = pd.read_csv(folder / "sample_embeddings.csv")
        emb = np.load(folder / "tile_embeddings.npy")
        rows = rows[rows["position"] == args.position].copy()
        has = pd.to_numeric(rows["sst_valid_fraction"], errors="coerce").fillna(0) > 0
        if args.best_tile == "prefer-sst":
            chosen = rows[has].assign(gap=rows["match_hour_gap"].abs()).sort_values("gap").drop_duplicates("sample_id")
        else:
            chosen = rows[(rows["is_best"] == 1) & has]
        chosen["embedding_row"] += offset
        frames.append(chosen)
        matrices.append(emb)
        offset += len(emb)
    return pd.concat(frames, ignore_index=True), np.vstack(matrices)


def metadata_features(unit: pd.DataFrame) -> np.ndarray:
    doy = pd.to_datetime(unit["sample_time"], utc=True).dt.dayofyear.to_numpy()
    lon = np.radians(unit["sample_lon"].to_numpy())
    return np.column_stack([unit["sample_lat"].to_numpy() / 90.0, np.sin(lon), np.cos(lon),
                            np.sin(2 * np.pi * doy / 365.25), np.cos(2 * np.pi * doy / 365.25)])  # fmt: skip


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings", REPO_ROOT / "data" / "faire" / "embeddings"])  # fmt: skip
    parser.add_argument("--position", choices=["fixed", "advected"], default="fixed")
    parser.add_argument("--best-tile", choices=["prefer-sst", "closest"], default="prefer-sst")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "results" / "ko_screen")
    parser.add_argument("--instrument", help="Only samples sequenced on this instrument (e.g. OBI_2631, OBI_3310).")
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    main()
