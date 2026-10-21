#!/usr/bin/env python3
"""Quick check that the satellite features carry signal: CatBoost, one target at a time.

Compares feature sets on the same samples, so differences come from the features:

  fixed      Granite embedding of the best tile at the sampling position
  advected   Granite embedding of the best tile at the back-tracked position
  raw_fixed  mean and s.d. of the 17 raw tile channels, best fixed tile
             (needs data/embeddings/tile_raw_features.csv from pool_tile_features.py)
  raw_advected  the same for the best advected tile

"Best" is the clear tile closest in time to sampling. Only samples that have every
feature set compared are used, one row per sample.

Cross-validation is grouped by spatial block (section x 5-degree latitude band), so
neighbouring underway samples never sit on both sides of a split; random splits
would let near-duplicate neighbours inflate the scores. A mean-only baseline gives
R2 = 0 by construction, and a "metadata" set (latitude, longitude, day of year)
shows how much a model gets from position and season alone.

Targets: in-situ measurements in biogoship_samples.csv (temperature, phosphate,
nitrate + nitrite) until NMF and KEGG tables exist; any numeric column of a CSV
keyed by sample_id can be added with --targets-csv.

Outputs in data/results/:
  feature_set_comparison.csv          R2 and RMSE per target, feature set and fold
  feature_set_comparison_summary.csv  mean and s.d. over folds

Usage:
    python scripts/compare_feature_sets.py
    python scripts/compare_feature_sets.py --cv section      # hold out one cruise at a time
    python scripts/compare_feature_sets.py --targets-csv nmf_scores.csv --targets NMF1 NMF2
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import MANIFEST_DIR, REPO_ROOT  # noqa: E402

DEFAULT_TARGETS = ["temperature_c", "phosphate_umol_kg", "nitrate_nitrite_umol_l"]


def main() -> None:
    args = parse_args()
    from catboost import CatBoostRegressor

    rows = pd.read_csv(args.embeddings_dir / "sample_embeddings.csv")
    embeddings = np.load(args.embeddings_dir / "tile_embeddings.npy")
    best = rows[rows["is_best"] == 1]
    samples = best.drop_duplicates("sample_id").set_index("sample_id")

    features: dict[str, pd.DataFrame] = {}
    for position in ("fixed", "advected"):
        chosen = best[best["position"] == position].drop_duplicates("sample_id").set_index("sample_id")
        features[position] = pd.DataFrame(
            embeddings[chosen["embedding_row"].to_numpy()], index=chosen.index, columns=[f"e{i}" for i in range(embeddings.shape[1])]
        )
        raw_path = args.embeddings_dir / "tile_raw_features.csv"
        if raw_path.exists():
            raw = pd.read_csv(raw_path).set_index("tile_path")
            raw = raw[[c for c in raw.columns if c.endswith(("_mean", "_std"))]]
            features[f"raw_{position}"] = raw.reindex(chosen["tile_path"]).set_axis(chosen.index)
    sampled = pd.to_datetime(samples["sample_time"], utc=True)
    features["metadata"] = pd.DataFrame(
        {
            "lat": samples["sample_lat"],
            "lon": samples["sample_lon"],
            "doy_sin": np.sin(2 * np.pi * sampled.dt.dayofyear / 365.25),
            "doy_cos": np.cos(2 * np.pi * sampled.dt.dayofyear / 365.25),
        }
    )
    if args.feature_sets:
        features = {name: frame for name, frame in features.items() if name in args.feature_sets}

    # SST availability of the tile behind each sample (best fixed tile; advected if there is no fixed one).
    # A quarter of tiles have no SST at all, concentrated in I07N, so folds are balanced on it.
    best_tile = best.sort_values("position", key=lambda p: p.map({"fixed": 0, "advected": 1})).drop_duplicates("sample_id")
    has_sst = best_tile.set_index("sample_id")["sst_valid_fraction"].reindex(samples.index) > 0

    targets = load_targets(args, samples)
    # A sample counts if it has a row in every feature set. Missing values inside a row
    # (e.g. a tile with no SST) are left to CatBoost, which handles them natively.
    common = sorted(set.intersection(*(set(frame.index[frame.notna().any(axis=1)]) for frame in features.values())))
    print(f"Samples with every feature set ({', '.join(features)}): {len(common):,}")

    if args.cv == "section":  # each cruise held out in turn: can the model transfer to another ocean basin?
        blocks = samples.loc[common, "section"].astype(str)
    else:
        blocks = samples.loc[common, "section"].astype(str) + "_" + (samples.loc[common, "sample_lat"] // 5).astype(int).astype(str)
    results, predictions = [], []
    for target in targets.columns:
        y = targets.loc[common, target].dropna()
        if len(y) < 50:
            print(f"{target}: only {len(y)} samples with a value; skipped")
            continue
        n_folds = blocks.loc[y.index].nunique() if args.cv == "section" else args.folds
        folds = assign_folds(blocks.loc[y.index], n_folds, args.seed, flag=~has_sst.loc[y.index])
        no_sst_share = [float((~has_sst.loc[y.index]).to_numpy()[folds == k].mean()) for k in range(n_folds)]
        print(f"\n{target}: {len(y):,} samples, {blocks.loc[y.index].nunique()} groups, {n_folds} folds ({args.cv} CV); "
              f"share without SST per test fold: {', '.join(f'{v:.0%}' for v in no_sst_share)}")  # fmt: skip
        for name, frame in features.items():
            X = frame.loc[y.index].to_numpy(dtype=float)
            for fold in range(n_folds):
                test = folds == fold
                model = CatBoostRegressor(
                    iterations=args.iterations, depth=6, learning_rate=0.05, loss_function="RMSE",
                    random_seed=args.seed, verbose=False, thread_count=args.threads, allow_writing_files=False,
                )  # fmt: skip
                model.fit(X[~test], y.to_numpy()[~test])
                pred = model.predict(X[test])
                truth = y.to_numpy()[test]
                predictions.append(pd.DataFrame({
                    "sample_id": y.index[test], "target": target, "feature_set": name, "cv": args.cv, "fold": fold,
                    "observed": truth, "predicted": pred, "has_sst": has_sst.loc[y.index].to_numpy()[test],
                }))  # fmt: skip
                ss_res = float(((truth - pred) ** 2).sum())
                ss_tot = float(((truth - y.to_numpy()[~test].mean()) ** 2).sum())  # baseline: training-fold mean
                results.append(
                    {"target": target, "feature_set": name, "cv": args.cv, "fold": fold,
                     "held_out": blocks.loc[y.index][test].iloc[0] if args.cv == "section" else "", "n_test": int(test.sum()),
                     "r2": 1 - ss_res / ss_tot, "rmse": float(np.sqrt(ss_res / test.sum()))}
                )  # fmt: skip
            scores = [r["r2"] for r in results if r["target"] == target and r["feature_set"] == name]
            pooled = pd.concat([p for p in predictions if p["target"].iat[0] == target and p["feature_set"].iat[0] == name])
            by_sst = {flag: pooled_r2(pooled[pooled["has_sst"] == flag]) for flag in (True, False)}
            print(f"  {name:13s} R2 {np.mean(scores):6.3f} +/- {np.std(scores):.3f} | pooled R2, samples with SST "
                  f"{by_sst[True]:6.3f}, without {by_sst[False]:6.3f}", flush=True)  # fmt: skip

    out = pd.DataFrame(results)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.cv == "block" else f"_{args.cv}"
    out.to_csv(args.out_dir / f"feature_set_comparison{suffix}.csv", index=False)
    pd.concat(predictions).to_csv(args.out_dir / f"feature_set_predictions{suffix}.csv", index=False)
    summary = out.groupby(["target", "feature_set"]).agg(
        r2_mean=("r2", "mean"), r2_sd=("r2", "std"), rmse_mean=("rmse", "mean"), n_samples=("n_test", "sum")
    )
    summary.round(4).to_csv(args.out_dir / f"feature_set_comparison{suffix}_summary.csv")
    print(f"\nWrote {args.out_dir / f'feature_set_comparison{suffix}_summary.csv'}")


def load_targets(args: argparse.Namespace, samples: pd.DataFrame) -> pd.DataFrame:
    if args.targets_csv:
        table = pd.read_csv(args.targets_csv).set_index("sample_id")
        columns = args.targets or [c for c in table.columns if pd.api.types.is_numeric_dtype(table[c])]
        return table[columns].reindex(samples.index)
    catalogue = pd.read_csv(args.manifest_dir / "biogoship_samples.csv").set_index("sample_id")
    return catalogue[args.targets or DEFAULT_TARGETS].reindex(samples.index)


def assign_folds(blocks: pd.Series, n_folds: int, seed: int, flag: pd.Series | None = None) -> np.ndarray:
    """Whole spatial blocks to folds, largest first, keeping folds balanced in size and in flagged samples.

    Each block goes to the fold where, after adding it, the larger of (share of all
    samples, share of flagged samples) is smallest. With flag = "no SST", samples
    lacking SST are spread across folds as evenly as whole blocks allow.
    """
    flag = pd.Series(False, index=blocks.index) if flag is None else flag.astype(bool)
    table = pd.DataFrame({"n": blocks.value_counts(), "k": flag.groupby(blocks).sum()}).fillna(0)
    table = table.sample(frac=1.0, random_state=seed).sort_values(["n", "k"], ascending=False, kind="stable")
    total_n, total_k = max(table["n"].sum(), 1), max(table["k"].sum(), 1)
    load_n, load_k, fold_of = np.zeros(n_folds), np.zeros(n_folds), {}
    for block, row in table.iterrows():
        cost = np.maximum((load_n + row["n"]) / total_n, (load_k + row["k"]) / total_k)
        fold = int(np.argmin(cost + 1e-9 * np.arange(n_folds)))
        fold_of[block] = fold
        load_n[fold] += row["n"]
        load_k[fold] += row["k"]
    return blocks.map(fold_of).to_numpy()


def pooled_r2(frame: pd.DataFrame) -> float:
    if len(frame) < 5:
        return float("nan")
    residual = ((frame["observed"] - frame["predicted"]) ** 2).sum()
    return float(1 - residual / ((frame["observed"] - frame["observed"].mean()) ** 2).sum())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--embeddings-dir", type=Path, default=REPO_ROOT / "data" / "embeddings")
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "results")
    parser.add_argument("--targets-csv", type=Path, help="CSV with sample_id and target columns (e.g. NMF or KEGG).")
    parser.add_argument("--targets", nargs="+", help="Target columns. Default: in-situ temperature and nutrients.")
    parser.add_argument("--feature-sets", nargs="+", help="Subset of: fixed advected raw_fixed raw_advected metadata.")
    parser.add_argument(
        "--cv",
        choices=["block", "section"],
        default="block",
        help="block: 5 folds of section x 5-degree latitude blocks. section: hold out one cruise at a time.",
    )
    parser.add_argument("--folds", type=int, default=5, help="Folds for block CV.")
    parser.add_argument("--iterations", type=int, default=600)
    parser.add_argument("--threads", type=int, default=-1)
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


if __name__ == "__main__":
    main()
