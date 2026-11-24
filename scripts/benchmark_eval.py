#!/usr/bin/env python3
"""Standalone benchmark of satellite-tile embeddings against Bio-GO-SHIP metagenome targets (Granite protocol, R59-R64).

Needs only numpy, pandas, scikit-learn and the exported bundle (scripts/export_benchmark.py). Protocol, identical to
the Granite main model:
  * unit = tile (samples sharing a tile averaged; rows of tiles.csv / targets/*.csv.gz)
  * technical removal in-sample, per target column: OLS on [1, NovaSeq 6000 group, cubic spline (8 knots) of log10
    KO reads per tile], fitted on all tiles of that target; models are fitted and scored on the adjusted targets
  * features: standardise -> PCA (64 components, fitted on the training tiles) -> first --pcs components ->
    standardise (all inside each fold). --pcs 0 uses the standardised features directly
  * ridge, penalty per target column chosen by leave-one-cruise-out inside the training cruises
    (alphas logspace(-1, 5, 25))
  * one cruise held out at a time (13 cruises)
Scores: pooled explained share of the non-technical variance (headline), per-column pooled R2, variance partition
(technical / explained / cruise offsets / within-tile / within-cruise, % of total), per-cruise level (bias) and
pattern (within-cruise r) with within-cruise SD.

Usage:
  python benchmark_eval.py --bundle . --features granite/tile_embeddings.npy --name granite      # reproduces Granite
  python benchmark_eval.py --bundle . --features my_tiles.csv --name newfm --reference
     my_tiles.csv: column tile_key + one column per feature (or .npy with rows in tiles.csv order).
     Tiles with missing features are dropped; --reference reruns Granite on exactly the same tiles for a paired
     comparison (always use it when tiles are dropped).
  --targets clr_all C A B D Dp De (default all); --pcs 32 (default)
Outputs: <bundle>/eval_<name>/{scores.csv, ko_r2_<target>.csv, partition.csv, per_cruise.csv, predictions_<target>.npz}
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import SplineTransformer, StandardScaler

ALPHAS = np.logspace(-1, 5, 25)
TARGETS = ["clr_all", "C", "D", "A", "Dp", "B", "De"]


def ridge_path(xtr, ytr, xte):
    xm, ym = xtr.mean(0), ytr.mean(0)
    u, s, vt = np.linalg.svd(xtr - xm, full_matrices=False)
    uty = u.T @ (ytr - ym)
    proj = (xte - xm) @ vt.T
    shrink = s[None, :] / (s[None, :] ** 2 + ALPHAS[:, None])
    return ym + np.einsum("tc,ac,ck->atk", proj, shrink, uty)


def ridge_loco_alpha(xtr, ytr, xte, groups):
    sse = np.zeros((len(ALPHAS), ytr.shape[1]))
    for g in np.unique(groups):
        a, b = groups != g, groups == g
        sse += ((ytr[b][None] - ridge_path(xtr[a], ytr[a], xtr[b])) ** 2).sum(axis=1)
    idx = sse.argmin(0)
    return ridge_path(xtr, ytr, xte)[idx, :, np.arange(ytr.shape[1])].T


def features_fold(x, tr, te, pcs):
    s = StandardScaler().fit(x[tr])
    a, b = s.transform(x[tr]), s.transform(x[te])
    if pcs:
        p = PCA(n_components=min(64, a.shape[1], len(tr)), random_state=0).fit(a)
        a, b = p.transform(a)[:, :pcs], p.transform(b)[:, :pcs]
        s2 = StandardScaler().fit(a)
        a, b = s2.transform(a), s2.transform(b)
    return a, b


def technical(y, plat, logtot):
    X = np.hstack([np.ones((len(y), 1)), plat[:, None].astype(float), SplineTransformer(n_knots=8, degree=3).fit_transform(logtot[:, None])])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ beta
    return y - (fit - y.mean(0)), fit - fit.mean(0)


def run(bundle, feats, name, targets, pcs, tiles_keep):
    tiles = pd.read_csv(bundle / "tiles.csv", index_col="tile_key")
    samples = pd.read_csv(bundle / "samples.csv", index_col="sample_key")
    pairs = pd.read_csv(bundle / "within_tile_pairs.csv")
    out_dir = bundle / f"eval_{name}"
    out_dir.mkdir(exist_ok=True)
    scores, parts, cruises = [], [], []
    for t in targets:
        Y = pd.read_csv(bundle / "targets" / f"{t}.csv.gz", index_col="tile_key")
        Y = Y.loc[[k for k in Y.index if k in tiles_keep]]
        tt = tiles.loc[Y.index]
        x = feats.loc[Y.index].to_numpy(float)
        y = Y.to_numpy(float)
        cr, plat = tt["cruise"].to_numpy(), (tt["platform_group"] == "NovaSeq 6000 group").to_numpy()
        ya, tech = technical(y, plat, tt["log10_ko_reads"].to_numpy())
        pred = np.full_like(ya, np.nan)
        for c in np.unique(cr):
            te, tr = np.where(cr == c)[0], np.where(cr != c)[0]
            a, b = features_fold(x, tr, te, pcs)
            pred[te] = ridge_loco_alpha(a, ya[tr], b, cr[tr])
        np.savez_compressed(out_dir / f"predictions_{t}.npz", predicted=pred.astype(np.float32), observed=ya.astype(np.float32),
                            tile_key=Y.index.to_numpy(), columns=np.array(Y.columns), cruise=cr)  # fmt: skip
        yac = ya - ya.mean(0)
        res = ya - pred
        r2 = 1 - (res**2).sum(0) / (yac**2).sum(0)
        pd.DataFrame({"column": Y.columns, "r2": r2}).to_csv(out_dir / f"ko_r2_{t}.csv", index=False)
        tot = ((y - y.mean(0)) ** 2).sum()
        expl = (yac**2).sum() - (res**2).sum()
        cm = pd.DataFrame(res).groupby(cr).transform("mean").to_numpy()
        within = ((res - cm) ** 2).sum()
        # within-tile variability from same-tile sample pairs (raw sample targets; the technical shift is per tile)
        st = pd.read_csv(bundle / "targets_sample" / f"{t}.csv.gz", index_col="sample_key")[Y.columns]
        pp = pairs[pairs["tile_key"].isin(Y.index) & pairs["sample_a"].isin(st.index) & pairs["sample_b"].isin(st.index)]
        diff = st.loc[pp["sample_a"]].to_numpy() - st.loc[pp["sample_b"]].to_numpy()
        sigma2 = (diff**2).mean(0) / 2
        noise = min((sigma2.sum()) * (1 / tt["k_samples"].to_numpy()).sum(), within)
        p = {"target": t, "features": name, "tiles": len(y), "columns": y.shape[1], "technical": (tech**2).sum() / tot,
             "explained": expl / tot, "cruise_offsets": (cm**2).sum() / tot, "within_tile": noise / tot,
             "within_cruise": (within - noise) / tot}  # fmt: skip
        p["explained_of_nontechnical"] = p["explained"] / (1 - p["technical"])
        p["explained_of_reproducible"] = p["explained"] / (1 - p["technical"] - p["within_tile"])
        p["median_column_r2"] = float(np.median(r2))
        parts.append(p)
        for c in np.unique(cr):
            m = cr == c
            o, q = ya[m], pred[m]
            so, sp = o.std(0), q.std(0)
            r = ((o - o.mean(0)) * (q - q.mean(0))).mean(0) / np.maximum(so * sp, 1e-12)
            cruises.append({"target": t, "features": name, "cruise": c, "tiles": int(m.sum()), "bias_mean": float((o - q).mean()),
                            "abs_bias_median": float(np.median(np.abs((o - q).mean(0)))), "sd_obs_median": float(np.median(so)),
                            "sd_pred_median": float(np.median(sp)), "r_median": float(np.median(r))})  # fmt: skip
        print(f"  {name:12s} {t:8s} {len(y):5,} tiles  explained {p['explained_of_nontechnical']:7.2%} of non-technical, "
              f"{p['explained_of_reproducible']:7.2%} of reproducible; median column R2 {p['median_column_r2']:.4f}", flush=True)
    pd.DataFrame(parts).round(5).to_csv(out_dir / "partition.csv", index=False)
    pd.DataFrame(cruises).round(5).to_csv(out_dir / "per_cruise.csv", index=False)
    return pd.DataFrame(parts)


def load_features(path, tiles):
    path = Path(path)
    if path.suffix == ".npy":
        x = np.load(path)
        assert len(x) == len(tiles), f"{path}: {len(x)} rows, tiles.csv has {len(tiles)}"
        return pd.DataFrame(x, index=tiles.index)
    f = pd.read_csv(path)
    return f.set_index("tile_key").reindex(tiles.index)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", type=Path, default=Path(__file__).resolve().parent)
    ap.add_argument("--features", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--targets", nargs="+", default=TARGETS)
    ap.add_argument("--pcs", type=int, default=32)
    ap.add_argument("--reference", action="store_true", help="also run Granite on the same tiles")
    a = ap.parse_args()
    tiles = pd.read_csv(a.bundle / "tiles.csv", index_col="tile_key")
    feats = load_features(a.features if Path(a.features).is_absolute() else a.bundle / a.features, tiles)
    keep = set(feats.index[feats.notna().all(1)])
    print(f"{len(keep):,} of {len(tiles):,} tiles have features ({a.name})")
    res = [run(a.bundle, feats, a.name, a.targets, a.pcs, keep)]
    if a.reference:
        g = load_features(a.bundle / "granite" / "tile_embeddings.npy", tiles)
        res.append(run(a.bundle, g, f"granite_on_{a.name}_tiles", a.targets, 32, keep))
    t = pd.concat(res)
    t.round(5).to_csv(a.bundle / f"eval_{a.name}" / "scores.csv", index=False)
    print(t.pivot_table(index="target", columns="features", values="explained_of_nontechnical", sort=False).mul(100).round(2).to_string())


if __name__ == "__main__":
    main()
