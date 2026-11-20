#!/usr/bin/env python3
"""Bloom phase: satellite chlorophyll history against the residual axes and as a KO predictor (RESULTS.md R37).

Inputs:
  transfer/results/chl_history/sample_chl_history.csv  (transfer/cluster/chl_history.sbatch on Morrill)
  the real-model targets (data/targets/ko_clr.csv, A16N Protocol A excluded) and the residual PCA
  refitted on them (transfer/results/residual_modules_noA/scores.csv)
Per tile: the mean over its samples of each feature (log10 for chlorophyll levels).

1. Spearman r of the axis scores (F1, F2, W1, W2, W3) with each feature, over all tiles and
   within cruise (R33 method), and how many cruises share the sign.
2. W2 from the features: within-cruise ridge, cruise held out (R33 loco_ridge), against the
   PISCES day chlorophyll alone and the satellite day chlorophyll alone.
3. As a KO predictor: the R33 baseline (32 PCs + temperature spline, cruise held out) plus one
   feature's spline (added_fields_all_kos.cv, as R27): median KO R2 and its change.
4. Spatial scale of the features along track (R34 correlogram).

Outputs: transfer/results/chl_history_test/ (drivers.csv, w2_ridge.csv, ko_r2.csv,
correlogram.csv, summary.txt).

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/chl_history_test.py
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
sys.path.insert(0, str(REPO / "transfer" / "scripts"))
from added_fields_all_kos import cv  # noqa: E402
from residual_axes_followup import BINS, correlogram  # noqa: E402
from residual_modules import loco_ridge, one_driver  # noqa: E402
from tune_figure_kos import load  # noqa: E402

OUT = REPO / "transfer" / "results" / "chl_history_test"
FEATURES = {
    "chlorophyll, day 0 (log)": ("chl_day0", True), "chlorophyll, 28-day mean (log)": ("chl_mean28", True),
    "chlorophyll, 28-day max (log)": ("chl_max28", True), "slope, 14 days (log10 per 10 d)": ("chl_slope14", False),
    "slope, 28 days (log10 per 10 d)": ("chl_slope28", False), "past the 28-day peak (log10 max/day0)": ("chl_peak_ratio28", False),
    "days since the 28-day max": ("chl_days_since_max", False), "CV over 28 days": ("chl_cv28", False),
}  # fmt: skip


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--targets", type=Path, default=REPO / "data" / "targets" / "ko_clr.csv")
    p.add_argument("--embeddings-dirs", type=Path, nargs="+",
                   default=[REPO / "data" / "embeddings_sst_l4", REPO / "data" / "faire" / "embeddings_sst_l4"])
    p.add_argument("--chl", type=Path, default=REPO / "transfer" / "results" / "chl_history" / "sample_chl_history.csv")
    p.add_argument("--scores", type=Path, default=REPO / "transfer" / "results" / "residual_modules_noA" / "scores.csv")
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, _, y, cruises, unit, rows = load(args, kos, return_rows=True)
    scores = pd.read_csv(args.scores, index_col="tile_key").loc[unit.index]
    chl = pd.read_csv(args.chl).drop_duplicates("sample_id").set_index("sample_id")

    def tile(v):
        return rows["sample_id"].map(v).astype(float).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()

    X = {}
    for label, (col, logged) in FEATURES.items():
        X[label] = tile(np.log10(chl[col].clip(lower=1e-3)) if logged else chl[col])
    fields = pd.read_csv(REPO / "data" / "model_fields" / "sample_model_fields.csv").set_index("sample_id")
    X["PISCES chlorophyll, day 0 (log; R25)"] = tile(np.log10(fields["pisces_chl"].clip(lower=0.001)))
    n_ok = np.isfinite(X["slope, 14 days (log10 per 10 d)"]).sum()
    say(f"{len(unit):,} tiles, {len(kos):,} KOs; chlorophyll history on {n_ok:,} tiles")

    # 1. drivers
    comps = ["F1", "F2", "W1", "W2", "W3"]
    rows_out = []
    say("\n== 1. Spearman r, all tiles / within cruise (cruises with the pooled sign)")
    for label, x in X.items():
        parts = []
        for c in comps:
            d = one_driver(label, scores[c].to_numpy(), x, cruises)
            if d:
                rows_out.append({**d, "component": c})
                parts.append(f"{c} {d['all_r']:+.2f}/{d['within_r']:+.2f} ({d['cruises_same_sign']}/{d['cruises_tested']})")
        say(f"   {label:40s} " + "  ".join(parts))
    pd.DataFrame(rows_out).round(4).to_csv(OUT / "drivers.csv", index=False)

    # 2. W2 ridge
    say("\n== 2. W2 from the features: within-cruise ridge, cruise held out")
    s = scores["W2"].to_numpy()
    hist = [k for k in FEATURES]
    sets = {"satellite chlorophyll, day 0": ["chlorophyll, day 0 (log)"], "PISCES chlorophyll, day 0": ["PISCES chlorophyll, day 0 (log; R25)"],
            "slopes + peak ratio + days since max": ["slope, 14 days (log10 per 10 d)", "slope, 28 days (log10 per 10 d)",
                                                     "past the 28-day peak (log10 max/day0)", "days since the 28-day max"],
            "all satellite history features": hist}  # fmt: skip
    r2rows = []
    for name, cols in sets.items():
        x = np.column_stack([X[c] for c in cols])
        keep = np.isfinite(x).all(1)
        xv = x[keep] - pd.DataFrame(x[keep]).groupby(cruises[keep]).transform("mean").to_numpy()
        r2 = loco_ridge(xv, s[keep], cruises[keep])
        r2rows.append({"predictors": name, "tiles": int(keep.sum()), "cv_r2": r2})
        say(f"   {name:40s} {int(keep.sum()):5d} tiles  CV R2 {r2:+.3f}")
    pd.DataFrame(r2rows).round(4).to_csv(OUT / "w2_ridge.csv", index=False)

    # 3. KO predictor
    say("\n== 3. Baseline (32 PCs + temperature spline) plus one feature's spline, cruise held out: median KO R2")
    keep = np.isfinite(np.column_stack([X[k] for k in hist])).all(1)
    idx = np.flatnonzero(keep)
    r2_0, _ = cv(emb, temp[:, 0], None, y, cruises, idx)
    say(f"   baseline on the same {len(idx):,} tiles: {np.median(r2_0):.3f}")
    korows = [{"added": "none", "median_r2": np.median(r2_0)}]
    for label in ["chlorophyll, day 0 (log)", "chlorophyll, 28-day mean (log)", "slope, 14 days (log10 per 10 d)",
                  "slope, 28 days (log10 per 10 d)", "past the 28-day peak (log10 max/day0)"]:
        r2, _ = cv(emb, temp[:, 0], X[label], y, cruises, idx)
        g = r2 - r2_0
        korows.append({"added": label, "median_r2": np.median(r2), "median_gain": np.median(g), "share_gain_0.02": (g > 0.02).mean(),
                       "share_loss_0.02": (g < -0.02).mean()})
        say(f"   + {label:40s} {np.median(r2):.3f} (median gain {np.median(g):+.4f}; gain > 0.02: {(g > 0.02).mean():.0%}, "
            f"loss > 0.02: {(g < -0.02).mean():.0%})")
    pd.DataFrame(korows).round(4).to_csv(OUT / "ko_r2.csv", index=False)

    # 4. spatial scale
    g = rows.groupby("tile_key")
    lat = g["sample_lat"].mean().loc[unit.index].to_numpy()
    lon = g["sample_lon"].mean().loc[unit.index].to_numpy()
    say("\n== 4. Within-cruise correlation of pairs by distance (km)")
    cg = []
    for label in ["slope, 14 days (log10 per 10 d)", "past the 28-day peak (log10 max/day0)", "chlorophyll, day 0 (log)"]:
        r = correlogram(X[label], cruises, lat, lon)
        cg.append({"series": label, **{f"{BINS[k]}-{BINS[k + 1]}": v[0] for k, v in r.items()}})
    r = correlogram(scores["W2"].to_numpy(), cruises, lat, lon)
    cg.append({"series": "W2 (real model)", **{f"{BINS[k]}-{BINS[k + 1]}": v[0] for k, v in r.items()}})
    cg = pd.DataFrame(cg).set_index("series")
    say(cg.round(2).to_string())
    cg.round(4).to_csv(OUT / "correlogram.csv")
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
