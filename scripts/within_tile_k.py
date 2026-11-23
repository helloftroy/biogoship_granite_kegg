#!/usr/bin/env python3
"""Does averaging samples in a tile bring the target closer to what the satellite can predict? (RESULTS.md R59)

Inside each tile with k >= 2 samples: held-out error of the main model against one sample at a time vs against the
tile mean (both technically adjusted with the tile's technical fit). Same tile, same prediction, so no cruise
confound. The gap equals sigma2 * (1 - 1/k) by construction (bias-variance identity), so this measures how
much error averaging removes, not whether it happens; the informative parts are sigma2 per k and the remaining error.
Output: data/results/main_model/within_tile_gain.csv
Usage: OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/within_tile_k.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import main_model as mm  # noqa: E402
from water_mass_test import features  # noqa: E402

fs = features()
ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
out = []
for tname, (m, p) in [("clr_all", ("granite", "cv")), ("C", ("sat_mld", "cv")), ("A", ("granite", "cv")),
                      ("B", ("granite", "cv")), ("D", ("granite", "cv"))]:  # fmt: skip
    d = mm.setup(tname, fs, ks)
    pred = np.load(mm.OUT / f"predictions_{tname}_{m}_{p}.npz")["predicted"].astype(float)
    st = pd.read_csv(mm.TARGETS[tname], index_col="sample_key")[d["kos"]]
    members = d["rows"].groupby("tile_key")["sample_id"].agg(list).reindex(d["unit"].index)
    for i, samples in enumerate(members):
        v = st.reindex([x for x in samples if x in st.index]).to_numpy()
        if len(v) < 2:
            continue
        v = v - d["tech"][i]
        out.append({"target": tname, "k": len(v), "cruise": d["cr"][i],
                    "err_single": float(((v - pred[i]) ** 2).sum(1).mean()), "err_mean": float(((v.mean(0) - pred[i]) ** 2).sum()),
                    "scatter": float(((v - v.mean(0)) ** 2).sum(1).sum() / (len(v) - 1))})  # fmt: skip
t = pd.DataFrame(out)
t["kg"] = np.where(t["k"] >= 4, "4-6", t["k"].astype(str))
g = t.groupby(["target", "kg"], sort=False).agg(tiles=("k", "size"), err_single=("err_single", "mean"), err_mean=("err_mean", "mean"),
                                                 sigma2=("scatter", "mean"), k=("k", "mean"))  # fmt: skip
g["gain"] = g["err_single"] - g["err_mean"]
g["gain_share"] = g["gain"] / g["err_single"]
g["expected_gain"] = g["sigma2"] * (1 - 1 / g["k"])
g.round(4).to_csv(mm.OUT / "within_tile_gain.csv")
print(g.round(3).to_string())
