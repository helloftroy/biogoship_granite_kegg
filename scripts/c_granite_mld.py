#!/usr/bin/env python3
"""R72: does Granite + 30-day mixed-layer depth beat OLCI + SST + MLD for C?

Same setup as main_model.py (C target, technical removal with read-pair depth, 13 cruises held out, ridge with the
leave-cruise-out penalty `cv`). Models: granite, granite_mld (32 Granite PCs + spline of 30-day mean MLD along the
back-trajectory, as sat_mld), sat, sat_mld, and granite_sat_mld (Granite PCs + OLCI/SST + MLD). Scored as % of
non-technical variance (pooled) and per cruise. Output: data/results/main_model/c_granite_mld.csv, *_per_cruise.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import main_model as mm  # noqa: E402
from tune_figure_kos import build  # noqa: E402
from two_stage_models import spl  # noqa: E402
from water_mass_test import features  # noqa: E402

MODELS = {"granite": (["P32"], False), "granite_mld": (["P32"], True), "sat": (["S~"], False), "sat_mld": (["S~"], True),
          "granite_sat": (["P32", "S~"], False), "granite_sat_mld": (["P32", "S~"], True)}  # fmt: skip


def main():
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    d = mm.setup("C", features(), ks)
    y, cr = d["ya"], d["cr"]
    pred = {m: np.full_like(y, np.nan) for m in MODELS}
    for c in np.unique(cr):
        te, tr = np.where(cr == c)[0], np.where(cr != c)[0]
        cache: dict = {}
        for m, (blocks, mld) in MODELS.items():
            xtr, xte = build(blocks, tr, te, d["emb"], d["temp"], d["sat"], cache)
            if mld:
                a, b = spl(d["mld"][tr], d["mld"][te])
                xtr, xte = np.hstack([xtr, a]), np.hstack([xte, b])
            pred[m][te], _ = mm.ridge(xtr, y[tr], xte, "cv", cr[tr])
    ss = ((y - y.mean()) ** 2).sum()
    rows, pc = [], []
    for m, p in pred.items():
        rows.append({"model": m, "explained_pct_nontechnical": 100 * (1 - ((y - p) ** 2).sum() / ss)})
        for c in np.unique(cr):
            k = cr == c
            mu_tr = y[~k].mean()
            pc.append({"model": m, "cruise": c, "tiles": int(k.sum()),
                       "r2": 1 - ((y[k] - p[k]) ** 2).sum() / ((y[k] - mu_tr) ** 2).sum(),
                       "r": np.corrcoef(y[k, 0], p[k, 0])[0, 1], "bias": float((y[k] - p[k]).mean())})  # fmt: skip
    t, pc = pd.DataFrame(rows), pd.DataFrame(pc)
    t.round(3).to_csv(mm.OUT / "c_granite_mld.csv", index=False)
    pc.round(3).to_csv(mm.OUT / "c_granite_mld_per_cruise.csv", index=False)
    print(t.round(2).to_string(index=False))
    w = pc.pivot(index="cruise", columns="model", values="r2")[list(MODELS)]
    print("\nper-cruise R2 (vs training-cruise mean):")
    print(w.round(2).to_string())
    print("\nper-cruise r:")
    print(pc.pivot(index="cruise", columns="model", values="r")[list(MODELS)].round(2).to_string())


if __name__ == "__main__":
    main()
