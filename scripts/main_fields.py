#!/usr/bin/env python3
"""Ocean-model and history fields added to the main model (Granite PCs), every KO, technical variance removed (R61).

Redo of R28 (scripts/added_fields_all_kos.py) on the R59/R60 main model: base = 32 Granite PCs, ridge, penalty by
leave-cruise-out inside the training cruises (cv), on the technically adjusted all-KO CLR targets (in-sample platform +
read-depth spline removed). Added = the same + a spline (8 knots, cubic) of one variable, transformed as in R28.
Variables missing on some tiles (the forecast-period PISCES fields, from 2021-11) run on their tiles with their own
base. Controls: three standard-normal noise columns and shuffled phosphate (all tiles and the 637 forecast tiles).
Per KO: pooled R2 base and added, gain, held-out cruises with lower squared error.

Output: data/results/main_model/added_fields.csv (+ added_fields_summary.csv)
Usage: OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/main_fields.py [--shard i n]   (shards: variables i, i+n, ...;
       then run with --merge n to combine added_fields_shard*.csv)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import main_model as mm  # noqa: E402
from model_field_first import FIELDS, HISTORY, HISTORY_FIELDS, VARIABLES  # noqa: E402
from tune_figure_kos import build  # noqa: E402
from two_stage_models import spl  # noqa: E402
from water_mass_test import features  # noqa: E402


def cv(d, v, idx):
    y, cr = d["ya"], d["cr"]
    pred = np.full_like(y, np.nan)
    for c in np.unique(cr[idx]):
        te, tr = idx[cr[idx] == c], idx[cr[idx] != c]
        xtr, xte = build(["P32"], tr, te, d["emb"], d["temp"], d["sat"], {})
        if v is not None:
            a, b = spl(v[tr, None], v[te, None])
            xtr, xte = np.hstack([xtr, a]), np.hstack([xte, b])
        pred[te], _ = mm.ridge(xtr, y[tr], xte, "cv", cr[tr])
    yi = y[idx]
    return 1 - ((yi - pred[idx]) ** 2).sum(0) / ((yi - yi.mean(0)) ** 2).sum(0), pred


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", nargs=2, type=int)
    ap.add_argument("--merge", type=int)
    a = ap.parse_args()
    if a.merge:
        t = pd.concat([pd.read_csv(mm.OUT / f"added_fields_shard{i}.csv") for i in range(a.merge)])
        return summarise(t)
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    d = mm.setup("clr_all", features(), ks)
    rows, unit = d["rows"], d["unit"]
    fields = pd.read_csv(FIELDS).set_index("sample_id").join(pd.read_csv(HISTORY_FIELDS).set_index("sample_id")
                                                             .drop(columns=["section", "sample_time", "sample_lat", "sample_lon"]))  # fmt: skip
    fields["speed"] = np.hypot(fields["glorys_uo"], fields["glorys_vo"])
    per_tile = lambda col: rows["sample_id"].map(fields[col]).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()  # noqa: E731
    rng = np.random.default_rng(0)
    variables = {}
    for name, (col, floor) in {**VARIABLES, **HISTORY}.items():
        v = per_tile(col)
        variables[name] = np.log10(np.clip(v, floor, None)) if floor is not None else v
    recent = ~np.isnan(variables["iron"])
    for k in (1, 2, 3):
        variables[f"CONTROL: noise {k}"] = rng.standard_normal(len(unit))
    variables["CONTROL: noise 1 (forecast tiles)"] = np.where(recent, rng.standard_normal(len(unit)), np.nan)
    variables["CONTROL: shuffled phosphate"] = rng.permutation(variables["phosphate"])
    sp = np.full(len(unit), np.nan)
    sp[recent] = rng.permutation(variables["phosphate"][recent])
    variables["CONTROL: shuffled phosphate (forecast tiles)"] = sp
    print(f"{len(unit):,} tiles, {len(d['kos']):,} KOs, {len(variables)} variables", flush=True)
    if a.shard:
        names = list(variables)[a.shard[0]::a.shard[1]]
        variables = {k: variables[k] for k in names}
    bases, out = {}, []
    for name, v in variables.items():
        idx = np.flatnonzero(~np.isnan(v))
        if len(idx) not in bases:
            bases[len(idx)] = cv(d, None, idx)
        base, bp = bases[len(idx)]
        added, ap = cv(d, v, idx)
        sse = lambda p: pd.DataFrame((d["ya"][idx] - p[idx]) ** 2).groupby(d["cr"][idx]).sum()  # noqa: E731
        improved = (sse(ap) < sse(bp)).sum().to_numpy()
        out.append(pd.DataFrame({"ko": d["kos"], "variable": name, "tiles": len(idx), "baseline": base, "added": added,
                                 "gain": added - base, "cruises": len(set(d["cr"][idx])), "cruises_improved": improved}))  # fmt: skip
        g = added - base
        print(f"  {name:45s} {len(idx):5,} tiles  median gain {np.median(g):+.4f}  KOs >= +0.05: {(g >= 0.05).sum():4d}  "
              f"of those R2 >= 0.3 and >= 10 cruises improved: {((g >= 0.05) & (added >= 0.3) & (improved >= 10)).sum():3d}", flush=True)
    t = pd.concat(out)
    if a.shard:
        t.round(4).to_csv(mm.OUT / f"added_fields_shard{a.shard[0]}.csv", index=False)
        return
    summarise(t)


def summarise(t):
    t.round(4).to_csv(mm.OUT / "added_fields.csv", index=False)
    s = t.assign(g05=t["gain"] >= 0.05, hit=(t["gain"] >= 0.05) & (t["added"] >= 0.3) & (t["cruises_improved"] >= 10)).groupby(
        ["variable", "tiles"], sort=False).agg(median_gain=("gain", "median"), kos_gain_005=("g05", "sum"), kos_heatmap_rule=("hit", "sum"))
    s.round(5).to_csv(mm.OUT / "added_fields_summary.csv")
    print(s.round(4).to_string())


if __name__ == "__main__":
    main()
