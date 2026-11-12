#!/usr/bin/env python3
"""Each GLORYS / PISCES variable first, then the current best Granite model (RESULTS.md R25).

For every variable V at the samples (data/model_fields/sample_model_fields.csv):
    V~ > pca32 + temp~   spline of V first (ridge), then 32 Granite PCs + in-situ temperature
                         spline on the residuals (penalty by leave-one-cruise-out), as R20's
                         "P~ > pca32 + T~" with V in place of phosphate
    pca32 + temp~ + V~   V's spline added to the 32 PCs + temperature spline, one ridge
against the current best model, pca32 + temp~ (R18), on the same tiles. One cruise held out.
Positive variables (nutrients, chlorophyll, production, mixed layer, current speed, iron,
phytoplankton carbon) enter as log10(max(V, small floor)); the rest as they are.
Variables only in the forecast runs (from 2021-11: fe, phyc, ph, dissic, talk, spco2) are
tested on the tiles that have them, with their own baseline.

Outputs in data/results/model_fields/: first_r2.csv (KO x variable: baseline, V first,
V added; pooled R2, and the tiles used).

Usage:
    python scripts/model_field_first.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship.tables import REPO_ROOT  # noqa: E402
from phosphate_models import blocks, fit  # noqa: E402
from screen_ko_predictability import ALPHAS  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402
from tune_figure_kos import load, pooled_r2  # noqa: E402

RESULTS = REPO_ROOT / "data" / "results"
FIELDS = REPO_ROOT / "data" / "model_fields" / "sample_model_fields.csv"
# variable: (column or "speed", log floor or None)
VARIABLES = {
    "salinity": ("glorys_so", None), "mixed layer": ("glorys_mlotst", 1.0), "sea level": ("glorys_zos", None),
    "current speed": ("speed", 0.005), "GLORYS temperature": ("glorys_thetao", None),
    "phosphate": ("pisces_po4", 0.001), "nitrate": ("pisces_no3", 0.001), "silicate": ("pisces_si", 0.01),
    "oxygen": ("pisces_o2", None), "chlorophyll": ("pisces_chl", 0.001), "primary production": ("pisces_nppv", 0.01),
    "iron": ("pisces_fe", 1e-5), "phytoplankton carbon": ("pisces_phyc", 0.001), "pH": ("pisces_ph", None),
    "DIC": ("pisces_dissic", None), "alkalinity": ("pisces_talk", None), "pCO2": ("pisces_spco2", None),
}  # fmt: skip


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kos", type=Path, default=RESULTS / "ko_temperature_all_ostia_v2" / "ko_heatmap_rows.csv")
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings_sst_l4", REPO_ROOT / "data" / "faire" / "embeddings_sst_l4"])
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "model_fields")
    args = parser.parse_args()

    fields = pd.read_csv(FIELDS).set_index("sample_id")
    fields["speed"] = np.hypot(fields["glorys_uo"], fields["glorys_vo"])
    kos = pd.read_csv(args.kos)["ko"].tolist()
    emb, temp, _, y, cruises, unit, rows_all = load(args, kos, return_rows=True)

    baselines, out = {}, []
    for name, (col, floor) in VARIABLES.items():
        v = rows_all["sample_id"].map(fields[col]).groupby(rows_all["tile_key"]).mean().reindex(unit.index).to_numpy()
        if floor is not None:
            v = np.log10(np.clip(v, floor, None))
        idx = np.flatnonzero(~np.isnan(v))
        key = tuple(idx)
        if key not in baselines:
            baselines[key] = cv(emb, temp[:, 0], v, y, cruises, idx, "pca32 + temp~")
        base = baselines[key]
        first = cv(emb, temp[:, 0], v, y, cruises, idx, "P~ > pca32 + T~")
        added = cv(emb, temp[:, 0], v, y, cruises, idx, "added")
        out.append(pd.DataFrame({"ko": kos, "variable": name, "tiles": len(idx), "cruises": len(set(cruises[idx])),
                                 "baseline": base, "first": first, "gain": first - base,
                                 "added": added, "gain_added": added - base}))  # fmt: skip
        print(f"  {name:22s} {len(idx):5,} tiles  first: median gain {np.median(first - base):+.3f}, "
              f"{(first - base >= 0.05).sum()} KOs >= +0.05;  added: median gain {np.median(added - base):+.3f}, "
              f"{(added - base >= 0.05).sum()} KOs >= +0.05", flush=True)
    table = pd.concat(out)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table.round(4).to_csv(args.out_dir / "first_r2.csv", index=False)

    info = pd.read_csv(args.kos).set_index("ko")
    table = table.join(info[["gene", "block"]], on="ko")
    pd.set_option("display.width", 220)
    print("\nKOs that gain >= 0.05 pooled R2 with a variable first:")
    hits = table[(table["gain"] >= 0.05) | (table["gain_added"] >= 0.05)].sort_values(["variable", "gain"], ascending=[True, False])
    print(hits[["variable", "gene", "block", "tiles", "baseline", "first", "gain", "added", "gain_added"]].round(2).to_string(index=False))


def cv(emb, temp, v, y, cruises, idx, model):
    pred = np.full_like(y, np.nan)
    for cruise in np.unique(cruises[idx]):
        test, train = idx[cruises[idx] == cruise], idx[cruises[idx] != cruise]
        b = blocks(emb, temp[:, None], v[:, None], train, test)
        if model == "added":
            tr = np.hstack([b["P32"][0], b["T~"][0], b["P~"][0]])
            te = np.hstack([b["P32"][1], b["T~"][1], b["P~"][1]])
            pred[test] = RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(tr, y[train]).predict(te)
        else:
            pred[test] = fit(model, b, y[train], cruises[train])
    return pooled_r2(y[idx], pred[idx])


if __name__ == "__main__":
    main()
