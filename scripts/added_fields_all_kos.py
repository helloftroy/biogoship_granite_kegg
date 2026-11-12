#!/usr/bin/env python3
"""Ocean-model and history variables added to the best Granite model, on every KO (RESULTS.md R28).

R25 and R27 tested the added variables on the 54 figure KOs. Those were chosen for being
well predicted (the top of each heatmap block), which could tilt the test against new
inputs. Here the same "added" test runs on all target KOs, split by how well the
baseline predicts each KO:

    baseline   32 Granite PCs + in-situ temperature spline, ridge (R18)
    added      the same + a spline of the variable, one ridge
One cruise held out; pooled R2. Penalty by generalised cross-validation per KO (as the
baseline); spline: 8 knots, cubic, as elsewhere.

Variables: the sampling-day GLORYS / PISCES fields (R25) and the recent-history
summaries (R27), transformed as in scripts/model_field_first.py. Forecast-only
variables (from 2021-11) run on the tiles that have them, with their own baseline.
Controls, to measure how often a useless input "helps" by chance, on all tiles and on
the 637 forecast-period tiles:
    noise 1-3          standard normal values (three draws)
    shuffled phosphate PISCES phosphate permuted across tiles
Per KO, besides pooled R2: the number of held-out cruises where the added model has the
lower squared error ("cruises improved").

Outputs in data/results/model_fields/: added_all_kos.csv (KO x variable: tiles,
baseline, added, gain, cruises, cruises_improved), added_all_kos_summary.csv (per
variable and baseline stratum).

Usage:
    python scripts/added_fields_all_kos.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import SplineTransformer, StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship.tables import REPO_ROOT  # noqa: E402
from model_field_first import FIELDS, HISTORY, HISTORY_FIELDS, VARIABLES  # noqa: E402
from screen_ko_predictability import ALPHAS  # noqa: E402
from tune_figure_kos import load, pooled_r2  # noqa: E402

RESULTS = REPO_ROOT / "data" / "results"
STRATA = [(-np.inf, 0.0, "below 0"), (0.0, 0.1, "0 to 0.1"), (0.1, 0.3, "0.1 to 0.3"), (0.3, np.inf, "0.3 and above")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings_sst_l4", REPO_ROOT / "data" / "faire" / "embeddings_sst_l4"])
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "model_fields")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, _, y, cruises, unit, rows_all = load(args, kos, return_rows=True)
    fields = pd.read_csv(FIELDS).set_index("sample_id").join(pd.read_csv(HISTORY_FIELDS).set_index("sample_id")
                                                             .drop(columns=["section", "sample_time", "sample_lat", "sample_lon"]))  # fmt: skip
    fields["speed"] = np.hypot(fields["glorys_uo"], fields["glorys_vo"])
    per_tile = lambda col: rows_all["sample_id"].map(fields[col]).groupby(rows_all["tile_key"]).mean().reindex(unit.index).to_numpy()  # noqa: E731
    rng = np.random.default_rng(args.seed)
    variables = {}
    for name, (col, floor) in {**VARIABLES, **HISTORY}.items():
        v = per_tile(col)
        variables[name] = np.log10(np.clip(v, floor, None)) if floor is not None else v
    recent = ~np.isnan(variables["iron"])  # the forecast-period tiles
    for k in (1, 2, 3):
        variables[f"CONTROL: noise {k}"] = rng.standard_normal(len(unit))
        variables[f"CONTROL: noise {k} (637)"] = np.where(recent, rng.standard_normal(len(unit)), np.nan)
    variables["CONTROL: shuffled phosphate"] = rng.permutation(variables["phosphate"])
    variables["CONTROL: shuffled phosphate (637)"] = np.where(recent, variables["phosphate"], np.nan)
    variables["CONTROL: shuffled phosphate (637)"][recent] = rng.permutation(variables["phosphate"][recent])
    print(f"{len(unit):,} tiles, {len(kos):,} KOs, {len(variables)} variables", flush=True)

    baselines, out = {}, []
    for name, v in variables.items():
        idx = np.flatnonzero(~np.isnan(v))
        key = len(idx), idx[:5].tobytes()
        if key not in baselines:
            baselines[key] = cv(emb, temp[:, 0], None, y, cruises, idx)
        base, base_pred = baselines[key]
        added, pred = cv(emb, temp[:, 0], v, y, cruises, idx)
        sse = lambda pr: pd.DataFrame((y[idx] - pr[idx]) ** 2).groupby(cruises[idx]).sum()  # noqa: E731
        improved = (sse(pred) < sse(base_pred)).sum().to_numpy()
        out.append(pd.DataFrame({"ko": kos, "variable": name, "tiles": len(idx), "baseline": base, "added": added,
                                 "gain": added - base, "cruises": len(set(cruises[idx])), "cruises_improved": improved}))  # fmt: skip
        print(f"  {name:32s} {len(idx):5,} tiles  median gain {np.median(added - base):+.4f}  "
              f"KOs >= +0.05: {(added - base >= 0.05).sum():4d}", flush=True)
    table = pd.concat(out)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table.round(4).to_csv(args.out_dir / "added_all_kos.csv", index=False)

    rows = []
    for (name, tiles), part in table.groupby(["variable", "tiles"], sort=False):
        for lo, hi, label in STRATA:
            s = part[(part["baseline"] >= lo) & (part["baseline"] < hi)]
            rows.append({"variable": name, "tiles": tiles, "baseline_r2": label, "kos": len(s),
                         "median_gain": s["gain"].median(), "share_gain_0.05": (s["gain"] >= 0.05).mean(),
                         "kos_gain_0.05": int((s["gain"] >= 0.05).sum()), "kos_gain_0.10": int((s["gain"] >= 0.10).sum()),
                         "kos_gain_0.05_most_cruises": int(((s["gain"] >= 0.05) & (s["cruises_improved"] >= 0.75 * s["cruises"])).sum()),
                         "kos_now_useful": int(((s["gain"] >= 0.05) & (s["added"] >= 0.3)).sum()),
                         "kos_crossing_0.3": int(((s["baseline"] < 0.3) & (s["added"] >= 0.3)).sum()),
                         "max_gain": s["gain"].max()})  # fmt: skip
    summary = pd.DataFrame(rows)
    summary.round(4).to_csv(args.out_dir / "added_all_kos_summary.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    wide = summary.pivot_table(index=["variable", "tiles"], columns="baseline_r2", values="kos_gain_0.05", sort=False)
    print("\nKOs gaining >= 0.05, by baseline R2 stratum:")
    print(wide[[s[2] for s in STRATA]].to_string())
    totals = summary.groupby(["variable", "tiles"], sort=False)[["kos_gain_0.05", "kos_gain_0.05_most_cruises", "kos_now_useful",
                                                                  "kos_crossing_0.3"]].sum()
    print("\nAll strata: KOs gaining >= 0.05; of those, improved in >= 75% of cruises; gaining and R2 >= 0.3 after; crossing 0.3:")
    print(totals.to_string())
    wide = summary.pivot_table(index=["variable", "tiles"], columns="baseline_r2", values="median_gain", sort=False)
    print("\nMedian gain, by baseline R2 stratum:")
    print(wide[[s[2] for s in STRATA]].round(4).to_string())
    print("\nKOs per stratum:", summary.drop_duplicates(["tiles", "baseline_r2"]).set_index(["tiles", "baseline_r2"])["kos"].to_dict())


def cv(emb, temp, v, y, cruises, idx):
    """Pooled R2 per KO of 32 PCs + temperature spline (+ the variable's spline), one cruise held out."""
    pred = np.full_like(y, np.nan)
    for cruise in np.unique(cruises[idx]):
        test, train = idx[cruises[idx] == cruise], idx[cruises[idx] != cruise]
        s = StandardScaler().fit(emb[train])
        pca = PCA(n_components=32, random_state=0).fit(s.transform(emb[train]))
        p_tr, p_te = pca.transform(s.transform(emb[train])), pca.transform(s.transform(emb[test]))
        s2 = StandardScaler().fit(p_tr)
        parts_tr, parts_te = [s2.transform(p_tr)], [s2.transform(p_te)]
        for x in ([temp] if v is None else [temp, v]):
            sx = StandardScaler().fit(x[train, None])
            sp = SplineTransformer(n_knots=8, degree=3, extrapolation="linear").fit(sx.transform(x[train, None]))
            parts_tr.append(sp.transform(sx.transform(x[train, None])))
            parts_te.append(sp.transform(sx.transform(x[test, None])))
        model = RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(np.hstack(parts_tr), y[train])
        pred[test] = model.predict(np.hstack(parts_te))
    return pooled_r2(y[idx], pred[idx]), pred


if __name__ == "__main__":
    main()
