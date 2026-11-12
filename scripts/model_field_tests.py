#!/usr/bin/env python3
"""Ocean-model fields (GLORYS12, PISCES) in place of in-situ measurements (RESULTS.md R24).

A. Validation: model value at the sample (data/model_fields/sample_model_fields.csv,
   scripts/extract_model_fields.py) against the in-situ value from the manifests and
   the FAIRe sheet: temperature, salinity, phosphate, nitrate, chlorophyll.
B. Models on the figure KOs (as scripts/phosphate_models.py, one cruise held out):
     pca32 + temp~          in-situ temperature (R18)
     pca32 + Tm~            GLORYS temperature instead
     Pm~ > pca32            PISCES phosphate first, then 32 Granite PCs on the residuals
     Pm~ > pca32 + Tm~      the same + GLORYS temperature in the second stage
     P~ > pca32             in-situ phosphate first (R20); tiles with in-situ phosphate only
   Run on all tiles, and on the tiles with in-situ phosphate, so the PISCES and in-situ
   versions can be compared on the same rows.

Outputs in data/results/model_fields/: validation.csv, model_r2.csv (pooled, per KO,
model and row set), model_r2_cruise.csv, predictions_all.npz.

Usage:
    python scripts/model_field_tests.py
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
from phosphate_models import blocks, fit, phosphate  # noqa: E402
from tune_figure_kos import load, pooled_r2  # noqa: E402

RESULTS = REPO_ROOT / "data" / "results"
FIELDS = REPO_ROOT / "data" / "model_fields" / "sample_model_fields.csv"
MODELS = {  # name: (phosphate-first variable or None, second-stage temperature or None, phosphate model name in fit())
    "pca32 + temp~": ("insitu_T",),
    "pca32 + Tm~": ("model_T",),
    "Pm~ > pca32": ("model_P", None),
    "Pm~ > pca32 + Tm~": ("model_P", "model_T"),
    "P~ > pca32": ("insitu_P", None),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kos", type=Path, default=RESULTS / "ko_temperature_all_ostia_v2" / "ko_heatmap_rows.csv")
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings_sst_l4", REPO_ROOT / "data" / "faire" / "embeddings_sst_l4"])
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "model_fields")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    fields = pd.read_csv(FIELDS).set_index("sample_id")
    validate(fields, args.out_dir)

    kos = pd.read_csv(args.kos)["ko"].tolist()
    emb, temp, _, y, cruises, unit, rows_all = load(args, kos, return_rows=True)
    per_tile = lambda col: rows_all["sample_id"].map(fields[col]).groupby(rows_all["tile_key"]).mean().reindex(unit.index).to_numpy()  # noqa: E731
    vars_ = {
        "insitu_T": temp[:, 0],
        "model_T": per_tile("glorys_thetao"),
        "model_P": np.log10(np.clip(per_tile("pisces_po4"), 0.02, None)),
        "insitu_P": np.log10(np.clip(phosphate(unit, rows_all), 0.02, None)),
    }
    have_model = ~np.isnan(vars_["model_T"]) & ~np.isnan(vars_["model_P"])
    have_insitu = have_model & ~np.isnan(vars_["insitu_P"])
    print(f"\nTiles: {len(unit):,}; with GLORYS and PISCES: {have_model.sum():,}; also in-situ phosphate: {have_insitu.sum():,}")

    tables, per_cruise, saved = [], [], {}
    for row_set, mask in [("all tiles", have_model), ("tiles with in-situ phosphate", have_insitu)]:
        idx = np.flatnonzero(mask)
        models = [m for m in MODELS if row_set != "all tiles" or m != "P~ > pca32"]
        preds = {m: np.full_like(y, np.nan) for m in models}
        for cruise in np.unique(cruises[idx]):
            test, train = idx[cruises[idx] == cruise], idx[cruises[idx] != cruise]
            for m in models:
                preds[m][test] = run(m, emb, vars_, y, train, test, cruises)
                ss_res = ((y[test] - preds[m][test]) ** 2).sum(axis=0)
                ss_tot = ((y[test] - y[train].mean(axis=0)) ** 2).sum(axis=0)
                per_cruise.append(pd.DataFrame({"ko": kos, "model": m, "rows": row_set, "held_out": cruise, "r2": 1 - ss_res / ss_tot}))
        for m in models:
            tables.append(pd.DataFrame({"ko": kos, "model": m, "rows": row_set, "r2_pooled": pooled_r2(y[idx], preds[m][idx])}))
        if row_set == "all tiles":
            saved = preds
        print(f"  {row_set}: done", flush=True)

    table = pd.concat(tables)
    table.round(4).to_csv(args.out_dir / "model_r2.csv", index=False)
    pd.concat(per_cruise).round(4).to_csv(args.out_dir / "model_r2_cruise.csv", index=False)
    np.savez_compressed(args.out_dir / "predictions_all.npz", observed=y.astype(np.float32), ko=np.array(kos),
                        models=np.array(list(saved)), predicted=np.stack(list(saved.values())).astype(np.float32),
                        tile_path=unit["tile_path"].to_numpy().astype(str), section=cruises)  # fmt: skip

    info = pd.read_csv(args.kos).set_index("ko")
    pd.set_option("display.width", 220)
    for row_set in table["rows"].unique():
        wide = table[table["rows"] == row_set].pivot(index="ko", columns="model", values="r2_pooled")
        wide = wide[[m for m in MODELS if m in wide.columns]]
        print(f"\nMedian pooled R2, {row_set}, by heatmap block:")
        print(wide.groupby(info["block"]).median().round(3).to_string())
        show = ["K02037", "K02036", "K02040", "K02038", "K05382", "K14835", "K16090"]
        print(wide.loc[show].rename(index=info["gene"]).round(2).to_string())


def run(model, emb, v, y, train, test, cruises):
    spec = MODELS[model]
    if len(spec) == 1:  # pca32 + a temperature spline
        b = blocks(emb, v[spec[0]][:, None], v["model_P"][:, None], train, test)
        return fit("pca32 + temp~", b, y[train], cruises[train])
    first, second_t = spec
    t = v[second_t] if second_t else v["insitu_T"]
    b = blocks(emb, t[:, None], v[first][:, None], train, test)
    return fit("P~ > pca32 + T~" if second_t else "P~ > pca32", b, y[train], cruises[train])


def validate(fields: pd.DataFrame, out: Path) -> None:
    man = pd.concat([pd.read_csv(m, low_memory=False) for m in
                     [REPO_ROOT / "data/manifests/biogoship_samples.csv", REPO_ROOT / "data/faire/manifests/biogoship_samples.csv"]
                     if m.exists()]).drop_duplicates("sample_id").set_index("sample_id")  # fmt: skip
    faire = pd.read_csv(REPO_ROOT / "data/faire/faire_samples.csv", low_memory=False)
    fa = faire.set_index("samp_name")
    fa_orig = faire.dropna(subset=["original_sample_id"]).drop_duplicates("original_sample_id").set_index("original_sample_id")

    def insitu(col_manifest, col_faire):
        s = pd.Series(np.nan, index=fields.index)
        if col_manifest and col_manifest in man:
            s = fields.index.to_series().map(pd.to_numeric(man[col_manifest], errors="coerce"))
        if col_faire:
            v = pd.to_numeric(fa[col_faire], errors="coerce")
            s = s.fillna(fields.index.to_series().map(v[~v.index.duplicated()]))
            v2 = pd.to_numeric(fa_orig[col_faire], errors="coerce")
            s = s.fillna(fields.index.to_series().map(v2))
        return s

    pairs = [("temperature", "glorys_thetao", insitu("temperature_c", "temp"), False),
             ("salinity", "glorys_so", insitu(None, "salinity"), False),
             ("phosphate", "pisces_po4", insitu("phosphate_umol_kg", "phosphate"), True),
             ("nitrate", "pisces_no3", insitu("nitrate_nitrite_umol_l", "nitrate_plus_nitrite"), True),
             ("chlorophyll", "pisces_chl", insitu(None, "chlorophyll"), True)]
    rows = []
    for name, col, obs, logs in pairs:
        model = fields[col]
        ok = obs.notna() & model.notna()
        o, m = obs[ok], model[ok]
        if logs:
            lo, lm = np.log10(np.clip(o, 0.02, None)), np.log10(np.clip(m, 0.02, None))
        else:
            lo, lm = o, m
        rows.append({"variable": name, "model": col, "n": int(ok.sum()), "cruises": fields.loc[ok.index[ok], "section"].nunique(),
                     "r": np.corrcoef(lo, lm)[0, 1], "r_log" if logs else "r_raw": np.corrcoef(lo, lm)[0, 1],
                     "median_insitu": o.median(), "median_model": m.median(), "median_diff": (m - o).median(),
                     "r_within_cruise_median": pd.Series({c: np.corrcoef(lo[g.index], lm[g.index])[0, 1]
                                                          for c, g in fields.loc[ok.index[ok]].groupby("section") if len(g) > 10}).median()})  # fmt: skip
    table = pd.DataFrame(rows)
    table.round(3).to_csv(out / "validation.csv", index=False)
    print("Model vs in-situ at the samples (r on log10 for nutrients and chlorophyll):")
    print(table[["variable", "model", "n", "cruises", "r", "r_within_cruise_median", "median_insitu", "median_model", "median_diff"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
