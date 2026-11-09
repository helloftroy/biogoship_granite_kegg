#!/usr/bin/env python3
"""Granite models with in-situ phosphate, on the figure KOs (RESULTS.md R19, R20).

Samples: the R17/R18 tiles that have an in-situ phosphate value (manifest
phosphate_umol_kg, else the FAIRe sheet; averaged over merged samples): 451 tiles,
12 cruises. One cruise held out at a time, as elsewhere.

Models (ridge; "~" = cubic spline of the input, 8 knots):
  emb               Granite embedding (penalty by leave-one-sample-out, as R17)
  pca32 + temp~     32 PCs + temperature spline (R18's leading model)
  P~ > pca32        phosphate spline first, then 32 PCs fitted to its residuals
  T~P~ > pca32      temperature + phosphate splines first, then 32 PCs on the residuals
  P~ > pca32 + T~   phosphate spline first, then 32 PCs + temperature spline on the residuals
The second stage chooses its penalty by leave-one-cruise-out within the training cruises.
Phosphate enters as log10(max(PO4, 0.02 umol/kg)).

Outputs in --out-dir (default data/results/phosphate_models/): model_r2.csv (pooled),
model_r2_cruise.csv, predictions.npz (rows = the R17 tiles; NaN for tiles without
phosphate), phosphate.csv (per tile).

Usage:
    python scripts/phosphate_models.py
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
from screen_ko_predictability import ALPHAS  # noqa: E402
from tune_figure_kos import load, pooled_r2, ridge_grouped  # noqa: E402

RESULTS = REPO_ROOT / "data" / "results"
MODELS = ["emb", "pca32 + temp~", "P~ > pca32", "T~P~ > pca32", "P~ > pca32 + T~"]


def phosphate(unit: pd.DataFrame, rows_all: pd.DataFrame) -> np.ndarray:
    """In-situ phosphate per tile (umol/kg): manifest value, else the FAIRe sheet; mean over merged samples."""
    manifests = [REPO_ROOT / "data" / "manifests" / "biogoship_samples.csv",
                 REPO_ROOT / "data" / "faire" / "manifests" / "biogoship_samples.csv"]  # fmt: skip
    man = pd.concat([pd.read_csv(m, low_memory=False) for m in manifests if m.exists()]).drop_duplicates("sample_id")
    value = pd.to_numeric(man.set_index("sample_id")["phosphate_umol_kg"], errors="coerce")
    faire = pd.read_csv(REPO_ROOT / "data" / "faire" / "faire_samples.csv", low_memory=False)
    po4 = pd.to_numeric(faire["phosphate"], errors="coerce")
    extra = pd.concat([po4.set_axis(faire["samp_name"]),
                       po4[faire["original_sample_id"].notna().to_numpy()].set_axis(faire["original_sample_id"].dropna())])  # fmt: skip
    extra = extra[~extra.index.duplicated()]
    per_sample = rows_all["sample_id"].map(value).fillna(rows_all["sample_id"].map(extra))
    return per_sample.groupby(rows_all["tile_key"]).mean().reindex(unit.index).to_numpy()


def blocks(emb, temp, logp, train, test):
    def spline(x):
        s = StandardScaler().fit(x[train])
        sp = SplineTransformer(n_knots=8, degree=3, extrapolation="linear").fit(s.transform(x[train]))
        return sp.transform(s.transform(x[train])), sp.transform(s.transform(x[test]))

    s = StandardScaler().fit(emb[train])
    e_tr, e_te = s.transform(emb[train]), s.transform(emb[test])
    pca = PCA(n_components=32, random_state=0).fit(e_tr)
    s2 = StandardScaler().fit(pca.transform(e_tr))
    return {"E": (e_tr, e_te), "P32": (s2.transform(pca.transform(e_tr)), s2.transform(pca.transform(e_te))),
            "T~": spline(temp), "P~": spline(logp)}  # fmt: skip


def fit(model, b, y_tr, groups):
    def cat(*names):
        return np.hstack([b[n][0] for n in names]), np.hstack([b[n][1] for n in names])

    def gcv(names, target):
        tr, te = cat(*names)
        m = RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(tr, target)
        return m.predict(tr), m.predict(te)

    if model == "emb":
        return gcv(["E"], y_tr)[1]
    if model == "pca32 + temp~":
        return gcv(["P32", "T~"], y_tr)[1]
    first, second = {"P~ > pca32": (["P~"], ["P32"]), "T~P~ > pca32": (["T~", "P~"], ["P32"]),
                      "P~ > pca32 + T~": (["P~"], ["P32", "T~"])}[model]  # fmt: skip
    fitted, pred = gcv(first, y_tr)
    tr, te = cat(*second)
    return pred + ridge_grouped(tr, y_tr - fitted, te, groups)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kos", type=Path, default=RESULTS / "ko_temperature_all_ostia_v2" / "ko_heatmap_rows.csv")
    parser.add_argument("--targets", type=Path, default=REPO_ROOT / "data" / "targets" / "ko_clr.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings_sst_l4", REPO_ROOT / "data" / "faire" / "embeddings_sst_l4"])
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "phosphate_models")
    args = parser.parse_args()

    kos = pd.read_csv(args.kos)["ko"].tolist()
    emb, temp, _, y, cruises, unit, rows_all = load(args, kos, return_rows=True)
    po4 = phosphate(unit, rows_all)
    logp = np.log10(np.clip(po4, 0.02, None))[:, None]
    have = np.flatnonzero(~np.isnan(po4))
    print(f"Tiles with phosphate: {len(have):,} of {len(unit):,}; cruises: {len(set(cruises[have]))}", flush=True)

    preds = {m: np.full_like(y, np.nan) for m in MODELS}
    per_cruise = []
    for cruise in np.unique(cruises[have]):
        test, train = have[cruises[have] == cruise], have[cruises[have] != cruise]
        b = blocks(emb, temp, logp, train, test)
        for m in MODELS:
            preds[m][test] = fit(m, b, y[train], cruises[train])
            ss_res = ((y[test] - preds[m][test]) ** 2).sum(axis=0)
            ss_tot = ((y[test] - y[train].mean(axis=0)) ** 2).sum(axis=0)
            per_cruise.append(pd.DataFrame({"ko": kos, "model": m, "held_out": cruise, "r2": 1 - ss_res / ss_tot}))
        print(f"  {cruise}", flush=True)

    table = pd.concat([pd.DataFrame({"ko": kos, "model": m, "r2_pooled": pooled_r2(y[have], preds[m][have])}) for m in MODELS])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table.round(4).to_csv(args.out_dir / "model_r2.csv", index=False)
    pd.concat(per_cruise).round(4).to_csv(args.out_dir / "model_r2_cruise.csv", index=False)
    pd.DataFrame({"tile_path": unit["tile_path"].to_numpy(), "section": cruises, "phosphate_umol_kg": po4}).to_csv(
        args.out_dir / "phosphate.csv", index=False)
    np.savez_compressed(args.out_dir / "predictions.npz", observed=y.astype(np.float32), ko=np.array(kos),
                        models=np.array(MODELS), predicted=np.stack([preds[m] for m in MODELS]).astype(np.float32),
                        tile_path=unit["tile_path"].to_numpy().astype(str), section=cruises)  # fmt: skip

    cat = pd.read_csv(args.kos).set_index("ko")["category"]
    wide = table.pivot(index="ko", columns="model", values="r2_pooled")[MODELS]
    pd.set_option("display.width", 200)
    print("\nMedian pooled R2 (tiles with phosphate) by heatmap category:")
    print(wide.groupby(cat).median().round(3).to_string())
    names = pd.read_csv(args.kos).set_index("ko")["gene"]
    show = ["K02037", "K02036", "K02040", "K02038", "K05382", "K14835", "K01425", "K16090"]
    print("\nPer KO:")
    print(wide.loc[show].rename(index=names).round(2).to_string())


if __name__ == "__main__":
    main()
