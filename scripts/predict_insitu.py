#!/usr/bin/env python3
"""Can Granite predict the in-situ measurements better than GLORYS / PISCES reproduce them? (RESULTS.md R26)

Targets (one at a time, samples that have the value; from the manifests and the FAIRe
sheet as in scripts/model_field_tests.py): temperature, salinity, log10 phosphate,
log10 nitrate + nitrite, log10 chlorophyll.

Granite and other inputs: scripts/test_ko_temperature.py with --targets set to the
variable (one cruise held out, pooled R2; learners ridge, spline, knn; the best learner
per input is reported). For temperature, inputs that contain in-situ temperature are
left out.

Ocean model, on the same tiles: the GLORYS / PISCES value at the sample
  as is          R2 of the model value itself (pooled, against the mean)
  calibrated     a straight line from model value to measurement, fitted on the
                 training cruises, one cruise held out (removes PISCES's offset)
  r2             squared correlation (the ceiling of any linear calibration)

Outputs in data/results/insitu/: targets_<variable>.csv, <variable>/ (test_ko_temperature
output), summary.csv.

Usage:
    python scripts/predict_insitu.py
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

OUT = REPO_ROOT / "data" / "results" / "insitu"
FIELDS = REPO_ROOT / "data" / "model_fields" / "sample_model_fields.csv"
# name: (manifest column, FAIRe column, model column, log10?)
VARIABLES = {
    "temperature": ("temperature_c", "temp", "glorys_thetao", False),
    "salinity": (None, "salinity", "glorys_so", False),
    "phosphate": ("phosphate_umol_kg", "phosphate", "pisces_po4", True),
    "nitrate": ("nitrate_nitrite_umol_l", "nitrate_plus_nitrite", "pisces_no3", True),
    "chlorophyll": (None, "chlorophyll", "pisces_chl", True),
}
FLOOR = 0.01  # log10(max(value, FLOOR)) for nutrients and chlorophyll


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variables", nargs="+", default=list(VARIABLES))
    parser.add_argument("--skip-fit", action="store_true", help="Only summarise existing runs.")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    fields = pd.read_csv(FIELDS).set_index("sample_id")
    insitu = measurements(fields.index)

    rows = []
    for name in args.variables:
        _, _, model_col, logs = VARIABLES[name]
        value = insitu[name]
        if logs:
            value = np.log10(value.clip(lower=FLOOR))
        target = pd.DataFrame({name: value, f"{name}_copy": value}).dropna()  # two columns: the script expects a table
        path = OUT / f"targets_{name}.csv"
        target.rename_axis("sample_key").to_csv(path)
        sets = ["embedding", "satellite", "metadata"] + ([] if name == "temperature" else ["temp", "pca32_temp"])
        if not args.skip_fit:
            subprocess.run([sys.executable, str(REPO_ROOT / "scripts" / "test_ko_temperature.py"), "--kos", "all",
                            "--targets", str(path), "--learners", "ridge", "spline", "knn", "--feature-sets", *sets,
                            "--embeddings-dirs", str(REPO_ROOT / "data/embeddings_sst_l4"),
                            str(REPO_ROOT / "data/faire/embeddings_sst_l4"), "--save-predictions",
                            "--out-dir", str(OUT / name)], check=True, stdout=subprocess.DEVNULL)  # fmt: skip
        scores = pd.read_csv(OUT / name / "ko_temperature.csv")
        scores = scores[scores["ko"] == name]
        best = scores.loc[scores.groupby("feature_set")["r2_pooled"].idxmax()].set_index("feature_set")

        tiles = pd.read_csv(OUT / name / "prediction_rows.csv")
        z = np.load(OUT / name / f"predictions_embedding_{best.loc['embedding', 'learner']}.npz", allow_pickle=True)
        y = z["observed"][:, list(z["ko"]).index(name)].astype(float)
        m = fields[model_col].reindex(tiles["sample_id"]).to_numpy()
        if logs:
            m = np.log10(np.clip(m, FLOOR, None))
        ok = ~np.isnan(m)
        cruise = tiles["section"].to_numpy()
        calibrated = np.full_like(y, np.nan)
        for c in np.unique(cruise[ok]):
            tr, te = ok & (cruise != c), ok & (cruise == c)
            slope, icept = np.polyfit(m[tr], y[tr], 1)
            calibrated[te] = icept + slope * m[te]
        sst = ((y[ok] - y[ok].mean()) ** 2).sum()
        row = {"variable": name, "tiles": len(y), "cruises": len(set(cruise)),
               "model_as_is": 1 - ((y[ok] - m[ok]) ** 2).sum() / sst,
               "model_calibrated": 1 - ((y[ok] - calibrated[ok]) ** 2).sum() / sst,
               "model_r2": np.corrcoef(y[ok], m[ok])[0, 1] ** 2}  # fmt: skip
        for f, r in best.iterrows():
            row[f"{f}"] = r["r2_pooled"]
            row[f"{f}_learner"] = r["learner"]
        rows.append(row)
        print(f"{name}: done", flush=True)

    summary = pd.DataFrame(rows)
    summary.round(3).to_csv(OUT / "summary.csv", index=False)
    show = ["variable", "tiles", "cruises", "model_as_is", "model_calibrated", "model_r2", "embedding", "pca32_temp",
            "satellite", "temp", "metadata"]  # fmt: skip
    pd.set_option("display.width", 220)
    print(summary[[c for c in show if c in summary]].round(2).to_string(index=False))


def measurements(index: pd.Index) -> pd.DataFrame:
    """In-situ values per sample: manifest column first, then the FAIRe sheet (samp_name, then original_sample_id)."""
    manifests = [REPO_ROOT / "data/manifests/biogoship_samples.csv", REPO_ROOT / "data/faire/manifests/biogoship_samples.csv"]
    man = pd.concat([pd.read_csv(m, low_memory=False) for m in manifests if m.exists()]).drop_duplicates("sample_id").set_index("sample_id")
    faire = pd.read_csv(REPO_ROOT / "data/faire/faire_samples.csv", low_memory=False)
    by_name = faire.drop_duplicates("samp_name").set_index("samp_name")
    by_orig = faire.dropna(subset=["original_sample_id"]).drop_duplicates("original_sample_id").set_index("original_sample_id")
    ids = index.to_series()
    out = {}
    for name, (col_m, col_f, _, _) in VARIABLES.items():
        s = pd.Series(np.nan, index=index)
        if col_m:
            s = ids.map(pd.to_numeric(man[col_m], errors="coerce"))
        s = s.fillna(ids.map(pd.to_numeric(by_name[col_f], errors="coerce")))
        s = s.fillna(ids.map(pd.to_numeric(by_orig[col_f], errors="coerce")))
        out[name] = s
    return pd.DataFrame(out)


if __name__ == "__main__":
    main()
