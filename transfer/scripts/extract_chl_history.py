#!/usr/bin/env python3
"""Satellite chlorophyll over the 4 weeks before each sample: is a bloom growing or decaying? (RESULTS.md R35)

Gap-free daily L4 chlorophyll (GlobColour multi-sensor, 4 km, Copernicus Marine
OCEANCOLOUR_GLO_BGC_L4_*_009_10x), the first dataset in DATASETS whose time range
covers the window: multi-year, then interim, then near-real-time. Per sample (sampling
position and UTC day, as scripts/extract_history_fields.py), the box mean of CHL within
+-BOX deg (about 7 x 7 cells) on each of the 29 days from day -28 to the sampling day:

  chl_day0          chlorophyll on the sampling day (mg m-3)
  chl_mean28        mean of days -28..-1
  chl_max28         maximum of days -28..0
  chl_slope14       slope of log10 chlorophyll over days -14..0 (log10 units per 10 days):
                    > 0 growing, < 0 decaying
  chl_slope28       the same over days -28..0
  chl_peak_ratio28  log10(max over days -28..0 / day 0): 0 = at the peak now; > 0 = past a peak
  chl_days_since_max  days since the 28-day maximum (0 = today)
  chl_cv28          coefficient of variation of the 29 daily values
  chl_n             days with a value (gap-free products: normally 29)
The 29 daily box means are kept too (chl_d-28 ... chl_d0), for other summaries later.

Output transfer/results/chl_history/sample_chl_history.csv and provenance_chl_history.json.
Resumable: values are appended to transfer/results/chl_history/partial_chl.csv as computed;
a re-run skips samples already there.

Needs internet and the Copernicus Marine login: run transfer/cluster/chl_history.sbatch.

Usage:
    python transfer/scripts/extract_chl_history.py --list-datasets
    python transfer/scripts/extract_chl_history.py --workers 4
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time as clock
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from extract_history_fields import describe, open_dataset, subset  # noqa: E402
from extract_model_fields import load_samples  # noqa: E402

DATASETS = ["cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D",
            "cmems_obs-oc_glo_bgc-plankton_myint_l4-gapfree-multi-4km_P1D",
            "cmems_obs-oc_glo_bgc-plankton_nrt_l4-gapfree-multi-4km_P1D"]  # fmt: skip
VARIABLE = "CHL"
BOX, DAYS = 0.125, 28
OUT = REPO / "transfer" / "results" / "chl_history"
LOCK = threading.Lock()


def main() -> None:
    args = parse_args()
    if args.list_datasets:
        for d in DATASETS:
            print(describe(d, [VARIABLE]), flush=True)
        return
    samples = load_samples(args.manifests)
    print(f"{len(samples):,} samples, {samples['sample_time'].min():%Y-%m-%d} to {samples['sample_time'].max():%Y-%m-%d}", flush=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    path = args.out_dir / "partial_chl.csv"
    done = set(pd.read_csv(path, on_bad_lines="skip")["sample_id"]) if path.exists() else set()
    todo = samples[~samples["sample_id"].isin(done)]
    print(f"{len(done):,} samples already done, {len(todo):,} to do", flush=True)
    opened = []
    for dataset_id in DATASETS:
        if todo.empty:
            break
        try:
            ds = open_dataset(dataset_id)
        except Exception as exc:
            print(f"  {dataset_id}: not opened ({str(exc).splitlines()[0][:160]})", flush=True)
            opened.append({"dataset_id": dataset_id, "error": str(exc)[:500]})
            continue
        t0, t1 = pd.to_datetime(ds["time"].values[[0, -1]])
        opened.append({"dataset_id": dataset_id, "time_start": str(t0), "time_end": str(t1)})
        day = todo["sample_time"].dt.floor("D")
        covered = todo[(day - pd.Timedelta(days=DAYS) >= t0) & (day <= t1)]
        print(f"  {dataset_id}: {t0:%Y-%m-%d} to {t1:%Y-%m-%d}; {len(covered):,} samples covered", flush=True)
        started, count = clock.time(), [0]

        def work(row):
            values = compute(ds, row)
            line = pd.DataFrame([{"sample_id": row["sample_id"], "dataset_id": dataset_id, **values}])
            with LOCK:
                line.to_csv(path, mode="a", header=not path.exists(), index=False)
                count[0] += 1
                if count[0] % 100 == 0:
                    print(f"    {count[0]:,}/{len(covered):,} ({clock.time() - started:5.0f} s)", flush=True)

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(work, [r for _, r in covered.iterrows()]))
        todo = todo[~todo["sample_id"].isin(covered["sample_id"])]
        ds.close()
    if not todo.empty:
        print(f"  {len(todo):,} samples not covered by any dataset", flush=True)
    write_outputs(samples, path, opened, args)


def compute(ds, row) -> dict:
    t = row["sample_time"].floor("D")
    empty = {"chl_n": 0}
    for attempt in range(1, 6):
        try:
            da = subset(ds[VARIABLE], row["sample_lat"], row["sample_lon"], t - pd.Timedelta(days=DAYS),
                        t + pd.Timedelta(hours=23, minutes=59), BOX)
            break
        except Exception as exc:  # network hiccup
            if attempt == 5:
                print(f"    {row['sample_id']}: failed ({str(exc)[:100]})", flush=True)
                return empty
            clock.sleep(5 * attempt)
    series = pd.Series(da.mean(dim=["latitude", "longitude"], skipna=True).values,
                       index=(pd.to_datetime(da["time"].values).floor("D") - t).days)
    series = series.groupby(level=0).mean().reindex(range(-DAYS, 1))
    out = {f"chl_d{d}": series[d] for d in series.index}
    ok = series.notna()
    out["chl_n"] = int(ok.sum())
    if ok.sum() < 10 or pd.isna(series[0]):
        return out
    logc = np.log10(series.clip(lower=1e-3))

    def slope(days):
        s = logc[logc.index >= -days].dropna()
        return np.polyfit(s.index, s.values, 1)[0] * 10 if len(s) >= days // 2 else np.nan

    mx = series.max()
    out.update({"chl_day0": series[0], "chl_mean28": series[series.index < 0].mean(), "chl_max28": mx,
                "chl_slope14": slope(14), "chl_slope28": slope(28), "chl_peak_ratio28": np.log10(mx / series[0]),
                "chl_days_since_max": int(-series.idxmax()), "chl_cv28": series.std() / series.mean()})  # fmt: skip
    return out


def write_outputs(samples, path, opened, args) -> None:
    wide = samples[["sample_id", "section", "sample_time", "sample_lat", "sample_lon"]].set_index("sample_id")
    if path.exists():
        part = pd.read_csv(path, on_bad_lines="skip").drop_duplicates("sample_id", keep="last").set_index("sample_id")
        wide = wide.join(part)
    wide.reset_index().to_csv(args.out_dir / "sample_chl_history.csv", index=False)
    (args.out_dir / "provenance_chl_history.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "box_deg": BOX, "days": DAYS,
        "variable": VARIABLE, "datasets": opened, "n_samples": len(wide)}, indent=2))  # fmt: skip
    print(f"\nWrote {args.out_dir / 'sample_chl_history.csv'}")
    for f in ["chl_day0", "chl_slope14", "chl_slope28", "chl_peak_ratio28", "chl_days_since_max"]:
        if f in wide:
            print(f"  {f:20s} {int(wide[f].notna().sum()):,} values, median {wide[f].median():.3g}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifests", type=Path, nargs="+",
                   default=[REPO / "data" / "manifests" / "biogoship_tile_manifest.csv",
                            REPO / "data" / "faire" / "manifests" / "biogoship_tile_manifest.csv"])  # fmt: skip
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--out-dir", type=Path, default=OUT)
    p.add_argument("--list-datasets", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    main()
