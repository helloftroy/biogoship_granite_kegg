#!/usr/bin/env python3
"""Recent ocean history around each sample, before the sampling day (RESULTS.md R27).

A single day's model value misses what happened just before sampling. Per sample (the
sampling position and UTC day, as scripts/extract_model_fields.py), from Copernicus
Marine:

  mld_max30      mixing: maximum mixed-layer depth (m) over the 30 days up to and including
                 the sampling day. GLORYS12 mlotst, mean of the cells within +-BOX deg.
  mld_now        the same on the sampling day (for reference)
  sst_change14   cooling or upwelling: SST on the sampling day minus the mean of the 14
                 days before (K). OSTIA L4 (reprocessed, then near-real-time), +-BOX deg.
  wind_max14     storms: peak wind speed (m/s) over the 14 days up to the sampling day,
                 hourly L4 winds (multi-year, then near-real-time), +-BOX deg.
  wind_mean14    mean wind speed over the same 14 days
  ssh_gradient   fronts: mean magnitude of the sea-level gradient (m per 100 km) on the
                 sampling day within +-FRONT_BOX deg. GLORYS12 zos.
  eke30          eddies: mean eddy kinetic energy 0.5 (u'^2 + v'^2) (m^2/s^2) of the
                 geostrophic velocity anomalies over the 30 days up to the sampling day,
                 within +-FRONT_BOX deg. DUACS L4 altimetry (multi-year, interim, then
                 near-real-time).

Each feature takes the first dataset in its list whose time range covers the window;
the dataset used is recorded per value. A dataset that cannot be opened is skipped and
reported (dataset ids change between catalogue versions; --list-datasets checks them).
Boxes that cross the date line are handled. The tile itself is 12.6 km; the boxes are
model-grid sized (BOX 0.125 deg = about 3 GLORYS cells across; FRONT_BOX 0.25 deg).

Outputs (data/model_fields/):
  sample_history_fields.csv   one row per sample: sample columns, each feature, and
                              <feature>_dataset, <feature>_n (time steps used)
  provenance_history.json     datasets opened (time ranges, errors), settings
Resumable: each value is appended to data/model_fields/partial_history/<feature>.csv as
it is computed; a re-run skips (sample, feature) pairs already there.

Needs the Copernicus Marine login and internet: run cluster/history_fields.sbatch.

Usage:
    python scripts/extract_history_fields.py --list-datasets
    python scripts/extract_history_fields.py --workers 4
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship.tables import REPO_ROOT  # noqa: E402
from extract_model_fields import load_samples  # noqa: E402

BOX, FRONT_BOX = 0.125, 0.25
GLORYS = ["cmems_mod_glo_phy_my_0.083deg_P1D-m", "cmems_mod_glo_phy_anfc_0.083deg_P1D-m"]
OSTIA = ["METOFFICE-GLO-SST-L4-REP-OBS-SST", "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2"]
WIND = ["cmems_obs-wind_glo_phy_my_l4_0.125deg_PT1H", "cmems_obs-wind_glo_phy_nrt_l4_0.125deg_PT1H"]
DUACS = ["cmems_obs-sl_glo_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D",
         "cmems_obs-sl_glo_phy-ssh_myint_allsat-l4-duacs-0.125deg_P1D",
         "cmems_obs-sl_glo_phy-ssh_nrt_allsat-l4-duacs-0.125deg_P1D"]  # fmt: skip
# feature group: (datasets, variables, days before the sampling day, box half-width)
GROUPS = {
    "mld": (GLORYS, ["mlotst"], 30, BOX),
    "sst": (OSTIA, ["analysed_sst"], 14, BOX),
    "wind": (WIND, ["eastward_wind", "northward_wind"], 14, BOX),
    "front": (GLORYS, ["zos"], 0, FRONT_BOX),
    "eke": (DUACS, ["ugosa", "vgosa"], 30, FRONT_BOX),
}
FEATURES = {"mld": ["mld_max30", "mld_now"], "sst": ["sst_change14"], "wind": ["wind_max14", "wind_mean14"],
            "front": ["ssh_gradient"], "eke": ["eke30"]}  # fmt: skip
LOCK = threading.Lock()


def main() -> None:
    args = parse_args()
    if args.list_datasets:
        for group, (datasets, variables, _, _) in GROUPS.items():
            for dataset_id in datasets:
                print(f"{group:6s} {describe(dataset_id, variables)}", flush=True)
        return
    samples = load_samples(args.manifests)
    print(f"{len(samples):,} samples, {samples['sample_time'].min():%Y-%m-%d} to {samples['sample_time'].max():%Y-%m-%d}", flush=True)
    partial = args.out_dir / "partial_history"
    partial.mkdir(parents=True, exist_ok=True)
    opened = []
    for group in args.groups:
        datasets, variables, days, box = GROUPS[group]
        path = partial / f"{group}.csv"
        done = set(pd.read_csv(path, on_bad_lines="skip")["sample_id"]) if path.exists() else set()
        todo = samples[~samples["sample_id"].isin(done)]
        print(f"\n{group}: {len(done):,} samples already done, {len(todo):,} to do", flush=True)
        for dataset_id in datasets:
            if todo.empty:
                break
            try:
                ds = open_dataset(dataset_id)
            except Exception as exc:
                print(f"  {dataset_id}: not opened ({str(exc).splitlines()[0][:160]})", flush=True)
                opened.append({"group": group, "dataset_id": dataset_id, "error": str(exc)[:500]})
                continue
            t0, t1 = pd.to_datetime(ds["time"].values[[0, -1]])
            opened.append({"group": group, "dataset_id": dataset_id, "time_start": str(t0), "time_end": str(t1)})
            day = todo["sample_time"].dt.floor("D")
            covered = todo[(day - pd.Timedelta(days=days) >= t0) & (day <= t1)]
            print(f"  {dataset_id}: {t0:%Y-%m-%d} to {t1:%Y-%m-%d}; {len(covered):,} samples covered", flush=True)
            started, count = clock.time(), [0]

            def work(row):
                values = compute(group, ds, variables, row, days, box)
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
    write_outputs(samples, partial, opened, args)


def compute(group, ds, variables, row, days, box) -> dict:
    """The group's features for one sample; NaN when the data are missing."""
    t = row["sample_time"].floor("D")
    for attempt in range(1, 6):
        try:
            fields = {v: subset(ds[v], row["sample_lat"], row["sample_lon"], t - pd.Timedelta(days=days),
                                t + pd.Timedelta(hours=23, minutes=59), box) for v in variables}  # fmt: skip
            break
        except Exception as exc:  # network hiccup
            if attempt == 5:
                print(f"    {row['sample_id']}: failed ({str(exc)[:100]})", flush=True)
                return {f: np.nan for f in FEATURES[group]}
            clock.sleep(5 * attempt)
    first = fields[variables[0]]
    n = int(first.sizes.get("time", 1))
    if group == "mld":
        series = first.mean(dim=["latitude", "longitude"], skipna=True).values
        return {"mld_max30": np.nanmax(series) if np.isfinite(series).any() else np.nan,
                "mld_now": series[-1] if len(series) else np.nan, "mld_n": n}  # fmt: skip
    if group == "sst":
        series = first.mean(dim=["latitude", "longitude"], skipna=True).values
        before = series[:-1]
        change = series[-1] - np.nanmean(before) if len(series) > 1 and np.isfinite(before).any() else np.nan
        return {"sst_change14": change, "sst_n": n}
    if group == "wind":
        speed = np.hypot(fields["eastward_wind"], fields["northward_wind"]).mean(dim=["latitude", "longitude"], skipna=True).values
        ok = np.isfinite(speed).any()
        return {"wind_max14": np.nanmax(speed) if ok else np.nan, "wind_mean14": np.nanmean(speed) if ok else np.nan, "wind_n": n}
    if group == "front":
        zos = first.isel(time=-1) if "time" in first.dims else first
        lat, lon = zos["latitude"].values, zos["longitude"].values
        if len(lat) < 2 or len(lon) < 2:
            return {"ssh_gradient": np.nan, "front_n": n}
        dy = np.gradient(lat) * 111.2  # km per step
        dx = np.gradient(np.unwrap(np.radians(lon))) * 6371.0 * np.cos(np.radians(row["sample_lat"]))
        gy, gx = np.gradient(zos.values, axis=0) / dy[:, None], np.gradient(zos.values, axis=1) / dx[None, :]
        return {"ssh_gradient": np.nanmean(np.hypot(gx, gy)) * 100.0, "front_n": n}  # m per 100 km
    if group == "eke":
        eke = 0.5 * (fields["ugosa"] ** 2 + fields["vgosa"] ** 2)
        return {"eke30": float(eke.mean(skipna=True)), "eke_n": n}
    raise ValueError(group)


def subset(da, lat, lon, start, end, box):
    """The variable within +-box degrees of (lat, lon) between start and end, loaded; date line handled."""
    if "depth" in da.dims:
        da = da.isel(depth=0)
    lats, lons = da["latitude"].values, da["longitude"].values
    lat_slice = slice(lat - box, lat + box) if lats[0] < lats[-1] else slice(lat + box, lat - box)
    east = lons.max() > 180  # grid on 0-360
    if east:
        lon = lon % 360
    parts = []
    for lo, hi in lon_ranges(lon, box, east):
        parts.append(da.sel(time=slice(start, end), latitude=lat_slice, longitude=slice(lo, hi)))
    out = parts[0] if len(parts) == 1 else __import__("xarray").concat(parts, dim="longitude")
    if out.sizes["latitude"] == 0 or out.sizes["longitude"] == 0:  # box smaller than one cell: take the nearest
        out = da.sel(time=slice(start, end)).sel(latitude=[lat], longitude=[lon], method="nearest")
    return out.load(scheduler="synchronous")


def lon_ranges(lon, box, east=False):
    low, high = (0, 360) if east else (-180, 180)
    lo, hi = lon - box, lon + box
    if lo < low:
        return [(lo + 360, high), (low, hi)]
    if hi > high:
        return [(lo, high), (low, hi - 360)]
    return [(lo, hi)]


def open_dataset(dataset_id: str):
    import copernicusmarine

    # time-series chunks: a long window at one place is a few chunks (the history features read 15-31 days per sample)
    try:
        return copernicusmarine.open_dataset(dataset_id=dataset_id, service="arco-time-series")
    except TypeError:
        return copernicusmarine.open_dataset(dataset_id=dataset_id)


def describe(dataset_id: str, variables: list[str]) -> str:
    try:
        ds = open_dataset(dataset_id)
    except Exception as exc:
        return f"{dataset_id}: NOT OPENED ({str(exc).splitlines()[0][:120]})"
    t = pd.to_datetime(ds["time"].values[[0, -1]])
    missing = [v for v in variables if v not in ds.data_vars]
    return (f"{dataset_id}: {t[0]:%Y-%m-%d} to {t[1]:%Y-%m-%d}; variables {'all present' if not missing else 'MISSING ' + ', '.join(missing)}"
            f"; has {', '.join(list(ds.data_vars)[:8])}")


def write_outputs(samples, partial, opened, args) -> None:
    wide = samples[["sample_id", "section", "sample_time", "sample_lat", "sample_lon"]].set_index("sample_id")
    for group, features in FEATURES.items():
        path = partial / f"{group}.csv"
        if not path.exists():
            continue
        part = pd.read_csv(path, on_bad_lines="skip").drop_duplicates("sample_id", keep="last").set_index("sample_id")
        for f in features:
            wide[f] = part[f].reindex(wide.index) if f in part else np.nan
        wide[f"{group}_dataset"] = part["dataset_id"].reindex(wide.index)
        n_col = f"{group}_n"
        if n_col in part:
            wide[n_col] = part[n_col].reindex(wide.index)
    wide.reset_index().to_csv(args.out_dir / "sample_history_fields.csv", index=False)
    (args.out_dir / "provenance_history.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "box_deg": BOX, "front_box_deg": FRONT_BOX,
        "windows_days": {g: GROUPS[g][2] for g in GROUPS}, "datasets": opened, "n_samples": len(wide)}, indent=2))  # fmt: skip
    print(f"\nWrote {args.out_dir / 'sample_history_fields.csv'}")
    summary = {f: f"{int(wide[f].notna().sum()):,} values, median {wide[f].median():.3g}" for g in FEATURES for f in FEATURES[g] if f in wide}
    for f, text in summary.items():
        print(f"  {f:14s} {text}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifests", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "manifests" / "biogoship_tile_manifest.csv",
                                 REPO_ROOT / "data" / "faire" / "manifests" / "biogoship_tile_manifest.csv"])  # fmt: skip
    parser.add_argument("--groups", nargs="+", choices=list(GROUPS), default=list(GROUPS))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "model_fields")
    parser.add_argument("--list-datasets", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    main()
