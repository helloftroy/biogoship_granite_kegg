#!/usr/bin/env python3
"""Ocean-model values (GLORYS12 physics, PISCES biogeochemistry) at each sample's own position, time and depth.

For every sample in the tile manifests (both sample sets) this reads the Copernicus
Marine reanalyses at the **sampling** latitude, longitude, UTC date and depth, not
at the backtracked or satellite position:

  glorys   GLOBAL_MULTIYEAR_PHY_001_030 (GLORYS12V1, 1/12 deg, daily means):
           thetao, so, mlotst, zos, uo, vo
  pisces   GLOBAL_MULTIYEAR_BGC_001_029 (PISCES hindcast, 1/4 deg, daily means):
           no3, po4, si, fe, o2, chl, phyc, nppv, ph, dissic, talk, spco2

Each product has datasets in priority order: the multi-year reanalysis ("my"), then its
interim continuation ("myint"), then (only with --allow-forecast) the analysis/forecast
system for dates the reanalysis does not reach yet. The dataset that supplied each
value is recorded, so values from different datasets can be told apart or left out.

Matching, per sample and variable:
  time    the model day that contains the sample's UTC timestamp (daily means;
          the label offset of the dataset, 00:00 or 12:00, is taken into account).
          Samples outside a dataset's time range go to the next dataset.
  depth   the model level closest to the sample's depth (depth_m; when that is
          missing, the shallowest level, flagged depth_assumed = 1). 2-D variables
          (mlotst, zos, spco2) have no depth.
  space   bilinear interpolation between the 4 surrounding cells when all 4 have a
          value ("bilinear"); otherwise the closest cell with a value within
          --window cells ("nearest_valid"); otherwise missing ("no_valid_cell").
          distance_km is from the sample to the closest of the cells used.
Missing stays missing. In-situ measurements are never used to fill anything.

Outputs (data/model_fields/):
  sample_model_fields_long.csv   one row per sample x variable: value, unit, product,
                                 dataset_id, model_time, model_depth_m, method,
                                 distance_km, model_lat, model_lon, depth_assumed
  sample_model_fields.csv        one row per sample: sample columns, then per
                                 variable <product>_<var> plus <product>_<var>_dataset
                                 and <product>_<var>_distance_km
  provenance.json                toolbox version, datasets opened with their time
                                 range and attributes, settings, extraction time

Needs a Copernicus Marine login (COPERNICUSMARINE_SERVICE_USERNAME / _PASSWORD, as for
the currents and OSTIA jobs) and internet: run cluster/model_fields.sbatch.

Usage:
    python scripts/extract_model_fields.py
    python scripts/extract_model_fields.py --products pisces --allow-forecast
    python scripts/extract_model_fields.py --list-datasets   # print what each dataset offers
"""

from __future__ import annotations

import argparse
import json
import sys
import time as clock
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

PHY_VARS = ["thetao", "so", "mlotst", "zos", "uo", "vo"]
BGC_VARS = ["no3", "po4", "si", "fe", "o2", "chl", "phyc", "nppv", "ph", "dissic", "talk", "spco2"]
# product -> datasets in priority order: (dataset_id, variables it may hold, is_forecast)
PRODUCTS = {
    "glorys": [
        ("cmems_mod_glo_phy_my_0.083deg_P1D-m", PHY_VARS, False),
        ("cmems_mod_glo_phy_myint_0.083deg_P1D-m", PHY_VARS, False),
        ("cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m", ["thetao"], True),
        ("cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m", ["so"], True),
        ("cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m", ["uo", "vo"], True),
        ("cmems_mod_glo_phy_anfc_0.083deg_P1D-m", ["mlotst", "zos"], True),
    ],
    "pisces": [
        ("cmems_mod_glo_bgc_my_0.25deg_P1D-m", BGC_VARS, False),
        ("cmems_mod_glo_bgc_myint_0.25deg_P1D-m", BGC_VARS, False),
        ("cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m", ["no3", "po4", "si", "fe"], True),
        ("cmems_mod_glo_bgc-bio_anfc_0.25deg_P1D-m", ["o2", "nppv"], True),
        ("cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m", ["chl", "phyc"], True),
        ("cmems_mod_glo_bgc-car_anfc_0.25deg_P1D-m", ["ph", "dissic", "talk"], True),
        ("cmems_mod_glo_bgc-co2_anfc_0.25deg_P1D-m", ["spco2"], True),
    ],
}
SAMPLE_COLUMNS = ["sample_id", "section", "sample_time", "sample_lat", "sample_lon", "depth_m"]


def main() -> None:
    args = parse_args()
    if args.list_datasets:
        list_datasets(args)
        return
    samples = load_samples(args.manifests)
    print(f"{len(samples):,} samples from {len(args.manifests)} manifests, "
          f"{samples['sample_time'].min():%Y-%m-%d} to {samples['sample_time'].max():%Y-%m-%d}", flush=True)  # fmt: skip
    args.out_dir.mkdir(parents=True, exist_ok=True)
    long_rows, opened = [], []
    for product in args.products:
        done: dict[str, set] = {}  # variable -> sample_ids already filled
        for dataset_id, variables, is_forecast in PRODUCTS[product]:
            if is_forecast and not args.allow_forecast:
                continue
            try:
                ds = open_dataset(dataset_id, variables, samples, args)
            except Exception as exc:  # unknown dataset id, no access, ...
                print(f"{product:7s} {dataset_id}: not opened ({str(exc).splitlines()[0][:160]})", flush=True)
                opened.append({"product": product, "dataset_id": dataset_id, "error": str(exc)[:500]})
                continue
            present = [v for v in variables if v in ds.data_vars]
            info = dataset_info(ds, product, dataset_id, present)
            opened.append(info)
            print(f"{product:7s} {dataset_id}: {info['time_start']} to {info['time_end']}, {', '.join(present)}", flush=True)
            for var in present:
                todo = samples[~samples["sample_id"].isin(done.get(var, set()))]
                if todo.empty:
                    continue
                started = clock.time()
                rows = extract(ds, var, todo, window=args.window)
                rows = rows.assign(product=product, variable=var, dataset_id=dataset_id, forecast=int(is_forecast),
                                   unit=ds[var].attrs.get("units", ""))  # fmt: skip
                found = rows[rows["method"] != "outside_time_range"]
                long_rows.append(found)
                done.setdefault(var, set()).update(found["sample_id"])
                n_val = int(found["value"].notna().sum())
                print(f"        {var:7s} {len(found):5,} samples in range, {n_val:5,} with a value ({clock.time() - started:5.0f} s)",
                      flush=True)  # fmt: skip
            ds.close()
    long = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame()
    write_outputs(samples, long, opened, args)


def load_samples(manifests: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(m, usecols=lambda c: c in SAMPLE_COLUMNS, low_memory=False) for m in manifests if m.exists()]
    samples = pd.concat(frames).drop_duplicates("sample_id").reset_index(drop=True)
    samples["sample_time"] = pd.to_datetime(samples["sample_time"], utc=True).dt.tz_localize(None)
    samples = samples.dropna(subset=["sample_time", "sample_lat", "sample_lon"])
    samples["sample_lon"] = ((samples["sample_lon"] + 180) % 360) - 180
    return samples


def open_dataset(dataset_id: str, variables: list[str], samples: pd.DataFrame, args: argparse.Namespace):
    import copernicusmarine

    depth_max = float(np.nanmax(samples["depth_m"])) if samples["depth_m"].notna().any() else 1.0
    return copernicusmarine.open_dataset(
        dataset_id=dataset_id, variables=variables,
        minimum_depth=0.0, maximum_depth=max(depth_max, 1.0) * 1.5 + 10.0,
        start_datetime=(samples["sample_time"].min() - pd.Timedelta(days=2)).strftime("%Y-%m-%dT%H:%M:%S"),
        end_datetime=(samples["sample_time"].max() + pd.Timedelta(days=2)).strftime("%Y-%m-%dT%H:%M:%S"),
    )  # fmt: skip


def dataset_info(ds, product: str, dataset_id: str, variables: list[str]) -> dict:
    times = pd.to_datetime(ds["time"].values)
    keep = ("title", "source", "institution", "references", "product_id", "product_version", "contact", "Conventions")
    return {"product": product, "dataset_id": dataset_id, "variables": variables,
            "time_start": f"{times.min():%Y-%m-%d %H:%M}", "time_end": f"{times.max():%Y-%m-%d %H:%M}",
            "depths_m": [round(float(d), 3) for d in ds["depth"].values] if "depth" in ds.dims else [],
            "attributes": {k: str(v) for k, v in ds.attrs.items() if k in keep},
            "units": {v: ds[v].attrs.get("units", "") for v in variables}}  # fmt: skip


def extract(ds, var: str, samples: pd.DataFrame, window: int = 4) -> pd.DataFrame:
    """Value of ds[var] at each sample's day, nearest depth level and position (see module docstring)."""
    da = ds[var]
    lat = ds["latitude"].values.astype(float)
    lon = ds["longitude"].values.astype(float)
    times = pd.to_datetime(ds["time"].values)
    offset = times[0] - times[0].floor("D")  # daily means labelled at 00:00 or 12:00
    target = samples["sample_time"].dt.floor("D") + offset
    t_index = pd.Index(times).get_indexer(target)  # exact day only; -1 when outside the range
    has_depth = "depth" in da.dims
    depths = ds["depth"].values.astype(float) if has_depth else None

    out = samples[["sample_id"]].copy().reset_index(drop=True)
    for col, default in (("value", np.nan), ("model_time", pd.NaT), ("model_depth_m", np.nan), ("method", "outside_time_range"),
                         ("distance_km", np.nan), ("model_lat", np.nan), ("model_lon", np.nan), ("depth_assumed", 0)):  # fmt: skip
        out[col] = default
    s_lat, s_lon = samples["sample_lat"].to_numpy(float), samples["sample_lon"].to_numpy(float)
    s_depth = samples["depth_m"].to_numpy(float)
    if has_depth:
        assumed = np.isnan(s_depth)
        d_index = np.where(assumed, 0, np.abs(depths[None, :] - np.nan_to_num(s_depth)[:, None]).argmin(axis=1))
        out["depth_assumed"] = assumed.astype(int)
    lat_ascending = lat[0] < lat[-1]
    for ti in np.unique(t_index[t_index >= 0]):
        members = np.flatnonzero(t_index == ti)
        # one box of cells around all of this day's samples
        lat_i = nearest_index(lat, s_lat[members])
        lon_i = nearest_index(lon, s_lon[members], periodic=True)
        rows = np.arange(max(lat_i.min() - window - 1, 0), min(lat_i.max() + window + 2, len(lat)))
        cols = lon_box(lon_i, window, len(lon))
        sel = {"time": int(ti), "latitude": rows, "longitude": cols}
        if has_depth:
            levels = np.unique(d_index[members])
            sel["depth"] = levels
        block = np.asarray(da.isel(sel).values, dtype=float)  # (depth?, lat, lon)
        if not has_depth:
            block = block[None]
            levels = np.array([0])
        for m in members:
            level = int(np.searchsorted(levels, d_index[m])) if has_depth else 0
            field = block[level]
            result = point_value(field, lat[rows], lon[cols], s_lat[m], s_lon[m], window, lat_ascending)
            out.loc[m, ["value", "method", "distance_km", "model_lat", "model_lon"]] = result
            out.loc[m, "model_time"] = times[ti]
            if has_depth:
                out.loc[m, "model_depth_m"] = depths[d_index[m]]
    return out


def nearest_index(grid: np.ndarray, values: np.ndarray, periodic: bool = False) -> np.ndarray:
    diff = values[:, None] - grid[None, :]
    if periodic:
        diff = (diff + 180) % 360 - 180
    return np.abs(diff).argmin(axis=1)


def lon_box(lon_i: np.ndarray, window: int, n: int) -> np.ndarray:
    """Longitude indices around the samples, wrapping at the date line, in grid order of the circle."""
    pad = window + 1
    needed = np.unique(np.concatenate([(i + np.arange(-pad, pad + 1)) % n for i in lon_i]))
    if needed.size == n:
        return needed
    # start after the largest gap so a box across the date line stays contiguous
    gaps = np.diff(np.concatenate([needed, [needed[0] + n]]))
    start = (np.argmax(gaps) + 1) % needed.size
    return np.roll(needed, -start)


def point_value(field, lat, lon, p_lat, p_lon, window, lat_ascending) -> list:
    """Bilinear if the 4 surrounding cells have values, else the nearest valid cell within `window` cells."""
    dlon = (lon - p_lon + 180) % 360 - 180  # signed longitude offset of every column
    # the 4 surrounding cells
    order = np.argsort(lat) if not lat_ascending else np.arange(len(lat))
    lat_sorted = lat[order]
    k = np.searchsorted(lat_sorted, p_lat)
    col_order = np.argsort(dlon)
    dlon_sorted = dlon[col_order]
    j = np.searchsorted(dlon_sorted, 0.0)
    if 0 < k < len(lat_sorted) and 0 < j < len(dlon_sorted):
        r0, r1 = order[k - 1], order[k]
        c0, c1 = col_order[j - 1], col_order[j]
        quad = field[np.ix_([r0, r1], [c0, c1])]
        if np.isfinite(quad).all() and lat_sorted[k] > lat_sorted[k - 1] and dlon_sorted[j] > dlon_sorted[j - 1]:
            wy = (p_lat - lat_sorted[k - 1]) / (lat_sorted[k] - lat_sorted[k - 1])
            wx = (0.0 - dlon_sorted[j - 1]) / (dlon_sorted[j] - dlon_sorted[j - 1])
            value = (quad[0, 0] * (1 - wy) * (1 - wx) + quad[0, 1] * (1 - wy) * wx
                     + quad[1, 0] * wy * (1 - wx) + quad[1, 1] * wy * wx)  # fmt: skip
            corners = [(lat[r], lon[c]) for r in (r0, r1) for c in (c0, c1)]
            dist = [haversine_km(p_lat, p_lon, a, b) for a, b in corners]
            best = int(np.argmin(dist))
            return [float(value), "bilinear", float(dist[best]), corners[best][0], corners[best][1]]
    # nearest valid cell within the window
    ri = int(np.abs(lat - p_lat).argmin())
    ci = int(np.abs(dlon).argmin())
    r_lo, r_hi = max(ri - window, 0), min(ri + window + 1, len(lat))
    c_sel = np.flatnonzero(np.abs(np.arange(len(lon)) - ci) <= window)
    sub = field[r_lo:r_hi][:, c_sel]
    valid = np.argwhere(np.isfinite(sub))
    if valid.size == 0:
        return [np.nan, "no_valid_cell", np.nan, np.nan, np.nan]
    cand_lat = lat[r_lo:r_hi][valid[:, 0]]
    cand_lon = lon[c_sel][valid[:, 1]]
    dist = haversine_km(p_lat, p_lon, cand_lat, cand_lon)
    best = int(np.argmin(dist))
    return [float(sub[valid[best, 0], valid[best, 1]]), "nearest_valid", float(dist[best]), float(cand_lat[best]), float(cand_lon[best])]


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, np.asarray(lat2, float), np.asarray(lon2, float)))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def write_outputs(samples: pd.DataFrame, long: pd.DataFrame, opened: list, args: argparse.Namespace) -> None:
    if not long.empty:
        long = long.merge(samples[["sample_id", "section"]], on="sample_id", how="left")
    long.to_csv(args.out_dir / "sample_model_fields_long.csv", index=False)
    wide = samples.copy()
    if not long.empty:
        for (product, var), part in long.groupby(["product", "variable"]):
            part = part.set_index("sample_id")
            name = f"{product}_{var}"
            wide[name] = wide["sample_id"].map(part["value"])
            wide[f"{name}_dataset"] = wide["sample_id"].map(part["dataset_id"])
            wide[f"{name}_distance_km"] = wide["sample_id"].map(part["distance_km"]).round(2)
    wide.to_csv(args.out_dir / "sample_model_fields.csv", index=False)
    import copernicusmarine

    provenance = {"extracted_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                  "copernicusmarine_version": getattr(copernicusmarine, "__version__", "unknown"),
                  "manifests": [str(m) for m in args.manifests], "window_cells": args.window,
                  "allow_forecast": args.allow_forecast, "n_samples": len(samples), "datasets": opened}  # fmt: skip
    (args.out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2))
    print(f"Wrote {args.out_dir / 'sample_model_fields.csv'} ({wide.shape[1]} columns), the long table and provenance.json")
    if not long.empty:
        summary = long.groupby(["product", "variable", "dataset_id"]).agg(
            samples=("sample_id", "size"), with_value=("value", lambda v: int(v.notna().sum())),
            bilinear=("method", lambda m: int((m == "bilinear").sum())),
            median_km=("distance_km", "median"))  # fmt: skip
        print(summary.round(1).to_string())


def list_datasets(args: argparse.Namespace) -> None:
    import copernicusmarine

    for product in args.products:
        for dataset_id, _, _ in PRODUCTS[product]:
            try:
                ds = copernicusmarine.open_dataset(dataset_id=dataset_id)
                times = pd.to_datetime(ds["time"].values)
                print(f"{product:7s} {dataset_id}: {times.min():%Y-%m-%d} to {times.max():%Y-%m-%d}; {', '.join(ds.data_vars)}")
                ds.close()
            except Exception as exc:
                print(f"{product:7s} {dataset_id}: {str(exc).splitlines()[0][:160]}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifests", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "manifests" / "biogoship_tile_manifest.csv",
                                 REPO_ROOT / "data" / "faire" / "manifests" / "biogoship_tile_manifest.csv"])  # fmt: skip
    parser.add_argument("--products", nargs="+", choices=list(PRODUCTS), default=list(PRODUCTS))
    parser.add_argument("--allow-forecast", action="store_true",
                        help="Use the analysis/forecast datasets for dates the reanalyses do not cover (recorded per value).")
    parser.add_argument("--window", type=int, default=4, help="Cells searched in each direction for a valid cell (default 4).")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "model_fields")
    parser.add_argument("--list-datasets", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    main()
