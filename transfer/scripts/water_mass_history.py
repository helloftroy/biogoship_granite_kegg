#!/usr/bin/env python3
"""Water-mass origin and age, and light history, along 60-day back-trajectories (RESULTS.md R42).

Axis 2 (phototrophs) is smooth over ~560 km (R34), about a month of surface drift. If the community
carries the imprint of where its water came from, the Lagrangian history should explain it. Per sample:

  step currents  Copernicus Marine GlobCurrent (MULTIOBS_GLO_PHY_MYNRT_015_003, 0.25 deg, daily, 0 and 15 m)
                 per section, the samples' box +- MARGIN deg, from DAYS + 2 days before the first sample to
                 1 day after the last (scripts/download_currents.py fetch / join_across_dateline;
                 multi-year dataset, near-real-time after its end). transfer/data/currents60/.
  step track     one particle per sample, released at the sampling position and time, advected BACKWARD
                 DAYS days with the 15 m currents (OceanParcels RK4, dt 30 min, no diffusion; as
                 biogoship/advection.py, backward only), positions every 6 h. A particle that leaves the
                 currents box is dropped from then on; one stuck on land stays put (velocity 0) and is
                 flagged. transfer/results/water_mass/trajectories.csv.
  step fields    at the daily positions (day 0, -1, ..., -DAYS; interpolated linearly between the 6-hourly
                 track points, see daily_positions), nearest grid value of:
                   sst      OSTIA L4 analysed_sst (reprocessed, then near-real-time), C
                   chl      GlobColour L4 gap-free CHL (multi-year, interim, near-real-time), mg m-3
                   mld      GLORYS12 mlotst (multi-year, interim, analysis-forecast), m
                   kd490    GlobColour L4 gap-free transparency KD490 if the dataset opens, m-1
                   depth    GLORYS12 static deptho (sea-floor depth), m
                 transfer/results/water_mass/along_track_daily.csv (long: sample_id, day, lat, lon, variable, value).
                 (along_track.csv from the first run is obsolete: whole-day points were selected, not
                 interpolated, which left only day 0 for almost every track.)
                 Resumable (values appended per sample and variable).
  --probe-light  lists candidate Copernicus Marine ocean-colour datasets with their variables, to find
                 a daily PAR product (for the light-history test) and KD490.

Features (time since shelf water, trajectory-mean SST / chlorophyll, origin region, mixed-layer light)
are computed afterwards on the Mac from along_track_daily.csv and trajectories.csv
(transfer/scripts/water_mass_test.py), so their definitions can change without re-running this.

Usage (Morrill, service partition: transfer/cluster/water_mass.sbatch):
    python transfer/scripts/water_mass_history.py --probe-light
    python transfer/scripts/water_mass_history.py --steps currents track fields --workers 4
"""

from __future__ import annotations

import argparse
import sys
import threading
import time as clock
from concurrent.futures import ThreadPoolExecutor
from datetime import timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from biogoship import advection  # noqa: E402
from download_currents import covers, fetch, join_across_dateline  # noqa: E402
from extract_history_fields import open_dataset  # noqa: E402
from extract_model_fields import load_samples  # noqa: E402

DAYS, MARGIN = 60, 15.0
CURRENTS_DIR = REPO / "transfer" / "data" / "currents60"
OUT = REPO / "transfer" / "results" / "water_mass"
FIELDS = {
    "sst": (["METOFFICE-GLO-SST-L4-REP-OBS-SST", "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2"], "analysed_sst"),
    "chl": (["cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D",
             "cmems_obs-oc_glo_bgc-plankton_myint_l4-gapfree-multi-4km_P1D",
             "cmems_obs-oc_glo_bgc-plankton_nrt_l4-gapfree-multi-4km_P1D"], "CHL"),
    "mld": (["cmems_mod_glo_phy_my_0.083deg_P1D-m", "cmems_mod_glo_phy_myint_0.083deg_P1D-m",
             "cmems_mod_glo_phy_anfc_0.083deg_P1D-m"], "mlotst"),
    "kd490": (["cmems_obs-oc_glo_bgc-transp_my_l4-gapfree-multi-4km_P1D",
               "cmems_obs-oc_glo_bgc-transp_myint_l4-gapfree-multi-4km_P1D",
               "cmems_obs-oc_glo_bgc-transp_nrt_l4-gapfree-multi-4km_P1D"], "KD490"),
}  # fmt: skip
STATIC = (["cmems_mod_glo_phy_my_0.083deg_static", "cmems_mod_glo_phy_anfc_0.083deg_static"], "deptho")
LIGHT_CANDIDATES = [
    "cmems_obs-oc_glo_bgc-transp_my_l4-gapfree-multi-4km_P1D", "cmems_obs-oc_glo_bgc-transp_my_l3-multi-4km_P1D",
    "cmems_obs-oc_glo_bgc-optics_my_l3-multi-4km_P1D", "cmems_obs-oc_glo_bgc-optics_my_l4-multi-4km_P1M",
    "cmems_obs-oc_glo_bgc-pp_my_l4-multi-4km_P1M", "cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D",
    "cmems_obs-oc_glo_bgc-reflectance_my_l3-multi-4km_P1D",
]  # fmt: skip
LOCK = threading.Lock()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifests", type=Path, nargs="+",
                   default=[REPO / "data" / "manifests" / "biogoship_tile_manifest.csv",
                            REPO / "data" / "faire" / "manifests" / "biogoship_tile_manifest.csv"])  # fmt: skip
    p.add_argument("--steps", nargs="+", default=["currents", "track", "fields"], choices=["currents", "track", "fields"])
    p.add_argument("--fields", nargs="+", default=list(FIELDS) + ["depth"])
    p.add_argument("--days", type=int, default=DAYS)
    p.add_argument("--margin", type=float, default=MARGIN)
    p.add_argument("--depth-m", type=float, default=15.0)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--probe-light", action="store_true")
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    CURRENTS_DIR.mkdir(parents=True, exist_ok=True)
    if args.probe_light:
        probe()
        return
    samples = load_samples(args.manifests)
    samples["sample_time"] = samples["sample_time"].dt.tz_localize("UTC")
    print(f"{len(samples):,} samples in {samples['section'].nunique()} sections", flush=True)
    if "currents" in args.steps:
        step_currents(samples, args)
    if "track" in args.steps:
        step_track(samples, args)
    if "fields" in args.steps:
        step_fields(samples, args)


def boxes(samples: pd.DataFrame, args) -> list[dict]:
    out = []
    for section, g in samples.groupby("section"):
        lats, lons = g["sample_lat"].tolist(), g["sample_lon"].tolist()
        wraps = max(lons) - min(lons) > 180.0 and advection.narrower_in_360(lons)
        if wraps:
            lons = [x % 360.0 for x in lons]
        lo, hi = min(lons) - args.margin, max(lons) + args.margin
        if not wraps and (lo < -180.0 or hi > 180.0):  # the margin crosses the date line (e.g. HLY2202, Bering Strait)
            lo, hi, wraps = lo % 360.0, hi % 360.0 if hi > 180.0 else hi + 360.0, True
        out.append({"section": section, "sites": g, "lat_min": max(min(lats) - args.margin, -89.5),
                    "lat_max": min(max(lats) + args.margin, 89.5), "lon_min": lo,
                    "lon_max": hi, "wraps": wraps,
                    "start": (g["sample_time"].min() - timedelta(days=args.days + 2)).to_pydatetime(),
                    "end": (g["sample_time"].max() + timedelta(days=1)).to_pydatetime()})  # fmt: skip
    return out


def dataset_for(box) -> str:
    return advection.CURRENTS_DATASET if box["end"] <= advection.CURRENTS_MY_END else advection.CURRENTS_DATASET_NRT


def step_currents(samples, args) -> None:
    for box in boxes(samples, args):
        dest = CURRENTS_DIR / advection.section_filename(box["section"])
        if dest.exists() and covers(dest, box):
            print(f"{box['section']:14s} currents on disk", flush=True)
            continue
        ds = dataset_for(box)
        print(f"{box['section']:14s} {len(box['sites']):4d} samples | lon {box['lon_min']:.1f}..{box['lon_max']:.1f} "
              f"lat {box['lat_min']:.1f}..{box['lat_max']:.1f} | {box['start']:%Y-%m-%d}..{box['end']:%Y-%m-%d} | {ds}"
              + (" | across 180" if box["wraps"] else ""), flush=True)  # fmt: skip
        if not box["wraps"]:
            fetch(ds, box, box["lon_min"], box["lon_max"], CURRENTS_DIR, dest.name)
        else:
            west, east = dest.with_suffix(".west.nc"), dest.with_suffix(".east.nc")
            fetch(ds, box, box["lon_min"], 180.0, CURRENTS_DIR, west.name)
            fetch(ds, box, -180.0, box["lon_max"] - 360.0, CURRENTS_DIR, east.name)
            join_across_dateline(west, east, dest)


def backward(currents, releases: list[dict], max_hours: float, dt_minutes=30.0, output_minutes=360.0) -> dict:
    """advection.advect, backward only: id -> rows (hours <= 0, lat, lon)."""
    import parcels

    fieldset = parcels.FieldSet.from_xarray_dataset(currents, variables={"U": "uo", "V": "vo"},
                                                    dimensions={"lon": "longitude", "lat": "latitude", "time": "time"},
                                                    mesh="spherical")  # fmt: skip
    origin = np.datetime64(currents["time"].values[0], "ns")
    release_s = np.array([(np.datetime64(r["time"].astimezone(timezone.utc).replace(tzinfo=None), "ns") - origin)
                          / np.timedelta64(1, "s") for r in releases])  # fmt: skip

    def delete_lost(particle, fieldset, time):  # noqa: ARG001
        if particle.state >= 50:
            particle.delete()

    step = timedelta(minutes=output_minutes)
    tracks = {i: [] for i in range(len(releases))}
    pset = parcels.ParticleSet(fieldset=fieldset, pclass=parcels.ScipyParticle, lon=[r["lon"] for r in releases],
                               lat=[r["lat"] for r in releases], time=release_s)  # fmt: skip
    index_of = {pt.id: i for i, pt in enumerate(pset)}
    span = (release_s.max() - release_s.min()) / 3600.0 + max_hours
    for _ in range(int(np.ceil(span * 60.0 / output_minutes)) + 1):
        if len(pset) == 0:
            break
        pset.execute([parcels.AdvectionRK4, delete_lost], runtime=step, dt=-timedelta(minutes=dt_minutes), verbose_progress=False)
        done = []
        for pt in pset:
            i = index_of[pt.id]
            h = (float(pt.time) - release_s[i]) / 3600.0
            if h <= 0:
                tracks[i].append((h, float(pt.lat), float(pt.lon)))
            done.append(h <= -max_hours)
        if any(done):
            pset.remove_booleanvector(np.array(done))
    out = {}
    for i, r in enumerate(releases):
        rows = np.array([(0.0, r["lat"], r["lon"]), *tracks[i]], dtype=float)
        rows = rows[np.unique(np.round(rows[:, 0], 6), return_index=True)[1]]
        out[r["id"]] = rows
    return out


def step_track(samples, args) -> None:
    import logging
    import warnings

    warnings.filterwarnings("ignore", module="parcels")
    logging.getLogger("parcels").setLevel(logging.ERROR)
    path = OUT / "trajectories.csv"
    done = set(pd.read_csv(path)["sample_id"]) if path.exists() else set()
    for box in boxes(samples, args):
        g = box["sites"][~box["sites"]["sample_id"].isin(done)]
        if g.empty:
            continue
        cur_path = CURRENTS_DIR / advection.section_filename(box["section"])
        if not cur_path.exists():
            print(f"{box['section']}: no currents file", flush=True)
            continue
        t0 = clock.time()
        cur = advection.load_currents(cur_path, args.depth_m)
        in360 = float(cur["longitude"].max()) > 180.0
        rel = [{"id": r.sample_id, "lat": r.sample_lat, "lon": r.sample_lon % 360.0 if in360 else r.sample_lon,
                "time": r.sample_time.to_pydatetime()} for r in g.itertuples()]  # fmt: skip
        tr = backward(cur, rel, args.days * 24.0)
        rows = []
        for sid, a in tr.items():
            stuck = len(a) > 4 and np.allclose(a[-4:, 1:], a[-1, 1:], atol=1e-6)
            for h, la, lo in a:
                rows.append({"sample_id": sid, "hours": round(h, 3), "lat": round(la, 5),
                             "lon": round(((lo + 180) % 360) - 180, 5), "stuck_at_end": stuck})  # fmt: skip
        pd.DataFrame(rows).to_csv(path, mode="a", header=not path.exists(), index=False)
        reach = [a[:, 0].min() / -24 for a in tr.values()]
        print(f"{box['section']:14s} {len(rel):4d} tracks in {clock.time() - t0:5.0f} s; reach median {np.median(reach):.1f} d, "
              f"full {DAYS} d: {np.mean(np.array(reach) >= args.days - 0.3):.0%}", flush=True)  # fmt: skip


def daily_positions(tr: pd.DataFrame) -> pd.DataFrame:
    """Positions at exactly 0, -24, -48, ... hours, linearly interpolated between the 6-hourly outputs.

    Parcels writes on a 6-hour clock that is not aligned to each sample's own time (offsets of a fraction
    of an hour), so whole-day outputs exist only at day 0; selecting them (the first version) left one
    point per track. Longitude is unwrapped before interpolating so a track across 180 deg stays continuous."""
    rows = []
    for sid, g in tr.groupby("sample_id", sort=False):
        g = g.sort_values("hours")
        h = g["hours"].to_numpy()
        days = np.arange(0, int(np.floor(h.min() / 24)) - 1, -1)
        days = days[24 * days >= h.min() - 1e-6]
        lon = np.degrees(np.unwrap(np.radians(g["lon"].to_numpy())))
        lo = np.interp(24.0 * days, h, lon)
        rows.append(pd.DataFrame({"sample_id": sid, "day": days, "lat": np.interp(24.0 * days, h, g["lat"].to_numpy()),
                                  "lon": (lo + 180) % 360 - 180}))  # fmt: skip
    return pd.concat(rows, ignore_index=True)


def step_fields(samples, args) -> None:
    import xarray as xr

    tr = daily_positions(pd.read_csv(OUT / "trajectories.csv"))
    st = samples.set_index("sample_id")["sample_time"]
    tr["time"] = tr["sample_id"].map(st).dt.tz_convert(None).dt.floor("D") + pd.to_timedelta(tr["day"], unit="D")
    path = OUT / "along_track_daily.csv"
    done = set(map(tuple, pd.read_csv(path, usecols=["sample_id", "variable"]).drop_duplicates().values)) if path.exists() else set()
    groups = {sid: g for sid, g in tr.groupby("sample_id")}
    print(f"{len(groups):,} tracks with daily points", flush=True)

    def write(sid, var, g, values):
        df = pd.DataFrame({"sample_id": sid, "day": g["day"].to_numpy(), "lat": g["lat"].to_numpy(), "lon": g["lon"].to_numpy(),
                           "variable": var, "value": values})  # fmt: skip
        with LOCK:
            df.to_csv(path, mode="a", header=not path.exists(), index=False, float_format="%.5g")

    for var in args.fields:
        if var == "depth":
            ds = None
            for d in STATIC[0]:
                try:
                    import copernicusmarine

                    ds = copernicusmarine.open_dataset(dataset_id=d)
                    break
                except Exception as exc:
                    print(f"depth: {d} not opened ({str(exc)[:120]})", flush=True)
            if ds is None:
                continue
            da = ds[STATIC[1]]
            todo = [s for s in groups if (s, "depth") not in done]
            for i, sid in enumerate(todo):
                g = groups[sid]
                v = da.sel(latitude=xr.DataArray(g["lat"].to_numpy(), dims="p"),
                           longitude=xr.DataArray(g["lon"].to_numpy(), dims="p"), method="nearest").values
                write(sid, "depth", g, v)
                if i % 200 == 0:
                    print(f"  depth {i:,}/{len(todo):,}", flush=True)
            continue
        ids, name = FIELDS[var]
        todo = {s: g for s, g in groups.items() if (s, var) not in done}
        for d in ids:
            if not todo:
                break
            try:
                ds = open_dataset(d)
            except Exception as exc:
                print(f"{var}: {d} not opened ({str(exc).splitlines()[0][:140]})", flush=True)
                continue
            t0d, t1d = pd.to_datetime(ds["time"].values[[0, -1]])
            cov = {s: g for s, g in todo.items() if g["time"].min() >= t0d.floor("D") and g["time"].max() <= t1d}
            print(f"{var}: {d} {t0d:%Y-%m-%d}..{t1d:%Y-%m-%d}: {len(cov):,} tracks", flush=True)
            da = ds[name]
            if "depth" in da.dims:
                da = da.isel(depth=0)
            started, cnt = clock.time(), [0]

            def work(item):
                sid, g = item
                for attempt in range(5):
                    try:
                        v = da.sel(time=xr.DataArray(g["time"].to_numpy(), dims="p"),
                                   latitude=xr.DataArray(g["lat"].to_numpy(), dims="p"),
                                   longitude=xr.DataArray(g["lon"].to_numpy(), dims="p"), method="nearest").values
                        break
                    except Exception as exc:
                        if attempt == 4:
                            print(f"    {sid} {var}: failed ({str(exc)[:80]})", flush=True)
                            return
                        clock.sleep(5 * (attempt + 1))
                if var == "sst":
                    v = v - 273.15
                write(sid, var, g, v)
                with LOCK:
                    cnt[0] += 1
                    if cnt[0] % 200 == 0:
                        print(f"    {var} {cnt[0]:,}/{len(cov):,} ({clock.time() - started:.0f} s)", flush=True)

            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                list(pool.map(work, list(cov.items())))
            todo = {s: g for s, g in todo.items() if s not in cov}
            ds.close()
        if todo:
            print(f"{var}: {len(todo):,} tracks not covered by any dataset", flush=True)


def probe() -> None:
    import copernicusmarine

    for d in LIGHT_CANDIDATES:
        try:
            ds = copernicusmarine.open_dataset(dataset_id=d)
            t = pd.to_datetime(ds["time"].values[[0, -1]]) if "time" in ds else ["?", "?"]
            print(f"OPENED {d}: {t[0]}..{t[1]}; variables {list(ds.data_vars)}", flush=True)
        except Exception as exc:
            print(f"not opened {d}: {str(exc).splitlines()[0][:140]}", flush=True)
    try:  # catalogue search for anything with PAR
        cat = copernicusmarine.describe(contains=["PAR"], disable_progress_bar=True)
        text = str(cat)
        print("describe(contains=['PAR']) mentions:", sorted(set(w for w in text.replace('"', " ").replace("'", " ").split()
                                                                if "par" in w.lower() and "oc" in w.lower()))[:40])  # fmt: skip
    except Exception as exc:
        print(f"describe() not available: {str(exc)[:140]}")


if __name__ == "__main__":
    main()
