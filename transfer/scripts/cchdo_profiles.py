#!/usr/bin/env python3
"""Nutricline and deep-chlorophyll-maximum per CCHDO cast, for the modelled cruises (RESULTS.md R35).

Sources (public, CCHDO; cached under transfer/data/cchdo/<expocode>/):
  bottle file (WHP exchange hy1.csv)   nitrate (NITRAT, else NO2+NO3), flag 2/6 kept
  CTD files (WHP exchange ct1 zip)     2-dbar temperature, salinity, fluorescence (CTDFLUOR,
                                        when a fluorometer was fitted)

Per cast:
  no3_surface        nitrate of the shallowest bottle at <= 15 m (umol/kg)
  nitracline_1um_m   depth where nitrate first reaches 1 umol/kg, linear interpolation between
                     bottles; 0 when the surface is already >= 1; NaN when 1 is never reached
                     above 500 m or the profile has < 4 bottles above 500 m
  nitracline_grad_m  depth of the steepest nitrate increase (largest dNO3/dz between successive
                     bottles, 0-400 m; midpoint of that interval)
  bottles_0_200      bottles above 200 m (resolution of the above)
  From the CTD profile (down-cast, binned to 2 m, 5-bin running median):
  dcm_m              depth of the fluorescence maximum in 0-300 m
  dcm_ratio          fluorescence at the maximum / mean fluorescence 0-10 m (unitless; the
                     fluorometers are not cross-calibrated, so only ratios and depths are used)
  mld_ctd_m          mixed layer: first depth where temperature is 0.2 C below its value at 10 m
  When a cruise has no CTD fluorescence but bottle CTDFLUOR (flag 1, uncalibrated, accepted), dcm_m /
  dcm_ratio come from the bottles (dcm_source says which). A16N-2023 and A13.5-2024 have no
  fluorescence in either file (no DCM); A13.5-2024 has no CTD file on CCHDO (no CTD mixed layer).

Output: transfer/results/nutricline/station_profiles.csv (one row per cast).

Usage:
    ~/miniforge3/bin/python transfer/scripts/cchdo_profiles.py
"""

from __future__ import annotations

import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from build_new_cruise_manifest import CCHDO, cchdo_cruise_record, cchdo_files, depth_from_pressure, fetch, utc  # noqa: E402

# tile-table section name -> CCHDO expocodes (legs)
EXPOCODES = {
    "AMT28": ["74JC20180923"], "I07N": ["33RO20180423"], "I09N": ["33RR20160321"], "P18": ["33RO20161119"],
    "C13.5/A13.5": ["33RO20200321"], "P02-2022": ["33RR20220430", "33RR20220613"],
    "A16N-2023": ["33RO20230306", "33RO20230413"], "I05-2023": ["33RR20230722"], "A13.5-2024": ["33H320240201"],
    "I08S-2024": ["325020240221"], "I09N-2025": ["325020250321"],
}  # fmt: skip
DATA = REPO / "transfer" / "data" / "cchdo"
OK_FLAGS = {"2", "6", ""}


def ctd_zip(expocode: str, folder: Path) -> Path | None:
    """The newest CTD exchange zip of a cruise, downloaded once."""
    index = folder / "ctd_file.json"
    if index.exists():
        meta = json.loads(index.read_text())
    else:
        cruise = cchdo_cruise_record(expocode, DATA / "cruises_all.json")
        meta = {}
        for file_id in cruise.get("files", []):
            try:
                m = json.loads(fetch(f"{CCHDO}/api/v1/file/{file_id}"))
            except Exception:
                continue
            if m.get("role") == "dataset" and m.get("data_type") == "ctd" and m.get("data_format") == "exchange":
                if not meta or file_id > meta["id"]:
                    meta = {"id": file_id, "name": m["file_name"], "path": m["file_path"]}
        index.write_text(json.dumps(meta))
    if not meta:
        return None
    dest = folder / meta["name"]
    if not dest.exists():
        fetch(CCHDO + meta["path"], dest)
    return dest


def read_exchange(text: str) -> tuple[dict, pd.DataFrame]:
    lines = [l for l in text.splitlines()[1:] if not l.startswith("#")]
    head = {}
    i = 0
    if lines[0].startswith("NUMBER_HEADERS"):
        n = int(lines[0].split("=")[1])
        for l in lines[:n]:
            k, v = l.split("=", 1)
            head[k.strip()] = v.strip()
        i = n
    cols = [c.strip() for c in lines[i].split(",")]
    body = lines[i + 2:]
    body = body[: next((j for j, l in enumerate(body) if l.strip() == "END_DATA"), len(body))]
    df = pd.read_csv(io.StringIO("\n".join(body)), header=None, names=cols, dtype=str, skipinitialspace=True)
    return head, df.apply(lambda c: c.str.strip())


def key(v) -> str:
    """Station / cast number without padding ('001' and '1' match)."""
    v = str(v).strip()
    return v.lstrip("0") or "0"


def num(s):
    v = pd.to_numeric(s, errors="coerce")
    return v.where(v > -998)


def bottle_casts(path: Path) -> pd.DataFrame:
    _, raw = read_exchange(path.read_text(errors="replace"))
    out = []
    no3_col = "NITRAT" if "NITRAT" in raw else "NO2+NO3"
    fl_col = next((c for c in ["CTDFLUOR", "CHLORA"] if c in raw), None)
    raw["p"] = num(raw["CTDPRS"])
    raw["lat"], raw["lon"] = num(raw["LATITUDE"]), num(raw["LONGITUDE"])
    raw["no3"] = num(raw[no3_col]).where(raw.get(f"{no3_col}_FLAG_W", pd.Series("", index=raw.index)).fillna("").isin(OK_FLAGS))
    if fl_col:
        # flag 1 (uncalibrated) kept for fluorescence: only depths and ratios are used
        raw["fl"] = num(raw[fl_col]).where(raw.get(f"{fl_col}_FLAG_W", pd.Series("", index=raw.index)).fillna("").isin(OK_FLAGS | {"1"}))
    raw["t"] = utc(raw["DATE"], raw["TIME"])
    for (stn, cast), g in raw.groupby(["STNNBR", "CASTNO"]):
        g = g.assign(z=depth_from_pressure(g["p"], g["lat"])).sort_values("z")
        r = {"station": key(stn), "cast": key(cast), "lat": g["lat"].median(), "lon": g["lon"].median(), "time": g["t"].dropna().min(),
             "bottles_0_200": int((g["z"] <= 200).sum())}
        n = g[(g["z"] <= 500) & g["no3"].notna()][["z", "no3"]].groupby("z").mean().reset_index()
        surf = n[n["z"] <= 15]
        r["no3_surface"] = surf["no3"].iloc[0] if len(surf) else np.nan
        if len(n) >= 4:
            if n["no3"].iloc[0] >= 1:
                r["nitracline_1um_m"] = 0.0
            else:
                above = np.flatnonzero(n["no3"].to_numpy() >= 1)
                if len(above):
                    j = above[0]
                    z0, z1, v0, v1 = n["z"].iloc[j - 1], n["z"].iloc[j], n["no3"].iloc[j - 1], n["no3"].iloc[j]
                    r["nitracline_1um_m"] = z0 + (1 - v0) * (z1 - z0) / (v1 - v0) if v1 > v0 else z1
                else:
                    r["nitracline_1um_m"] = np.nan
            m = n[n["z"] <= 400]
            if len(m) >= 4:
                dz = np.diff(m["z"].to_numpy())
                gr = np.diff(m["no3"].to_numpy()) / np.where(dz > 0, dz, np.nan)
                k = np.nanargmax(gr)
                r["nitracline_grad_m"] = (m["z"].iloc[k] + m["z"].iloc[k + 1]) / 2
        if fl_col:
            f = g[(g["z"] <= 300) & g["fl"].notna()]
            if len(f) >= 4 and (f["z"] <= 15).any():
                top = f[f["z"] <= 15]["fl"].mean()
                k = f["fl"].idxmax()
                r["dcm_bottle_m"], r["dcm_bottle_ratio"] = f.loc[k, "z"], f.loc[k, "fl"] / top if top > 0 else np.nan
        out.append(r)
    return pd.DataFrame(out)


def ctd_casts(path: Path) -> pd.DataFrame:
    out = []
    with zipfile.ZipFile(path) as z:
        for name in z.namelist():
            if not name.lower().endswith(".csv"):
                continue
            try:
                head, d = read_exchange(z.read(name).decode(errors="replace"))
            except Exception:
                continue
            if "CTDPRS" not in d or "CTDTMP" not in d:
                continue
            lat = float(head.get("LATITUDE", "nan"))
            p = num(d["CTDPRS"])
            z_ = depth_from_pressure(p, pd.Series(lat, index=p.index))
            tmp = num(d["CTDTMP"])
            fl = num(d["CTDFLUOR"]) if "CTDFLUOR" in d else None
            if fl is not None and "CTDFLUOR_FLAG_W" in d:
                fl = fl.where(d["CTDFLUOR_FLAG_W"].fillna("").isin(OK_FLAGS | {"1"}))
            prof = pd.DataFrame({"z": z_, "t": tmp, "fl": fl if fl is not None else np.nan})
            prof = prof[prof["z"] <= 400].assign(b=lambda x: (x["z"] // 2) * 2 + 1).groupby("b").mean()
            r = {"station": key(head.get("STNNBR", "")), "cast": key(head.get("CASTNO", "")),
                 "ctd_lat": lat, "ctd_lon": float(head.get("LONGITUDE", "nan"))}
            t10 = prof["t"][(prof.index >= 8) & (prof.index <= 12)].mean()
            below = prof.index[(prof.index > 10) & (prof["t"] < t10 - 0.2)]
            r["mld_ctd_m"] = below[0] if len(below) else np.nan
            if prof["fl"].notna().sum() >= 20:
                f = prof["fl"].rolling(5, center=True, min_periods=3).median()
                f = f[f.index <= 300]
                top = f[f.index <= 10].mean()
                r["dcm_m"] = f.idxmax()
                r["dcm_ratio"] = f.max() / top if top and top > 0 else np.nan
            out.append(r)
    return pd.DataFrame(out)


def main() -> None:
    rows = []
    for section, codes in EXPOCODES.items():
        for code in codes:
            folder = DATA / code
            files = cchdo_files(code, folder)
            b = bottle_casts(files["bottle"]).assign(section=section, expocode=code)
            try:
                zp = ctd_zip(code, folder)
            except Exception as exc:
                print(f"{code}: CTD zip not fetched ({exc})")
                zp = None
            if zp is not None:
                c = ctd_casts(zp)
                if len(c):
                    b = b.merge(c, on=["station", "cast"], how="left")
            for col in ["dcm_m", "dcm_ratio", "mld_ctd_m", "dcm_bottle_m", "dcm_bottle_ratio"]:
                if col not in b:
                    b[col] = np.nan
            b["dcm_source"] = np.where(b["dcm_m"].notna(), "ctd", np.where(b["dcm_bottle_m"].notna(), "bottle", ""))
            b["dcm_m"] = b["dcm_m"].fillna(b["dcm_bottle_m"])
            b["dcm_ratio"] = b["dcm_ratio"].fillna(b["dcm_bottle_ratio"])
            print(f"{section:12s} {code}: {len(b):4d} casts; nitracline {b['nitracline_1um_m'].notna().sum():4d}; "
                  f"DCM {b['dcm_m'].notna().sum():4d} ({', '.join(f'{k} {v}' for k, v in b['dcm_source'].value_counts().items() if k)}); "
                  f"CTD MLD {b['mld_ctd_m'].notna().sum()}; bottles <= 200 m median {b['bottles_0_200'].median():.0f}", flush=True)
            rows.append(b)
    out = pd.concat(rows, ignore_index=True)
    dest = REPO / "transfer" / "results" / "nutricline"
    dest.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest / "station_profiles.csv", index=False)
    print(out.groupby("section")[["no3_surface", "nitracline_1um_m", "nitracline_grad_m", "dcm_m", "dcm_ratio", "mld_ctd_m"]].median().round(1).to_string())


if __name__ == "__main__":
    main()
