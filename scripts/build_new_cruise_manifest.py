#!/usr/bin/env python3
"""Candidate surface sampling events for newer Bio-GO-SHIP cruises (2022-2026).

The 2011-2020 cruises have sequencing runs with positions and times in ENA. For the
newer cruises sequencing may not be public yet, so this builds *candidate*
locations and times from the public hydrographic and bio-optical records, to be
matched to sequencing accessions later. Nothing here is a confirmed metagenomic
sample.

Sources (all public):
  CCHDO (cchdo.ucsd.edu): the latest bottle file (*_hy1.csv, WHP exchange) and
      cruise summary (*su.txt) for each leg. No login.
  SeaBASS, experiment BIO-GO-SHIP: biological sampling files (HPLC pigments; on
      I09N 2025 also IFCB imaging, particle backscatter and VSF). Downloads need a
      free NASA Earthdata login (EARTHDATA_USERNAME / EARTHDATA_PASSWORD; see
      cdse_credentials.sh.example). Without it, SeaBASS is skipped and flagged.

Outputs in data/new_cruises/ (sources kept under sources/ with their identifiers):
  cchdo_bottles.csv        every bottle, selected variables, all depths
  surface_candidates.csv   events at or above --max-depth (default 10 m):
                           hydrographic bottles (repeat trips at one depth merged)
                           and SeaBASS biological samples, each with its
                           timestamp source and missing-data flags
  cruise_coverage.csv      per cruise: legs, dates, counts, missing-data rates
  manifests/biogoship_sites.csv   one row per surface sampling event, in the
                           format search_products.py and download_tiles.py read;
                           cruises with satellite=no (A16S) are left out

Usage:
    python scripts/build_new_cruise_manifest.py
    python scripts/build_new_cruise_manifest.py --max-depth 10 --ifcb-files
"""

from __future__ import annotations

import argparse
import csv
import http.cookiejar
import io
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

CCHDO = "https://cchdo.ucsd.edu"
SEABASS = "https://seabass.gsfc.nasa.gov"
GETFILE = "https://oceandata.sci.gsfc.nasa.gov/getfile/"
URS = "https://urs.earthdata.nasa.gov"

# (section, year, CCHDO expocodes per leg, SeaBASS cruise names, fetch satellite data?)
CRUISES = [
    ("P02", 2022, ["33RR20220430", "33RR20220613"], ["P02"], True),
    ("A16N", 2023, ["33RO20230306", "33RO20230413"], ["A16N_L1_2023", "A16N_L2_2023"], True),
    ("I05", 2023, ["33RR20230722"], ["I05_2023"], True),
    ("A13.5", 2024, ["33H320240201"], [], True),
    ("I08S", 2024, ["325020240221"], ["I08S_2024"], True),
    ("I09N", 2025, ["325020250321"], ["I09N"], True),
    ("A16S", 2026, ["33RR20260309"], [], False),  # too recent: no satellite matching yet
]

# Exchange column -> output column. Bottle values preferred over CTD where both exist.
VARIABLES = {
    "temperature_c": ["CTDTMP"],
    "salinity_pss": ["SALNTY", "CTDSAL"],
    "oxygen_umol_kg": ["OXYGEN", "CTDOXY"],
    "nitrate_umol_kg": ["NITRAT", "NO2+NO3"],
    "phosphate_umol_kg": ["PHSPHT"],
    "silicate_umol_kg": ["SILCAT"],
}
BAD_FLAGS = {"4", "5", "9"}  # WOCE: 4 bad, 5 not reported, 9 not sampled


def main() -> None:
    args = parse_args()
    out = args.out_dir
    (out / "sources").mkdir(parents=True, exist_ok=True)
    opener = earthdata_opener()

    bottles, summaries, bio, coverage = [], [], [], []
    for section, year, expocodes, seabass_cruises, satellite in CRUISES:
        cruise = f"{section}_{year}"
        notes = []
        for leg, expocode in enumerate(expocodes, start=1):
            files = cchdo_files(expocode, out / "sources" / "cchdo" / expocode)
            if "bottle" in files:
                table = read_exchange_bottles(files["bottle"])
                bottles.append(table.assign(cruise=cruise, section=section, leg=leg, source_file=files["bottle"].name))
            else:
                notes.append(f"{expocode}: no bottle file on CCHDO yet")
            if "summary" in files:
                summaries.append(read_sumfile(files["summary"]).assign(expocode=expocode))
            else:
                notes.append(f"{expocode}: no summary file")
        for name in seabass_cruises:
            found, note = seabass_cruise(name, out / "sources" / "seabass" / name, opener, args)
            bio += [{**row, "cruise": cruise, "section": section} for row in found]
            notes += [note] if note else []
        coverage.append({"cruise": cruise, "section": section, "year": year, "legs": len(expocodes),
                         "expocodes": ";".join(expocodes), "seabass_cruises": ";".join(seabass_cruises),
                         "satellite": "yes" if satellite else "no", "notes": " | ".join(notes)})  # fmt: skip
        print(f"{cruise}: " + ("; ".join(notes) if notes else "ok"), flush=True)

    bottles_df = pd.concat(bottles, ignore_index=True) if bottles else pd.DataFrame()
    summary_df = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    bottles_df = attach_cast_events(bottles_df, summary_df)
    bottles_df.to_csv(out / "cchdo_bottles.csv", index=False)

    surface = surface_bottle_events(bottles_df, args.max_depth)
    bio_df = pd.DataFrame(bio)
    if len(bio_df):
        bio_df = bio_df[(bio_df["depth_m"].isna()) | (bio_df["depth_m"] <= args.max_depth)]
    candidates = pd.concat([surface, bio_df], ignore_index=True, sort=False)
    candidates["confirmed_metagenome"] = "no"  # never assumed; to be matched to sequencing accessions
    candidates.to_csv(out / "surface_candidates.csv", index=False)

    coverage_df = add_coverage_stats(pd.DataFrame(coverage), bottles_df, candidates, args.max_depth)
    coverage_df.to_csv(out / "cruise_coverage.csv", index=False)
    write_sites(candidates, coverage_df, out / "manifests" / "biogoship_sites.csv")
    print(coverage_df.drop(columns=["notes", "expocodes", "seabass_cruises"]).to_string(index=False))
    print(f"\nWrote {out / 'surface_candidates.csv'} ({len(candidates):,} rows)")


# ----------------------------------------------------------------------------
# CCHDO
# ----------------------------------------------------------------------------
def fetch(url: str, dest: Path | None = None, opener=None, retries: int = 5) -> bytes:
    for attempt in range(retries):
        try:
            with (opener or urllib.request.build_opener()).open(url, timeout=300) as resp:
                data = resp.read()
            if dest is not None:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
            return data
        except urllib.error.HTTPError as exc:
            if exc.code < 500 and exc.code != 429 or attempt == retries - 1:  # a missing file will not appear on retry
                raise
            time.sleep(10 * 2**attempt)
        except Exception:  # noqa: BLE001 - network errors are retried; the last one is raised
            if attempt == retries - 1:
                raise
            time.sleep(10 * 2**attempt)
    raise AssertionError("unreachable")


def cchdo_files(expocode: str, folder: Path) -> dict[str, Path]:
    """Latest bottle (exchange CSV) and summary files for a cruise, downloaded once."""
    folder.mkdir(parents=True, exist_ok=True)
    index = folder / "files.json"
    if index.exists():
        chosen = json.loads(index.read_text())
    else:
        cruise = cchdo_cruise_record(expocode, folder.parent / "cruises_all.json")
        chosen = {}
        for file_id in cruise.get("files", []):
            meta = json.loads(fetch(f"{CCHDO}/api/v1/file/{file_id}"))
            if meta.get("role") != "dataset":
                continue
            kind = {("bottle", "exchange"): "bottle", ("summary", "woce"): "summary"}.get(
                (meta.get("data_type"), meta.get("data_format"))
            )
            if kind and (kind not in chosen or file_id > chosen[kind]["id"]):  # newest upload wins
                chosen[kind] = {"id": file_id, "name": meta["file_name"], "path": meta["file_path"], "hash": meta.get("file_hash")}
        index.write_text(json.dumps(chosen, indent=2))
    paths = {}
    for kind, meta in chosen.items():
        dest = folder / meta["name"]
        if not dest.exists():
            fetch(CCHDO + meta["path"], dest)
        paths[kind] = dest
    return paths


def cchdo_cruise_record(expocode: str, cache: Path) -> dict:
    """The CCHDO record of one cruise, from the full cruise list (fetched once, about 12 MB)."""
    if not cache.exists():
        fetch(f"{CCHDO}/api/v1/cruise/all", cache)
    for cruise in json.loads(cache.read_text()):
        if cruise["expocode"] == expocode:
            return cruise
    raise RuntimeError(f"{expocode} not found on CCHDO")


def read_exchange_bottles(path: Path) -> pd.DataFrame:
    """WHP exchange bottle file -> one row per bottle with the variables of interest."""
    lines = [l for l in path.read_text(errors="replace").splitlines()[1:] if not l.startswith("#")]
    header, body = lines[0], lines[2:]  # line 1 holds units
    body = body[: next((i for i, l in enumerate(body) if l.strip() == "END_DATA"), len(body))]
    raw = pd.read_csv(io.StringIO("\n".join([header, *body])), dtype=str, skipinitialspace=True)
    raw.columns = [c.strip() for c in raw.columns]
    raw = raw.apply(lambda col: col.str.strip())
    out = pd.DataFrame({
        "expocode": raw["EXPOCODE"], "station": raw["STNNBR"], "cast": raw["CASTNO"],
        "sample": raw.get("SAMPNO"), "bottle": raw.get("BTLNBR"),
        "lat": num(raw["LATITUDE"]), "lon": num(raw["LONGITUDE"]),
        "pressure_dbar": num(raw.get("CTDPRS")),
        "cast_time_utc": utc(raw["DATE"], raw["TIME"]),
        "bottle_time_utc": utc(raw.get("BTL_DATE"), raw.get("BTL_TIME")),
    })  # fmt: skip
    out["depth_m"] = depth_from_pressure(out["pressure_dbar"], out["lat"])
    for name, columns in VARIABLES.items():
        value, source = pd.Series(np.nan, index=raw.index), pd.Series("", index=raw.index)
        for column in columns:
            if column not in raw:
                continue
            v = num(raw[column])
            flag = raw.get(f"{column}_FLAG_W", pd.Series("", index=raw.index)).fillna("")
            v[flag.isin(BAD_FLAGS)] = np.nan
            fill = value.isna() & v.notna()
            value[fill], source[fill] = v[fill], column
        out[name], out[f"{name}_source"] = value, source
    return out


def num(series) -> pd.Series:
    if series is None:
        return pd.Series(dtype=float)
    values = pd.to_numeric(series, errors="coerce")
    return values.where(values > -998)  # -999 marks missing in WHP exchange


def utc(dates, times) -> pd.Series:
    if dates is None:
        return pd.Series(dtype=object)
    times = (times if times is not None else pd.Series("0000", index=dates.index)).fillna("0000").str.zfill(4)
    stamp = pd.to_datetime(dates.fillna("") + times, format="%Y%m%d%H%M", errors="coerce", utc=True)
    return stamp.dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def depth_from_pressure(pressure: pd.Series, lat: pd.Series) -> pd.Series:
    """Depth (m) from pressure (dbar), UNESCO 1983 (Fofonoff and Millard)."""
    x = np.sin(np.radians(lat)) ** 2
    g = 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-6 * pressure
    return ((((-1.82e-15 * pressure + 2.279e-10) * pressure - 2.2512e-5) * pressure + 9.72659) * pressure) / g


SUM_LINE = re.compile(
    r"^(?P<expo>\S+)\s+(?P<sect>\S+)\s+(?P<stn>\S+)\s+(?P<cast>\d+)\s+(?P<type>\S+)\s+(?P<date>\d{6})\s+(?P<time>\d{4})\s+"
    r"(?P<code>BE|BO|EN|UN|MR|DE)\s+(?P<latd>\d+)\s+(?P<latm>[\d.]+)\s+(?P<lath>[NS])\s+(?P<lond>\d+)\s+(?P<lonm>[\d.]+)\s+(?P<lonh>[EW])"
)


def read_sumfile(path: Path) -> pd.DataFrame:
    """Cast events from a WOCE summary file: begin (BE), bottom (BO), end (EN), with time and position."""
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        m = SUM_LINE.match(line.strip())
        if not m:
            continue
        lat = (int(m["latd"]) + float(m["latm"]) / 60) * (1 if m["lath"] == "N" else -1)
        lon = (int(m["lond"]) + float(m["lonm"]) / 60) * (1 if m["lonh"] == "E" else -1)
        when = pd.to_datetime(m["date"] + m["time"], format="%m%d%y%H%M", utc=True)
        rows.append({"station": m["stn"], "cast": str(int(m["cast"])), "cast_type": m["type"], "code": m["code"],
                     "time": when.strftime("%Y-%m-%dT%H:%M:%SZ"), "lat": round(lat, 5), "lon": round(lon, 5)})  # fmt: skip
    return pd.DataFrame(rows)


def attach_cast_events(bottles: pd.DataFrame, summaries: pd.DataFrame) -> pd.DataFrame:
    if bottles.empty or summaries.empty:
        return bottles
    bottles = bottles.assign(cast=bottles["cast"].astype(str).str.lstrip("0").replace("", "0"))
    wide = summaries.pivot_table(index=["expocode", "station", "cast"], columns="code", values="time", aggfunc="first")
    wide = wide.reindex(columns=["BE", "BO", "EN"]).rename(columns={"BE": "cast_begin_utc", "BO": "cast_bottom_utc", "EN": "cast_end_utc"})
    types = summaries.groupby(["expocode", "station", "cast"])["cast_type"].first()
    return bottles.merge(wide.join(types).reset_index(), on=["expocode", "station", "cast"], how="left")


def surface_bottle_events(bottles: pd.DataFrame, max_depth: float) -> pd.DataFrame:
    """Bottles at or above max_depth; repeat trips at the same cast and depth (1 m bins) merged into one row."""
    if bottles.empty:
        return bottles
    shallow = bottles[bottles["depth_m"] <= max_depth].copy()
    shallow["depth_bin"] = shallow["depth_m"].round(0)
    keys = ["cruise", "section", "leg", "expocode", "station", "cast", "depth_bin"]
    numeric = ["lat", "lon", "pressure_dbar", "depth_m", *VARIABLES]
    agg = {c: "mean" for c in numeric}
    agg.update({c: "first" for c in ["cast_time_utc", "bottle_time_utc", "cast_begin_utc", "cast_bottom_utc",
                                      "cast_end_utc", "cast_type", "source_file"] if c in shallow})  # fmt: skip
    agg.update({f"{v}_source": "first" for v in VARIABLES})
    agg.update({"sample": lambda s: ";".join(s.dropna().astype(str)), "bottle": lambda s: ";".join(s.dropna().astype(str))})
    events = shallow.groupby(keys, dropna=False).agg(agg).reset_index()
    events["n_bottles"] = shallow.groupby(keys, dropna=False).size().to_numpy()
    events["event_type"] = "hydrographic_bottle"
    events["source"] = "CCHDO"
    events["time_utc"] = events["cast_time_utc"]
    events["time_source"] = "CTD cast (DATE/TIME in the bottle file; see cast_begin/bottom/end_utc from the summary)"
    events["event_id"] = events["expocode"] + "_s" + events["station"].astype(str) + "_c" + events["cast"].astype(str) + "_" + events["depth_bin"].astype(int).astype(str) + "m"
    for v in VARIABLES:
        events[f"missing_{v}"] = events[v].isna().astype(int)
    return events.drop(columns="depth_bin")


# ----------------------------------------------------------------------------
# SeaBASS
# ----------------------------------------------------------------------------
def earthdata_opener():
    user, password = os.getenv("EARTHDATA_USERNAME"), os.getenv("EARTHDATA_PASSWORD")
    if not user or not password:
        return None
    manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    manager.add_password(None, URS, user, password)
    return urllib.request.build_opener(
        urllib.request.HTTPBasicAuthHandler(manager), urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar())
    )


def seabass_cruise(name: str, folder: Path, opener, args: argparse.Namespace) -> tuple[list[dict], str]:
    """Biological sampling events from one SeaBASS cruise. Returns (rows, note)."""
    folder.mkdir(parents=True, exist_ok=True)
    page = fetch(f"{SEABASS}/cruise/{urllib.parse.quote(name)}").decode("utf-8", "replace")
    paths = sorted(set(re.findall(r"[A-Z_]+/[A-Za-z_]+/BIO-GO-SHIP/[^/\s<\"']+/archive/[^\s<\"']+\.sb", page)))
    (folder / "file_list.txt").write_text("\n".join(paths) + "\n")
    ifcb = [p for p in paths if p.endswith("_ifcb.sb")]
    others = [p for p in paths if p not in ifcb]
    rows = [ifcb_event(p, name) for p in ifcb]  # timestamps from file names; positions need the file
    wanted = others + (ifcb if args.ifcb_files else [])
    if opener is None:
        return rows, f"SeaBASS {name}: {len(paths)} files listed, not downloaded (set EARTHDATA_USERNAME/PASSWORD)"
    failed = 0
    for path in wanted:
        dest = folder / Path(path).name
        try:
            if not dest.exists():
                data = fetch(GETFILE + Path(path).name, opener=opener)
                if data.lstrip()[:1] == b"<":
                    raise RuntimeError("got an HTML page instead of the file (login failed?)")
                dest.write_bytes(data)
            rows += read_seabass(dest, name, path)
        except Exception as exc:  # noqa: BLE001 - keep going; report the count
            failed += 1
            print(f"  SeaBASS {Path(path).name}: {exc}", flush=True)
    return rows, f"SeaBASS {name}: {failed} of {len(wanted)} downloads failed" if failed else ""


def ifcb_event(path: str, cruise: str) -> dict:
    m = re.search(r"D(\d{8}T\d{6})", path)
    when = pd.to_datetime(m.group(1), format="%Y%m%dT%H%M%S", utc=True).strftime("%Y-%m-%dT%H:%M:%SZ") if m else ""
    return {"event_type": "biological_sample_ifcb", "source": "SeaBASS", "seabass_cruise": cruise, "event_id": Path(path).stem,
            "time_utc": when, "time_source": "IFCB sample time from the file name (UTC)", "source_file": path,
            "lat": np.nan, "lon": np.nan, "depth_m": np.nan, "position_note": "position in the file header; download with --ifcb-files"}  # fmt: skip


def read_seabass(path: Path, cruise: str, archive_path: str) -> list[dict]:
    """Rows of a SeaBASS .sb file: identifiers, time, position and depth of each record."""
    header, data = {}, []
    in_header = True
    for line in path.read_text(errors="replace").splitlines():
        if in_header:
            if line.lower().startswith("/end_header"):
                in_header = False
            elif line.startswith("/") and "=" in line:
                key, _, value = line[1:].partition("=")
                header[key.strip().lower()] = value.strip()
            continue
        if line.strip() and not line.startswith("!"):
            data.append(line)
    fields = [f.strip().lower() for f in header.get("fields", "").split(",")]
    delimiter = {"comma": ",", "tab": "\t"}.get(header.get("delimiter", "").lower())
    missing = header.get("missing", "-9999")
    kind = "biological_sample_hplc" if "hplc" in path.name.lower() else "biological_sample_" + path.stem.split("_")[-1].lower()
    rows = []
    for line in data:
        values = line.split(delimiter) if delimiter else line.split()
        rec = dict(zip(fields, [v.strip() for v in values]))
        rec = {k: (None if v == missing else v) for k, v in rec.items()}
        date = rec.get("date") or header.get("start_date")
        clock = rec.get("time") or header.get("start_time", "").replace("[gmt]", "").strip()
        when = pd.to_datetime(f"{date} {clock}", errors="coerce", utc=True) if date else pd.NaT
        identifiers = {k: rec[k] for k in ("station", "sample", "cast", "bottle", "sample_id", "cruise") if rec.get(k)}
        rows.append({
            "event_type": kind, "source": "SeaBASS", "seabass_cruise": cruise, "source_file": archive_path,
            "event_id": f"{path.stem}_{identifiers.get('sample', identifiers.get('station', len(rows)))}",
            "time_utc": when.strftime("%Y-%m-%dT%H:%M:%SZ") if when is not pd.NaT and not pd.isna(when) else "",
            "time_source": "biological sample time recorded in the SeaBASS file (UTC)",
            "lat": to_float(rec.get("lat")) if rec.get("lat") else to_float(header.get("north_latitude", "").split("[")[0]),
            "lon": to_float(rec.get("lon")) if rec.get("lon") else to_float(header.get("east_longitude", "").split("[")[0]),
            "depth_m": to_float(rec.get("depth")),
            "identifiers": json.dumps(identifiers),
        })  # fmt: skip
    return rows


def to_float(text) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


# ----------------------------------------------------------------------------
def add_coverage_stats(coverage: pd.DataFrame, bottles: pd.DataFrame, candidates: pd.DataFrame, max_depth: float) -> pd.DataFrame:
    rows = []
    for _, c in coverage.iterrows():
        b = bottles[bottles["cruise"] == c["cruise"]] if len(bottles) else bottles
        s = candidates[(candidates["cruise"] == c["cruise"]) & (candidates["event_type"] == "hydrographic_bottle")]
        bio = candidates[(candidates["cruise"] == c["cruise"]) & candidates["event_type"].astype(str).str.startswith("biological")]
        times = pd.to_datetime(b["cast_time_utc"], utc=True) if len(b) else pd.Series(dtype="datetime64[ns, UTC]")
        row = {**c.to_dict(),
               "start": times.min().strftime("%Y-%m-%d") if len(times) else "", "end": times.max().strftime("%Y-%m-%d") if len(times) else "",
               "stations": b[["expocode", "station"]].drop_duplicates().shape[0] if len(b) else 0,
               "casts": b[["expocode", "station", "cast"]].drop_duplicates().shape[0] if len(b) else 0,
               "bottles": len(b), f"surface_bottle_events_le_{max_depth:g}m": len(s),
               "surface_casts": s[["expocode", "station", "cast"]].drop_duplicates().shape[0] if len(s) else 0,
               "seabass_events": len(bio),
               "has_bottle_file": "yes" if len(b) else "no"}  # fmt: skip
        for v in VARIABLES:
            row[f"missing_{v}"] = round(float(s[v].isna().mean()), 3) if len(s) else math.nan
        rows.append(row)
    return pd.DataFrame(rows)


def write_sites(candidates: pd.DataFrame, coverage: pd.DataFrame, path: Path) -> None:
    """One site per surface sampling event with a position and time, for the satellite search."""
    allowed = set(coverage.loc[coverage["satellite"] == "yes", "cruise"])
    c = candidates[candidates["cruise"].isin(allowed) & candidates["lat"].notna() & candidates["time_utc"].astype(str).str.len().gt(0)].copy()
    hydro = c[c["event_type"] == "hydrographic_bottle"].sort_values("depth_m").drop_duplicates(["expocode", "station", "cast"])
    hydro = hydro.assign(site_id="new_" + hydro["expocode"] + "_s" + hydro["station"].astype(str) + "_c" + hydro["cast"].astype(str))
    bio = c[c["event_type"] != "hydrographic_bottle"].copy()
    bio = bio.assign(site_id="new_" + bio["event_id"].astype(str))
    sites = pd.concat([hydro, bio])[["site_id", "time_utc", "lat", "lon", "section", "cruise", "event_type"]]
    sites = sites.drop_duplicates("site_id").assign(n_samples=1, n_metagenome_samples=0, sample_ids="")
    path.parent.mkdir(parents=True, exist_ok=True)
    sites.sort_values("time_utc").to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"Sites for satellite matching: {len(sites):,} ({path})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "new_cruises")
    parser.add_argument("--max-depth", type=float, default=10.0, help="Deepest sample kept as surface, metres. Default 10.")
    parser.add_argument("--ifcb-files", action="store_true", help="Also download every IFCB file (thousands) for positions.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
