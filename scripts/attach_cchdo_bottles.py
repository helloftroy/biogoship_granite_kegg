#!/usr/bin/env python3
"""Attach CCHDO bottle chemistry to the sequenced samples of the original cruises.

The sample manifest has temperature, phosphate and nitrate+nitrite from the curated
metadata table, with gaps (C13.5/A13.5 has no nutrients at all). The cruises'
CCHDO bottle files carry much more: salinity, oxygen, nitrate, nitrite, phosphate,
silicate, DIC, alkalinity and pH. Each sample is matched to bottles in one of two
ways, and the match is recorded so it can be filtered on:

  station       CTD samples with a station number (curated ctd_station_no_pub):
                the cast at that station closest in time; per variable, the
                bottle closest in depth within --depth-tolerance.
  nearest_cast  underway samples (or CTD samples without a station number): the
                cast closest in space-time, if within --max-km and --max-hours;
                per variable, the shallowest valid bottle at or above --max-depth.
                Underway water is drawn from about 5 m, so this is the surface
                water at the nearest station, not the sample's own water.

Output data/manifests/biogoship_samples_cchdo.csv: one row per sample of the
matched cruises with match_method, station/cast, distance_km, gap_hours,
depth_diff_m and the CCHDO variables (prefixed cchdo_). Sources are kept under
data/cchdo/<expocode>/.

Usage:
    python scripts/attach_cchdo_bottles.py
    python scripts/attach_cchdo_bottles.py --max-km 30 --max-hours 12
"""

from __future__ import annotations

import argparse
import io
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, haversine_km  # noqa: E402
from build_new_cruise_manifest import VARIABLES, cchdo_files, read_exchange_bottles  # noqa: E402
from build_sample_manifest import MERENLAB_URL  # noqa: E402

# Section -> CCHDO expocode, for the cruises with Sentinel-3 era samples.
EXPOCODES = {
    "AMT28": "74JC20180923",
    "I07N": "33RO20180423",
    "I09N": "33RR20160321",
    "P18": "33RO20161119",
    "C13.5/A13.5": "33RO20200321",
}


def main() -> None:
    args = parse_args()
    samples = pd.read_csv(args.manifest_dir / "biogoship_samples.csv")
    samples = samples[samples["section"].isin(EXPOCODES)].copy()
    samples["ctd_station"] = samples["sample_id"].map(station_numbers())
    rows = []
    for section, expocode in EXPOCODES.items():
        files = cchdo_files(expocode, args.data_dir / "cchdo" / expocode)
        if "bottle" not in files:
            print(f"{section}: no bottle file on CCHDO")
            continue
        bottles = read_exchange_bottles(files["bottle"]).assign(source_file=files["bottle"].name)
        bottles["station_key"] = bottles["station"].map(station_key)
        for _, sample in samples[samples["section"] == section].iterrows():
            rows.append({"sample_id": sample["sample_id"], "section": section, "expocode": expocode,
                         "sample_method": sample["sample_method"], **match_sample(sample, bottles, args)})  # fmt: skip
    out = pd.DataFrame(rows)
    path = args.manifest_dir / "biogoship_samples_cchdo.csv"
    out.to_csv(path, index=False)
    summarise(out, samples)
    print(f"\nWrote {path} ({len(out):,} samples)")


def station_numbers() -> dict[str, str]:
    """sample_id -> CTD station number from the curated metadata table (blank for underway samples)."""
    with urllib.request.urlopen(MERENLAB_URL, timeout=180) as resp:
        curated = pd.read_csv(io.BytesIO(resp.read()), low_memory=False)
    curated = curated.dropna(subset=["ctd_station_no_pub"]).drop_duplicates("sample_alias")
    return dict(zip(curated["sample_alias"], curated["ctd_station_no_pub"].map(station_key)))


def station_key(value) -> str:
    """'005', '5', 5.0 -> '5'; non-numeric station names are kept as written."""
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value).strip()


def match_sample(sample: pd.Series, bottles: pd.DataFrame, args: argparse.Namespace) -> dict:
    when = pd.Timestamp(sample["collection_time_utc"])
    casts = bottles.groupby(["station", "cast"]).agg(
        station_key=("station_key", "first"), lat=("lat", "mean"), lon=("lon", "mean"), time=("cast_time_utc", "first")
    ).reset_index()
    casts["gap_hours"] = (pd.to_datetime(casts["time"], utc=True) - when).dt.total_seconds() / 3600
    casts["distance_km"] = [haversine_km(sample["lat"], sample["lon"], lat, lon) for lat, lon in zip(casts["lat"], casts["lon"])]
    station = sample.get("ctd_station")
    if isinstance(station, str) and station and sample.get("sample_method") != "underway" and (casts["station_key"] == station).any():
        method = "station"
        here = casts[casts["station_key"] == station]
        cast = here.loc[here["gap_hours"].abs().idxmin()]
        depth = sample["depth_m"] if pd.notna(sample["depth_m"]) else 5.0
        tolerance = max(args.depth_tolerance, 0.2 * depth)
    else:
        method = "nearest_cast"
        # Space-time distance: 1 h of time counts like 2 km (a ship steams ~20 km/h; casts are hours apart).
        score = casts["distance_km"] + 2.0 * casts["gap_hours"].abs()
        cast = casts.loc[score.idxmin()]
        if cast["distance_km"] > args.max_km or abs(cast["gap_hours"]) > args.max_hours:
            return {"match_method": "none", "station": cast["station"], "cast": cast["cast"],
                    "distance_km": round(cast["distance_km"], 2), "gap_hours": round(cast["gap_hours"], 2)}  # fmt: skip
        depth, tolerance = None, None
    in_cast = bottles[(bottles["station"] == cast["station"]) & (bottles["cast"] == cast["cast"])]
    row = {"match_method": method, "station": cast["station"], "cast": cast["cast"],
           "distance_km": round(cast["distance_km"], 2), "gap_hours": round(cast["gap_hours"], 2)}  # fmt: skip
    diffs = []
    for name in VARIABLES:
        valid = in_cast[in_cast[name].notna()]
        if depth is None:  # surface water: shallowest valid bottle above the cut-off
            valid = valid[valid["depth_m"] <= args.max_depth]
            pick = valid.loc[valid["depth_m"].idxmin()] if len(valid) else None
        else:
            valid = valid[(valid["depth_m"] - depth).abs() <= tolerance]
            pick = valid.loc[(valid["depth_m"] - depth).abs().idxmin()] if len(valid) else None
        row[f"cchdo_{name}"] = pick[name] if pick is not None else np.nan
        if pick is not None:
            diffs.append(pick["depth_m"] - (depth if depth is not None else 0.0))
    row["bottle_depth_max_m"] = round(max(diffs, default=np.nan) + (depth or 0.0), 1) if diffs else np.nan
    return row


def summarise(out: pd.DataFrame, samples: pd.DataFrame) -> None:
    print("\nSamples matched per cruise (match_method):")
    print(out.groupby(["section", "match_method"]).size().unstack(fill_value=0).to_string())
    matched = out[out["match_method"] != "none"]
    if matched.empty:
        return
    print("\nNearest-cast distance and time gap (matched underway samples), median / 90th percentile:")
    near = matched[matched["match_method"] == "nearest_cast"]
    for column in ("distance_km", "gap_hours"):
        values = near[column].abs()
        print(f"  {column}: {values.median():.1f} / {values.quantile(0.9):.1f}" if len(values) else f"  {column}: -")
    print("\nSamples with a value, per cruise (all matched samples):")
    have = matched.groupby("section")[[f"cchdo_{v}" for v in VARIABLES]].count()
    print(have.rename(columns=lambda c: c.removeprefix("cchdo_")).to_string())
    check = matched.merge(samples[["sample_id", "temperature_c", "phosphate_umol_kg"]], on="sample_id")
    for ours, curated in (("cchdo_temperature_c", "temperature_c"), ("cchdo_phosphate_umol_kg", "phosphate_umol_kg")):
        both = check[[ours, curated, "match_method"]].dropna()
        for method, part in both.groupby("match_method"):
            diff = (part[ours] - part[curated]).abs()
            print(f"Check vs curated {curated} ({method}, n={len(part)}): median |diff| {diff.median():.3f}, 90th pct {diff.quantile(0.9):.3f}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data")
    parser.add_argument("--max-km", type=float, default=30.0, help="Underway samples: largest distance to the cast. Default 30 km.")
    parser.add_argument("--max-hours", type=float, default=12.0, help="Underway samples: largest time gap to the cast. Default 12 h.")
    parser.add_argument("--max-depth", type=float, default=10.0, help="Underway samples: deepest bottle used as surface water.")
    parser.add_argument("--depth-tolerance", type=float, default=5.0, help="CTD samples: depth mismatch allowed (or 20%% of depth).")
    return parser.parse_args()


if __name__ == "__main__":
    main()
