#!/usr/bin/env python3
"""Sample and site manifests from the Bio-GO-SHIP FAIRe-NOAA metadata workbook.

FAIRe-NOAA_uci-ess-biogoship.xlsx, tab sampleMetadata, lists every sequenced
sample (metagenomes _MG, metatranscriptomes _MT) of every cruise with its time,
position, depth, collection method, station/cast/bottle and the environmental
measurements (temperature, salinity, oxygen, nutrients, carbonate system,
chlorophyll). It replaces the candidate events built from CCHDO casts for the
newer cruises: about a third to half of their samples were taken from the
underway line between stations, where no cast is.

Outputs in data/faire/:
  sampleMetadata.csv            the tab as it is, for reference
  faire_samples.csv             one row per sample: identifiers, time, position,
                                depth, method, numeric environmental values (units
                                in <name>_unit), assay (MG/MT), site_id, and the
                                matching sample_id of data/manifests/biogoship_samples.csv
                                for the original cruises (same site within 1 h, 2 km)
  manifests/biogoship_sites.csv one row per unique time and position on the cruises
                                not yet matched to Sentinel-3, in the format
                                search_products.py and download_tiles.py read

Reading .xlsx needs openpyxl (pip install openpyxl). Run it where that is
installed and copy data/faire/ to the cluster.

Usage:
    python scripts/build_faire_manifest.py FAIRe-NOAA_uci-ess-biogoship.xlsx
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, haversine_km  # noqa: E402

# Cruises already in data/manifests (ENA-based) -> their section name there.
ORIGINAL = {
    "AE1319-2013": "AE1319", "AMT28-2018": "AMT28", "BVAL46-2011": "BVAL46", "C13.5/A13.5-2020": "C13.5/A13.5",
    "I07N-2018": "I07N", "I09N-2016": "I09N", "NH1418-2014": "NH1418", "P18-2016": "P18",
}  # fmt: skip
SENTINEL3_START = pd.Timestamp("2016-04-25", tz="UTC")  # first OLCI level-2 data
ENVIRONMENT = [
    "temp", "salinity", "chlorophyll", "light_intensity", "ph", "diss_oxygen", "diss_inorg_carb", "tot_alkalinity",
    "carbonate", "hydrogen_ion", "omega_arag", "pco2", "nitrate", "nitrite", "nitrate_plus_nitrite", "ammonium",
    "phosphate", "silicate", "pressure", "transmittance", "tot_depth_water_col",
]  # fmt: skip
IDENTIFIERS = [
    "samp_name", "short_name", "line_id", "eventDate", "verbatimEventDate", "collection_method", "samp_collect_device",
    "station_id", "ctd_cast_number", "ctd_bottle_number", "replicate_number", "size_frac_low", "size_frac",
    "biosample_accession", "materialsampleID", "samp_collect_notes",
]  # fmt: skip
MISSING = {"not applicable", "not provided", "missing", "not collected", "nan", ""}


def main() -> None:
    args = parse_args()
    raw = pd.read_excel(args.workbook, sheet_name="sampleMetadata", header=2, dtype=str)
    raw = raw[~raw["samp_name"].astype(str).str.startswith("#")]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw.to_csv(args.out_dir / "sampleMetadata.csv", index=False)

    samples = tidy(raw[raw["samp_category"] == "sample"])
    samples = link_original(samples, pd.read_csv(args.original_manifest / "biogoship_samples.csv"))
    samples, sites = build_sites(samples)
    samples.to_csv(args.out_dir / "faire_samples.csv", index=False)
    (args.out_dir / "manifests").mkdir(exist_ok=True)
    sites.to_csv(args.out_dir / "manifests" / "biogoship_sites.csv", index=False)
    summarise(samples, sites)
    print(f"\nWrote {args.out_dir / 'faire_samples.csv'} and {args.out_dir / 'manifests' / 'biogoship_sites.csv'}")


def tidy(raw: pd.DataFrame) -> pd.DataFrame:
    out = raw[[c for c in IDENTIFIERS if c in raw]].copy()
    out = out.rename(columns={"short_name": "cruise"})
    out["assay"] = raw["samp_name"].str.extract(r"_(MG|MT)$", expand=False).fillna("other")
    out["time_utc"] = pd.to_datetime(raw["eventDate"], errors="coerce", utc=True)
    out["lat"] = pd.to_numeric(raw["decimalLatitude"], errors="coerce")
    out["lon"] = pd.to_numeric(raw["decimalLongitude"], errors="coerce")
    out["depth_min_m"] = pd.to_numeric(raw["minimumDepthInMeters"], errors="coerce")
    out["depth_max_m"] = pd.to_numeric(raw["maximumDepthInMeters"], errors="coerce")
    for name in ENVIRONMENT:
        if name in raw:
            out[name] = pd.to_numeric(raw[name].where(~raw[name].astype(str).str.lower().isin(MISSING)), errors="coerce")
            if f"{name}_unit" in raw:
                out[f"{name}_unit"] = raw[f"{name}_unit"].where(out[name].notna())
    for column in IDENTIFIERS:
        if column in out:
            out[column] = out[column].where(~out[column].astype(str).str.lower().isin(MISSING))
    return out


def link_original(samples: pd.DataFrame, original: pd.DataFrame) -> pd.DataFrame:
    """sample_id of the ENA-based manifest for the same sample: same section, within 1 h and 2 km."""
    original = original.assign(t=pd.to_datetime(original["collection_time_utc"], utc=True))
    samples["original_sample_id"] = None
    for cruise, section in ORIGINAL.items():
        ours = original[original["section"] == section]
        for index, row in samples[(samples["cruise"] == cruise) & samples["time_utc"].notna()].iterrows():
            gap = (ours["t"] - row["time_utc"]).dt.total_seconds().abs().to_numpy() / 3600
            dist = np.array([haversine_km(row["lat"], row["lon"], a, b) for a, b in zip(ours["lat"], ours["lon"])])
            ok = (gap <= 1) & (dist <= 2)
            if ok.any():
                samples.at[index, "original_sample_id"] = ours["sample_id"].iloc[int(np.argmin(np.where(ok, gap + dist, np.inf)))]
    return samples


def build_sites(samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One site per unique time (minute) and position (4 decimals) on the cruises not yet in data/manifests."""
    usable = samples["time_utc"].notna() & samples["lat"].notna() & samples["lon"].notna()
    new = ~samples["cruise"].isin(list(ORIGINAL)) & samples["cruise"].notna()
    in_era = samples["time_utc"] >= SENTINEL3_START
    pick = samples[usable & new & in_era].copy()
    key = pick["time_utc"].dt.strftime("%Y-%m-%dT%H:%MZ") + "|" + pick["lat"].round(4).astype(str) + "|" + pick["lon"].round(4).astype(str)
    pick["site_key"] = key
    first = pick.groupby("site_key")["samp_name"].transform("first")
    pick["site_id"] = "site_" + first.str.replace(r"_(MG|MT)$", "", regex=True).map(lambda s: re.sub(r"[^A-Za-z0-9_.-]", "_", s))
    samples = samples.join(pick["site_id"])
    sites = pick.groupby("site_id").agg(
        time_utc=("time_utc", "first"), lat=("lat", "first"), lon=("lon", "first"), section=("cruise", "first"),
        n_samples=("samp_name", "size"), n_metagenome_samples=("assay", lambda s: int((s == "MG").sum())),
        sample_ids=("samp_name", lambda s: ";".join(s)),
    ).reset_index()  # fmt: skip
    sites["time_utc"] = sites["time_utc"].dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    samples["time_utc"] = samples["time_utc"].dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    sites = sites[["site_id", "time_utc", "lat", "lon", "section", "n_samples", "n_metagenome_samples", "sample_ids"]]
    return samples, sites.sort_values("time_utc")


def summarise(samples: pd.DataFrame, sites: pd.DataFrame) -> None:
    table = samples.groupby("cruise").agg(
        samples=("samp_name", "size"), MG=("assay", lambda s: int((s == "MG").sum())), MT=("assay", lambda s: int((s == "MT").sum())),
        no_time_or_position=("lat", lambda s: int(s.isna().sum())), at_or_above_10m=("depth_min_m", lambda s: int((s <= 10).sum())),
        start=("time_utc", lambda s: s.dropna().min()[:10] if s.notna().any() else ""),
        end=("time_utc", lambda s: s.dropna().max()[:10] if s.notna().any() else ""),
        in_original_manifest=("original_sample_id", lambda s: int(s.notna().sum())),
    )  # fmt: skip
    table = table.join(sites.groupby("section").size().rename("sites"))
    print(table.fillna({"sites": 0}).astype({"sites": int}).to_string())
    print(f"\nSites for satellite matching: {len(sites):,} (cruises not yet in data/manifests, sampled since {SENTINEL3_START:%Y-%m-%d})")
    have = samples[ENVIRONMENT].notna().sum()
    print("Samples with a value: " + ", ".join(f"{k} {v:,}" for k, v in have[have > 0].items()))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("workbook", type=Path, help="FAIRe-NOAA_uci-ess-biogoship.xlsx")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "faire")
    parser.add_argument("--original-manifest", type=Path, default=MANIFEST_DIR, help="Folder with the ENA-based biogoship_samples.csv.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
