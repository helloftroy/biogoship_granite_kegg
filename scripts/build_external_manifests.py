#!/usr/bin/env python3
"""Sample and site manifests for the external test sets (RESULTS.md R13, R14).

Three sets, each written to data/external/<set>/manifests/ in the format of
data/manifests/ (biogoship_samples.csv for build_tile_manifest.py, biogoship_sites.csv
for search_products.py and download_tiles.py), plus data/external/<set>/samples.csv
with every candidate sample, why it was or was not kept, and the source metadata.

  jgi_flinkstrom  Flinkstrom et al. 2024 (figshare 187989), IMG KO table already
                  available. Marine seawater, collected on or after the first OLCI
                  level-2 data, at most --max-depth m deep. Input: obs of
                  adata_raw_filtered_v2.h5ad (needs h5py).
  tara_pacific    Tara Pacific 0.2-3 um metagenomes (ENA PRJEB64422). Open-ocean
                  (O-SRF) and island-vicinity (I-SRF) surface samples only; reef and
                  coral-adjacent water (S-SRF, C-CSW) is left out. Input: the ENA run
                  table (source/ena_runs_PRJEB64422.tsv). UTC time from the sample title.
  emo_bon         EMO BON water samples with metaGOflow results on GitHub
                  (analysis-results-cluster-01/02-crate), any size fraction. Input: the
                  observatory logsheets (source/logsheets/<obs>_water_{sampling,observatory}.csv)
                  and the crate listings (source/analysis_results_cluster_0*.txt).

Times: Tara Pacific gives UTC to the minute. JGI and EMO BON give the date only; the
sample time is set to local solar noon (12:00 - lon/15 h UTC) and time_precision is
"day". Sites: one per set/date/position (4 decimals), so replicates and size fractions
of one sampling share a site.

Usage:
    python scripts/build_external_manifests.py            # all three
    python scripts/build_external_manifests.py --sets tara_pacific
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

EXTERNAL = REPO_ROOT / "data" / "external"
SENTINEL3_START = pd.Timestamp("2016-04-25", tz="UTC")  # first OLCI level-2 data
SETS = ["jgi_flinkstrom", "tara_pacific", "emo_bon"]


def solar_noon(dates: pd.Series, lon: pd.Series) -> pd.Series:
    return dates.dt.normalize() + pd.to_timedelta(12.0 - lon / 15.0, unit="h")


# ---------------------------------------------------------------- JGI

def read_h5ad_obs(path: Path) -> pd.DataFrame:
    """obs of an anndata file written by an old anndata (categoricals in obs/__categories)."""
    import h5py

    def dec(a):
        return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)

    with h5py.File(path, "r") as f:
        obs = f["obs"]
        cats = obs["__categories"] if "__categories" in obs else {}
        cols = {"img_genome_id": dec(obs[obs.attrs["_index"]][:])}
        for name in obs.attrs["column-order"]:
            v = obs[name][:]  # names with "/" are nested groups in HDF5; indexing by path still works
            if v.dtype.kind in "SO":
                v = dec(v)
            if name in cats:
                c = dec(cats[name][:])
                v = np.where(v >= 0, c[np.clip(v, 0, None)], None)
            cols[name] = v
    return pd.DataFrame(cols)


def jgi(max_depth: float) -> pd.DataFrame:
    d = read_h5ad_obs(EXTERNAL / "jgi_flinkstrom" / "adata_raw_filtered_v2.h5ad")
    date = pd.to_datetime(d["Sample Collection Date"], format="%Y/%m/%d", errors="coerce").dt.tz_localize("UTC")
    lat = pd.to_numeric(d["Latitude"], errors="coerce")
    lon = pd.to_numeric(d["Longitude"], errors="coerce")
    depth = pd.to_numeric(d["Depth In Meters"], errors="coerce")
    seawater = (d["Ecosystem Type"] == "Marine") & d["Habitat"].astype(str).str.lower().isin(["seawater", "marine"])
    study = d["Study Name"].astype(str)
    out = pd.DataFrame({
        "sample_id": "JGI_" + d["img_genome_id"].astype(str), "biosample": "", "time_utc": solar_noon(date, lon),
        "time_precision": "day", "lat": lat, "lon": lon, "depth_m": depth, "section": study.map(short_study),
        "size_fraction": "", "habitat": d["Habitat"], "study": study, "sample_name": d["Genome Name / Sample Name"],
    })  # fmt: skip
    out["reason"] = np.select(
        [~seawater, date.isna() | lat.isna() | lon.isna(), date < SENTINEL3_START, ~(depth <= max_depth)],
        ["not marine seawater", "no date or position", "before Sentinel-3", f"deeper than {max_depth:g} m"], "kept")
    return out


def short_study(name: str) -> str:
    known = {"Monterey Bay": "JGI_MontereyBay", "Marineland": "JGI_Marineland", "Southern Atlantic": "JGI_SAtlantic",
             "North Sea": "JGI_NorthSea", "Arctic Ocean": "JGI_Arctic"}  # fmt: skip
    for key, value in known.items():
        if key in name:
            return value
    return "JGI_" + re.sub(r"[^A-Za-z0-9]", "", name.title())[:24]


# ---------------------------------------------------------------- Tara Pacific

def tara_pacific(max_depth: float) -> pd.DataFrame:
    runs = pd.read_csv(EXTERNAL / "tara_pacific" / "source" / "ena_runs_PRJEB64422.tsv", sep="\t")
    agg = runs.groupby("sample_accession").agg(
        sample_title=("sample_title", "first"), lat=("lat", "first"), lon=("lon", "first"),
        runs=("run_accession", ";".join), base_count=("base_count", "sum"),
    ).reset_index()  # fmt: skip
    parts = agg["sample_title"].str.extract(r"TARA_SAMPLE_(?P<when>\d{8}T\d{4}Z)_(?P<station>[^_]+)_(?P<kind>[A-Z]-[A-Z]+)_")
    t = pd.to_datetime(parts["when"], format="%Y%m%dT%H%MZ", utc=True, errors="coerce")
    no_time = parts["when"].str.endswith("9999Z")  # time unknown: solar noon of that date
    day = pd.to_datetime(parts["when"].str[:8], format="%Y%m%d", utc=True, errors="coerce")
    t = t.where(~no_time, solar_noon(day, agg["lon"]))
    out = pd.DataFrame({
        "sample_id": "TARA_" + agg["sample_accession"], "biosample": agg["sample_accession"], "time_utc": t,
        "time_precision": np.where(no_time, "day", "minute"), "lat": agg["lat"], "lon": agg["lon"], "depth_m": np.nan,
        "section": "TARA_" + parts["station"].str[6:9],  # Inn: island nn, I00 = open ocean
        "size_fraction": "0.2-3", "habitat": parts["kind"], "station": parts["station"], "runs": agg["runs"],
        "base_count": agg["base_count"], "sample_name": agg["sample_title"],
    })  # fmt: skip
    # Surface plankton was pumped from about 2 m; no per-sample depth in the run table.
    out["reason"] = np.select(
        [~parts["kind"].isin(["O-SRF", "I-SRF"]), t.isna() | out["lat"].isna() | out["lon"].isna()],
        ["reef or coral-adjacent water", "no time or position"], "kept")
    return out


# ---------------------------------------------------------------- EMO BON

def emo_bon(max_depth: float) -> pd.DataFrame:
    src = EXTERNAL / "emo_bon" / "source"
    analysed = set()
    for path in sorted(src.glob("analysis_results_cluster_*.txt")):
        analysed |= set(path.read_text().split())
    rows = []
    for path in sorted((src / "logsheets").glob("*_water_sampling.csv")):
        obs_name = path.name.replace("_water_sampling.csv", "")
        s = pd.read_csv(path, dtype=str)
        site = pd.read_csv(src / "logsheets" / f"{obs_name}_water_observatory.csv", dtype=str)
        site = site[pd.to_numeric(site["latitude"], errors="coerce").notna()].iloc[0]
        s = s[s["source_mat_id"].notna()].copy()
        s["obs"], s["lat"], s["lon"] = obs_name.upper(), float(site["latitude"]), float(site["longitude"])
        rows.append(s)
    d = pd.concat(rows, ignore_index=True)
    date = pd.to_datetime(d["collection_date"], errors="coerce").dt.tz_localize("UTC")
    depth = pd.to_numeric(d["depth"], errors="coerce")
    out = pd.DataFrame({
        "sample_id": d["source_mat_id"], "biosample": d["ENA_accession_number_sample"], "time_utc": solar_noon(date, d["lon"]),
        "time_precision": "day", "lat": d["lat"], "lon": d["lon"], "depth_m": depth, "section": "EMOBON_" + d["obs"],
        "size_fraction": [size_fraction(a, b) for a, b in zip(d["size_frac_low"], d["size_frac_up"])], "habitat": "coastal water",
        "replicate": d["replicate"], "sample_name": d["source_mat_id_orig"],
    })  # fmt: skip
    out["reason"] = np.select(
        [~out["sample_id"].isin(analysed), date.isna(), date < SENTINEL3_START, ~(depth <= max_depth)],
        ["no metaGOflow results", "no date", "before Sentinel-3", f"deeper than {max_depth:g} m"], "kept")
    return out.drop_duplicates("sample_id")


def size_fraction(low, up) -> str:
    """'0.2-3' or '3-200'; some logsheets swap the bounds, a few lack one."""
    bounds = sorted(float(v) for v in (low, up) if pd.notna(v) and re.fullmatch(r"[0-9.]+", str(v)))
    return "-".join(f"{v:g}" for v in bounds) if len(bounds) == 2 else "unknown"


# ---------------------------------------------------------------- shared

def write(name: str, samples: pd.DataFrame) -> None:
    root = EXTERNAL / name
    (root / "manifests").mkdir(parents=True, exist_ok=True)
    samples = samples.copy()
    kept = samples["reason"] == "kept"
    key = (samples["time_utc"].dt.strftime("%Y-%m-%dT%H:%MZ") + "|" + samples["lat"].round(4).astype(str) + "|"
           + samples["lon"].round(4).astype(str))  # fmt: skip
    first = samples[kept].groupby(key[kept])["sample_id"].transform("first")
    samples["site_id"] = None
    samples.loc[kept, "site_id"] = "site_" + first.map(lambda s: re.sub(r"[^A-Za-z0-9_.-]", "_", s))
    samples["time_utc"] = samples["time_utc"].dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    samples.to_csv(root / "samples.csv", index=False)

    k = samples[kept]
    pipeline = pd.DataFrame({
        "sample_id": k["sample_id"], "biosample": k["biosample"], "collection_time_utc": k["time_utc"], "lat": k["lat"],
        "lon": k["lon"], "assay": "MG", "cruise_code": k["section"], "section": k["section"], "depth_m": k["depth_m"],
        "sample_method": "", "expocode": "", "ocean": "", "temperature_c": "", "phosphate_umol_kg": "",
        "nitrate_nitrite_umol_l": "", "site_id": k["site_id"],
    })  # fmt: skip
    pipeline.to_csv(root / "manifests" / "biogoship_samples.csv", index=False)
    sites = k.groupby("site_id").agg(
        time_utc=("time_utc", "first"), lat=("lat", "first"), lon=("lon", "first"), section=("section", "first"),
        n_samples=("sample_id", "size"), sample_ids=("sample_id", ";".join),
    ).reset_index()  # fmt: skip
    sites["n_metagenome_samples"] = sites["n_samples"]
    sites = sites[["site_id", "time_utc", "lat", "lon", "section", "n_samples", "n_metagenome_samples", "sample_ids"]]
    sites.sort_values("time_utc").to_csv(root / "manifests" / "biogoship_sites.csv", index=False)

    print(f"\n== {name}: {kept.sum()} of {len(samples)} samples kept, {len(sites)} sites, "
          f"{k['section'].nunique()} sections, {k['time_utc'].min()[:10]} to {k['time_utc'].max()[:10]}")
    print(samples["reason"].value_counts().to_string())
    print(k.groupby("section").size().to_string())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sets", nargs="+", choices=SETS, default=SETS)
    parser.add_argument("--max-depth", type=float, default=20.0, help="Deepest sample kept (m). Default 20.")
    args = parser.parse_args()
    builders = {"jgi_flinkstrom": jgi, "tara_pacific": tara_pacific, "emo_bon": emo_bon}
    for name in args.sets:
        write(name, builders[name](args.max_depth))


if __name__ == "__main__":
    main()
