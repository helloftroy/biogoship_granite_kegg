#!/usr/bin/env python3
"""Build the Bio-GO-SHIP sample and site manifests from public metadata.

Sources (both public, fetched on each run):
  * ENA run report for BioProject PRJNA656268: every sequencing run with its
    BioSample, position and collection time. This is the authoritative sample list.
  * merenlab/public-marine-omics-metadata `BGS_patch2.csv`: curated extras for the
    metagenomes (section, depth, sampling method, temperature, nutrients).

Outputs in data/manifests/:
  biogoship_samples.csv   one row per BioSample (a filter of seawater)
  biogoship_sites.csv     one row per unique position and time; the unit that gets
                          matched to satellite data, since replicates share a site

Usage:
    python scripts/build_sample_manifest.py
"""

from __future__ import annotations

import argparse
import io
import sys
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import MANIFEST_DIR  # noqa: E402

BIOPROJECT = "PRJNA656268"
ENA_FIELDS = [
    "run_accession",
    "sample_accession",
    "experiment_accession",
    "sample_alias",
    "library_strategy",
    "library_source",
    "instrument_model",
    "collection_date",
    "lat",
    "lon",
    "read_count",
    "base_count",
    "first_public",
]
ENA_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}"
    "&result=read_run&fields={fields}&format=tsv&limit=0"
)
MERENLAB_URL = "https://raw.githubusercontent.com/merenlab/public-marine-omics-metadata/main/data/BGS_patch2.csv"

# Sample aliases start with a cruise code; map it to the hydrographic section name.
CRUISE_SECTIONS = {
    "A392013": "AE1319",
    "AMT2018": "AMT28",
    "B462011": "BVAL46",
    "C132020": "C13.5/A13.5",
    "I7N2018": "I07N",
    "I9N2016": "I09N",
    "N482014": "NH1418",
    "P182016": "P18",
}
CURATED_COLUMNS = {
    "depth": "depth_m",
    "sample_method_pub": "sample_method",
    "ccdho_expocode_pub": "expocode",
    "geo_loc_name": "ocean",
    "temperature_degC": "temperature_c",
    "phosphate_umolKg": "phosphate_umol_kg",
    "nitrate+nitrite_bottle_umol_l_1_pub": "nitrate_nitrite_umol_l",
}


def main() -> None:
    args = parse_args()
    runs = fetch_table(ENA_URL.format(accession=BIOPROJECT, fields=",".join(ENA_FIELDS)), sep="\t")
    curated = fetch_table(MERENLAB_URL, sep=",")
    samples = build_samples(runs, curated)
    sites = build_sites(samples)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    samples.to_csv(args.out_dir / "biogoship_samples.csv", index=False)
    sites.to_csv(args.out_dir / "biogoship_sites.csv", index=False)
    summarise(samples, sites, args.out_dir)


def fetch_table(url: str, sep: str) -> pd.DataFrame:
    with urllib.request.urlopen(url, timeout=180) as resp:
        return pd.read_csv(io.BytesIO(resp.read()), sep=sep, low_memory=False)


def build_samples(runs: pd.DataFrame, curated: pd.DataFrame) -> pd.DataFrame:
    runs = runs.copy()
    runs["kind"] = runs["library_strategy"].map({"WGS": "metagenome", "AMPLICON": "amplicon"}).fillna("other")

    def join_runs(kind: str) -> pd.Series:
        subset = runs[runs["kind"] == kind]
        return subset.groupby("sample_accession")["run_accession"].agg(lambda s: ";".join(sorted(s)))

    first = runs.sort_values("run_accession").groupby("sample_accession").first()
    for column in ("collection_date", "lat", "lon"):
        if (runs.groupby("sample_accession")[column].nunique() > 1).any():
            raise SystemExit(f"Runs of one BioSample disagree on {column}; inspect the ENA report.")

    samples = pd.DataFrame(
        {
            "sample_id": first["sample_alias"],
            "biosample": first.index,
            "collection_time_utc": pd.to_datetime(first["collection_date"], utc=True).dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "lat": first["lat"].astype(float),
            "lon": first["lon"].astype(float),
            "metagenome_runs": join_runs("metagenome"),
            "amplicon_runs": join_runs("amplicon"),
        }
    ).reset_index(drop=True)
    samples[["metagenome_runs", "amplicon_runs"]] = samples[["metagenome_runs", "amplicon_runs"]].fillna("")
    samples["cruise_code"] = samples["sample_id"].str.split("_").str[0]
    samples["section"] = samples["cruise_code"].map(CRUISE_SECTIONS)
    if samples["section"].isna().any():
        unknown = sorted(samples.loc[samples["section"].isna(), "cruise_code"].unique())
        raise SystemExit(f"Unknown cruise codes (add them to CRUISE_SECTIONS): {unknown}")

    extras = curated.rename(columns=CURATED_COLUMNS)[["biosample", *CURATED_COLUMNS.values()]]
    extras = extras.drop_duplicates("biosample")
    samples = samples.merge(extras, on="biosample", how="left")
    return samples.sort_values(["collection_time_utc", "sample_id"]).reset_index(drop=True)


def build_sites(samples: pd.DataFrame) -> pd.DataFrame:
    """Collapse samples that share a position and time, and tag each sample with its site."""
    keys = ["collection_time_utc", "lat", "lon"]
    first_sample = samples.groupby(keys)["sample_id"].transform("min")
    samples["site_id"] = "site_" + first_sample
    sites = (
        samples.groupby("site_id")
        .agg(
            time_utc=("collection_time_utc", "first"),
            lat=("lat", "first"),
            lon=("lon", "first"),
            section=("section", "first"),
            n_samples=("sample_id", "size"),
            n_metagenome_samples=("metagenome_runs", lambda s: int((s != "").sum())),
            sample_ids=("sample_id", lambda s: ";".join(sorted(s))),
        )
        .reset_index()
    )
    return sites.sort_values(["time_utc", "site_id"]).reset_index(drop=True)


def summarise(samples: pd.DataFrame, sites: pd.DataFrame, out_dir: Path) -> None:
    has_mg = samples["metagenome_runs"] != ""
    has_amp = samples["amplicon_runs"] != ""
    print(f"BioSamples: {len(samples):,} ({int(has_mg.sum()):,} with a metagenome, {int(has_amp.sum()):,} with amplicons)")
    print(f"Sites (unique position and time): {len(sites):,}")
    table = samples.groupby("section").agg(
        samples=("sample_id", "size"),
        sites=("site_id", "nunique"),
        start=("collection_time_utc", "min"),
        end=("collection_time_utc", "max"),
    )
    print(table.to_string())
    print(f"Wrote {out_dir / 'biogoship_samples.csv'}")
    print(f"Wrote {out_dir / 'biogoship_sites.csv'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-dir", type=Path, default=MANIFEST_DIR)
    return parser.parse_args()


if __name__ == "__main__":
    main()
