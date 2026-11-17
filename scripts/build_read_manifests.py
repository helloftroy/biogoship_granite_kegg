#!/usr/bin/env python3
"""Raw-read download lists for the external test sets (RESULTS.md R31).

For every external sample that has a Sentinel-3 tile (sample_embeddings.csv), find its
sequencing runs on ENA and write data/external/<set>/reads/runs.csv: one row per run,
with the fastq URLs, sizes and MD5s from the ENA file report. cluster/download_reads.sbatch
downloads from these lists (no lookups on the cluster).

How runs are found:
    tara_pacific  the run accessions in samples.csv ("runs", ;-separated)
    emo_bon       ENA read runs of the sample's biosample accession
    jgi_flinkstrom
                  1. ENA read runs whose sample_title equals the JGI sample name
                  2. else the NCBI BioSample carrying IMG_genome_ID = the IMG id (often a
                     MAG record) -> its SRA_experiment_ID -> ENA read runs
Only metagenomic runs (strategy WGS, or OTHER as Tara Pacific labels them) are kept. Samples with no runs are listed in reads/unresolved.csv.

Usage (needs internet; run locally or on the service partition):
    python scripts/build_read_manifests.py [--sets tara_pacific emo_bon jgi_flinkstrom]
"""

from __future__ import annotations

import argparse
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship.tables import REPO_ROOT  # noqa: E402

EXT = REPO_ROOT / "data" / "external"
ENA = "https://www.ebi.ac.uk/ena/portal/api"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FIELDS = ("run_accession,experiment_accession,sample_accession,secondary_sample_accession,sample_title,fastq_ftp,"
          "fastq_bytes,fastq_md5,read_count,base_count,library_layout,library_strategy,library_source,instrument_model")  # fmt: skip


def get(url: str, tries: int = 5) -> str:
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "biogoship-read-manifest"}), timeout=60) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code == 400:  # unknown or malformed accession
                print(f"    ENA: bad request for {url[:120]}", flush=True)
                return ""
            if k == tries - 1:
                raise
            time.sleep(3 * (k + 1))
        except Exception as e:  # noqa: BLE001
            if k == tries - 1:
                raise
            print(f"    retry {k + 1}: {e}", flush=True)
            time.sleep(3 * (k + 1))
    return ""


def ena_runs(accession: str) -> pd.DataFrame:
    text = get(f"{ENA}/filereport?accession={accession}&result=read_run&fields={FIELDS}&format=tsv")
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str) if text.strip() else pd.DataFrame()


def ena_by_title(title: str) -> pd.DataFrame:
    q = urllib.parse.quote(f'sample_title="{title}"')
    text = get(f"{ENA}/search?result=read_run&query={q}&fields={FIELDS}&format=tsv&limit=0")
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str) if text.strip() else pd.DataFrame()


def jgi_via_biosample(img_id: str) -> pd.DataFrame:
    ids = json.loads(get(f"{EUTILS}/esearch.fcgi?db=biosample&term={img_id}&retmode=json&retmax=50"))["esearchresult"]["idlist"]
    time.sleep(0.4)
    if not ids:
        return pd.DataFrame()
    xml = get(f"{EUTILS}/efetch.fcgi?db=biosample&id={','.join(ids)}")
    time.sleep(0.4)
    blocks = [b for b in xml.split("<BioSample ")[1:] if re.search(rf'IMG_genome_ID">\s*{img_id}\s*<', b)]
    exps = sorted({e for b in blocks for e in re.findall(r'SRA_experiment_ID">\s*([SED]RX\d+)', b)})
    return pd.concat([ena_runs(e) for e in exps]) if exps else pd.DataFrame()


def resolve(set_name: str, s: pd.Series) -> tuple[pd.DataFrame, str]:
    if set_name == "tara_pacific":
        return pd.concat([ena_runs(r) for r in str(s["runs"]).split(";")]), "runs in samples.csv"
    if set_name == "emo_bon":
        return ena_runs(s["biosample"]), "biosample"
    runs = ena_by_title(s["sample_name"])
    if len(runs):
        return runs, "ENA sample_title"
    return jgi_via_biosample(s["sample_id"].removeprefix("JGI_")), "IMG id -> BioSample -> SRA experiment"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sets", nargs="+", default=["tara_pacific", "emo_bon", "jgi_flinkstrom"])
    args = parser.parse_args()
    for set_name in args.sets:
        folder = EXT / set_name
        samples = pd.read_csv(folder / "samples.csv", dtype=str)
        tiled = set(pd.read_csv(folder / "embeddings_sst_l4" / "sample_embeddings.csv", usecols=["sample_id"])["sample_id"])
        samples = samples[(samples["reason"] == "kept") & samples["sample_id"].isin(tiled)]
        print(f"{set_name}: {len(samples)} samples with a tile", flush=True)
        rows, missing = [], []
        for _, s in samples.iterrows():
            runs, how = resolve(set_name, s)
            if len(runs):
                runs = runs[runs["library_strategy"].isin(["WGS", "OTHER"]) & runs["library_source"].eq("METAGENOMIC")
                            & runs["fastq_ftp"].notna()]
            if not len(runs):
                missing.append({"sample_id": s["sample_id"], "sample_name": s.get("sample_name"), "biosample": s.get("biosample")})
                print(f"  {s['sample_id']}: no metagenomic runs found", flush=True)
                continue
            rows.append(runs.drop_duplicates("run_accession").assign(sample_id=s["sample_id"], section=s["section"], resolved_by=how))
        out = pd.concat(rows, ignore_index=True)
        files = out["fastq_ftp"].str.split(";")
        out["n_files"] = files.str.len()
        out["bytes_total"] = out["fastq_bytes"].str.split(";").map(lambda v: sum(int(x) for x in v))
        cols = ["sample_id", "section", "run_accession", "experiment_accession", "sample_accession", "library_layout",
                "instrument_model", "read_count", "base_count", "n_files", "bytes_total", "fastq_ftp", "fastq_bytes",
                "fastq_md5", "resolved_by", "sample_title"]  # fmt: skip
        (folder / "reads").mkdir(exist_ok=True)
        out[cols].to_csv(folder / "reads" / "runs.csv", index=False)
        pd.DataFrame(missing, columns=["sample_id", "sample_name", "biosample"]).to_csv(folder / "reads" / "unresolved.csv", index=False)
        print(f"  {out['sample_id'].nunique()} samples, {len(out)} runs, {out['bytes_total'].sum() / 1e12:.2f} TB, "
              f"{out['base_count'].astype(float).sum() / 1e12:.2f} Tbp; layouts {out['library_layout'].value_counts().to_dict()}, "
              f"files per run {out['n_files'].value_counts().to_dict()}; unresolved {len(missing)}", flush=True)


if __name__ == "__main__":
    main()
