#!/usr/bin/env python3
"""One row per sequenced metagenome library, for the axis-1 batch diagnosis (RESULTS.md R39).

Rows: every MetaG_LibPrep_Illumina library in the FAIRe workbook (experimentRunMetadata), plus any
library in ko_counts_cov50_all.tsv without a workbook row. Columns, with where each comes from:

  lib_id                 KO-table library name (or the FAIRe lib_id if not in the KO table)
  faire_lib_id           FAIRe experimentRunMetadata lib_id
  sample_id              our sample key (as in data/targets), via build_ko_targets.library_samples
  samp_name, cruise      FAIRe
  in_ko_table, in_real_targets, excluded_protocol_a
  seq_run_id             FAIRe experimentRunMetadata
  sequencing_date        date part of seq_run_id ("20189999" = 2018, day unknown)
  flowcell_id            flowcell part of seq_run_id
  read_length            "PE150" etc. from seq_run_id; checked against ENA base_count / read_count
  lane                   from the library or file name (_L00n, -Ln-) where encoded
  instrument_id          FAIRe instrument term (OBI id)
  sequencer_model        ENA/SRA run record where an accession exists; else the published table;
                         else inferred from the flowcell-ID format ("...LT4" = NovaSeq X series,
                         unverified); see sequencer_source
  library_prep_kit       published table (2011-2020 cruises; HLY2202), else "not verified"
  extraction             published table, else "not verified"
  extraction_batch       FAIRe sampleMetadata extract_id plate (+ date_ext) where recorded
  preservation_protocol  cruise level (methods table, R35), refined per sample: A16N Protocol A /
                         Protocol B (cruise report), P02 "shield" notes
  total_reads            FAIRe input_read_count; ena_read_count (spots) where an accession exists
  ko_assigned_reads      sum of the library's row in ko_counts_cov50_all.tsv (reads assigned to KOs
                         at >= 50% reference coverage, split 1/n). Mapped-read counts (flagstat) are
                         not available to us: mapped_reads is left empty.
  sra_run                SRA run accession, where one exists
  source_of_information  per-field sources, compact

ENA run records are fetched once (study PRJNA656268) into transfer/data/ena_PRJNA656268_runs.tsv.
FASTQ headers: the ENA copies are rewritten (@SRRxxxx.n n/1), so instrument / run / flowcell / lane
cannot be read from them; the original files would be needed.

Output: transfer/results/library_metadata/library_metadata.csv and cruise_summary.csv.

Usage:
    ~/miniforge3/bin/python transfer/scripts/library_metadata.py
"""

from __future__ import annotations

import argparse
import io
import re
import subprocess
import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from biogoship.tables import MANIFEST_DIR  # noqa: E402
from build_ko_targets import library_samples  # noqa: E402

warnings.filterwarnings("ignore")
ENA = ("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession=PRJNA656268"
       "&fields=run_accession,sample_accession,instrument_model,library_name,library_strategy,library_layout,"
       "read_count,base_count,run_alias&limit=0&format=tsv")  # fmt: skip

# cruise (FAIRe short_name) -> (extraction, library prep, sequencer, confidence), from the tables supplied 2026-11-19/20
DOCUMENTED = {
    "AE1319-2013": ("enzymatic lysis + Zymo purification", "modified Nextera tagmentation, Illumina 20034197, custom dual indexes", "HiSeq 4000 or NovaSeq 6000", "published"),
    "BVAL46-2011": ("enzymatic lysis + Zymo purification", "modified Nextera tagmentation, Illumina 20034197, custom dual indexes", "HiSeq 4000 or NovaSeq 6000", "published"),
    "NH1418-2014": ("enzymatic lysis + Zymo purification", "modified Nextera tagmentation, Illumina 20034197, custom dual indexes", "HiSeq 4000 or NovaSeq 6000", "published"),
    "I09N-2016": ("enzymatic lysis + Zymo purification", "modified Nextera tagmentation, Illumina 20034197, custom dual indexes", "HiSeq 4000 or NovaSeq 6000", "published"),
    "P18-2016": ("enzymatic lysis + Zymo purification", "modified Nextera tagmentation, Illumina 20034197, custom dual indexes", "NovaSeq 6000", "published"),
    "I07N-2018": ("enzymatic lysis + Zymo purification", "modified Nextera tagmentation, Illumina 20034197, custom dual indexes", "NovaSeq 6000", "published"),
    "AMT28-2018": ("enzymatic lysis + Zymo purification", "modified Nextera tagmentation, Illumina 20034197, custom dual indexes", "NovaSeq 6000", "published"),
    "C13.5/A13.5-2020": ("enzymatic lysis + Zymo purification", "modified Nextera tagmentation, Illumina 20034197, custom dual indexes", "NovaSeq 6000", "published"),
    "HLY2202-2022": ("original-style enzymatic lysis + Zymo purification", "Nextera DNA Flex, low-volume, custom dual indexes", "NovaSeq 6000", "published"),
    "P02-2022": ("P02W shipped to UC Irvine; final extraction not confirmed", "not verified", "not verified", "partial"),
    "A16N-2023": ("NOAA AOML / Northern Gulf Institute; workflow not verified", "not verified", "not verified", "collection documented"),
    "I05-2023": ("not verified", "not verified", "not verified", "collection documented"),
    "A13.5-2024": ("not verified", "not verified", "not verified", "collection documented"),
    "I08S-2024": ("not verified", "not verified", "not verified", "collection documented"),
    "I09N-2025": ("not verified", "not verified", "not verified", "collection documented"),
    "GOMECC4": ("ZymoBIOMICS 96 MagBead DNA kit, KingFisher (FAIRe workbook)", "not verified", "not verified", "workbook"),
}  # fmt: skip
PRESERVATION = {
    "AE1319-2013": "1.62 mL lysis buffer, -20 C", "BVAL46-2011": "1.62 mL lysis buffer, -20 C", "NH1418-2014": "1.62 mL lysis buffer, -20 C",
    "I09N-2016": "1.62 mL lysis buffer, -20 C", "P18-2016": "1.62 mL lysis buffer, -20 C", "I07N-2018": "1.62 mL lysis buffer, -20 C",
    "AMT28-2018": "1.62 mL lysis buffer, -20 C", "C13.5/A13.5-2020": "1.62 mL lysis buffer, -20 C", "HLY2202-2022": "not verified",
    "P02-2022": "lysis buffer, -80 C (some DNA/RNA Shield)", "A16N-2023": "Protocol B: bashing beads + 1 mL DNA/RNA Shield, -80 C",
    "I05-2023": "1 mL DNA/RNA Shield, -80 C", "A13.5-2024": "1 mL DNA/RNA Shield, -80 C", "I08S-2024": "1 mL DNA/RNA Shield, -80 C",
    "I09N-2025": "1 mL DNA/RNA Shield, -80 C", "GOMECC4": "none, -80 C",
}  # fmt: skip


def strip(token: str) -> str:
    t = re.sub(r"(_R[12])?(_001)?\.fastq\.gz$", "", str(token).strip())
    return re.sub(r"-READ[12]-Sequences\.txt\.gz$", "", t)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", type=Path, default=REPO / "ko_counts_cov50_all.tsv")
    p.add_argument("--workbook", type=Path, default=REPO / "FAIRe-NOAA_uci-ess-biogoship.xlsx")
    p.add_argument("--faire-samples", type=Path, default=REPO / "data" / "faire" / "faire_samples.csv")
    p.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    p.add_argument("--out-dir", type=Path, default=REPO / "transfer" / "results" / "library_metadata")
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # KO table: library names and KO-assigned reads (row sums), read without loading all 18,836 columns at once
    sums = {}
    with open(args.counts) as fh:
        fh.readline()
        for line in fh:
            name, rest = line.split("\t", 1)
            sums[name] = float(np.fromstring(rest, sep="\t").sum())
    ko_libs = pd.Index(list(sums))
    link = library_samples(ko_libs, args)

    er = pd.read_excel(args.workbook, sheet_name="experimentRunMetadata", header=2)
    er = er[er["assay_name"] == "MetaG_LibPrep_Illumina"].reset_index(drop=True)
    er["sra_run"] = er["associatedSequences"].astype(str).str.extract(r"sra/(SRR\d+)")[0]
    token = {}
    for i, r in er.iterrows():
        for c in ["lib_id", "filename", "filename2", "sra_run"]:
            if pd.notna(r[c]):
                token.setdefault(strip(r[c]), i)
    def er_row(lib):
        for t in (lib, re.sub(r"_meta.*$", "_meta", lib), re.sub(r"_SA_L\d+$", "", lib)):
            if t in token:
                return token[t]
        return None
    link["er_row"] = [er_row(l) for l in link.index]
    say = print
    say(f"KO-table libraries: {len(ko_libs):,}; matched to a FAIRe run row: {link['er_row'].notna().sum():,}")

    # ENA
    ena_path = REPO / "transfer" / "data" / "ena_PRJNA656268_runs.tsv"
    if not ena_path.exists():
        ena_path.parent.mkdir(parents=True, exist_ok=True)
        ena_path.write_bytes(urllib.request.urlopen(ENA, timeout=300).read())
    ena = pd.read_csv(ena_path, sep="\t")
    ena = ena[ena["library_strategy"] == "WGS"]
    by_run = ena.set_index("run_accession")
    by_name = ena.drop_duplicates("library_name").set_index("library_name")

    # rows: FAIRe libraries + KO-table libraries without a FAIRe row
    rows = []
    matched = set(link["er_row"].dropna().astype(int))
    ko_by_er = {int(v): k for k, v in link["er_row"].dropna().items()}
    for i, r in er.iterrows():
        rows.append({"lib_id": ko_by_er.get(i, r["lib_id"]), "faire_lib_id": r["lib_id"], "er": i, "in_ko_table": i in matched})
    for lib in link.index[link["er_row"].isna()]:
        rows.append({"lib_id": lib, "faire_lib_id": None, "er": None, "in_ko_table": True})
    t = pd.DataFrame(rows)

    sm = pd.read_excel(args.workbook, sheet_name="sampleMetadata", header=2)
    sm = sm[~sm["samp_name"].astype(str).str.startswith("#")].set_index("samp_name")
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    excl = set(pd.read_csv(REPO / "data" / "manifests" / "excluded_samples.csv")["sample_key"])

    out = []
    for _, r in t.iterrows():
        e = er.loc[int(r["er"])] if pd.notna(r["er"]) else None
        lk = link.loc[r["lib_id"]] if r["lib_id"] in link.index else None
        samp = (lk["samp_name"] if lk is not None and pd.notna(lk["samp_name"]) else (e["samp_name"] if e is not None else None))
        sample_id = lk["sample_key"] if lk is not None else samp
        s = sm.loc[samp] if samp in sm.index else None
        cruise = (lk["cruise"] if lk is not None and pd.notna(lk["cruise"]) else (s["short_name"] if s is not None else None))
        run = e["seq_run_id"] if e is not None else None
        m = re.match(r"(\d{8})_([A-Z0-9]+)_MetaG_([A-Z]+\d+)", str(run)) if run is not None else None
        date, flowcell, layout = (m.group(1), m.group(2), m.group(3)) if m else (None, None, None)
        names = " ".join(str(x) for x in [r["lib_id"], e["filename"] if e is not None else ""])
        lane = re.search(r"_L00(\d)|-L(\d)-", names)
        lane = (lane.group(1) or lane.group(2)) if lane else None
        sra = e["sra_run"] if e is not None and pd.notna(e["sra_run"]) else (r["lib_id"] if str(r["lib_id"]).startswith("SRR") else None)
        en = by_run.loc[sra] if sra is not None and sra in by_run.index else (by_name.loc[r["lib_id"]] if r["lib_id"] in by_name.index else None)
        if en is not None and sra is None:
            sra = en.name
        doc = DOCUMENTED.get(cruise, ("not verified", "not verified", "not verified", "none"))
        if en is not None:
            seq, seq_src = en["instrument_model"], "ENA/SRA run record"
        elif doc[2] != "not verified":
            seq, seq_src = doc[2], "published table (cruise level)"
        elif flowcell and flowcell.endswith("LT4"):
            seq, seq_src = "NovaSeq X series? (flowcell ID format ...LT4; unverified)", "inferred from flowcell ID"
        elif flowcell and re.search(r"DS[XY][A-Z0-9]$|DRXX$", flowcell):
            seq, seq_src = "NovaSeq 6000? (flowcell ID format; unverified)", "inferred from flowcell ID"
        else:
            seq, seq_src = "unknown", ""
        read_len = layout
        if en is not None and en["read_count"] > 0:
            read_len = f"{layout or ''} (ENA {en['base_count'] / en['read_count']:.0f} bases/spot)"
        note = str(s["samp_collect_notes"]) if s is not None else ""
        pres = PRESERVATION.get(cruise, "not verified")
        pres_src = "methods table (cruise level)"
        if sample_id in excl or "Protocol A" in note:
            pres, pres_src = "Protocol A: no bashing beads, 1.6 mL lysis buffer, -80 C", "FAIRe note + A16N cruise report"
        elif note.strip().lower() == "shield":
            pres, pres_src = "DNA/RNA Shield (P02 buffer comparison)", "FAIRe note"
        ext_id = str(s["extract_id"]) if s is not None else ""
        plate = re.match(r"(Plate\d+)", ext_id)
        batch = (plate.group(1) + (f" ({s['date_ext']})" if pd.notna(s.get("date_ext")) and str(s["date_ext"]) != "not applicable" else "")) if plate else None
        out.append({
            "lib_id": r["lib_id"], "faire_lib_id": r["faire_lib_id"], "sample_id": sample_id, "samp_name": samp, "cruise": cruise,
            "in_ko_table": bool(r["in_ko_table"]), "in_real_targets": sample_id in ks.index, "excluded_protocol_a": sample_id in excl,
            "seq_run_id": run, "sequencing_date": date, "flowcell_id": flowcell, "read_length": read_len, "lane": lane,
            "instrument_id": re.search(r"OBI_\d+", str(e["instrument"])).group(0) if e is not None and re.search(r"OBI_\d+", str(e["instrument"])) else None,
            "sequencer_model": seq, "sequencer_source": seq_src, "library_prep_kit": doc[1], "extraction": doc[0],
            "documentation": doc[3], "extraction_batch": batch, "preservation_protocol": pres,
            "total_reads": pd.to_numeric(e["input_read_count"], errors="coerce") if e is not None else np.nan,
            "ena_read_count": en["read_count"] if en is not None else np.nan, "mapped_reads": np.nan,
            "ko_assigned_reads": sums.get(r["lib_id"], np.nan), "sra_run": sra,
            "source_of_information": "; ".join(x for x in [
                "lib/sample: KO table + FAIRe" if r["in_ko_table"] else "lib/sample: FAIRe (not in KO table)",
                "run/date/flowcell/read length: FAIRe seq_run_id" if run else "",
                f"sequencer: {seq_src}" if seq_src else "", f"prep/extraction: {doc[3]}",
                f"preservation: {pres_src}", "extraction batch: FAIRe extract_id" if batch else "",
                "total_reads: FAIRe input_read_count" if e is not None else "", "lane: library/file name" if lane else ""] if x),
        })  # fmt: skip
    df = pd.DataFrame(out)
    df.to_csv(args.out_dir / "library_metadata.csv", index=False)

    summ = df.groupby("cruise", dropna=False).agg(
        libraries=("lib_id", "size"), in_ko_table=("in_ko_table", "sum"), in_real_targets=("in_real_targets", "sum"),
        flowcells=("flowcell_id", lambda v: "; ".join(f"{k} ({n})" for k, n in v.value_counts().items())),
        sequencer=("sequencer_model", lambda v: "; ".join(f"{k} ({n})" for k, n in v.value_counts().items())),
        lanes=("lane", lambda v: ", ".join(sorted(v.dropna().unique()))), sra_runs=("sra_run", lambda v: v.notna().sum()),
        extraction_batches=("extraction_batch", lambda v: v.nunique()), prep=("library_prep_kit", "first"))
    summ.to_csv(args.out_dir / "cruise_summary.csv")
    pd.set_option("display.width", 250, "display.max_colwidth", 90)
    print(summ.to_string())
    print(f"\nWrote {len(df):,} libraries; unmatched KO-table libraries: {int((df['faire_lib_id'].isna()).sum())}")


if __name__ == "__main__":
    main()
