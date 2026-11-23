#!/usr/bin/env python3
"""Per-cruise metadata status: what is confirmed, from where, and what is still missing ("?") (RESULTS.md R57).

Sources, combined:
  transfer/results/protocol/cruise_methods.csv      methods table supplied 2026-11-19 ("?" / "ND" kept as "?")
  transfer/results/library_metadata/library_metadata.csv   per library: flowcell, sequencer (ENA instrument_model
                                                    or original FASTQ header, R39-R41), library prep, extraction
  transfer/results/protocol/sample_protocol.csv     FAIRe workbook fields per sample (fill rate per cruise)
  data/results/mld_platform/mld_per_cruise.csv      modelled tiles per cruise (R43 tiles, C target)
The sequencer column replaces cruise_methods.csv's stale "NovaSeq X? (unverified)" (confirmed in R40/R41).

Output: data/results/metadata_status/cruise_metadata_status.csv (one row per cruise),
faire_fill_rate.csv (share of samples with each FAIRe field filled, per cruise), and
biogoship_metadata_check.xlsx, for the Bio-GO-SHIP team: sheet "Questions" (the open asks) and sheet
"Cruise methods" (our values, gaps highlighted, a column for corrections).

Usage:
    ~/miniforge3/bin/python3 scripts/metadata_status.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "results" / "metadata_status"
Q = "?"

# Hand-checked cruise-level values (sources in RESULTS.md R36-R41, R56). "?" = not known to us.
# cruise: (lab, sequencer, sequencer evidence, library prep, extraction, preservation, prefilter, filter, volume)
CRUISES = {
    "BVAL46-2011": ("whoi", "HiSeq 4000 + NovaSeq 6000", "ENA", "Nextera-style tagmentation (Illumina 20034197)", "enzymatic + Zymo", "lysis buffer, -20 C", "none reported", "0.22 um Sterivex", "2 L"),
    "AE1319-2013": ("whoi", "HiSeq 4000 + NovaSeq 6000", "ENA", "Nextera-style tagmentation (Illumina 20034197)", "enzymatic + Zymo", "lysis buffer, -20 C", "none reported", "0.22 um Sterivex", "2 L"),
    "NH1418-2014": ("whoi", "HiSeq 4000 + NovaSeq 6000", "ENA", "Nextera-style tagmentation (Illumina 20034197)", "enzymatic + Zymo", "lysis buffer, -20 C", "none reported", "0.22 um Sterivex", "2 L"),
    "I09N-2016": ("whoi", "NovaSeq 6000 (215 libs); HiSeq 4000 (27 libs)", "ENA", "Nextera-style tagmentation (Illumina 20034197)", "enzymatic + Zymo", "lysis buffer, -20 C", "none reported", "0.22 um Sterivex", "10 L"),
    "P18-2016": ("whoi", "NovaSeq 6000", "ENA", "Nextera-style tagmentation (Illumina 20034197)", "enzymatic + Zymo", "lysis buffer, -20 C", "none reported", "0.22 um Sterivex", "2 L"),
    "I07N-2018": ("whoi", "NovaSeq 6000", "ENA", "Nextera-style tagmentation (Illumina 20034197)", "enzymatic + Zymo", "lysis buffer, -20 C", "none reported", "0.22 um Sterivex", "2-4 L CTD; 4 L underway"),
    "AMT28-2018": ("whoi", "NovaSeq 6000", "ENA", "Nextera-style tagmentation (Illumina 20034197)", "enzymatic + Zymo", "lysis buffer, -20 C", "none reported", "0.22 um Sterivex", "2 L"),
    "C13.5/A13.5-2020": ("whoi", "NovaSeq 6000", "ENA", "Nextera-style tagmentation (Illumina 20034197)", "enzymatic + Zymo", "lysis buffer, -20 C", "none for DNA (30 um mesh for POM only)", "0.22 um Sterivex", "5-10 L"),
    "HLY2202-2022": ("uci", "NovaSeq 6000", "publication + flowcell format (no ENA, no header)", "Nextera DNA Flex, low-volume", "enzymatic + Zymo (original style)", Q, Q, "0.22 um Sterivex", "2-8 L"),
    "P02-2022": ("uci", "NovaSeq X", "FASTQ header (22YJMLLT4)", Q, Q + " (P02W shipped to UCI; final extraction not confirmed)", "lysis buffer, -80 C (some DNA/RNA Shield)", "30 um underway (P02W); ? (P02E, CTD)", "0.22 um Sterivex", "4-8 L"),
    "GOMECC4": ("noaa", "NovaSeq X", "FASTQ header (23LT2GLT4)", Q, "ZymoBIOMICS 96 MagBead, KingFisher (bead beating)", "none, -80 C", "none", "0.22 um Sterivex-GP PES", "~2 L (Niskin)"),
    "A16N-2023": ("whoi", "NovaSeq X", "FASTQ header (22NCJYLT4, 23LT2GLT4)", Q, "beads (kit ?)", "DNA/RNA Shield, -80 C (Protocol B; A excluded) - FAIRe says lysis buffer", Q, "0.22 um Sterivex", "~8 L (leg 2: ?)"),
    "I05-2023": ("noaa", "NovaSeq X", "FASTQ header (23LT2GLT4)", Q, "beads (kit ?)", "1 mL DNA/RNA Shield, -80 C - FAIRe says lysis buffer", Q, "0.22 um Sterivex", "~8 L"),
    "A13.5-2024": ("whoi", "NovaSeq X", "FASTQ header (23LT2GLT4)", Q, "beads (kit ?)", "1 mL DNA/RNA Shield, -80 C - FAIRe says lysis buffer", Q, "0.22 um Sterivex", "~8 L"),
    "I08S-2024": ("uci", "NovaSeq X", "FASTQ header (22YJMLLT4)", Q, "beads (kit ?)", "1 mL DNA/RNA Shield, -80 C", Q, "Sterivex, pore size ?", "~8 L"),
    "I09N-2025": ("uci", "NovaSeq X", "FASTQ header (22YJMLLT4)", Q, Q + " (beads?)", "1 mL DNA/RNA Shield, -80 C", Q, "0.22 um Sterivex", "~9 L"),
}  # fmt: skip
COLS = ["seq_lab", "sequencer", "sequencer_evidence", "library_prep", "extraction", "preservation", "prefilter", "filter", "volume"]
FAIRE = ["collection_method", "volume_filtered_l", "date_ext", "extract_plate", "dna_yield_ng", "nucl_acid_ext_kit", "size_frac",
         "samp_store_sol", "input_read_count"]  # fmt: skip
TILE_NAME = {"I09N-2016": "I09N", "P18-2016": "P18", "I07N-2018": "I07N", "AMT28-2018": "AMT28", "C13.5/A13.5-2020": "C13.5/A13.5"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.DataFrame.from_dict(CRUISES, orient="index", columns=COLS)
    lib = pd.read_csv(REPO / "transfer" / "results" / "library_metadata" / "library_metadata.csv", low_memory=False)
    real = lib[lib["in_real_targets"]]
    t["libraries_in_targets"] = real.groupby("cruise").size().reindex(t.index).fillna(0).astype(int)
    t["flowcells"] = lib.groupby("cruise")["flowcell_id"].agg(lambda v: "; ".join(f"{k} ({n})" for k, n in v.value_counts().items())).reindex(t.index)
    t["lane_recorded"] = real.groupby("cruise")["lane"].agg(lambda v: f"{v.notna().mean():.0%}").reindex(t.index).fillna(Q)
    tiles = pd.read_csv(REPO / "data" / "results" / "mld_platform" / "mld_per_cruise.csv", index_col="cruise")["tiles"]
    t["modelled_tiles"] = [int(tiles.get(TILE_NAME.get(c, c), 0)) for c in t.index]
    sp = pd.read_csv(REPO / "transfer" / "results" / "protocol" / "sample_protocol.csv", low_memory=False)
    fill = sp.groupby("cruise")[FAIRE].agg(lambda v: v.notna().mean()).reindex(t.index)
    fill.round(3).to_csv(OUT / "faire_fill_rate.csv")
    for f in ["collection_method", "volume_filtered_l", "extract_plate", "dna_yield_ng"]:
        t[f"{f}_per_sample"] = fill[f].map(lambda x: Q if pd.isna(x) or x == 0 else f"{x:.0%}")
    t.index.name = "cruise"
    order = ["modelled_tiles", "libraries_in_targets", "seq_lab", "flowcells", "sequencer", "sequencer_evidence", "library_prep",
             "extraction", "preservation", "prefilter", "filter", "volume", "collection_method_per_sample",
             "volume_filtered_l_per_sample", "extract_plate_per_sample", "dna_yield_ng_per_sample", "lane_recorded"]  # fmt: skip
    t[order].to_csv(OUT / "cruise_metadata_status.csv")
    print(t[order].to_string())
    workbook(t)


# Library prep (2023-2025) and same-DNA-on-both-sequencers asks dropped 2026-11-23: the sequencer term already
# captures that shift; the paper notes that prep and sequencer cannot be separated.
QUESTIONS = [
    ("Extraction, 2022-2025 cruises",
     "Which extraction kit and lysis (bead beating or enzymatic) was used for P02-2022, A16N-2023, I05-2023, "
     "A13.5-2024, I08S-2024 and I09N-2025? Was it the same protocol for all of them?",
     "Bead beating vs enzymatic lysis changes how well eukaryotes and Gram-positive cells are captured."),
    ("Preservation buffer conflict",
     "For A16N-2023 (Protocol B), I05-2023 and A13.5-2024, the methods table says 1 mL DNA/RNA Shield, but the FAIRe "
     "workbook (samp_store_sol) says lysis buffer. Which is right? And what was the buffer and storage temperature "
     "for HLY2202-2022?",
     "Preservation is one of the candidate causes of the shift between the older and newer cruises."),
    ("Pre-filtration",
     "Was any pre-filter or mesh used before the Sterivex on HLY2202-2022, P02-2022 (CTD / P02E), A16N-2023, "
     "I05-2023, A13.5-2024, I08S-2024 or I09N-2025? And what was the Sterivex pore size on I08S-2024?",
     "Pre-filtration removes larger eukaryotes and directly changes the eukaryote : prokaryote balance, which is "
     "the variable the satellite predicts best."),
    ("Per-sample records (lower priority)",
     "Are volume filtered, extraction plate / well and DNA yield recorded per sample for the cruises where the "
     "FAIRe workbook leaves them empty (see the sheet)?",
     "These let us test batch effects within a cruise."),
]


def workbook(t) -> None:
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill

    yellow = PatternFill("solid", fgColor="FFF2A8")
    orange = PatternFill("solid", fgColor="FFD3A8")
    bold = Font(name="Arial", bold=True)
    norm = Font(name="Arial")
    wrap = Alignment(wrap_text=True, vertical="top")
    wb = Workbook()
    q = wb.active
    q.title = "Questions"
    q.append(["Bio-GO-SHIP metagenome methods: open questions"])
    q.append(["Please answer in column D; the 'Cruise methods' sheet lists what we have now, with gaps highlighted."])
    q.append([])
    q.append(["#", "Topic", "Question", "Answer", "Why it matters"])
    for i, (topic, ques, why) in enumerate(QUESTIONS, 1):
        q.append([i, topic, ques, "", why])
    for row in q.iter_rows():
        for c in row:
            c.font, c.alignment = norm, wrap
    q["A1"].font = Font(name="Arial", bold=True, size=12)
    for c in q[4]:
        c.font = bold
    for col, w in zip("ABCDE", (4, 24, 70, 40, 50)):
        q.column_dimensions[col].width = w

    m = wb.create_sheet("Cruise methods")
    lab = {"seq_lab": "Lab (sequencing run list)", "libraries_in_targets": "Libraries in our analysis",
           "flowcells": "Flowcell (libraries)", "sequencer": "Sequencer", "sequencer_evidence": "How we know the sequencer",
           "library_prep": "Library prep", "extraction": "Extraction", "preservation": "Preservation / storage",
           "prefilter": "Pre-filter", "filter": "Filter", "volume": "Volume filtered",
           "volume_filtered_l_per_sample": "FAIRe: volume per sample", "extract_plate_per_sample": "FAIRe: extraction plate per sample",
           "dna_yield_ng_per_sample": "FAIRe: DNA yield per sample"}  # fmt: skip
    cols = list(lab)
    m.append(["Cruise", *[lab[c] for c in cols], "All correct? (Y / N)", "Corrections / notes"])
    for cruise, r in t.iterrows():
        m.append([cruise, *[r[c] for c in cols], "", ""])
    for row in m.iter_rows():
        for c in row:
            c.font, c.alignment = norm, wrap
            v = str(c.value or "")
            if c.row > 1 and ("?" in v):
                c.fill = yellow
            if c.row > 1 and "FAIRe says" in v:
                c.fill = orange
    for c in m[1]:
        c.font = bold
    m.freeze_panes = "B2"
    widths = [16, 10, 9, 26, 16, 24, 26, 26, 30, 22, 18, 16, 12, 12, 12, 12, 40]
    for i, w in enumerate(widths):
        m.column_dimensions[m.cell(1, i + 1).column_letter].width = w
    n = m.max_row + 2
    m.cell(n, 1, "Yellow: not known to us (?). Orange: our two sources disagree. FAIRe columns: share of samples with "
                 "the field filled in the FAIRe workbook ('?' = none). Sequencer 'FASTQ header': instrument ID LH... "
                 "(NovaSeq X series); 'ENA': instrument_model of the archived runs.").font = Font(name="Arial", italic=True)  # fmt: skip
    wb.save(OUT / "biogoship_metadata_check.xlsx")
    print("wrote", OUT / "biogoship_metadata_check.xlsx")


if __name__ == "__main__":
    main()
