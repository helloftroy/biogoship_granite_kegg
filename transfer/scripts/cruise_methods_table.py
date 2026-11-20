#!/usr/bin/env python3
"""Master cruise-methods table for the Bio-GO-SHIP DNA samples (RESULTS.md R35).

Transcribed from the methods table supplied on 2026-11-19, compiled from the cruise papers and
reports. Conventions there: "?" = not confirmed; "ND" = not described in the source (not
proof that the step was absent). Volumes refer to DNA samples. GOMECC4 is not in that
table; its row comes from the FAIRe workbook (sampleMetadata: 0.22 um Sterivex-GP, no
pre-filter, about 2 L from Niskin bottles, no preservative, -80 C, ZymoBIOMICS 96
MagBead kit on a KingFisher).

Coded columns added for testing (cruise level, as best known; "unknown" where the table
has ? or ND):
  preservation  lysis buffer | DNA/RNA Shield | mixed (paired protocols) | none | unknown
  extraction    enzymatic (lysozyme, proteinase K/SDS, Zymo clean-up) | beads | mixed | unknown
  prefilter     none reported | 30 um (underway) | unknown
  volume_l      typical DNA volume, litres (midpoint of the range given)
  store_c       storage temperature, C

Output: transfer/results/protocol/cruise_methods.csv (one row per cruise or leg), with
`cruise` matching the section names of the tile tables and `legs_by_date` giving the
split used for P02 and A16N (first day of the later leg).

Usage:
    ~/miniforge3/bin/python transfer/scripts/cruise_methods_table.py
"""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
COLS = ["row", "cruise", "year", "dna_prefilter", "dna_filter", "volume", "preservation_storage", "beads",
        "extraction_purification", "differences", "preservation", "extraction", "prefilter", "volume_l", "store_c",
        "legs_by_date", "source"]  # fmt: skip
T = "user table 2026-11-19"
ROWS = [
    ("BVAL46", "BVAL46-2011", "2011", "None reported", "0.22 um Sterivex", "2 L", "1.62 mL lysis buffer; -20 C", "No, enzymatic",
     "Lysozyme, proteinase K/SDS; Zymo Genomic DNA Clean & Concentrator", "Original protocol", "lysis buffer", "enzymatic", "none reported", 2, -20, "", T),
    ("AE1319", "AE1319-2013", "2013", "None reported", "0.22 um Sterivex", "2 L", "1.62 mL lysis buffer; -20 C", "No, enzymatic",
     "Same original extraction", "Original protocol", "lysis buffer", "enzymatic", "none reported", 2, -20, "", T),
    ("NH1418", "NH1418-2014", "2014", "None reported", "0.22 um Sterivex", "2 L", "1.62 mL lysis buffer; -20 C", "No, enzymatic",
     "Same original extraction", "Original protocol", "lysis buffer", "enzymatic", "none reported", 2, -20, "", T),
    ("I09N", "I09N", "2016", "None reported", "0.22 um Sterivex", "10 L", "1.62 mL lysis buffer; -20 C", "No, enzymatic",
     "Same original extraction", "Underway sampling", "lysis buffer", "enzymatic", "none reported", 10, -20, "", T),
    ("P18", "P18", "2016-17", "None reported", "0.22 um Sterivex", "2 L", "1.62 mL lysis buffer; -20 C", "No, enzymatic",
     "Same original extraction", "CTD sampling", "lysis buffer", "enzymatic", "none reported", 2, -20, "", T),
    ("I07N", "I07N", "2018", "None reported", "0.22 um Sterivex", "2-4 L CTD; 4 L underway", "1.62 mL lysis buffer; -20 C", "No, enzymatic",
     "Same original extraction", "Mixed underway/CTD", "lysis buffer", "enzymatic", "none reported", 3.5, -20, "", T),
    ("AMT28", "AMT28", "2018", "None reported", "0.22 um Sterivex", "2 L", "1.62 mL lysis buffer; -20 C", "No, enzymatic",
     "Same original extraction", "CTD sampling", "lysis buffer", "enzymatic", "none reported", 2, -20, "", T),
    ("C13.5/A13.5", "C13.5/A13.5", "2020", "None reported for DNA", "0.22 um Sterivex", "5-10 L", "1.62 mL lysis buffer; -20 C", "No, enzymatic",
     "Same original extraction", "Underway; POM sampling separately used 30 um mesh", "lysis buffer", "enzymatic", "none reported", 7.5, -20, "", T),
    ("P02W", "P02-2022", "2022", "30 um underway; CTD ?", "0.22 um Sterivex", "Up to 8 L", "1.62 mL lysis buffer; -80 C", "ND",
     "?", "Source-dependent pre-filtration", "lysis buffer", "unknown", "30 um (underway)", 8, -80, "leg 1 before 2022-06-13", T),
    ("P02E", "P02-2022", "2022", "?", "0.22 um Sterivex", "4-8 L", "Lysis buffer normally; some DNA/RNA Shield; frozen", "ND",
     "?", "10 paired buffer-comparison samples", "mixed (paired protocols)", "unknown", "unknown", 6, None, "leg 2 from 2022-06-13", T),
    ("HLY2202", "HLY2202-2022", "2022", "?", "0.22 um Sterivex", "2-8 L", "?", "No, enzymatic extraction",
     "Lysozyme, proteinase K/SDS; Zymo purification", "Old extraction method on newer cruise", "unknown", "enzymatic", "unknown", 5, None, "", T),
    ("A16N, leg 1", "A16N-2023", "2023", "ND", "0.22 um Sterivex", "~8 L", "Predominantly DNA/RNA Shield; -80 C", "Yes, usual protocol",
     "?", "Paired old/new protocols; underway pump changed", "DNA/RNA Shield", "beads", "unknown", 8, -80, "leg 1 before 2023-04-10", T),
    ("A16N, leg 2", "A16N-2023", "2023", "ND", "0.22 um Sterivex", "?", "DNA/RNA Shield or lysis buffer by protocol", "Protocol-dependent",
     "?", "Paired old/new protocols", "mixed (paired protocols)", "mixed", "unknown", None, -80, "leg 2 from 2023-04-10", T),
    ("I05", "I05-2023", "2023", "ND", "0.22 um Sterivex", "~8 L", "1 mL DNA/RNA Shield; -80 C", "Yes",
     "?", "Newer bead/Shield workflow", "DNA/RNA Shield", "beads", "unknown", 8, -80, "", T),
    ("A13.5", "A13.5-2024", "2024", "ND", "0.22 um Sterivex", "~8 L", "1 mL DNA/RNA Shield; -80 C", "Yes",
     "?", "Newer bead/Shield workflow", "DNA/RNA Shield", "beads", "unknown", 8, -80, "", T),
    ("I08S", "I08S-2024", "2024", "ND", "Sterivex; pore size ?", "~8 L", "1 mL DNA/RNA Shield; -80 C", "Yes",
     "?", "Newer bead/Shield workflow", "DNA/RNA Shield", "beads", "unknown", 8, -80, "", T),
    ("I09N 2025", "I09N-2025", "2025", "ND", "0.22 um Sterivex", "~9 L", "1 mL DNA/RNA Shield; -80 C", "?",
     "?", "Reoccupation of 2016 line", "DNA/RNA Shield", "unknown", "unknown", 9, -80, "", T),
    ("GOMECC4", "GOMECC4", "2021", "no pre-filter", "0.22 um Sterivex-GP PES", "~2 L (Niskin)", "none; -80 C", "Yes",
     "ZymoBIOMICS 96 MagBead DNA Kit (KingFisher)", "AOML Niskin samples, A/B/C replicates", "none", "beads", "none reported", 2, -80, "",
     "FAIRe workbook"),
]


# Library preparation and sequencer (second table supplied 2026-11-19, from the published run table:
# all six original runs used the same custom Nextera-style prep, Illumina Tagment DNA Enzyme and
# Buffer kit 20034197, custom 8-bp dual indexes, modified low-input tagmentation; runs 1-2 HiSeq,
# 3-6 NovaSeq). Flowcells from the FAIRe workbook (experimentRunMetadata seq_run_id), libraries
# in the KO table. Sequencer of the 2023-2025 flowcells ("...LT4") inferred from the flowcell-ID
# format of the NovaSeq X series: NOT verified.
LIBRARY = {
    "BVAL46": ("Nextera-style tagmentation (Illumina 20034197), custom 8-bp dual indexes", "HiSeq (run 2) + NovaSeq (run 4)", "HWKN2BBXX; HMMFVDSXX"),
    "AE1319": ("Nextera-style tagmentation (Illumina 20034197), custom 8-bp dual indexes", "HiSeq (run 2) + NovaSeq (run 4)", "HWKN2BBXX; HMMFVDSXX"),
    "NH1418": ("Nextera-style tagmentation (Illumina 20034197), custom 8-bp dual indexes", "HiSeq (run 2) + NovaSeq (run 4)", "HWKN2BBXX; HMMFVDSXX"),
    "I09N": ("Nextera-style tagmentation (Illumina 20034197), custom 8-bp dual indexes", "HiSeq (runs 1-2, 26 libraries) + NovaSeq (run 3, 213)", "HKMFTBBXX; HWKN2BBXX; HHGJFDSXX"),
    "P18": ("Nextera-style tagmentation (Illumina 20034197), custom 8-bp dual indexes", "NovaSeq (run 4)", "HMMFVDSXX"),
    "I07N": ("Nextera-style tagmentation (Illumina 20034197), custom 8-bp dual indexes", "NovaSeq (run 5)", "H357JDSXY"),
    "AMT28": ("Nextera-style tagmentation (Illumina 20034197), custom 8-bp dual indexes", "NovaSeq (run 4)", "HMMFVDSXX"),
    "C13.5/A13.5": ("Nextera-style tagmentation (Illumina 20034197), custom 8-bp dual indexes", "NovaSeq (run 6)", "HKNJLDRXX"),
    "HLY2202": ("Nextera DNA Flex, low-volume, custom dual indexes", "NovaSeq 6000", "HFJLFDSX7"),
    "P02W": ("not verified", "KO-table libraries: 22YJMLLT4 (NovaSeq X?, 2025); Shield pairs: HWKVFDSX5 (2023, not in KO table)", "22YJMLLT4; HWKVFDSX5"),
    "P02E": ("not verified", "as P02W", "22YJMLLT4; HWKVFDSX5"),
    "A16N, leg 1": ("not verified", "NovaSeq X? (not verified)", "22NCJYLT4; 23LT2GLT4"),
    "A16N, leg 2": ("not verified", "NovaSeq X? (not verified)", "22NCJYLT4; 23LT2GLT4"),
    "I05": ("not verified", "NovaSeq X? (not verified)", "23LT2GLT4"),
    "A13.5": ("not verified", "NovaSeq X? (not verified)", "23LT2GLT4"),
    "I08S": ("not verified", "NovaSeq X? (not verified)", "22YJMLLT4"),
    "I09N 2025": ("not verified", "NovaSeq X? (not verified)", "22YJMLLT4"),
    "GOMECC4": ("not verified", "NovaSeq X? (not verified)", "23LT2GLT4"),
}
# A16N Protocol A vs B (A16N 2023 leg 1 cruise report, section 6.13.3, p. 57): both 0.22 um Sterivex,
# -80 C. A = duplicate underway DNA sample, no bashing beads, 1,600 uL lysis buffer, about every
# 3 days; B = standard underway sample, Zymo ZR BashingBeads, 1,000 uL DNA/RNA Shield.


def main() -> None:
    out = REPO / "transfer" / "results" / "protocol"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(ROWS, columns=COLS)
    df["library_prep"] = df["row"].map(lambda r: LIBRARY[r][0])
    df["sequencer"] = df["row"].map(lambda r: LIBRARY[r][1])
    df["flowcells"] = df["row"].map(lambda r: LIBRARY[r][2])
    # coded for the libraries in the KO table (P02: 22YJMLLT4; its NovaSeq 6000 run HWKVFDSX5 is not in the table)
    df["sequencer_coded"] = df["row"].map(lambda r: "NovaSeq X? (unverified)" if LIBRARY[r][2].startswith("2")
                                          else ("HiSeq + NovaSeq 6000" if "BBXX" in LIBRARY[r][2] else "NovaSeq 6000"))
    df.to_csv(out / "cruise_methods.csv", index=False)
    print(df[["row", "cruise", "preservation", "extraction", "prefilter", "volume_l", "sequencer_coded", "flowcells"]].to_string(index=False))


if __name__ == "__main__":
    main()
