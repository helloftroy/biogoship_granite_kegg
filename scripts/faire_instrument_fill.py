#!/usr/bin/env python3
"""Fill the FAIRe `instrument` field (experimentRunMetadata) from each library's flowcell (RESULTS.md R58).

Flowcell -> sequencer, with the evidence of R39-R41:
  ...BBXX (HKMFTBBXX, HWKN2BBXX)                          Illumina HiSeq 4000     ENA instrument_model
  HHGJFDSXX, HMMFVDSXX, H357JDSXY, HKNJLDRXX              Illumina NovaSeq 6000   ENA instrument_model
  HFJLFDSX7, HMY7WDSX7 (HLY2202)                          Illumina NovaSeq 6000   publication + flowcell format
  HWKVFDSX5 (2023 mixed run)                              Illumina NovaSeq 6000   flowcell format only (unverified)
  22NCJYLT4, 23LT2GLT4, 22YJMLLT4                         Illumina NovaSeq X      original FASTQ headers (LH...)
OBI terms (OLS4 / Ontobee, checked 2026-11-23): NovaSeq 6000 OBI_0002630, NovaSeq X OBI_0003663,
NovaSeq X Plus OBI_0003664, HiSeq 4000 OBI_0002049. "LH" instrument IDs are used by both NovaSeq X and X Plus,
so X vs X Plus is not known; the X term is used and flagged.
The workbook's current values ("INSTRUMENT NAME [OBI_3310]", "[OBI_2631]", "[OBI_2131]") point to non-sequencer
OBI terms (OBI_0003310 is an assay, OBI_0002631 an scRNA-seq assay, OBI_0002131 a methylation BeadChip).

Rows without a recognised flowcell in seq_run_id (amplicon and most MetaT runs) are left as they are.
The original workbook is not modified. Outputs in data/results/metadata_status/:
  FAIRe-NOAA_uci-ess-biogoship_instrument.xlsx   copy with the column filled; changed cells highlighted
  instrument_column.txt                          the new column, one line per data row, in sheet order (for pasting)
  instrument_changes.csv                         lib_id, seq_run_id, old value, new value, evidence

Usage:
    ~/miniforge3/bin/python3 scripts/faire_instrument_fill.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "FAIRe-NOAA_uci-ess-biogoship.xlsx"
OUT = REPO / "data" / "results" / "metadata_status"
H4000 = ("Illumina HiSeq 4000 [OBI_0002049]", "ENA instrument_model")
N6000 = ("Illumina NovaSeq 6000 [OBI_0002630]", "ENA instrument_model")
NX = ("Illumina NovaSeq X [OBI_0003663]", "original FASTQ header, instrument ID LH...; X vs X Plus not known")
FLOWCELL = {
    "HKMFTBBXX": H4000, "HWKN2BBXX": H4000,
    "HHGJFDSXX": N6000, "HMMFVDSXX": N6000, "H357JDSXY": N6000, "HKNJLDRXX": N6000,
    "HFJLFDSX7": (N6000[0], "publication (HLY2202) + flowcell format"),
    "HMY7WDSX7": (N6000[0], "publication (HLY2202) + flowcell format"),
    "HWKVFDSX5": (N6000[0], "flowcell format only, unverified"),
    "22NCJYLT4": NX, "23LT2GLT4": NX, "22YJMLLT4": NX,
}  # fmt: skip


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wb = load_workbook(SRC)
    ws = wb["experimentRunMetadata"]
    head = {c.value: c.column for c in ws[3]}
    ci, cr, cl = head["instrument"], head["seq_run_id"], head["lib_id"]
    fill = PatternFill("solid", fgColor="FFF2A8")
    lines, changes = [], []
    last = max(r for r in range(4, ws.max_row + 1) if any(ws.cell(r, c).value not in (None, "") for c in (cl, cr)))
    for r in range(4, last + 1):
        run = ws.cell(r, cr).value
        old = ws.cell(r, ci).value
        fc = str(run).split("_")[1] if run and str(run).count("_") >= 2 else None
        if fc in FLOWCELL:
            new, ev = FLOWCELL[fc]
            if new != old:
                ws.cell(r, ci).value = new
                ws.cell(r, ci).fill = fill
                changes.append({"excel_row": r, "lib_id": ws.cell(r, cl).value, "seq_run_id": run, "old": old, "new": new, "evidence": ev})
            lines.append(new)
        else:
            lines.append("" if old is None else str(old))
    wb.save(OUT / "FAIRe-NOAA_uci-ess-biogoship_instrument.xlsx")
    (OUT / "instrument_column.txt").write_text("\n".join(lines) + "\n")
    ch = pd.DataFrame(changes)
    ch.to_csv(OUT / "instrument_changes.csv", index=False)
    print(f"data rows {len(lines)}; changed {len(ch)}; left as is {len(lines) - len(ch)}")
    print(ch.groupby(["new", "evidence"]).size().to_string())
    print(ch.groupby(["old", "new"], dropna=False).size().to_string())


if __name__ == "__main__":
    main()
