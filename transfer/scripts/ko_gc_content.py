#!/usr/bin/env python3
"""GC content of the Bio-GO-SHIP reference genes per KO (and per domain), for the platform-bias check (RESULTS.md R40).

If the axis-1 loadings follow the GC content of the KOs' genes, a sequencing-chemistry GC bias
(HiSeq 4000 four-colour vs NovaSeq 6000 / NovaSeq X two-colour) is the likely mechanism of the
axis-1 platform offset (R39).

Inputs (data/reference/biogoship_2026/ on Morrill):
  biogoship_eggnogmapper.tsv          eggNOG-mapper table, #query (col 1) and KEGG_ko (col 12)
  biogoship_2026_taxonomy.tsv         transcript (col 1), domain (col 3)
  biogoship_2026_database_v2.fasta    the reference sequences (212 GB; download with
                                      cluster/download_reference.sbatch FILES=...)
A deterministic 1-in-SAMPLE subset of the genes (zlib.crc32 of the sequence id) keeps the memory
small; 1 in 10 of a few hundred million genes leaves tens of thousands per common KO.

Pass 1, eggNOG: sampled id -> its KOs (genes with KEGG_ko only).
Pass 2, taxonomy: sampled id -> domain.
Pass 3, FASTA: GC fraction and length of each sampled gene with a KO; added to every KO it carries
(1/n weight, as the counting script splits a gene's reads), overall and per domain.

Output: transfer/results/ko_gc/ko_gc.csv: ko, genes (weighted), gc_mean, gc_sd, length_mean, and
gc_mean_<domain> / genes_<domain> for Bacteria, Archaea, Eukaryota, Viruses. Also progress lines.

Usage (Morrill, transfer/cluster/ko_gc.sbatch):
    python transfer/scripts/ko_gc_content.py --sample 10
"""

from __future__ import annotations

import argparse
import sys
import time
import zlib
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
REF = REPO / "data" / "reference" / "biogoship_2026"
DOMAINS = ["Bacteria", "Archaea", "Eukaryota", "Viruses"]


def keep(seq_id: bytes, m: int) -> bool:
    return zlib.crc32(seq_id) % m == 0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eggnog", type=Path, default=REF / "biogoship_eggnogmapper.tsv")
    p.add_argument("--taxonomy", type=Path, default=REF / "biogoship_2026_taxonomy.tsv")
    p.add_argument("--fasta", type=Path, default=REF / "biogoship_2026_database_v2.fasta")
    p.add_argument("--sample", type=int, default=10, help="keep 1 gene in this many")
    p.add_argument("--out-dir", type=Path, default=REPO / "transfer" / "results" / "ko_gc")
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    def say(s):
        print(f"[{(time.time() - t0) / 60:6.1f} min] {s}", flush=True)

    # pass 1: eggNOG
    kos: dict[bytes, tuple] = {}
    n = 0
    with open(args.eggnog, "rb") as fh:
        for line in fh:
            if line.startswith(b"#"):
                continue
            n += 1
            tab = line.find(b"\t")
            sid = line[:tab]
            if not keep(sid, args.sample):
                continue
            f = line.split(b"\t", 12)
            if len(f) < 12 or f[11] in (b"-", b""):
                continue
            kos[sid] = tuple(k.replace(b"ko:", b"").decode() for k in f[11].split(b","))
            if n % 50_000_000 == 0:
                say(f"eggNOG: {n:,} lines, {len(kos):,} sampled genes with a KO")
    say(f"eggNOG done: {n:,} lines; sampled genes with a KO: {len(kos):,}")

    # pass 2: taxonomy (domain)
    dom: dict[bytes, int] = {}
    n = 0
    with open(args.taxonomy, "rb") as fh:
        fh.readline()
        for line in fh:
            n += 1
            tab = line.find(b"\t")
            sid = line[:tab]
            if sid in kos:
                f = line.split(b"\t", 3)
                d = f[2].strip().decode(errors="replace") if len(f) > 2 else ""
                dom[sid] = DOMAINS.index(d) if d in DOMAINS else -1
    say(f"taxonomy done: {n:,} lines; domain for {len(dom):,} sampled genes")

    # pass 3: FASTA
    acc = defaultdict(lambda: np.zeros(4 + 2 * len(DOMAINS)))  # w, w*gc, w*gc^2, w*len, then per domain (w, w*gc)
    n = hit = 0

    def add(sid, seq):
        s = b"".join(seq).upper()
        acgt = s.count(b"A") + s.count(b"C") + s.count(b"G") + s.count(b"T")
        if acgt == 0:
            return
        gc = (s.count(b"G") + s.count(b"C")) / acgt
        ks = kos[sid]
        w = 1.0 / len(ks)
        d = dom.get(sid, -1)
        for k in ks:
            a = acc[k]
            a[0] += w
            a[1] += w * gc
            a[2] += w * gc * gc
            a[3] += w * len(s)
            if d >= 0:
                a[4 + 2 * d] += w
                a[5 + 2 * d] += w * gc

    with open(args.fasta, "rb") as fh:
        sid, seq = None, None
        for line in fh:
            if line.startswith(b">"):
                if seq is not None:
                    add(sid, seq)
                    hit += 1
                n += 1
                cand = line[1:].split(None, 1)[0]
                if cand in kos:
                    sid, seq = cand, []
                else:
                    sid, seq = None, None
                if n % 50_000_000 == 0:
                    say(f"FASTA: {n:,} sequences, {hit:,} sampled genes with a KO measured")
            elif seq is not None:
                seq.append(line.rstrip())
        if seq is not None:
            add(sid, seq)
            hit += 1
    say(f"FASTA done: {n:,} sequences; measured {hit:,} of {len(kos):,} sampled genes with a KO")

    rows = []
    for k, a in acc.items():
        w = a[0]
        mean = a[1] / w
        r = {"ko": k, "genes": w, "gc_mean": mean, "gc_sd": np.sqrt(max(a[2] / w - mean ** 2, 0)), "length_mean": a[3] / w}
        for i, d in enumerate(DOMAINS):
            wd = a[4 + 2 * i]
            r[f"genes_{d}"] = wd
            r[f"gc_mean_{d}"] = a[5 + 2 * i] / wd if wd > 0 else np.nan
        rows.append(r)
    out = pd.DataFrame(rows).sort_values("ko")
    out.round(5).to_csv(args.out_dir / "ko_gc.csv", index=False)
    say(f"wrote {args.out_dir / 'ko_gc.csv'}: {len(out):,} KOs; median genes per KO {out['genes'].median():.0f}")


if __name__ == "__main__":
    main()
