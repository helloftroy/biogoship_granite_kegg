#!/usr/bin/env python3
"""Merge per-library KO counts (one folder per library) into one wide table like ko_counts_cov50_all.tsv.

The counting pipeline (biogoship-metagenomics count_biogoship_annotations.py) writes, per library folder,
level_<n>_<rank>_ko_counts.tsv (taxon, ko, length, count) for every taxonomic level. Every counted gene has a
taxonomy row, so the sum over taxa of any one level is the library's KO count; checked against
ko_counts_cov50_all.tsv for C132020_098_meta and C132020_099_meta (identical to < 1e-10). The lowest level
present is used (level 1, domain, if there). Level 2 also gives the full count; deeper levels lose genes
without that rank (order level: about 3 % fewer counts, tested on 52 libraries), and are flagged.

Also collects level_1_domain_abundance.tsv read counts (if present) into <out stem>_domain_counts.tsv
(libraries x Archaea / Bacteria / Eukaryota / Viruses, as data/targets/taxa/domain_counts.tsv).

Output layout = ko_counts_cov50_all.tsv: header "count.ko.Kxxxxx ..." with NO row-name column; each row
"<library>\\t<counts...>". KO columns: those of --template (ko_counts_cov50_all.tsv) in its order, then any
other KOs. Files may be plain or .gz. A library found twice: the first is kept and the repeat is reported.

Roots may also be Google Drive zips (drive-download-*.zip, members <library>/level_*.tsv): with --stash DIR,
the level-1 KO and abundance files of every library are first copied gzipped to DIR/<library>/ (small; the
zips can then be deleted), and the merge reads DIR.

Usage:
    python scripts/merge_ko_counts.py <root> [<root> ...] --out ko_counts_cov0_all.tsv
    python scripts/merge_ko_counts.py drive-download-*.zip --stash data/targets/cov0_stash --out ko_counts_cov0_all.tsv
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
KO_FILE = re.compile(r"^level_(\d)_[a-z_]+_ko_counts\.tsv(\.gz)?$")
DOMAINS = ["Archaea", "Bacteria", "Eukaryota", "Viruses"]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip" if path.suffix == ".gz" else None)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("roots", type=Path, nargs="+", help="folders holding <library>/level_*_ko_counts.tsv (searched recursively)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--template", type=Path, default=REPO / "ko_counts_cov50_all.tsv")
    p.add_argument("--stash", type=Path, help="extract level-1 files from zip roots into this folder first")
    args = p.parse_args()
    zips = [r for r in args.roots if r.suffix == ".zip"]
    if zips:
        if args.stash is None:
            raise SystemExit("zip roots need --stash")
        import zipfile

        n = 0
        for z in zips:
            with zipfile.ZipFile(z) as zf:
                for name in zf.namelist():
                    m = re.search(r"(?:^|/)([^/]+)/(level_1_[a-z_]+_(ko_counts|abundance))\.tsv$", name)
                    if not m:
                        continue
                    dest = args.stash / m.group(1) / (m.group(2) + ".tsv.gz")
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    with zf.open(name) as src, gzip.open(dest, "wb") as out:
                        out.write(src.read())
                    n += 1
        print(f"stashed {n:,} level-1 files from {len(zips)} zips into {args.stash}")
        args.roots = [r for r in args.roots if r.suffix != ".zip"] + [args.stash]

    libs: dict[str, Path] = {}
    for root in args.roots:
        for f in sorted(root.rglob("level_*_ko_counts.tsv*")):
            m = KO_FILE.match(f.name)
            if not m:
                continue
            lib, lvl = f.parent.name, int(m.group(1))
            cur = libs.get(lib)
            if cur is None:
                libs[lib] = f
            elif cur.parent != f.parent:
                if int(KO_FILE.match(cur.name).group(1)) == lvl:
                    print(f"duplicate library {lib}: keeping {cur.parent}, ignoring {f.parent}")
            elif lvl < int(KO_FILE.match(cur.name).group(1)):
                libs[lib] = f
    print(f"libraries found: {len(libs):,}")
    deep = sorted(lib for lib, f in libs.items() if int(KO_FILE.match(f.name).group(1)) > 2)
    if deep:
        print(f"WARNING: {len(deep)} libraries have no level 1/2 KO file; a deeper level drops genes without that rank "
              f"(about 3 % of counts at order level): {', '.join(deep[:10])}{' ...' if len(deep) > 10 else ''}")
    if not libs:
        raise SystemExit("no level_*_ko_counts.tsv files found under " + ", ".join(map(str, args.roots)))

    template = [c.replace("count.ko.", "") for c in open(args.template).readline().rstrip("\n").split("\t")] if args.template.exists() else []
    series, doms = {}, {}
    for i, (lib, f) in enumerate(sorted(libs.items())):
        k = read_tsv(f)
        s = k.groupby("ko")["count"].sum()
        s.index = s.index.str.replace("ko:", "", regex=False)
        series[lib] = s
        a = next(iter(sorted(f.parent.glob("level_1_*_abundance.tsv*"))), None)
        if a is not None:
            d = read_tsv(a)
            doms[lib] = d.groupby("taxon")["read_count"].sum().reindex(DOMAINS).fillna(0)
        if (i + 1) % 250 == 0:
            print(f"  {i + 1:,}/{len(libs):,}", flush=True)
    allk = set().union(*[s.index for s in series.values()])
    cols = template + sorted(allk - set(template))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join("count.ko." + k for k in cols) + "\n")
        for lib, s in series.items():
            v = s.reindex(cols).fillna(0).to_numpy()
            fh.write(lib + "\t" + "\t".join(f"{x:.10g}" for x in v) + "\n")
    tot = np.array([s.sum() for s in series.values()])
    print(f"wrote {args.out}: {len(series):,} libraries x {len(cols):,} KOs; KO reads per library median {np.median(tot):,.0f}")
    if doms:
        dp = args.out.with_name(args.out.name.replace(".tsv", "") + "_domain_counts.tsv")
        pd.DataFrame(doms).T.to_csv(dp, sep="\t")
        print(f"wrote {dp}: {len(doms):,} libraries")


if __name__ == "__main__":
    main()
