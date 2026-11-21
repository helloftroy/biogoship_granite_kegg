#!/usr/bin/env python3
"""Merge per-library KO counts (one folder per library) into one wide table like ko_counts_cov50_all.tsv.

The counting pipeline (biogoship-metagenomics count_biogoship_annotations.py) writes, per library folder,
level_<n>_<rank>_ko_counts.tsv (taxon, ko, length, count) for every taxonomic level. The official
ko_counts_cov50_all.tsv is the sum over taxa of the level-2 (phylum) file: identical for all 3,081 libraries
(checked 2026-11-20). The level-1 (domain) file holds MORE counts (median +8 %, up to +44 %): an
"Unclassified" taxon plus genes with a domain but no phylum. Deeper levels hold fewer (order: about -3 %).
So --level (default 2) picks the level; a library without that level falls back to the nearest one present
and is flagged. Compare tables only when built from the same level.

Also collects level_1_domain_abundance.tsv read counts (if present) into <out stem>_domain_counts.tsv
(libraries x Archaea / Bacteria / Eukaryota / Viruses, as data/targets/taxa/domain_counts.tsv).

Output layout = ko_counts_cov50_all.tsv: header "count.ko.Kxxxxx ..." with NO row-name column; each row
"<library>\\t<counts...>". KO columns: those of --template (ko_counts_cov50_all.tsv) in its order, then any
other KOs. Files may be plain or .gz. A library found twice: the first is kept and the repeat is reported.

Roots may also be Google Drive zips (drive-download-*.zip, members <library>/level_*.tsv): with --stash DIR,
the KO and abundance files of the --stash-levels levels (default 1) of every library are first copied gzipped
to DIR/<library>/ (small; the zips can then be deleted; files already there are kept), and the merge reads DIR.
Zips are extracted in parallel (--jobs).

Usage:
    python scripts/merge_ko_counts.py <root> [<root> ...] --out ko_counts_cov0_all.tsv
    python scripts/merge_ko_counts.py drive-download-*.zip --stash data/targets/cov0_stash --level 1 --out ko_counts_cov0_all.tsv
    python scripts/merge_ko_counts.py coverage_50-*.zip --stash data/targets/taxa/ko_by_taxon --stash-levels 1 2 \
        --out ko_counts_cov50_rebuilt.tsv
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


def stash_zip(z: Path, stash: Path, levels: tuple[int, ...]) -> tuple[int, int]:
    import zipfile

    pat = re.compile(r"(?:^|/)([^/]+)/(level_([%s])_[a-z_]+_(ko_counts|abundance))\.tsv$" % "".join(map(str, levels)))
    n = kept = 0
    with zipfile.ZipFile(z) as zf:
        for name in zf.namelist():
            m = pat.search(name)
            if not m:
                continue
            dest = stash / m.group(1) / (m.group(2) + ".tsv.gz")
            if dest.exists():
                kept += 1
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_suffix(".tmp")
            with zf.open(name) as src, gzip.open(tmp, "wb", compresslevel=6) as out:
                out.write(src.read())
            tmp.rename(dest)
            n += 1
    return n, kept


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("roots", type=Path, nargs="+", help="folders holding <library>/level_*_ko_counts.tsv (searched recursively)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--template", type=Path, default=REPO / "ko_counts_cov50_all.tsv")
    p.add_argument("--stash", type=Path, help="extract KO/abundance files from zip roots into this folder first")
    p.add_argument("--stash-levels", type=int, nargs="+", default=[1], help="taxonomic levels to extract (default 1)")
    p.add_argument("--level", type=int, default=2, help="taxonomic level summed over taxa (default 2 = official table)")
    p.add_argument("--jobs", type=int, default=6, help="zips extracted in parallel")
    args = p.parse_args()
    zips = [r for r in args.roots if r.suffix == ".zip"]
    if zips:
        if args.stash is None:
            raise SystemExit("zip roots need --stash")
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(args.jobs) as ex:
            got = list(ex.map(stash_zip, zips, [args.stash] * len(zips), [tuple(args.stash_levels)] * len(zips)))
        print(f"stashed {sum(g[0] for g in got):,} files ({sum(g[1] for g in got):,} already present, kept) "
              f"from {len(zips)} zips into {args.stash}")
        args.roots = [r for r in args.roots if r.suffix != ".zip"] + [args.stash]

    found: dict[str, dict[int, Path]] = {}
    for root in args.roots:
        for f in sorted(root.rglob("level_*_ko_counts.tsv*")):
            m = KO_FILE.match(f.name)
            if not m:
                continue
            lib, lvl = f.parent.name, int(m.group(1))
            have = found.setdefault(lib, {})
            if lvl in have and have[lvl].parent != f.parent:
                print(f"duplicate library {lib}: keeping {have[lvl].parent}, ignoring {f.parent}")
                continue
            have.setdefault(lvl, f)
    libs = {lib: have[min(have, key=lambda v: (abs(v - args.level), v))] for lib, have in found.items()}
    print(f"libraries found: {len(libs):,}")
    off = sorted(lib for lib, f in libs.items() if int(KO_FILE.match(f.name).group(1)) != args.level)
    if off:
        print(f"WARNING: {len(off)} libraries have no level-{args.level} KO file and use another level (level 1 holds "
              f"about 8 % more counts, order level about 3 % fewer): {', '.join(off[:10])}{' ...' if len(off) > 10 else ''}")
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
