#!/usr/bin/env python3
"""Summarise every tile's raw input channels: the "raw satellite" baseline for models.

For each tile in data/manifests/biogoship_tile_manifest.csv, writes the mean and
standard deviation of each of the 17 normalised channels (16 OLCI bands, SST) over
the tile's clear pixels, plus the clear and SST-valid fractions:

    data/embeddings/tile_raw_features.csv    one row per tile, keyed by tile_path

Small enough to copy off the cluster (a few MB). CPU only, a minute or two.

Usage:
    python scripts/pool_tile_features.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from biogoship import tiles  # noqa: E402
from biogoship.tables import MANIFEST_DIR, REPO_ROOT, from_repo, read_csv, write_csv  # noqa: E402

CHANNELS = [*tiles.OLCI_BANDS, "SST"]


def main() -> None:
    args = parse_args()
    paths = sorted({row["tile_path"] for row in read_csv(args.manifest_csv)})
    rows, missing = [], 0
    for i, path in enumerate(paths, start=1):
        tile_file = from_repo(path)
        if not tile_file.exists():
            missing += 1
            continue
        with np.load(tile_file, allow_pickle=False) as tile:
            data, clear = tile["data"].astype(np.float64), tile["valid_mask"].astype(bool)
        row = {"tile_path": path, "valid_fraction": f"{clear.mean():.4f}", "sst_valid_fraction": f"{np.isfinite(data[-1]).mean():.4f}"}
        for channel, values in zip(CHANNELS, data):
            pixels = values[clear & np.isfinite(values)]
            row[f"{channel}_mean"] = f"{pixels.mean():.6g}" if pixels.size else ""
            row[f"{channel}_std"] = f"{pixels.std():.6g}" if pixels.size else ""
        rows.append(row)
        if i % 500 == 0:
            print(f"  {i:,}/{len(paths):,}", flush=True)
    write_csv(args.out_csv, rows)
    print(f"Wrote {len(rows):,} tiles to {args.out_csv}" + (f"; {missing:,} tiles not on disk" if missing else ""))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest-csv", type=Path, default=MANIFEST_DIR / "biogoship_tile_manifest.csv")
    parser.add_argument("--out-csv", type=Path, default=REPO_ROOT / "data" / "embeddings" / "tile_raw_features.csv")
    return parser.parse_args()


if __name__ == "__main__":
    main()
