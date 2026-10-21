#!/usr/bin/env python3
"""Natural-colour images of the tiles listed in a figure's *_tiles.csv.

Converts each tile's normalised OLCI bands back to water-leaving reflectance and
combines them with the EUMETSAT broad-band recipe used in IBM's Granite helper
(red from Oa08-Oa11, green from Oa05-Oa07, blue from Oa02-Oa04), log-stretched.
Flagged pixels (cloud, glint and so on) are shown in grey. One shared stretch is
used for all tiles, so their colours can be compared.

The tiles live on the cluster (data/tiles/); copy the ones listed in the CSV to
the same relative paths before running.

Usage:
    python figures/make_tile_rgb.py figures/output/spectral_recovery_tiles.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from biogoship import tiles  # noqa: E402

BAND = {name: i for i, name in enumerate(tiles.OLCI_BANDS)}
RECIPE = {
    "red": {"Oa08": 0.16666, "Oa09": 0.66666, "Oa10": 0.08333, "Oa11": 0.08333},
    "green": {"Oa05": 0.16666, "Oa06": 0.66666, "Oa07": 0.16666},
    "blue": {"Oa02": 0.16666, "Oa03": 0.66666, "Oa04": 0.16666},
}
FLAGGED = np.array([0.82, 0.82, 0.80])


def main() -> None:
    args = parse_args()
    plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 8})
    listing = pd.read_csv(args.tiles_csv)
    images = []
    for _, row in listing.iterrows():
        path = ROOT / row["tile_path"]
        if not path.exists():
            raise SystemExit(f"Missing {path}: copy it from the cluster (same relative path) first.")
        with np.load(path) as tile:
            data, clear = tile["data"][:16].astype(float), tile["valid_mask"].astype(bool)
        reflectance = tiles.OLCI_MEANS[:, None, None] + data * tiles.OLCI_STDS[:, None, None]
        rgb = np.stack([sum(w * reflectance[BAND[b]] for b, w in RECIPE[c].items()) for c in ("red", "green", "blue")], -1)
        images.append((row, np.log10(np.clip(rgb, 1e-4, None)), clear))
    values = np.concatenate([img[clear].ravel() for _, img, clear in images])
    low, high = np.percentile(values, [1, 99.5])

    fig, axes = plt.subplots(1, len(images), figsize=(1.8 * len(images), 2.1))
    for ax, (row, img, clear) in zip(np.atleast_1d(axes), images):
        shown = np.clip((img - low) / (high - low), 0, 1)
        shown[~clear] = FLAGGED
        ax.imshow(shown, interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(f"{row['panel']}  {row['label']}", loc="left", fontsize=8, fontweight="bold", pad=3)
    out = Path(args.tiles_csv).with_name(Path(args.tiles_csv).stem.replace("_tiles", "") + "_tile_rgb.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.03)
    print(f"Wrote {out}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tiles_csv", type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    main()
