#!/usr/bin/env python3
"""Figure: three KOs along one cruise drawn on the map as coloured tracks.

One map per KO. Each map has two parallel copies of the cruise track, side by side
(in reality they lie on top of each other): observed KO abundance (CLR) and the
held-out prediction of a Granite model, on one colour scale per KO. R2 on this cruise
is given under the map (against the mean of the training cruises, as elsewhere), with
temperature alone for reference.

The Granite model can differ per KO (--models, one per KO):
  pca32 + temp~   32 PCs of the Granite embedding + temperature spline, ridge
                  (RESULTS.md R18; scripts/tune_figure_kos.py, all tiles)
  emb             the Granite embedding alone, ridge (same run)
  P~ > pca32 ...  in-situ phosphate first, then Granite PCs on the residuals (R20;
                  scripts/phosphate_models.py, tiles with phosphate only)
Defaults: the model with the highest pooled R2 per KO in R20 (pstC: phosphate first).

Input:   data/results/figure_ko_models/, data/results/phosphate_models/
         data/results/ko_temperature_all_ostia_v2/prediction_rows.csv (positions)
Outputs: figures/output/cruise_map.png and .pdf

Usage:
    python figures/make_cruise_map.py
    python figures/make_cruise_map.py --kos K05382 K02037 K14835 --models "pca32 + temp~" "P~ > pca32" "pca32 + temp~"
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
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_cruise_line import LAND_FILL, MUTED, INK, RULE, land_rings  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TUNED = ROOT / "data" / "results" / "figure_ko_models"
PHOSPHATE = ROOT / "data" / "results" / "phosphate_models"
ROWS = ROOT / "data" / "results" / "ko_temperature_all_ostia_v2" / "prediction_rows.csv"
OUT = ROOT / "figures" / "output"
CMAP = plt.get_cmap("viridis")
LABELS = {
    "K05382": ("cpeS", "phycoerythrin linker", "Synechococcus"),
    "K02037": ("pstC", "phosphate transporter", "heterotrophic bacteria"),
    "K14835": ("NOP2", "rRNA methyltransferase", "eukaryotes"),
    "K02036": ("pstB", "phosphate transporter", "heterotrophic bacteria"),
    "K18292": ("E4.1.3.25", "citramalyl-CoA lyase", "heterotrophic bacteria"),
    "K13991": ("puhA", "reaction centre H", "purple bacteria"),
}
OFFSET = 5.0  # degrees of longitude from the cruise line to each track
TRACKS = [("Observed", -OFFSET), ("Granite\nmodel", OFFSET)]
MODEL_TEXT = {
    "pca32 + temp~": "32 Granite PCs + temperature",
    "emb": "Granite embedding",
    "P~ > pca32": "phosphate, then 32 Granite PCs",
    "T~P~ > pca32": "temperature + phosphate, then 32 Granite PCs",
    "P~ > pca32 + T~": "phosphate, then 32 Granite PCs + temperature",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cruise", default="P18")
    parser.add_argument("--kos", nargs=3, default=["K05382", "K02037", "K14835"])
    parser.add_argument("--models", nargs=3, default=["pca32 + temp~", "P~ > pca32", "pca32 + temp~"],
                        help="Granite model per KO (see the module docstring).")
    args = parser.parse_args()
    plt.rcParams.update({
        "font.family": ["Arial", "DejaVu Sans"], "font.size": 8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "axes.edgecolor": RULE, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.bf": "Arial:bold",
    })  # fmt: skip

    rows = pd.read_csv(ROWS)
    runs = {}
    for folder in (TUNED, PHOSPHATE):
        z = np.load(folder / "predictions.npz", allow_pickle=True)
        assert (z["tile_path"] == rows["tile_path"].to_numpy().astype(str)).all()
        r2 = pd.read_csv(folder / "model_r2_cruise.csv")
        runs[folder] = (z, r2[r2["held_out"] == args.cruise].set_index(["ko", "model"])["r2"])

    on = (rows["section"] == args.cruise).to_numpy()
    line = rows[on].assign(i=np.flatnonzero(on)).sort_values("sample_lat")
    lat, lon = line["sample_lat"].to_numpy(), line["sample_lon"].to_numpy()
    joined = np.abs(np.diff(lat)) <= 8  # segments only between neighbouring samples

    rings = land_rings()
    lon0, lon1 = lon.min() - OFFSET - 7, lon.max() + OFFSET + 7
    lat0, lat1 = lat.min() - 4, lat.max() + 4
    map_w = 2.25  # inches
    map_h = map_w * (lat1 - lat0) / (lon1 - lon0)
    bottom, top = 1.35, 0.45  # inches below (track labels, colour key, R2) and above (title) the maps
    fig_h = bottom + map_h + top
    fig = plt.figure(figsize=(7.5, fig_h))
    for k, (ko, model) in enumerate(zip(args.kos, args.models)):
        z, r2 = runs[PHOSPHATE if "P~" in model else TUNED]
        j = list(z["ko"]).index(ko)
        values = {
            "Observed": z["observed"][line["i"], j],
            "Granite\nmodel": z["predicted"][list(z["models"]).index(model)][line["i"], j],
        }
        lo, hi = np.percentile(values["Observed"], [2, 98])
        norm = Normalize(lo, hi)
        x0 = 0.045 + k * 0.325
        ax = fig.add_axes([x0, bottom / fig_h, map_w / 7.5, map_h / fig_h])
        for ring in rings:
            ax.fill(ring[:, 0], ring[:, 1], color=LAND_FILL, lw=0)
        for name, shift in TRACKS:
            v = values[name]
            ok = ~np.isnan(v)  # the phosphate model has no prediction where phosphate is missing
            pts = np.column_stack([lon + shift, lat])
            link = joined & ok[:-1] & ok[1:]
            segs = np.stack([pts[:-1], pts[1:]], axis=1)[link]
            colour = ((v[:-1] + v[1:]) / 2)[link]
            ax.add_collection(LineCollection(segs, colors=CMAP(norm(colour)), linewidths=5.0, capstyle="round"))
            alone = ok & ~(np.r_[False, link] | np.r_[link, False])  # samples with no neighbour within 8 degrees
            ax.scatter(pts[alone, 0], pts[alone, 1], s=26, c=CMAP(norm(v[alone])), linewidths=0, zorder=3)
            south = np.argmin(lat)
            ax.text(lon[south] + shift, lat0 - 1.5, name, ha="center", va="top", fontsize=6.3, color=INK, clip_on=False)
        ax.set_xlim(lon0, lon1)
        ax.set_ylim(lat0, lat1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks(np.arange(-60, 40, 20))
        ax.set_yticklabels([f"{abs(v)}°{'S' if v < 0 else 'N' if v > 0 else ''}" for v in range(-60, 40, 20)]
                           if k == 0 else [])  # fmt: skip
        ax.tick_params(length=2)
        gene, what, who = LABELS.get(ko, (ko, "", ""))
        ax.set_title(rf"$\bf{{{gene}}}$ · {what}" + (f"\n{who}" if who else ""), fontsize=7.5, color=INK, pad=4,
                     linespacing=1.15)  # fmt: skip

        cax = fig.add_axes([x0 + 0.04, 0.98 / fig_h, map_w / 7.5 - 0.08, 0.07 / fig_h])
        cax.imshow(np.linspace(lo, hi, 200)[None, :], aspect="auto", cmap=CMAP, extent=[lo, hi, 0, 1])
        cax.set_yticks([])
        cax.tick_params(length=2, labelsize=6.3, pad=1.5)
        for side in cax.spines.values():
            side.set_visible(False)
        cax.set_xlabel("KO abundance (CLR)", fontsize=6.5, color=MUTED, labelpad=1.5)

        r2_model = r2.get((ko, model), np.nan)
        r2_temp = runs[TUNED][1].get((ko, "temp~"), np.nan)
        fig.text(x0 + map_w / 15, 0.36 / fig_h, f"Granite model R² {r2_model:.2f}", ha="center", va="center",
                 fontsize=8, color=INK)  # fmt: skip
        fig.text(x0 + map_w / 15, 0.16 / fig_h, f"{MODEL_TEXT.get(model, model)}\n(temperature alone: R² {r2_temp:.2f})",
                 ha="center", va="center", fontsize=6.3, color=MUTED, linespacing=1.3)  # fmt: skip

    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"cruise_map.{suffix}", dpi=300)
    print(f"Wrote {OUT / 'cruise_map'}.png and .pdf ({len(line)} tiles on {args.cruise})")


if __name__ == "__main__":
    main()
