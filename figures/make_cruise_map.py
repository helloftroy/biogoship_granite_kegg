#!/usr/bin/env python3
"""Figure: three KOs along one cruise drawn on the map as coloured tracks.

One map per KO. Each map has three parallel copies of the cruise track, side by side
(in reality they lie on top of each other): observed KO abundance (CLR), the held-out
prediction of the best Granite model for that KO, and the held-out prediction from
in-situ temperature (spline). All three share one colour scale per KO. R2 on this
cruise for the two predictions is given under the map (against the mean of the
training cruises, as elsewhere).

Granite model = 32 principal components of the Granite embedding + a spline of in-situ
temperature, ridge (the leading single model of RESULTS.md R18, --model); temperature =
temp~ of the same run (scripts/tune_figure_kos.py). The per-KO choice ("best per KO") is
not used by default: for some KOs it picks temperature alone.

Input:   data/results/figure_ko_models/ (scripts/tune_figure_kos.py)
         data/results/ko_temperature_all_ostia_v2/prediction_rows.csv (positions)
Outputs: figures/output/cruise_map.png and .pdf

Usage:
    python figures/make_cruise_map.py
    python figures/make_cruise_map.py --kos K05382 K02037 K14835 --cruise P18
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
from make_cruise_line import LAND_FILL, MUTED, INK, RULE, land_rings, readable  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TUNED = ROOT / "data" / "results" / "figure_ko_models"
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
OFFSET = 9.0  # degrees of longitude between the three tracks
TRACKS = [("Observed", -OFFSET), ("Granite", 0.0), ("Temp.", OFFSET)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cruise", default="P18")
    parser.add_argument("--kos", nargs=3, default=["K05382", "K02037", "K14835"])
    parser.add_argument("--model", default="pca32 + temp~", help="Model of tune_figure_kos.py drawn as the Granite track.")
    args = parser.parse_args()
    plt.rcParams.update({
        "font.family": ["Arial", "DejaVu Sans"], "font.size": 8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "axes.edgecolor": RULE, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.bf": "Arial:bold",
    })  # fmt: skip

    rows = pd.read_csv(ROWS)
    tuned = np.load(TUNED / "predictions.npz", allow_pickle=True)
    assert (tuned["tile_path"] == rows["tile_path"].to_numpy().astype(str)).all()
    models, kos = list(tuned["models"]), list(tuned["ko"])
    r2 = pd.read_csv(TUNED / "model_r2_cruise.csv")
    r2 = r2[r2["held_out"] == args.cruise].set_index(["ko", "model"])["r2"]
    picks = pd.read_csv(TUNED / "best_per_ko.csv")
    picks = picks[picks["held_out"] == args.cruise].set_index("ko")["model"]

    on = (rows["section"] == args.cruise).to_numpy()
    line = rows[on].assign(i=np.flatnonzero(on)).sort_values("sample_lat")
    lat, lon = line["sample_lat"].to_numpy(), line["sample_lon"].to_numpy()
    joined = np.abs(np.diff(lat)) <= 8  # segments only between neighbouring samples

    rings = land_rings()
    lon0, lon1 = lon.min() - OFFSET - 7, lon.max() + OFFSET + 7
    lat0, lat1 = lat.min() - 4, lat.max() + 4
    map_w = 2.25  # inches
    map_h = map_w * (lat1 - lat0) / (lon1 - lon0)
    bottom, top = 1.25, 0.45  # inches below (track labels, colour key, R2) and above (title) the maps
    fig_h = bottom + map_h + top
    fig = plt.figure(figsize=(7.5, fig_h))
    for k, ko in enumerate(args.kos):
        j = kos.index(ko)
        values = {
            "Observed": tuned["observed"][line["i"], j],
            "Granite": tuned["predicted"][models.index(args.model)][line["i"], j],
            "Temp.": tuned["predicted"][models.index("temp~")][line["i"], j],
        }
        lo, hi = np.percentile(values["Observed"], [2, 98])
        norm = Normalize(lo, hi)
        x0 = 0.045 + k * 0.325
        ax = fig.add_axes([x0, bottom / fig_h, map_w / 7.5, map_h / fig_h])
        for ring in rings:
            ax.fill(ring[:, 0], ring[:, 1], color=LAND_FILL, lw=0)
        for name, shift in TRACKS:
            v = values[name]
            pts = np.column_stack([lon + shift, lat])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)[joined]
            colour = ((v[:-1] + v[1:]) / 2)[joined]
            ax.add_collection(LineCollection(segs, colors=CMAP(norm(colour)), linewidths=4.2, capstyle="round"))
            alone = ~(np.r_[False, joined] | np.r_[joined, False])  # samples with no neighbour within 8 degrees
            ax.scatter(pts[alone, 0], pts[alone, 1], s=22, c=CMAP(norm(v[alone])), linewidths=0, zorder=3)
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

        cax = fig.add_axes([x0 + 0.04, 0.92 / fig_h, map_w / 7.5 - 0.08, 0.07 / fig_h])
        cax.imshow(np.linspace(lo, hi, 200)[None, :], aspect="auto", cmap=CMAP, extent=[lo, hi, 0, 1])
        cax.set_yticks([])
        cax.tick_params(length=2, labelsize=6.3, pad=1.5)
        for side in cax.spines.values():
            side.set_visible(False)
        cax.set_xlabel("KO abundance (CLR)", fontsize=6.5, color=MUTED, labelpad=1.5)

        r2_model, r2_temp = r2.get((ko, args.model), np.nan), r2.get((ko, "temp~"), np.nan)
        fig.text(x0 + map_w / 15, 0.25 / fig_h, f"Granite model R² {r2_model:.2f}\nTemperature R² {r2_temp:.2f}",
                 ha="center", va="center", fontsize=7.5, color=INK, linespacing=1.5)  # fmt: skip
        if args.model == "best per KO" and ko in picks:
            fig.text(x0 + map_w / 15, 0.06 / fig_h, f"({readable(picks[ko])})", ha="center", va="center",
                     fontsize=6.0, color=MUTED)  # fmt: skip

    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"cruise_map.{suffix}", dpi=300)
    print(f"Wrote {OUT / 'cruise_map'}.png and .pdf ({len(line)} tiles on {args.cruise})")


if __name__ == "__main__":
    main()
