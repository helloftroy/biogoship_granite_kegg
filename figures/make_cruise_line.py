#!/usr/bin/env python3
"""Figure: three KOs along one cruise line, observed and predicted with the cruise held out.

a  Track of the cruise (sampling positions).
b  In-situ temperature along the line.
c-e  One KO each: observed CLR (black dots), held-out prediction from the Granite
     embedding (ridge; blue, solid) and from in-situ temperature (spline + ridge;
     magenta, dashed). R2 for this cruise at the right (the per-cruise R2 of
     scripts/test_ko_temperature.py: against the mean of the training cruises).

The three KOs show the three cases of the heatmap figure: both about equal and
high, Granite clearly better, temperature better. Chosen on P18 (2016-17, the
eastern Pacific from 33 N to 70 S), where each category has a KO that varies
strongly along the line.

Input:   data/results/ko_temperature_all_ostia/ (scripts/test_ko_temperature.py --kos all
         --learners ridge spline --feature-sets temp embedding --embeddings-dirs
         data/embeddings_sst_l4 data/faire/embeddings_sst_l4 --save-predictions)
         Land outline: Natural Earth 110 m land (cartopy's cache), drawn if present.
Outputs: figures/output/cruise_line.png and .pdf

Usage:
    python figures/make_cruise_line.py
    python figures/make_cruise_line.py --cruise A16N-2023 --kos K05382 K04518 K02327
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results" / "ko_temperature_all_ostia"
OUT = ROOT / "figures" / "output"
LAND = Path.home() / ".local/share/cartopy/shapefiles/natural_earth/physical/ne_110m_land.shp"

INK, MUTED, RULE, LAND_FILL = "#0b0b0b", "#52514e", "#d6d5d0", "#e4e3de"
GRANITE, TEMP = "#2a78d6", "#c0457e"
LABELS = {  # short labels: gene, what it is, and who carries it
    "K05382": ("cpeS", "phycoerythrin linker protein", "Synechococcus light harvesting"),
    "K01425": ("glsA", "glutaminase", "heterotrophic bacteria, N metabolism"),
    "K16090": ("fiu", "catecholate siderophore receptor", "iron uptake"),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cruise", default="P18")
    parser.add_argument("--kos", nargs=3, default=["K05382", "K01425", "K16090"])
    parser.add_argument("--out-name", default="cruise_line")
    args = parser.parse_args()
    plt.rcParams.update({
        "font.family": ["Arial", "DejaVu Sans"], "font.size": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": RULE, "axes.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": MUTED,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.bf": "Arial:bold",
    })  # fmt: skip

    rows = pd.read_csv(RESULTS / "prediction_rows.csv")
    granite = np.load(RESULTS / "predictions_embedding_ridge.npz", allow_pickle=True)
    temp = np.load(RESULTS / "predictions_temp_spline.npz", allow_pickle=True)
    kos = list(granite["ko"])
    per_cruise = pd.read_csv(RESULTS / "ko_temperature_per_cruise.csv")
    per_cruise = per_cruise[per_cruise["held_out"] == args.cruise].set_index(["ko", "feature_set", "learner"])["r2"]
    names = pd.read_csv(ROOT / "data/reference/kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "name"]).set_index("ko")["name"]

    on = (rows["section"] == args.cruise).to_numpy()
    line = rows[on].assign(i=np.flatnonzero(on)).sort_values("sample_lat")
    lat = line["sample_lat"].to_numpy()
    gap = np.r_[False, np.diff(lat) > 8]  # do not draw lines across stretches with no samples
    def broken(values):
        return np.insert(np.asarray(values, dtype=float), np.flatnonzero(gap), np.nan)
    lat_line = broken(lat)

    fig = plt.figure(figsize=(7.5, 6.6))
    left, width = 0.30, 0.56
    heights = [0.13, 0.22, 0.22, 0.22]
    tops = [0.93]
    for h in heights[:-1]:
        tops.append(tops[-1] - h - 0.025)
    axes = [fig.add_axes([left, top - h, width, h]) for top, h in zip(tops, heights)]
    panel_map(fig.add_axes([0.06, 0.08, 0.15, 0.85]), rows, line)

    ax = axes[0]
    ax.plot(lat_line, broken(line["temperature_c"]), color=INK, lw=1.0)
    ax.scatter(lat, line["temperature_c"], s=8, color=INK, zorder=3, linewidths=0)
    ax.set_ylabel("In-situ\ntemperature (°C)")
    tag(ax, "b")

    for k, (ax, ko) in enumerate(zip(axes[1:], args.kos)):
        j = kos.index(ko)
        observed = granite["observed"][line["i"], j]
        ax.scatter(lat, observed, s=9, color=INK, zorder=4, linewidths=0, label="Observed")
        ax.plot(lat_line, broken(granite["predicted"][line["i"], j]), color=GRANITE, lw=1.6, marker="o", ms=2.2, zorder=3, label="Granite embedding")
        ax.plot(lat_line, broken(temp["predicted"][line["i"], j]), color=TEMP, lw=1.6, ls=(0, (4, 2)), marker="o", ms=2.2, zorder=2, label="In-situ temperature")
        gene, what, who = LABELS.get(ko, (names.get(ko, ko).split(";")[0], names.get(ko, "").split(";")[-1].strip(), ""))
        ax.text(0.01, 0.97, rf"$\bf{{{gene}}}$   {ko} · {what}" + (f" ({who})" if who else ""),
                transform=ax.transAxes, ha="left", va="top", fontsize=7.5, color=INK)  # fmt: skip
        r2g = per_cruise.get((ko, "embedding", "ridge"), np.nan)
        r2t = per_cruise.get((ko, "temp", "spline"), np.nan)
        x0 = 1.015
        ax.plot([x0, x0 + 0.035], [0.62, 0.62], transform=ax.transAxes, color=GRANITE, lw=1.6, clip_on=False)
        ax.text(x0 + 0.045, 0.62, f"R² {r2g:.2f}", transform=ax.transAxes, va="center", fontsize=7.5, color=INK)
        ax.plot([x0, x0 + 0.035], [0.42, 0.42], transform=ax.transAxes, color=TEMP, lw=1.6, ls=(0, (4, 2)), clip_on=False)
        ax.text(x0 + 0.045, 0.42, f"R² {r2t:.2f}", transform=ax.transAxes, va="center", fontsize=7.5, color=INK)
        ax.set_ylabel("KO abundance\n(CLR)")
        lo, hi = np.nanmin(observed), np.nanmax(observed)
        ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.30 * (hi - lo))  # room for the label
        tag(ax, "cde"[k])
        if k == 0:
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(left - 0.01, 0.995), ncol=3, frameon=False,
                       fontsize=7.5, handlelength=2.6, columnspacing=1.6, labelcolor=INK)  # fmt: skip

    for ax in axes:
        ax.set_xlim(lat.min() - 2, lat.max() + 2)
        ax.grid(axis="y", color=RULE, lw=0.4)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(length=2)
        if ax is not axes[-1]:
            ax.tick_params(labelbottom=False)
    axes[-1].set_xticks(np.arange(-60, 40, 20))
    axes[-1].set_xticklabels([f"{abs(v)}°{'S' if v < 0 else 'N' if v > 0 else ''}" for v in range(-60, 40, 20)])
    axes[-1].set_xlabel("Latitude")

    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"{args.out_name}.{suffix}", dpi=300)
    print(f"Wrote {OUT / args.out_name}.png and .pdf ({len(line)} tiles on {args.cruise})")


def panel_map(ax, rows: pd.DataFrame, line: pd.DataFrame) -> None:
    lon_pad, lat_pad = 22, 6
    lon0, lon1 = line["sample_lon"].min() - lon_pad, line["sample_lon"].max() + lon_pad
    lat0, lat1 = line["sample_lat"].min() - lat_pad, line["sample_lat"].max() + lat_pad
    for ring in land_rings():
        ax.fill(ring[:, 0], ring[:, 1], color=LAND_FILL, lw=0)
    others = rows[rows["section"] != line["section"].iloc[0]]
    ax.scatter(others["sample_lon"], others["sample_lat"], s=2, color=RULE, linewidths=0)
    ax.plot(line["sample_lon"], line["sample_lat"], color=INK, lw=0.6)
    ax.scatter(line["sample_lon"], line["sample_lat"], s=9, color=INK, linewidths=0, zorder=3)
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(lat0, lat1)
    ax.set_aspect("equal")
    ax.set_xticks([round(lon0 / 10) * 10 + 5, round(lon1 / 10) * 10 - 5])
    ax.set_xticklabels([f"{abs(v)}°W" for v in ax.get_xticks()])
    ax.set_yticks(np.arange(-60, 40, 20))
    ax.set_yticklabels([f"{abs(v)}°{'S' if v < 0 else 'N' if v > 0 else ''}" for v in range(-60, 40, 20)])
    ax.tick_params(length=2)
    ax.text(0.5, 1.01, f"{line['section'].iloc[0]}, {line['sample_time'].min()[:7]} to {line['sample_time'].max()[:7]}",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=7, color=MUTED)  # fmt: skip
    tag(ax, "a", x=-0.40)


def land_rings() -> list[np.ndarray]:
    """Polygon rings of a Natural Earth land shapefile (shape type 5), without extra packages."""
    if not LAND.exists():
        return []
    data, rings, pos = LAND.read_bytes(), [], 100
    while pos < len(data):
        _, length = struct.unpack(">ii", data[pos:pos + 8])
        rec = data[pos + 8:pos + 8 + 2 * length]
        pos += 8 + 2 * length
        if struct.unpack("<i", rec[:4])[0] != 5:
            continue
        n_parts, n_points = struct.unpack("<ii", rec[36:44])
        parts = list(struct.unpack(f"<{n_parts}i", rec[44:44 + 4 * n_parts])) + [n_points]
        points = np.frombuffer(rec[44 + 4 * n_parts:44 + 4 * n_parts + 16 * n_points], dtype="<f8").reshape(-1, 2)
        rings += [points[a:b] for a, b in zip(parts[:-1], parts[1:])]
    return rings


def tag(ax, letter: str, x: float = -0.13) -> None:
    ax.text(x, 1.0, letter, transform=ax.transAxes, ha="left", va="top", fontsize=10, fontweight="bold", color=INK)


if __name__ == "__main__":
    main()
