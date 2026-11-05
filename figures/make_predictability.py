#!/usr/bin/env python3
"""Figure: what Granite embeddings predict, one cruise held out at a time.

a  In-situ temperature, phosphate and nitrate + nitrite (original 5 cruises).
   R2 per held-out cruise (small dots) and mean +/- s.d. over cruises (large
   marker, line) for the Granite embedding and the raw band averages of the
   fixed (filled) and advected (hollow) tile, and for position + season.
   Tile rule: closest-in-time tile with SST; tiles without SST are not used.
b  The 20 KOs best predicted by the Granite embedding (13 cruises), with their
   KEGG gene names: embedding R2 against position + season R2 (pooled over all
   held-out cruises).
c  Pooled R2 of all screened KOs: Granite embedding, position + season, and the
   embedding shuffled between samples (chance).

Inputs:  data/results/sst_only_prefer/feature_set_comparison_section_unique_tiles.csv
         (compare_feature_sets.py --require-sst --best-tile prefer-sst)
         data/results/ko_screen/ko_screen.csv (screen_ko_predictability.py)
         data/reference/kegg_ko_list.tsv (KEGG KO names; downloaded if missing)
Outputs: figures/output/predictability.png and .pdf

Usage:
    python figures/make_predictability.py
"""

from __future__ import annotations

import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INSITU = ROOT / "data" / "results" / "sst_only_prefer" / "feature_set_comparison_section_unique_tiles.csv"
KO_SCREEN = ROOT / "data" / "results" / "ko_screen" / "ko_screen.csv"
KO_NAMES = ROOT / "data" / "reference" / "kegg_ko_list.tsv"
OUT = ROOT / "figures" / "output"

INK, MUTED, RULE = "#0b0b0b", "#52514e", "#d6d5d0"
GRANITE, RAW, META, NULL = "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
N_TOP_KOS = 20
R2_FLOOR = -1.0  # panel a: per-cruise values below this are drawn at the floor as a left-pointing marker

# Panel a rows: (feature set, label, colour, filled?)
FEATURE_SETS = [
    ("fixed", "Granite, fixed", GRANITE, True),
    ("advected", "Granite, advected", GRANITE, False),
    ("raw_fixed", "Raw bands, fixed", RAW, True),
    ("raw_advected", "Raw bands, advected", RAW, False),
    ("metadata", "Position + season", META, True),
]
TARGETS = [("temperature_c", "Temperature"), ("phosphate_umol_kg", "Phosphate"), ("nitrate_nitrite_umol_l", "Nitrate + nitrite")]


def main() -> None:
    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": RULE, "axes.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": MUTED,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })  # fmt: skip
    fig = plt.figure(figsize=(7.5, 8.4))
    panel_insitu(fig.add_axes([0.17, 0.555, 0.33, 0.43]))
    panel_ko_distribution(fig.add_axes([0.66, 0.60, 0.32, 0.36]))
    panel_top_kos(fig.add_axes([0.47, 0.055, 0.50, 0.44]))
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"predictability.{suffix}", dpi=300)
    print(f"Wrote {OUT / 'predictability.png'} and .pdf")


def style(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=2, width=0.6)


def panel_insitu(ax) -> None:
    folds = pd.read_csv(INSITU)
    style(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    yticks, ylabels, y = [], [], 0.0
    for target, target_label in TARGETS:
        rows = folds[folds["target"] == target]
        for name, label, colour, filled in FEATURE_SETS:
            r2 = rows.loc[rows["feature_set"] == name, "r2"].to_numpy()
            below = r2 < R2_FLOOR
            jitter = np.linspace(-0.12, 0.12, len(r2))
            ax.scatter(np.clip(r2[~below], R2_FLOOR, None), y + jitter[~below], s=7, color=colour, alpha=0.45, lw=0, zorder=2)
            ax.scatter(np.full(below.sum(), R2_FLOOR), y + jitter[below], s=12, marker="<", color=colour, alpha=0.7, lw=0, zorder=2)
            mean, sd = r2.mean(), r2.std(ddof=1)  # sample s.d., as in the comparison tables
            lo, hi = max(mean - sd, R2_FLOOR), min(mean + sd, 1.0)
            ax.plot([lo, hi], [y, y], color=colour, lw=1.2, zorder=3, solid_capstyle="round")
            ax.scatter([max(mean, R2_FLOOR)], [y], s=34, marker="o" if mean >= R2_FLOOR else "<",
                       facecolor=colour if filled else "white", edgecolor=colour, lw=1.2, zorder=4)  # fmt: skip
            ax.text(1.02, y, f"{mean:.2f}", va="center", ha="left", fontsize=6.5, color=MUTED, transform=ax.get_yaxis_transform())
            yticks.append(y)
            ylabels.append(label)
            y += 1
        ax.text(-0.02, y - len(FEATURE_SETS) - 0.75, target_label, ha="right", va="center", fontsize=8, color=INK,
                fontweight="bold", transform=ax.get_yaxis_transform())  # fmt: skip
        y += 0.9
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_ylim(y - 0.6, -1.4)
    ax.set_xlim(R2_FLOOR - 0.05, 1.0)
    ax.axvline(0, color=RULE, lw=0.6, zorder=1)
    ax.set_xlabel("R² on held-out cruise")
    ax.text(-0.5, 1.0, "a", transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="top")


def panel_ko_distribution(ax) -> None:
    screen = pd.read_csv(KO_SCREEN)
    style(ax)
    bins = np.linspace(-0.6, 0.8, 57)
    for column, label, colour in (("r2_pooled_embedding_shuffled", "Shuffled (chance)", NULL),
                                  ("r2_pooled_metadata", "Position + season", META),
                                  ("r2_pooled_embedding", "Granite, fixed", GRANITE)):  # fmt: skip
        values = np.clip(screen[column].to_numpy(), bins[0], bins[-1])
        ax.hist(values, bins=bins, histtype="step", lw=1.4, color=colour, label=label)
    ax.axvline(0, color=RULE, lw=0.6, zorder=0)
    ax.set_xlabel("Pooled R² per KO, cruises held out")
    ax.text(0.0, -0.17, "values below −0.6 drawn at −0.6", transform=ax.transAxes, fontsize=6, color=MUTED)
    ax.set_ylabel(f"KOs (of {len(screen):,})")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right", handlelength=1.2)
    ax.text(-0.3, 1.07, "c", transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="top")


def panel_top_kos(ax) -> None:
    screen = pd.read_csv(KO_SCREEN).sort_values("r2_pooled_embedding", ascending=False).head(N_TOP_KOS)
    names = ko_names()
    style(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    y = np.arange(len(screen))
    granite, meta = screen["r2_pooled_embedding"].to_numpy(), screen["r2_pooled_metadata"].to_numpy()
    for yi, g, m in zip(y, granite, meta):
        ax.plot([m, g], [yi, yi], color=RULE, lw=1.2, zorder=1)
    ax.scatter(meta, y, s=26, color=META, zorder=2, label="Position + season", edgecolor="white", lw=0.8)
    ax.scatter(granite, y, s=26, color=GRANITE, zorder=3, label="Granite, fixed", edgecolor="white", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([label_ko(ko, names) for ko in screen["ko"]], fontsize=6.5)
    ax.set_ylim(len(screen) - 0.4, -0.8)
    ax.set_xlim(-0.45, 0.8)
    ax.axvline(0, color=RULE, lw=0.6, zorder=0)
    ax.set_xlabel("Pooled R², 13 cruises held out in turn")
    ax.tick_params(axis="y", labelcolor=INK)
    ax.legend(frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(0.0, 0.99), ncol=2, handletextpad=0.3)
    ax.text(-0.92, 1.03, "b", transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="top")


def ko_names() -> dict[str, str]:
    if not KO_NAMES.exists():
        KO_NAMES.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen("https://rest.kegg.jp/list/ko", timeout=120) as resp:
            KO_NAMES.write_bytes(resp.read())
    table = pd.read_csv(KO_NAMES, sep="\t", header=None, names=["ko", "name"])
    return dict(zip(table["ko"], table["name"]))


def label_ko(ko: str, names: dict[str, str]) -> str:
    """'K05377' -> 'cpeB  phycoerythrin beta chain  K05377'."""
    full = names.get(ko)
    if not isinstance(full, str):
        return f"{ko} (no KEGG name)"
    symbol, _, description = full.partition("; ")
    symbol = symbol.split(",")[0].strip()
    description = re.sub(r"\s*\[EC:[^\]]*\]", "", description or symbol)
    description = description if len(description) <= 42 else description[:40].rstrip(" ,") + "…"
    if symbol == ko or symbol.startswith("E") and "." in symbol:
        symbol = ""
    return f"{symbol + '  ' if symbol else ''}{description}  {ko}"


if __name__ == "__main__":
    main()
