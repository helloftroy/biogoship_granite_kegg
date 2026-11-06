#!/usr/bin/env python3
"""Figure: how much the Granite embedding adds beyond in-situ temperature, over every KO.

a  Per KO: pooled R2 from in-situ temperature alone (x) against the Granite
   embedding (y), each cruise held out in turn, best of the learners per KO
   (temperature: ridge, spline + ridge, kNN; embedding: ridge, kNN). KOs where the
   embedding clearly wins are blue, where temperature clearly wins magenta.
b  Distribution of R2(embedding) - R2(temperature) over all KOs.

"Clearly wins": a difference of at least WIN_MARGIN in pooled R2 and the same
winner in at least WIN_CRUISES of the 13 held-out cruises (per-cruise R2, same
learners).

Input:   data/results/ko_temperature_all/ (scripts/test_ko_temperature.py --kos all
         --learners ridge spline knn --feature-sets temp embedding)
Outputs: figures/output/temperature_gap.png and .pdf;
         data/results/ko_temperature_all/ko_temperature_gap.csv (per KO: both R2,
         difference, cruises won by the embedding, verdict)

Usage:
    python figures/make_temperature_gap.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results" / "ko_temperature_all"
OUT = ROOT / "figures" / "output"

INK, MUTED, RULE = "#0b0b0b", "#52514e", "#d6d5d0"
GRANITE, TEMP, NEUTRAL = "#2a78d6", "#c0457e", "#b5b4af"
WIN_MARGIN = 0.1
WIN_CRUISES = 9
FLOOR = -0.6  # values below are drawn at the floor


def main() -> None:
    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": RULE, "axes.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": MUTED,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })  # fmt: skip
    gap = per_ko()
    gap.round(4).to_csv(RESULTS / "ko_temperature_gap.csv")
    fig = plt.figure(figsize=(7.5, 3.4))
    panel_scatter(fig.add_axes([0.08, 0.14, 0.36, 0.80]), gap)
    panel_histogram(fig.add_axes([0.58, 0.14, 0.40, 0.80]), gap)
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"temperature_gap.{suffix}", dpi=300)
    counts = gap["verdict"].value_counts()
    print(f"KOs: {len(gap):,}; " + ", ".join(f"{k}: {v:,}" for k, v in counts.items()))
    print(f"Wrote {OUT / 'temperature_gap.png'} and .pdf, {RESULTS / 'ko_temperature_gap.csv'}")


def per_ko() -> pd.DataFrame:
    scores = pd.read_csv(RESULTS / "ko_temperature.csv")
    best = scores.loc[scores.groupby(["ko", "feature_set"])["r2_pooled"].idxmax()]
    pooled = best.pivot(index="ko", columns="feature_set", values="r2_pooled")
    learner = best.set_index(["ko", "feature_set"])["learner"]
    cruise = pd.read_csv(RESULTS / "ko_temperature_per_cruise.csv")
    cruise = cruise[cruise["learner"].to_numpy() == learner.reindex(pd.MultiIndex.from_frame(cruise[["ko", "feature_set"]])).to_numpy()]
    wide = cruise.pivot_table(index=["ko", "held_out"], columns="feature_set", values="r2")
    won = (wide["embedding"] > wide["temp"]).groupby(level="ko").sum()
    cruises = wide.index.get_level_values("held_out").nunique()
    gap = pd.DataFrame({"r2_temp": pooled["temp"], "r2_embedding": pooled["embedding"],
                        "learner_temp": learner.xs("temp", level=1), "learner_embedding": learner.xs("embedding", level=1)})  # fmt: skip
    gap["difference"] = gap["r2_embedding"] - gap["r2_temp"]
    gap["cruises_won_by_embedding"] = won.reindex(gap.index)
    gap["cruises"] = cruises
    embedding_wins = (gap["difference"] >= WIN_MARGIN) & (gap["cruises_won_by_embedding"] >= WIN_CRUISES)
    temp_wins = (gap["difference"] <= -WIN_MARGIN) & (cruises - gap["cruises_won_by_embedding"] >= WIN_CRUISES)
    gap["verdict"] = np.select([embedding_wins, temp_wins], ["embedding clearly better", "temperature clearly better"], "neither clearly")
    return gap.sort_values("difference", ascending=False)


def style(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=2, width=0.6)


def panel_scatter(ax, gap: pd.DataFrame) -> None:
    style(ax)
    x, y = gap["r2_temp"].clip(lower=FLOOR), gap["r2_embedding"].clip(lower=FLOOR)
    for verdict, colour, size, z in (("neither clearly", NEUTRAL, 3, 1), ("temperature clearly better", TEMP, 5, 2),
                                     ("embedding clearly better", GRANITE, 5, 3)):  # fmt: skip
        pick = gap["verdict"] == verdict
        ax.scatter(x[pick], y[pick], s=size, color=colour, lw=0, alpha=0.75, zorder=z,
                   label=f"{verdict[0].upper() + verdict[1:]} ({int(pick.sum()):,})")  # fmt: skip
    ax.plot([FLOOR, 0.85], [FLOOR, 0.85], color=MUTED, lw=0.7, ls=(0, (3, 2)), zorder=0)
    ax.axhline(0, color=RULE, lw=0.6, zorder=0)
    ax.axvline(0, color=RULE, lw=0.6, zorder=0)
    ax.set_xlim(FLOOR - 0.03, 0.85)
    ax.set_ylim(FLOOR - 0.03, 0.85)
    ax.set_aspect("equal")
    ax.set_xlabel("R², in-situ temperature alone")
    ax.set_ylabel("R², Granite embedding")
    ax.legend(frameon=False, fontsize=6.5, loc="upper left", handletextpad=0.2, markerscale=2.2, borderaxespad=0.1)
    ax.text(-0.22, 1.0, "a", transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="top")


def panel_histogram(ax, gap: pd.DataFrame) -> None:
    style(ax)
    bins = np.arange(-0.4, 0.46, 0.02)
    diff = gap["difference"].clip(bins[0], bins[-1] - 1e-9)
    counts, edges = np.histogram(diff, bins=bins)
    centres = (edges[:-1] + edges[1:]) / 2
    colours = np.where(centres >= WIN_MARGIN, GRANITE, np.where(centres <= -WIN_MARGIN, TEMP, NEUTRAL))
    ax.bar(edges[:-1], counts, width=np.diff(edges) * 0.9, align="edge", color=colours, lw=0)
    ax.axvline(0, color=MUTED, lw=0.7)
    ax.set_xlabel("R² embedding − R² temperature, per KO")
    ax.set_ylabel(f"KOs (of {len(gap):,})")
    top = counts.max()
    better, worse = int((gap["difference"] >= WIN_MARGIN).sum()), int((gap["difference"] <= -WIN_MARGIN).sum())
    ax.text(0.30, top * 0.55, f"embedding ahead\nby ≥ {WIN_MARGIN:g}: {better:,}", color=INK, fontsize=6.5, ha="center")
    ax.text(-0.30, top * 0.55, f"temperature ahead\nby ≥ {WIN_MARGIN:g}: {worse:,}", color=INK, fontsize=6.5, ha="center")
    ax.text(0.01, top * 0.98, f"median {gap['difference'].median():+.2f}", color=MUTED, fontsize=6.5, ha="left", va="top")
    ax.text(-0.17, 1.0, "b", transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="top")


if __name__ == "__main__":
    main()
