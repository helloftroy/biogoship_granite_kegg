#!/usr/bin/env python3
"""Figure: how well each sample is predicted with its cruise held out, on the globe, and why.

a  Every sample (distinct tile) coloured by its prediction skill when its cruise is
   held out: 1 - mean over KOs of (error^2 / variance of that KO over all samples).
   Averaged over samples it equals the mean pooled R2 of the KOs, so it is "R2 at a
   point". KOs: those with pooled R2 >= 0.25 for the model (920 of 3,921), so the
   skill reflects KOs the model can predict at all. Model: 32 PCs of the Granite
   embedding + in-situ temperature, ridge (R17 pca32_temp).
b  Skill against the distance to the nearest sample of another cruise (does a cruise
   that crosses another one's track get predicted better?). Median and quartiles per bin.
c  Skill against how unusual the sample's conditions are for the other cruises: mean
   distance to its 5 nearest samples from other cruises in (first 8 PCs of the
   embedding, in-situ temperature), all standardised. Quintiles.

Input:   data/results/ko_temperature_all_ostia_v2/ (predictions_pca32_temp_ridge.npz,
         prediction_rows.csv, ko_temperature.csv), the OSTIA embeddings (for c)
Outputs: figures/output/skill_map.png and .pdf; data/results/ko_temperature_all_ostia_v2/sample_skill.csv

Usage:
    python figures/make_skill_map.py
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
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_cruise_line import INK, LAND_FILL, MUTED, RULE, land_rings  # noqa: E402

RESULTS = ROOT / "data" / "results" / "ko_temperature_all_ostia_v2"
OUT = ROOT / "figures" / "output"
SKILL_CMAP = LinearSegmentedColormap.from_list("skill", ["#b2182b", "#e9a38f", "#e8e6e1", "#8fb7dc", "#1f5fa6"])
SKILL_NORM = Normalize(-1, 1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--min-r2", type=float, default=0.25)
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[ROOT / "data" / "embeddings_sst_l4", ROOT / "data" / "faire" / "embeddings_sst_l4"])
    args = parser.parse_args()
    plt.rcParams.update({
        "font.family": ["Arial", "DejaVu Sans"], "font.size": 8, "axes.labelsize": 7.5, "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5, "axes.edgecolor": RULE, "axes.linewidth": 0.6, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.labelcolor": MUTED, "pdf.fonttype": 42, "svg.fonttype": "none",
    })  # fmt: skip

    rows = pd.read_csv(RESULTS / "prediction_rows.csv")
    z = np.load(RESULTS / "predictions_pca32_temp_ridge.npz", allow_pickle=True)
    scores = pd.read_csv(RESULTS / "ko_temperature.csv")
    scores = scores[(scores["feature_set"] == "pca32_temp") & (scores["learner"] == "ridge")]
    good = set(scores.loc[scores["r2_pooled"] >= args.min_r2, "ko"])
    keep = np.array([k in good for k in z["ko"]])
    y, p = z["observed"][:, keep], z["predicted"][:, keep]
    skill = 1 - ((y - p) ** 2 / y.var(axis=0)).mean(axis=1)
    cruise = rows["section"].to_numpy()
    other = cruise[:, None] != cruise[None, :]

    lat, lon = np.radians(rows["sample_lat"].to_numpy()), np.radians(rows["sample_lon"].to_numpy())
    km = 6371 * np.arccos(np.clip(np.sin(lat)[:, None] * np.sin(lat) + np.cos(lat)[:, None] * np.cos(lat)
                                  * np.cos(lon[:, None] - lon), -1, 1))  # fmt: skip
    nearest_km = np.where(other, km, np.inf).min(axis=1)

    from tune_figure_kos import load

    ns = argparse.Namespace(embeddings_dirs=args.embeddings_dirs, targets=ROOT / "data" / "targets" / "ko_clr.csv")
    emb, temp, _, _, _, unit = load(ns, [str(z["ko"][0])])
    assert (unit["tile_path"].to_numpy() == rows["tile_path"].to_numpy()).all()
    x = StandardScaler().fit_transform(np.hstack([PCA(8, random_state=0).fit_transform(StandardScaler().fit_transform(emb)), temp]))
    dist = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(-1))
    novelty = np.sort(np.where(other, dist, np.inf), axis=1)[:, :5].mean(axis=1)

    table = rows[["sample_id", "section", "sample_lat", "sample_lon"]].assign(skill=skill.round(4),
        nearest_other_cruise_km=nearest_km.round(1), novelty=novelty.round(4))  # fmt: skip
    table.to_csv(RESULTS / "sample_skill.csv", index=False)
    rho_km, rho_nov = spearmanr(nearest_km, skill)[0], spearmanr(novelty, skill)[0]
    print(f"{keep.sum()} KOs; mean skill {skill.mean():.3f}; Spearman skill vs distance {rho_km:.2f}, vs novelty {rho_nov:.2f}")

    fig = plt.figure(figsize=(7.5, 6.4))
    ax = fig.add_axes([0.05, 0.445, 0.92, 0.545])
    for ring in land_rings():
        ax.fill(ring[:, 0], ring[:, 1], color=LAND_FILL, lw=0)
    order = np.argsort(-np.abs(skill))[::-1]  # draw the extreme ones last
    ax.scatter(rows["sample_lon"].to_numpy()[order], rows["sample_lat"].to_numpy()[order], s=9,
               c=SKILL_CMAP(SKILL_NORM(skill[order])), linewidths=0.25, edgecolors="white", zorder=3)  # fmt: skip
    ax.set_xlim(-180, 180)
    ax.set_ylim(-75, 82)
    ax.set_aspect("equal")
    ax.set_xticks(np.arange(-150, 181, 50))
    ax.set_xticklabels([f"{abs(v)}°{'W' if v < 0 else 'E' if v > 0 else ''}" for v in np.arange(-150, 181, 50)])
    ax.set_yticks(np.arange(-60, 81, 30))
    ax.set_yticklabels([f"{abs(v)}°{'S' if v < 0 else 'N' if v > 0 else ''}" for v in np.arange(-60, 81, 30)])
    ax.tick_params(length=2)
    tag(ax, "a", -0.045)

    kax = fig.add_axes([0.06, 0.375, 0.22, 0.012])
    kax.imshow(np.linspace(-1, 1, 200)[None, :], aspect="auto", cmap=SKILL_CMAP, extent=[-1, 1, 0, 1])
    kax.set_yticks([])
    kax.set_xticks([-1, 0, 1])
    kax.set_xticklabels(["≤ −1\nworse than the mean", "0", "1\nperfect"], fontsize=6)
    kax.tick_params(length=2, pad=1.5)
    for side in kax.spines.values():
        side.set_visible(False)
    fig.text(0.17, 0.393, "Skill, cruise held out", ha="center", va="bottom", fontsize=6.5, color=MUTED)

    # per-cruise table: median skill, share of samples within 200 km of another cruise
    per = table.groupby("section").agg(skill=("skill", "median"), near=("nearest_other_cruise_km", lambda v: (v < 200).mean()))
    per = per.sort_values("skill", ascending=False)
    top, step = 0.295, 0.0205
    fig.text(0.075, top + 0.012, "Cruise", fontsize=6.3, color=MUTED, va="bottom")
    fig.text(0.225, top + 0.012, "Median\nskill", fontsize=6.3, color=MUTED, va="bottom", ha="right", linespacing=1.0)
    fig.text(0.315, top + 0.012, "Within\n200 km", fontsize=6.3, color=MUTED, va="bottom", ha="right", linespacing=1.0)
    for r, (name, row) in enumerate(per.iterrows()):
        yy = top - r * step
        fig.patches.append(plt.Circle((0.065, yy), 0.0045, transform=fig.transFigure, color=SKILL_CMAP(SKILL_NORM(row["skill"])),
                                      figure=fig))  # fmt: skip
        fig.text(0.075, yy, name, fontsize=6.3, color=INK, va="center")
        fig.text(0.225, yy, f"{row['skill']:.2f}", fontsize=6.3, color=INK, va="center", ha="right")
        fig.text(0.315, yy, f"{row['near']:.0%}", fontsize=6.3, color=INK, va="center", ha="right")

    bins = [0, 100, 200, 500, 1000, 2000, np.inf]
    labels = ["<100", "100–\n200", "200–\n500", "500–\n1000", "1000–\n2000", ">2000"]
    bax = fig.add_axes([0.43, 0.115, 0.24, 0.20])
    binned(bax, pd.cut(nearest_km, bins, labels=labels), skill)
    bax.set_xlabel("Distance to the nearest sample\nof another cruise (km)")
    bax.set_ylabel("Skill")
    bax.text(0.98, 0.98, f"ρ = {rho_km:.2f}", transform=bax.transAxes, ha="right", va="top", fontsize=6.5, color=INK)
    tag(bax, "b", -0.30)

    cax = fig.add_axes([0.73, 0.115, 0.24, 0.20], sharey=bax)
    binned(cax, pd.qcut(novelty, 5, labels=["most\ntypical", "2", "3", "4", "most\nunusual"]), skill)
    cax.set_xlabel("Conditions compared with other cruises\n(embedding + temperature, quintiles)")
    cax.tick_params(labelleft=False)
    cax.text(0.98, 0.98, f"ρ = {rho_nov:.2f}", transform=cax.transAxes, ha="right", va="top", fontsize=6.5, color=INK)
    tag(cax, "c", -0.10)

    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"skill_map.{suffix}", dpi=300)
    print(f"Wrote {OUT / 'skill_map'}.png and .pdf")


def binned(ax, groups: pd.Categorical, skill: np.ndarray) -> None:
    frame = pd.DataFrame({"g": groups, "s": skill})
    stats = frame.groupby("g", observed=False)["s"].describe()
    x = np.arange(len(stats))
    ax.vlines(x, stats["25%"], stats["75%"], color="#8fb7dc", lw=3)
    ax.scatter(x, stats["50%"], s=18, color="#1f5fa6", zorder=3, linewidths=0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{name}\nn={int(n)}" for name, n in zip(stats.index, stats["count"])], fontsize=5.8)
    ax.set_ylim(-0.1, 1.0)
    ax.axhline(0, color=RULE, lw=0.6)
    ax.grid(axis="y", color=RULE, lw=0.4)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=2)


def tag(ax, letter: str, x: float) -> None:
    ax.text(x, 1.0, letter, transform=ax.transAxes, ha="left", va="top", fontsize=10, fontweight="bold", color=INK)


if __name__ == "__main__":
    main()
