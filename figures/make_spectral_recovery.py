#!/usr/bin/env python3
"""Figure: observed vs embedding-recovered OLCI spectra for four contrasting samples.

For every tile, a linear model (ridge regression) reads the tile's mean value in
each channel back out of its 512-number Granite embedding. Predictions are
cross-validated with whole cruises held out, so each recovered spectrum comes
from a model that never saw that cruise. Four samples are shown, chosen to be as
different as possible: the clearest water, the greenest water, the warmest and
the coldest (best fixed-position tile, at least 95% clear, with SST).

Inputs:  data/embeddings/tile_embeddings.npy, tile_embeddings.csv,
         tile_raw_features.csv, sample_embeddings.csv
Outputs: figures/output/spectral_recovery.png and .pdf,
         figures/output/spectral_recovery_tiles.csv (the four tiles, for imagery)

Usage:
    python figures/make_spectral_recovery.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, cross_val_predict

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from biogoship import tiles  # noqa: E402

EMB = ROOT / "data" / "embeddings"
OUT = ROOT / "figures" / "output"
# OLCI band centre wavelengths (nm), in tiles.OLCI_BANDS order.
WAVELENGTH_NM = [400, 412.5, 442.5, 490, 510, 560, 620, 665, 673.75, 681.25, 708.75, 753.75, 778.75, 865, 885, 1020]
INK, MUTED, RULE, RECOVERED = "#0b0b0b", "#52514e", "#d6d5d0", "#2a78d6"
SST_K = (tiles.SST_MEAN_K, tiles.SST_STD_K)


def main() -> None:
    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": RULE, "axes.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": MUTED,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })  # fmt: skip
    embeddings = np.load(EMB / "tile_embeddings.npy")
    meta = pd.read_csv(EMB / "tile_embeddings.csv")
    raw = pd.read_csv(EMB / "tile_raw_features.csv").set_index("tile_path").reindex(meta["tile_path"]).reset_index()
    samples = pd.read_csv(EMB / "sample_embeddings.csv")

    channels = [f"{band}_mean" for band in tiles.OLCI_BANDS] + ["SST_mean"]
    cruise = meta["site_id"].str.extract(r"site_([A-Z0-9]+)_")[0]
    recovered = pd.DataFrame(index=meta.index, columns=channels, dtype=float)
    for column in channels:
        y = raw[column].to_numpy()
        ok = np.isfinite(y)
        recovered.loc[ok, column] = cross_val_predict(
            RidgeCV(alphas=np.logspace(-2, 4, 13)), embeddings[ok], y[ok], groups=cruise[ok], cv=GroupKFold(5)
        )

    chosen = choose_tiles(meta, raw, samples)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.1, top=0.9, wspace=0.08, hspace=0.32)
    for ax, (letter, row) in zip(axes.ravel(), chosen.iterrows()):
        i = row["index"]
        observed = to_reflectance(raw.loc[i, channels[:16]].to_numpy(float))
        predicted = to_reflectance(recovered.loc[i, channels[:16]].to_numpy(float))
        ax.plot(WAVELENGTH_NM, observed, color=INK, lw=1.6, marker="o", ms=3.2, zorder=3)
        ax.plot(WAVELENGTH_NM, predicted, color=RECOVERED, lw=1.6, ls=(0, (4, 2)), marker="o", ms=3.2,
                markerfacecolor="white", markeredgewidth=1.0, zorder=4)  # fmt: skip
        ax.axhline(0, color=RULE, lw=0.6, zorder=1)
        sst_obs = raw.loc[i, "SST_mean"] * SST_K[1] + SST_K[0] - 273.15
        sst_rec = recovered.loc[i, "SST_mean"] * SST_K[1] + SST_K[0] - 273.15
        ax.text(0.02, 1.03, f"{letter}  {row['label']}", transform=ax.transAxes, ha="left", va="bottom",
                fontsize=8, fontweight="bold", color=INK)  # fmt: skip
        ax.text(0.98, 0.95, f"{row['section']}  {fmt_latlon(row['tile_lat'], row['tile_lon'])}\n{row['tile_date'][:10]}\n"
                f"SST {sst_obs:.1f} °C observed, {sst_rec:.1f} °C recovered",
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color=MUTED, linespacing=1.4)  # fmt: skip
        ax.grid(axis="y", color=RULE, lw=0.4)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Wavelength (nm)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Water-leaving reflectance")
    handles = [
        plt.Line2D([], [], color=INK, lw=1.6, marker="o", ms=3.2, label="Observed (tile mean)"),
        plt.Line2D([], [], color=RECOVERED, lw=1.6, ls=(0, (4, 2)), marker="o", ms=3.2, markerfacecolor="white",
                   label="Recovered from the Granite embedding (cruise held out)"),
    ]  # fmt: skip
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=7.5, bbox_to_anchor=(0.54, 1.0))

    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"spectral_recovery.{suffix}", dpi=300, bbox_inches="tight", pad_inches=0.03)
    keep = ["label", "sample_id", "section", "tile_id", "tile_path", "granule_id", "tile_date", "tile_lat", "tile_lon",
            "valid_fraction", "blue_green", "sst_c"]  # fmt: skip
    chosen.rename_axis("panel").reset_index()[["panel", *keep]].to_csv(OUT / "spectral_recovery_tiles.csv", index=False)
    print(chosen[["label", "sample_id", "tile_path", "tile_date", "tile_lat", "tile_lon"]].to_string())
    print(f"Wrote {OUT / 'spectral_recovery.png'}")


def to_reflectance(normalised: np.ndarray) -> np.ndarray:
    return tiles.OLCI_MEANS + normalised * tiles.OLCI_STDS


def choose_tiles(meta: pd.DataFrame, raw: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """Clearest, greenest, warmest and coldest among the best fixed tiles that are >=95% clear with SST.

    Clear vs green is ranked by the 443/560 nm reflectance ratio (the core of standard
    chlorophyll algorithms), only among tiles whose visible reflectances are all positive
    and whose 560 nm reflectance is above 0.002, so the ratio is not inflated by a tiny or
    negative denominator (a sign of poor atmospheric correction).
    """
    best = samples[(samples["is_best"] == 1) & (samples["position"] == "fixed")].drop_duplicates("tile_path")
    pool = meta.reset_index().merge(best[["tile_path", "sample_id", "section"]], on="tile_path")
    pool = pool[(pool["valid_fraction"] >= 0.95) & (pool["sst_valid_fraction"] > 0)].copy()
    visible = np.column_stack([to_reflectance_band(raw.loc[pool["index"], f"{b}_mean"], k)
                               for k, b in enumerate(tiles.OLCI_BANDS[:8])])  # fmt: skip
    blue, green = visible[:, 2], visible[:, 5]  # Oa03 (442.5 nm), Oa06 (560 nm)
    pool["blue_green"] = blue / green
    pool = pool[(visible > 0).all(axis=1) & (green > 0.002)].copy()
    pool["sst_c"] = (raw.loc[pool["index"], "SST_mean"] * SST_K[1] + SST_K[0] - 273.15).to_numpy()
    picks, used = [], set()
    for label, column, ascending in (("Clearest water", "blue_green", False), ("Greenest water", "blue_green", True),
                                     ("Warmest water", "sst_c", False), ("Coldest water", "sst_c", True)):  # fmt: skip
        ranked = pool.sort_values(column, ascending=ascending)
        row = ranked[~ranked["site_id"].isin(used)].iloc[0]
        used.add(row["site_id"])
        picks.append({**row.to_dict(), "label": label})
    return pd.DataFrame(picks, index=list("abcd"))


def to_reflectance_band(normalised: pd.Series, band: int) -> pd.Series:
    return tiles.OLCI_MEANS[band] + normalised * tiles.OLCI_STDS[band]


def fmt_latlon(lat: float, lon: float) -> str:
    return f"{abs(lat):.1f}°{'N' if lat >= 0 else 'S'} {abs(lon):.1f}°{'E' if lon >= 0 else 'W'}"


if __name__ == "__main__":
    main()
