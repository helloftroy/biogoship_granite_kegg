#!/usr/bin/env python3
"""Variance budget of the KO CLR targets and a scree plot of the residual components (RESULTS.md R43).

A first draft: parts still unknown are labelled unknown. Real model (data/targets/ko_clr.csv, A16N Protocol A
excluded; 1,215 tiles x 3,933 KOs); baseline = 32 Granite PCs + temperature spline, ridge, one cruise held out.

Budget (sum of squares over tiles and KOs, raw CLR units; the same with every KO standardised is also written).
Sequential, each step on what the previous left:
  1  explained by the satellite model           total - residual (cross-validated)
  2  sequencing platform                        residual: means of the NovaSeq 6000 and NovaSeq X cruise groups
  3  other cruise offsets (unknown)             residual: cruise means minus the platform means
  4  sequencing depth                           within-cruise residual on a spline of log10 total KO counts
                                                (8 knots, cubic; basis centred per cruise; in-sample OLS per KO)
  5  axis 1 within cruise (W1)                  regression on the R36 W1 score (biology; weak links to
                                                stratification and bloom phase)
  6  axis 2 within cruise (W2)                  regression on the R36 W2 score (biology; driver unknown)
  7  replicate noise (floor)                    pooled same-tile sample-pair variance per KO, sigma2 = mean (a-b)^2/2
                                                over pairs, scaled by sum over tiles of 1/k (k samples averaged per
                                                tile); capped at what is left
  8  remaining within cruise (unknown)          the rest
Caveats: steps 2-3 are an upper bound for platform (the platform groups are also years and regions); W1 / W2
scores are from a PCA of the standardised within-cruise residuals, so their raw-unit shares differ from the
R33 percentages; the replicate pairs are GOMECC4 A/B/C, C13.5 v1/v2 and co-located samples (same water,
separate filters / libraries).

Scree: share of each component of the full (F) and within-cruise (W) standardised residual PCA
(transfer/results/residual_modules_noA/components.csv), with the permutation-null 95th percentile and the
split-half congruence.

Outputs: transfer/results/variance_budget/budget.csv; transfer/figures/variance_budget.png, residual_scree.png.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/variance_budget.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from added_fields_all_kos import cv  # noqa: E402
from tune_figure_kos import load  # noqa: E402

OUT = REPO / "transfer" / "results" / "variance_budget"
FIG = REPO / "transfer" / "figures"
NOVASEQ_6000 = ["I09N", "P18", "I07N", "AMT28", "C13.5/A13.5", "HLY2202-2022"]


def centre(m, g):
    return m - pd.DataFrame(m).groupby(g).transform("mean").to_numpy()


def fit_out(m, x):
    """Least-squares fit of every column of m on the columns of x (both centred); returns fitted values."""
    beta, *_ = np.linalg.lstsq(x, m, rcond=None)
    return x @ beta


def main() -> None:
    from sklearn.preprocessing import SplineTransformer

    OUT.mkdir(parents=True, exist_ok=True)
    args = argparse.Namespace(targets=REPO / "data" / "targets" / "ko_clr.csv",
                              embeddings_dirs=[REPO / "data" / "embeddings_sst_l4", REPO / "data" / "faire" / "embeddings_sst_l4"])
    kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, _, y, cr, unit, rows = load(args, kos, return_rows=True)
    _, pred = cv(emb, temp[:, 0], None, y, cr, np.arange(len(unit)))
    yc = y - y.mean(0)
    res = y - pred
    res = res - res.mean(0)
    sd = yc.std(0)
    out = {}
    for scale_name, s in (("raw", np.ones(y.shape[1])), ("standardised", 1 / np.where(sd > 0, sd, 1))):
        tot = ((yc * s) ** 2).sum()
        ss = lambda m: ((m * s) ** 2).sum()  # noqa: E731
        b = {"explained by the satellite model (CV)": tot - ss(res)}
        plat = np.isin(cr, NOVASEQ_6000)
        pm = pd.DataFrame(res).groupby(plat).transform("mean").to_numpy()
        cm = pd.DataFrame(res).groupby(cr).transform("mean").to_numpy()
        b["sequencing platform (NovaSeq 6000 vs X; upper bound)"] = ss(pm)
        b["other cruise offsets (unknown)"] = ss(cm - pm)
        w = res - cm
        ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
        lt = np.log10(ks["total"])
        logtot = rows["sample_id"].map(lt).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()
        basis = SplineTransformer(n_knots=8, degree=3).fit_transform(logtot[:, None])
        basis = centre(basis, cr)
        d = fit_out(w, basis)
        b["sequencing depth"] = ss(d)
        w2 = w - d
        sc = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "scores.csv", index_col="tile_key").loc[unit.index]
        for comp, label in (("W1", "axis 1 within cruise (bacteria vs eukaryotes; biology)"),
                            ("W2", "axis 2 within cruise (phototrophs; biology, driver unknown)")):
            x = centre(sc[[comp]].to_numpy(), cr)
            f = fit_out(w2, x)
            b[label] = ss(f)
            w2 = w2 - f
        # replicate noise
        multi = rows.groupby("tile_key")["sample_id"].agg(list)
        clr = pd.read_csv(args.targets, index_col="sample_key")
        diffs = []
        for samples in multi[multi.str.len() >= 2]:
            v = clr.loc[[x for x in samples if x in clr.index]].to_numpy()
            for i in range(len(v)):
                for j in range(i + 1, len(v)):
                    diffs.append(v[i] - v[j])
        diffs = np.array(diffs)
        sigma2 = (diffs ** 2).mean(0) / 2
        k = multi.reindex(unit.index).str.len().fillna(1).to_numpy()
        noise = (sigma2 * s ** 2).sum() * (1 / k).sum()
        left = ss(w2)
        b["replicate noise (unpredictable floor)"] = min(noise, left)
        b["remaining within cruise (unknown)"] = left - min(noise, left)
        out[scale_name] = pd.Series(b) / tot
        if scale_name == "raw":
            print(f"replicate pairs: {len(diffs):,} from {int((multi.str.len() >= 2).sum())} tiles")
    budget = pd.DataFrame(out)
    budget.index.name = "component"
    budget.round(4).to_csv(OUT / "budget.csv")
    print((budget * 100).round(1).to_string())
    figures(budget["raw"])


def figures(b: pd.Series) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    kind = {"explained by the satellite model (CV)": ("model", "#0072B2", None),
            "sequencing platform (NovaSeq 6000 vs X; upper bound)": ("technical", "#E69F00", None),
            "other cruise offsets (unknown)": ("unknown", "#BBBBBB", "////"),
            "sequencing depth": ("technical", "#E69F00", None),
            "axis 1 within cruise (bacteria vs eukaryotes; biology)": ("biology", "#009E73", None),
            "axis 2 within cruise (phototrophs; biology, driver unknown)": ("biology", "#009E73", "...."),
            "replicate noise (unpredictable floor)": ("floor", "#555555", None),
            "remaining within cruise (unknown)": ("unknown", "#BBBBBB", "////")}  # fmt: skip
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    left = 0.0
    labels = list(b.index)
    for i, name in enumerate(labels):
        v = b[name] * 100
        _, col, hatch = kind[name]
        ax.barh(i, v, left=left, height=0.62, color=col if hatch is None else "white", edgecolor=col,
                hatch=hatch, linewidth=0.8)
        ax.text(left + v + 0.8, i, f"{v:.1f}%", va="center", fontsize=7)
        if i < len(labels) - 1:
            ax.plot([left + v, left + v], [i + 0.31, i + 0.69], color="#999999", lw=0.5)
        left += v
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels([n[0].upper() + n[1:] for n in labels])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of KO CLR variance (%), cumulative from left")
    from matplotlib.patches import Patch

    ax.legend(handles=[Patch(facecolor="#0072B2", label="Model"), Patch(facecolor="#E69F00", label="Technical"),
                       Patch(facecolor="#009E73", label="Biology"), Patch(facecolor="#555555", label="Noise floor"),
                       Patch(facecolor="white", edgecolor="#BBBBBB", hatch="////", label="Unknown")],
              frameon=False, fontsize=6.5, loc="upper right")  # fmt: skip
    fig.tight_layout()
    fig.savefig(FIG / "variance_budget.png", dpi=300)
    plt.close(fig)

    comp = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "components.csv")
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    for res, col, lab in (("full", "#0072B2", "Full residual (F)"), ("within cruise", "#D55E00", "Within cruise (W)")):
        c = comp[comp["residuals"] == res].sort_values("rank")
        ax.plot(c["rank"], c["share"] * 100, color=col, lw=1.2, label=lab, zorder=2)
        stable = c["split_half_congruence"] >= 0.75
        ax.scatter(c["rank"][stable], c["share"][stable] * 100, s=22, color=col, zorder=3, edgecolor="white", linewidth=0.8)
        ax.scatter(c["rank"][~stable], c["share"][~stable] * 100, s=22, facecolor="white", edgecolor=col, zorder=3, linewidth=0.9)
    null = comp[comp["residuals"] == "full"].sort_values("rank")
    ax.plot(null["rank"], null["null95_share"] * 100, color="#999999", lw=0.8, ls="--", label="Permutation null (95th pct)")
    ax.set_xlabel("Component")
    ax.set_ylabel("Share of standardised residual variance (%)")
    ax.set_xticks(range(1, 11))
    ax.scatter([], [], s=22, color="#555555", label="Reproducible (split-half >= 0.75)")
    ax.scatter([], [], s=22, facecolor="white", edgecolor="#555555", label="Not reproducible")
    ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(FIG / "residual_scree.png", dpi=300)
    plt.close(fig)
    print(f"Figures: {FIG / 'variance_budget.png'}, {FIG / 'residual_scree.png'}")


if __name__ == "__main__":
    main()
