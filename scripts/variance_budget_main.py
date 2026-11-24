#!/usr/bin/env python3
"""Variance budget of all KOs (CLR) under the main model, as a waterfall (RESULTS.md R63; update of R43).

Main model (R59/R60): 32 Granite PCs, ridge, penalty by leave-cruise-out inside the training cruises, one cruise held
out; technical variance removed in-sample. 1,215 tiles x 3,933 KOs. Sums of squares over tiles and KOs (raw CLR units;
the same with every KO standardised is written too). The parts add to 100 % of the total variance across tiles:

  1  explained by the satellite model      SS(adjusted) - SS(held-out residual)
  2  sequencing platform (upper bound)     OLS of the targets on platform group (NovaSeq 6000 vs NovaSeq X cruises)
  3  sequencing depth                      the rest of the in-sample technical fit (platform + 8-knot cubic spline of
                                           log10 KO reads per tile) beyond platform
  4  other cruise offsets                  cruise means of the held-out residual
  5  within-cruise axis 1, 6  axis 2       the two leading principal components of the within-cruise residual (KOs
                                           standardised), regressed out (in-sample); described by their correlation
                                           with the R36 axes and with log(eukaryote / prokaryote reads)
  7  within-tile variability               sigma2 per KO from 305 same-tile sample pairs x sum over tiles of 1/k
                                           (capped at what is left)
  8  remaining within cruise               the rest

Outputs: data/results/main_model/variance_budget.csv; figures/output/variance_budget.png and .pdf
Usage: OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/variance_budget_main.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import main_model as mm  # noqa: E402
from water_mass_test import features  # noqa: E402

LABELS = {
    "explained": "Explained by Granite PCs (cruise held out)",
    "platform": "Sequencing platform (NovaSeq 6000 vs X; upper bound)",
    "depth": "Sequencing depth",
    "offsets": "Other cruise offsets (unknown)",
    "W1": "Axis 1 within cruise (eukaryotes vs prokaryotes)",
    "W2": "Axis 2 within cruise (phototrophs; driver unknown)",
    "within_tile": "Within-tile variability (unpredictable floor)",
    "rest": "Remaining within cruise (unknown)",
}


def centre(m, g):
    return m - pd.DataFrame(m).groupby(g).transform("mean").to_numpy()


def fit_out(m, x):
    beta, *_ = np.linalg.lstsq(x, m, rcond=None)
    return x @ beta


def main():
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    fs = features()
    d = mm.setup("clr_all", fs, ks)
    y, ya, cr, unit, rows = d["y"], d["ya"], d["cr"], d["unit"], d["rows"]
    pred = np.load(mm.OUT / "predictions_clr_all_granite_cv.npz")["predicted"].astype(float)
    yc = y - y.mean(0)
    plat = d["plat"].astype(float)[:, None]
    pf = fit_out(yc, np.hstack([np.ones((len(y), 1)), plat]))
    res = ya - pred  # not centred: its overall mean goes into the cruise offsets
    cm = pd.DataFrame(res).groupby(cr).transform("mean").to_numpy()
    w = res - cm
    sd = np.where(yc.std(0) > 0, yc.std(0), 1)
    u, s, vt = np.linalg.svd(w / sd, full_matrices=False)
    scores = u[:, :2] * s[:2]
    old = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "scores.csv", index_col="tile_key").reindex(unit.index)
    cdat = mm.setup("C", fs, ks)
    c_w = centre(cdat["ya"], cr)[:, 0]
    desc = {}
    for i in range(2):
        sc = scores[:, i]
        r_old = {k: np.corrcoef(sc, old[k])[0, 1] for k in ("W1", "W2") if old[k].notna().all()}
        desc[f"W{i + 1}"] = {"r_with_R36": r_old, "r_with_C_within_cruise": np.corrcoef(sc, c_w)[0, 1]}
        load = vt[i] * s[i]
        top = pd.Series(load, index=d["kos"])
        desc[f"W{i + 1}"]["top_pos"] = top.nlargest(8).index.tolist()
        desc[f"W{i + 1}"]["top_neg"] = top.nsmallest(8).index.tolist()
    st = pd.read_csv(mm.TARGETS["clr_all"], index_col="sample_key")[d["kos"]]
    members = rows.groupby("tile_key")["sample_id"].agg(list).reindex(unit.index)
    diffs = []
    for samples in members[members.str.len() >= 2]:
        v = st.reindex([x for x in samples if x in st.index]).to_numpy()
        diffs += [v[i] - v[j] for i in range(len(v)) for j in range(i + 1, len(v))]
    sigma2 = (np.array(diffs) ** 2).mean(0) / 2
    out = {}
    for scale, sw in (("raw", np.ones(y.shape[1])), ("standardised", 1 / sd)):
        ss = lambda m: ((m * sw) ** 2).sum()  # noqa: E731
        tot = ss(yc)
        b = {"explained": ss(ya - ya.mean(0)) - ss(res), "platform": ss(pf), "depth": ss(d["tech"]) - ss(pf),
             "offsets": ss(cm)}  # fmt: skip
        left = w
        for i in range(2):
            f = fit_out(left, centre(scores[:, [i]], cr))
            b[f"W{i + 1}"] = ss(f)
            left = left - f
        noise = (sigma2 * sw**2).sum() * (1 / d["k"]).sum()
        b["within_tile"] = min(noise, ss(left))
        b["rest"] = ss(left) - b["within_tile"]
        out[scale] = pd.Series(b) / tot
    t = pd.DataFrame(out)
    t["label"] = [LABELS[k] for k in t.index]
    t.round(5).to_csv(mm.OUT / "variance_budget.csv")
    print((t[["raw", "standardised"]] * 100).round(2).assign(label=t["label"]).to_string())
    print("sums:", t["raw"].sum().round(6), t["standardised"].sum().round(6), f"; pairs {len(diffs)}")
    names = pd.read_csv(REPO / "data/reference/kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "n"]).set_index("ko")["n"]
    for k, v in desc.items():
        print(k, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items() if not kk.startswith("top")} | {
            "r_with_R36": {a: round(b, 3) for a, b in v["r_with_R36"].items()}})
        for side in ("top_pos", "top_neg"):
            print("  ", side, "; ".join(str(names.get(x, x))[:40] for x in v[side]))
    figure(t)


def figure(t):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})  # fmt: skip
    kind = {"explained": ("#0072B2", None), "platform": ("#E69F00", None), "depth": ("#E69F00", None),
            "offsets": ("#BBBBBB", "////"), "W1": ("#009E73", None), "W2": ("#009E73", "...."),
            "within_tile": ("#555555", None), "rest": ("#BBBBBB", "////")}  # fmt: skip
    b = t["raw"]
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    left = 0.0
    for i, k in enumerate(b.index):
        v = b[k] * 100
        col, hatch = kind[k]
        ax.barh(i, v, left=left, height=0.62, color=col if hatch is None else "white", edgecolor=col, hatch=hatch, linewidth=0.8)
        ax.text(left + v + 0.8, i, f"{v:.1f}%", va="center", fontsize=7)
        if i < len(b) - 1:
            ax.plot([left + v] * 2, [i + 0.31, i + 0.69], color="#999999", lw=0.5)
        left += v
    ax.set_yticks(range(len(b)))
    ax.set_yticklabels([LABELS[k] for k in b.index])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of KO CLR variance (%), cumulative from left")
    ax.legend(handles=[Patch(facecolor="#0072B2", label="Model"), Patch(facecolor="#E69F00", label="Technical"),
                       Patch(facecolor="#009E73", label="Biology"), Patch(facecolor="#555555", label="Within-tile floor"),
                       Patch(facecolor="white", edgecolor="#BBBBBB", hatch="////", label="Unknown")],
              frameon=False, fontsize=6.5, loc="upper right")  # fmt: skip
    fig.tight_layout()
    for suf in ("png", "pdf"):
        fig.savefig(mm.FIG / f"variance_budget.{suf}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
