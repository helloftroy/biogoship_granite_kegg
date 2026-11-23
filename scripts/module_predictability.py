#!/usr/bin/env python3
"""Pathway-level predictability: which kinds of function are visible from space (RESULTS.md R61).

Per KO (main model, R59/R60: 32 Granite PCs, ridge, leave-cruise-out penalty, technical variance removed):
  R2        pooled held-out R2 on the technically adjusted CLR (data/results/main_model/ko_temperature.csv)
  ceiling   reproducible share of the KO's non-technical variance = 1 - within-tile / SS(adjusted KO over tiles),
            within-tile = sigma2 (single-sample within-tile variance from the 305 sample pairs) * sum over tiles of 1/k
  fraction  R2 / ceiling: skill as a share of what is predictable at all given the within-tile scatter
Per KEGG module (KEGG REST link/module/ko and br:ko00002, downloaded 2026-11-23 into data/reference/kegg_modules/):
  median fraction over the module's modelled KOs (modules with >= MIN_KOS of them), median R2, median ceiling, and the
  pooled fraction (sum of explained variance / sum of reproducible variance), which is variance-weighted.
Ceilings are floored at CEIL_FLOOR so very noisy KOs do not explode the ratio (count reported).

Outputs: data/results/main_model/ko_ceiling.csv, module_predictability.csv; figures/output/module_predictability.png/pdf
Usage: OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/module_predictability.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import main_model as mm  # noqa: E402
from water_mass_test import features  # noqa: E402

KD = REPO / "data" / "reference" / "kegg_modules"
MIN_KOS, CEIL_FLOOR = 4, 0.1


def ceilings():
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    d = mm.setup("clr_all", features(), ks)
    st = pd.read_csv(mm.TARGETS["clr_all"], index_col="sample_key")[d["kos"]]
    members = d["rows"].groupby("tile_key")["sample_id"].agg(list).reindex(d["unit"].index)
    diffs = []
    for samples in members[members.str.len() >= 2]:
        v = st.reindex([x for x in samples if x in st.index]).to_numpy()
        diffs += [v[i] - v[j] for i in range(len(v)) for j in range(i + 1, len(v))]
    sigma2 = (np.array(diffs) ** 2).mean(0) / 2
    within = sigma2 * (1 / d["k"]).sum()
    ss = ((d["ya"] - d["ya"].mean(0)) ** 2).sum(0)
    return pd.DataFrame({"ko": d["kos"], "within_tile": within, "ss_nontech": ss, "ceiling": 1 - within / ss}).set_index("ko")


def modules():
    link = pd.read_csv(KD / "link_module_ko.tsv", sep="\t", header=None, names=["ko", "module"])
    link["ko"], link["module"] = link["ko"].str[3:], link["module"].str[3:]
    tree = json.loads((KD / "ko00002.json").read_text())
    rows = []
    for top in tree["children"]:  # "Pathway modules" and "Signature modules"
        for l1 in top["children"]:
            for l2 in l1["children"]:
                for l3 in l2["children"]:
                    mid, name = l3["name"].split("  ", 1)
                    rows.append({"module": mid, "name": name.split(" [")[0], "type": top["name"], "class1": l1["name"],
                                 "class2": l2["name"]})
    return link, pd.DataFrame(rows).drop_duplicates("module").set_index("module")


def main():
    path = mm.OUT / "ko_ceiling.csv"
    c = pd.read_csv(path, index_col="ko") if path.exists() else ceilings()
    s = pd.read_csv(mm.OUT / "ko_temperature.csv")
    for fs, lab in (("embedding", "r2"), ("temp", "r2_temp")):
        c[lab] = s[(s["feature_set"] == fs) & (s["learner"] == ("ridge" if fs == "embedding" else "spline"))].set_index("ko")["r2_pooled"]
    c["ceiling_used"] = c["ceiling"].clip(lower=CEIL_FLOOR)
    c["fraction"] = c["r2"] / c["ceiling_used"]
    c["fraction_temp"] = c["r2_temp"] / c["ceiling_used"]
    c.round(5).to_csv(path)
    print(f"KOs {len(c):,}; ceiling median {c['ceiling'].median():.3f} (10-90% {c['ceiling'].quantile(.1):.3f}-"
          f"{c['ceiling'].quantile(.9):.3f}); below floor {CEIL_FLOOR}: {(c['ceiling'] < CEIL_FLOOR).sum()}; "
          f"fraction median {c['fraction'].median():.3f}")  # fmt: skip
    link, meta = modules()
    link = link[link["ko"].isin(c.index)]
    out = []
    for m, g in link.groupby("module"):
        k = c.loc[g["ko"].unique()]
        if len(k) < MIN_KOS or m not in meta.index:
            continue
        expl = (k["r2"] * k["ss_nontech"]).sum()
        repro = (k["ss_nontech"] - k["within_tile"]).clip(lower=0).sum()
        out.append({"module": m, **meta.loc[m].to_dict(), "kos": len(k), "median_fraction": k["fraction"].median(),
                    "q25": k["fraction"].quantile(0.25), "q75": k["fraction"].quantile(0.75), "median_r2": k["r2"].median(),
                    "median_ceiling": k["ceiling"].median(), "pooled_fraction": expl / repro,
                    "median_fraction_temp": k["fraction_temp"].median()})  # fmt: skip
    t = pd.DataFrame(out).sort_values("median_fraction", ascending=False)
    t.round(4).to_csv(mm.OUT / "module_predictability.csv", index=False)
    pd.set_option("display.width", 250)
    print(f"{len(t)} modules with >= {MIN_KOS} modelled KOs")
    cls = t.groupby("class2").agg(modules=("module", "size"), median=("median_fraction", "median")).sort_values("median", ascending=False)
    print(cls.round(3).to_string())
    print(t[["module", "name", "class2", "kos", "median_fraction", "median_r2", "median_ceiling", "pooled_fraction",
             "median_fraction_temp"]].round(3).to_string(index=False))  # fmt: skip
    print(f"median over modules: Granite {t['median_fraction'].median():.3f}, temperature {t['median_fraction_temp'].median():.3f}; "
          f"KO-level median fraction Granite {c['fraction'].median():.3f}, temperature {c['fraction_temp'].median():.3f}")
    figure(c, t)


BRITE_A = ["Metabolism", "Genetic Information Processing", "Environmental Information Processing", "Cellular Processes"]


def figure(c, t):
    """a: KEGG modules (dots), one row per module class, sorted by the class median; best module named at the right.
    b: KEGG BRITE level-B categories (all KOs, incl. transport and genetic information processing, which modules do
    not cover): median R2 / ceiling with interquartile range, Granite PCs (filled) and temperature (open)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})  # fmt: skip
    BLUE, GREY, MAG = "#0072B2", "#8c8c8c", "#c0457e"
    cls = t.groupby("class2")["median_fraction"].median().sort_values()
    b = pd.read_csv(REPO / "data" / "reference" / "kegg_brite_ko00001.csv")
    b = b[b["ko"].isin(c.index) & b["A"].isin(BRITE_A + ["Not Included in Pathway or Brite"])].drop_duplicates(["ko", "B"])
    b = b[~b["B"].str.startswith("Unclassified")]
    rows_b = []
    for B, g in b.groupby("B"):
        f, ft = c.loc[g["ko"], "fraction"], c.loc[g["ko"], "fraction_temp"]
        if len(f) >= 15:
            rows_b.append({"category": B, "kos": len(f), "median": f.median(), "q25": f.quantile(.25), "q75": f.quantile(.75),
                           "median_temp": ft.median()})  # fmt: skip
    rb = pd.DataFrame(rows_b).sort_values("median")
    rb.round(4).to_csv(mm.OUT / "category_predictability.csv", index=False)
    na, nb = len(cls), len(rb)
    h = (na + nb) * 0.16 + 1.4
    fig = plt.figure(figsize=(7.2, h))
    ya0 = (nb * 0.16 + 0.95) / h
    a = fig.add_axes([0.30, ya0, 0.38, (na * 0.16) / h])
    bx = fig.add_axes([0.30, 0.45 / h, 0.38, (nb * 0.16) / h])
    lo, hi = -0.2, 0.8
    for i, k in enumerate(cls.index):
        g = t[t["class2"] == k]
        a.scatter(g["median_fraction"].clip(lo, hi), np.full(len(g), i) + np.linspace(-0.18, 0.18, len(g)) * (len(g) > 1),
                  s=np.clip(g["kos"], 4, 30) * 0.9, color=BLUE, alpha=0.55, lw=0, zorder=3)  # fmt: skip
        a.plot([cls[k]] * 2, [i - 0.35, i + 0.35], color="black", lw=1.2, zorder=4)
        best = g.iloc[g["median_fraction"].argmax()]
        nm = best["name"].split(",")[0]
        a.text(hi + 0.02, i, (nm[:46] + "…") if len(nm) > 47 else nm, va="center", fontsize=5.6, color="#444444", clip_on=False)
    a.set_yticks(range(na))
    a.set_yticklabels([f"{k} ({(t['class2'] == k).sum()})" for k in cls.index], fontsize=6)
    a.set_ylim(-0.6, na - 0.4)
    a.set_xlim(lo, hi)
    a.axvline(0, color="#cccccc", lw=0.6, zorder=1)
    a.text(hi + 0.02, na - 0.2, "Best module in the class", fontsize=6, fontweight="bold", va="bottom")
    a.tick_params(axis="x", labelbottom=False)
    a.text(-0.72, 1.0, "a", transform=a.transAxes, fontsize=9, fontweight="bold", va="bottom")
    y = np.arange(nb)
    bx.hlines(y, rb["q25"].clip(lo, hi), rb["q75"].clip(lo, hi), color=BLUE, alpha=0.35, lw=3, zorder=2)
    bx.scatter(rb["median"].clip(lo, hi), y, s=16, color=BLUE, zorder=3, lw=0, label="Granite PCs (median, IQR)")
    bx.scatter(rb["median_temp"].clip(lo, hi), y, s=14, facecolor="white", edgecolor=MAG, lw=0.8, zorder=3,
               label="In-situ temperature (median)")  # fmt: skip
    bx.set_yticks(y)
    bx.set_yticklabels([f"{k} ({n})" for k, n in zip(rb["category"], rb["kos"])], fontsize=6)
    bx.set_ylim(-0.6, nb - 0.4)
    bx.set_xlim(lo, hi)
    bx.axvline(0, color="#cccccc", lw=0.6, zorder=1)
    bx.set_xlabel("Held-out R² as a share of the predictable variance (R² / within-tile ceiling)")
    bx.legend(frameon=False, fontsize=6, loc="lower left", bbox_to_anchor=(1.02, 0.0), handletextpad=0.2)
    bx.text(-0.72, 1.0, "b", transform=bx.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for suf in ("png", "pdf"):
        fig.savefig(REPO / "figures" / "output" / f"module_predictability.{suf}", dpi=300)
    plt.close(fig)
    print(rb.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
