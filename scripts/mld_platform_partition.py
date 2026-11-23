#!/usr/bin/env python3
"""Within-platform check of the 30-day MLD gain on C, and the 100 % variance partition per target (RESULTS.md R56).

Part "mld" (C, 1,215 R43 tiles, one cruise held out, ridge with per-target GCV alpha as R51/R53):
  a. raw_sat vs raw_sat + MLD30 spline, with and without platform (NovaSeq 6000 cruise group, R43) as a covariate.
  b. The same comparison with the cruise hold-out run inside each platform group only (train and test on one
     platform), so MLD cannot stand in for the platform difference.
  c. How MLD30 differs between the platform groups (tile means, per-cruise means).
  d. Within-cruise Spearman of MLD30 with the raw_sat C residual (centred per cruise), all tiles and per group.
  e. Per-cruise squared-error change from adding MLD30, with the cruise's platform.
Part "partition": each target's total variance split into parts that sum to 100 %:
  technical      in-sample OLS on platform + cubic spline (8 knots) of log10 KO reads (R54 covariates)
  explained      held-out (cruise) skill on the technically adjusted targets: for C, raw_sat + MLD30; for KO
                 targets, the R54 plug-in linear two-stage (32 Granite PCs + standardised C; observed C on
                 the training tiles, C predicted from raw_sat + MLD30 for the test cruise)
  cruise offsets per-cruise means of the held-out residual (not explained)
  replicate noise within-tile replicate variance / replicates per tile (as the R51 budget), capped at what is left
  within-cruise  the rest (biological, reproducible, not explained)
The technical fit here uses all tiles (a description: 10 coefficients per KO); the strict fold-wise number,
technical fit on the training cruises only, is R54's and is printed beside it.
Headlines: explained / (100 - technical) and explained / (100 - technical - replicate noise).

Outputs: data/results/mld_platform/{mld_check.csv, mld_per_cruise.csv, partition.csv, summary.txt};
figures/output/variance_partition.png.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/mld_platform_partition.py --parts mld partition figure
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import SplineTransformer, StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from domain_split_models import EMB, FSETS, NOVASEQ_6000  # noqa: E402
from tune_figure_kos import build, load, ridge_gcv  # noqa: E402
from two_stage_models import TARGETS, spl  # noqa: E402
from water_mass_test import features  # noqa: E402

OUT = REPO / "data" / "results" / "mld_platform"
FIG = REPO / "figures" / "output"
TD = REPO / "data" / "targets" / "domain_split"
LOG: list[str] = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def setup(tname, fs, ks):
    path = TARGETS[tname]
    kos = pd.read_csv(path, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, sat, y, cr, unit, rows = load(argparse.Namespace(targets=path, embeddings_dirs=EMB), kos, return_rows=True)
    tile = lambda s: s.reindex(rows["sample_id"]).set_axis(rows.index).groupby(rows["tile_key"]).mean().reindex(unit.index)  # noqa: E731
    c = tile(pd.read_csv(TD / "C.csv", index_col="sample_key")["euk_prok_log_ratio"]).to_numpy()[:, None]
    mld = tile(fs["mld_mean30"])
    mld = mld.fillna(mld.median()).to_numpy()[:, None]
    logtot = tile(np.log10(ks["total"])).to_numpy()
    return dict(emb=emb, temp=temp, sat=sat, y=y, cr=cr, unit=unit, rows=rows, c=c, mld=mld, logtot=logtot, kos=kos,
                plat=np.isin(cr, NOVASEQ_6000))  # fmt: skip


def predict(d, tr, te, y, use_mld, use_plat, blocks=("S~",)):
    xtr, xte = build(list(blocks), tr, te, d["emb"], d["temp"], d["sat"], {})
    if use_mld:
        a, b = spl(d["mld"][tr], d["mld"][te])
        xtr, xte = np.hstack([xtr, a]), np.hstack([xte, b])
    if use_plat:
        p = d["plat"].astype(float)[:, None]
        xtr, xte = np.hstack([xtr, p[tr]]), np.hstack([xte, p[te]])
    return np.asarray(ridge_gcv(xtr, y[tr], xte)).reshape(len(te), -1)


def loco(d, y, idx, **kw):
    pred = np.full_like(y, np.nan)
    cr = d["cr"]
    for k in np.unique(cr[idx]):
        te, tr = idx[cr[idx] == k], idx[cr[idx] != k]
        pred[te] = predict(d, tr, te, y, **kw)
    return pred


def expl(y, p, idx):
    yc = y[idx] - y[idx].mean(0)
    return 1 - ((y[idx] - p[idx]) ** 2).sum() / (yc**2).sum()


def mld_part(fs, ks):
    d = setup("C", fs, ks)
    y, cr, plat, mld = d["c"], d["cr"], d["plat"], d["mld"][:, 0]
    allidx = np.arange(len(y))
    groups = {"all tiles": allidx, "NovaSeq 6000 cruises": np.where(plat)[0], "other cruises": np.where(~plat)[0]}
    say(f"C: {len(y):,} tiles, {len(np.unique(cr))} cruises; NovaSeq 6000 group {plat.sum()} tiles in "
        f"{len(np.unique(cr[plat]))} cruises, other {(~plat).sum()} tiles in {len(np.unique(cr[~plat]))} cruises")  # fmt: skip
    rows, preds = [], {}
    say("\na/b. Variance of C explained, cruise held out (scored on the tiles of each group):")
    for scheme in ("pooled training", "within-platform training"):
        for usep in (False, True):
            if scheme == "within-platform training" and usep:
                continue
            for usem in (False, True):
                name = f"raw_sat{' + MLD30' if usem else ''}{' + platform' if usep else ''}"
                if scheme == "pooled training":
                    p = loco(d, y, allidx, use_mld=usem, use_plat=usep)
                else:
                    p = np.full_like(y, np.nan)
                    for g in (np.where(plat)[0], np.where(~plat)[0]):
                        p[g] = loco(d, y, g, use_mld=usem, use_plat=False)[g]
                preds[(scheme, name)] = p
                for gname, idx in groups.items():
                    rows.append({"scheme": scheme, "model": name, "scored_on": gname, "tiles": len(idx), "explained": expl(y, p, idx)})
    t = pd.DataFrame(rows)
    piv = t.pivot_table(index=["scheme", "model"], columns="scored_on", values="explained", sort=False)[list(groups)]
    say((piv * 100).round(2).to_string())
    say("\n   MLD30 gain (percentage points):")
    for scheme, base, plus in [("pooled training", "raw_sat", "raw_sat + MLD30"),
                               ("pooled training", "raw_sat + platform", "raw_sat + MLD30 + platform"),
                               ("within-platform training", "raw_sat", "raw_sat + MLD30")]:  # fmt: skip
        g = (piv.loc[(scheme, plus)] - piv.loc[(scheme, base)]) * 100
        say(f"   {scheme:25s} {base:22s} -> +MLD30: " + "  ".join(f"{k} {v:+.2f}" for k, v in g.items()))
    t.round(5).to_csv(OUT / "mld_check.csv", index=False)

    say("\nc. MLD30 (log10 m, along-track 30-day mean) by platform group:")
    say(f"   tiles: NovaSeq 6000 mean {mld[plat].mean():.3f} (SD {mld[plat].std():.3f}); other {mld[~plat].mean():.3f} "
        f"(SD {mld[~plat].std():.3f}); point-biserial r {np.corrcoef(mld, plat)[0, 1]:+.3f}")  # fmt: skip
    cm = pd.DataFrame({"cruise": cr, "plat": plat, "mld": mld, "C": y[:, 0]}).groupby("cruise").agg(
        platform=("plat", "first"), tiles=("mld", "size"), mld30=("mld", "mean"), mld30_sd=("mld", "std"), C=("C", "mean"))
    share = 1 - ((mld - pd.Series(mld).groupby(cr).transform("mean").to_numpy()) ** 2).sum() / ((mld - mld.mean()) ** 2).sum()
    shp = 1 - ((mld - pd.Series(mld).groupby(plat).transform("mean").to_numpy()) ** 2).sum() / ((mld - mld.mean()) ** 2).sum()
    say(f"   share of MLD30 variance between cruises {share:.1%}; between platform groups {shp:.1%}")

    say("\nd. Within-cruise Spearman of MLD30 with the raw_sat C residual (both centred per cruise):")
    res = y[:, 0] - preds[("pooled training", "raw_sat")][:, 0]
    cen = lambda v, idx: v[idx] - pd.Series(v[idx]).groupby(cr[idx]).transform("mean").to_numpy()  # noqa: E731
    for gname, idx in groups.items():
        r, pv = spearmanr(cen(mld, idx), cen(res, idx))
        r2, pv2 = spearmanr(cen(mld, idx), cen(y[:, 0], idx))
        say(f"   {gname:22s} residual rho {r:+.3f} (p {pv:.1e});  C itself rho {r2:+.3f} (p {pv2:.1e})")

    say("\ne. Per cruise, pooled training: squared error of C, raw_sat vs raw_sat + MLD30:")
    p0, p1 = preds[("pooled training", "raw_sat")][:, 0], preds[("pooled training", "raw_sat + MLD30")][:, 0]
    pw0, pw1 = preds[("within-platform training", "raw_sat")][:, 0], preds[("within-platform training", "raw_sat + MLD30")][:, 0]
    pc = pd.DataFrame({"cruise": cr, "e0": (y[:, 0] - p0) ** 2, "e1": (y[:, 0] - p1) ** 2, "w0": (y[:, 0] - pw0) ** 2,
                       "w1": (y[:, 0] - pw1) ** 2, "b0": y[:, 0] - p0, "b1": y[:, 0] - p1}).groupby("cruise").mean()  # fmt: skip
    pc = cm.join(pc)
    pc["change_pooled"] = pc["e1"] / pc["e0"] - 1
    pc["change_within_platform"] = pc["w1"] / pc["w0"] - 1
    pc = pc.sort_values(["platform", "change_pooled"])
    say(pc[["platform", "tiles", "mld30", "C", "b0", "b1", "e0", "e1", "change_pooled", "change_within_platform"]].round(3).to_string())
    for lab, m in (("NovaSeq 6000", pc["platform"]), ("other", ~pc["platform"])):
        say(f"   {lab}: error lower with MLD30 in {(pc.loc[m, 'change_pooled'] < 0).sum()} of {m.sum()} cruises (pooled), "
            f"{(pc.loc[m, 'change_within_platform'] < 0).sum()} of {m.sum()} (within-platform)")  # fmt: skip
    pc.round(5).to_csv(OUT / "mld_per_cruise.csv")


def halves_part(fs, ks, n=20):
    """Control for b: the within-group gain when the 13 cruises are split into random 6 / 7 groups instead of by
    platform. Fewer training cruises alone may erase a between-cruise effect such as MLD30."""
    d = setup("C", fs, ks)
    y, cr = d["c"], d["cr"]
    cruises = np.unique(cr)
    rng = np.random.default_rng(0)
    out = []
    for i in range(n):
        g6 = rng.choice(cruises, 6, replace=False)
        same = sorted(g6) == sorted(np.unique(cr[d["plat"]]))
        mask = np.isin(cr, g6)
        r = {"split": i, "is_platform_split": same}
        for usem in (False, True):
            p = np.full_like(y, np.nan)
            for g in (np.where(mask)[0], np.where(~mask)[0]):
                p[g] = loco(d, y, g, use_mld=usem, use_plat=False)[g]
            r["mld" if usem else "base"] = expl(y, p, np.arange(len(y)))
        r["gain"] = r["mld"] - r["base"]
        out.append(r)
        say(f"   random split {i:2d}: raw_sat {r['base']:.2%}  + MLD30 {r['mld']:.2%}  gain {r['gain'] * 100:+.2f} pp")
    t = pd.DataFrame(out)
    t.round(5).to_csv(OUT / "mld_random_halves.csv", index=False)
    say(f"   random 6/7 cruise splits (n {n}): MLD30 gain median {t['gain'].median() * 100:+.2f} pp, "
        f"10-90% {t['gain'].quantile(0.1) * 100:+.2f} to {t['gain'].quantile(0.9) * 100:+.2f}; "
        f"gain > 0 in {(t['gain'] > 0).sum()} of {n}; platform split +0.53 pp")  # fmt: skip


def replicate_noise(y_tile, d, sample_targets):
    rows, unit = d["rows"], d["unit"]
    multi = rows.groupby("tile_key")["sample_id"].agg(list)
    diffs = []
    for samples in multi[multi.str.len() >= 2]:
        v = sample_targets.reindex([x for x in samples if x in sample_targets.index]).to_numpy()
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                diffs.append(v[i] - v[j])
    sigma2 = (np.array(diffs) ** 2).mean(0) / 2
    k = multi.reindex(unit.index).str.len().fillna(1).to_numpy()
    return sigma2.sum() * (1 / k).sum(), len(diffs)


def partition_part(fs, ks):
    out = []
    strict = pd.read_csv(REPO / "data" / "results" / "two_stage" / "plugin_headline.csv")
    for tname in ["clr_all", "C", "D", "A", "Dp", "B", "De"]:
        d = setup(tname, fs, ks)
        y = d["c"] if tname == "C" else d["y"]
        cr, plat, logtot = d["cr"], d["plat"], d["logtot"]
        yc = y - y.mean(0)
        tot = (yc**2).sum()
        X = np.hstack([np.ones((len(y), 1)), plat[:, None].astype(float), SplineTransformer(n_knots=8, degree=3).fit_transform(logtot[:, None])])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        fit = X @ beta
        tech = ((fit - fit.mean(0)) ** 2).sum()
        ya = y - (fit - y.mean(0))
        pred = np.full_like(ya, np.nan)
        pg = np.full_like(ya, np.nan)
        for k in np.unique(cr):
            te, tr = np.where(cr == k)[0], np.where(cr != k)[0]
            c1 = predict(d, tr, te, d["c"], use_mld=True, use_plat=False)
            if tname == "C":
                pred[te] = c1
                continue
            xtr, xte = build(FSETS["granite"], tr, te, d["emb"], d["temp"], d["sat"], {})
            pg[te] = np.asarray(ridge_gcv(xtr, ya[tr], xte)).reshape(len(te), -1)
            sc = StandardScaler().fit(d["c"][tr])
            pred[te] = np.asarray(ridge_gcv(np.hstack([xtr, sc.transform(d["c"][tr])]), ya[tr],
                                            np.hstack([xte, sc.transform(c1)]))).reshape(len(te), -1)  # fmt: skip
        yac = ya - ya.mean(0)
        res = ya - pred
        explained = (yac**2).sum() - (res**2).sum()
        res = res - res.mean(0)
        cm = pd.DataFrame(res).groupby(cr).transform("mean").to_numpy()
        offs = (cm**2).sum()
        within = ((res - cm) ** 2).sum()
        st = pd.read_csv(TARGETS[tname], index_col="sample_key")[d["kos"]]
        noise, npairs = replicate_noise(y, d, st)
        noise = min(noise, within)
        r = {"target": tname, "tiles": len(y), "kos": y.shape[1], "technical": tech / tot, "explained": explained / tot,
             "cruise_offsets": offs / tot, "replicate_noise": noise / tot, "within_cruise": (within - noise) / tot,
             "replicate_pairs": npairs}  # fmt: skip
        if tname != "C":
            r["granite_alone_explained"] = 1 - ((ya - pg) ** 2).sum() / (yac**2).sum()
        r["sum"] = r["technical"] + r["explained"] + r["cruise_offsets"] + r["replicate_noise"] + r["within_cruise"]
        r["of_nontechnical"] = r["explained"] / (1 - r["technical"])
        r["of_reproducible_nontechnical"] = r["explained"] / (1 - r["technical"] - r["replicate_noise"])
        m = "two_stage" if tname != "C" else "granite"
        s = strict[(strict["target"] == tname) & (strict["mode"] == "tech") & (strict["model"] == m)]
        r["strict_foldwise_of_nontechnical"] = s["explained_raw"].iloc[0] if (len(s) and tname != "C") else np.nan
        out.append(r)
        say(f"{tname:8s} tech {r['technical']:6.2%}  explained {r['explained']:6.2%}  cruise offsets {r['cruise_offsets']:6.2%}  "
            f"noise {r['replicate_noise']:6.2%}  within-cruise {r['within_cruise']:6.2%}  (sum {r['sum']:.4f})  -> "
            f"{r['of_nontechnical']:6.2%} of non-technical, {r['of_reproducible_nontechnical']:6.2%} of reproducible non-technical")  # fmt: skip
    p = pd.DataFrame(out)
    p.round(5).to_csv(OUT / "partition.csv", index=False)
    say("\n" + (p.set_index("target").drop(columns=["tiles", "kos", "replicate_pairs"]) * 100).round(2).to_string())


def figure():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    p = pd.read_csv(OUT / "partition.csv").set_index("target")
    tg = ["clr_all", "C", "D", "A", "Dp", "B", "De"]
    tl = {"clr_all": "All KOs\n(CLR)", "C": "C\nEuk : prok", "D": "D\nShared", "A": "A\nProk-specific", "Dp": "Dp\nShared, prok",
          "B": "B\nEuk-specific", "De": "De\nShared, euk"}  # fmt: skip
    parts = [("explained", "Explained by satellite (cruise held out)", "#E69F00"),
             ("cruise_offsets", "Cruise offsets, not explained", "#56B4E9"),
             ("within_cruise", "Within-cruise, not explained", "#0072B2"),
             ("replicate_noise", "Replicate noise", "#BBBBBB"),
             ("technical", "Technical (platform, depth)", "#666666")]  # fmt: skip
    fig, ax = plt.subplots(figsize=(7.2, 3.3))
    x = np.arange(len(tg))
    bottom = np.zeros(len(tg))
    for col, lab, colr in parts:
        v = p.loc[tg, col].clip(lower=0).to_numpy() * 100
        ax.bar(x, v, 0.62, bottom=bottom, color=colr, label=lab, edgecolor="white", linewidth=0.8)
        if col == "explained":
            for xi, vi in zip(x, v):
                ax.text(xi, vi / 2, f"{vi:.0f}", ha="center", va="center", fontsize=7, color="black")
        bottom += v
    for xi, t in zip(x, tg):
        ax.text(xi, 101.5, f"{p.loc[t, 'of_nontechnical'] * 100:.0f}% | {p.loc[t, 'of_reproducible_nontechnical'] * 100:.0f}%",
                ha="center", va="bottom", fontsize=6.5)  # fmt: skip
    ax.set_xticks(x)
    ax.set_xticklabels([tl[t] for t in tg], fontsize=7)
    ax.set_ylim(0, 108)
    ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel("Share of total variance (%)")
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[::-1], l[::-1], frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig(FIG / "variance_partition.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", default=["mld", "partition", "figure"], choices=["mld", "halves", "partition", "figure"])
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    fs = features() if {"mld", "halves", "partition"} & set(a.parts) else None
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    if "mld" in a.parts:
        mld_part(fs, ks)
    if "halves" in a.parts:
        say("\nControl: within-group training on random 6 / 7 cruise splits:")
        halves_part(fs, ks)
    if "partition" in a.parts:
        say("\nVariance partition (share of total variance; technical fit in-sample, model cruise held out):")
        partition_part(fs, ks)
    if "figure" in a.parts:
        figure()
    if LOG:
        (OUT / f"summary_{'_'.join(p for p in a.parts if p != 'figure')}.txt").write_text("\n".join(LOG) + "\n")
