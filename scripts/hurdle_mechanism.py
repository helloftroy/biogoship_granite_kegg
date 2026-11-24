#!/usr/bin/env python3
"""R66: why are hurdle KOs "satellite" or "technical"? Test the detection-limit hypothesis.

Hypothesis: technical-class KOs are low-abundance genes near the detection limit. Their absences are sampling
dropouts, and when present their counts are small, so the abundance-where-present is mostly counting noise.
Satellite-class KOs are mid-prevalence genes whose presence and abundance follow real biogeography.

Per KO (the 5,808 hurdle KOs; samples and tiles as in hurdle_kos.py build, G = 1e6):
  expected_count_G   mean share over kept samples x G (expected count in a library of G genes)
  count_present_med  median count in present samples
  noise_frac         Poisson (counting) variance of the tile abundance-where-present / its observed variance.
                     Per present sample var(log c) ~ 1/c; the tile value is a mean over n present samples,
                     so var = sum(1/c_i) / n^2. noise_frac ~ 1 means nothing beyond counting noise to predict.
  r_depth_within     correlation of tile presence with log10 depth after removing cruise means
  r_temp_within      same, with in-situ temperature
  eta2_cruise        share of tile-presence variance between cruises
  length_mean        mean reference gene length (nt; R44 ko_gc.csv), euk_share = eukaryote share of its genes

Outputs: data/results/hurdle/mechanism.csv, mechanism_summary.txt, figures/output/hurdle_mechanism.{png,pdf}
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "results" / "hurdle"
BX = REPO / "data" / "results" / "benchmark_export"
FIG = REPO / "figures" / "output"
G = 1e6
LOG: list = []


def say(s=""):
    print(s)
    LOG.append(str(s))


def per_ko():
    f = pd.read_csv(OUT / "hurdle_all.csv").set_index("ko")
    kos = f.index
    S = pd.read_pickle(OUT / "sample_counts.pkl")[kos].astype(np.float64)
    T = pd.read_pickle(OUT / "sample_counts.pkl").astype(np.float64).sum(1)
    keep = T >= G
    S, T = S[keep], T[keep]
    P = S.div(T, axis=0)
    pres = P >= 1 / G
    tiles = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    tp = pd.read_csv(OUT / "presence_tile.csv.gz", index_col=0)[kos]
    ta = pd.read_csv(OUT / "abundance_tile.csv.gz", index_col=0)[kos]
    members = tiles.loc[tp.index, "sample_keys"].str.split(";").explode()
    members = members[members.isin(S.index)]
    tile_of = members.index.values
    Sm, Pm = S.loc[members.values], pres.loc[members.values]
    inv = (1 / Sm.where(Pm)).fillna(0)
    inv.index = tile_of
    npres = Pm.astype(float)
    npres.index = tile_of
    noise_tile = inv.groupby(level=0).sum() / npres.groupby(level=0).sum().pow(2)  # var of mean log over present samples
    noise_tile = noise_tile.reindex(tp.index).where(tp == 1)
    ta_p = ta.where(tp == 1)
    r = pd.DataFrame(index=kos)
    r["expected_count_G"] = P.mean() * G
    r["count_present_med"] = S.where(pres).median()
    r["noise_var"] = noise_tile.mean()
    r["obs_var"] = ta_p.var()
    r["noise_frac"] = r["noise_var"] / r["obs_var"]
    # within-cruise correlations of tile presence
    t = tiles.loc[tp.index]
    cr = t["cruise"].values
    depth = pd.read_csv(OUT / "tile_depth.csv", index_col=0)["log10_total_genes"].loc[tp.index]
    Y = tp.astype(float)
    Yd = Y - Y.groupby(cr).transform("mean")

    def r_within(x):
        x = pd.Series(x, index=tp.index)
        ok = x.notna().values
        xd = (x - x.groupby(cr).transform("mean"))[ok].values
        yd = Yd[ok].values
        num = (yd * xd[:, None]).sum(0)
        den = np.sqrt((yd**2).sum(0) * (xd**2).sum())
        return np.where(den > 0, num / den, np.nan)

    r["r_depth_within"] = r_within(depth.values)
    r["r_temp_within"] = r_within(t["temperature_c"].values)
    # Eukaryote : prokaryote balance (target C). Observed: within-cruise correlation with presence.
    # Predicted (Granite, cruise held out; R60): within-cruise AUC of +/- predicted C, the sign taken from the
    # training cruises' correlation, so this is a one-number model comparable with auc_within_granite.
    Cobs = pd.read_csv(BX / "targets" / "C.csv.gz", index_col="tile_key")["euk_prok_log_ratio"]
    z = np.load(REPO / "data/results/main_model/predictions_C_granite_cv.npz", allow_pickle=True)
    Cpred = pd.Series(z["predicted"][:, 0], index=z["tile"]).loc[tp.index].values
    r["r_C_within"] = r_within(Cobs.reindex(tp.index).values)
    Yv = Y.values
    aucs = np.zeros(Yv.shape[1])
    wsum = np.zeros(Yv.shape[1])
    from sklearn.metrics import roc_auc_score

    for c in np.unique(cr):
        te, trn = cr == c, cr != c
        xc = Cpred[trn] - Cpred[trn].mean()
        sign = np.sign(((Yv[trn] - Yv[trn].mean(0)) * xc[:, None]).sum(0))
        for j in range(Yv.shape[1]):
            yj = Yv[te, j]
            if 0 < yj.mean() < 1:
                aucs[j] += roc_auc_score(yj, sign[j] * Cpred[te]) * te.sum()
                wsum[j] += te.sum()
    r["auc_within_Cpred"] = np.where(wsum > 0, aucs / np.maximum(wsum, 1), np.nan)
    tot = ((Y - Y.mean()) ** 2).sum()
    r["eta2_cruise"] = (((Y.groupby(cr).transform("mean") - Y.mean()) ** 2).sum() / tot).values
    gc = pd.read_csv(REPO / "transfer/results/ko_gc/ko_gc.csv", index_col="ko")
    r["length_mean"] = gc["length_mean"].reindex(kos)
    r["euk_share"] = (gc["genes_Eukaryota"] / gc["genes"]).reindex(kos)
    cols = ["prevalence", "presence_class", "in_core", "auc_within_tech", "auc_within_granite", "auc_within_both",
            "auc_tech", "auc_granite", "auc_both", "abund_r2", "abund_technical_share", "name"]  # fmt: skip
    r = f[cols].join(r)
    r.round(5).to_csv(OUT / "mechanism.csv")
    return r


CLS = ["satellite", "weak", "technical"]
PBINS = [0.2, 0.5, 0.7, 0.9]


def summarise(r):
    r = r.copy()
    r["pbin"] = pd.cut(r["prevalence"], PBINS, include_lowest=True)
    r["log10_expected"] = np.log10(r["expected_count_G"])
    say("## Q1. Does adding Granite to the technical covariates improve presence AUC?")
    for sc in ("auc_within", "auc"):
        d = r[f"{sc}_both"] - r[f"{sc}_tech"]
        say(f"  {sc}: both - technical, median {d.median():+.3f} (10-90% {d.quantile(.1):+.3f} to {d.quantile(.9):+.3f}); "
            f">= +0.02 for {(d >= .02).sum():,}, <= -0.02 for {(d <= -.02).sum():,} of {len(d):,}")  # fmt: skip
        for c in CLS:
            dc = d[r["presence_class"] == c]
            say(f"    {c:9s} median {dc.median():+.3f}  (both {r.loc[dc.index, sc + '_both'].median():.3f} vs technical "
                f"{r.loc[dc.index, sc + '_tech'].median():.3f}; >= +0.02 for {(dc >= .02).mean():.0%})")  # fmt: skip
    say("\n## Per-class medians (all prevalences)")
    vars_ = ["prevalence", "expected_count_G", "count_present_med", "noise_frac", "abund_r2", "r_depth_within",
             "r_temp_within", "eta2_cruise", "length_mean", "euk_share"]  # fmt: skip
    say(r.groupby("presence_class")[vars_].median().reindex(CLS).round(3).T.to_string())
    say("\n## Per-class medians within prevalence bins (class is confounded with prevalence)")
    for v in ("expected_count_G", "noise_frac", "abund_r2", "r_depth_within", "r_temp_within", "length_mean", "euk_share"):
        tab = r.pivot_table(index="pbin", columns="presence_class", values=v, aggfunc="median", observed=True).reindex(columns=CLS)
        say(f"  {v}:\n" + tab.round(3).to_string())
    n = r.pivot_table(index="pbin", columns="presence_class", values="abund_r2", aggfunc="size", observed=True).reindex(columns=CLS)
    say("  n KOs:\n" + n.to_string())
    say("\n## Is low abundance R2 explained by counting noise?")
    ok = r[["noise_frac", "abund_r2"]].notna().all(1)
    rho = spearmanr(r.loc[ok, "noise_frac"], r.loc[ok, "abund_r2"])[0]
    say(f"  Spearman(noise_frac, abund_r2) = {rho:.3f} (n {ok.sum():,})")
    say(f"  noise_frac quantiles: " + ", ".join(f"{q:.0%} {r['noise_frac'].quantile(q):.3f}" for q in (.1, .25, .5, .75, .9)))
    say(f"  abund_r2 above its noise ceiling (1 - noise_frac): {(r['abund_r2'] > 1 - r['noise_frac']).sum()} KOs")
    r["nbin"] = pd.qcut(r["noise_frac"], 5)
    tab = r.pivot_table(index="nbin", columns="presence_class", values="abund_r2", aggfunc="median", observed=True).reindex(columns=CLS)
    cnt = r.pivot_table(index="nbin", columns="presence_class", values="abund_r2", aggfunc="size", observed=True).reindex(columns=CLS)
    say("  median abund_r2 by noise_frac quintile x class:\n" + tab.round(3).to_string() + "\n  n:\n" + cnt.to_string())
    say("\n## Is presence class explained by abundance / gene length / domain? (satellite vs technical, ranks)")
    st = r[r["presence_class"].isin(["satellite", "technical"])].dropna(subset=["length_mean", "euk_share"])
    y = (st["presence_class"] == "satellite").astype(int).values
    Xn = ["log10_expected", "length_mean", "euk_share", "prevalence"]
    X = st[Xn].rank().values
    X = (X - X.mean(0)) / X.std(0)
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import cross_val_predict

    for cols in (["prevalence"], ["log10_expected"], ["length_mean"], ["euk_share"], ["log10_expected", "prevalence"], Xn):
        idx = [Xn.index(c) for c in cols]
        p = cross_val_predict(LogisticRegression(C=1e3), X[:, idx], y, cv=10, method="predict_proba")[:, 1]
        m = LogisticRegression(C=1e3).fit(X[:, idx], y)
        say(f"  {' + '.join(cols):55s} 10-fold AUC {roc_auc_score(y, p):.3f}; std. coef "
            + ", ".join(f"{c} {b:+.2f}" for c, b in zip(cols, m.coef_[0])))  # fmt: skip
    say(f"  n = {len(st):,} ({y.sum():,} satellite)")
    say("\n## Is the satellite signal the eukaryote : prokaryote balance (target C)?")
    r["euk_bin"] = pd.cut(r["euk_share"], [-0.01, 0.1, 0.5, 0.9, 1.0])
    say("  class by eukaryote share of the KO's reference genes:\n"
        + pd.crosstab(r["euk_bin"], r["presence_class"]).reindex(columns=CLS).to_string())  # fmt: skip
    say("  medians by class:\n" + r.groupby("presence_class")[["r_C_within", "auc_within_Cpred", "auc_within_granite"]]
        .median().reindex(CLS).round(3).to_string())  # fmt: skip
    gap = r["auc_within_granite"] - r["auc_within_Cpred"]
    for c in CLS:
        g = gap[r["presence_class"] == c]
        say(f"  {c:9s} Granite - predicted-C within-cruise AUC: median {g.median():+.3f}; Granite ahead by >= 0.05 for "
            f"{(g >= .05).mean():.0%}")  # fmt: skip
    for lo, hi in ((-0.01, 0.1), (0.9, 1.0)):
        m = (r["euk_share"] > lo) & (r["euk_share"] <= hi)
        g = r[m]
        say(f"  euk_share ({lo:.1f}, {hi:.1f}]: n {m.sum():,}; median Granite {g['auc_within_granite'].median():.3f}, "
            f"predicted C {g['auc_within_Cpred'].median():.3f}, technical {g['auc_within_tech'].median():.3f}; "
            f"abund R2 {g['abund_r2'].median():.3f}; r_C_within {g['r_C_within'].median():+.3f}")  # fmt: skip
    say("  prokaryote-dominated KOs (euk_share <= 0.1), class by prevalence bin:\n"
        + pd.crosstab(r.loc[r["euk_share"] <= 0.1, "pbin"], r.loc[r["euk_share"] <= 0.1, "presence_class"]).reindex(columns=CLS).to_string())  # fmt: skip
    pk = r[(r["euk_share"] <= 0.1) & (r["presence_class"] == "satellite")].sort_values("auc_within_granite", ascending=False)
    say(f"  prokaryote-dominated satellite KOs: {len(pk)}; top 25:")
    say(pk[["name", "prevalence", "auc_within_tech", "auc_within_granite", "auc_within_Cpred", "r_temp_within", "abund_r2"]]
        .head(25).assign(name=lambda x: x["name"].str[:60]).round(3).to_string())  # fmt: skip
    say("\n## Split by domain (euk_share > 0.9 vs <= 0.1): length, depth, abundance")
    for nm, m in (("prokaryotic (<= 0.1)", r["euk_share"] <= 0.1), ("eukaryotic (> 0.9)", r["euk_share"] > 0.9)):
        g = r[m]
        say(f"  {nm}: n {len(g):,}; median length {g['length_mean'].median():.0f} nt; median r_depth {g['r_depth_within'].median():+.3f}; "
            f"rho(length, r_depth) {spearmanr(g['length_mean'], g['r_depth_within'])[0]:+.3f}; "
            f"rho(length, auc_within_tech) {spearmanr(g['length_mean'], g['auc_within_tech'])[0]:+.3f}")  # fmt: skip
        for cl in ("satellite", "weak", "technical"):
            h = g[g["presence_class"] == cl]
            say(f"    {cl:9s} n {len(h):5,}  length {h['length_mean'].median():5.0f}  r_depth {h['r_depth_within'].median():+.3f}  "
                f"expected count {h['expected_count_G'].median():5.2f}  noise_frac {h['noise_frac'].median():.3f}  "
                f"abund_r2 {h['abund_r2'].median():+.3f}")  # fmt: skip
        if nm.startswith("prok"):
            say("    rho(length, r_depth) within prevalence bins: " + ", ".join(
                f"{b}: {spearmanr(h['length_mean'], h['r_depth_within'])[0]:+.3f} (n {len(h)})"
                for b, h in g.groupby("pbin", observed=True)))  # fmt: skip
    t = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    dep = pd.read_csv(OUT / "tile_depth.csv", index_col=0).iloc[:, 0]
    Cc = pd.read_csv(BX / "targets" / "C.csv.gz", index_col="tile_key").iloc[:, 0].loc[dep.index]
    crr = t.loc[dep.index, "cruise"]
    rr = np.corrcoef(dep - dep.groupby(crr).transform("mean"), Cc - Cc.groupby(crr).transform("mean"))[0, 1]
    say(f"  within-cruise r(eukaryote : prokaryote log ratio C, log10 depth) = {rr:+.3f}")
    say("\n## Spearman across all 5,808 KOs")
    for a, b in (("length_mean", "r_depth_within"), ("length_mean", "auc_within_tech"), ("log10_expected", "auc_within_tech"),
                 ("log10_expected", "r_depth_within"), ("noise_frac", "log10_expected"), ("euk_share", "auc_within_granite"),
                 ("r_temp_within", "auc_within_granite"), ("abs_r_temp", "auc_within_granite")):  # fmt: skip
        if a == "abs_r_temp":
            r["abs_r_temp"] = r["r_temp_within"].abs()
        ok = r[[a, b]].notna().all(1)
        say(f"  {a:16s} vs {b:20s} rho {spearmanr(r.loc[ok, a], r.loc[ok, b])[0]:+.3f}")
    # length partial: within prevalence bins
    say("  length_mean vs r_depth_within, within prevalence bins: " + ", ".join(
        f"{b}: {spearmanr(g['length_mean'], g['r_depth_within'], nan_policy='omit')[0]:+.3f}"
        for b, g in r.groupby("pbin", observed=True)))  # fmt: skip
    focus = {"K02588": "nifH", "K02586": "nifD", "K02591": "nifK", "K10944": "amoA", "K10945": "amoB", "K10946": "amoC",
             "K00372": "nasA", "K06163": "phnJ", "K06164": "phnI", "K06165": "phnH", "K06166": "phnG", "K05780": "phnL",
             "K05781": "phnK", "K06162": "phnM"}  # fmt: skip
    say("\n## Focus genes")
    fr = r.loc[[k for k in focus if k in r.index]].copy()
    fr.insert(0, "gene", fr.index.map(focus))
    say(fr[["gene", "presence_class", "prevalence", "expected_count_G", "count_present_med", "noise_frac", "abund_r2",
            "r_depth_within", "r_temp_within", "length_mean", "euk_share"]].round(3).to_string())  # fmt: skip
    (OUT / "mechanism_summary.txt").write_text("\n".join(LOG) + "\n")
    return r


def figure(r):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})  # fmt: skip
    col = {"satellite": "#0072B2", "technical": "#E69F00", "weak": "#aaaaaa"}
    lab = {"satellite": "Satellite", "technical": "Technical", "weak": "Weak"}
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))
    plt.subplots_adjust(left=0.11, right=0.97, top=0.93, bottom=0.08, wspace=0.28, hspace=0.38)
    a, b, c, d = ax.ravel()
    # a: class share by eukaryote share of the KO's reference genes
    eb = pd.cut(r["euk_share"], [-0.01, 0.1, 0.5, 0.9, 1.0], labels=["\u2264 0.1", "0.1\u20130.5", "0.5\u20130.9", "> 0.9"])
    ct = pd.crosstab(eb, r["presence_class"]).reindex(columns=["satellite", "weak", "technical"])
    sh = ct.div(ct.sum(1), axis=0)
    left = np.zeros(len(sh))
    for cl in sh.columns:
        a.barh(range(len(sh)), sh[cl], left=left, color=col[cl], height=0.7, label=lab[cl])
        left += sh[cl].values
    for i_, n_ in enumerate(ct.sum(1)):
        a.text(1.01, i_, f"{n_:,}", va="center", fontsize=6)
    a.set_yticks(range(len(sh)))
    a.set_yticklabels(sh.index)
    a.set_xlim(0, 1)
    a.set_xlabel("Share of KOs in each presence class")
    a.set_ylabel("Eukaryote share of the KO's reference genes")
    a.legend(frameon=False, fontsize=6, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=1)
    # b: Granite vs predicted eukaryote : prokaryote ratio, coloured by eukaryote share
    sc = b.scatter(r["auc_within_Cpred"], r["auc_within_granite"], c=r["euk_share"], cmap="cividis", s=2, lw=0, alpha=0.7,
                   rasterized=True, vmin=0, vmax=1)  # fmt: skip
    b.plot([0.2, 1], [0.2, 1], color="#555555", lw=0.6, ls=(0, (3, 2)))
    b.set_xlim(0.2, 1)
    b.set_ylim(0.3, 1)
    b.set_xlabel("Presence AUC within cruise, predicted eukaryote : prokaryote ratio")
    b.set_ylabel("Presence AUC within cruise, Granite PCs")
    cb = fig.colorbar(sc, ax=b, fraction=0.04, pad=0.02)
    cb.set_label("Eukaryote share of genes", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    # c: median abundance R2 by counting-noise quintile, class x domain (euk share > 0.9 vs <= 0.1)
    EUK, PROK = "#6a3d9a", "#33a02c"
    nb = pd.qcut(r["noise_frac"], 5)
    mids = r.groupby(nb, observed=True)["noise_frac"].median().values
    dom = np.where(r["euk_share"] > 0.9, "euk", np.where(r["euk_share"] <= 0.1, "prok", ""))
    for dm, dc, dl in (("euk", EUK, "Eukaryotic KOs"), ("prok", PROK, "Prokaryotic KOs")):
        for cl, ls in (("satellite", "-"), ("technical", "--")):
            m = (dom == dm) & (r["presence_class"] == cl).values
            v = r[m].groupby(nb[m], observed=True)["abund_r2"].median().reindex(nb.cat.categories)
            c.plot(mids, v.values, marker="o" if cl == "satellite" else "s", ms=3, lw=1, ls=ls, color=dc,
                   label=f"{dl}, {lab[cl].lower()} ({m.sum():,})")  # fmt: skip
    c.set_xlim(0, mids.max() * 1.08)
    c.axhline(0, color="#cccccc", lw=0.6)
    c.set_xlabel("Counting-noise share of abundance variance (quintile median)")
    c.set_ylabel("Median abundance-where-present R$^2$")
    c.legend(frameon=False, fontsize=6)
    # d: within-cruise depth correlation vs gene length, by domain
    for dm, dc, dl in (("euk", EUK, "Eukaryotic KOs"), ("prok", PROK, "Prokaryotic KOs")):
        g = r[dom == dm]
        d.scatter(g["length_mean"], g["r_depth_within"], s=2, lw=0, color=dc, alpha=0.3, rasterized=True)
        lb = pd.qcut(g["length_mean"], 10)
        d.plot(g.groupby(lb, observed=True)["length_mean"].median(), g.groupby(lb, observed=True)["r_depth_within"].median(),
               color=dc, lw=1.4, marker="o", ms=3, mec="white", mew=0.4, label=f"{dl} (decile median)")  # fmt: skip
    d.set_xscale("log")
    d.set_xticks([200, 500, 1000, 2000, 5000])
    d.set_xticklabels(["200", "500", "1,000", "2,000", "5,000"])
    d.minorticks_off()
    d.axhline(0, color="#cccccc", lw=0.6)
    d.set_xlabel("Mean reference gene length (nt)")
    d.set_ylabel("Presence vs sequencing depth, within cruise (r)")
    d.legend(frameon=False, fontsize=6, loc="upper left")
    for axx, letter in zip((a, b, c, d), "abcd"):
        axx.text(-0.22 if axx is a else -0.14, 1.02, letter, transform=axx.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for suf in ("png", "pdf"):
        fig.savefig(FIG / f"hurdle_mechanism.{suf}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    if "--figure-only" in sys.argv:
        r = pd.read_csv(OUT / "mechanism.csv", index_col=0)
    else:
        r = per_ko()
        summarise(r)
    figure(r)
