#!/usr/bin/env python3
"""R69: patchy-KO presence with one base model (Granite PCs + depth + platform), +/- in-situ temperature, against the
eukaryote : prokaryote ratio (C) and temperature as references. Replaces the satellite / technical / weak classes of
R65 / R66.

Same tiles, presence and cruise-held-out scheme as scripts/hurdle_kos.py (1,111 tiles, 13 cruises, 5,808 KOs; rarefied
presence at 1e6 KO reads). Models per KO, each cruise held out, L2 logistic with the penalty (and, for Granite, the
number of PCs: 4 / 8 / 16 / 32) chosen by leave-cruise-out inside the training cruises:

  base_temp   Granite PCs + technical (platform + depth spline) + in-situ temperature spline        [new]
  c_tech      Granite-predicted C + technical                                                      [new]
  temp_tech   in-situ temperature spline + technical                                               [new]
  (R65, reused from hurdle_all.csv: tech = technical only; granite = Granite PCs; base = auc_*_both = Granite +
   technical. Identical folds, tiles and presence.)

Predicted C is nested, so no information from the held-out cruise enters: for held-out cruise c, a training tile in
cruise t gets C predicted by a ridge on 32 Granite PCs fitted on the cruises other than c and t; the tiles of c get C
predicted from all cruises but c. C is the technically adjusted target of the main model (R68 definition), ridge
penalty by leave-cruise-out (benchmark_eval.ridge_loco_alpha), fitted on all 1,215 benchmark tiles.

Parts:
  cpred    data/results/hurdle/c_nested.npz (per outer fold: predicted C for its training and test tiles)
  fit      --shard i n; appends to base_shard<i>.csv.partial (fixed columns, resumable)
  summary  base_all.csv (joined with hurdle_all.csv and ko_domain_reads.csv), base_summary.txt
  figure   figures/output/hurdle_base.{png,pdf}; depth_detection.csv; figure numbers appended to base_summary.txt
Usage:
  OMP_NUM_THREADS=2 ~/miniforge3/bin/python3 scripts/hurdle_base.py --parts cpred
  bash cluster/run_hurdle_base_local.sh            (shards in parallel, resumable)
  OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/hurdle_base.py --parts summary
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import SplineTransformer, StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import hurdle_kos as h  # noqa: E402
from benchmark_eval import ridge_loco_alpha  # noqa: E402

OUT, BX = h.OUT, h.BX
MODELS = ("base_temp", "c_tech", "temp_tech")
COLS = ["ko", "tiles", "prevalence"] + [f"{s}_{m}" for m in MODELS for s in ("auc", "devexp", "auc_within")] + ["pcs_median_base_temp"]


def data():
    tiles = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    tp = pd.read_csv(OUT / "presence_tile.csv.gz", index_col=0)
    depth = pd.read_csv(OUT / "tile_depth.csv", index_col=0)["log10_total_genes"]
    idx = tp.index
    t = tiles.loc[idx]
    cr, plat = t["cruise"].to_numpy(), (t["platform_group"] == "NovaSeq 6000 group").to_numpy()
    emb = pd.DataFrame(np.load(BX / "granite" / "tile_embeddings.npy"), index=tiles.index)
    return tiles, tp, idx, cr, plat, emb, depth.loc[idx].to_numpy(), t[["temperature_c"]].to_numpy()


def spline_fold(x, tr, te):
    s = StandardScaler().fit(x[tr])
    sp = SplineTransformer(n_knots=5, degree=3, extrapolation="linear").fit(s.transform(x[tr]))
    return sp.transform(s.transform(x[tr])), sp.transform(s.transform(x[te]))


def cpred():
    """Nested Granite predictions of C for every outer fold (hurdle tiles), fitted on all 1,215 benchmark tiles."""
    tiles, tp, idx, cr, *_ = data()
    emb = np.load(BX / "granite" / "tile_embeddings.npy")
    z = np.load(BX / "granite" / "predictions_C.npz", allow_pickle=True)
    yC = pd.Series(z["observed"][:, 0].astype(float), index=z["tile"]).reindex(tiles.index)  # technically adjusted C
    allcr = tiles["cruise"].to_numpy()
    rowpos = pd.Series(np.arange(len(tiles)), index=tiles.index)
    hpos = rowpos.loc[idx].to_numpy()

    def fit_predict(train_mask, pred_rows):
        trn = np.where(train_mask)[0]
        a, b = h.pcs(emb, trn, pred_rows)
        return ridge_loco_alpha(a, yC.to_numpy()[trn][:, None], b, allcr[trn])[:, 0]

    out = {}
    for c in np.unique(cr):
        t0 = time.time()
        test_rows = hpos[cr == c]
        pc_test = fit_predict(allcr != c, test_rows)
        pc_train = np.empty((cr != c).sum())
        tr_idx = np.where(cr != c)[0]
        for t in np.unique(cr[tr_idx]):
            sel = cr[tr_idx] == t
            pc_train[sel] = fit_predict((allcr != c) & (allcr != t), hpos[tr_idx[sel]])
        out[c] = (pc_train, pc_test)
        r_te = np.corrcoef(pc_test, yC.loc[idx[cr == c]])[0, 1]
        print(f"  {c:14s} test r(pred, obs C) {r_te:+.3f}  ({time.time() - t0:.0f} s)", flush=True)
    full = np.empty(len(idx))
    for c, (_, te) in out.items():
        full[cr == c] = te
    obs = yC.loc[idx].to_numpy()
    print(f"held-out C on hurdle tiles: r {np.corrcoef(full, obs)[0, 1]:.3f}, R2 {1 - ((obs - full) ** 2).sum() / ((obs - obs.mean()) ** 2).sum():.3f}")
    np.savez_compressed(OUT / "c_nested.npz", cruises=np.array(list(out)), observed=obs, heldout=full,
                        **{f"train_{c}": v[0] for c, v in out.items()}, **{f"test_{c}": v[1] for c, v in out.items()})  # fmt: skip


def fit(shard):
    tiles, tp, idx, cr, plat, emb, logd, temp = data()
    X = emb.loc[idx].to_numpy()
    tech = h.tech_matrix(plat, logd)
    zc = np.load(OUT / "c_nested.npz", allow_pickle=True)
    folds = []
    for c in np.unique(cr):
        te, tr = np.where(cr == c)[0], np.where(cr != c)[0]
        a, b = h.pcs(X, tr, te)
        st = StandardScaler().fit(tech[tr])
        sa, sb = st.transform(tech[tr]), st.transform(tech[te])
        ta, tb = spline_fold(temp, tr, te)
        ca, cb = zc[f"train_{c}"][:, None], zc[f"test_{c}"][:, None]
        cs = StandardScaler().fit(ca)
        folds.append(dict(te=te, tr=tr, a=a, b=b, ext=(np.hstack([sa, ta]), np.hstack([sb, tb])),
                          c=(np.hstack([cs.transform(ca), sa]), np.hstack([cs.transform(cb), sb])),
                          t=(np.hstack([ta, sa]), np.hstack([tb, sb]))))  # fmt: skip
    kos = list(tp.columns)[shard[0]::shard[1]] if shard else list(tp.columns)
    part = OUT / (f"base_shard{shard[0]}.csv.partial" if shard else "base.csv.partial")
    done = set(pd.read_csv(part)["ko"]) if part.exists() and part.stat().st_size else set()
    if done:
        print(f"resuming: {len(done)} KOs already in {part.name}", flush=True)
    t0 = time.time()
    for n, k in enumerate(kos):
        if k in done:
            continue
        y = tp[k].to_numpy()
        P = {m: np.empty(len(y)) for m in MODELS}
        kch = []
        for f in folds:
            g, ytr = cr[f["tr"]], y[f["tr"]]
            pb = h.logit_loco_k(f["a"], ytr, f["b"], g, extra=f["ext"])
            P["base_temp"][f["te"]] = pb[0] if isinstance(pb, tuple) else pb
            if isinstance(pb, tuple):
                kch.append(pb[1])
            P["c_tech"][f["te"]] = h.logit_loco(f["c"][0], ytr, f["c"][1], g)
            P["temp_tech"][f["te"]] = h.logit_loco(f["t"][0], ytr, f["t"][1], g)
        r = {"ko": k, "tiles": len(y), "prevalence": y.mean(), "pcs_median_base_temp": float(np.median(kch)) if kch else np.nan}
        for m, p in P.items():
            r[f"auc_{m}"], r[f"devexp_{m}"], r[f"auc_within_{m}"] = h.scores(y, p, cr)
        pd.DataFrame([r], columns=COLS).round(5).to_csv(part, mode="a", header=not part.exists() or part.stat().st_size == 0, index=False)
        if n % 25 == 0 or k in h.FOCUS:
            print(f"  {n:5d}/{len(kos)} {k} {h.FOCUS.get(k, ''):8s} within-cruise AUC base+temp {r['auc_within_base_temp']:.3f} "
                  f"C+tech {r['auc_within_c_tech']:.3f} temp+tech {r['auc_within_temp_tech']:.3f}  "
                  f"({(time.time() - t0) / (n + 1 - len(done & set(kos[:n + 1]))):.1f} s/KO)", flush=True)  # fmt: skip
    pd.read_csv(part).round(5).to_csv(OUT / part.name.replace(".partial", ""), index=False)


def summary():
    parts = {p.name.split(".")[0]: p for p in OUT.glob("base_shard*.csv.partial")}
    parts.update({p.stem: p for p in OUT.glob("base_shard*.csv")})
    f = pd.concat([pd.read_csv(p) for _, p in sorted(parts.items())]).drop_duplicates("ko").set_index("ko")
    r65 = pd.read_csv(OUT / "hurdle_all.csv", index_col="ko")
    dom = pd.read_csv(REPO / "data/results/ko_domain/ko_domain_reads.csv", index_col="ko")
    f = f.join(r65[["name", "in_core"] + [f"{s}_{m}" for m in ("tech", "granite", "both") for s in ("auc", "devexp", "auc_within")]
                   + ["abund_tiles", "abund_r2", "abund_technical_share"]]).join(dom[["euk_read_share", "reads"]])  # fmt: skip
    f = f.rename(columns=lambda c: c.replace("_both", "_base"))
    f["domain"] = np.select([f["euk_read_share"] > 0.9, f["euk_read_share"] <= 0.1], ["eukaryotic", "prokaryotic"], "mixed")
    f.round(5).to_csv(OUT / "base_all.csv")
    log = [f"KOs {len(f):,} of {len(r65):,}"]
    show = ["tech", "temp_tech", "c_tech", "granite", "base", "base_temp"]
    for sc in ("auc_within", "auc"):
        log.append(f"\nmedian {sc} (10-90%) by read domain:")
        rows = []
        for d, g in [("all", f)] + list(f.groupby("domain")):
            rows.append({"domain": d, "n": len(g), **{m: f"{g[f'{sc}_{m}'].median():.3f} ({g[f'{sc}_{m}'].quantile(.1):.2f}-"
                                                         f"{g[f'{sc}_{m}'].quantile(.9):.2f})" for m in show}})  # fmt: skip
        log.append(pd.DataFrame(rows).to_string(index=False))
    for a, b in (("base", "tech"), ("base_temp", "base"), ("base", "c_tech"), ("base", "temp_tech"), ("c_tech", "tech")):
        log.append(f"\nwithin-cruise AUC {a} - {b}:")
        for d, g in [("all", f)] + list(f.groupby("domain")):
            x = g[f"auc_within_{a}"] - g[f"auc_within_{b}"]
            log.append(f"  {d:12s} n {len(g):5,}  median {x.median():+.3f} (10-90% {x.quantile(.1):+.3f} to {x.quantile(.9):+.3f}); "
                       f">= +0.02 {(x >= .02).mean():.0%}, >= +0.05 {(x >= .05).mean():.0%}, <= -0.02 {(x <= -.02).mean():.0%}")  # fmt: skip
    (OUT / "base_summary.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))


def depth_detection():
    """Relative detection by within-cruise depth quintile, per read domain and gene-length tercile (within domain)."""
    tp = pd.read_csv(OUT / "presence_tile.csv.gz", index_col=0)
    cr = pd.read_csv(BX / "tiles.csv", index_col="tile_key").loc[tp.index, "cruise"]
    dep = pd.read_csv(OUT / "tile_depth.csv", index_col=0)["log10_total_genes"].loc[tp.index]
    q = dep.groupby(cr).transform(lambda s: pd.qcut(s.rank(method="first"), 5, labels=False)).to_numpy()
    rate = tp.groupby(q).mean()
    rel = rate / rate.mean()  # 5 quintiles x KOs, relative to the KO's mean detection
    f = pd.read_csv(OUT / "base_all.csv", index_col="ko")
    L = pd.read_csv(REPO / "transfer/results/ko_gc/ko_gc.csv", index_col="ko")["length_mean"].reindex(f.index)
    rows, cuts = [], {}
    for dm in ("prokaryotic", "mixed", "eukaryotic"):
        m = f["domain"] == dm
        cuts[dm] = [round(x) for x in L[m].quantile([1 / 3, 2 / 3])]
        terc = pd.qcut(L[m], 3, labels=["short", "medium", "long"])
        for tt in ("short", "medium", "long"):
            kk = terc.index[terc == tt]
            for i in range(5):
                rows.append({"domain": dm, "length": tt, "n": len(kk), "quintile": i + 1, "rel_median": rel.loc[i, kk].median(),
                             "rel_q25": rel.loc[i, kk].quantile(.25), "rel_q75": rel.loc[i, kk].quantile(.75),
                             "rate_median": rate.loc[i, kk].median()})  # fmt: skip
    dd = dep - dep.groupby(cr).transform("mean")
    C = pd.read_csv(BX / "targets" / "C.csv.gz", index_col="tile_key")["euk_prok_log_ratio"].reindex(tp.index)
    cc = C - C.groupby(cr).transform("mean")
    ok = cc.notna()
    return pd.DataFrame(rows), cuts, float(np.corrcoef(dd[ok], cc[ok])[0, 1])


def figure():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})  # fmt: skip
    FIG = REPO / "figures" / "output"
    f = pd.read_csv(OUT / "base_all.csv", index_col="ko")
    PROK, EUK, MIX = "#33a02c", "#6a3d9a", "#8c8c8c"
    dcol = {"prokaryotic": PROK, "mixed": MIX, "eukaryotic": EUK}
    cmap = LinearSegmentedColormap.from_list("dom", [PROK, "#c8c8c8", EUK])
    fig = plt.figure(figsize=(7.2, 5.8))
    a = fig.add_axes([0.09, 0.57, 0.38, 0.38])
    b = fig.add_axes([0.58, 0.57, 0.33, 0.38])
    cb = fig.add_axes([0.925, 0.62, 0.012, 0.28])
    c = fig.add_axes([0.09, 0.08, 0.38, 0.36])
    d = fig.add_axes([0.58, 0.08, 0.38, 0.36])
    log = []

    # a: within-cruise presence AUC by model and read domain (median, IQR)
    models = (("tech", "Depth + platform", "#E69F00", "s"), ("temp_tech", "Temperature + depth + platform", "#c0457e", "D"),
              ("c_tech", "Predicted eukaryote : prokaryote + depth + platform", "#56B4E9", "^"),
              ("base", "Granite + depth + platform", "#0072B2", "o"))  # fmt: skip
    doms = ("prokaryotic", "mixed", "eukaryotic")
    for j, (m, lb, cc, mk) in enumerate(models):
        for i, dm in enumerate(doms):
            v = f.loc[f["domain"] == dm, f"auc_within_{m}"]
            x = i + (j - 1.5) * 0.17
            a.plot([x, x], [v.quantile(.25), v.quantile(.75)], color=cc, lw=1.2, solid_capstyle="butt")
            a.scatter(x, v.median(), s=18, marker=mk, color=cc, lw=0, zorder=3, label=lb if i == 0 else None)
    a.axhline(0.5, color="#cccccc", lw=0.6, zorder=0)
    a.set_xticks(range(3))
    a.set_xticklabels([f"{dm.capitalize()}\n({(f['domain'] == dm).sum():,} KOs)" for dm in doms])
    a.set_ylabel("Presence AUC within cruise (median, IQR)")
    a.set_xlim(-0.5, 2.5)
    a.legend(frameon=False, fontsize=5.8, loc="upper left", handletextpad=0.2, bbox_to_anchor=(0, 1.04))
    a.set_ylim(0.44, 0.86)

    # b: Granite + technical vs predicted C + technical, coloured by eukaryotic read share
    o = f.sample(frac=1, random_state=0)  # random draw order so neither domain hides the other
    sc = b.scatter(o["auc_within_c_tech"], o["auc_within_base"], c=o["euk_read_share"], cmap=cmap, vmin=0, vmax=1, s=2.2,
                   lw=0, alpha=0.75, rasterized=True)  # fmt: skip
    b.plot([0.25, 1], [0.25, 1], color="#555555", lw=0.6, ls=(0, (3, 2)))
    LABEL_OFF = {"nifH": (-22, -5), "nasA": (5, -4), "nifD": (-23, 2), "amoB": (5, 2), "nifK": (5, -5), "phnJ": (4, -7)}
    for k, gname in h.FOCUS.items():
        if k in f.index and gname in ("nifH", "nifD", "nifK", "amoB / pmoB", "nasA", "phnJ"):
            r = f.loc[k]
            gn = gname.split(" /")[0]
            b.scatter(r["auc_within_c_tech"], r["auc_within_base"], s=16, facecolor="none", edgecolor="black", lw=0.8, zorder=4)
            b.annotate(gn, (r["auc_within_c_tech"], r["auc_within_base"]), xytext=LABEL_OFF.get(gn, (4, 3)),
                       textcoords="offset points", fontsize=6.5, style="italic")  # fmt: skip
            log.append(f"  focus {gn:5s} {k}: domain {r['domain']}, euk read share {r['euk_read_share']:.3f}, AUC within: tech "
                       f"{r['auc_within_tech']:.3f}, temp+tech {r['auc_within_temp_tech']:.3f}, C+tech {r['auc_within_c_tech']:.3f}, "
                       f"Granite+tech {r['auc_within_base']:.3f}, +temp {r['auc_within_base_temp']:.3f}; abundance R2 {r['abund_r2']:.3f}")  # fmt: skip
    b.set_xlim(0.25, 1)
    b.set_ylim(0.25, 1)
    b.set_xlabel("Presence AUC within cruise,\npredicted eukaryote : prokaryote + depth + platform")
    b.set_ylabel("Presence AUC within cruise, Granite + depth + platform")
    cbar = fig.colorbar(sc, cax=cb)
    cbar.set_ticks([0, 0.5, 1])
    cbar.outline.set_linewidth(0.5)
    cbar.set_label("Eukaryotic share of the KO's reads", fontsize=6.5)
    cbar.solids.set_rasterized(True)

    # c: abundance-where-present R2 by eukaryotic read share (deciles of share)
    edges = np.linspace(0, 1, 11)
    bins = pd.cut(f["euk_read_share"], edges, include_lowest=True)
    g = f.groupby(bins, observed=True)["abund_r2"]
    mid = (edges[:-1] + edges[1:]) / 2
    med, q1, q3, n = g.median().to_numpy(), g.quantile(.25).to_numpy(), g.quantile(.75).to_numpy(), g.size().to_numpy()
    c.fill_between(mid, q1, q3, color="#0072B2", alpha=0.18, lw=0)
    c.plot(mid, med, color="#0072B2", lw=1.2, marker="o", ms=3)
    for x, y, k in zip(mid, q3, n):
        c.text(x, y + 0.012, f"{k:,}", ha="center", va="bottom", fontsize=5.3, color="#555555")
    c.axhline(0, color="#cccccc", lw=0.6, zorder=0)
    c.set_xlim(0, 1)
    c.set_xlabel("Eukaryotic share of the KO's reads")
    c.set_ylabel("Abundance where present, R² (Granite, cruise held out)")
    log.append("\nabundance-where-present R2 by eukaryotic read share bin (n, median, IQR):")
    for e0, e1, k, m_, a_, b_ in zip(edges[:-1], edges[1:], n, med, q1, q3):
        log.append(f"  {e0:.1f}-{e1:.1f}  n {k:5,}  median {m_:+.3f}  IQR {a_:+.3f} to {b_:+.3f}")
    for dm in doms:
        v = f.loc[f["domain"] == dm, "abund_r2"]
        log.append(f"  {dm:12s} n {len(v):5,}  median {v.median():+.3f}  IQR {v.quantile(.25):+.3f} to {v.quantile(.75):+.3f}  "
                   f"R2 > 0.1: {(v > 0.1).mean():.0%}")  # fmt: skip

    # d: relative detection vs within-cruise depth quintile, by gene length, prokaryotic vs eukaryotic
    dd, cuts, r_dc = depth_detection()
    dd.round(4).to_csv(OUT / "depth_detection.csv", index=False)
    shade = {"short": 0.45, "medium": 0.72, "long": 1.0}
    for dm, cc in (("prokaryotic", PROK), ("eukaryotic", EUK)):
        for tt in ("short", "medium", "long"):
            s_ = dd[(dd["domain"] == dm) & (dd["length"] == tt)]
            rgb = np.array(matplotlib.colors.to_rgb(cc))
            col = 1 - shade[tt] * (1 - rgb)
            d.plot(s_["quintile"], s_["rel_median"], color=col, lw=1.3 if tt == "long" else 1.0, marker="o", ms=2.6,
                   label=f"{dm.capitalize()}, {tt}")  # fmt: skip
    d.axhline(1, color="#cccccc", lw=0.6, zorder=0)
    d.set_xticks(range(1, 6))
    d.set_xticklabels(["1\nshallowest", "2", "3", "4", "5\ndeepest"])
    d.set_xlabel("Sequencing depth, fifth of the tiles within each cruise")
    d.set_ylabel("Detection rate relative to the KO's average")
    d.legend(frameon=False, fontsize=5.8, ncol=2, loc="upper left", handletextpad=0.3, columnspacing=0.8)
    d.set_ylim(0.75, 1.32)
    log.append(f"\ndepth detection: length tercile cuts (nt) {cuts}; within-cruise r(log10 depth, observed C) {r_dc:+.3f}")
    for (dm, tt), s_ in dd.groupby(["domain", "length"], sort=False):
        log.append(f"  {dm:12s} {tt:6s} n {s_['n'].iloc[0]:4d}  relative detection Q1..Q5 "
                   + " ".join(f"{v:.2f}" for v in s_["rel_median"]) + "  | absolute rate " + " ".join(f"{v:.2f}" for v in s_["rate_median"]))  # fmt: skip

    for ax, letter, dx in ((a, "a", -0.17), (b, "b", -0.22), (c, "c", -0.17), (d, "d", -0.17)):
        ax.text(dx, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for suf in ("png", "pdf"):
        fig.savefig(FIG / f"hurdle_base.{suf}", dpi=300)
    plt.close(fig)
    with open(OUT / "base_summary.txt", "a") as fh:
        fh.write("\nfigure numbers:\n" + "\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", default=["fit"], choices=["cpred", "fit", "summary", "figure"])
    ap.add_argument("--shard", nargs=2, type=int)
    a = ap.parse_args()
    if "cpred" in a.parts:
        cpred()
    if "fit" in a.parts:
        fit(a.shard)
    if "summary" in a.parts:
        summary()
    if "figure" in a.parts:
        figure()
