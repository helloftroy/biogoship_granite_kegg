#!/usr/bin/env python3
"""R69: patchy-KO presence with one base model (Granite PCs + depth + platform), +/- in-situ temperature, against the
eukaryote : prokaryote ratio (C) and temperature as references. Replaces the satellite / technical / weak classes of
R65 / R66.

Same tiles, presence and cruise-held-out scheme as scripts/hurdle_kos.py (1,111 tiles, 13 cruises, 5,808 KOs; rarefied
presence at 1e6 KO reads). Models per KO, each cruise held out, L2 logistic with the penalty (and, for Granite, the
number of PCs: 4 / 8 / 16 / 32) chosen by leave-cruise-out inside the training cruises:

  c_tech      Granite-predicted C + technical
  (from hurdle_kos.py, hurdle_all.csv: tech = technical only; granite = Granite PCs; base = auc_*_both = Granite +
   technical. Identical folds, tiles and presence.)
Technical = platform group + spline of log10 read pairs sequenced (R70; before R70: log10 KO-assigned reads).
R69 also fitted base_temp (base + temperature spline) and temp_tech; temperature added nothing (+0.001) and was dropped.

Predicted C is nested, so no information from the held-out cruise enters: for held-out cruise c, a training tile in
cruise t gets C predicted by a ridge on 32 Granite PCs fitted on the cruises other than c and t; the tiles of c get C
predicted from all cruises but c. C is the technically adjusted target of the main model (R68 definition), ridge
penalty by leave-cruise-out (benchmark_eval.ridge_loco_alpha), fitted on all 1,215 benchmark tiles.

Parts:
  cpred    data/results/hurdle/c_nested.npz (per outer fold: predicted C for its training and test tiles)
  fit      --shard i n; appends to base_shard<i>.csv.partial (fixed columns, resumable)
  summary  base_all.csv (joined with hurdle_all.csv and ko_domain_reads.csv), base_summary.txt
  figure   figures/output/hurdle_base.{png,pdf}; depth_detection.csv, core_domain_depth.csv; figure numbers appended
           to base_summary.txt. Core KOs = the 3,933 KOs of the main (full-KO) model, R2 from main_model/ko_ceiling.csv.
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
MODELS = ("c_tech",)  # R70: temperature dropped (base_temp - base +0.001 in R69); base / tech / granite come from hurdle_kos
COLS = ["ko", "tiles", "prevalence"] + [f"{s}_{m}" for m in MODELS for s in ("auc", "devexp", "auc_within")]


def data():
    tiles = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    tp = pd.read_csv(OUT / "presence_tile.csv.gz", index_col=0)
    depth = pd.read_csv(OUT / "tile_depth.csv", index_col=0)[h.DEPTH_COL]  # read pairs since R70
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
        st = StandardScaler().fit(tech[tr])
        sa, sb = st.transform(tech[tr]), st.transform(tech[te])
        ca, cb = zc[f"train_{c}"][:, None], zc[f"test_{c}"][:, None]
        cs = StandardScaler().fit(ca)
        folds.append(dict(te=te, tr=tr, c=(np.hstack([cs.transform(ca), sa]), np.hstack([cs.transform(cb), sb]))))
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
        for f in folds:
            g, ytr = cr[f["tr"]], y[f["tr"]]
            P["c_tech"][f["te"]] = h.logit_loco(f["c"][0], ytr, f["c"][1], g)
        r = {"ko": k, "tiles": len(y), "prevalence": y.mean()}
        for m, p in P.items():
            r[f"auc_{m}"], r[f"devexp_{m}"], r[f"auc_within_{m}"] = h.scores(y, p, cr)
        pd.DataFrame([r], columns=COLS).round(5).to_csv(part, mode="a", header=not part.exists() or part.stat().st_size == 0, index=False)
        if n % 25 == 0 or k in h.FOCUS:
            print(f"  {n:5d}/{len(kos)} {k} {h.FOCUS.get(k, ''):8s} within-cruise AUC C+tech {r['auc_within_c_tech']:.3f}  "
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
    show = ["tech", "c_tech", "granite", "base"]
    for sc in ("auc_within", "auc"):
        log.append(f"\nmedian {sc} (10-90%) by read domain:")
        rows = []
        for d, g in [("all", f)] + list(f.groupby("domain")):
            rows.append({"domain": d, "n": len(g), **{m: f"{g[f'{sc}_{m}'].median():.3f} ({g[f'{sc}_{m}'].quantile(.1):.2f}-"
                                                         f"{g[f'{sc}_{m}'].quantile(.9):.2f})" for m in show}})  # fmt: skip
        log.append(pd.DataFrame(rows).to_string(index=False))
    for a, b in (("base", "tech"), ("base", "c_tech"), ("c_tech", "tech"), ("base", "granite")):
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
    eff = pd.read_csv(REPO / "data/results/sequencing_effort/tile_effort.csv", index_col="tile_key").loc[tp.index]
    dep = eff["log10_read_pairs"]  # raw sequencing effort (R70), not KO-assigned reads
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
    cc = eff["C"] - eff["C"].groupby(cr).transform("mean")
    ok = cc.notna()
    return pd.DataFrame(rows), cuts, float(np.corrcoef(dd[ok], cc[ok])[0, 1])


def core_table():
    """Core (full-KO model) KOs: held-out Granite R2 (main model, ko_ceiling.csv), read domain, gene length and the
    within-cruise correlation of abundance (CLR, before technical adjustment) with log10 read pairs sequenced (raw
    effort, R70), over all 1,215 tiles."""
    tiles = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    X = pd.read_csv(BX / "targets" / "clr_all.csv.gz", index_col="tile_key")
    eff = pd.read_csv(REPO / "data/results/sequencing_effort/tile_effort.csv", index_col="tile_key").loc[X.index]
    cr, dp = tiles.loc[X.index, "cruise"], eff["log10_read_pairs"]  # raw sequencing effort (R70)
    Xd, dd = X - X.groupby(cr).transform("mean"), dp - dp.groupby(cr).transform("mean")
    k = pd.read_csv(REPO / "data/results/main_model/ko_ceiling.csv", index_col="ko")[["r2"]]
    k["r_depth_within"] = Xd.mul(dd, axis=0).sum() / np.sqrt((Xd**2).sum() * (dd**2).sum())
    k = k.join(pd.read_csv(REPO / "data/results/ko_domain/ko_domain_reads.csv", index_col="ko")[["euk_read_share"]])
    k["domain"] = np.select([k["euk_read_share"] > 0.9, k["euk_read_share"] <= 0.1], ["eukaryotic", "prokaryotic"], "mixed")
    k["length_mean"] = pd.read_csv(REPO / "transfer/results/ko_gc/ko_gc.csv", index_col="ko")["length_mean"].reindex(k.index)
    k["length_tercile"] = ""
    for dm in ("prokaryotic", "mixed", "eukaryotic"):
        m = k["domain"] == dm
        k.loc[m, "length_tercile"] = pd.qcut(k.loc[m, "length_mean"], 3, labels=["short", "medium", "long"]).astype(str)
    k.round(5).to_csv(OUT / "core_domain_depth.csv")
    return k


def figure():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    from scipy.stats import spearmanr

    plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})  # fmt: skip
    FIG = REPO / "figures" / "output"
    f = pd.read_csv(OUT / "base_all.csv", index_col="ko")
    core = core_table()
    PROK, EUK, MIX = "#33a02c", "#6a3d9a", "#7f7f7f"
    dcol = {"prokaryotic": PROK, "mixed": MIX, "eukaryotic": EUK}
    doms = ("prokaryotic", "mixed", "eukaryotic")
    CORE_C, PATCHY_C = "#D55E00", "#0072B2"
    cmap = LinearSegmentedColormap.from_list("dom", [PROK, "#c8c8c8", EUK])

    def tint(c, s):
        rgb = np.array(to_rgb(c))
        return 1 - s * (1 - rgb)

    SHADE = {"short": 0.42, "medium": 0.7, "long": 1.0}
    fig = plt.figure(figsize=(7.2, 8.6))
    a = fig.add_axes([0.09, 0.715, 0.36, 0.255])
    b = fig.add_axes([0.60, 0.715, 0.37, 0.255])
    c = fig.add_axes([0.09, 0.375, 0.32, 0.255])
    cb = fig.add_axes([0.425, 0.40, 0.010, 0.20])
    d = fig.add_axes([0.60, 0.375, 0.37, 0.255])
    e = [fig.add_axes([0.09 + i * 0.125, 0.06, 0.105, 0.225]) for i in range(3)]
    g6 = fig.add_axes([0.60, 0.06, 0.37, 0.225])
    log = []

    def dotrange(ax, x, v, cc, mk, lb):
        ax.plot([x, x], [v.quantile(.25), v.quantile(.75)], color=cc, lw=1.2, solid_capstyle="butt")
        ax.scatter(x, v.median(), s=18, marker=mk, color=cc, lw=0, zorder=3, label=lb)

    # a: presence AUC (patchy KOs)
    models = (("tech", "Depth + platform", "#E69F00", "s"),
              ("c_tech", "Predicted eukaryote : prokaryote\n+ depth + platform", "#56B4E9", "^"),
              ("base", "Granite + depth + platform", PATCHY_C, "o"))  # fmt: skip
    for j, (m, lb, cc, mk) in enumerate(models):
        for i, dm in enumerate(doms):
            dotrange(a, i + (j - 1) * 0.22, f.loc[f["domain"] == dm, f"auc_within_{m}"], cc, mk, lb if i == 0 else None)
    a.axhline(0.5, color="#cccccc", lw=0.6, zorder=0)
    a.set_xticks(range(3))
    a.set_xticklabels([f"{dm.capitalize()}\n({(f['domain'] == dm).sum():,})" for dm in doms])
    a.set_xlim(-0.5, 2.5)
    a.set_ylim(0.44, 0.84)
    a.set_ylabel("Presence AUC within cruise\n(patchy KOs; median, IQR)")
    a.legend(frameon=False, fontsize=5.8, loc="upper left", handletextpad=0.2, bbox_to_anchor=(0, 1.06))

    # b: abundance R2, core (main model) vs patchy (where present), side by side, scored separately
    for j, (src, lb, cc, mk) in enumerate(((core, "Core KOs (main model)", CORE_C, "D"),
                                            (f.rename(columns={"abund_r2": "r2"}), "Patchy KOs, where present", PATCHY_C, "o"))):  # fmt: skip
        for i, dm in enumerate(doms):
            v = src.loc[src["domain"] == dm, "r2"]
            dotrange(b, i + (j - 0.5) * 0.28, v, cc, mk, lb if i == 0 else None)
            b.text(i + (j - 0.5) * 0.28, v.quantile(.75) + 0.012, f"{len(v):,}", ha="center", va="bottom", fontsize=5.3, color="#555555")
    b.axhline(0, color="#cccccc", lw=0.6, zorder=0)
    b.set_xticks(range(3))
    b.set_xticklabels([dm.capitalize() for dm in doms])
    b.set_xlim(-0.5, 2.5)
    b.set_ylabel("Abundance R², Granite, cruise held out\n(median, IQR)")
    b.legend(frameon=False, fontsize=5.8, loc="upper left", handletextpad=0.2, bbox_to_anchor=(0, 1.06))
    log.append("abundance R2 by read domain (median, IQR, share > 0.1):")
    for nm, src in (("core", core), ("patchy where present", f.rename(columns={"abund_r2": "r2"}))):
        for dm in doms:
            v = src.loc[src["domain"] == dm, "r2"]
            log.append(f"  {nm:21s} {dm:12s} n {len(v):5,}  median {v.median():+.3f}  IQR {v.quantile(.25):+.3f} to "
                       f"{v.quantile(.75):+.3f}  > 0.1 {(v > 0.1).mean():.0%}")  # fmt: skip
    log.append(f"  core KOs also in the patchy set: {len(core.index.intersection(f.index))}")

    # c: per-KO presence AUC, Granite + tech vs predicted C + tech
    o = f.sample(frac=1, random_state=0)
    sc = c.scatter(o["auc_within_c_tech"], o["auc_within_base"], c=o["euk_read_share"], cmap=cmap, vmin=0, vmax=1, s=2,
                   lw=0, alpha=0.75, rasterized=True)  # fmt: skip
    c.plot([0.25, 1], [0.25, 1], color="#555555", lw=0.6, ls=(0, (3, 2)))
    LABEL_AT = {"amoB": (0.50, 0.97), "nifD": (0.30, 0.91), "nasA": (0.28, 0.66), "nifH": (0.28, 0.60),
                "nifK": (0.78, 0.42), "phnJ": (0.88, 0.62)}  # fmt: skip
    for k, gname in h.FOCUS.items():
        gn = gname.split(" /")[0]
        if k in f.index and gn in LABEL_AT:
            r = f.loc[k]
            xy = (r["auc_within_c_tech"], r["auc_within_base"])
            c.scatter(*xy, s=16, facecolor="none", edgecolor="black", lw=0.8, zorder=4)
            c.annotate(gn, xy, xytext=LABEL_AT[gn], textcoords="data", fontsize=6.8, style="italic", fontweight="bold",
                       va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="#333333", shrinkA=1, shrinkB=3))  # fmt: skip
            log.append(f"  focus {gn:5s} {k}: euk read share {r['euk_read_share']:.3f}, AUC within: tech {r['auc_within_tech']:.3f}, "
                       f"C+tech {r['auc_within_c_tech']:.3f}, Granite+tech {r['auc_within_base']:.3f}")  # fmt: skip
    c.set_xlim(0.25, 1)
    c.set_ylim(0.25, 1)
    c.set_xlabel("Presence AUC, predicted eukaryote :\nprokaryote + depth + platform")
    c.set_ylabel("Presence AUC, Granite + depth + platform")
    cbar = fig.colorbar(sc, cax=cb)
    cbar.set_ticks([0, 0.5, 1])
    cbar.outline.set_linewidth(0.5)
    cbar.set_label("Eukaryotic share of reads", fontsize=6.3)
    cbar.solids.set_rasterized(True)

    # d: abundance R2 vs eukaryotic read share, core and patchy
    edges = np.linspace(0, 1, 11)
    mid = (edges[:-1] + edges[1:]) / 2
    log.append("\nabundance R2 by eukaryotic read share bin (n, median, IQR):")
    for nm, src, cc, dy in (("core", core, CORE_C, 0.012), ("patchy", f.rename(columns={"abund_r2": "r2"}), PATCHY_C, -0.012)):
        g = src.groupby(pd.cut(src["euk_read_share"], edges, include_lowest=True), observed=False)["r2"]
        med, q1, q3, n = g.median().to_numpy(), g.quantile(.25).to_numpy(), g.quantile(.75).to_numpy(), g.size().to_numpy()
        d.fill_between(mid, q1, q3, color=cc, alpha=0.15, lw=0)
        d.plot(mid, med, color=cc, lw=1.2, marker="D" if nm == "core" else "o", ms=2.8,
               label="Core KOs (main model)" if nm == "core" else "Patchy KOs, where present")  # fmt: skip
        for e0, e1, k_, m_, a_, b_ in zip(edges[:-1], edges[1:], n, med, q1, q3):
            log.append(f"  {nm:6s} {e0:.1f}-{e1:.1f}  n {k_:5,}  median {m_:+.3f}  IQR {a_:+.3f} to {b_:+.3f}")
        log.append(f"  {nm:6s} Spearman(R2, euk read share) {spearmanr(src['r2'], src['euk_read_share'], nan_policy='omit')[0]:.3f}")
    d.axhline(0, color="#cccccc", lw=0.6, zorder=0)
    d.set_xlim(0, 1)
    d.set_xlabel("Eukaryotic share of the KO's reads")
    d.set_ylabel("Abundance R², Granite, cruise held out")
    d.legend(frameon=False, fontsize=5.8, loc="upper left", handletextpad=0.3)

    # e: patchy KOs, detection by within-cruise depth fifth, per domain, by gene length
    dd, cuts, r_dc = depth_detection()
    dd.round(4).to_csv(OUT / "depth_detection.csv", index=False)
    for ax, dm in zip(e, doms):
        for tt in ("short", "medium", "long"):
            s_ = dd[(dd["domain"] == dm) & (dd["length"] == tt)]
            ax.plot(s_["quintile"], 100 * (s_["rel_median"] - 1), color=tint(dcol[dm], SHADE[tt]), lw=1.3 if tt == "long" else 1.0,
                    marker="o", ms=2.3, label=f"{tt.capitalize()} genes")  # fmt: skip
        ax.axhline(0, color="#cccccc", lw=0.6, zorder=0)
        ax.set_ylim(-22, 22)
        ax.set_xticks(range(1, 6))
        ax.text(0.5, 1.0, dm.capitalize(), transform=ax.transAxes, ha="center", va="bottom", fontsize=6.8, color=dcol[dm])
        if ax is not e[0]:
            ax.set_yticklabels([])
    from matplotlib.lines import Line2D

    LEN_HANDLES = [Line2D([], [], color=tint("#444444", SHADE[tt]), marker="o", ms=2.3, lw=1.3 if tt == "long" else 1.0,
                          label=f"{tt.capitalize()} genes") for tt in ("short", "medium", "long")]  # fmt: skip
    e[0].legend(handles=LEN_HANDLES, frameon=False, fontsize=5.4, loc="upper left", handlelength=1.4, handletextpad=0.3,
                borderaxespad=0.1)  # fmt: skip
    e[0].set_ylabel("Detection vs the KO's average (%)\n(patchy KOs)")
    e[1].set_xlabel("Read pairs sequenced, within cruise\n(fifths of the tiles; 1 = fewest, 5 = most)")

    # f: core KOs, within-cruise correlation of abundance with depth, by gene length
    for j, tt in enumerate(("short", "medium", "long")):
        for i, dm in enumerate(doms):
            v = core.loc[(core["domain"] == dm) & (core["length_tercile"] == tt), "r_depth_within"]
            dotrange(g6, i + (j - 1) * 0.22, v, tint(dcol[dm], SHADE[tt]), "o", None)
    g6.axhline(0, color="#cccccc", lw=0.6, zorder=0)
    g6.set_xticks(range(3))
    g6.set_xticklabels([f"{dm.capitalize()}\n({(core['domain'] == dm).sum():,})" for dm in doms])
    g6.set_xlim(-0.5, 2.5)
    g6.set_ylabel("Within-cruise r, abundance vs\nread pairs sequenced (core KOs)")
    g6.legend(handles=[Line2D([], [], color=tint("#444444", SHADE[tt]), marker="o", ms=3, lw=1.2, label=f"{tt.capitalize()} genes")
                       for tt in ("short", "medium", "long")], frameon=False, fontsize=5.8, loc="upper right", handlelength=1.2)
    log.append(f"\npatchy depth detection: length tercile cuts (nt) {cuts}; within-cruise r(log10 read pairs, observed C) {r_dc:+.3f}")
    for (dm, tt), s_ in dd.groupby(["domain", "length"], sort=False):
        log.append(f"  {dm:12s} {tt:6s} n {s_['n'].iloc[0]:4d}  % vs KO average Q1..Q5 "
                   + " ".join(f"{100 * (v - 1):+.0f}" for v in s_["rel_median"]))  # fmt: skip
    log.append("\ncore KOs: within-cruise r(CLR abundance, log10 read pairs), median (IQR) by domain x length tercile:")
    for dm in doms:
        g = core[core["domain"] == dm]
        cut = [round(x) for x in g["length_mean"].quantile([1 / 3, 2 / 3])]
        parts = [f"{tt} {g.loc[g['length_tercile'] == tt, 'r_depth_within'].median():+.3f} "
                 f"({g.loc[g['length_tercile'] == tt, 'r_depth_within'].quantile(.25):+.3f} to "
                 f"{g.loc[g['length_tercile'] == tt, 'r_depth_within'].quantile(.75):+.3f})" for tt in ("short", "medium", "long")]  # fmt: skip
        log.append(f"  {dm:12s} n {len(g):5,} cuts {cut}: " + "; ".join(parts)
                   + f"; Spearman(r, length) {spearmanr(g['r_depth_within'], g['length_mean'])[0]:.3f}")  # fmt: skip

    for ax, letter, dx in ((a, "a", -0.2), (b, "b", -0.2), (c, "c", -0.25), (d, "d", -0.2), (e[0], "e", -0.6), (g6, "f", -0.2)):
        ax.text(dx, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
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
