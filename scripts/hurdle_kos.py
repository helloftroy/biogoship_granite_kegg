#!/usr/bin/env python3
"""Hurdle models for patchy KOs (nifH, amoA, nasA, phnJ, ...) that the 90 % prevalence filter drops (RESULTS.md R65).

Unit and folds as the main model (benchmark tiles, data/results/benchmark_export/tiles.csv; one cruise held out).

build  Per sample: KO gene counts summed over libraries (data/results/hurdle/sample_counts.pkl, from
       ko_counts_cov50_all.tsv with the build_ko_targets.py library -> sample link and filters).
       Rarefied presence: a KO is present if it would be expected at least once in a library of G genes,
       i.e. count / total >= 1 / G (G = --depth, default 1,000,000; samples with fewer than G genes are dropped).
       This fixes the detection threshold in relative terms, so deeper libraries do not detect more by depth alone.
       Abundance where present: log(count / total) minus the mean log share of the 3,933 core KOs (the CLR reference
       of the main targets), so it is on the CLR scale.
       Tile level: presence = majority of the tile's kept samples present; abundance = mean over present samples.
       Hurdle set: KOs present in 20-90 % of tiles (--min-prev / --max-prev), plus the focus KOs whatever their
       prevalence.
fit    Presence (logistic, L2; penalty C chosen by leave-one-cruise-out inside the training cruises), three models:
         tech      platform group + cubic spline (5 knots) of log10 total genes per tile (technical only)
         granite   32 Granite PCs only
         both      Granite PCs + the technical columns
       Scores pooled over held-out cruises: AUC, deviance explained (1 - deviance / deviance of the training-prevalence
       null). Granite adds information about presence if both beats tech; presence is "technical" if tech alone is
       about as good as Granite (AUC_tech >= AUC_granite - 0.02).
       Abundance where present (tiles with the KO present): technical removal in-sample (platform + depth spline), ridge
       on 32 Granite PCs, penalty by leave-cruise-out (as the main model); pooled R2 of the non-technical variance.
       Also, for the focus KOs: tech AUC with unrarefied presence (any count > 0), to show what rarefaction removes.
       --shard i n splits the KOs.
summary / figure: per-KO table, focus-KO table, figures/output/hurdle_kos.png.

Outputs: data/results/hurdle/
Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/hurdle_kos.py --parts build
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/hurdle_kos.py --parts fit --shard 0 3   (x3)
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/hurdle_kos.py --parts summary figure
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import SplineTransformer, StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from benchmark_eval import ridge_loco_alpha  # noqa: E402

BX = REPO / "data" / "results" / "benchmark_export"
OUT = REPO / "data" / "results" / "hurdle"
FIG = REPO / "figures" / "output"
CS = np.logspace(-3, 1, 5)
KS = (4, 8, 16, 32)  # Granite PCs tried for presence; chosen with C inside the training cruises
FOCUS = {"K02588": "nifH", "K02586": "nifD", "K02591": "nifK", "K10944": "amoA / pmoA", "K10945": "amoB / pmoB",
         "K10946": "amoC / pmoC", "K00372": "nasA", "K06163": "phnJ", "K06164": "phnI", "K06165": "phnH",
         "K06166": "phnG", "K05780": "phnL", "K05781": "phnK", "K06162": "phnM"}  # fmt: skip
warnings.filterwarnings("ignore")

# Columns of a fit row. Before 2026-11-24 the .partial writer took the header from the first KO and appended each later
# row with only the blocks that ran for it (13 / 16 / 17 / 20 fields), so the file is not a rectangular CSV.
# read_partial() maps each row by its field count; new rows are always written with all 20 columns.
_BASE = ["ko", "tiles", "prevalence", "granite_pcs_median"] + [
    f"{s}_{m}" for m in ("tech", "granite", "both") for s in ("auc", "devexp", "auc_within")]
_UNR = ["prevalence_unrarefied", "auc_tech_unrarefied", "devexp_tech_unrarefied", "auc_within_tech_unrarefied"]
_AB = ["abund_tiles", "abund_r2", "abund_technical_share"]
FIT_COLS = _BASE + _UNR + _AB
_BY_N = {len(c): c for c in (_BASE, _BASE + _AB, _BASE + _UNR, FIT_COLS)}


def read_partial(path):
    """Read a fit .partial file whose rows may have 13 / 16 / 17 / 20 fields (see FIT_COLS)."""
    rows = []
    with open(path) as fh:
        next(fh, None)  # header (from whichever KO came first)
        for line in fh:
            v = line.rstrip("\n").split(",")
            if v == [""]:
                continue
            cols = _BY_N.get(len(v))
            if cols is None:
                raise ValueError(f"{path}: row with {len(v)} fields: {line[:80]}")
            rows.append(dict(zip(cols, v)))
    f = pd.DataFrame(rows, columns=FIT_COLS).replace("", np.nan)
    return f.astype({c: float for c in FIT_COLS if c != "ko"})


def build(G, lo, hi):
    S = pd.read_pickle(OUT / "sample_counts.pkl").astype(float)
    T = S.sum(1)
    core = pd.read_csv(BX / "targets" / "clr_all.csv.gz", nrows=0).columns.drop("tile_key")
    keep = T >= G
    P = S[keep].div(T[keep], axis=0)
    pres = P >= 1 / G
    raw = S[keep] > 0
    logcore = np.log(P[core].where(P[core] > 0)).mean(1)  # core KOs are > 0 almost everywhere
    tiles = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    members = tiles["sample_keys"].str.split(";").explode()
    members = members[members.isin(pres.index)]
    tp = pres.loc[members.values].groupby(members.index.values).mean() >= 0.5
    tr = raw.loc[members.values].groupby(members.index.values).mean() >= 0.5
    ab = np.log(P.loc[members.values].where(pres.loc[members.values])).sub(logcore.loc[members.values].values, axis=0)
    ta = ab.groupby(members.index.values).mean()
    tlog = np.log10(T.loc[members.values]).groupby(members.index.values).mean()
    prev = tp.mean()
    sel = prev[((prev >= lo) & (prev <= hi)) | prev.index.isin(list(FOCUS))].index
    sel = sel[prev[sel] > 0]
    tp, tr, ta = tp[sel], tr[sel], ta[sel]
    info = pd.DataFrame({"tiles_kept": len(tp), "samples_dropped": int((~keep).sum())}, index=[0])
    tp.astype(np.int8).to_csv(OUT / "presence_tile.csv.gz")
    tr.loc[:, tr.columns.isin(list(FOCUS))].astype(np.int8).to_csv(OUT / "presence_unrarefied_focus.csv.gz")
    ta.astype(np.float32).to_csv(OUT / "abundance_tile.csv.gz", float_format="%.5g")
    tlog.rename("log10_total_genes").to_csv(OUT / "tile_depth.csv")
    print(f"G {G:,}: samples kept {int(keep.sum()):,} of {len(T):,} (dropped {int((~keep).sum())}); benchmark tiles with a kept "
          f"sample {len(tp):,} of {len(tiles):,}; hurdle KOs {len(sel):,} (tile prevalence {lo:.0%}-{hi:.0%}, plus focus); "
          f"of them in the core 3,933: {int(sel.isin(core).sum())}")  # fmt: skip
    names = pd.read_csv(REPO / "data/reference/kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "n"]).set_index("ko")["n"]
    for k, g in FOCUS.items():
        if k in tp:
            print(f"  {k} {g:12s} tile prevalence rarefied {tp[k].mean():.3f}, unrarefied {tr[k].mean():.3f}   {str(names.get(k, ''))[:50]}")
    info.to_csv(OUT / "build_info.csv", index=False)


def tech_matrix(plat, logd):
    return np.hstack([plat[:, None].astype(float), SplineTransformer(n_knots=5, degree=3).fit_transform(logd[:, None])])


def pcs(x, tr, te, k=32):
    s = StandardScaler().fit(x[tr])
    p = PCA(n_components=64, random_state=0).fit(s.transform(x[tr]))
    a, b = p.transform(s.transform(x[tr]))[:, :k], p.transform(s.transform(x[te]))[:, :k]
    s2 = StandardScaler().fit(a)
    return s2.transform(a), s2.transform(b)


def logit_loco(xtr, ytr, xte, groups):
    """L2 logistic, C chosen by leave-one-cruise-out log loss inside the training rows."""
    if ytr.min() == ytr.max():
        return np.full(len(xte), ytr.mean())
    loss = np.zeros(len(CS))
    for g in np.unique(groups):
        a, b = groups != g, groups == g
        if ytr[a].min() == ytr[a].max():
            continue
        for i, c in enumerate(CS):
            p = LogisticRegression(C=c, max_iter=500).fit(xtr[a], ytr[a]).predict_proba(xtr[b])[:, 1]
            p = np.clip(p, 1e-6, 1 - 1e-6)
            loss[i] -= (ytr[b] * np.log(p) + (1 - ytr[b]) * np.log(1 - p)).sum()
    return LogisticRegression(C=CS[loss.argmin()], max_iter=500).fit(xtr, ytr).predict_proba(xte)[:, 1]


def logit_loco_k(a, ytr, b, groups, extra=None):
    """As logit_loco, but the number of Granite PCs (KS) is chosen together with C by inner leave-cruise-out.
    extra = (train, test) technical columns appended to every variant."""
    if ytr.min() == ytr.max():
        return np.full(len(b), ytr.mean())
    cat = (lambda x, k, e: np.hstack([x[:, :k], e])) if extra is not None else (lambda x, k, e: x[:, :k])  # noqa: E731
    ea, eb = extra if extra is not None else (None, None)
    loss = np.zeros((len(KS), len(CS)))
    for g in np.unique(groups):
        i, j = groups != g, groups == g
        if ytr[i].min() == ytr[i].max():
            continue
        for ki, k in enumerate(KS):
            xa, xb = cat(a[i], k, None if ea is None else ea[i]), cat(a[j], k, None if ea is None else ea[j])
            for ci, c in enumerate(CS):
                p = np.clip(LogisticRegression(C=c, max_iter=500).fit(xa, ytr[i]).predict_proba(xb)[:, 1], 1e-6, 1 - 1e-6)
                loss[ki, ci] -= (ytr[j] * np.log(p) + (1 - ytr[j]) * np.log(1 - p)).sum()
    ki, ci = np.unravel_index(loss.argmin(), loss.shape)
    k = KS[ki]
    m = LogisticRegression(C=CS[ci], max_iter=500).fit(cat(a, k, ea), ytr)
    return m.predict_proba(cat(b, k, eb))[:, 1], k


def scores(y, p, cr):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    dev = -2 * (y * np.log(p) + (1 - y) * np.log(1 - p)).sum()
    null = np.empty(len(y))
    for c in np.unique(cr):
        null[cr == c] = y[cr != c].mean()
    null = np.clip(null, 1e-6, 1 - 1e-6)
    dev0 = -2 * (y * np.log(null) + (1 - y) * np.log(1 - null)).sum()
    auc = roc_auc_score(y, p) if 0 < y.mean() < 1 else np.nan
    w = [(roc_auc_score(y[cr == c], p[cr == c]), (cr == c).sum()) for c in np.unique(cr) if 0 < y[cr == c].mean() < 1]
    within = np.average([a for a, _ in w], weights=[n for _, n in w]) if w else np.nan
    return auc, 1 - dev / dev0, within


FIXED_C: dict = {}  # model -> C; empty = nested selection


def fit(shard, only=None, name=None):
    tiles = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    emb = pd.DataFrame(np.load(BX / "granite" / "tile_embeddings.npy"), index=tiles.index)
    tp = pd.read_csv(OUT / "presence_tile.csv.gz", index_col=0)
    tr_un = pd.read_csv(OUT / "presence_unrarefied_focus.csv.gz", index_col=0)
    ta = pd.read_csv(OUT / "abundance_tile.csv.gz", index_col=0)
    depth = pd.read_csv(OUT / "tile_depth.csv", index_col=0)["log10_total_genes"]
    idx = tp.index
    t = tiles.loc[idx]
    cr, plat = t["cruise"].to_numpy(), (t["platform_group"] == "NovaSeq 6000 group").to_numpy()
    X = emb.loc[idx].to_numpy()
    tech = tech_matrix(plat, depth.loc[idx].to_numpy())
    kos = list(tp.columns) if only is None else [k for k in only if k in tp.columns]
    if shard:
        kos = kos[shard[0]::shard[1]]
    folds = []
    for c in np.unique(cr):
        te, trn = np.where(cr == c)[0], np.where(cr != c)[0]
        a, b = pcs(X, trn, te)
        st = StandardScaler().fit(tech[trn])
        folds.append((te, trn, a, b, st.transform(tech[trn]), st.transform(tech[te])))
    name = name or (f"fit_shard{shard[0]}.csv" if shard else "fit.csv")
    part = OUT / (name + ".partial")  # one row appended per KO; a rerun skips KOs already there
    done = set(read_partial(part)["ko"]) if part.exists() and part.stat().st_size else set()
    if done:
        print(f"resuming: {len(done)} KOs already in {part.name}", flush=True)
    out = []
    for n, k in enumerate(kos):
        if k in done:
            continue
        y = tp[k].to_numpy()
        r = {"ko": k, "tiles": len(y), "prevalence": y.mean()}
        preds = {m: np.empty(len(y)) for m in ("tech", "granite", "both")}
        kchosen = []
        for te, trn, a, b, sa, sb in folds:
            g = cr[trn]
            preds["tech"][te] = logit_loco(sa, y[trn], sb, g)
            pg = logit_loco_k(a, y[trn], b, g)
            pb = logit_loco_k(a, y[trn], b, g, extra=(sa, sb))
            preds["granite"][te] = pg[0] if isinstance(pg, tuple) else pg
            preds["both"][te] = pb[0] if isinstance(pb, tuple) else pb
            if isinstance(pg, tuple):
                kchosen.append(pg[1])
        r["granite_pcs_median"] = float(np.median(kchosen)) if kchosen else np.nan
        for m, p in preds.items():
            r[f"auc_{m}"], r[f"devexp_{m}"], r[f"auc_within_{m}"] = scores(y, p, cr)
        if k in tr_un:
            yu = tr_un.loc[idx, k].to_numpy()
            pu = np.empty(len(yu))
            for te, trn, a, b, sa, sb in folds:
                pu[te] = logit_loco(sa, yu[trn], sb, cr[trn])
            r["prevalence_unrarefied"] = yu.mean()
            r["auc_tech_unrarefied"], r["devexp_tech_unrarefied"], r["auc_within_tech_unrarefied"] = scores(yu, pu, cr)
        # abundance where present
        m = y == 1
        if m.sum() >= 50 and len(np.unique(cr[m])) >= 4:
            ya = ta.loc[idx[m], [k]].to_numpy()
            Tm = tech_matrix(plat[m], depth.loc[idx[m]].to_numpy())
            Xd = np.hstack([np.ones((m.sum(), 1)), Tm])
            beta, *_ = np.linalg.lstsq(Xd, ya, rcond=None)
            fitv = Xd @ beta
            yadj = ya - (fitv - ya.mean(0))
            crm, Xm = cr[m], X[m]
            pa = np.full_like(yadj, np.nan)
            for c in np.unique(crm):
                te, trn = np.where(crm == c)[0], np.where(crm != c)[0]
                if len(np.unique(crm[trn])) < 2:
                    continue
                a, b = pcs(Xm, trn, te)
                pa[te] = ridge_loco_alpha(a, yadj[trn], b, crm[trn])
            ok = ~np.isnan(pa[:, 0])
            yc = yadj[ok] - yadj[ok].mean()
            r["abund_tiles"] = int(m.sum())
            r["abund_r2"] = 1 - ((yadj[ok] - pa[ok]) ** 2).sum() / (yc**2).sum()
            r["abund_technical_share"] = ((fitv - fitv.mean()) ** 2).sum() / ((ya - ya.mean()) ** 2).sum()
        out.append(r)
        pd.DataFrame([r], columns=FIT_COLS).round(5).to_csv(
            part, mode="a", header=not part.exists() or part.stat().st_size == 0, index=False)
        if n % 50 == 0 or k in FOCUS:
            print(f"  {n:5d}/{len(kos)} {k} {FOCUS.get(k, ''):10s} prev {r['prevalence']:.2f}  AUC tech {r['auc_tech']:.3f} "
                  f"granite {r['auc_granite']:.3f} both {r['auc_both']:.3f} | within-cruise AUC tech {r['auc_within_tech']:.3f} "
                  f"granite {r['auc_within_granite']:.3f} both {r['auc_within_both']:.3f} | abund R2 {r.get('abund_r2', np.nan):.3f}", flush=True)
    read_partial(part).round(5).to_csv(OUT / name, index=False)


def focus_models():
    """Focus KOs: held-out presence probabilities from technical-only, Granite PCs (nested PC count), in-situ
    temperature (spline) and raw satellite (OLCI band means + SST, splines); saved for the figure."""
    tiles = pd.read_csv(BX / "tiles.csv", index_col="tile_key")
    tp = pd.read_csv(OUT / "presence_tile.csv.gz", index_col=0)
    depth = pd.read_csv(OUT / "tile_depth.csv", index_col=0)["log10_total_genes"]
    idx = tp.index
    t = tiles.loc[idx]
    cr, plat = t["cruise"].to_numpy(), (t["platform_group"] == "NovaSeq 6000 group").to_numpy()
    X = pd.DataFrame(np.load(BX / "granite" / "tile_embeddings.npy"), index=tiles.index).loc[idx].to_numpy()
    sat = pd.DataFrame(np.load(BX / "granite" / "raw_satellite.npy"), index=tiles.index).loc[idx].to_numpy()
    sat = np.where(np.isnan(sat), np.nanmedian(sat, 0), sat)
    temp = t[["temperature_c"]].to_numpy()
    tech = tech_matrix(plat, depth.loc[idx].to_numpy())

    def spline_fold(x, tr, te):
        s = StandardScaler().fit(x[tr])
        sp = SplineTransformer(n_knots=5, degree=3, extrapolation="linear").fit(s.transform(x[tr]))
        return sp.transform(s.transform(x[tr])), sp.transform(s.transform(x[te]))

    folds = []
    for c in np.unique(cr):
        te, tr = np.where(cr == c)[0], np.where(cr != c)[0]
        st = StandardScaler().fit(tech[tr])
        folds.append((te, tr, pcs(X, tr, te), (st.transform(tech[tr]), st.transform(tech[te])), spline_fold(temp, tr, te),
                      spline_fold(sat, tr, te)))  # fmt: skip
    rows, preds = [], {}
    for k, gname in FOCUS.items():
        if k not in tp:
            continue
        y = tp[k].to_numpy()
        P = {m: np.empty(len(y)) for m in ("technical", "granite", "temperature", "raw_satellite")}
        for te, tr, (a, b), (sa, sb), (ta_, tb_), (ra, rb) in folds:
            g = cr[tr]
            P["technical"][te] = logit_loco(sa, y[tr], sb, g)
            pg = logit_loco_k(a, y[tr], b, g)
            P["granite"][te] = pg[0] if isinstance(pg, tuple) else pg
            P["temperature"][te] = logit_loco(ta_, y[tr], tb_, g)
            P["raw_satellite"][te] = logit_loco(ra, y[tr], rb, g)
        for m, pr in P.items():
            auc, dev, w = scores(y, pr, cr)
            rows.append({"ko": k, "gene": gname, "model": m, "prevalence": y.mean(), "auc": auc, "devexp": dev, "auc_within": w})
        preds[k] = np.vstack([y] + [P[m] for m in P])
    out = pd.DataFrame(rows)
    out.round(4).to_csv(OUT / "focus_models.csv", index=False)
    np.savez_compressed(OUT / "focus_predictions.npz", tile_key=idx.to_numpy(), cruise=cr, lat=t["sample_lat"].to_numpy(),
                        models=np.array(["observed", "technical", "granite", "temperature", "raw_satellite"]), **preds)  # fmt: skip
    print(out.pivot_table(index=["gene"], columns="model", values=["auc", "auc_within"], sort=False).round(3).to_string())


def summary():
    # Final fit_shard<i>.csv where present, else that shard's .partial (cluster jobs that died on the old final read).
    parts = {p.name.split(".")[0]: p for p in OUT.glob("fit_shard*.csv.partial")}
    parts.update({p.stem: p for p in OUT.glob("fit_shard*.csv")})
    if parts:
        f = pd.concat([read_partial(p) if p.suffix == ".partial" else pd.read_csv(p) for _, p in sorted(parts.items())])
    else:
        f = pd.read_csv(OUT / "fit.csv")
    expect = list(pd.read_csv(OUT / "presence_tile.csv.gz", index_col=0, nrows=0).columns)
    dup = f["ko"].duplicated().sum()
    f = f.drop_duplicates("ko")
    miss = len(set(expect) - set(f["ko"]))
    LOG.append(f"read {len(parts)} shard files: {len(f):,} KOs of {len(expect):,} expected; {miss:,} missing; {dup} duplicates dropped")
    print(LOG[-1])
    if miss:
        print("WARNING: summary is on an incomplete fit", file=sys.stderr)
    names = pd.read_csv(REPO / "data/reference/kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "n"]).set_index("ko")["n"]
    core = set(pd.read_csv(BX / "targets" / "clr_all.csv.gz", nrows=0).columns)
    f["name"] = f["ko"].map(names).fillna("")
    f["in_core"] = f["ko"].isin(core)
    gain = f["auc_within_granite"] - f["auc_within_tech"]
    f["presence_class"] = np.select(
        [(f["auc_within_granite"] >= 0.6) & (gain >= 0.05), f["auc_within_tech"] >= f["auc_within_granite"] - 0.02],
        ["satellite", "technical"], "weak")  # fmt: skip
    f.round(4).to_csv(OUT / "hurdle_all.csv", index=False)
    n = len(f)
    say = lambda s: (print(s), LOG.append(s))  # noqa: E731
    say(f"KOs {n:,}; presence class: " + ", ".join(f"{k} {v:,}" for k, v in f["presence_class"].value_counts().items()))
    for col in ("auc_granite", "auc_tech", "auc_both", "auc_within_granite", "auc_within_tech", "auc_within_both", "abund_r2"):
        say(f"  median {col}: {f[col].median():.3f} (10-90% {f[col].quantile(.1):.3f}-{f[col].quantile(.9):.3f})")
    say(f"  abundance R2 > 0.1: {(f['abund_r2'] > 0.1).sum():,} of {f['abund_r2'].notna().sum():,} KOs with >= 50 present tiles")
    top = f[f["presence_class"] == "satellite"].sort_values("auc_within_granite", ascending=False).head(30)
    say(top[["ko", "name", "prevalence", "auc_within_tech", "auc_within_granite", "auc_granite", "abund_r2"]].round(3).to_string(index=False))
    (OUT / "summary.txt").write_text("\n".join(LOG) + "\n")


def figure():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "pdf.fonttype": 42})  # fmt: skip
    f = pd.read_csv(OUT / "hurdle_all.csv")
    fm = pd.read_csv(OUT / "focus_models.csv")
    z = np.load(OUT / "focus_predictions.npz", allow_pickle=True)
    BLUE, GREY, MAG, ORANGE, LBLUE = "#0072B2", "#8c8c8c", "#c0457e", "#E69F00", "#56B4E9"
    fig = plt.figure(figsize=(7.2, 5.6))
    a = fig.add_axes([0.08, 0.55, 0.36, 0.41])
    b = fig.add_axes([0.60, 0.55, 0.37, 0.41])
    c = fig.add_axes([0.08, 0.08, 0.40, 0.33])
    d = fig.add_axes([0.57, 0.08, 0.40, 0.33])
    col = {"satellite": BLUE, "technical": ORANGE, "weak": "#cccccc"}
    lab = {"satellite": "Satellite signal", "technical": "Technical (depth / platform)", "weak": "Weak / neither"}
    for cls in ("weak", "technical", "satellite"):
        g = f[f["presence_class"] == cls]
        a.scatter(g["auc_within_tech"], g["auc_within_granite"], s=2.5, lw=0, color=col[cls], alpha=0.6,
                  label=f"{lab[cls]} ({len(g):,})", rasterized=True)  # fmt: skip
    a.plot([0.3, 1], [0.3, 1], color="#555555", lw=0.6, ls=(0, (3, 2)))
    for k, gname in FOCUS.items():
        if k in set(f["ko"]) and gname in ("nifH", "nifD", "amoB / pmoB", "nasA", "phnJ"):
            r = f[f["ko"] == k].iloc[0]
            a.scatter(r["auc_within_tech"], r["auc_within_granite"], s=16, facecolor="none", edgecolor="black", lw=0.8, zorder=4)
            a.annotate(gname.split(" /")[0], (r["auc_within_tech"], r["auc_within_granite"]), xytext=(4, 3),
                       textcoords="offset points", fontsize=6.5, style="italic")  # fmt: skip
    a.set_xlim(0.3, 1)
    a.set_ylim(0.3, 1)
    a.set_xlabel("Presence AUC within cruise, depth + platform only")
    a.set_ylabel("Presence AUC within cruise, Granite PCs")
    a.legend(frameon=False, fontsize=5.8, loc="upper left", markerscale=3, handletextpad=0.2)
    genes = [g for g in FOCUS.values() if g in set(fm["gene"])]
    w = fm.pivot_table(index="gene", columns="model", values="auc_within").reindex(genes)
    y = np.arange(len(genes))[::-1]
    for m, cc, mk, lb in (("technical", ORANGE, "s", "Depth + platform"), ("temperature", MAG, "D", "In-situ temperature"),
                          ("raw_satellite", LBLUE, "^", "OLCI + SST"), ("granite", BLUE, "o", "Granite PCs")):  # fmt: skip
        b.scatter(w[m], y, s=18 if m == "granite" else 11, marker=mk, color=cc, lw=0, label=lb, zorder=3)
    b.axvline(0.5, color="#cccccc", lw=0.6)
    b.set_yticks(y)
    b.set_yticklabels(genes, fontsize=6.5, style="italic")
    b.set_xlabel("Presence AUC within cruise")
    b.legend(frameon=False, fontsize=5.8, loc="lower right", handletextpad=0.1)
    models = list(z["models"])
    lat = z["lat"].astype(float)
    bins = np.arange(-75, 85, 10)
    for ax, k, gname in ((c, "K02588", "nifH"), (d, "K10945", "amoB")):
        v = z[k]
        mid = (bins[:-1] + bins[1:]) / 2
        cut = np.digitize(lat, bins) - 1
        for row, cc, lb, ls in ((models.index("observed"), "black", "Observed (share of tiles)", "-"),
                                (models.index("granite"), BLUE, "Granite PCs", "-"),
                                (models.index("temperature"), MAG, "In-situ temperature", "--"),
                                (models.index("technical"), ORANGE, "Depth + platform", ":")):  # fmt: skip
            mval = np.array([v[row][cut == i].mean() if (cut == i).sum() >= 8 else np.nan for i in range(len(mid))])
            ax.plot(mid, mval, color=cc, lw=1.2 if row == 0 else 1.0, ls=ls, marker="o" if row == 0 else None, ms=2.5, label=lb)
        ax.set_xlabel("Latitude (°)")
        ax.set_ylabel(f"$\\it{{{gname}}}$ present (share of tiles / mean probability)")
        ax.set_ylim(0, None)
    c.legend(frameon=False, fontsize=5.8, loc="upper left", handletextpad=0.3)
    for ax, letter in ((a, "a"), (b, "b"), (c, "c"), (d, "d")):
        ax.text(-0.17 if ax is not b else -0.3, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    for suf in ("png", "pdf"):
        fig.savefig(FIG / f"hurdle_kos.{suf}", dpi=300)
    plt.close(fig)


LOG: list = []


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", default=["build"], choices=["build", "fit", "focus", "summary", "figure"])
    ap.add_argument("--depth", type=float, default=1e6)
    ap.add_argument("--min-prev", type=float, default=0.2)
    ap.add_argument("--max-prev", type=float, default=0.9)
    ap.add_argument("--shard", nargs=2, type=int)
    ap.add_argument("--kos", nargs="+", help="fit only these KOs")
    ap.add_argument("--fixed-c", nargs=3, type=float, help="C for tech, granite, both (skips nested selection)")
    ap.add_argument("--out-name")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if "build" in a.parts:
        build(a.depth, a.min_prev, a.max_prev)
    if a.fixed_c:
        FIXED_C.update(dict(zip(["tech", "granite", "both"], a.fixed_c)))
    if "fit" in a.parts:
        fit(a.shard, a.kos, a.out_name)
    if "focus" in a.parts:
        focus_models()
    if "summary" in a.parts:
        summary()
    if "figure" in a.parts:
        figure()
