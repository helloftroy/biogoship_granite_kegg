#!/usr/bin/env python3
"""R73: which combination of inputs predicts C (and each KO) best? Nested search over input groups.

Question (user): R72 showed Granite + 30-day MLD beats OLCI + SST + MLD for C, while putting everything in one ridge
did not help. Could some combination of Granite PCs, OLCI bands, SST, MLD and ocean-model fields do better, per cruise
or overall? And does the same hold for the KOs?

Design. Inputs are split into groups. Every non-empty subset of groups is a candidate model (ridge on the stacked
group features; penalty chosen by leave-cruise-out inside the training cruises, exactly as main_model `cv`).
For each held-out cruise, the subset is chosen from the training cruises only: the leave-cruise-out error inside
the 12 training cruises at the chosen penalty ("inner SSE"). The chosen subset is then scored once on the held-out
cruise. This nested score is the honest one. The "oracle" (the best subset per held-out cruise, chosen on that
cruise) is reported only as an optimistic upper bound.

Groups for C (15 -> 32,767 subsets):
  PC1-4, PC5-8, PC9-16, PC17-32   Granite PCs (PCA 64 on the standardised embedding of the training tiles, as P32)
  OLCI blue (Oa01-05), OLCI green-red (Oa06-08), OLCI red-edge (Oa09-12), OLCI NIR (Oa16-18, Oa21), SST
                                  quantile-normalised on the training tiles, 8-knot splines (as S~)
  MLD30                           30-day mean mixed-layer depth along the back-trajectory (spline, as sat_mld)
  nutrients                       PISCES nitrate, phosphate, silicate, oxygen (log as in R61; splines)
  physics                         GLORYS salinity, temperature, same-day mixed layer, sea level, current speed
  history                         30-day max MLD, 14-day SST change, 14-day wind mean / peak, sea-level gradient,
                                  30-day eddy kinetic energy
  PISCES biology                  chlorophyll, primary production
  time & place                    as the main-model "place" model (linear)
Groups for KOs (7 -> 127 subsets): Granite PCs (32), OLCI (16 bands), SST, MLD30, nutrients + PISCES biology,
physics + history, time & place.
No in-situ (ship) measurement is used. Fields only on the forecast-period tiles (iron, pH, DIC, ...) are left out.

Parts:
  bundle   (local) data/results/combo/inputs.npz from main_model.setup (read-pair technical adjustment, R71)
  fit      --target C|ko --shard i n   per subset: inner SSE and held-out SSE per (outer cruise, column);
           data/results/combo/<target>_shard<i>.npz (written when the shard finishes)
  summary  (local) nested selection, oracle, group frequencies -> combo_<target>_summary.txt, csv tables
Usage:
  ~/miniforge3/bin/python3 scripts/combo_search.py --parts bundle
  sbatch --account=191001-364393 cluster/combo_search.sbatch          (arrays for C and ko)
  ~/miniforge3/bin/python3 scripts/combo_search.py --parts summary --target C   (and --target ko)
"""

import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "results" / "combo"
ALPHAS = np.logspace(-1, 5, 25)  # as screen_ko_predictability.ALPHAS
OLCI = ["Oa01", "Oa02", "Oa03", "Oa04", "Oa05", "Oa06", "Oa07", "Oa08", "Oa09", "Oa10", "Oa11", "Oa12", "Oa16", "Oa17", "Oa18", "Oa21"]
NUTRIENTS = ["nitrate", "phosphate", "silicate", "oxygen"]
PHYSICS = ["salinity", "GLORYS temperature", "mixed layer", "sea level", "current speed"]
HISTORY = ["mixed layer, 30-day max", "SST change, 14 days", "wind, 14-day mean", "wind, 14-day peak", "sea-level gradient",
           "eddy kinetic energy, 30 days"]  # fmt: skip
BIO = ["chlorophyll", "primary production"]
GROUPS = {
    "C": {"PC1-4": ("pc", 0, 4), "PC5-8": ("pc", 4, 8), "PC9-16": ("pc", 8, 16), "PC17-32": ("pc", 16, 32),
          "OLCI blue": ("sat", 0, 5), "OLCI green-red": ("sat", 5, 8), "OLCI red-edge": ("sat", 8, 12),
          "OLCI NIR": ("sat", 12, 16), "SST": ("sat", 16, 17), "MLD30": ("mld",), "nutrients": ("field", NUTRIENTS),
          "physics": ("field", PHYSICS), "history": ("field", HISTORY), "PISCES biology": ("field", BIO),
          "time & place": ("place",)},
    "ko": {"Granite PCs": ("pc", 0, 32), "OLCI": ("sat", 0, 16), "SST": ("sat", 16, 17), "MLD30": ("mld",),
           "nutrients + biology": ("field", NUTRIENTS + BIO), "physics + history": ("field", PHYSICS + HISTORY),
           "time & place": ("place",)},
}  # fmt: skip
REFERENCE = {"C": ["PC1-4", "PC5-8", "PC9-16", "PC17-32"], "ko": ["Granite PCs"]}  # = main-model Granite (P32)


# ----------------------------------------------------------------------------------------------------------- bundle
def bundle():
    sys.path[:0] = [str(REPO), str(REPO / "scripts")]
    import main_model as mm
    from main_fields import HISTORY_FIELDS, FIELDS, VARIABLES
    from model_field_first import HISTORY as HIST
    from screen_ko_predictability import metadata_features
    from water_mass_test import features

    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    fs = features()
    dc, dk = mm.setup("C", fs, ks), mm.setup("clr_all", fs, ks)
    assert (dc["unit"].index == dk["unit"].index).all()
    rows, unit = dk["rows"], dk["unit"]
    fields = pd.read_csv(FIELDS).set_index("sample_id").join(pd.read_csv(HISTORY_FIELDS).set_index("sample_id")
                                                             .drop(columns=["section", "sample_time", "sample_lat", "sample_lon"]))  # fmt: skip
    fields["speed"] = np.hypot(fields["glorys_uo"], fields["glorys_vo"])
    per_tile = lambda col: rows["sample_id"].map(fields[col]).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()  # noqa: E731
    names = NUTRIENTS + PHYSICS + HISTORY + BIO
    F = []
    for n in names:
        col, floor = {**VARIABLES, **HIST}[n]
        v = per_tile(col)
        F.append(np.log10(np.clip(v, floor, None)) if floor is not None else v)
    F = np.column_stack(F)
    miss = np.isnan(F).mean(0)
    print("field missing share:", dict(zip(names, miss.round(4))))
    F = np.where(np.isnan(F), np.nanmedian(F, 0), F)  # few gaps; median (all tiles) as in the R61 per-tile means
    place = np.hstack([metadata_features(unit), np.abs(unit["sample_lat"].to_numpy())[:, None] / 90])
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "inputs.npz", emb=dk["emb"].astype(np.float32), sat=dk["sat"].astype(np.float32),
                        mld=dk["mld"][:, 0].astype(np.float32), fields=F.astype(np.float32), field_names=np.array(names),
                        place=place.astype(np.float32), cruise=dk["cr"], tile=unit.index.to_numpy(),
                        yC=dc["ya"].astype(np.float64), yK=dk["ya"].astype(np.float32), kos=np.array(dk["kos"]))  # fmt: skip
    print(f"inputs.npz: tiles {len(unit):,}, KOs {len(dk['kos']):,}, fields {F.shape[1]}, sat {dk['sat'].shape[1]}, "
          f"place {place.shape[1]}; missing field values imputed {int(np.isnan(np.column_stack([per_tile({**VARIABLES, **HIST}[n][0]) for n in names])).sum())}")  # fmt: skip


# -------------------------------------------------------------------------------------------------------------- fit
def ridge_path(xtr, ytr, xte):
    xm, ym = xtr.mean(axis=0), ytr.mean(axis=0)
    u, s, vt = np.linalg.svd(xtr - xm, full_matrices=False)
    uty = u.T @ (ytr - ym)
    proj = (xte - xm) @ vt.T
    shrink = s[None, :] / (s[None, :] ** 2 + ALPHAS[:, None])
    return ym + np.einsum("tc,ac,ck->atk", proj, shrink, uty)


def ridge_cv(xtr, ytr, xte, groups):
    """main_model.ridge(..., 'cv'): per-column penalty by leave-cruise-out in the training rows.
    Returns held-out predictions and the inner leave-cruise-out SSE at the chosen penalty (per column)."""
    sse = np.zeros((len(ALPHAS), ytr.shape[1]))
    for g in np.unique(groups):
        a, b = groups != g, groups == g
        sse += ((ytr[b][None] - ridge_path(xtr[a], ytr[a], xtr[b])) ** 2).sum(axis=1)
    idx = sse.argmin(0)
    path = ridge_path(xtr, ytr, xte)
    return path[idx, :, np.arange(ytr.shape[1])].T, sse[idx, np.arange(ytr.shape[1])]


def group_features(z, tr, te, groups):
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import QuantileTransformer, SplineTransformer, StandardScaler

    def scale(a, b):
        s = StandardScaler().fit(a)
        return s.transform(a), s.transform(b)

    def spline(a, b):
        s = SplineTransformer(n_knots=8, degree=3, extrapolation="linear").fit(a)
        return s.transform(a), s.transform(b)

    e_tr, e_te = scale(z["emb"][tr], z["emb"][te])
    pca = PCA(n_components=64, random_state=0).fit(e_tr)
    p_tr, p_te = pca.transform(e_tr), pca.transform(e_te)
    sat = z["sat"].astype(float)
    med = np.nanmedian(sat[tr], axis=0)
    s_tr, s_te = np.where(np.isnan(sat[tr]), med, sat[tr]), np.where(np.isnan(sat[te]), med, sat[te])
    q = QuantileTransformer(n_quantiles=min(1000, len(tr)), output_distribution="normal", random_state=0).fit(s_tr)
    s_tr, s_te = scale(q.transform(s_tr), q.transform(s_te))
    fn = list(z["field_names"])
    out = {}
    for name, spec in groups.items():
        if spec[0] == "pc":
            out[name] = scale(p_tr[:, spec[1]:spec[2]], p_te[:, spec[1]:spec[2]])  # P32 scales its 32 columns the same way
        elif spec[0] == "sat":
            out[name] = spline(s_tr[:, spec[1]:spec[2]], s_te[:, spec[1]:spec[2]])
        elif spec[0] == "mld":
            a, b = scale(z["mld"][tr, None].astype(float), z["mld"][te, None].astype(float))
            out[name] = spline(a, b)
        elif spec[0] == "field":
            cols = [fn.index(v) for v in spec[1]]
            a, b = scale(z["fields"][tr][:, cols].astype(float), z["fields"][te][:, cols].astype(float))
            out[name] = spline(a, b)
        elif spec[0] == "place":
            out[name] = scale(z["place"][tr].astype(float), z["place"][te].astype(float))
    return out


def subsets(target):
    names = list(GROUPS[target])
    return [c for k in range(1, len(names) + 1) for c in itertools.combinations(range(len(names)), k)]


def fit(target, shard):
    z = np.load(OUT / "inputs.npz", allow_pickle=True)
    y = z["yC"] if target == "C" else z["yK"].astype(np.float64)
    cr = z["cruise"]
    names = list(GROUPS[target])
    subs = subsets(target)
    mine = list(range(shard[0], len(subs), shard[1]))
    part = OUT / f"{target}_shard{shard[0]}.npz"
    cruises = np.unique(cr)
    t0 = time.time()
    folds = []
    for c in cruises:
        te, tr = np.where(cr == c)[0], np.where(cr != c)[0]
        folds.append((te, tr, group_features(z, tr, te, GROUPS[target])))
    print(f"{target} shard {shard[0]}/{shard[1]}: {len(mine):,} of {len(subs):,} subsets; features built in {time.time() - t0:.0f} s", flush=True)
    inner = np.full((len(mine), len(cruises), y.shape[1]), np.nan, np.float32)
    held = np.full_like(inner, np.nan)
    for j, s in enumerate(mine):
        for i, (te, tr, feats) in enumerate(folds):
            xtr = np.hstack([feats[names[g]][0] for g in subs[s]])
            xte = np.hstack([feats[names[g]][1] for g in subs[s]])
            p, sse_in = ridge_cv(xtr, y[tr], xte, cr[tr])
            inner[j, i], held[j, i] = sse_in, ((y[te] - p) ** 2).sum(0)
        if j % 200 == 0 or j == len(mine) - 1:
            print(f"  {j + 1:6,}/{len(mine):,} subsets  ({(time.time() - t0) / (j + 1):.2f} s/subset)", flush=True)
    np.savez_compressed(part, subsets=np.array(mine), inner=inner, held=held, cruises=cruises.astype(str))
    print(f"wrote {part.name} ({time.time() - t0:.0f} s)")


# ---------------------------------------------------------------------------------------------------------- summary
def summary(target):
    z = np.load(OUT / "inputs.npz", allow_pickle=True)
    y = z["yC"] if target == "C" else z["yK"].astype(np.float64)
    cr = z["cruise"]
    names = list(GROUPS[target])
    subs = subsets(target)
    parts = sorted(OUT.glob(f"{target}_shard*.npz"))
    idx = np.concatenate([np.load(p)["subsets"] for p in parts])
    inner = np.concatenate([np.load(p)["inner"] for p in parts]).astype(np.float64)
    held = np.concatenate([np.load(p)["held"] for p in parts]).astype(np.float64)
    cruises = np.unique(cr)  # the fold order used in fit (shards store it too; read it from inputs.npz to avoid
    # unpickling object arrays written by a different numpy version on the cluster)
    order = np.argsort(idx)
    idx, inner, held = idx[order], inner[order], held[order]
    log = [f"{target}: subsets {len(idx):,} of {len(subs):,}; shards {len(parts)}; columns {y.shape[1]:,}"]
    assert len(idx) == len(subs), "missing shards"
    ss = ((y - y.mean(0)) ** 2).sum(0)  # per column, pooled over tiles (main-model 'explained')
    ref = subs.index(tuple(names.index(g) for g in REFERENCE[target]))
    lab = lambda s: " + ".join(names[g] for g in subs[s])  # noqa: E731

    def pooled(sse_cols):  # explained, % of non-technical, pooled over columns
        return 100 * (1 - sse_cols.sum() / ss.sum())

    # nested: per column, per outer cruise, the subset with the lowest inner SSE
    choice = inner.argmin(0)  # cruises x columns
    nested_sse = np.take_along_axis(held, choice[None], 0)[0]  # cruises x columns
    ref_sse = held[ref]
    oracle_sse = held.min(0)
    log.append(f"\npooled explained, % of non-technical variance:")
    log.append(f"  Granite PCs (reference, = main model P32): {pooled(ref_sse.sum(0)):.2f}")
    log.append(f"  nested selection (chosen in the training cruises): {pooled(nested_sse.sum(0)):.2f}")
    log.append(f"  oracle, best subset per held-out cruise and column (optimistic): {pooled(oracle_sse.sum(0)):.2f}")
    fixed = held.sum(1)  # subsets x columns
    best_fixed = fixed.argmin(0)
    log.append(f"  best single subset per column chosen on held-out error (optimistic): {pooled(fixed.min(0)):.2f}")
    if target == "C":
        tot = fixed[:, 0]
        top = np.argsort(tot)[:15]
        log.append("\ntop 15 subsets by held-out error over all cruises (optimistic ranking; explained %):")
        for s in top:
            log.append(f"  {pooled(np.array([tot[s]])):6.2f}  {lab(idx[s])}")
        ins = inner.sum(1)[:, 0]
        log.append("\ntop 10 subsets by inner (training-cruise) error summed over outer folds (no held-out information):")
        for s in np.argsort(ins)[:10]:
            log.append(f"  inner rank; held-out explained {pooled(np.array([tot[s]])):6.2f}  {lab(idx[s])}")
        log.append("\nper held-out cruise: chosen subset (from the training cruises) and R2 vs the training-cruise mean")
        rows = []
        for i, c in enumerate(cruises):
            k = cr == c
            mu = y[~k].mean()
            den = ((y[k] - mu) ** 2).sum()
            rows.append({"cruise": c, "tiles": int(k.sum()), "R2_granite": 1 - held[ref, i, 0] / den,
                         "R2_nested": 1 - nested_sse[i, 0] / den, "R2_oracle": 1 - oracle_sse[i, 0] / den,
                         "chosen": lab(idx[choice[i, 0]]), "oracle_subset": lab(idx[held[:, i, 0].argmin()])})  # fmt: skip
        pc = pd.DataFrame(rows)
        pc.round(3).to_csv(OUT / "combo_C_per_cruise.csv", index=False)
        log.append(pc.round(2).to_string(index=False))
        # group usage among the 100 best subsets by inner error, per outer fold
        use = np.zeros(len(names))
        for i in range(len(cruises)):
            for s in np.argsort(inner[:, i, 0])[:100]:
                for g in subs[idx[s]]:
                    use[g] += 1
        log.append("\nshare of the 100 best subsets (inner error, per outer fold) that contain each group:")
        for g, u in sorted(zip(names, use / (100 * len(cruises))), key=lambda x: -x[1]):
            log.append(f"  {g:16s} {u:.2f}")
    else:
        kos = z["kos"]
        r2 = lambda sse: 1 - sse / ss  # noqa: E731
        per = pd.DataFrame({"ko": kos, "r2_granite": r2(ref_sse.sum(0)), "r2_nested": r2(nested_sse.sum(0)),
                            "r2_oracle": r2(oracle_sse.sum(0)), "best_fixed_subset": [lab(idx[s]) for s in best_fixed],
                            "r2_best_fixed": r2(fixed.min(0))})  # fmt: skip
        per["gain_nested"] = per["r2_nested"] - per["r2_granite"]
        better = (nested_sse < ref_sse).sum(0)
        per["cruises_nested_better"] = better
        per.round(4).to_csv(OUT / "combo_ko_per_ko.csv", index=False)
        g = per["gain_nested"]
        log.append(f"\nper-KO median R2: Granite {per['r2_granite'].median():.4f}, nested {per['r2_nested'].median():.4f}, "
                   f"oracle {per['r2_oracle'].median():.4f}")
        log.append(f"nested - Granite per KO: median {g.median():+.4f} (10-90% {g.quantile(.1):+.4f} to {g.quantile(.9):+.4f}); "
                   f">= +0.05: {(g >= .05).sum():,}; >= +0.05 and better in >= 10 of 13 cruises: "
                   f"{((g >= .05) & (better >= 10)).sum():,}; <= -0.05: {(g <= -.05).sum():,}")  # fmt: skip
        # one subset shared by all KOs: chosen per outer fold by inner SSE summed over KOs (127 options per fold)
        shared = inner.sum(2).argmin(0)  # per outer fold
        sh_sse = np.stack([held[shared[i], i] for i in range(len(cruises))])
        log.append(f"\none subset for all KOs, chosen in the training cruises (nested): {pooled(sh_sse.sum(0)):.2f}; "
                   f"chosen per fold: " + "; ".join(f"{c}: {lab(idx[shared[i]])}" for i, c in enumerate(cruises)))  # fmt: skip
        tab = pd.DataFrame({"subset": [lab(s) for s in idx], "explained_heldout": [pooled(fixed[j]) for j in range(len(idx))],
                            "explained_inner_rank": pd.Series(inner.sum((1, 2))).rank().to_numpy()})  # fmt: skip
        tab["groups"] = [len(subs[s]) for s in idx]
        tab = tab.sort_values("explained_heldout", ascending=False)
        tab.round(3).to_csv(OUT / "combo_ko_subsets.csv", index=False)
        log.append("\nall-KO pooled explained per subset, held out (top 15; one model for every KO; optimistic ranking over 127):")
        log.append(tab.head(15).round(2).to_string(index=False))
        for nm in ("Granite PCs", "Granite PCs + MLD30", "Granite PCs + OLCI", "Granite PCs + SST + MLD30",
                   "Granite PCs + nutrients + biology", "Granite PCs + SST + MLD30 + nutrients + biology",
                   "OLCI + SST", "OLCI + SST + MLD30", "SST + MLD30 + nutrients + biology"):
            if nm in set(tab["subset"]):
                log.append(f"  {nm:50s} {tab.loc[tab['subset'] == nm, 'explained_heldout'].iloc[0]:6.2f}")
        names_ko = pd.read_csv(REPO / "data/reference/kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "name"]).set_index("ko")["name"]
        win = per[(per["gain_nested"] >= .05) & (per["cruises_nested_better"] >= 10)].sort_values("gain_nested", ascending=False)
        log.append(f"\nKOs gaining >= 0.05 with >= 10 of 13 cruises better (nested; {len(win)}), top 25:")
        for _, r in win.head(25).iterrows():
            log.append(f"  {r['ko']}  {r['r2_granite']:+.3f} -> {r['r2_nested']:+.3f} ({r['gain_nested']:+.3f}, {int(r['cruises_nested_better'])} cruises)  "
                       f"{str(names_ko.get(r['ko'], ''))[:60]}")  # fmt: skip
        use = np.zeros(len(names))
        for s in choice.ravel():
            for gg in subs[idx[s]]:
                use[gg] += 1
        log.append("\nshare of (KO, outer fold) choices that include each group:")
        for nm, u in sorted(zip(names, use / choice.size), key=lambda x: -x[1]):
            log.append(f"  {nm:22s} {u:.2f}")
        top = pd.Series([lab(idx[s]) for s in choice.ravel()]).value_counts().head(12) / choice.size
        log.append("\nmost often chosen subsets (share of KO x fold choices):")
        log.append(top.round(3).to_string())
    (OUT / f"combo_{target}_summary.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", nargs="+", required=True, choices=["bundle", "fit", "summary"])
    ap.add_argument("--target", choices=["C", "ko"])
    ap.add_argument("--shard", nargs=2, type=int)
    a = ap.parse_args()
    if "bundle" in a.parts:
        bundle()
    if "fit" in a.parts:
        fit(a.target, a.shard or (0, 1))
    if "summary" in a.parts:
        summary(a.target)
