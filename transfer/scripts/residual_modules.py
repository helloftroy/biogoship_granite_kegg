#!/usr/bin/env python3
"""What the satellite misses: shared structure in the cross-validated residuals (RESULTS.md R33).

1. Baseline: 32 Granite PCs + in-situ temperature spline, ridge, one cruise held out
   (R18; the model of R28), on every target KO. Residual = observed CLR - prediction.
2. Residual variance split into between-cruise (per-cruise mean offsets, which a
   held-out-cruise model cannot know) and within-cruise parts.
3. PCA of the residuals (each KO standardised, so a module is KOs that co-vary, whatever
   their size). Components kept: eigenvalue above the 95th percentile of the same rank
   in 20 null PCAs where each KO's residuals are permuted across tiles independently.
   Run on the full residuals and on the within-cruise residuals (cruise means removed).
   PCA, not NMF: residuals are signed, and a module's two poles are both of interest.
4. Each component is described by:
   - cruise eta2: share of its score variance between cruises (offset-like: batch,
     region or year)
   - correlations with candidate drivers, Spearman, over all tiles and within cruise
     (ranks demeaned per cruise), and in how many cruises the within-cruise sign holds:
       sequencing  log total counts, libraries, instrument, year 2021 on, share unclassified
       sampling    depth, year, season, local solar hour (diel), latitude
       in situ     phosphate, nitrate (manifest), CCHDO bottle salinity/oxygen/nutrients/carbonate
       model       GLORYS / PISCES day fields and history (R25, R27), transformed as there
       taxa        domain CLRs, unclassified logit, best phylum CLRs
   - KEGG: KOs with the largest loadings; BRITE level B and C (pathway) categories whose
     mean loading departs most from the rest (z of the mean).

Outputs in data/results/residual_modules/: components.csv, drivers.csv, categories.csv,
top_kos.csv, scores.csv (tile x component), loadings.csv.gz (KO x component), summary.txt.

Usage:
    OMP_NUM_THREADS=4 python scripts/residual_modules.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biogoship.tables import REPO_ROOT  # noqa: E402
from added_fields_all_kos import cv  # noqa: E402
from model_field_first import FIELDS, HISTORY, HISTORY_FIELDS, VARIABLES  # noqa: E402
from tune_figure_kos import load  # noqa: E402

RESULTS = REPO_ROOT / "data" / "results"
TARGETS = REPO_ROOT / "data" / "targets"
PATHWAY_A = {"Metabolism", "Genetic Information Processing", "Environmental Information Processing", "Cellular Processes"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", type=Path, default=TARGETS / "ko_clr.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings_sst_l4", REPO_ROOT / "data" / "faire" / "embeddings_sst_l4"])
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "residual_modules")
    parser.add_argument("--components", type=int, default=10, help="components described (at least)")
    parser.add_argument("--null-draws", type=int, default=20)
    parser.add_argument("--splits", type=int, default=20, help="random halvings of the cruises")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, _, y, cruises, unit, rows = load(args, kos, return_rows=True)
    n = len(unit)
    say(f"{n:,} tiles, {len(kos):,} KOs, {len(set(cruises))} cruises")
    r2, pred = cv(emb, temp[:, 0], None, y, cruises, np.arange(n))
    res = y - pred
    say(f"baseline (32 PCs + temperature spline, cruise held out): median KO R2 {np.median(r2):.3f}")

    # 2. between / within cruise
    cm = pd.DataFrame(res).groupby(cruises).transform("mean").to_numpy()
    within = res - cm
    ss_tot = ((res - res.mean(0)) ** 2).sum()
    ss_between = ((cm - res.mean(0)) ** 2).sum()
    obs_tot = ((y - y.mean(0)) ** 2).sum()
    say(f"residual variance = {((res - res.mean(0)) ** 2).sum() / obs_tot:.1%} of the KO CLR variance; "
        f"of it, between cruises {ss_between / ss_tot:.1%}, within cruises {1 - ss_between / ss_tot:.1%}")

    # 3. PCA with permutation null
    drivers, phyla = driver_table(unit, rows)
    brite = pd.read_csv(REPO_ROOT / "data" / "reference" / "kegg_brite_ko00001.csv")
    names = pd.read_csv(REPO_ROOT / "data" / "reference" / "kegg_ko_list.tsv", sep="\t", header=None, index_col=0)[1]
    ko_ids = [k.replace("count.ko.", "").replace("ko:", "") for k in kos]
    comp_rows, drv_rows, cat_rows, top_rows, score_cols, load_cols = [], [], [], [], {}, {}
    for variant, mat in (("full", res), ("within cruise", within)):
        z = standardise(mat)
        ev, vec = top_pca(z)
        null = np.array([top_pca(standardise(np.column_stack([rng.permutation(c) for c in mat.T])))[0][:30]
                         for _ in range(args.null_draws)])
        null95 = np.percentile(null, 95, axis=0)
        k_sig = int(np.argmax(ev[:30] <= null95)) if (ev[:30] <= null95).any() else 30
        total = ev.sum()
        K = args.components
        rep_med, rep_share = split_half(mat, cruises, vec[:, :K], rng, args.splits)
        say(f"\n== {variant} residuals: {k_sig} of the first 30 components above the permutation null "
            f"(each KO shuffled on its own: a weak bar). Split-half reproducibility (cruises split in two, "
            f"{args.splits} splits; congruence of the matched loadings, >= 0.85 = same component):")
        say("   " + "  ".join(f"{j + 1}: {rep_med[j]:.2f}" for j in range(K)))
        for j in range(K):
            score = z @ vec[:, j] / np.sqrt(ev[j])  # unit variance
            loading = vec[:, j] * np.sqrt(ev[j] / (len(z) - 1))  # correlation of KO residual with score
            if loading.sum() < 0:  # sign: most KOs load positively
                score, loading = -score, -loading
            name = f"{'F' if variant == 'full' else 'W'}{j + 1}"
            eta2 = pd.Series(score).groupby(cruises).transform("mean").var() / score.var()
            pr = group_r2(score, cruises, drivers, phyla, within=variant != "full")
            comp_rows.append({"component": name, "residuals": variant, "rank": j + 1, "share": ev[j] / total,
                              "null95_share": null95[j] / total if j < 30 else np.nan, "above_null": j < k_sig,
                              "split_half_congruence": rep_med[j], "share_splits_0.85": rep_share[j],
                              "cruise_eta2": eta2, "share_kos_positive": (loading > 0).mean(),
                              **{f"r2 {g}": v for g, v in pr.items()}})
            score_cols[name], load_cols[name] = score, loading
            d = describe_drivers(score, cruises, drivers, phyla)
            drv_rows.append(d.assign(component=name))
            c = categories(loading, ko_ids, brite)
            cat_rows.append(c.assign(component=name))
            order = np.argsort(loading)
            for side, sel in (("+", order[::-1][:10]), ("-", order[:10])):
                for i in sel:
                    top_rows.append({"component": name, "side": side, "ko": ko_ids[i], "loading": loading[i],
                                     "name": str(names.get(ko_ids[i], ""))[:70]})
            best = d.reindex(d["within_r"].abs().sort_values(ascending=False).index).head(4)
            say(f"{name}: {ev[j] / total:5.1%} of residual variance  split-half {rep_med[j]:.2f}  cruise eta2 {eta2:.2f}  "
                f"KOs + {(loading > 0).mean():.0%}")
            say("      CV R2 of the score from: " + "  ".join(f"{g} {v:+.2f}" for g, v in pr.items()))
            for _, b in best.iterrows():
                say(f"      {b['driver']:38s} all r {b['all_r']:+.2f}  within r {b['within_r']:+.2f}  "
                    f"(n {int(b['n'])}, same sign in {int(b['cruises_same_sign'])}/{int(b['cruises_tested'])} cruises)")
            pos = c.sort_values("z", ascending=False).head(2)
            neg = c.sort_values("z").head(2)
            say("      KEGG +: " + "; ".join(f"{r.category} (z {r.z:+.1f})" for r in pos.itertuples())
                + " | -: " + "; ".join(f"{r.category} (z {r.z:+.1f})" for r in neg.itertuples()))

    out = args.out_dir
    pd.DataFrame(comp_rows).round(4).to_csv(out / "components.csv", index=False)
    pd.concat(drv_rows).round(4).to_csv(out / "drivers.csv", index=False)
    pd.concat(cat_rows).round(4).to_csv(out / "categories.csv", index=False)
    pd.DataFrame(top_rows).round(4).to_csv(out / "top_kos.csv", index=False)
    pd.DataFrame(score_cols, index=unit.index).assign(section=cruises, sample_id=unit["sample_id"].to_numpy()).round(4) \
        .to_csv(out / "scores.csv")
    pd.DataFrame(load_cols, index=ko_ids).round(4).to_csv(out / "loadings.csv.gz")
    (out / "summary.txt").write_text("\n".join(log) + "\n")


def standardise(m: np.ndarray) -> np.ndarray:
    m = m - m.mean(0)
    sd = m.std(0)
    return m / np.where(sd > 0, sd, 1)


def top_pca(z: np.ndarray):
    """Eigenvalues (descending) and KO-space eigenvectors of z'z, through the tile Gram matrix."""
    g = z @ z.T
    ev, u = np.linalg.eigh(g)
    ev, u = ev[::-1], u[:, ::-1]
    keep = ev > 1e-8 * ev[0]
    ev, u = ev[keep], u[:, keep]
    vec = z.T @ u[:, :60] / np.sqrt(ev[:60])
    return ev, vec


def split_half(mat, cruises, full_vec, rng, n_splits):
    """Per full-data component: median congruence of its best matches in two halves of the cruises."""
    names = np.unique(cruises)
    cong = np.zeros((n_splits, full_vec.shape[1]))
    for s in range(n_splits):
        half = set(rng.permutation(names)[: len(names) // 2])
        a = np.array([c in half for c in cruises])
        va, vb = top_pca(standardise(mat[a]))[1][:, :30], top_pca(standardise(mat[~a]))[1][:, :30]
        for j in range(full_vec.shape[1]):
            ia = np.abs(va.T @ full_vec[:, j]).argmax()
            ib = np.abs(vb.T @ full_vec[:, j]).argmax()
            cong[s, j] = abs(va[:, ia] @ vb[:, ib])
    return np.median(cong, axis=0), (cong >= 0.85).mean(axis=0)


GROUPS = ["sequencing", "sampling", "GLORYS+history", "all environment", "taxa (carriers)",
          "GLORYS+history (PISCES tiles)", "+ PISCES (PISCES tiles)"]


def group_r2(score, cruises, drivers: pd.DataFrame, phyla: pd.DataFrame, within: bool = False) -> dict:
    """Pooled R2 of the component score from each group of candidate drivers; linear ridge, one cruise held out.

    For within-cruise components the drivers are centred per cruise too (a within-cruise regression)."""
    cover = drivers.notna().mean()
    model_all = [c for c in drivers if c.startswith("model:")]
    full_model = [c for c in model_all if cover[c] > 0.9]
    seq = [c for c in drivers if c.startswith("sequencing:")]
    samp = [c for c in drivers if c.startswith("sampling:")]
    taxa = drivers[[c for c in drivers if c.startswith("taxa:")]].join(phyla)
    pisces_tiles = drivers[model_all].notna().all(axis=1).to_numpy()
    sets = {"sequencing": (drivers[seq], None), "sampling": (drivers[samp], None),
            "GLORYS+history": (drivers[full_model], None), "all environment": (drivers[seq + samp + full_model], None),
            "taxa (carriers)": (taxa, None), "GLORYS+history (PISCES tiles)": (drivers[full_model], pisces_tiles),
            "+ PISCES (PISCES tiles)": (drivers[model_all], pisces_tiles)}  # fmt: skip
    out = {}
    for g in GROUPS:
        x, keep = sets[g]
        idx = np.flatnonzero(keep) if keep is not None else np.arange(len(score))
        xv = x.to_numpy()[idx]
        if within:
            xv = xv - pd.DataFrame(xv).groupby(cruises[idx]).transform("mean").to_numpy()
        out[g] = loco_ridge(xv, score[idx], cruises[idx])
    return out


def loco_ridge(x, s, groups):
    pred = np.full(len(s), np.nan)
    for g in np.unique(groups):
        te, tr = groups == g, groups != g
        med = np.nanmedian(x[tr], axis=0)
        xtr, xte = np.where(np.isnan(x[tr]), med, x[tr]), np.where(np.isnan(x[te]), med, x[te])
        sc = StandardScaler().fit(xtr)
        m = RidgeCV(alphas=np.logspace(-2, 4, 25)).fit(sc.transform(xtr), s[tr])
        pred[te] = m.predict(sc.transform(xte))
    return 1 - ((s - pred) ** 2).sum() / ((s - s.mean()) ** 2).sum()


def tile_mean(rows: pd.DataFrame, unit: pd.DataFrame, by_sample: pd.Series) -> np.ndarray:
    return rows["sample_id"].map(by_sample).astype(float).groupby(rows["tile_key"]).mean().reindex(unit.index).to_numpy()


def driver_table(unit: pd.DataFrame, rows: pd.DataFrame):
    d = {}
    fields = pd.read_csv(FIELDS).set_index("sample_id").join(
        pd.read_csv(HISTORY_FIELDS).set_index("sample_id").drop(columns=["section", "sample_time", "sample_lat", "sample_lon"]))
    fields["speed"] = np.hypot(fields["glorys_uo"], fields["glorys_vo"])
    t = pd.to_datetime(fields["sample_time"], utc=True, errors="coerce")
    lat, lon = fields["sample_lat"], fields["sample_lon"]
    solar = (t.dt.hour + t.dt.minute / 60 + lon / 15) % 24
    doy = t.dt.dayofyear
    hemi = np.sign(lat).replace(0, 1)
    man = pd.read_csv(REPO_ROOT / "data" / "manifests" / "biogoship_samples.csv").drop_duplicates("sample_id").set_index("sample_id")
    ks = pd.read_csv(TARGETS / "ko_samples.csv").set_index("sample_key")
    for name, v in {
        "sampling: depth (m)": ks["depth_m"], "sampling: year": t.dt.year + doy / 365.25,
        "sampling: season (cos, local summer +)": np.cos(2 * np.pi * (doy - 172) / 365.25) * hemi,
        "sampling: season (sin)": np.sin(2 * np.pi * (doy - 172) / 365.25) * hemi,
        "sampling: solar hour (cos, noon +)": np.cos(2 * np.pi * (solar - 12) / 24),
        "sampling: solar hour (sin, morning +)": -np.sin(2 * np.pi * (solar - 12) / 24),
        "sampling: latitude": lat, "sampling: |latitude|": lat.abs(),
        "sampling: longitude (cos)": np.cos(np.radians(lon)), "sampling: longitude (sin)": np.sin(np.radians(lon)),
        "in situ: temperature": man["temperature_c"],
        "in situ: log phosphate": np.log10(man["phosphate_umol_kg"].clip(lower=0.005)),
        "in situ: log nitrate+nitrite": np.log10(man["nitrate_nitrite_umol_l"].clip(lower=0.01)),
    }.items():
        d[name] = tile_mean(rows, unit, v)
    cch = pd.read_csv(REPO_ROOT / "data" / "manifests" / "biogoship_samples_cchdo.csv").drop_duplicates("sample_id").set_index("sample_id")
    for col, (label, logf) in {"cchdo_salinity_pss": ("salinity", None), "cchdo_oxygen_umol_kg": ("oxygen", None),
                               "cchdo_silicate_umol_kg": ("log silicate", 0.05), "cchdo_nitrite_umol_kg": ("log nitrite", 0.001),
                               "cchdo_dic_umol_kg": ("DIC", None), "cchdo_alkalinity_umol_kg": ("alkalinity", None),
                               "cchdo_ph_total": ("pH", None)}.items():  # fmt: skip
        v = pd.to_numeric(cch[col], errors="coerce")
        d[f"in situ (CCHDO): {label}"] = tile_mean(rows, unit, np.log10(v.clip(lower=logf)) if logf else v)
    for name, (col, floor) in {**VARIABLES, **HISTORY}.items():
        v = fields[col]
        d[f"model: {name}"] = tile_mean(rows, unit, np.log10(v.clip(lower=floor)) if floor is not None else v)
    d["sequencing: log total counts"] = tile_mean(rows, unit, np.log10(ks["total"]))
    d["sequencing: libraries"] = tile_mean(rows, unit, ks["libraries"])
    inst = pd.read_csv(TARGETS / "taxa" / "instrument_target.csv").set_index("sample_key")
    d["sequencing: instrument 3310"] = tile_mean(rows, unit, inst["instrument_3310"])
    d["sequencing: year 2021 on"] = tile_mean(rows, unit, inst["year_2021_on"])
    unc = pd.read_csv(TARGETS / "taxa" / "unclassified_targets.csv").set_index("sample_key")
    d["taxa: unclassified (logit)"] = tile_mean(rows, unit, unc["unclassified_logit"])
    dom = pd.read_csv(TARGETS / "taxa" / "taxa_domain_clr.csv").set_index("sample_key")
    for c in dom.columns:
        d[f"taxa: {c} (CLR)"] = tile_mean(rows, unit, dom[c])
    phy = pd.read_csv(TARGETS / "taxa" / "taxa_phylum_clr.csv").set_index("sample_key")
    phyla = pd.DataFrame({c: tile_mean(rows, unit, phy[c]) for c in phy.columns}, index=unit.index)
    return pd.DataFrame(d, index=unit.index), phyla


def within_rank(x: np.ndarray, groups: np.ndarray) -> np.ndarray:
    r = pd.Series(x).rank().to_numpy()
    return r - pd.Series(r).groupby(groups).transform("mean").to_numpy()


def corr(a, b):
    a, b = a - a.mean(), b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return (a * b).sum() / den if den > 0 else np.nan


def describe_drivers(score, cruises, drivers: pd.DataFrame, phyla: pd.DataFrame) -> pd.DataFrame:
    out = []
    for name, x in drivers.items():
        out.append(one_driver(name, score, x.to_numpy(), cruises))
    ph = [one_driver(f"taxa: phylum {c}", score, x.to_numpy(), cruises) for c, x in phyla.items()]
    ph = sorted([p for p in ph if p], key=lambda p: -abs(p["within_r"]))[:3]
    return pd.DataFrame([o for o in out if o] + ph)


def one_driver(name, score, x, cruises):
    ok = ~np.isnan(x)
    if ok.sum() < 100 or len(set(cruises[ok])) < 3 or np.nanstd(x) == 0:
        return None
    s, v, g = score[ok], x[ok], cruises[ok]
    all_r = corr(pd.Series(s).rank().to_numpy(), pd.Series(v).rank().to_numpy())
    w_r = corr(within_rank(s, g), within_rank(v, g))
    same, tested = 0, 0
    for c in np.unique(g):
        m = g == c
        if m.sum() >= 15 and np.std(v[m]) > 0:
            rc = corr(pd.Series(s[m]).rank().to_numpy(), pd.Series(v[m]).rank().to_numpy())
            tested += 1
            same += int(np.sign(rc) == np.sign(w_r))
    return {"driver": name, "n": int(ok.sum()), "cruises": len(set(g)), "all_r": all_r, "within_r": w_r,
            "cruises_tested": tested, "cruises_same_sign": same}


def categories(loading, ko_ids, brite: pd.DataFrame) -> pd.DataFrame:
    lo = pd.Series(loading, index=ko_ids)
    b = brite[brite["A"].isin(PATHWAY_A) & brite["ko"].isin(lo.index)]
    mu, sd = lo.mean(), lo.std()
    out = []
    for level in ("B", "C"):
        for cat, kk in b.groupby(level)["ko"]:
            kk = kk.drop_duplicates()
            if len(kk) < 10:
                continue
            m = lo[kk].mean()
            out.append({"level": level, "category": cat, "kos": len(kk), "mean_loading": m,
                        "z": (m - mu) / (sd / np.sqrt(len(kk)))})
    return pd.DataFrame(out)


if __name__ == "__main__":
    main()
