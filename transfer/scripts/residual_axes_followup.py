#!/usr/bin/env python3
"""Follow-up tests of the R33 residual axes (RESULTS.md R34).

Parts (run one or several with --parts):
  protocol  Do sampling / extraction / sequencing protocol differences explain the axes?
            a) cruise level: per-cruise mean of each F score against the sequencing
               flowcell, sequencing year, sampling year, storage time (sequencing date
               minus sampling date), lysis buffer, instrument (13 cruises)
            b) within cruise, categorical (collection method, sequencing run, extraction
               plate): eta2 of the W score by label inside each cruise, pooled; null =
               labels rotated along the sampling order inside each cruise (keeps runs of
               consecutive samples together, so a label that only tracks position along
               the track is not mistaken for batch), 999 rotations
            c) within cruise, numeric (volume filtered, log input reads, log KO yield =
               KO counts per input read pair): Spearman, as in R33
  depth     Axis 3: refit the baseline with log10 total KO counts as a covariate
            (its own spline, as added_fields_all_kos.cv does for any variable), redo the
            residual PCA, and match the new components to the old ones.
  region    Temperature x region: the temperature spline fitted per ocean basin.
  spatial   Axis scores along each transect: correlation of score pairs as a function of
            along-track distance (correlogram), and a figure of the scores along track.
  methods   (R35) The cruise-methods table (transfer/results/protocol/cruise_methods.csv):
            cruise mean scores by extraction / preservation / pre-filter, and within-cruise
            contrasts of legs (P02 W/E, A16N 1/2) and of P02 flow-through vs cast.
  pairs     (R35) The paired protocol comparisons: A16N "Protocol A" samples and P02
            "shield" samples, each with partner samples taken at the same time and place.
            Per pair, the sample-level differences in domain CLRs, the R33 axis scores,
            and all KO CLRs; the KO-level protocol effect is compared with the axis loadings.
  nutricline (R35) Nutricline and deep-chlorophyll-maximum depth from CCHDO bottle and
            CTD profiles (transfer/results/nutricline/station_profiles.csv) against the axes.

Inputs: as scripts/residual_modules.py, plus transfer/results/protocol/sample_protocol.csv
(transfer/scripts/protocol_metadata.py). Outputs: transfer/results/residual_followup/.
The baseline residuals are cached in residuals.npz there (delete it to recompute).

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/residual_axes_followup.py --parts protocol depth
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
from added_fields_all_kos import cv  # noqa: E402
from residual_modules import corr, loco_ridge, one_driver, split_half, standardise, top_pca  # noqa: E402
from tune_figure_kos import load  # noqa: E402

OUT = REPO / "transfer" / "results" / "residual_followup"
R33 = REPO / "data" / "results" / "residual_modules"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parts", nargs="+", default=["protocol", "depth"], choices=["protocol", "depth", "region", "spatial", "methods", "pairs", "nutricline"])
    p.add_argument("--targets", type=Path, default=REPO / "data" / "targets" / "ko_clr.csv")
    p.add_argument("--embeddings-dirs", type=Path, nargs="+",
                   default=[REPO / "data" / "embeddings_sst_l4", REPO / "data" / "faire" / "embeddings_sst_l4"])
    p.add_argument("--protocol", type=Path, default=REPO / "transfer" / "results" / "protocol" / "sample_protocol.csv")
    p.add_argument("--out-dir", type=Path, default=OUT)
    p.add_argument("--rotations", type=int, default=999)
    p.add_argument("--splits", type=int, default=20)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ctx = Context(args)
    for part in args.parts:
        {"protocol": part_protocol, "depth": part_depth, "region": part_region, "spatial": part_spatial,
         "methods": part_methods, "pairs": part_pairs, "nutricline": part_nutricline}[part](ctx)
    (args.out_dir / f"summary_{'_'.join(args.parts)}.txt").write_text("\n".join(ctx.log) + "\n")


class Context:
    """Tiles, targets, baseline residuals and the R33 component scores."""

    def __init__(self, args):
        self.args, self.log = args, []
        self.rng = np.random.default_rng(args.seed)
        kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()
        self.emb, self.temp, _, self.y, self.cruises, self.unit, self.rows = load(args, kos, return_rows=True)
        cache = args.out_dir / "residuals.npz"
        if cache.exists():
            self.pred = np.load(cache)["pred"]
        else:
            _, self.pred = cv(self.emb, self.temp[:, 0], None, self.y, self.cruises, np.arange(len(self.unit)))
            np.savez_compressed(cache, pred=self.pred)
        self.res = self.y - self.pred
        self.within = self.res - pd.DataFrame(self.res).groupby(self.cruises).transform("mean").to_numpy()
        self.scores = pd.read_csv(R33 / "scores.csv", index_col="tile_key").loc[self.unit.index]
        self.loadings = pd.read_csv(R33 / "loadings.csv.gz", index_col=0)
        # sampling time and position per tile (mean over its samples)
        r = self.rows.assign(t=pd.to_datetime(self.rows["sample_time"], utc=True, errors="coerce"))
        g = r.groupby("tile_key")
        self.time = g["t"].min().loc[self.unit.index]
        self.lat = g["sample_lat"].mean().loc[self.unit.index].to_numpy()
        self.lon = g["sample_lon"].mean().loc[self.unit.index].to_numpy()
        self.say(f"{len(self.unit):,} tiles, {self.y.shape[1]:,} KOs, {len(set(self.cruises))} cruises")

    def say(self, s=""):
        print(s, flush=True)
        self.log.append(s)

    def tile_values(self, by_sample: pd.Series, how="mean") -> pd.Series:
        v = self.rows["sample_id"].map(by_sample)
        g = v.groupby(self.rows["tile_key"])
        out = g.mean() if how == "mean" else g.agg(lambda x: x.dropna().mode().iloc[0] if x.notna().any() else np.nan)
        return out.reindex(self.unit.index)


# ---------------------------------------------------------------- protocol
def part_protocol(ctx: Context) -> None:
    say = ctx.say
    pr = pd.read_csv(ctx.args.protocol, index_col="sample_key")
    vol = pd.to_numeric(pr["samp_size"], errors="coerce")
    unit = pr["samp_size_unit"].astype(str).str.lower()
    # "L" entries of 830-10,000 are mL (an 8,000 L Sterivex filtration is impossible)
    pr["volume_l"] = np.where(unit == "ml", vol / 1000, np.where((unit == "l") & (vol > 100), vol / 1000, np.where(unit == "l", vol, np.nan)))
    pr["log_input_reads"] = np.log10(pr["input_read_count"])
    pr["log_ko_yield"] = np.log10(pr["total_counts"] / pr["input_read_count"])  # KO counts per input read pair
    method = pr["collection_method"].astype(str).str.lower()
    pr["collection"] = np.select([method.isin(["underway", "inline"]), method.isin(["ctd rosette", "ctd", "bio-cast", "biocast"])],
                                 ["flow-through", "cast"], default=None)  # 'T', 'TOP' (I09N-2025) unknown  # fmt: skip
    pr["seq_date"] = pd.to_datetime(pr["seq_run_date"].astype("Int64").astype(str), format="%Y%m%d", errors="coerce")
    # 20189999 = '2018, date unknown': mid-year
    pr.loc[pr["seq_run_date"].astype("Int64").astype(str).str.endswith("9999"), "seq_date"] = pd.Timestamp("2018-07-01")

    tile = pd.DataFrame(index=ctx.unit.index)
    for c in ["volume_l", "log_input_reads", "log_ko_yield"]:
        tile[c] = ctx.tile_values(pr[c])
    for c in ["collection", "seq_flowcell", "extract_plate", "lysis_buffer", "instrument"]:
        tile[c] = ctx.tile_values(pr[c], how="mode")
    seq_date = ctx.tile_values(pr["seq_date"].astype("int64").where(pr["seq_date"].notna()) / 1e9)
    tile["storage_years"] = (seq_date - ctx.time.astype("int64") / 1e9) / (365.25 * 86400)
    tile["seq_year"] = 1970 + seq_date / (365.25 * 86400)
    tile["sample_year"] = ctx.time.dt.year + ctx.time.dt.dayofyear / 365.25
    tile["cruise"] = ctx.cruises
    tile.join(ctx.scores[["F1", "F2", "F3", "W1", "W2", "W3"]]).to_csv(ctx.args.out_dir / "tile_protocol.csv")

    # a) cruise level
    F = ["F1", "F2", "F3"]
    cm = ctx.scores[F].groupby(ctx.cruises).mean() / ctx.scores[F].std()
    cr = tile.groupby("cruise").agg(tiles=("seq_year", "size"), flowcell=("seq_flowcell", lambda v: v.mode().iloc[0]),
                                    seq_year=("seq_year", "median"), sample_year=("sample_year", "median"),
                                    storage_years=("storage_years", "median"), lysis=("lysis_buffer", lambda v: v.mode().iloc[0] if v.notna().any() else "not recorded"),
                                    instrument=("instrument", lambda v: v.mode().iloc[0]), log_ko_yield=("log_ko_yield", "median"),
                                    log_input_reads=("log_input_reads", "median"))  # fmt: skip
    cr = cr.join(cm.round(2))
    cr.round(3).to_csv(ctx.args.out_dir / "cruise_protocol.csv")
    say("\n== protocol a) cruise level: per-cruise mean score (in SDs of the tile scores)")
    say(cr.round(2).to_string())
    say("\nSpearman r across the 13 cruises (cruise mean score vs cruise median):")
    rows = []
    for v in ["sample_year", "seq_year", "storage_years", "log_ko_yield", "log_input_reads"]:
        rr = {f: corr(cr[f].rank().to_numpy(), cr[v].rank().to_numpy()) for f in F}
        rows.append({"variable": v, **rr})
        say(f"   {v:16s} " + "  ".join(f"{f} {rr[f]:+.2f}" for f in F))
    # flowcell sharing: eta2 of cruise means by flowcell, permutation over cruises
    say("Cruise means grouped by flowcell (cruises sharing a flowcell; eta2 over cruises, permutation p, 9999 draws):")
    fc = cr["flowcell"].to_numpy()
    for f in F:
        e = eta2(cr[f].to_numpy(), fc)
        null = np.array([eta2(cr[f].to_numpy(), ctx.rng.permutation(fc)) for _ in range(9999)])
        p = (1 + (null >= e).sum()) / 10000
        rows.append({"variable": "flowcell eta2 (p)", f: e, f"{f} p": p})
        say(f"   {f}: eta2 {e:.2f}  p {p:.3f}")
    for v in ["instrument", "lysis"]:
        say(f"   by {v}: " + "; ".join(f"{k}: " + ", ".join(f"{f} {g[f].mean():+.2f}" for f in F) + f" ({len(g)} cruises)"
                                     for k, g in cr.groupby(v)))
    pd.DataFrame(rows).round(3).to_csv(ctx.args.out_dir / "cruise_level_tests.csv", index=False)

    # b) within cruise, categorical
    W = ["W1", "W2", "W3"]
    order = np.argsort(ctx.time.astype("int64").to_numpy(), kind="stable")
    say(f"\n== protocol b) within cruise, categorical: pooled within-cruise eta2 of the W score; null = labels rotated "
        f"along the sampling order inside each cruise ({ctx.args.rotations} rotations)")
    rows = []
    for c in ["collection", "seq_flowcell", "extract_plate"]:
        lab = tile[c].to_numpy(dtype=object)
        usable = []
        for cruise in np.unique(ctx.cruises):
            m = (ctx.cruises == cruise) & pd.notna(lab)
            vc = pd.Series(lab[m]).value_counts()
            if (vc >= 5).sum() >= 2:
                usable.append(cruise)
        if not usable:
            say(f"   {c}: no cruise with two labels of >= 5 tiles")
            continue
        keep = np.isin(ctx.cruises, usable) & pd.notna(lab)
        say(f"   {c}: {len(usable)} cruises ({', '.join(usable)}), {keep.sum()} tiles")
        for w in W:
            e, p, per = rotation_test(ctx.scores[w].to_numpy(), lab, ctx.cruises, keep, order, ctx.rng, ctx.args.rotations)
            rows.append({"label": c, "component": w, "cruises": len(usable), "tiles": int(keep.sum()), "eta2": e, "p_rotation": p,
                         **{f"eta2 {k}": v for k, v in per.items()}})
            say(f"      {w}: eta2 {e:.3f}  p {p:.3f}   per cruise: " + ", ".join(f"{k} {v:.2f}" for k, v in per.items()))
    pd.DataFrame(rows).round(4).to_csv(ctx.args.out_dir / "within_cruise_labels.csv", index=False)

    # c) within cruise, numeric
    say("\n== protocol c) within cruise, numeric: Spearman all / within cruise (cruises with the R33 sign)")
    rows = []
    for c in ["volume_l", "log_input_reads", "log_ko_yield", "storage_years"]:
        for comp in F + W:
            d = one_driver(c, ctx.scores[comp].to_numpy(), tile[c].to_numpy(float), ctx.cruises)
            if d:
                rows.append({**d, "component": comp})
        sub = [r for r in rows if r["driver"] == c]
        if sub:
            say(f"   {c:16s} (n {sub[0]['n']}, {sub[0]['cruises']} cruises) " + "  ".join(
                f"{r['component']} {r['all_r']:+.2f}/{r['within_r']:+.2f} ({r['cruises_same_sign']}/{r['cruises_tested']})" for r in sub))
        else:
            say(f"   {c}: too few tiles")
    pd.DataFrame(rows).round(4).to_csv(ctx.args.out_dir / "within_cruise_numeric.csv", index=False)


def eta2(x, lab):
    s = pd.Series(x)
    m = s.groupby(lab).transform("mean")
    return ((m - s.mean()) ** 2).sum() / ((s - s.mean()) ** 2).sum()


def rotation_test(score, lab, cruises, keep, order, rng, n_rot):
    """Pooled within-cruise eta2 (score and labels inside cruises) and a rotation null."""
    idx = {c: [i for i in order if keep[i] and cruises[i] == c] for c in np.unique(cruises[keep])}

    def stat(labels):
        num = den = 0.0
        per = {}
        for c, ii in idx.items():
            s = pd.Series(score[ii])
            m = s.groupby(np.asarray(labels[c], dtype=object)).transform("mean")
            a, b = ((m - s.mean()) ** 2).sum(), ((s - s.mean()) ** 2).sum()
            num, den = num + a, den + b
            per[c] = a / b
        return num / den, per

    obs = {c: lab[ii] for c, ii in idx.items()}
    e, per = stat(obs)
    null = np.array([stat({c: np.roll(v, rng.integers(1, len(v))) for c, v in obs.items()})[0] for _ in range(n_rot)])
    return e, (1 + (null >= e).sum()) / (n_rot + 1), per


# ---------------------------------------------------------------- depth
def part_depth(ctx: Context) -> None:
    say = ctx.say
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    logtot = ctx.tile_values(np.log10(ks["total"])).to_numpy()
    n = len(ctx.unit)
    r2_0 = pooled_r2_each(ctx.y, ctx.pred)
    r2_1, pred1 = cv(ctx.emb, ctx.temp[:, 0], logtot, ctx.y, ctx.cruises, np.arange(n))
    np.savez_compressed(ctx.args.out_dir / "residuals_depth.npz", pred=pred1)
    gain = r2_1 - r2_0
    say("\n== depth: baseline + log10 total KO counts (spline), cruise held out")
    say(f"   median KO R2 {np.median(r2_0):.3f} -> {np.median(r2_1):.3f}; median gain {np.median(gain):+.3f}; "
        f"KOs gaining > 0.02: {(gain > 0.02).mean():.0%}; losing > 0.02: {(gain < -0.02).mean():.0%}")
    res1 = ctx.y - pred1
    obs_tot = ((ctx.y - ctx.y.mean(0)) ** 2).sum()
    say(f"   residual share of KO CLR variance {((ctx.res - ctx.res.mean(0)) ** 2).sum() / obs_tot:.1%} -> "
        f"{((res1 - res1.mean(0)) ** 2).sum() / obs_tot:.1%}")
    old = ctx.loadings
    rows = []
    for variant, mat in (("full", res1), ("within cruise", res1 - pd.DataFrame(res1).groupby(ctx.cruises).transform("mean").to_numpy())):
        z = standardise(mat)
        ev, vec = top_pca(z)
        K = 6
        rep, share = split_half(mat, ctx.cruises, vec[:, :K], ctx.rng, ctx.args.splits)
        pre = "F" if variant == "full" else "W"
        say(f"   {variant} residuals after the refit:")
        for j in range(K):
            score = z @ vec[:, j] / np.sqrt(ev[j])
            loading = vec[:, j] * np.sqrt(ev[j] / (len(z) - 1))
            if loading.sum() < 0:
                score, loading = -score, -loading
            cos = {c: float(np.dot(loading, old[c]) / np.linalg.norm(loading) / np.linalg.norm(old[c]))
                   for c in old.columns if c.startswith(pre)}
            best = max(cos, key=lambda c: abs(cos[c]))
            d = one_driver("log total", score, logtot, ctx.cruises)
            rows.append({"residuals": variant, "rank": j + 1, "share": ev[j] / ev.sum(), "split_half": rep[j],
                         "share_splits_0.85": share[j], "best_old": best, "cos_best_old": cos[best],
                         f"cos_old_{pre}3": cos[f"{pre}3"], "r_log_total_all": d["all_r"], "r_log_total_within": d["within_r"]})
            say(f"      {pre}{j + 1}': {ev[j] / ev.sum():5.1%}  split-half {rep[j]:.2f}  matches old {best} (cos {cos[best]:+.2f}; "
                f"cos with old {pre}3 {cos[f'{pre}3']:+.2f})  log total r {d['all_r']:+.2f} / within {d['within_r']:+.2f}")
    pd.DataFrame(rows).round(4).to_csv(ctx.args.out_dir / "depth_refit_components.csv", index=False)
    pd.DataFrame({"ko": old.index, "r2_baseline": r2_0, "r2_with_depth": r2_1}).round(4).to_csv(
        ctx.args.out_dir / "depth_refit_ko_r2.csv", index=False)


def pooled_r2_each(y, pred):
    return 1 - ((y - pred) ** 2).sum(0) / ((y - y.mean(0)) ** 2).sum(0)


# ---------------------------------------------------------------- region
BASIN = {"A13.5-2024": "Atlantic", "A16N-2023": "Atlantic", "AMT28": "Atlantic", "C13.5/A13.5": "Atlantic",
         "GOMECC4": "Atlantic", "HLY2202-2022": "Arctic", "I05-2023": "Indian", "I07N": "Indian", "I08S-2024": "Indian",
         "I09N": "Indian", "I09N-2025": "Indian", "P02-2022": "Pacific", "P18": "Pacific"}  # fmt: skip


def cv_interaction(emb, temp, region, y, cruises, shrink=1.0):
    """As added_fields_all_kos.cv (32 PCs + temperature spline, ridge, one cruise held out), plus the
    temperature spline multiplied by each region's indicator. A region with no training tiles gets zero
    columns, so the held-out cruise then falls back to the shared temperature curve. shrink < 1 scales the
    interaction columns down, which is the same as a ridge penalty 1/shrink^2 times larger on them."""
    from sklearn.decomposition import PCA
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import SplineTransformer, StandardScaler
    from screen_ko_predictability import ALPHAS

    pred = np.full_like(y, np.nan)
    levels = sorted(set(region))
    for cruise in np.unique(cruises):
        te, tr = cruises == cruise, cruises != cruise
        s = StandardScaler().fit(emb[tr])
        pca = PCA(n_components=32, random_state=0).fit(s.transform(emb[tr]))
        p_tr, p_te = pca.transform(s.transform(emb[tr])), pca.transform(s.transform(emb[te]))
        s2 = StandardScaler().fit(p_tr)
        sx = StandardScaler().fit(temp[tr, None])
        sp = SplineTransformer(n_knots=8, degree=3, extrapolation="linear").fit(sx.transform(temp[tr, None]))
        b_tr, b_te = sp.transform(sx.transform(temp[tr, None])), sp.transform(sx.transform(temp[te, None]))
        parts_tr, parts_te = [s2.transform(p_tr), b_tr], [s2.transform(p_te), b_te]
        for r in levels:
            parts_tr.append(shrink * b_tr * (region[tr] == r)[:, None])
            parts_te.append(shrink * b_te * (region[te] == r)[:, None])
        model = RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(np.hstack(parts_tr), y[tr])
        pred[te] = model.predict(np.hstack(parts_te))
    return pooled_r2_each(y, pred), pred


def part_region(ctx: Context) -> None:
    say = ctx.say
    lat = np.abs(ctx.lat)
    basin = np.array([BASIN[c] for c in ctx.cruises])
    band = np.select([lat < 20, lat < 35], ["tropical", "subtropical"], "temperate/polar")
    inst = pd.read_csv(REPO / "data" / "targets" / "taxa" / "instrument_target.csv").set_index("sample_key")["instrument_3310"]
    era = np.where(ctx.tile_values(inst).to_numpy() > 0.5, "OBI_3310 (2022-2025 cruises)", "OBI_2131/2631 (2011-2020 cruises)")
    variants = {f"{n}, shrink {k:g}": (lab, k) for n, lab in (("ocean basin", basin), ("latitude band (<20, 20-35, >35)", band),
                                                             ("sequencing era", era)) for k in (1.0, 0.3, 0.1)}
    r2_0 = pooled_r2_each(ctx.y, ctx.pred)
    t = ctx.temp[:, 0]
    w_old = ctx.loadings["W1"].to_numpy()
    say("\n== region: temperature spline x region, cruise held out")
    say("   regions: " + "; ".join(f"{b}: {', '.join(sorted({c for c in ctx.cruises if BASIN[c] == b}))}" for b in sorted(set(BASIN.values()))))
    d0 = one_driver("T", ctx.scores["W1"].to_numpy(), t, ctx.cruises)
    say(f"   before: median KO R2 {np.median(r2_0):.3f}; W1 vs in-situ temperature (all 1,216 tiles): "
        f"all r {d0['all_r']:+.2f}, within r {d0['within_r']:+.2f} ({d0['cruises_same_sign']}/{d0['cruises_tested']})")
    rows = [{"variant": "baseline", "median_r2": np.median(r2_0), "W1_T_within_r": d0["within_r"]}]
    for name, (region, k) in variants.items():
        r2, pred = cv_interaction(ctx.emb, t, region, ctx.y, ctx.cruises, shrink=k)
        res = ctx.y - pred
        within = res - pd.DataFrame(res).groupby(ctx.cruises).transform("mean").to_numpy()
        s = standardise(within) @ (w_old / np.linalg.norm(w_old))  # the old W1 direction in the new residuals
        d = one_driver("T", s, t, ctx.cruises)
        gain = r2 - r2_0
        say(f"   {name}: median KO R2 {np.median(r2):.3f} (median gain {np.median(gain):+.3f}; KOs gaining > 0.02: "
            f"{(gain > 0.02).mean():.0%}, losing > 0.02: {(gain < -0.02).mean():.0%}); old-W1 direction vs temperature: "
            f"within r {d['within_r']:+.2f} ({d['cruises_same_sign']}/{d['cruises_tested']})")
        per = pd.Series(pooled_r2_each(ctx.y, pred) - r2_0)
        rows.append({"variant": name, "median_r2": np.median(r2), "median_gain": np.median(gain), "share_gain_0.02": (gain > 0.02).mean(),
                     "share_loss_0.02": (gain < -0.02).mean(), "W1_T_within_r": d["within_r"], "best_gain": per.max()})
        for b in sorted(set(region)):
            m = region == b
            g = (1 - ((ctx.y[m] - pred[m]) ** 2).sum(0) / ((ctx.y[m] - ctx.y[m].mean(0)) ** 2).sum(0)) - \
                (1 - ((ctx.y[m] - ctx.pred[m]) ** 2).sum(0) / ((ctx.y[m] - ctx.y[m].mean(0)) ** 2).sum(0))
            say(f"      {b:16s} {m.sum():4d} tiles  median KO R2 change within the region {np.median(g):+.3f}")
    pd.DataFrame(rows).round(4).to_csv(ctx.args.out_dir / "region_interaction.csv", index=False)


# ---------------------------------------------------------------- spatial
BINS = [0, 15, 40, 75, 150, 300, 600, 1200, 2400, 5000]


def haversine(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def correlogram(values: np.ndarray, cruises, lat, lon):
    """Pooled within-cruise Moran-type correlation of pairs, per distance bin: sum(x_i x_j) / sqrt(sum x_i^2 sum x_j^2),
    x centred per cruise; pairs only within the same cruise."""
    out = {b: [0.0, 0.0, 0.0, 0] for b in range(len(BINS) - 1)}
    for c in np.unique(cruises):
        m = np.flatnonzero((cruises == c) & ~np.isnan(values))
        if len(m) < 10:
            continue
        x = values[m] - values[m].mean()
        i, j = np.triu_indices(len(m), 1)
        d = haversine(lat[m][i], lon[m][i], lat[m][j], lon[m][j])
        b = np.digitize(d, BINS) - 1
        for k in range(len(BINS) - 1):
            s = b == k
            if s.any():
                o = out[k]
                o[0] += (x[i[s]] * x[j[s]]).sum()
                o[1] += (x[i[s]] ** 2).sum()
                o[2] += (x[j[s]] ** 2).sum()
                o[3] += int(s.sum())
    return {k: (v[0] / np.sqrt(v[1] * v[2]) if v[3] >= 30 else np.nan, v[3]) for k, v in out.items()}


def part_spatial(ctx: Context) -> None:
    from scipy.optimize import curve_fit

    say = ctx.say
    lat, lon, cr = ctx.lat, ctx.lon, ctx.cruises
    # sample spacing along track
    order = np.argsort(ctx.time.astype("int64").to_numpy(), kind="stable")
    along = np.zeros(len(cr))
    nn = []
    for c in np.unique(cr):
        ii = [i for i in order if cr[i] == c]
        step = np.r_[0, haversine(lat[ii][:-1], lon[ii][:-1], lat[ii][1:], lon[ii][1:])]
        along[ii] = np.cumsum(step)
        m = np.array(ii)
        d = haversine(lat[m][:, None], lon[m][:, None], lat[m][None], lon[m][None])
        np.fill_diagonal(d, np.inf)
        nn.extend(d.min(1))
    say(f"\n== spatial: tiles' nearest neighbour in the same cruise: median {np.median(nn):.0f} km "
        f"(quartiles {np.percentile(nn, 25):.0f}-{np.percentile(nn, 75):.0f} km)")

    fields = pd.read_csv(REPO / "data" / "model_fields" / "sample_model_fields.csv").set_index("sample_id")
    comps = {"W1 (bacteria vs eukaryotes)": ctx.scores["W1"].to_numpy(), "W2 (phototrophs)": ctx.scores["W2"].to_numpy(),
             "W3 (sequencing depth)": ctx.scores["W3"].to_numpy()}
    refs = {"in-situ temperature": ctx.temp[:, 0],
            "log PISCES chlorophyll": np.log10(ctx.tile_values(fields["pisces_chl"].clip(lower=0.001)).to_numpy()),
            "log GLORYS mixed layer": np.log10(ctx.tile_values(fields["glorys_mlotst"].clip(lower=1)).to_numpy())}
    # single-KO residuals: median correlogram over 500 random KOs
    ko_pick = ctx.rng.choice(ctx.within.shape[1], 500, replace=False)
    rows = []
    for name, v in {**comps, **refs}.items():
        cg = correlogram(v, cr, lat, lon)
        for k, (r, n) in cg.items():
            rows.append({"series": name, "bin": f"{BINS[k]}-{BINS[k + 1]}", "mid_km": np.sqrt(max(BINS[k], 5) * BINS[k + 1]), "r": r, "pairs": n})
    ko_r = np.array([[correlogram(ctx.within[:, j], cr, lat, lon)[k][0] for k in range(len(BINS) - 1)] for j in ko_pick])
    for k in range(len(BINS) - 1):
        rows.append({"series": "single KO residuals (median of 500)", "bin": f"{BINS[k]}-{BINS[k + 1]}",
                     "mid_km": np.sqrt(max(BINS[k], 5) * BINS[k + 1]), "r": np.nanmedian(ko_r[:, k]), "pairs": np.nan})
    cg = pd.DataFrame(rows)
    cg.round(4).to_csv(ctx.args.out_dir / "correlogram.csv", index=False)
    wide = cg.pivot(index="series", columns="bin", values="r")[[f"{BINS[k]}-{BINS[k + 1]}" for k in range(len(BINS) - 1)]]
    say("Within-cruise correlation of pairs of tiles by distance (km); pooled over cruises, >= 30 pairs per bin:")
    say(wide.round(2).to_string())
    say("pairs per bin: " + ", ".join(f"{b}: {int(n)}" for b, n in cg[cg["series"] == "W2 (phototrophs)"][["bin", "pairs"]].values))

    # exponential fit r(d) = c exp(-d / L): c = share of variance that is spatially structured at the smallest scale
    fit_rows = []
    say("Exponential fit r(d) = c * exp(-d / L) over the bins (weights = pairs): c = correlation at d -> 0 between "
        "neighbouring tiles (1 - c = share that is sample-to-sample), L = decorrelation length:")
    for name in list(comps) + list(refs):
        g = cg[(cg["series"] == name) & cg["r"].notna()]
        try:
            (c, L), _ = curve_fit(lambda d, c, L: c * np.exp(-d / L), g["mid_km"], g["r"], p0=[0.5, 300],
                                  sigma=1 / np.sqrt(g["pairs"]), bounds=([-1, 1], [1.5, 20000]))
        except RuntimeError:
            c, L = np.nan, np.nan
        fit_rows.append({"series": name, "c": c, "L_km": L})
        say(f"   {name:34s} c {c:.2f}  L {L:5.0f} km")
    pd.DataFrame(fit_rows).round(3).to_csv(ctx.args.out_dir / "correlogram_fit.csv", index=False)

    # replicates / co-located samples: samples that share one tile
    replicate_check(ctx)
    along_track_figure(ctx, along)
    correlogram_figure(ctx, cg)


def replicate_check(ctx: Context) -> None:
    """Samples merged into one tile (same image; GOMECC4 A/B/C replicates and co-located samples): how well do two
    samples of the same tile agree on each W score? Sample residual = sample CLR - the tile's prediction - the
    cruise mean residual; standardised with the tile matrix's KO SDs; projected on the component's loadings."""
    say = ctx.say
    multi = ctx.rows.groupby("tile_key")["sample_id"].agg(list)
    multi = multi[multi.str.len() >= 2].reindex(ctx.unit.index).dropna()
    if multi.empty:
        say("Replicates: no tile with two samples")
        return
    need = sorted({s for v in multi for s in v})
    clr = pd.read_csv(ctx.args.targets, index_col="sample_key")
    clr = clr.loc[need, [c for c in clr.columns if c != "sample_key"]]
    pos = {k: i for i, k in enumerate(ctx.unit.index)}
    cm = pd.DataFrame(ctx.res).groupby(ctx.cruises).transform("mean").to_numpy()
    mu, sd = ctx.within.mean(0), ctx.within.std(0)
    sd = np.where(sd > 0, sd, 1)
    rows = []
    for comp in ["W1", "W2", "W3"]:
        w = ctx.loadings[comp].to_numpy()
        w = w / np.linalg.norm(w)
        a, b, dist, cru = [], [], [], []
        for tk, samples in multi.items():
            i = pos[tk]
            z = (clr.loc[samples].to_numpy() - ctx.pred[i] - cm[i] - mu) / sd
            s = z @ w
            a.append(s[0])
            b.append(s[1])
            cru.append(ctx.cruises[i])
            ss = ctx.rows.set_index("sample_id").loc[samples[:2]]
            dist.append(haversine(ss["sample_lat"].iloc[0], ss["sample_lon"].iloc[0], ss["sample_lat"].iloc[1], ss["sample_lon"].iloc[1]))
        a, b = np.array(a), np.array(b)
        tile_s = standardise(ctx.within) @ w
        r = corr(a, b)
        # share of the tile-score variance that is between-sample noise: var(a - b) / 2 over var of tile scores
        noise = np.var(a - b) / 2 / np.var(tile_s)
        rows.append({"component": comp, "tiles": len(a), "r_pair": r, "noise_share": noise,
                     "median_pair_distance_km": float(np.median(dist)), "cruises": "; ".join(f"{k} {v}" for k, v in pd.Series(cru).value_counts().items())})
    out = pd.DataFrame(rows)
    out.round(4).to_csv(ctx.args.out_dir / "replicates.csv", index=False)
    say(f"\nTwo samples in one tile ({len(multi)} tiles; {out['cruises'].iloc[0]}; median pair distance "
        f"{out['median_pair_distance_km'].iloc[0]:.1f} km): agreement of the sample scores")
    for r in out.itertuples():
        say(f"   {r.component}: r between the two samples {r.r_pair:+.2f}; sample-to-sample variance = {r.noise_share:.0%} "
            f"of the tile-score variance")


def _arial():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    return plt


def along_track_figure(ctx: Context, along: np.ndarray) -> None:
    plt = _arial()
    cruises = sorted(np.unique(ctx.cruises), key=lambda c: -np.sum(ctx.cruises == c))
    fig, axes = plt.subplots(5, 3, figsize=(7.2, 8.4), sharey=True)
    s = ctx.scores["W2"].to_numpy()
    s = s / s.std()
    for ax, c in zip(axes.flat, cruises):
        m = ctx.cruises == c
        o = np.argsort(along[m])
        x, y = along[m][o], s[m][o]
        ax.plot(x, y, color="#0072B2", lw=0.6, alpha=0.5)
        ax.scatter(x, y, s=6, color="#0072B2", lw=0)
        k = max(3, int(round(len(y) / 15)) | 1)
        ax.plot(x, pd.Series(y).rolling(k, center=True, min_periods=2).median(), color="#D55E00", lw=1.4)
        ax.axhline(0, color="#999999", lw=0.5)
        ax.text(0.02, 0.95, f"{c}  (n {m.sum()})", transform=ax.transAxes, va="top", fontsize=7)
        ax.set_xlim(0, max(x.max(), 1))
    for ax in axes.flat[len(cruises):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("W2 score (SD)")
    for ax in axes[-1, :]:
        ax.set_xlabel("Distance along track (km)")
    axes.flat[len(cruises) - 1].set_xlabel("Distance along track (km)")
    fig.tight_layout()
    out = REPO / "transfer" / "figures" / "w2_along_track.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    ctx.say(f"Figure: {out.relative_to(REPO)}")


def correlogram_figure(ctx: Context, cg: pd.DataFrame) -> None:
    plt = _arial()
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    style = {"W1 (bacteria vs eukaryotes)": ("#009E73", "-", "o"), "W2 (phototrophs)": ("#0072B2", "-", "o"),
             "W3 (sequencing depth)": ("#CC79A7", "-", "o"), "in-situ temperature": ("#555555", "--", "s"),
             "log PISCES chlorophyll": ("#E69F00", "--", "s"), "single KO residuals (median of 500)": ("#999999", ":", "^")}
    for name, (col, ls, mk) in style.items():
        g = cg[(cg["series"] == name) & cg["r"].notna()]
        ax.plot(g["mid_km"], g["r"], ls=ls, marker=mk, ms=3.5, lw=1.2, color=col, label=name)
    ax.set_xscale("log")
    ax.axhline(0, color="#999999", lw=0.5)
    ax.set_xlabel("Distance between tiles, same cruise (km)")
    ax.set_ylabel("Correlation of pairs")
    ax.legend(frameon=False, fontsize=6.5, loc="lower left")
    fig.tight_layout()
    out = REPO / "transfer" / "figures" / "residual_axes_correlogram.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    ctx.say(f"Figure: {out.relative_to(REPO)}")


# ---------------------------------------------------------------- methods (R35)
LEG_SPLIT = {"P02-2022": ("2022-06-13", "P02W (leg 1)", "P02E (leg 2)"), "A16N-2023": ("2023-04-10", "A16N leg 1", "A16N leg 2")}


def part_methods(ctx: Context) -> None:
    say = ctx.say
    meth = pd.read_csv(REPO / "transfer" / "results" / "protocol" / "cruise_methods.csv")
    # one coding per cruise for the cruise-level comparison (legs: the majority / standard protocol)
    per = meth.drop_duplicates("cruise", keep="first").set_index("cruise")
    per.loc["A16N-2023", "preservation"] = "DNA/RNA Shield"  # predominantly Shield; the Protocol A samples are paired (part pairs)
    per.loc["P02-2022", "preservation"] = "lysis buffer"  # lysis buffer normally; 11 Shield samples are paired
    F = ["F1", "F2", "F3"]
    cm = ctx.scores[F].groupby(ctx.cruises).mean() / ctx.scores[F].std()
    inst = pd.read_csv(ctx.args.out_dir / "cruise_protocol.csv", index_col=0)[["instrument", "sample_year"]]
    cr = cm.join(per[["preservation", "extraction", "prefilter", "volume_l", "store_c"]]).join(inst)
    cr.round(3).to_csv(ctx.args.out_dir / "cruise_methods_scores.csv")
    say("\n== methods: cruise mean scores (SD) with the coded methods table")
    say(cr.round(2).sort_values(["extraction", "sample_year"]).to_string())
    for col in ["extraction", "preservation", "prefilter", "instrument"]:
        say(f"   by {col}: " + "; ".join(f"{k} ({len(g)}): " + ", ".join(f"{f} {g[f].mean():+.2f}" for f in F) for k, g in cr.groupby(col)))
    say("   Spearman r with sampling year inside each extraction group (does the era trend survive within a protocol?):")
    for k, g in cr.groupby("extraction"):
        if len(g) >= 4:
            say(f"      {k} ({len(g)} cruises): " + ", ".join(f"{f} {corr(g[f].rank().to_numpy(), g['sample_year'].rank().to_numpy()):+.2f}" for f in F))
    say("   Spearman r with typical volume filtered (13 cruises): " + ", ".join(
        f"{f} {corr(cr[f].rank().to_numpy(), cr['volume_l'].rank().to_numpy()):+.2f}" for f in F))

    # within-cruise contrasts: legs, and P02 flow-through (30 um pre-filtered on leg 1) vs cast
    tp = pd.read_csv(ctx.args.out_dir / "tile_protocol.csv", index_col=0).loc[ctx.unit.index]
    order = np.argsort(ctx.time.astype("int64").to_numpy(), kind="stable")
    rows = []
    say("   Within-cruise contrasts (W scores, SD of all tile scores; rotation p as in R34):")
    for cruise, (cut, a, b) in LEG_SPLIT.items():
        m = ctx.cruises == cruise
        lab = np.where(ctx.time.dt.tz_convert(None) < pd.Timestamp(cut), a, b).astype(object)
        lab[~m] = None
        for w in ["W1", "W2", "W3"]:
            s = ctx.scores[w].to_numpy() / ctx.scores[w].std()
            e, p, _ = rotation_test(s, lab, ctx.cruises, m, order, ctx.rng, ctx.args.rotations)
            d = s[m & (lab == b)].mean() - s[m & (lab == a)].mean()
            rows.append({"contrast": f"{b} - {a}", "component": w, "n_a": int((m & (lab == a)).sum()), "n_b": int((m & (lab == b)).sum()),
                         "difference_sd": d, "eta2": e, "p_rotation": p})
            say(f"      {b} - {a} ({rows[-1]['n_b']} vs {rows[-1]['n_a']} tiles), {w}: {d:+.2f} SD, eta2 {e:.3f}, p {p:.3f}")
    m = (ctx.cruises == "P02-2022") & (ctx.time.dt.tz_convert(None) < pd.Timestamp("2022-06-13")).to_numpy()
    lab = tp["collection"].to_numpy(dtype=object).copy()
    lab[~m] = None
    keep = m & pd.notna(lab)
    for w in ["W1", "W2", "W3"]:
        s = ctx.scores[w].to_numpy() / ctx.scores[w].std()
        e, p, _ = rotation_test(s, lab, ctx.cruises, keep, order, ctx.rng, ctx.args.rotations)
        d = s[keep & (lab == "flow-through")].mean() - s[keep & (lab == "cast")].mean()
        rows.append({"contrast": "P02W flow-through (30 um pre-filter) - cast", "component": w, "n_a": int((keep & (lab == "cast")).sum()),
                     "n_b": int((keep & (lab == "flow-through")).sum()), "difference_sd": d, "eta2": e, "p_rotation": p})
        say(f"      P02W flow-through (30 um pre-filtered) - cast ({rows[-1]['n_b']} vs {rows[-1]['n_a']} tiles), {w}: {d:+.2f} SD, "
            f"eta2 {e:.3f}, p {p:.3f}")
    pd.DataFrame(rows).round(4).to_csv(ctx.args.out_dir / "methods_within_cruise.csv", index=False)


# ---------------------------------------------------------------- pairs (R35)
def find_pairs(workbook: Path) -> pd.DataFrame:
    """Treated sample + its partners (same cruise, same eventDate, <= 10 m, other MG samples)."""
    import warnings

    warnings.filterwarnings("ignore")
    sm = pd.read_excel(workbook, sheet_name="sampleMetadata", header=2)
    sm = sm[sm["samp_name"].astype(str).str.endswith("_MG")]
    sm["depth"] = pd.to_numeric(sm["minimumDepthInMeters"], errors="coerce")
    sm = sm[sm["depth"] <= 10]
    notes = sm["samp_collect_notes"].astype(str)
    treated = {"A16N-2023": ("Protocol A", notes.str.contains("Protocol A")),
               "P02-2022": ("DNA/RNA Shield", notes.str.strip().str.lower().eq("shield"))}
    rows = []
    for cruise, (label, flag) in treated.items():
        for _, t in sm[flag & (sm["short_name"] == cruise)].iterrows():
            partners = sm[(sm["short_name"] == cruise) & (sm["eventDate"] == t["eventDate"]) & (sm["samp_name"] != t["samp_name"])]
            for _, q in partners.iterrows():
                rows.append({"cruise": cruise, "treatment": label, "treated": t["samp_name"], "partner": q["samp_name"],
                             "event": t["eventDate"], "partners": len(partners), "treated_note": t["samp_collect_notes"],
                             "partner_note": q["samp_collect_notes"]})
    return pd.DataFrame(rows)


def part_pairs(ctx: Context) -> None:
    from scipy.stats import wilcoxon

    say = ctx.say
    pairs = find_pairs(REPO / "FAIRe-NOAA_uci-ess-biogoship.xlsx")
    clr = pd.read_csv(ctx.args.targets, index_col="sample_key")
    kos = list(clr.columns)
    ks = pd.read_csv(REPO / "data" / "targets" / "ko_samples.csv").set_index("sample_key")
    unc = pd.read_csv(REPO / "data" / "targets" / "taxa" / "unclassified_targets.csv").set_index("sample_key")
    dom = pd.read_csv(REPO / "data" / "targets" / "taxa" / "taxa_domain_clr.csv").set_index("sample_key")
    tile_of = ctx.rows.set_index("sample_id")["tile_key"]
    pos = {k: i for i, k in enumerate(ctx.unit.index)}
    have = pairs["treated"].isin(clr.index) & pairs["partner"].isin(clr.index)
    say(f"\n== pairs: {len(pairs)} treated-partner pairs found ({pairs.groupby('cruise').size().to_dict()}); "
        f"{int(have.sum())} with both samples among the modelled KO samples")
    pairs = pairs[have].copy()
    pairs["same_tile"] = [tile_of.get(a) == tile_of.get(b) and tile_of.get(a) is not None for a, b in zip(pairs["treated"], pairs["partner"])]
    say(f"   both samples in the same modelled tile: {int(pairs['same_tile'].sum())} of {len(pairs)} (the R33 tiles averaged them)")

    # sample-level axis scores (as replicate_check): residual vs the tile prediction, minus the cruise mean residual
    cm = pd.DataFrame(ctx.res).groupby(ctx.cruises).transform("mean").to_numpy()
    mu_w, sd_w = ctx.within.mean(0), np.where(ctx.within.std(0) > 0, ctx.within.std(0), 1)
    mu_f, sd_f = ctx.res.mean(0), np.where(ctx.res.std(0) > 0, ctx.res.std(0), 1)
    comps = ["F1", "F2", "F3", "W1", "W2", "W3"]
    unit_w = {c: ctx.loadings[c].to_numpy() / np.linalg.norm(ctx.loadings[c].to_numpy()) for c in comps}
    tile_sd = {c: (standardise(ctx.within if c[0] == "W" else ctx.res) @ unit_w[c]).std() for c in comps}

    def sample_scores(sample):
        t = tile_of.get(sample)
        if t is None or t not in pos:
            return {c: np.nan for c in comps}
        i = pos[t]
        y = clr.loc[sample].to_numpy()
        out = {}
        for c in comps:
            z = (y - ctx.pred[i] - cm[i] - mu_w) / sd_w if c[0] == "W" else (y - ctx.pred[i] - mu_f) / sd_f
            out[c] = z @ unit_w[c] / tile_sd[c]
        return out

    feats = {}
    for s in set(pairs["treated"]) | set(pairs["partner"]):
        feats[s] = {"Bacteria CLR": dom.loc[s, "Bacteria"] if s in dom.index else np.nan,
                    "Eukaryota CLR": dom.loc[s, "Eukaryota"] if s in dom.index else np.nan,
                    "Archaea CLR": dom.loc[s, "Archaea"] if s in dom.index else np.nan,
                    "bacteria vs eukaryota": unc.loc[s, "bacteria_vs_eukaryota"] if s in unc.index else np.nan,
                    "unclassified logit": unc.loc[s, "unclassified_logit"] if s in unc.index else np.nan,
                    "log10 total KO counts": np.log10(ks.loc[s, "total"]) if s in ks.index else np.nan,
                    **sample_scores(s)}
    feats = pd.DataFrame(feats).T
    for c in feats.columns:
        pairs[f"d {c}"] = feats.loc[pairs["treated"], c].to_numpy() - feats.loc[pairs["partner"], c].to_numpy()

    # replicate baseline: partner vs partner at the same event (P02 triples)
    reps = []
    for (cruise, ev), g in pairs.groupby(["cruise", "event"]):
        p = list(g["partner"])
        if len(p) >= 2:
            reps.append({"cruise": cruise, "a": p[0], "b": p[1]})
    reps = pd.DataFrame(reps)

    rows = []
    for cruise, g in pairs.groupby("cruise"):
        # one difference per treated sample: treated minus the mean of its partners
        per = g.groupby("treated")[[c for c in g if c.startswith("d ")]].mean()
        rb = reps[reps["cruise"] == cruise] if len(reps) else reps
        say(f"   {cruise}: {g['treatment'].iloc[0]} minus partner, {len(per)} treated samples "
            f"({len(g)} pairs; partner-partner baseline from {len(rb)} events)")
        for c in feats.columns:
            d = per[f"d {c}"].dropna()
            if len(d) < 4:
                continue
            p = wilcoxon(d).pvalue if (d != 0).any() else np.nan
            base = (feats.loc[rb["a"], c].to_numpy() - feats.loc[rb["b"], c].to_numpy()) if len(rb) else np.array([])
            base_sd = np.nanstd(base, ddof=1) / np.sqrt(2) * np.sqrt(2) if len(base) >= 3 else np.nan
            rows.append({"cruise": cruise, "treatment": g["treatment"].iloc[0], "variable": c, "n": len(d), "mean_diff": d.mean(),
                         "sd_diff": d.std(ddof=1), "same_sign": f"{int((np.sign(d) == np.sign(d.mean())).sum())}/{len(d)}",
                         "wilcoxon_p": p, "replicate_sd_diff": base_sd, "replicate_n": len(base)})
            say(f"      {c:24s} mean {d.mean():+.3f} (SD {d.std(ddof=1):.3f}; {rows[-1]['same_sign']} same sign; Wilcoxon p {p:.3f})"
                + (f"; partner-partner SD of difference {base_sd:.3f} (n {len(base)})" if np.isfinite(base_sd) else ""))
    pd.DataFrame(rows).round(4).to_csv(ctx.args.out_dir / "pairs_effects.csv", index=False)
    pairs.round(4).to_csv(ctx.args.out_dir / "pairs.csv", index=False)

    # KO-level protocol effect vs the axis loadings
    sd_ko = ctx.within.std(0)
    say("   KO-level effect (mean treated - partner difference per KO / the KO's within-cruise residual SD): "
        "cosine with the R33 loadings; and its own strength")
    eff_rows = []
    for cruise, g in pairs.groupby("cruise"):
        diff = clr.loc[g["treated"]].to_numpy() - clr.loc[g["partner"]].to_numpy()
        per = pd.DataFrame(diff).groupby(g["treated"].to_numpy()).mean().to_numpy()
        eff = per.mean(0) / sd_ko
        t = per.mean(0) / (per.std(0, ddof=1) / np.sqrt(len(per)))
        cos = {c: float(eff @ ctx.loadings[c].to_numpy() / np.linalg.norm(eff) / np.linalg.norm(ctx.loadings[c].to_numpy()))
               for c in ["F1", "F2", "F3", "F4", "F5", "W1", "W2", "W3", "W4"]}
        # null for the cosine: random sign flips of whole pairs
        null = []
        for _ in range(999):
            sg = ctx.rng.choice([-1, 1], len(per))[:, None]
            e0 = (per * sg).mean(0) / sd_ko
            null.append(abs(e0 @ ctx.loadings["F1"].to_numpy()) / np.linalg.norm(e0) / np.linalg.norm(ctx.loadings["F1"].to_numpy()))
        p_f1 = (1 + (np.array(null) >= abs(cos["F1"])).sum()) / 1000
        frac = (np.abs(t) > 3).mean()
        eff_rows.append({"cruise": cruise, "treated": len(per), "rms_effect_sd": float(np.sqrt((eff ** 2).mean())),
                         "share_kos_abs_t_gt_3": frac, **{f"cos {c}": v for c, v in cos.items()}, "p_cos_F1_signflip": p_f1})
        top = pd.Series(eff, index=[k.replace("count.ko.", "") for k in kos]).sort_values()
        say(f"      {cruise}: RMS effect {eff_rows[-1]['rms_effect_sd']:.2f} SD per KO; KOs with |t| > 3: {frac:.1%}; cosine "
            + ", ".join(f"{c} {v:+.2f}" for c, v in cos.items()) + f"; p(|cos F1|, sign flips) {p_f1:.3f}")
        say(f"         most up in the treated samples: {', '.join(top.index[-6:][::-1])}; most down: {', '.join(top.index[:6])}")
    pd.DataFrame(eff_rows).round(4).to_csv(ctx.args.out_dir / "pairs_ko_effect.csv", index=False)


PROFILE = {"nitracline depth (1 umol/kg, m)": "nitracline_1um_m", "nitracline, steepest gradient (m)": "nitracline_grad_m",
           "DCM depth (m)": "dcm_m", "DCM / surface fluorescence (log)": "log_dcm_ratio", "CTD mixed layer (m)": "mld_ctd_m",
           "nitracline - mixed layer (m)": "nitracline_minus_mld", "surface nitrate (log)": "log_no3_surface"}  # fmt: skip


def match_casts(ctx: Context, casts: pd.DataFrame, max_km: float, max_days: float) -> pd.DataFrame:
    """Per tile: the nearest cast of its own cruise within max_km and max_days."""
    casts = casts.assign(t=pd.to_datetime(casts["time"], utc=True, errors="coerce"))
    out = []
    for i, (c, la, lo, t) in enumerate(zip(ctx.cruises, ctx.lat, ctx.lon, ctx.time)):
        g = casts[casts["section"] == c]
        if g.empty or pd.isna(t):
            out.append({})
            continue
        d = haversine(la, lo, g["lat"].to_numpy(), g["lon"].to_numpy())
        dt = np.abs((g["t"] - t).dt.total_seconds().to_numpy()) / 86400
        ok = (d <= max_km) & (dt <= max_days)
        if not ok.any():
            out.append({})
            continue
        k = np.flatnonzero(ok)[np.argmin(d[ok])]
        out.append({**g.iloc[k][list(PROFILE.values())].to_dict(), "cast_km": d[k], "cast_days": dt[k]})
    return pd.DataFrame(out, index=ctx.unit.index)


def part_nutricline(ctx: Context) -> None:
    say = ctx.say
    casts = pd.read_csv(REPO / "transfer" / "results" / "nutricline" / "station_profiles.csv")
    casts["log_dcm_ratio"] = np.log10(casts["dcm_ratio"].clip(lower=0.1))
    casts["nitracline_minus_mld"] = casts["nitracline_1um_m"] - casts["mld_ctd_m"]
    casts["log_no3_surface"] = np.log10(casts["no3_surface"].clip(lower=0.01))
    fields = pd.read_csv(REPO / "data" / "model_fields" / "sample_model_fields.csv").set_index("sample_id")
    chl = np.log10(ctx.tile_values(fields["pisces_chl"].clip(lower=0.001)).to_numpy())
    comps = ["F1", "F2", "W1", "W2", "W3"]
    rows = []
    for max_km, max_days in ((100, 3), (30, 1)):
        m = match_casts(ctx, casts, max_km, max_days)
        n = m["cast_km"].notna()
        say(f"\n== nutricline: tiles matched to a cast of their cruise within {max_km} km and {max_days} d: {int(n.sum())} of {len(n)} "
            f"(median {m['cast_km'].median():.0f} km, {m['cast_days'].median():.2f} d); per cruise: "
            + ", ".join(f"{c} {k}" for c, k in pd.Series(ctx.cruises[n.to_numpy()]).value_counts().items()))
        m.round(3).assign(cruise=ctx.cruises).to_csv(ctx.args.out_dir / f"tile_profiles_{max_km}km.csv")
        for label, col in PROFILE.items():
            x = m[col].to_numpy(float)
            parts = []
            for comp in comps:
                d = one_driver(label, ctx.scores[comp].to_numpy(), x, ctx.cruises)
                if d:
                    rows.append({**d, "component": comp, "match": f"{max_km} km / {max_days} d"})
                    parts.append(f"{comp} {d['all_r']:+.2f}/{d['within_r']:+.2f} ({d['cruises_same_sign']}/{d['cruises_tested']})")
            if parts:
                say(f"   {label:34s} n {int(np.isfinite(x).sum()):4d}  " + "  ".join(parts))
        # does the profile add to W2 beyond surface chlorophyll? within-cruise ridge, cruise held out, same tiles
        keep = np.isfinite(m["nitracline_1um_m"].to_numpy(float)) & np.isfinite(chl)
        if keep.sum() > 100:
            s = ctx.scores["W2"].to_numpy()[keep]
            cr = ctx.cruises[keep]

            def r2(cols):
                x = np.column_stack(cols)[keep]
                x = x - pd.DataFrame(x).groupby(cr).transform("mean").to_numpy()
                return loco_ridge(x, s, cr)

            nit = m["nitracline_1um_m"].to_numpy(float)
            mld = m["mld_ctd_m"].to_numpy(float)
            say(f"   W2, within-cruise ridge R2 (cruise held out, {int(keep.sum())} tiles): chlorophyll {r2([chl]):+.3f}; "
                f"nitracline {r2([nit]):+.3f}; both {r2([chl, nit]):+.3f}; nitracline + CTD mixed layer + chlorophyll "
                f"{r2([chl, nit, np.where(np.isnan(mld), np.nanmedian(mld), mld)]):+.3f}")
    pd.DataFrame(rows).round(4).to_csv(ctx.args.out_dir / "nutricline_drivers.csv", index=False)



if __name__ == "__main__":
    main()
