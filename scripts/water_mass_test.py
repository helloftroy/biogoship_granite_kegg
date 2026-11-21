#!/usr/bin/env python3
"""Water-mass history along 60-day back-trajectories as predictors (RESULTS.md R53; data R42 / R52).

Input: transfer/results/water_mass/along_track_daily.csv (sample_id, day 0 ... -60, lat, lon, variable, value;
variables sst, chl, kd490, mld, depth). A track whose particle got stuck at the coast has every variable
missing from some day on (R52); its history is used only up to the day before the first missing day and the
sample is flagged (features needing missing days become missing, then the training median).

Features per sample (then averaged over the samples of a tile, as the targets):
  shelf_days        days since the water was last over the shelf (sea floor < 200 m); 0 = on the shelf now,
                    61 = not within 60 days
  sst_change30      SST day 0 - mean SST days -30..-1 (C)
  sst_sd60          SD of SST over the track (C)
  chl_mean60        mean log10 chlorophyll over days -60..0
  chl_change30      log10 chl day 0 - mean log10 chl days -30..-1
  chl_peak60        max log10 chl days -60..-1 - log10 chl day 0 (a past bloom the water has left behind)
  mld_mean30        mean log10 MLD days -30..0
  mld_change30      log10 MLD day 0 - mean log10 MLD days -30..-1
  light_ml0         log10 mean light fraction in the mixed layer at day 0, (1 - exp(-K MLD)) / (K MLD),
                    K = KD490
  light_ml30        the same, mean over days -30..0
  displacement_km   great-circle distance day 0 -> oldest valid day
  dlat_origin       latitude day 0 - latitude of the oldest valid day (water from the pole side > 0 in N)
  control_shuffled  displacement_km permuted within cruise (null for the gain of one added spline)

Tests (R43 tiles, one cruise held out, ridge, per-KO GCV alpha; scripts/domain_split_models.cv_predict):
  1  added value: baseline (32 Granite PCs + temperature spline) + one feature (spline) or + all features
     (splines); explained share (raw, KOs standardised) on the current CLR (3,933 KOs), C (eukaryote /
     prokaryote amount) and A (prokaryote-specific KOs per genome equivalent) (R51); gain over baseline.
     Also Granite PCs alone + all features (Granite alone is the better R51 model).
  2  residual axes: Spearman, within cruise (features and scores centred per cruise), of each feature with
     the R36 within-cruise residual axis scores W1 (bacteria vs eukaryotes) and W2 (phototrophs), and with
     the out-of-fold baseline residual of C.

Outputs: data/results/water_mass/{features_sample,features_tile,added_value,axis_correlations}.csv,
summary.txt; figures/output/water_mass_test.png.

Usage:
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python3 scripts/water_mass_test.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from domain_split_models import EMB, FSETS, cv_predict, scores  # noqa: E402
from tune_figure_kos import load  # noqa: E402

OUT = REPO / "data" / "results" / "water_mass"
FIG = REPO / "figures" / "output"
TRACK = REPO / "transfer" / "results" / "water_mass" / "along_track_daily.csv"
FEATS = ["shelf_days", "sst_change30", "sst_sd60", "chl_mean60", "chl_change30", "chl_peak60", "mld_mean30", "mld_change30",
         "light_ml0", "light_ml30", "displacement_km", "dlat_origin"]  # fmt: skip


def features() -> pd.DataFrame:
    d = pd.read_csv(TRACK)
    w = d.pivot_table(index=["sample_id", "day"], columns="variable", values="value", dropna=False)
    pos = d.drop_duplicates(["sample_id", "day"]).set_index(["sample_id", "day"])[["lat", "lon"]]
    w = w.join(pos).reset_index()
    allnan = w[["sst", "chl", "kd490", "mld", "depth"]].isna().all(1)
    first_bad = w[allnan].groupby("sample_id")["day"].max()  # most recent all-missing day
    w = w[~(w["day"] <= w["sample_id"].map(first_bad).fillna(-999))]
    V = {v: w.pivot(index="sample_id", columns="day", values=v) for v in ["sst", "chl", "kd490", "mld", "depth", "lat", "lon"]}
    days = V["sst"].columns  # -60 .. 0
    prior30 = [c for c in days if -30 <= c <= -1]
    prior60 = [c for c in days if c <= -1]
    last30 = [c for c in days if c >= -30]
    lchl, lmld = np.log10(V["chl"]), np.log10(V["mld"])
    k, m = V["kd490"], V["mld"]
    light = np.log10((1 - np.exp(-k * m)) / (k * m))
    shelf = V["depth"] < 200
    on = shelf.where(V["depth"].notna())
    sd = pd.Series(61.0, index=on.index)
    for i in on.index:
        r = on.loc[i].dropna()
        hit = r[r.astype(bool)]
        if len(hit):
            sd[i] = -hit.index.max()
    lat, lon = V["lat"], V["lon"]
    oldest = lat.notna().iloc[:, ::1].idxmax(1)  # smallest valid day (columns ascending)
    lat0, lon0 = lat[0], lon[0]
    lato = pd.Series([lat.at[s, oldest[s]] for s in lat.index], index=lat.index)
    lono = pd.Series([lon.at[s, oldest[s]] for s in lon.index], index=lon.index)
    p1, p2, dl = np.radians(lat0), np.radians(lato), np.radians(lon0 - lono)
    disp = 6371 * 2 * np.arcsin(np.sqrt(np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2))
    f = pd.DataFrame({
        "shelf_days": sd,
        "sst_change30": V["sst"][0] - V["sst"][prior30].mean(1),
        "sst_sd60": V["sst"].std(1),
        "chl_mean60": lchl.mean(1),
        "chl_change30": lchl[0] - lchl[prior30].mean(1),
        "chl_peak60": lchl[prior60].max(1) - lchl[0],
        "mld_mean30": lmld[last30].mean(1),
        "mld_change30": lmld[0] - lmld[prior30].mean(1),
        "light_ml0": light[0],
        "light_ml30": light[last30].mean(1),
        "displacement_km": disp,
        "dlat_origin": lat0 - lato,
    })  # fmt: skip
    f["valid_days"] = V["sst"].notna().sum(1)
    f["stuck"] = f.index.isin(first_bad.index)
    return f


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    f = features()
    f.round(5).rename_axis("sample_id").to_csv(OUT / "features_sample.csv")
    say(f"samples with tracks: {len(f):,}; stuck (history truncated): {int(f['stuck'].sum())}; valid days median "
        f"{f['valid_days'].median():.0f} (min {f['valid_days'].min()})")
    say((f[FEATS].describe(percentiles=[0.1, 0.5, 0.9]).T[["count", "10%", "50%", "90%"]]).round(3).to_string())

    TD = REPO / "data" / "targets" / "domain_split"
    targets = {"clr_all": REPO / "data" / "targets" / "ko_clr.csv", "C": TD / "C.csv", "A": TD / "A.csv"}
    rng = np.random.default_rng(0)
    rows_av, rows_cor = [], []
    resid_C = None
    for name, path in targets.items():
        kos = pd.read_csv(path, nrows=0).columns.drop("sample_key").tolist()
        emb, temp, sat, y, cr, unit, rows = load(argparse.Namespace(targets=path, embeddings_dirs=EMB), kos, return_rows=True)
        ft = f[FEATS].reindex(rows["sample_id"]).set_axis(rows.index).groupby(rows["tile_key"]).mean().reindex(unit.index)
        cover = ft.notna().all(1).mean()
        ctl = ft["displacement_km"].copy()
        for c in np.unique(cr):
            ix = np.where(cr == c)[0]
            ctl.iloc[ix] = rng.permutation(ctl.iloc[ix].to_numpy())
        ft["control_shuffled"] = ctl
        ftf = ft.fillna(ft.median())
        if name == "clr_all":
            ft.round(5).to_csv(OUT / "features_tile.csv")
        say(f"\n{name}: {y.shape[0]:,} tiles x {y.shape[1]:,}; tiles with every feature {cover:.1%}")
        base = cv_predict(FSETS["baseline"], emb, temp, sat, y, cr)
        sb, _ = scores(y, base)
        gran = cv_predict(FSETS["granite"], emb, temp, sat, y, cr)
        sg, _ = scores(y, gran)
        if name == "C":
            resid_C = pd.Series((y - base)[:, 0], index=unit.index)
        rows_av.append({"target": name, "model": "baseline", **sb, "gain_raw": 0.0, "gain_std": 0.0})
        rows_av.append({"target": name, "model": "granite", **sg, "gain_raw": sg["explained_raw"] - sb["explained_raw"],
                        "gain_std": sg["explained_std"] - sb["explained_std"]})  # fmt: skip
        say(f"   baseline        explained raw {sb['explained_raw']:6.2%}  std {sb['explained_std']:6.2%}")
        say(f"   granite         explained raw {sg['explained_raw']:6.2%}  std {sg['explained_std']:6.2%}")
        for feat in FEATS + ["control_shuffled", "ALL"]:
            extra = ftf[FEATS].to_numpy() if feat == "ALL" else ftf[[feat]].to_numpy()
            p = cv_predict(FSETS["baseline"], emb, temp, sat, y, cr, extra=extra)
            s, _ = scores(y, p)
            rows_av.append({"target": name, "model": f"baseline + {feat}", **s, "gain_raw": s["explained_raw"] - sb["explained_raw"],
                            "gain_std": s["explained_std"] - sb["explained_std"]})  # fmt: skip
            say(f"   + {feat:16s} explained raw {s['explained_raw']:6.2%} ({s['explained_raw'] - sb['explained_raw']:+.2%})  "
                f"std {s['explained_std']:6.2%} ({s['explained_std'] - sb['explained_std']:+.2%})")
        p = cv_predict(FSETS["granite"], emb, temp, sat, y, cr, extra=ftf[FEATS].to_numpy())
        s, _ = scores(y, p)
        rows_av.append({"target": name, "model": "granite + ALL", **s, "gain_raw": s["explained_raw"] - sg["explained_raw"],
                        "gain_std": s["explained_std"] - sg["explained_std"]})  # fmt: skip
        say(f"   granite + ALL    explained raw {s['explained_raw']:6.2%} ({s['explained_raw'] - sg['explained_raw']:+.2%} over granite)  "
            f"std {s['explained_std']:6.2%} ({s['explained_std'] - sg['explained_std']:+.2%})")
        if name == "clr_all":
            sc = pd.read_csv(REPO / "transfer" / "results" / "residual_modules_noA" / "scores.csv", index_col="tile_key").reindex(unit.index)
            axes = sc[["W1", "W2"]].copy()
            tiles_ft, tiles_cr = ft, cr

    # residual-axis correlations, within cruise
    axes["C_residual"] = resid_C.reindex(axes.index)
    say("\nWithin-cruise Spearman with residual axes (tiles; centred per cruise):")
    for feat in FEATS + ["control_shuffled"]:
        r = {"feature": feat}
        for a in ["W1", "W2", "C_residual"]:
            ok = tiles_ft[feat].notna() & axes[a].notna()
            x = tiles_ft.loc[ok, feat] - tiles_ft.loc[ok, feat].groupby(tiles_cr[ok.to_numpy()]).transform("mean")
            yv = axes.loc[ok, a] - axes.loc[ok, a].groupby(tiles_cr[ok.to_numpy()]).transform("mean")
            rho, pv = spearmanr(x, yv)
            r[f"{a}_rho"], r[f"{a}_p"], r[f"{a}_n"] = rho, pv, int(ok.sum())
        rows_cor.append(r)
        say(f"   {feat:16s} W1 {r['W1_rho']:+.3f} (p {r['W1_p']:.1e})  W2 {r['W2_rho']:+.3f} (p {r['W2_p']:.1e})  "
            f"C residual {r['C_residual_rho']:+.3f} (p {r['C_residual_p']:.1e})")
    av = pd.DataFrame(rows_av)
    av.round(5).to_csv(OUT / "added_value.csv", index=False)
    cor = pd.DataFrame(rows_cor)
    cor.round(5).to_csv(OUT / "axis_correlations.csv", index=False)
    (OUT / "summary.txt").write_text("\n".join(log) + "\n")
    figure(av, cor)


def figure(av, cor):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    labels = {"shelf_days": "Days since shelf", "sst_change30": "SST change, 30 d", "sst_sd60": "SST variability, 60 d",
              "chl_mean60": "Chlorophyll mean, 60 d", "chl_change30": "Chlorophyll change, 30 d", "chl_peak60": "Past chlorophyll peak",
              "mld_mean30": "MLD mean, 30 d", "mld_change30": "MLD change, 30 d", "light_ml0": "Mixed-layer light, day 0",
              "light_ml30": "Mixed-layer light, 30 d", "displacement_km": "Displacement, 60 d", "dlat_origin": "Latitude of origin shift",
              "control_shuffled": "Control (shuffled)", "ALL": "All features"}  # fmt: skip
    feats = FEATS + ["control_shuffled", "ALL"]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.0), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axs[0]
    cols = {"clr_all": "#000000", "C": "#E69F00", "A": "#0072B2"}
    names = {"clr_all": "Current CLR, all KOs", "C": "C  eukaryote / prokaryote", "A": "A  prokaryote KOs per GE"}
    for j, t in enumerate(["clr_all", "C", "A"]):
        g = av[av["target"] == t].set_index("model")["gain_raw"]
        v = [g.get(f"baseline + {x}", np.nan) * 100 for x in feats]
        ax.barh(np.arange(len(feats)) + (j - 1) * 0.26, v, 0.24, color=cols[t], label=names[t])
    ax.axvline(0, color="#666666", lw=0.5)
    ax.set_yticks(range(len(feats)))
    ax.set_yticklabels([labels[x] for x in feats])
    ax.invert_yaxis()
    ax.set_xlabel("Gain in variance explained over baseline\n(percentage points, cruise held out)")
    ax.legend(frameon=False, fontsize=6.5, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=1)
    ax = axs[1]
    c = cor.set_index("feature")
    fs = FEATS + ["control_shuffled"]
    mat = c.loc[fs, ["W1_rho", "W2_rho", "C_residual_rho"]].to_numpy()
    lim = max(0.15, np.nanmax(np.abs(mat)))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=6.5, color="white" if abs(mat[i, j]) > lim * 0.6 else "black")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Axis 1\n(W1)", "Axis 2\n(W2)", "C\nresidual"], fontsize=7)
    ax.set_yticks(range(len(fs)))
    ax.set_yticklabels([labels[x] for x in fs])
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("Spearman ρ within cruise")
    fig.tight_layout()
    fig.savefig(FIG / "water_mass_test.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
