#!/usr/bin/env python3
"""External test on the JGI marine metagenomes (Flinkstrom et al. 2024; RESULTS.md R29).

Models are trained on all Bio-GO-SHIP tiles (R17 units: fixed position, OSTIA-filled
embeddings, samples on identical tiles merged) and predict the JGI tiles. Targets: the
2,832 KOs shared with JGI, CLR over that set in both datasets (scripts/build_external_targets.py).

The JGI samples have no in-situ temperature, so "temperature" here is the tile's own SST
(SLSTR, or OSTIA where SLSTR is missing; Granite-normalised). Models:
    mean          the Bio-GO-SHIP mean per KO (no input)
    granite       ridge on the 512-d embedding
    pca32         ridge on 32 Granite PCs
    sst           spline of tile SST
    pca32 + sst   the best model (R18) with tile SST in place of in-situ temperature
    raw satellite spline of 16 OLCI band means + SST (rank-normalised)
    time+place    ridge on latitude, sin/cos longitude, sin/cos day of year
Reference: the same models, cruise-held-out on Bio-GO-SHIP, plus pca32 + in-situ temperature,
so the cost of using tile SST is visible.

Per KO on JGI: pooled R2 against the JGI mean; "shift R2" after removing the mean error
(one constant per KO, fitted on JGI, so slightly optimistic: it forgives a dataset offset
such as annotation or extraction differences); Pearson r. Chance level: the expected number
of KOs with shift R2 >= 0.3 when the same predictions are given to permuted tiles (200 draws). Also on Monterey Bay alone (one
position, Sep-Nov 2016: change over time) and on the other JGI tiles (spatial contrast).

Outputs in data/results/external_jgi/: external_kos.csv (KO x model), summary.csv,
units.csv (the JGI tiles used).

Usage:
    python scripts/external_test_jgi.py
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
from biogoship import tiles  # noqa: E402
from biogoship.tables import REPO_ROOT  # noqa: E402
from screen_ko_predictability import ALPHAS, load_embeddings, metadata_features  # noqa: E402
from test_ko_temperature import satellite_features  # noqa: E402
from tune_figure_kos import build, load, pooled_r2  # noqa: E402

EXT = REPO_ROOT / "data" / "external" / "jgi_flinkstrom"
RESULTS = REPO_ROOT / "data" / "results"
MODELS = {"mean": [], "granite": ["E"], "pca32": ["P32"], "sst": ["SST~"], "pca32 + sst": ["P32", "SST~"],
          "raw satellite": ["S~"], "time+place": ["M"]}  # fmt: skip
CV_ONLY = {"pca32 + in-situ temperature": ["P32", "T~"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", type=Path, default=EXT / "targets" / "ko_clr_shared_biogoship.csv")
    parser.add_argument("--external-targets", type=Path, default=EXT / "targets" / "ko_clr_shared.csv")
    parser.add_argument("--embeddings-dirs", type=Path, nargs="+",
                        default=[REPO_ROOT / "data" / "embeddings_sst_l4", REPO_ROOT / "data" / "faire" / "embeddings_sst_l4"])
    parser.add_argument("--external-embeddings", type=Path, default=EXT / "embeddings_sst_l4")
    parser.add_argument("--categories", type=Path, default=RESULTS / "ko_temperature_all_ostia_v2" / "ko_categories.csv")
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "external_jgi")
    args = parser.parse_args()

    kos = pd.read_csv(args.targets, nrows=0).columns.drop("sample_key").tolist()
    emb, temp, sat, y, cruises, unit = load(args, kos)
    ok = ~np.isnan(sat[:, -1])  # tile SST present
    emb, temp, sat, y, cruises, unit = emb[ok], temp[ok], sat[ok], y[ok], cruises[ok], unit[ok]
    meta = metadata_features(unit)

    # JGI tiles: same tile choice as training (fixed, closest tile with SST), identical tiles merged
    ext_targets = pd.read_csv(args.external_targets, index_col="sample_id")[kos]
    rows, x_ext = load_embeddings(argparse.Namespace(embeddings_dirs=[args.external_embeddings], position="fixed", best_tile="prefer-sst"))
    rows = rows[rows["sample_id"].isin(ext_targets.index)]
    content = pd.Series([v.tobytes() for v in x_ext[rows["embedding_row"].to_numpy()].round(5)], index=rows.index)
    rows = rows.assign(tile_key=content.map({k: i for i, k in enumerate(content.unique())}))
    y_ext = ext_targets.loc[rows["sample_id"]].set_axis(rows.index).groupby(rows["tile_key"]).mean()
    u_ext = rows.drop_duplicates("tile_key").set_index("tile_key").loc[y_ext.index]
    sat_ext = satellite_features([args.external_embeddings], u_ext["tile_path"])
    keep = ~np.isnan(sat_ext[:, -1])
    u_ext, y_ext, sat_ext = u_ext[keep], y_ext.to_numpy()[keep], sat_ext[keep]
    e_ext, m_ext = x_ext[u_ext["embedding_row"].to_numpy()], metadata_features(u_ext)
    sst_c = lambda s: s * tiles.SST_STD_K + tiles.SST_MEAN_K - 273.15  # noqa: E731
    print(f"Bio-GO-SHIP: {len(y):,} tiles, {len(set(cruises))} cruises; JGI: {len(y_ext)} tiles from "
          f"{rows['sample_id'].nunique()} samples; {len(kos):,} KOs")
    print(f"Tile SST, deg C: Bio-GO-SHIP {np.nanpercentile(sst_c(sat[:, -1]), [0, 50, 100]).round(1)}, "
          f"JGI {np.percentile(sst_c(sat_ext[:, -1]), [0, 50, 100]).round(1)}")
    print("JGI tiles per section:", u_ext["section"].value_counts().to_dict())

    # Cruise-held-out reference on Bio-GO-SHIP
    cv = {}
    for name, blocks in {**MODELS, **CV_ONLY}.items():
        pred = np.full_like(y, np.nan)
        for c in np.unique(cruises):
            test, train = np.flatnonzero(cruises == c), np.flatnonzero(cruises != c)
            pred[test] = fit_predict(blocks, train, test, emb, temp, sat, meta, y)
        cv[name] = pooled_r2(y, pred)
        print(f"  CV {name:30s} median R2 {np.median(cv[name]):+.3f}  KOs >= 0.3: {(cv[name] >= 0.3).sum():5,}", flush=True)

    # Train on all Bio-GO-SHIP, predict JGI
    n = len(y)
    E, T = np.vstack([emb, e_ext]), np.vstack([temp, np.full((len(y_ext), 1), np.nan)])
    S, M = np.vstack([sat, sat_ext]), np.vstack([meta, m_ext])
    train, test = np.arange(n), np.arange(n, n + len(y_ext))
    groups = {"all": np.ones(len(y_ext), bool), "Monterey Bay": (u_ext["section"] == "JGI_MontereyBay").to_numpy()}
    groups["other JGI"] = ~groups["Monterey Bay"]
    out, rng = [], np.random.default_rng(0)
    for name, blocks in MODELS.items():
        pred = fit_predict(blocks, train, test, E, T, S, M, y)
        frame = {"ko": kos, "model": name, "cv_r2": cv[name]}
        for g, mask in groups.items():
            yt, pt = y_ext[mask], pred[mask]
            frame[f"r2_{g}"] = pooled_r2(yt, pt)
            frame[f"shift_r2_{g}"] = pooled_r2(yt - yt.mean(0), pt - pt.mean(0))
            frame[f"r_{g}"] = corr(yt, pt)
            # chance: the same predictions given to the wrong tiles (200 permutations)
            perm = [pooled_r2(yt - yt.mean(0), (q := pt[rng.permutation(len(pt))]) - q.mean(0)) >= 0.3 for _ in range(200)]
            frame[f"chance_p_{g}"] = np.mean(perm, axis=0)
        out.append(pd.DataFrame(frame))
    table = pd.concat(out)
    table["cv_r2_best_insitu"] = table["ko"].map(dict(zip(kos, cv["pca32 + in-situ temperature"])))
    table["bias"] = table["model"].map(lambda _: np.nan)
    table.loc[table["model"] == "mean", "bias"] = (y_ext.mean(0) - y.mean(0))  # JGI mean minus Bio-GO-SHIP mean
    cats = pd.read_csv(args.categories).set_index("ko")["category"]
    table["category"] = table["ko"].map(cats)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table.round(4).to_csv(args.out_dir / "external_kos.csv", index=False)
    u_ext.assign(sst_c=sst_c(sat_ext[:, -1]).round(2))[["sample_id", "section", "sample_time", "sample_lat", "sample_lon",
                                                        "tile_path", "match_hour_gap", "sst_source", "sst_c"]].to_csv(args.out_dir / "units.csv")  # fmt: skip

    # Summaries
    strata = {"all KOs": np.ones(len(kos), bool), "Bio-GO-SHIP CV R2 >= 0.3 (best model)": cv["pca32 + in-situ temperature"] >= 0.3}
    for cat in ("granite", "same", "temp"):
        strata[f"heatmap: {cat}"] = cats.reindex(kos).eq(cat).to_numpy()
    rows_out = []
    for name in MODELS:
        part = table[table["model"] == name]
        for s, mask in strata.items():
            p = part[mask]
            row = {"model": name, "kos": s, "n": int(mask.sum()), "cv_median_r2": p["cv_r2"].median()}
            for g in groups:
                row |= {f"{g}: median r2": p[f"r2_{g}"].median(), f"{g}: median shift r2": p[f"shift_r2_{g}"].median(),
                        f"{g}: median r": p[f"r_{g}"].median(), f"{g}: kos shift r2 >= 0.3": int((p[f"shift_r2_{g}"] >= 0.3).sum()),
                        f"{g}: chance kos": p[f"chance_p_{g}"].sum()}  # fmt: skip
            rows_out.append(row)
    summary = pd.DataFrame(rows_out)
    summary.round(4).to_csv(args.out_dir / "summary.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 50)
    pd.set_option("display.max_rows", 200)
    print("\nAll JGI tiles:")
    print(summary[["model", "kos", "n", "cv_median_r2", "all: median r2", "all: median shift r2", "all: median r",
                   "all: kos shift r2 >= 0.3", "all: chance kos"]].round(3).to_string(index=False))  # fmt: skip
    for g in ("Monterey Bay", "other JGI"):
        print(f"\n{g}:")
        print(summary[["model", "kos", f"{g}: median r2", f"{g}: median shift r2", f"{g}: median r",
                       f"{g}: kos shift r2 >= 0.3", f"{g}: chance kos"]].round(3).to_string(index=False))  # fmt: skip
    b = table.loc[table["model"] == "mean", "bias"]
    print(f"\nJGI mean - Bio-GO-SHIP mean per KO (CLR): median {b.median():+.3f}, |bias| median {b.abs().median():.3f}, "
          f"> 1 in {(b.abs() > 1).sum()} KOs; JGI SD median {np.median(y_ext.std(0)):.3f}, Bio-GO-SHIP SD median {np.median(y.std(0)):.3f}")


def fit_predict(blocks, train, test, emb, temp, sat, meta, y):
    if not blocks:
        return np.tile(y[train].mean(0), (len(test), 1))
    if blocks == ["M"]:
        s = StandardScaler().fit(meta[train])
        x_tr, x_te = s.transform(meta[train]), s.transform(meta[test])
    else:
        x_tr, x_te = build(blocks, train, test, emb, temp, sat, {})
    return RidgeCV(alphas=ALPHAS, alpha_per_target=True).fit(x_tr, y[train]).predict(x_te)


def corr(a, b):
    a, b = a - a.mean(0), b - b.mean(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a * b).sum(0) / np.sqrt((a**2).sum(0) * (b**2).sum(0))


if __name__ == "__main__":
    main()
