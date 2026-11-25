#!/usr/bin/env python3
"""R74: follow-ups to R73. (1) The C winner (Granite PC5-32 + SST + MLD30) and its parts applied to all KOs, with
per-KO tables by read domain. (2) What Granite PC1-4 encode: correlation (pooled and within cruise) of each PC
(PCA 64 on the standardised embedding of all 1,215 tiles, as in the models) with SST, MLD30, |latitude|, C, OLCI bands,
ocean-model fields and the platform group.
Same setup as combo_search.py (inputs.npz; ridge, cv penalty, 13 cruises held out).
Outputs: data/results/combo/followup_ko_models.csv, followup_ko_per_ko.csv, followup_pc_correlations.csv, followup.txt
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "scripts")]
import combo_search as cs  # noqa: E402

GROUPS = {"PC1-4": ("pc", 0, 4), "PC5-32": ("pc", 4, 32), "PC1-32": ("pc", 0, 32), "SST": ("sat", 16, 17), "MLD30": ("mld",),
          "OLCI": ("sat", 0, 16), "nutrients + biology": ("field", cs.NUTRIENTS + cs.BIO),
          "physics + history": ("field", cs.PHYSICS + cs.HISTORY)}  # fmt: skip
MODELS = {"Granite (PC1-32)": ["PC1-32"], "PC1-4": ["PC1-4"], "PC5-32": ["PC5-32"], "SST + MLD30": ["SST", "MLD30"],
          "Granite + SST + MLD30": ["PC1-32", "SST", "MLD30"], "PC5-32 + SST + MLD30": ["PC5-32", "SST", "MLD30"],
          "PC5-32 + SST + MLD30 + nutrients + physics": ["PC5-32", "SST", "MLD30", "nutrients + biology", "physics + history"],
          "Granite + SST + MLD30 + nutrients + physics": ["PC1-32", "SST", "MLD30", "nutrients + biology", "physics + history"]}  # fmt: skip
LOG: list = []


def say(s=""):
    print(s, flush=True)
    LOG.append(str(s))


def ko_models(z):
    y, cr = z["yK"].astype(np.float64), z["cruise"]
    ss = ((y - y.mean(0)) ** 2).sum(0)
    folds = []
    for c in np.unique(cr):
        te, tr = np.where(cr == c)[0], np.where(cr != c)[0]
        folds.append((te, tr, cs.group_features(z, tr, te, GROUPS)))
    held = {}
    for m, gs in MODELS.items():
        h = np.zeros((len(folds), y.shape[1]))
        for i, (te, tr, f) in enumerate(folds):
            p, _ = cs.ridge_cv(np.hstack([f[g][0] for g in gs]), y[tr], np.hstack([f[g][1] for g in gs]), cr[tr])
            h[i] = ((y[te] - p) ** 2).sum(0)
        held[m] = h
        say(f"  {m:45s} all KOs {100 * (1 - h.sum() / ss.sum()):6.2f}%  per-KO median R2 {np.median(1 - h.sum(0) / ss):.4f}")
    return y, ss, held


def main():
    z = np.load(cs.OUT / "inputs.npz", allow_pickle=True)
    say("All KOs (3,933), % of non-technical variance, cruise held out, fixed models (no selection):")
    y, ss, held = ko_models(z)
    tab = pd.DataFrame({m: [100 * (1 - h.sum() / ss.sum())] for m, h in held.items()}).T.rename(columns={0: "explained"})
    tab.round(3).to_csv(cs.OUT / "followup_ko_models.csv")
    dom = pd.read_csv(REPO / "data/results/ko_domain/ko_domain_reads.csv", index_col="ko")["euk_read_share"].reindex(z["kos"])
    d = np.select([dom > 0.9, dom <= 0.1], ["eukaryotic", "prokaryotic"], "mixed")
    per = pd.DataFrame({"ko": z["kos"], "domain": d})
    for m, h in held.items():
        per[m] = 1 - h.sum(0) / ss
    ref, alt = "Granite (PC1-32)", "PC5-32 + SST + MLD30"
    per["gain"] = per[alt] - per[ref]
    per["cruises_better"] = (held[alt] < held[ref]).sum(0)
    per.round(4).to_csv(cs.OUT / "followup_ko_per_ko.csv", index=False)
    say(f"\n{alt} vs {ref}, per KO by read domain (pooled explained of the domain's KOs; median per-KO R2; gains):")
    for dm in ("prokaryotic", "mixed", "eukaryotic"):
        k = (d == dm)
        a = {m: 100 * (1 - held[m][:, k].sum() / ss[k].sum()) for m in (ref, "Granite + SST + MLD30", alt)}
        g = per.loc[k, "gain"]
        say(f"  {dm:12s} n {k.sum():5,}  pooled: Granite {a[ref]:.2f}, Granite + SST + MLD30 {a['Granite + SST + MLD30']:.2f}, "
            f"{alt} {a[alt]:.2f};  per-KO gain median {g.median():+.4f}, >= +0.05 {(g >= .05).sum():,} "
            f"(and >= 10 cruises {((g >= .05) & (per.loc[k, 'cruises_better'] >= 10)).sum():,}), <= -0.05 {(g <= -.05).sum():,}")  # fmt: skip

    # what the first PCs encode (descriptive; PCA on all tiles)
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    e = StandardScaler().fit_transform(z["emb"].astype(float))
    pca = PCA(n_components=64, random_state=0).fit(e)
    P = pca.transform(e)
    say(f"\nPCA on the embedding (all tiles): variance share PC1-4 {pca.explained_variance_ratio_[:4].sum():.3f} "
        f"({', '.join(f'{v:.3f}' for v in pca.explained_variance_ratio_[:4])}); PC5-32 {pca.explained_variance_ratio_[4:32].sum():.3f}")  # fmt: skip
    tiles = pd.read_csv(REPO / "data/results/benchmark_export/tiles.csv", index_col="tile_key").loc[z["tile"]]
    sat = pd.DataFrame(z["sat"].astype(float), columns=cs.OLCI + ["SST"])
    sat = sat.fillna(sat.median())
    X = pd.concat([sat[["SST"]], pd.DataFrame({"MLD30": z["mld"].astype(float), "|latitude|": np.abs(tiles["sample_lat"].to_numpy()),
                   "C (adjusted)": z["yC"][:, 0], "platform NovaSeq 6000": (tiles["platform_group"] == "NovaSeq 6000 group").astype(float).to_numpy()}),
                   pd.DataFrame(z["fields"].astype(float), columns=list(z["field_names"])), sat[cs.OLCI]], axis=1)  # fmt: skip
    cr = z["cruise"]

    def corr(a, b, within):
        a, b = pd.Series(a), pd.Series(b)
        if within:
            a, b = a - a.groupby(cr).transform("mean"), b - b.groupby(cr).transform("mean")
        return np.corrcoef(a, b)[0, 1]

    rows = []
    for j in range(8):
        for v in X.columns:
            rows.append({"pc": f"PC{j + 1}", "variable": v, "r_pooled": corr(P[:, j], X[v].to_numpy(), False),
                         "r_within_cruise": corr(P[:, j], X[v].to_numpy(), True)})  # fmt: skip
    t = pd.DataFrame(rows)
    t.round(3).to_csv(cs.OUT / "followup_pc_correlations.csv", index=False)
    say("\nPC correlations (pooled / within cruise), strongest 6 per PC:")
    for pc, g in t.groupby("pc", sort=False):
        g = g.reindex(g["r_pooled"].abs().sort_values(ascending=False).index).head(6)
        say(f"  {pc}: " + "; ".join(f"{r.variable} {r.r_pooled:+.2f} / {r.r_within_cruise:+.2f}" for r in g.itertuples()))
    # how much of SST, MLD30, |lat|, C each PC block reproduces linearly (in-sample R2)
    say("\nIn-sample linear R2 of each variable from PC blocks (descriptive):")
    for v in ("SST", "MLD30", "|latitude|", "C (adjusted)", "platform NovaSeq 6000", "nitrate", "GLORYS temperature"):
        yv = X[v].to_numpy()
        out = []
        for nm, sl in (("PC1-4", slice(0, 4)), ("PC5-32", slice(4, 32)), ("PC1-32", slice(0, 32))):
            A = np.hstack([np.ones((len(yv), 1)), P[:, sl]])
            b, *_ = np.linalg.lstsq(A, yv, rcond=None)
            out.append(f"{nm} {1 - ((yv - A @ b) ** 2).sum() / ((yv - yv.mean()) ** 2).sum():.3f}")
        say(f"  {v:22s} " + ", ".join(out))
    (cs.OUT / "followup.txt").write_text("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
