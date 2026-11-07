#!/usr/bin/env python3
"""Figure: KOs where the Granite embedding beats in-situ temperature, ties with it, or loses.

Heatmap rows are KOs in three blocks; columns are pooled R2 from the Granite
embedding (ridge), pooled R2 from in-situ temperature (spline + ridge), and their
difference. Each cruise is held out in turn (13 cruises, OSTIA-filled embeddings,
1,216 tiles). Within a block, rows are grouped by KEGG function (BRITE ko00001,
level B, merged into a few groups; a few KOs without a pathway are placed by
hand, see OVERRIDES). The right-hand column names the phylum that carries most
of the KO's genes in one example metagenome (A392013_008, "coverage0" output of
the KO pipeline), with its share.

Blocks:
  Granite clearly better   difference >= 0.15 and Granite ahead in >= 10 of 13
                           held-out cruises; the N_GRANITE with highest Granite R2
  Both high, about equal   both R2 >= 0.6 and |difference| < 0.05; top N_SAME
  Temperature clearly better  difference <= -0.1 and temperature ahead in >= 9 of
                           13 cruises; the N_TEMP with highest temperature R2

Inputs:  data/results/ko_temperature_all_ostia/ (scripts/test_ko_temperature.py --kos all
         --learners ridge spline --feature-sets temp embedding --embeddings-dirs
         data/embeddings_sst_l4 data/faire/embeddings_sst_l4 --save-predictions),
         data/reference/kegg_brite_ko00001.csv, data/reference/kegg_ko_list.tsv,
         --taxa: the example pipeline output's level_2_phylum_division_ko_counts.tsv
Outputs: figures/output/ko_heatmap.png and .pdf; data/results/ko_temperature_all_ostia/
         ko_heatmap_rows.csv (the rows shown) and ko_categories.csv (every KO)

Usage:
    python figures/make_ko_heatmap.py --taxa <unzipped A392013_008_meta>/level_2_phylum_division_ko_counts.tsv
"""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "data" / "results" / "ko_temperature_all_ostia_v2"
OUT = ROOT / "figures" / "output"

INK, MUTED, RULE, SURFACE = "#0b0b0b", "#52514e", "#d6d5d0", "#ffffff"
GRANITE, TEMP = "#2a78d6", "#c0457e"
N_GRANITE, N_SAME, N_TEMP = 28, 14, 12
R2_CMAP = LinearSegmentedColormap.from_list("r2", ["#f6f5f2", "#bdbab2", "#6d6a63", "#232220"])
DIFF_CMAP = LinearSegmentedColormap.from_list("diff", [TEMP, "#f6f5f2", GRANITE])
BLUE_CMAP = LinearSegmentedColormap.from_list("blue", ["#ffffff", "#9cc3ee", "#2a78d6", "#0c2f5c"])
# --columns inputs: four predictors side by side (feature set, learner, header); the learner is
# temperature and raw satellite get a smooth curve per input (spline), as their responses bend;
# time & place ridge (best median over all KOs; its spline overfits).
INPUTS = [("embedding", "ridge", "Granite\nembedding"), ("temp", "spline", "In-situ\ntemperature"),
          ("metadata", "ridge", "Time &\nplace"), ("satellite", "spline", "Raw\nsatellite")]  # fmt: skip

# KEGG level-B categories -> display group, in priority order for KOs listed under several.
GROUPS = [
    ("Membrane transport", "Transport"),
    ("Energy metabolism", "Energy & photosynthesis"),
    ("Amino acid metabolism", "Amino acids"),
    ("Metabolism of other amino acids", "Amino acids"),
    ("Nucleotide metabolism", "Nucleotides"),
    ("Carbohydrate metabolism", "Carbohydrate, lipid & cofactors"),
    ("Lipid metabolism", "Carbohydrate, lipid & cofactors"),
    ("Metabolism of cofactors and vitamins", "Carbohydrate, lipid & cofactors"),
    ("Metabolism of terpenoids and polyketides", "Carbohydrate, lipid & cofactors"),
    ("Xenobiotics biodegradation and metabolism", "Carbohydrate, lipid & cofactors"),
    ("Translation", "DNA, RNA & protein synthesis"),
    ("Replication and repair", "DNA, RNA & protein synthesis"),
    ("Protein families: genetic information processing", "DNA, RNA & protein synthesis"),
]
GROUP_ORDER = ["Energy & photosynthesis", "Transport", "Amino acids", "Nucleotides",
               "Carbohydrate, lipid & cofactors", "DNA, RNA & protein synthesis", "Other / uncharacterised"]  # fmt: skip
OVERRIDES = {  # KOs whose KEGG placement is missing or misleading
    "K00436": "Energy & photosynthesis",  # hoxH, bidirectional hydrogenase (cyanobacterial)
    "K13991": "Energy & photosynthesis",  # puhA, bacterial photosynthetic reaction centre
    "K08928": "Energy & photosynthesis",  # pufL, bacterial photosynthetic reaction centre
    "K16090": "Transport",  # fiu, siderophore receptor (TonB-dependent uptake)
    "K14330": "Carbohydrate, lipid & cofactors",  # acyl-ACP reductase, alkane synthesis
    "K01689": "Carbohydrate, lipid & cofactors",  # enolase (glycolysis)
    "K10808": "Nucleotides",  # RRM2, ribonucleotide reductase
}
SHORT = {  # shorter descriptions for row labels
    "K07015": "Rab family protein", "K14330": "fatty aldehyde-generating acyl-ACP reductase",
    "K06893": "uncharacterised protein", "K03810": "virulence factor MviM",
}  # fmt: skip


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--taxa", type=Path, help="level_2_phylum_division_ko_counts.tsv of the example pipeline output")
    parser.add_argument("--columns", choices=["difference", "inputs"], default="difference",
                        help="difference: Granite, temperature, their difference (ko_heatmap); inputs: Granite, temperature, "
                             "time & place, raw satellite, white to dark blue (ko_heatmap_inputs)")  # fmt: skip
    args = parser.parse_args()
    plt.rcParams.update({
        "font.family": ["Arial", "DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42, "svg.fonttype": "none",
        "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.bf": "Arial:bold",
    })  # fmt: skip

    ko = categories()
    taxa = phylum_share(args.taxa) if args.taxa else pd.DataFrame(columns=["phylum", "share"])
    ko = ko.join(taxa)
    scores = pd.read_csv(RESULTS / "ko_temperature.csv")
    for feature_set, learner, _ in INPUTS:
        ko[f"r2_{feature_set}"] = scores[(scores["feature_set"] == feature_set) & (scores["learner"] == learner)].set_index("ko")["r2_pooled"]
    ko.round(4).to_csv(RESULTS / "ko_categories.csv")
    blocks = [
        ("Granite clearly better", ko[ko["category"] == "granite"].nlargest(N_GRANITE, "granite"), "granite"),
        ("Both high, about equal", ko[ko["category"] == "same"].assign(m=lambda d: d[["granite", "temp"]].mean(axis=1))
         .nlargest(N_SAME, "m"), "same"),
        ("Temperature clearly better", ko[ko["category"] == "temp"].nlargest(N_TEMP, "temp"), "temp"),
    ]  # fmt: skip
    counts = ko["category"].value_counts()

    rows = []
    for label, block, key in blocks:
        block = block.assign(group=pd.Categorical(block["group"], GROUP_ORDER, ordered=True))
        block = block.sort_values(["group", "granite" if key != "temp" else "temp"], ascending=[True, False])
        rows.append((label, key, block))
    pd.concat([b.assign(block=k) for _, k, b in rows]).round(4).to_csv(RESULTS / "ko_heatmap_rows.csv")

    # Layout in inches: x positions of the columns, y = one unit per row.
    W = 7.5
    X_GROUP, X_BRACKET, X_LABEL = 0.10, 1.30, 3.92  # group text (left-aligned), bracket, row labels (right-aligned)
    if args.columns == "difference":
        X_CELLS, CELL = [3.98, 4.50, 5.17], 0.50  # Granite, temperature, difference
        columns = [("granite", R2_CMAP, "{:.2f}"), ("temp", R2_CMAP, "{:.2f}"), ("diff", DIFF_CMAP, "{:+.2f}")]
        headers = ["Granite\nembedding", "In-situ\ntemperature", "Difference\n(G − T)"]
    else:
        X_CELLS, CELL = [3.98, 4.50, 5.02, 5.54], 0.50
        columns = [(f"r2_{f}", BLUE_CMAP, "{:.2f}") for f, _, _ in INPUTS]
        headers = [h for _, _, h in INPUTS]
    X_TAXON = X_CELLS[-1] + CELL + 0.08
    n_rows = sum(len(b) for _, _, b in rows)
    gap = 2.4
    total = n_rows + gap * len(rows)
    row_h, head, foot = 0.155, 0.55, 1.3
    fig_h = total * row_h + head + foot
    fig = plt.figure(figsize=(W, fig_h))
    ax = fig.add_axes([0, foot / fig_h, 1, (total * row_h) / fig_h])
    ax.set_xlim(0, W)
    ax.set_ylim(total, 0)
    ax.axis("off")

    r2_norm, diff_norm = Normalize(0, 0.8), Normalize(-0.4, 0.4)
    y = 0.0
    for label, key, block in rows:
        y += gap
        ax.text(X_GROUP, y - 0.6, label, ha="left", va="bottom", fontsize=8.5, fontweight="bold", color=INK)
        ax.text(X_CELLS[-1] + CELL, y - 0.6, f"{len(block)} shown of {int(counts.get(key, 0)):,} KOs", ha="right",
                va="bottom", fontsize=6.5, color=MUTED)  # fmt: skip
        ax.plot([X_GROUP, W - 0.1], [y - 0.1, y - 0.1], color=RULE, lw=0.5)
        for group, members in block.groupby("group", observed=True, sort=False):
            g0, g1 = y, y + len(members)
            ax.plot([X_BRACKET, X_BRACKET], [g0 + 0.18, g1 - 0.18], color=RULE, lw=1.2)
            text = textwrap.fill(str(group), 16) if len(members) > 1 else str(group)
            ax.text(X_GROUP, (g0 + g1) / 2, text, ha="left", va="center", fontsize=6.5, color=MUTED, linespacing=1.1)
            for _, row in members.iterrows():
                for x, (column, cmap, fmt) in zip(X_CELLS, columns):
                    value = row[column]
                    colour = cmap((diff_norm if column == "diff" else r2_norm)(value))
                    ax.add_patch(plt.Rectangle((x + 0.015, y + 0.07), CELL - 0.03, 0.86, color=colour, lw=0))
                    dark = (0.299 * colour[0] + 0.587 * colour[1] + 0.114 * colour[2]) < 0.55
                    ax.text(x + CELL / 2, y + 0.5, fmt.format(value), ha="center", va="center", fontsize=6.3,
                            color=SURFACE if dark else INK)  # fmt: skip
                ax.text(X_LABEL, y + 0.5, rf"$\bf{{{row['gene']}}}$  {row['what']}", ha="right", va="center",
                        fontsize=6.5, color=INK)  # fmt: skip
                if isinstance(row.get("phylum"), str):
                    ax.text(X_TAXON, y + 0.5, f"{row['phylum']} {row['share']:.0%}", ha="left", va="center",
                            fontsize=6.3, color=MUTED)  # fmt: skip
                y += 1

    for x, text in zip(X_CELLS, headers):
        ax.text(x + CELL / 2, 0.3, text, ha="center", va="bottom", fontsize=6.8, color=INK)
    ax.text(X_CELLS[0] + CELL if args.columns == "difference" else X_CELLS[0] + 2 * CELL, -1.7, "Pooled R², cruise held out", ha="center", va="bottom", fontsize=6.5, color=MUTED)
    ax.text(X_TAXON, 0.3, "Main carrier, share of genes\n(one example metagenome)", ha="left", va="bottom", fontsize=7, color=INK)

    key_y, key_h = 0.88 / fig_h, 0.09 / fig_h
    key_axes = [fig.add_axes([X_CELLS[0] / W, key_y, 0.85 / W, key_h])]
    key_axes[0].imshow(np.linspace(0, 0.8, 200)[None, :], aspect="auto", norm=r2_norm, extent=[0, 0.8, 0, 1],
                       cmap=R2_CMAP if args.columns == "difference" else BLUE_CMAP)  # fmt: skip
    key_axes[0].set_xticks([0, 0.4, 0.8])
    key_axes[0].set_xticklabels(["≤0", "0.4", "0.8"])
    key_axes[0].set_xlabel("R²", fontsize=6.5, color=MUTED, labelpad=1)
    if args.columns == "difference":
        dax = fig.add_axes([(X_CELLS[2] - 0.05) / W, key_y, (CELL + 0.4) / W, key_h])
        dax.imshow(np.linspace(-0.4, 0.4, 200)[None, :], aspect="auto", cmap=DIFF_CMAP, norm=diff_norm, extent=[-0.4, 0.4, 0, 1])
        dax.set_xticks([-0.4, 0, 0.4])
        dax.set_xticklabels(["T", "0", "Granite"])
        dax.set_xlabel("ΔR², better:", fontsize=6.5, color=MUTED, labelpad=1)
        key_axes.append(dax)
    for a in key_axes:
        a.set_yticks([])
        a.tick_params(length=2, labelsize=6.3, colors=MUTED)
        for spine in a.spines.values():
            spine.set_visible(False)

    compared = ("the embedding alone (ridge) with temperature alone (spline)." if args.columns == "difference" else
                "the embedding alone (ridge), temperature alone (spline), time & place\n(latitude, longitude, day of year; ridge) "
                "and the raw satellite values (16 OLCI band means + SST mean of the tile, rank-normalised; spline).")  # fmt: skip
    fig.text(X_GROUP / W, 0.10 / fig_h,
             "Best model so far: ridge regression on 32 principal components of the OSTIA-filled Granite embedding + in-situ "
             "temperature\n(median pooled R² 0.34 on the 20-KO benchmark and 0.12 over all 3,921 KOs, cruise held out; "
             "RESULTS.md R12, R17).\nThe heatmap compares " + compared,
             ha="left", va="bottom", fontsize=6.3, color=MUTED, linespacing=1.25)  # fmt: skip
    name = "ko_heatmap" if args.columns == "difference" else "ko_heatmap_inputs"
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{suffix}", dpi=300)
    print(f"Wrote {OUT / name}.png and .pdf; categories: " + ", ".join(f"{k} {v:,}" for k, v in counts.items()))


def categories() -> pd.DataFrame:
    scores = pd.read_csv(RESULTS / "ko_temperature.csv")
    per_cruise = pd.read_csv(RESULTS / "ko_temperature_per_cruise.csv")
    pick = {"granite": ("embedding", "ridge"), "temp": ("temp", "spline")}
    out = pd.DataFrame({k: scores[(scores["feature_set"] == f) & (scores["learner"] == l)].set_index("ko")["r2_pooled"]
                        for k, (f, l) in pick.items()})  # fmt: skip
    cruise = {k: per_cruise[(per_cruise["feature_set"] == f) & (per_cruise["learner"] == l)].set_index(["ko", "held_out"])["r2"]
              for k, (f, l) in pick.items()}  # fmt: skip
    out["cruises_granite_ahead"] = (cruise["granite"] > cruise["temp"]).groupby(level=0).sum()
    out["diff"] = out["granite"] - out["temp"]
    n = per_cruise["held_out"].nunique()
    out["category"] = np.select(
        [(out["diff"] >= 0.15) & (out["cruises_granite_ahead"] >= 10),
         (out["granite"] >= 0.6) & (out["temp"] >= 0.6) & (out["diff"].abs() < 0.05),
         (out["diff"] <= -0.1) & (out["cruises_granite_ahead"] <= n - 9)],
        ["granite", "same", "temp"], "")  # fmt: skip

    names = pd.read_csv(ROOT / "data/reference/kegg_ko_list.tsv", sep="\t", header=None, names=["ko", "name"]).set_index("ko")["name"]
    brite = pd.read_csv(ROOT / "data/reference/kegg_brite_ko00001.csv")
    brite_names = brite.drop_duplicates("ko").set_index("ko")
    out["gene"], out["what"] = zip(*[label(k, names, brite_names) for k in out.index])
    rank = {b: i for i, (b, _) in enumerate(GROUPS)}
    mapping = dict(GROUPS)
    by_ko = brite[brite["B"].isin(mapping)].assign(r=lambda d: d["B"].map(rank)).sort_values("r").drop_duplicates("ko")
    out["group"] = by_ko.set_index("ko")["B"].map(mapping).reindex(out.index).fillna("Other / uncharacterised")
    for k, g in OVERRIDES.items():
        if k in out.index:
            out.loc[k, "group"] = g
    return out


def label(ko: str, names: pd.Series, brite: pd.DataFrame) -> tuple[str, str]:
    name = names.get(ko)
    if not isinstance(name, str):
        name = ko + "; " + SHORT.get(ko, "")
    gene, _, what = name.partition(";")
    gene = gene.split(",")[0].strip().replace("_", r"\_")
    what = SHORT.get(ko, what.split("[EC")[0].strip())
    if len(what) > 46:
        what = what[:44].rstrip(" ,") + "…"
    return gene, what


def phylum_share(path: Path) -> pd.DataFrame:
    t = pd.read_csv(path, sep="\t")
    t["ko"] = t["ko"].str.replace("ko:", "", regex=False)
    t["phylum"] = t["taxon"].str.split(";").str[-1].replace({"Cyanobacteriota": "Cyanobacteria", "dsDNA viruses, no RNA stage": "dsDNA viruses"})
    t.loc[t["taxon"].str.startswith("Eukaryota"), "phylum"] = "Eukaryotes"
    g = t.groupby(["ko", "phylum"])["count"].sum()
    share = (g / g.groupby(level=0).sum()).rename("share").reset_index().sort_values("share").groupby("ko").tail(1)
    return share.set_index("ko")[["phylum", "share"]]


if __name__ == "__main__":
    main()
