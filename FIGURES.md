# Figures

One entry per figure: what it shows, the code that makes it, its inputs, and how
to regenerate it. Outputs go to `figures/output/` (committed). Inputs under
`data/` are not in git; they live on the cluster or are copied back as noted.

Shared conventions: Arial, no titles baked into the image (captions go in the
text), observed values in black, model-derived values in blue (#2a78d6) and dashed.

The figure scripts need numpy, pandas, matplotlib and scikit-learn. Run them with
a Python that has these; on the project owner's laptop that is
`/Users/hellpark/miniforge3/bin/python` (the default `python3` there lacks
matplotlib).

## Contents

| Figure | Script | Output | Result entry |
|---|---|---|---|
| Predictability: in-situ variables and KOs | `figures/make_predictability.py` | `figures/output/predictability.png`, `.pdf` | RESULTS.md R1.3, R3.1 |
| Spectral recovery from Granite embeddings | `figures/make_spectral_recovery.py` | `figures/output/spectral_recovery.png`, `.pdf` | METHODS_RESULTS.md section 7 |
| Natural-colour tiles (companion) | `figures/make_tile_rgb.py` | `figures/output/spectral_recovery_tile_rgb.png` | - |
| Embedding vs in-situ temperature, all KOs | `figures/make_temperature_gap.py` | `figures/output/temperature_gap.png`, `.pdf` | RESULTS.md R6 |
| Three KOs along a cruise line (P18) | `figures/make_cruise_line.py` | `figures/output/cruise_line.png`, `.pdf` | RESULTS.md R16 |
| KO heatmap: Granite vs temperature by function | `figures/make_ko_heatmap.py` | `figures/output/ko_heatmap.png`, `.pdf` | RESULTS.md R16 |
| Cruise line: temperature vs best model (P18) | `figures/make_cruise_line.py --compare best` | `figures/output/cruise_line_best.png`, `.pdf` | RESULTS.md R17 |
| KO heatmap: four predictors | `figures/make_ko_heatmap.py --columns inputs` | `figures/output/ko_heatmap_inputs.png`, `.pdf` (`--cmap turbo`: `ko_heatmap_inputs_turbo`) | RESULTS.md R17 |
| Cruise map: three KOs as coloured tracks (P18) | `figures/make_cruise_map.py` | `figures/output/cruise_map.png`, `.pdf` | RESULTS.md R18, R20 |
| Prediction skill per sample on the globe; overlap vs similar conditions | `figures/make_skill_map.py` | `figures/output/skill_map.png`, `.pdf` | RESULTS.md R20 |
| Cruise line: temperature vs best model per gene (P18) | `figures/make_cruise_line.py --compare perko` | `figures/output/cruise_line_perko.png`, `.pdf` | RESULTS.md R18 |
| Presence of patchy genes (hurdle models; superseded by hurdle_base) | `scripts/hurdle_kos.py --parts summary figure` | `figures/output/hurdle_kos.png`, `.pdf` | RESULTS.md R65 |
| Why patchy genes are satellite or technical (superseded by hurdle_base) | `scripts/hurdle_mechanism.py` | `figures/output/hurdle_mechanism.png`, `.pdf` | RESULTS.md R66 |
| Patchy genes, base model (Granite + depth + platform) | `scripts/hurdle_base.py --parts summary figure` | `figures/output/hurdle_base.png`, `.pdf` | RESULTS.md R69 |

Colours shared by the prediction figures: Granite embedding blue `#2a78d6`, raw
band averages orange `#eb6834`, position + season aqua `#1baf7a`, chance grey
`#9a9993`, in-situ temperature magenta `#c0457e`; fixed position filled markers, advected hollow. Set at the top of
`figures/make_predictability.py` (`GRANITE`, `RAW`, `META`, `NULL`).

## Predictability: in-situ variables and KOs

![Predictability](figures/output/predictability.png)

- **Shows:**
  - **a** In-situ temperature, phosphate and nitrate + nitrite, original 5
    cruises, one cruise held out at a time: R2 per held-out cruise (small dots)
    and mean +/- s.d. over cruises (large marker and line; value at right) for the
    Granite embedding and the raw band averages at the fixed (filled) and advected
    (hollow) position, and for position + season. Values below -1 are drawn at -1
    as left-pointing markers (raw bands on I07N, about -3).
  - **b** The 20 KOs best predicted by the embedding over 13 cruises, labelled
    with KEGG gene symbol, name and KO: pooled R2 of the embedding (blue) against
    position + season (aqua).
  - **c** Pooled R2 of all 3,921 screened KOs for the embedding, position + season
    and the shuffled embedding (chance). Values below -0.6 drawn at -0.6.
- **Settings:** option A (closest-in-time tile with SST; tiles without SST not
  used), unique tiles, fixed position for b and c. Full settings: RESULTS.md R1.3
  (panel a) and R3.1 (panels b, c).
- **Code:** `figures/make_predictability.py`. Where to edit:
  - panel a: `panel_insitu()`; rows in `FEATURE_SETS` and `TARGETS`; floor `R2_FLOOR`.
  - panel b: `panel_top_kos()`; number of KOs `N_TOP_KOS`; label text `label_ko()`.
  - panel c: `panel_ko_distribution()`; bins and clipping in its first lines.
  - layout (panel positions, figure size): the `fig.add_axes` calls in `main()`.
- **Inputs:**
  - `data/results/sst_only_prefer/feature_set_comparison_section_unique_tiles.csv`
    from `python scripts/compare_feature_sets.py --require-sst --best-tile prefer-sst --out-dir data/results/sst_only_prefer`
  - `data/results/ko_screen/ko_screen.csv` from `python scripts/screen_ko_predictability.py`
    (needs `data/targets/ko_clr.csv` from `python scripts/build_ko_targets.py`)
  - `data/reference/kegg_ko_list.tsv`, KEGG KO names, downloaded from
    https://rest.kegg.jp/list/ko on first run (copy of 2026-11-04 in place).
- **Regenerate:** `python figures/make_predictability.py`
- **To switch to option B** (closest tile, SST-only): point `INSITU` at
  `data/results/sst_only_closest/...`; for the KO panels re-run the screen with
  `--best-tile closest --out-dir data/results/ko_screen_closest` and point `KO_SCREEN` there.

## Embedding vs in-situ temperature, all KOs

![Embedding vs temperature](figures/output/temperature_gap.png)

- **Shows:**
  - **a** Per KO (3,921): pooled R2 from in-situ temperature alone (x) against the
    Granite embedding (y), cruises held out in turn; dashed line = equal. Blue:
    embedding clearly better; magenta: temperature clearly better; grey: neither.
    "Clearly" = pooled R2 differs by at least 0.1 and the same side wins in at
    least 9 of the 13 held-out cruises. Values below -0.6 drawn at -0.6.
  - **b** Distribution of R2(embedding) - R2(temperature) over all KOs; bars
    beyond +/-0.1 coloured, counts beside them (by difference only, without the
    9-of-13 rule).
- **Settings:** RESULTS.md R6 (fixed position, option A tiles, samples with an
  in-situ temperature; temperature: best of ridge, spline + ridge, kNN;
  embedding: best of ridge, kNN).
- **Code:** `figures/make_temperature_gap.py`. Where to edit:
  - what counts as a clear win: `WIN_MARGIN`, `WIN_CRUISES`.
  - colours: `GRANITE`, `TEMP` (`#c0457e`), `NEUTRAL` at the top.
  - panel a: `panel_scatter()`; panel b: `panel_histogram()` (bins in its first line).
  - layout: the `fig.add_axes` calls in `main()`.
- **Inputs:** `data/results/ko_temperature_all/ko_temperature.csv` and
  `ko_temperature_per_cruise.csv` from
  `python scripts/test_ko_temperature.py --kos all --learners ridge spline knn --feature-sets temp embedding --out-dir data/results/ko_temperature_all`.
  Also writes `data/results/ko_temperature_all/ko_temperature_gap.csv` (per KO:
  both R2, difference, cruises won, verdict).
- **Regenerate:** `python figures/make_temperature_gap.py`

## Spectral recovery from Granite embeddings

![Observed vs recovered spectra](figures/output/spectral_recovery.png)

- **Shows:** for four contrasting samples, the tile's observed mean water-leaving
  reflectance at each OLCI wavelength (black) and the same spectrum recovered
  linearly from the tile's 512-number Granite embedding (blue, dashed), with
  observed and recovered SST. Evidence that the embeddings retain the absolute
  colour and temperature of the water.
- **Method:** ridge regression from embedding to each channel's tile mean, over
  all 3,503 tiles, cross-validated with whole cruises held out (5 folds grouped by
  cruise); the plotted recovery for each tile comes from a model that never saw
  its cruise. Samples: among best fixed tiles at least 95% clear with SST and
  well-behaved visible reflectances (all positive, 560 nm above 0.002), the
  clearest and greenest (highest and lowest 443/560 nm ratio), the warmest and the
  coldest. Panel d shows the weakest case: at 0 C the recovered SST is 6.6 C,
  because very cold tiles are rare (Granite also excluded polar water from
  pretraining).
- **Code:** `figures/make_spectral_recovery.py`
- **Inputs:** `data/embeddings/` tile_embeddings.npy, tile_embeddings.csv,
  tile_raw_features.csv, sample_embeddings.csv (copied back from the cluster).
- **Regenerate:** `python figures/make_spectral_recovery.py`
- **Outputs:** `figures/output/spectral_recovery.png`, `.pdf`, and
  `spectral_recovery_tiles.csv` (the four tiles: sample, scene, date, position,
  tile path).

### Companion: natural-colour images of those four tiles

![Natural-colour tiles](figures/output/spectral_recovery_tile_rgb.png)

- **Shows:** the four 42 x 42-pixel tiles (12.6 km across) behind panels a-d, in
  natural colour with one shared stretch, so differences in colour are real:
  deep blue for the South Pacific gyre (a), blue-green off southern California
  (b), lighter blue for the Bay of Bengal (c) and the Antarctic margin (d). Flagged
  pixels are grey (a few at the left edge of d). The pixel-to-pixel speckle is
  sensor noise at 300 m; open-ocean tiles have little real structure at this
  scale.

- **Code:** `figures/make_tile_rgb.py figures/output/spectral_recovery_tiles.csv`
- **Inputs:** the four tile files listed in `spectral_recovery_tiles.csv`, copied
  from the cluster to the same relative paths under `data/tiles/`.
- **Output:** `figures/output/spectral_recovery_tile_rgb.png` (1 x 4 strip, shared
  colour stretch, flagged pixels grey).
- The full scenes (with a quick-look JPEG inside each archive) are the products
  named in the `granule_id` column, in `data/products/olci/` on the cluster.

## Three KOs along a cruise line (P18)

![Cruise line](figures/output/cruise_line.png)

- **Shows:** P18 (eastern Pacific, 2016-11 to 2017-01, 64 tiles) held out, models
  trained on the other 12 cruises (OSTIA-filled embeddings). **a** track; **b** in-situ
  temperature; **c-e** observed CLR (black) with predictions from the Granite
  embedding (ridge, blue) and from temperature (spline, magenta dashed), and the R2
  on this cruise.
  - cpeS (*Synechococcus* phycoerythrin linker): both about 0.9.
  - glsA (glutaminase): Granite 0.50, temperature -0.15. It picks up the subtropical
    gyre maximum near 25 S, which temperature misses.
  - fiu (siderophore receptor): temperature 0.52, Granite 0.32.
  Lines are broken where samples are more than 8 degrees of latitude apart.
- **Inputs:** `data/results/ko_temperature_all_ostia_v2/` (RESULTS.md R16, R17; needs `--save-predictions`).
- **Regenerate:** `python figures/make_cruise_line.py` (`--cruise`, `--kos` for others).

## KO heatmap: Granite vs temperature by function

![KO heatmap](figures/output/ko_heatmap.png)

- **Shows:** pooled R2 (13 cruises, each held out) from the Granite embedding (ridge)
  and from in-situ temperature (spline), and their difference, for three blocks of KOs.
  - Granite clearly better: the 28 best of 87 (Granite R2 > 0.3, difference >= 0.15, ahead in >= 10 cruises).
  - Both high, about equal: the 14 best of 37 (both >= 0.6, |difference| < 0.05).
  - Temperature clearly better: the 12 best of 26 (temperature R2 > 0.3, difference <= -0.1, behind in >= 9 cruises).
  The R2 > 0.3 floor (2026-11-11) leaves out KOs neither input predicts; without it the
  blocks held 360 and 280 KOs, and most of the 280 had temperature R2 near 0.
  Rows are grouped by KEGG BRITE level B (merged; a few KOs placed by hand).
  The right column gives the pooled R2 of the best model so far, 32 Granite PCs + a
  temperature spline (ridge; the R28 baseline). Until 2026-11-12 it gave the main
  carrier phylum in one example metagenome; that was dropped as misleading after R22
  (the KO split follows the bacterial share of the community, not one carrier).
  The footnote names the best model so far (R12).
- **Inputs:** as above, plus `data/reference/kegg_brite_ko00001.csv` (from
  https://rest.kegg.jp/get/br:ko00001/json) and `data/results/model_fields/added_all_kos.csv`
  (`scripts/added_fields_all_kos.py`, R28).
- **Regenerate:** `python figures/make_ko_heatmap.py`

## Cruise line: temperature vs best model (P18)

![Cruise line, best model](figures/output/cruise_line_best.png)

- **Shows:** as the cruise line above, but the blue line is the best model so far
  (ridge on 32 PCs of the embedding + in-situ temperature, PCA fitted on the
  training cruises). It shows what Granite adds on top of temperature. P18 R2,
  best model vs temperature alone:
  - cpeS: 0.89 vs 0.92.
  - glsA: 0.47 vs -0.15.
  - fiu: 0.30 vs 0.52. The best model takes temperature as a straight line, so it
    misses fiu's curved response, which the temperature spline catches.
- **Regenerate:** `python figures/make_cruise_line.py --compare best`

## KO heatmap: four predictors

![KO heatmap, four predictors](figures/output/ko_heatmap_inputs.png)

- **Shows:** the same rows as the KO heatmap. Columns are pooled R2 from the Granite
  embedding (ridge), in-situ temperature (spline), time & place (latitude, sin/cos
  longitude, sin/cos day of year; ridge) and the raw satellite values (16 OLCI band
  means + SST mean of the tile, rank-normalised; spline). The colour scale runs 0 to 1
  through several hues (pale yellow, green, blue, navy), because a single blue made
  0.4 and 0.8 hard to tell apart; `--cmap turbo` gives a rainbow version and `--cmap
  blue` the first one. Values at or below 0 take the lowest colour.
- **Fourth block (2026-11-12), "Improved by adding an ocean-model field":** KOs where
  adding one GLORYS / PISCES field or recent-history summary to 32 PCs + temperature
  gains >= 0.05 R2, reaches R2 > 0.3, and is better in >= 10 of 13 held-out cruises (R28):
  the 16 largest gains of 93 such KOs, leaving out KOs already shown. Rows are grouped by
  the field. In this block the second column is 32 PCs + temperature + that field (not
  temperature), and the right column gives baseline -> with the field (gain). The 16 are
  8 phosphate (phoU +0.38, adeR +0.31, urtD, ureA, ureC, phoB, E3.1.3.70, vicR), 7 nitrate
  (ureF +0.34, ureD, urtA, urtB, urtC, ureB, gpgS) and 1 salinity (ttuC +0.19). Random
  inputs pass the same filter for 0 KOs. Rows: `data/results/ko_temperature_all_ostia_v2/ko_heatmap_field_rows.csv`
  (kept apart from `ko_heatmap_rows.csv`, which other scripts use as the figure KOs).
- **Regenerate:** `python figures/make_ko_heatmap.py --columns inputs`


## Cruise line: temperature vs best model per gene (P18)

- **Shows:** as the cruise line, but the blue line is the model chosen for each gene by
  `scripts/tune_figure_kos.py` (RESULTS.md R18), named under each gene. The choice is
  nested: made from the other 12 cruises only, so P18 does not pick its own model.
- **Regenerate:** `python scripts/tune_figure_kos.py` (about 15 minutes), then
  `python figures/make_cruise_line.py --compare perko`

## Cruise map: three KOs as coloured tracks (P18)

- **Shows:** one map per KO (cpeS, pstC, NOP2). Two copies of the cruise track side by
  side (really they overlap): observed CLR and the held-out prediction of a Granite
  model chosen per KO by pooled R2 (R20): 32 PCs + temperature spline for cpeS and NOP2;
  phosphate spline first, then 32 PCs on the residuals for pstC (tiles with in-situ
  phosphate only, so one P18 sample has no prediction). One colour scale per KO
  (viridis, 2nd-98th percentile of observed). R2 on P18 under each map, with
  temperature alone for reference.
- **Regenerate:** `python scripts/tune_figure_kos.py`, `python scripts/phosphate_models.py`,
  then `python figures/make_cruise_map.py`

## Prediction skill per sample on the globe

- **Shows:** a, every tile coloured by its skill with its cruise held out (1 - mean over
  the 920 KOs with pooled R2 >= 0.25 of error^2 / the KO's variance over all tiles;
  32 PCs + temperature, ridge), diverging red-grey-blue around 0; table of cruises with
  median skill and the share of samples within 200 km of another cruise. b, skill vs
  distance to the nearest sample of another cruise. c, skill vs how unusual the
  sample's conditions are for the other cruises (5 nearest in 8 embedding PCs +
  temperature). Median and quartiles.
- **Regenerate:** `python figures/make_skill_map.py`

## Presence of patchy genes (hurdle models)

- **Shows:** a, within-cruise presence AUC for 5,808 patchy KOs (20-90% of 1,111 tiles, rarefied to 10^6 genes):
  depth + platform logistic (x) vs Granite PCs logistic (y), coloured by class (satellite 3,047 / technical 1,617 /
  weak 1,144); focus genes circled. b, focus genes (nif, amo, nasA, phn) under four inputs: depth + platform,
  in-situ temperature, OLCI + SST, Granite PCs. c, d, nifH and amoB presence by 10 deg latitude bin: observed share
  vs mean held-out probability. Legend draft in RESULTS.md R65.
- **Regenerate:** fit on Morrill (`cluster/hurdle_fit.sbatch`, 24 shards), copy back
  `data/results/hurdle/fit_shard*.csv.partial`, then `python scripts/hurdle_kos.py --parts summary figure`.

## Why patchy genes are satellite or technical

- **Shows:** a, presence class by eukaryote share of each KO's reference genes; b, within-cruise presence AUC for
  Granite PCs vs the Granite-predicted eukaryote : prokaryote ratio alone, coloured by eukaryote share; c, median
  abundance-where-present R2 by counting-noise quintile, eukaryotic vs prokaryotic KOs x satellite vs technical;
  d, within-cruise presence-depth correlation vs mean reference gene length, by domain, with decile medians.
- **Regenerate:** `python scripts/hurdle_mechanism.py` (needs R65 outputs, R60 C predictions and R44 `ko_gc.csv`);
  `--figure-only` redraws from `mechanism.csv`.

## Patchy genes, base model (supersedes the two hurdle figures above)

- **Shows:** KOs grouped by read-based domain (R67). a, patchy KOs (5,808): within-cruise presence AUC (median, IQR) for
  depth + platform, predicted eukaryote : prokaryote + technical, Granite + technical. b, abundance R2 by domain, core
  KOs (main model, 3,933) and patchy KOs where present, side by side. c, per-KO presence AUC, Granite + technical vs
  predicted C + technical, coloured by eukaryotic read share; focus genes labelled. d, abundance R2 by eukaryotic read
  share, core and patchy. e, patchy KOs: detection (% vs the KO's average) by within-cruise depth fifth, per domain, by
  gene-length tercile. f, core KOs: within-cruise r of abundance with depth, by domain x length. Legend draft in
  RESULTS.md R69 addendum.
- **Regenerate:** fits: `cluster/run_hurdle_base_local.sh` (next time on Morrill), then
  `python scripts/hurdle_base.py --parts summary figure` (figure needs main_model/ko_ceiling.csv and
  targets/clr_all.csv.gz).
