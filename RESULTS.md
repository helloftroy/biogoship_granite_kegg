# Results

Every result is recorded with what it assumed and what went into it, so later
runs can be compared without re-running old ones. Each entry has: **inputs**
(data snapshot and sample set), **settings** (tile rule, SST handling, features,
model, validation), the **result tables**, **files**, and the **command** that
made it. Shared definitions are in section 0; an entry lists only how it differs.

Narrative, methods and earlier history (2026-10-20 to 10-23) are in
[METHODS_RESULTS.md](METHODS_RESULTS.md). Figures are indexed in [FIGURES.md](FIGURES.md).

Results that will change, and why:
- **SST for the samples that still have none** (about 105 per position, see R4):
  moving the window does not recover them; a gap-free SST product (e.g. OSTIA)
  would, and R1 and R3 would then gain samples, mostly on I07N.
- **In-situ temperature as a baseline**: done for 20 KOs (R5) and all KOs (R6).
- **Temperature + ocean colour** (R7, R10): done on 13 cruises.
- **SST variants** (`cluster/sst_variants.sbatch`, then
  `cluster/embed_sst_variants.sbatch`): embeddings with the SST channel blanked,
  and with SST gaps filled from OSTIA; to compare with the current embeddings.
- **CatBoost on a KO shortlist**: R3 uses ridge and nearest-neighbour models only.
- **Model choice** (R8): reductions, other learners and metadata tested on 20 KOs
  with the current embeddings; to be rerun on the SST variants.

---

## 0. Shared definitions

### Data snapshots

| ID | What | Date | Contents |
|---|---|---|---|
| E-orig-1102 | Granite embeddings, original cruises (`data/embeddings/`) | job finished 2026-11-02 09:11 | 3,503 tiles (1,692 fixed, 1,811 advected); SST: 2,642 paired, 305 filled (daytime overpasses, up to 120 h), 556 none |
| E-faire-1103 | Granite embeddings, FAIRe cruises (`data/faire/embeddings/`) | 2026-11-03 | 5,553 tiles (2,811 fixed, 2,742 advected) |
| RAW-1103 | Raw band summaries, original cruises (`data/embeddings/tile_raw_features.csv`) | 2026-11-03 | mean and s.d. of the 17 tile channels (16 OLCI bands + SST) over clear pixels; 3,503 tiles |
| E-orig-1020 | Earlier embeddings, before SST fill (backup in scratchpad only) | 2026-10-20 | 3,503 tiles; 861 without SST |
| KO-1104 | KO counts `ko_counts_cov50_all.tsv` (repo root, not in git) | received 2026-11-04 | 3,081 libraries x 18,835 KOs, coverage-weighted counts |
| PF-1104 | Pfam counts `pfam_counts_cov50_all.tsv` (repo root, not in git) | received 2026-11-04 | 3,084 libraries x 13,161 families (not yet analysed) |
| SSTNEAR-1105 | SST check for moved windows (`data/sst_nearby_check.csv`, log `data/bgs_sst_nearby_83336.out`; not in git) | job 83336 finished 2026-11-05 | 2,195 rows: samples whose best tile has no SST x scene x SST source (paired + 2 fill candidates) |

Embedding model: `ibm-granite/granite-geospatial-ocean` encoder, mean over output
tokens, 512 numbers per tile; missing SST pixels filled with Granite's mean
(normalised 0) before encoding. Tiles are 42 x 42 OLCI pixels (12.6 km), at
least 80% clear, centred on the position or the nearest clear window within 20 km.

### Positions
- **fixed**: tile at the sampling position.
- **advected**: tile at where the sampled water was at the overpass time, from
  back-tracking with daily surface currents (CMEMS MULTIOBS 0.25 deg; near-real-time
  version for P04E 2026).

### Tile rule (which tile represents a sample)
- **closest**: the clear tile closest in time to sampling (scenes within +/- 3 days).
- **prefer-SST**: the closest-in-time tile that has SST; closest overall if none has.
- **SST-only** (`--require-sst`): tiles without SST are never used; a sample whose
  chosen tile has none drops out. Option A = prefer-SST + SST-only; option B =
  closest + SST-only.

### Feature sets
| Name | Numbers | What |
|---|---:|---|
| fixed / advected | 512 | Granite embedding of the chosen tile at that position |
| raw_fixed / raw_advected | 34 | mean and s.d. of each of the 17 channels over the tile's clear pixels |
| fixed+raw / advected+raw | 546 | both side by side |
| metadata ("position + season") | 4-5 | latitude, longitude, day of year as sine and cosine (the seasonal cycle, not the calendar date: with whole cruises held out, an absolute date is never seen in training and cannot be extrapolated) |

### Validation and scores
- **Cruise held out** (default): train on all other cruises, test on every sample
  of one cruise, for each cruise in turn. Each cruise is exactly one fold (fixed
  2026-11-04; before that, fold balancing could merge two cruises).
- **Unique tiles** (default): samples whose chosen fixed tile is identical (same
  pixels: replicates, or neighbours snapped to one clear window) become one row,
  targets averaged.
- **R2** against the training-fold mean. "Mean +/- s.d." is over held-out cruises
  (sample s.d.). **Pooled R2** pools all held-out predictions; it includes
  differences between cruises, so it is higher than the typical within-cruise R2.

---

## R1. In-situ temperature, phosphate, nitrate + nitrite (original cruises)

**Common settings:** targets from the curated sample metadata
(`data/manifests/biogoship_samples.csv`); original 5 Sentinel-3-era cruises
(AMT28, C13.5/A13.5, I07N, I09N, P18); C13.5 has no nutrients, so nutrients use 4
cruises. CatBoost, 600 trees, depth 6, learning rate 0.05, seed 7. Samples must
have a fixed and an advected tile. Script `scripts/compare_feature_sets.py`.

### R1.0 Baseline before SST filling (2026-10-21)
Inputs E-orig-1020; tile rule closest; tiles without SST included (SST filled
with Granite's mean). Rows: temperature 535, phosphate 274, nitrate 177.

| Target | Fixed | Raw fixed | Advected | Raw advected | Position + season |
|---|---:|---:|---:|---:|---:|
| Temperature | 0.71 +/- 0.06 | 0.60 +/- 0.20 | 0.71 +/- 0.08 | 0.61 +/- 0.17 | 0.69 +/- 0.17 |
| Phosphate | 0.46 +/- 0.28 | -0.02 +/- 0.88 | 0.35 +/- 0.40 | 0.03 +/- 0.74 | 0.48 +/- 0.44 |
| Nitrate + nitrite | 0.51 +/- 0.44 | -0.31 +/- 1.16 | 0.47 +/- 0.53 | -0.22 +/- 1.15 | 0.80 +/- 0.20 |

(s.d. here are population s.d., as printed by the script at the time.)

### R1.1 SST filled, closest tile, tiles without SST included (2026-11-04)
Inputs E-orig-1102, RAW-1103. Rows: temperature 535, phosphate 274, nitrate 177.
Files `data/results/sst_filled/`. Command:
`python scripts/compare_feature_sets.py --out-dir data/results/sst_filled`

| Target | Fixed | Raw fixed | Advected | Raw advected | Position + season | Fixed + raw | Advected + raw |
|---|---:|---:|---:|---:|---:|---:|---:|
| Temperature | 0.73 +/- 0.08 | 0.66 +/- 0.20 | 0.73 +/- 0.09 | 0.63 +/- 0.20 | 0.69 +/- 0.17 | 0.71 +/- 0.07 | 0.73 +/- 0.10 |
| Phosphate | 0.47 +/- 0.35 | -0.03 +/- 0.91 | 0.41 +/- 0.44 | -0.15 +/- 0.79 | 0.48 +/- 0.44 | 0.45 +/- 0.38 | 0.43 +/- 0.43 |
| Nitrate + nitrite | 0.57 +/- 0.40 | -0.21 +/- 1.04 | 0.53 +/- 0.50 | -0.47 +/- 1.29 | 0.80 +/- 0.20 | 0.58 +/- 0.38 | 0.58 +/- 0.44 |

Pooled temperature R2 (fixed): samples whose tile has SST 0.78, without SST 0.00.

### R1.2 Prefer-SST tile, tiles without SST still included (2026-11-04)
As R1.1 with `--best-tile prefer-sst`. Rows: 533 / 272 / 177. About 45 samples per
position switch to a tile with SST, at a median 44 h from sampling instead of 16 h.
Files `data/results/prefer_sst/`.

| Target | Fixed | Raw fixed | Advected | Raw advected | Position + season | Fixed + raw | Advected + raw |
|---|---:|---:|---:|---:|---:|---:|---:|
| Temperature | 0.78 +/- 0.05 | 0.63 +/- 0.33 | 0.79 +/- 0.05 | 0.59 +/- 0.36 | 0.68 +/- 0.15 | 0.78 +/- 0.03 | 0.77 +/- 0.06 |
| Phosphate | 0.48 +/- 0.35 | -0.00 +/- 0.87 | 0.39 +/- 0.46 | -0.18 +/- 0.83 | 0.46 +/- 0.52 | 0.49 +/- 0.27 | 0.46 +/- 0.40 |
| Nitrate + nitrite | 0.61 +/- 0.37 | -0.11 +/- 1.05 | 0.49 +/- 0.58 | -0.29 +/- 1.13 | 0.80 +/- 0.20 | 0.64 +/- 0.35 | 0.52 +/- 0.54 |

### R1.3 Option A: prefer-SST, SST-only (2026-11-04) - current standard
As R1.1 with `--best-tile prefer-sst --require-sst`. Rows: 478 / 237 / 156 (493
distinct tiles). Files `data/results/sst_only_prefer/`. Figure: FIGURES.md,
"Predictability", panel a.

| Target | Fixed | Raw fixed | Advected | Raw advected | Position + season | Fixed + raw | Advected + raw |
|---|---:|---:|---:|---:|---:|---:|---:|
| Temperature | **0.87 +/- 0.09** | 0.02 +/- 1.69 | **0.87 +/- 0.07** | 0.04 +/- 1.69 | 0.68 +/- 0.11 | 0.86 +/- 0.07 | 0.87 +/- 0.05 |
| Phosphate | 0.43 +/- 0.45 | -0.53 +/- 1.08 | 0.42 +/- 0.50 | -0.53 +/- 1.12 | 0.54 +/- 0.49 | 0.48 +/- 0.37 | 0.45 +/- 0.47 |
| Nitrate + nitrite | 0.64 +/- 0.39 | -0.48 +/- 1.18 | 0.62 +/- 0.43 | -0.26 +/- 0.97 | 0.78 +/- 0.21 | 0.59 +/- 0.46 | 0.56 +/- 0.53 |

Temperature per held-out cruise:

| Held out | Fixed | Raw fixed | Advected | Raw advected | Position + season |
|---|---:|---:|---:|---:|---:|
| AMT28 | 0.93 | 0.90 | 0.92 | 0.92 | 0.60 |
| C13.5/A13.5 | 0.87 | 0.88 | 0.85 | 0.89 | 0.52 |
| I07N | 0.92 | -2.98 | 0.93 | -2.96 | 0.71 |
| I09N | 0.90 | 0.87 | 0.89 | 0.87 | 0.79 |
| P18 | 0.71 | 0.46 | 0.77 | 0.49 | 0.77 |

### R1.4 Option B: closest, SST-only (2026-11-04)
As R1.1 with `--best-tile closest --require-sst`. Rows: 441 / 222 / 149 (455
distinct tiles). Files `data/results/sst_only_closest/`.

| Target | Fixed | Raw fixed | Advected | Raw advected | Position + season | Fixed + raw | Advected + raw |
|---|---:|---:|---:|---:|---:|---:|---:|
| Temperature | **0.87 +/- 0.07** | 0.07 +/- 1.60 | 0.85 +/- 0.07 | -0.01 +/- 1.75 | 0.66 +/- 0.11 | 0.86 +/- 0.08 | 0.85 +/- 0.06 |
| Phosphate | 0.43 +/- 0.48 | -0.62 +/- 1.21 | 0.42 +/- 0.48 | -0.60 +/- 1.19 | 0.59 +/- 0.42 | 0.45 +/- 0.44 | 0.47 +/- 0.39 |
| Nitrate + nitrite | 0.59 +/- 0.45 | -0.64 +/- 1.36 | 0.59 +/- 0.44 | -0.53 +/- 1.26 | 0.80 +/- 0.19 | 0.59 +/- 0.45 | 0.59 +/- 0.46 |

Temperature per held-out cruise: AMT28 0.93, C13.5 0.88, I07N 0.93, I09N 0.84, P18
0.75 (fixed); raw fixed I07N -2.77.

### Reading R1
- Missing SST was the main limit on temperature: dropping tiles without SST
  takes the embeddings from 0.73 to 0.87, well above position + season (0.66-0.68).
  Options A and B agree; A keeps 37 more rows.
- Nutrients: embeddings 0.43-0.64, position + season 0.54-0.80; differences sit
  inside the cruise-to-cruise spread.
- Raw band averages match the embeddings on AMT28, C13.5 and I09N but collapse on
  I07N (-2.98) and are weak on P18. Not the SST channel (R2 below): the other 33
  band summaries push I07N predictions out of range (predicted s.d. 6.4 C against
  2.3 C observed). The embedding handles I07N (0.92).
- Fixed and advected remain indistinguishable.

## R2. Satellite SST against in-situ temperature (2026-11-04)

Inputs E-orig-1102, RAW-1103; best fixed tile (closest); tile SST = `SST_mean`
de-normalised (mean 293.908 K, s.d. 2.676 K). Ad hoc check, not a script.

| Cruise | Tiles with SST | Median SST - in situ (C) | Correlation |
|---|---:|---:|---:|
| AMT28 | 46 | -0.13 | 0.99 |
| C13.5/A13.5 | 251 | -0.21 | 0.99 |
| I07N | 194 | -0.05 | 0.92 |
| I09N | 56 | +0.16 | 0.44 (temperature range under 3 C) |
| P18 | 67 | +0.02 | 0.99 |

Filled SST (other overpass) reads warm: median +0.81 C on I07N, +1.21 C on P18,
+0.28 C on C13.5, against -0.06 to -0.22 C for paired SST.

## R3. KOs predicted from Granite embeddings (2026-11-04)

**Targets** (`scripts/build_ko_targets.py`, inputs KO-1104 and the FAIRe workbook):
- Library to sample: original-cruise library names (e.g. `C132020_021v2_meta`)
  by name to `biogoship_samples.csv`; SRR accessions via `metagenome_runs`; all
  others through the workbook's experimentRunMetadata sheet (3,048 of 3,084
  libraries matched; 35 ERR + 1 SRR unplaced). The FAIRe `original_sample_id`
  link is by time and place and can pick a neighbouring replicate, so names win.
- Dropped: 129 libraries without a sample, controls (blanks, mocks, references)
  or another project (`NP_*`); 58 libraries under 100,000 counts; 859 libraries
  deeper than 10 m (no recorded depth = kept).
- Proportions per library; libraries of one sample averaged (28 samples).
- KOs present in at least 90% of samples: 3,921 of 18,835. CLR with half a count
  as pseudocount.
- 2,007 surface samples from 16 cruises (`data/targets/ko_clr.csv`,
  `ko_samples.csv`). No P04E 2026 metagenomes yet.
- Sequencing: cruises sampled 2021-2024 were sequenced in 2025-26 on instrument
  OBI_3310 at 2-4x depth; older cruises on OBI_2631/2131. Batch and cruise are
  confounded.

**Settings** (`scripts/screen_ko_predictability.py`): inputs E-orig-1102 +
E-faire-1103; fixed position; option A (prefer-SST, SST-only); unique tiles;
13 cruises held out in turn (AMT28, C13.5/A13.5, I07N, I09N, P18, A13.5-2024,
A16N-2023, GOMECC4, HLY2202-2022, I05-2023, I08S-2024, I09N-2025, P02-2022).
1,339 samples, 1,127 distinct tiles. Two learners per feature set, the better one
per KO reported: ridge (alpha per KO by leave-one-out on training cruises) and
15-nearest-neighbour (distance-weighted); features standardised on training
cruises. Feature sets: embedding (512), position + season (5: latitude/90,
sin/cos longitude, sin/cos day of year), embedding shuffled between samples
(chance). Raw band summaries not included (none yet for FAIRe tiles).

### R3.1 All 13 cruises

| | KOs with pooled R2 >= 0.1 | >= 0.25 | >= 0.5 | Embedding > position + season |
|---|---:|---:|---:|---:|
| Granite embedding | 1,748 | 693 | 69 | 3,224 of 3,921 |
| Position + season | 225 | - | - | |
| Shuffled embedding | 0 (99th percentile -0.015) | | | |

Ridge was the better learner for 2,817 KOs, nearest-neighbour for 1,104.
Files `data/results/ko_screen/ko_screen.csv` (per KO), `ko_screen_per_cruise.csv`.
Command: `python scripts/screen_ko_predictability.py`

### R3.2 Within one sequencing batch
Same settings, targets restricted by instrument (`--instrument`).

| Batch | Cruises | Samples | KOs R2 >= 0.1 | >= 0.25 | >= 0.5 | Embedding > position + season | KO ranking vs R3.1 (Spearman) |
|---|---:|---:|---:|---:|---:|---:|---:|
| OBI_2631 (older) | 6 | 549 | 1,593 | 643 | 103 | 84% | 0.74 |
| OBI_3310 (newer) | 7 | 790 | 1,337 | 640 | 19 | 69% | 0.73 |

Files `data/results/ko_screen_OBI_2631/`, `ko_screen_OBI_3310/`.

### Top 20 KOs (R3.1), with the within-batch scores

| KO | Pooled R2 | Median per cruise | Position + season | Within OBI_2631 | Within OBI_3310 |
|---|---:|---:|---:|---:|---:|
| hoxH · NAD-reducing hydrogenase large subunit, K00436 | 0.70 | 0.36 | -0.26 | 0.74 | 0.27 |
| ndhK · NAD(P)H-quinone oxidoreductase subunit K, K05582 | 0.69 | 0.36 | -0.24 | 0.70 | 0.26 |
| ndhJ · NAD(P)H-quinone oxidoreductase subunit J, K05581 | 0.68 | 0.31 | -0.23 | 0.75 | 0.22 |
| mntB · manganese transport permease, K11602 | 0.68 | 0.54 | -0.09 | 0.78 | 0.19 |
| ndhD · NAD(P)H-quinone oxidoreductase subunit 4, K05575 | 0.68 | 0.28 | -0.23 | 0.71 | 0.29 |
| fatty aldehyde-generating acyl-ACP reductase, K14330 | 0.65 | 0.45 | -0.17 | 0.75 | 0.21 |
| cpeS · phycoerythrin-associated linker protein, K05382 | 0.65 | 0.43 | -0.19 | 0.74 | 0.18 |
| mntC · manganese transport substrate-binding protein, K11601 | 0.64 | 0.37 | -0.13 | 0.76 | 0.08 |
| K07015 (no KEGG name) | 0.64 | 0.32 | -0.23 | 0.68 | 0.24 |
| nrdJ · ribonucleotide reductase class II, K00524 | 0.64 | 0.39 | -0.19 | 0.74 | 0.15 |
| aroH · chorismate mutase, K06208 | 0.64 | 0.59 | -0.13 | 0.59 | 0.30 |
| cpeB · phycoerythrin beta chain, K05377 | 0.64 | 0.40 | -0.12 | 0.73 | 0.05 |
| 1,4-dihydroxy-2-naphthoyl-CoA hydrolase, K12073 | 0.63 | 0.22 | -0.27 | 0.71 | 0.13 |
| ecfT · energy-coupling factor transporter permease, K16785 | 0.63 | 0.47 | -0.10 | 0.72 | 0.06 |
| stpA · glucosylglycerol 3-phosphatase, K05978 | 0.62 | 0.21 | -0.27 | 0.70 | 0.16 |
| gpgS · glucosyl-3-phosphoglycerate synthase, K13693 | 0.62 | 0.49 | -0.13 | 0.46 | 0.59 |
| pebA · 15,16-dihydrobiliverdin:ferredoxin oxidoreductase, K05369 | 0.62 | 0.32 | -0.13 | 0.79 | 0.06 |
| fatty aldehyde decarbonylase, K14331 | 0.62 | 0.17 | -0.24 | 0.69 | 0.12 |
| tyrAa · arogenate dehydrogenase, K15226 | 0.62 | 0.29 | -0.20 | 0.71 | 0.07 |
| gap2 · glyceraldehyde-3-phosphate dehydrogenase (NADP), K00150 | 0.62 | 0.13 | -0.18 | 0.72 | 0.15 |

Also: pstS (K02040, phosphate uptake) 0.38 vs 0.13; pufL (K08928, aerobic
anoxygenic phototrophs) 0.25 vs -0.37. Classic markers nifH, nifD, proteorhodopsin,
amoA, narG and phnX are not screened: present in under 90% of samples.

### Reading R3
- About 1,750 KOs are predictable from the embedding in new cruises at pooled
  R2 >= 0.1, against about 225 from position + season, and none by chance.
- The best KOs are cyanobacterial, mostly Synechococcus (phycoerythrin, NDH-1,
  glucosylglycerol, alkane biosynthesis, ntcA); plausible, since phycoerythrin and
  chlorophyll colour the water.
- Not a sequencing-batch shortcut overall: within each batch the number of
  predictable KOs holds (about 640 at R2 >= 0.25 in each) and the ranking matches
  (0.73-0.74). But the top cyanobacterial KOs are much weaker within the newer
  batch (0.05-0.30 vs 0.7-0.8 within the older one). Unresolved: the newer batch's
  cruises may span less Synechococcus variation, or the newer sequencing may record
  it differently. To check: per-cruise ranges of these KOs in each batch.
- Pooled R2 includes between-cruise differences; median R2 within a cruise is
  0.13-0.59 for the top 20.
- Position + season does not transfer to unseen cruises for KOs (median pooled R2
  below 0), unlike for temperature and nutrients.

## R4. Does moving the tile recover SST? (2026-11-05)

Inputs SSTNEAR-1105, E-orig-1102 tile manifest. Samples: original cruises whose
best tile (closest) has no SST: 103 fixed, 105 advected (site level). For every
scene of the sample (tile_log valid, cloudy or edge) the check searched window
positions up to 40 km from the current one (3-pixel steps) that are at least 80%
clear in OLCI **and** have SST of quality 4 or higher (the level the pipeline uses)
anywhere in the window's SST box (as in `SstSwath.resample`). SST sources: the
paired SLSTR overpass and up to 2 daytime fill candidates. Compared with option A
(prefer-SST: same position, tile from another scene that has SST).
`python scripts/check_sst_nearby.py` (sbatch `cluster/check_sst_nearby.sbatch`);
comparison against option A ad hoc.

Samples recovered (counts; gaps are median hours from sampling to the scene used):

| Position | Cruise | No SST | Option A (other scene) | Move <= 10 km | <= 20 km | <= 40 km | Only by moving <= 20 km | Option A or move <= 20 km | Gap option A (h) | Gap move <= 20 km (h) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| fixed | AMT28 | 11 | 3 | 3 | 4 | 7 | 1 | 4 | 28 | 36 |
| fixed | C13.5/A13.5 | 39 | 21 | 18 | 23 | 32 | 3 | 24 | 43 | 28 |
| fixed | I07N | 52 | 20 | 10 | 19 | 33 | 1 | 21 | 60 | 57 |
| fixed | I09N | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 177 | 177 |
| fixed | all | 103 | 45 | 32 | 47 | 73 | 5 | 50 | 46 | 38 |
| advected | AMT28 | 11 | 2 | 2 | 3 | 8 | 1 | 3 | 47 | 3 |
| advected | C13.5/A13.5 | 34 | 20 | 14 | 20 | 27 | 2 | 22 | 45 | 32 |
| advected | I07N | 60 | 20 | 20 | 23 | 36 | 5 | 25 | 51 | 57 |
| advected | all | 105 | 42 | 36 | 46 | 71 | 8 | 50 | 51 | 44 |

Other numbers:
- In the closest scene itself (the one that lacks SST), a move of <= 20 km finds
  SST for only 25 of 208; quality 3 SST at the current window would add 17.
- Fill candidates add nothing: candidate 1 is the paired overpass in 441 of 455
  rows, and candidate 2 (about 2 h later) has SST near the sample once.
- Moving gives a scene closer in time than option A for 10 fixed and 9 advected
  samples.

### Reading R4

- SST gaps are not near misses: where SLSTR has no good SST at the sample, it
  usually has none for tens of km (cloud masking by SLSTR is stricter than the OLCI
  clear test). Moving <= 20 km recovers about the same samples as option A, adding
  5 fixed and 8 advected; 40 km recovers more but the tile is then 3 tile-widths
  from the sample.
- Not worth building into `download_tiles`. Option A stays the standard; about
  half of these samples (53 fixed, 55 advected) have no SST by either route and
  drop out under SST-only. Recovering them needs a gap-free SST product.

## R5. Position + season + in-situ temperature against Granite, 20 KOs (2026-11-05)

Inputs E-orig-1102, E-faire-1103, KO-1104 (targets as in R3). KOs from the R3
screen: the 10 best by embedding pooled R2, and 10 middling ones (embedding pooled
R2 0.25-0.40, at least 0.15 above position + season), evenly spaced in rank.
Samples as in R3 (fixed, prefer-SST, SST-only, identical tiles merged) **and with
an in-situ temperature** (manifest `temperature_c`): 1,320 samples, 1,110 tiles,
13 cruises; every feature set sees the same rows. Leave-one-cruise-out.
Feature sets: metadata (position + season, as in R3); metadata + temperature;
temperature alone; embedding; embedding + temperature. Learners: ridge, kNN
(k = 15), CatBoost (500 trees, depth 6, lr 0.05, MultiRMSE); score = best learner
per KO and feature set. Files `data/results/ko_temperature/`. Command:
`python scripts/test_ko_temperature.py` (67 min locally; CatBoost on the
embedding is most of it).

Median over the 10 KOs of each group (pooled R2; in brackets median R2 within a
held-out cruise):

| Group | Position + season | + temperature | Temperature alone | Embedding | Embedding + temperature |
|---|---:|---:|---:|---:|---:|
| Top 10 | 0.31 (0.29) | 0.34 (0.24) | 0.64 (0.37) | 0.67 (0.36) | 0.67 (0.35) |
| Middle 10 | -0.01 (0.03) | 0.19 (0.07) | 0.15 (0.24) | 0.33 (0.30) | 0.33 (0.31) |

Per KO (pooled R2):

| KO | Gene | Group | Pos + season | + temp | Temp | Embedding | Emb + temp |
|---|---|---|---:|---:|---:|---:|---:|
| K00436 | hoxH | top | 0.29 | 0.35 | 0.70 | 0.71 | 0.71 |
| K05582 | ndhK | top | 0.32 | 0.34 | 0.65 | 0.69 | 0.69 |
| K05581 | ndhJ | top | 0.31 | 0.34 | 0.70 | 0.70 | 0.70 |
| K11602 | mntB | top | 0.43 | 0.44 | 0.61 | 0.69 | 0.68 |
| K05575 | ndhD | top | 0.29 | 0.32 | 0.71 | 0.69 | 0.69 |
| K14330 | acyl-ACP reductase | top | 0.32 | 0.34 | 0.67 | 0.66 | 0.66 |
| K05382 | cpeS | top | 0.32 | 0.35 | 0.62 | 0.65 | 0.65 |
| K11601 | mntC | top | 0.40 | 0.40 | 0.57 | 0.65 | 0.65 |
| K07015 | (no KEGG name) | top | 0.25 | 0.28 | 0.63 | 0.65 | 0.65 |
| K00524 | nrdJ | top | 0.30 | 0.31 | 0.63 | 0.65 | 0.65 |
| K21147 | moeZR | middle | 0.28 | 0.24 | 0.39 | 0.41 | 0.41 |
| K03523 | bioY | middle | 0.04 | 0.18 | 0.15 | 0.38 | 0.37 |
| K06181 | rluE | middle | -0.04 | 0.11 | 0.28 | 0.36 | 0.36 |
| K20222 | importin-5 | middle | -0.09 | 0.22 | 0.26 | 0.34 | 0.34 |
| K21929 | udg | middle | 0.03 | 0.26 | 0.13 | 0.32 | 0.32 |
| K01092 | suhB | middle | 0.10 | 0.17 | 0.10 | 0.30 | 0.30 |
| K05380 | cpeE | middle | -0.16 | -0.07 | 0.15 | 0.34 | 0.37 |
| K03526 | ispG | middle | -0.01 | 0.25 | 0.08 | 0.28 | 0.28 |
| K01662 | dxs | middle | -0.19 | 0.03 | 0.06 | 0.27 | 0.27 |
| K22125 | myoferlin | middle | -0.00 | 0.21 | 0.23 | 0.26 | 0.26 |

All values in both tables are the **best of the three learners** per KO and
feature set, so the columns mix learners: CatBoost wins for position + season
(ridge and kNN are below 0 there) and for temperature alone on the top KOs; ridge
wins for the embedding (9 of 10 in each group). (R6 later found spline + ridge a
little better than CatBoost on temperature alone: median 0.53 vs 0.48 over these
20 KOs.)

### Reading R5

- For the top (cyanobacterial) KOs, in-situ temperature alone is as good as the
  embedding (0.64 vs 0.67 pooled, 0.37 vs 0.36 within cruise): these KOs are
  mostly a temperature signal.
- For the middling KOs the embedding stays ahead (0.33 vs 0.15-0.19 pooled, 0.30
  vs 0.24 within cruise): it carries something beyond temperature.
- Adding position + season to temperature makes the top KOs worse (0.64 -> 0.34):
  the models learn cruise-specific geography that does not carry to a new cruise.
  Temperature alone transfers better.
- Adding temperature to the embedding changes nothing: the embedding already
  holds the temperature information (SST channel and/or ocean colour).
- Open: whether the embedding without its SST channel keeps this (that decides
  whether SST gaps matter); and temperature on all KOs, not just 20.

## R6. In-situ temperature alone against Granite, all KOs (2026-11-05)

As R5 (same 1,320 samples, 1,110 tiles, 13 cruises; leave-one-cruise-out), but
all 3,921 KOs and two feature sets: in-situ temperature alone, and the embedding.
Learners without CatBoost: temperature with ridge, spline + ridge (cubic
B-spline, 8 knots) and kNN; embedding with ridge and kNN. Best learner per KO and
feature set. Spline + ridge is a fair stand-in for CatBoost on temperature: on
R5's 20 KOs it scores a median 0.53 against CatBoost's 0.48 (kNN 0.41, ridge 0.22).
Files `data/results/ko_temperature_all/`. Command:
`python scripts/test_ko_temperature.py --kos all --learners ridge spline knn --feature-sets temp embedding --out-dir data/results/ko_temperature_all`
(a few minutes). Figure: FIGURES.md, "Embedding vs in-situ temperature, all KOs".

"Clearly better": pooled R2 at least 0.1 higher **and** better in at least 9 of
the 13 held-out cruises.

| | KOs |
|---|---:|
| All | 3,921 |
| Pooled R2 >= 0.1: embedding / temperature | 1,772 / 1,246 |
| Embedding ahead by >= 0.1 (pooled) | 1,067 |
| Temperature ahead by >= 0.1 (pooled) | 469 |
| Embedding clearly better (also >= 9 of 13 cruises) | 299 (111 of them with embedding R2 >= 0.25) |
| Temperature clearly better | 366 (38 with temperature R2 >= 0.25) |
| Neither clearly | 3,256 |

Median difference (embedding - temperature) +0.02; 10th to 90th percentile -0.11
to +0.20. Among KOs with embedding R2 >= 0.25 (722) the median difference is +0.06.

Embedding clearly better, highest embedding R2 (temperature R2 -> embedding R2):
mannosyl-3-phosphoglycerate phosphatase K07026 (0.32 -> 0.59), aguB K12251
(0.18 -> 0.51), pct K01026 (0.39 -> 0.50), ureB K01429 (0.22 -> 0.50), aliphatic
nitrilase K01502, sucrose phosphorylase K00690, glsA K01425 (0.16 -> 0.47), minC,
purL K01952 (0.17 -> 0.43, better in all 13 cruises), pheA2, citramalyl-CoA lyase.
Temperature clearly better: gap2/gapB K00150 (0.72 vs 0.62), ndhF K05577 (0.64 vs
0.50), fiu K16090, AGPAT, NOP2, TOP2, citT, and the purple-bacterial reaction
centre genes puhA, pufL; most of the rest have temperature R2 near 0 and the
embedding below 0.

### Reading R6

- Across all KOs the embedding adds little beyond temperature on average (median
  +0.02); for most KOs whatever is predictable is mostly temperature.
- It does add for a sizeable minority: about 1,070 KOs are 0.1 or more ahead, and
  299 hold up in at least 9 of 13 cruises. Their median is temperature R2 -0.01 ->
  embedding 0.20: KOs temperature cannot predict at all.
- Where temperature clearly wins, the embedding is usually worse than the mean
  (below 0): it overfits those KOs rather than missing a signal.
- Satellite SST is in the embedding; whether that is what carries the temperature
  signal is the SST-variant test (blanked vs OSTIA-filled).

## R7. Temperature + ocean colour, 20 KOs (2026-11-05; rerun on 13 cruises 2026-11-06)

Question: is the embedding's "beyond temperature" just chlorophyll? Feature sets:
temperature; band means (16 OLCI bands, mean over clear pixels, RAW-1103);
temperature + band means; temperature + a chlorophyll index
log10(max(Rrs443, Rrs490, Rrs510)/Rrs560) (the OC4ME band ratio); embedding.
Learners ridge, kNN, CatBoost (not on the embedding). R5's 20 KOs. **Only samples
with band means: the 5 original cruises** (517 samples, 471 tiles), since the FAIRe
tiles are not pooled yet; so 5 folds, and the scores are not comparable to R5.
Files `data/results/ko_colour_orig5/`. Command:
`python scripts/test_ko_temperature.py --feature-sets temp colour temp_colour temp_chl embedding --out-dir data/results/ko_colour_orig5`
(run before the spline learner was added).

Median over the 10 KOs (pooled R2; in brackets median within a held-out cruise):

| Group | Temperature | Band means | Temp + band means | Temp + chl index | Embedding |
|---|---:|---:|---:|---:|---:|
| Top 10 | 0.37 (0.27) | -0.15 (-0.12) | 0.14 (0.21) | 0.37 (0.33) | 0.26 (0.16) |
| Middle 10 | 0.03 (0.12) | 0.05 (0.12) | 0.09 (0.20) | 0.05 (0.10) | 0.17 (0.19) |

Reading (provisional, 5 cruises): band means on their own do not transfer between
cruises, and adding all 16 to temperature hurts (too many inputs for 5 cruises);
the chlorophyll index adds little to temperature. On the middling KOs the
embedding stays ahead of temperature + colour. To be rerun on 13 cruises.

**Rerun on all 13 cruises (2026-11-06)**, after `cluster/sst_variants.sbatch`
pooled the FAIRe tiles (`data/faire/embeddings/tile_raw_features.csv`, 5,553
tiles). It used the same feature sets, with learners ridge, spline, kNN and
CatBoost, best per KO. Samples: 1,320 / 1,110 tiles / 13 cruises (every R5
sample's tile has band means), so the scores are comparable to R5. Files
`data/results/ko_colour/`. Command:
`python scripts/test_ko_temperature.py --feature-sets temp colour temp_colour temp_chl embedding --out-dir data/results/ko_colour`

| Group | Temperature | Band means | Temp + band means | Temp + chl index | Embedding |
|---|---:|---:|---:|---:|---:|
| Top 10 | 0.67 (0.37) | 0.05 (0.19) | 0.60 (0.37) | 0.69 (0.38) | 0.67 (0.35) |
| Middle 10 | 0.16 (0.22) | 0.20 (0.25) | 0.24 (0.28) | 0.24 (0.26) | 0.31 (0.30) |

Winning learners: spline for temperature and temperature + chl (15 and 16 of 20
KOs), kNN or CatBoost for the band means, ridge for the embedding. Ridge and the
spline fail on the band means pooled (median R2 far below 0), because some FAIRe
tiles have extreme band means (up to 50 s.d. in Oa12-Oa21, probably glint, ice or
cloud edges). Only the tree and neighbour learners cope with them.

Reading:
- On the top KOs, temperature already does everything. Colour adds nothing
  (temperature + chl 0.69 against 0.67).
- On the middling KOs, simple colour closes about half the gap between
  temperature and the embedding pooled (0.16 to 0.24, against 0.31), and most of
  it within cruises (0.22 to 0.26-0.28, against 0.30).
- **So part of "beyond temperature" is chlorophyll / ocean colour, but not all
  of it.** The embedding still leads on 8 of the 10 middling KOs.


## R8. Reductions, other learners and metadata with the embedding, 20 KOs (2026-11-06)

Question: does a different handling of the 512 embedding numbers (PCA or PLS,
then kernel ridge, a small MLP or CatBoost), or adding metadata (position,
season, in-situ temperature, salinity, and in a subset nitrate + nitrite and
phosphate) to the embedding, make a better model than ridge on the embedding?

Inputs: E-orig-1102 + E-faire-1103, KO targets as in R5/R6, option A, fixed
position. Sample sets: **all** = samples with in-situ temperature and salinity,
1,281 samples / 1,093 tiles / 13 cruises; **np** = also measured nitrate + nitrite
and phosphate (manifest, filled from `faire_samples.csv`), 456 samples / 331
tiles / 12 cruises. Settings: one fold per held-out cruise; PCA and PLS fitted on
the standardized embedding of the training cruises only (k = 16, 32, 64 for PCA;
32 for PLS). Kernel ridge: RBF, alpha and gamma picked by an inner split of the
training cruises into 3 groups. MLP: one hidden layer of 64, weight decay picked
the same way, early stopping. CatBoost: 500 trees, depth 6, rate 0.05, one
multi-output model. "meta" = position + season features + temperature +
salinity. KOs: 5 per category from R6's `ko_temperature_gap.csv`, spread in rank:
**clear** (embedding clearly better, embedding R2 >= 0.2), **ahead** (embedding
ahead by >= 0.1 pooled but not on 9 cruises, embedding R2 >= 0.2), **similar**
(|difference| < 0.03, both >= 0.3), **worse** (temperature clearly better,
temperature R2 >= 0.2).
Files `data/results/model_experiments/` (scores.csv, per_cruise.csv). Command:
`python scripts/model_experiments.py` (about 2.5 h on the laptop, mostly kernel ridge).

Median over the KOs, pooled R2 (in brackets median R2 within a held-out cruise):

| Features | Learner | Clear | Ahead | Similar | Worse | All 20 |
|---|---|---:|---:|---:|---:|---:|
| Embedding (512) | ridge | 0.25 (0.24) | 0.29 (0.37) | 0.37 (0.36) | 0.06 (0.07) | 0.29 (0.29) |
| PCA 16 | ridge | 0.29 (0.26) | 0.29 (0.37) | 0.40 (0.33) | 0.18 (0.02) | 0.31 (0.24) |
| PCA 32 | ridge | 0.27 (0.24) | 0.29 (0.37) | 0.40 (0.29) | 0.14 (-0.01) | **0.32** (0.28) |
| PCA 64 | ridge | 0.29 (0.25) | 0.24 (0.39) | 0.34 (0.32) | 0.10 (0.04) | 0.28 (0.26) |
| PCA 16 | kernel ridge | 0.17 (0.24) | 0.12 (0.19) | 0.36 (0.30) | 0.14 (0.17) | 0.22 (0.22) |
| PCA 32 | kernel ridge | 0.21 (0.28) | 0.07 (0.05) | 0.33 (0.33) | 0.10 (0.18) | 0.18 (0.23) |
| PCA 64 | kernel ridge | 0.20 (0.23) | -0.02 (-0.15) | 0.33 (0.21) | 0.05 (0.17) | 0.20 (0.18) |
| PCA 16 | MLP | 0.15 (0.22) | 0.04 (0.04) | 0.37 (0.28) | 0.12 (0.03) | 0.23 (0.16) |
| PCA 32 | MLP | 0.18 (0.23) | 0.10 (-0.02) | 0.36 (0.36) | 0.11 (0.07) | 0.20 (0.21) |
| PCA 64 | MLP | 0.08 (0.09) | -0.89 (-0.01) | 0.32 (0.25) | 0.04 (0.10) | 0.03 (0.12) |
| PCA 32 | CatBoost | 0.28 (0.26) | 0.21 (0.31) | 0.35 (0.31) | 0.17 (0.23) | 0.27 (**0.31**) |
| PLS 32 | ridge | 0.27 (0.21) | 0.26 (0.32) | 0.37 (0.29) | 0.07 (0.02) | 0.26 (0.24) |
| PLS 32 | kernel ridge | 0.16 (0.19) | 0.07 (0.08) | 0.35 (0.28) | 0.10 (0.16) | 0.16 (0.19) |
| PLS 32 | CatBoost | 0.28 (0.20) | 0.20 (0.31) | 0.35 (0.31) | 0.14 (0.21) | 0.26 (0.30) |
| Temperature | spline | 0.15 (0.11) | 0.10 (0.09) | 0.35 (0.33) | 0.22 (0.24) | 0.22 (0.25) |
| Meta | ridge | 0.06 (0.09) | -0.10 (-0.03) | 0.31 (0.26) | **0.25** (0.22) | 0.22 (0.16) |
| Meta | kernel ridge | -0.11 (0.08) | -0.35 (0.08) | 0.26 (0.27) | 0.12 (0.15) | 0.07 (0.15) |
| Meta | CatBoost | 0.16 (0.25) | -0.02 (0.02) | 0.29 (0.28) | 0.03 (0.09) | 0.11 (0.19) |
| Temp + PCA 4 | kernel ridge | 0.18 (0.18) | 0.13 (0.32) | 0.37 (0.32) | 0.21 (0.22) | 0.25 (0.26) |
| Temp + PCA 8 | kernel ridge | 0.21 (0.22) | 0.17 (0.26) | 0.39 (0.31) | 0.19 (0.23) | 0.25 (0.26) |
| Temp + PCA 8 | CatBoost | 0.30 (0.22) | 0.24 (0.26) | 0.35 (0.32) | 0.17 (0.16) | 0.26 (0.26) |
| Embedding + meta | ridge | 0.19 (0.22) | 0.19 (0.21) | 0.35 (0.30) | 0.10 (0.11) | 0.21 (0.21) |
| PCA 32 + meta | ridge | 0.32 (0.24) | 0.15 (0.19) | 0.37 (0.21) | 0.11 (0.12) | 0.26 (0.20) |
| PCA 32 + meta | kernel ridge | 0.16 (0.14) | 0.09 (-0.16) | 0.37 (0.30) | 0.06 (0.09) | 0.16 (0.18) |
| PCA 32 + meta | MLP | 0.01 (0.15) | -0.18 (-0.03) | 0.38 (0.29) | -0.03 (0.11) | 0.11 (0.17) |
| PCA 32 + meta | CatBoost | 0.25 (0.28) | 0.07 (0.13) | 0.35 (0.31) | 0.10 (0.18) | 0.23 (0.22) |
| PLS 32 + meta | kernel ridge | 0.14 (0.12) | 0.08 (0.02) | 0.37 (0.29) | 0.04 (0.12) | 0.15 (0.16) |
| PLS 32 + meta | CatBoost | 0.26 (0.31) | 0.09 (0.15) | 0.37 (0.31) | 0.12 (0.26) | 0.24 (0.25) |

Sample set np (456 samples, 12 cruises):

| Features | Learner | Clear | Ahead | Similar | Worse | All 20 |
|---|---|---:|---:|---:|---:|---:|
| Embedding | ridge | 0.24 (0.31) | 0.21 (0.23) | 0.42 (0.40) | 0.18 (0.11) | 0.27 (0.27) |
| Temperature | spline | 0.16 (0.16) | 0.14 (0.13) | 0.45 (0.38) | 0.27 (0.18) | 0.23 (0.19) |
| Meta | CatBoost | 0.06 (0.04) | -0.07 (0.09) | 0.30 (0.30) | -0.03 (0.09) | 0.12 (0.11) |
| Meta + N, P | ridge | -0.19 (-0.24) | -0.27 (-0.14) | 0.41 (0.39) | -0.00 (-0.13) | -0.02 (-0.02) |
| Meta + N, P | kernel ridge | 0.05 (0.00) | -0.48 (-0.15) | 0.32 (0.40) | 0.11 (-0.04) | 0.09 (0.14) |
| Meta + N, P | CatBoost | 0.19 (0.22) | 0.02 (0.07) | 0.33 (0.34) | 0.00 (0.09) | 0.14 (0.16) |
| Embedding + meta + N, P | ridge | 0.26 (0.31) | 0.06 (0.23) | 0.43 (0.38) | 0.13 (0.13) | 0.29 (0.30) |
| PCA 32 + meta + N, P | kernel ridge | 0.24 (0.24) | -0.49 (-0.03) | 0.29 (0.33) | 0.13 (-0.16) | 0.19 (0.18) |
| PCA 32 + meta + N, P | CatBoost | 0.28 (0.32) | 0.04 (0.21) | 0.40 (0.39) | 0.16 (0.24) | 0.26 (**0.33**) |

Per KO against embedding + ridge (set all): PCA 32 + ridge is ahead on 13 of 20
KOs, median gain +0.006; PCA 16 + ridge 12 of 20, +0.003; PCA 32 + CatBoost 8 of
20, -0.015; embedding + meta ridge 7 of 20, -0.016; PCA 32 + meta ridge 8 of 20,
-0.030.

Reading:
- **No model is clearly better than ridge on the embedding.** PCA 16 or 32 +
  ridge ties it (the median gain per KO is under 0.01). PCA 32 + CatBoost is a
  little better within cruises (0.31 vs 0.29) and a little worse pooled.
- **Kernel ridge and the MLP are worse**, pooled and mostly within cruises, at
  every k; 64 PCs + MLP fails. With about 1,100 tiles and 12 training cruises the
  extra flexibility overfits.
- **PLS is no better than PCA** at k = 32.
- **Adding metadata to the embedding lowers R2** (0.29 to 0.21 with ridge; 0.32 to
  0.26 with PCA 32). Position and season do not carry over between cruises, as
  in R5. The "worse" KOs are the only place a model with temperature helps.
  There, temperature alone, or meta + ridge (0.25), is still best, and no
  embedding model reaches it.
- **N and P add little.** Meta + N, P is weak on its own (at best 0.14).
  Embedding + meta + N, P (ridge, 0.29) and PCA 32 + meta + N, P (CatBoost,
  within-cruise 0.33) are slightly ahead of the embedding on the same samples
  (0.27). This is on 456 samples, and it loses on the "ahead" KOs.
- The categories prefer different models: temperature-led for "worse", embedding
  ridge for "ahead". So a per-KO choice or an ensemble could gain, but it was not
  tested.

Not tried yet:
- Temperature + a chlorophyll index + PCs. This needs FAIRe band means from
  `cluster/sst_variants.sbatch`.
- PLS at k = 16 or 64. Reduced-rank regression or other supervised reductions.
- Kernel ridge on the full embedding, or other kernels.
- Deeper or wider MLPs, an MLP on the full 512 numbers, Gaussian processes.
- Tuned CatBoost (default settings only).
- Per-KO hyperparameters for kernel ridge. Here alpha and gamma are shared across
  the KOs.
- Stacking or averaging models, or a per-KO choice of model (needs a nested
  choice to stay fair).
- PCA fitted on all tiles of both sample sets. This leaks no targets, but these
  runs used a fold-fitted PCA.
- The advected position; the in-situ targets; all 3,921 KOs.
- **Everything here to be rerun on the blank-SST and OSTIA-filled embeddings**
  (the script takes `--embeddings-dirs`).

## R9. Which metadata helps the embedding? Ridge, 20 KOs (2026-11-06)

As R8 (sample set all: 1,281 samples, 1,093 tiles, 13 cruises; the same 20 KOs;
PCA fitted on the training cruises only), ridge only. Metadata blocks: **temp** =
in-situ temperature; **ts** = temperature + salinity; **season** = sine and cosine
of day of year; **meta** = ts + season + position (lat, sine and cosine of
longitude). "ts + season" is the R8 metadata without latitude and longitude, and
"ts" is without position and time. Each block is standardized and shares one
ridge penalty with the embedding inputs.
Files `data/results/model_experiments_meta/`. Command:
`python scripts/model_experiments.py --only all:emb+ts:ridge ... --out-dir data/results/model_experiments_meta`
(the R9 entries of EXPERIMENTS, plus the matching baselines).

Median over the KOs, pooled R2 (in brackets median R2 within a held-out cruise):

| Features | Clear | Ahead | Similar | Worse | All 20 |
|---|---:|---:|---:|---:|---:|
| Embedding | 0.25 (0.24) | 0.29 (0.37) | 0.37 (0.36) | 0.06 (0.07) | 0.29 (0.29) |
| + temp | 0.25 (0.24) | 0.29 (0.39) | 0.39 (0.34) | 0.06 (0.07) | 0.29 (0.29) |
| + ts | 0.24 (0.19) | 0.23 (0.23) | 0.39 (0.34) | 0.05 (0.04) | 0.27 (0.22) |
| + temp + season | 0.23 (0.27) | 0.17 (0.29) | 0.37 (0.31) | 0.03 (0.02) | 0.24 (0.28) |
| + ts + season | 0.21 (0.22) | 0.15 (0.10) | 0.37 (0.35) | 0.01 (-0.02) | 0.23 (0.23) |
| + meta | 0.19 (0.22) | 0.19 (0.21) | 0.35 (0.30) | 0.10 (0.11) | 0.21 (0.21) |
| PCA 16 | 0.29 (0.26) | 0.29 (0.37) | 0.40 (0.33) | 0.18 (0.02) | 0.31 (0.24) |
| + temp | 0.30 (0.24) | 0.29 (0.33) | 0.43 (0.37) | 0.20 (0.12) | 0.30 (0.28) |
| + ts | 0.33 (0.21) | 0.28 (0.28) | 0.42 (0.34) | 0.18 (0.14) | 0.32 (0.28) |
| + temp + season | 0.27 (0.20) | 0.19 (0.23) | 0.41 (0.31) | 0.19 (0.03) | 0.29 (0.26) |
| + ts + season | 0.30 (0.20) | 0.16 (0.15) | 0.40 (0.33) | 0.18 (0.20) | 0.30 (0.23) |
| + meta | 0.32 (0.20) | 0.13 (0.07) | 0.36 (0.26) | 0.14 (0.05) | 0.27 (0.22) |
| PCA 32 | 0.27 (0.24) | 0.29 (0.37) | 0.40 (0.29) | 0.14 (-0.01) | 0.32 (0.28) |
| + temp | 0.27 (0.24) | 0.29 (0.37) | 0.42 (0.33) | 0.15 (0.06) | 0.32 (0.31) |
| + ts | 0.32 (0.24) | 0.28 (0.20) | 0.42 (0.33) | 0.14 (0.15) | 0.33 (0.23) |
| + temp + season | 0.26 (0.23) | 0.19 (0.25) | 0.41 (0.26) | 0.14 (0.09) | 0.28 (0.24) |
| + ts + season | 0.30 (0.22) | 0.17 (0.21) | 0.40 (0.27) | 0.13 (0.11) | 0.28 (0.22) |
| + meta | 0.32 (0.24) | 0.15 (0.19) | 0.37 (0.21) | 0.11 (0.12) | 0.26 (0.20) |
| ts alone | 0.13 (0.19) | 0.05 (0.28) | 0.34 (0.28) | 0.24 (0.33) | 0.23 (0.27) |
| ts + season alone | -0.04 (0.06) | -0.09 (0.12) | 0.33 (0.28) | 0.20 (0.25) | 0.16 (0.22) |
| meta alone | 0.06 (0.09) | -0.10 (-0.03) | 0.31 (0.26) | 0.25 (0.22) | 0.22 (0.16) |
| temp alone (spline) | 0.15 (0.11) | 0.10 (0.09) | 0.35 (0.33) | 0.22 (0.24) | 0.22 (0.25) |

Per KO, change from adding a block, pooled R2 (KOs improved of 20, median change):

| Added | Embedding | PCA 16 | PCA 32 |
|---|---|---|---|
| temp | 13, +0.000 | 15, +0.010 | 16, +0.008 |
| ts | 8, -0.005 | 12, +0.004 | 12, +0.011 |
| ts + season | 4, -0.032 | 5, -0.025 | 6, -0.038 |
| meta | 7, -0.016 | 8, -0.050 | 6, -0.033 |

Reading:
- **Season is what hurts.** Adding season lowers R2 with every base, both pooled
  and within cruises, and mostly on the "ahead" KOs (for example PCA 32 drops from
  0.29 to 0.17-0.19). Each cruise samples one season, so day of year acts as a
  cruise label and does not carry over to a held-out cruise. Position does the
  same, but with season already in it adds little more.
- **Temperature helps a little on the PCs** (16 of 20 KOs with PCA 32, median
  +0.01; within cruises 0.28 to 0.31). It does nothing on the full embedding,
  where it is one input among 513 under a shared penalty.
- **Salinity is mixed.** PCA 32 + T + salinity is the best pooled result so far
  (0.33), but within cruises it falls (0.31 with T alone to 0.23), mostly on the
  "ahead" KOs. So salinity seems to help between cruises and not within them.
- The "worse" KOs are still best served by temperature + salinity alone (within
  cruises 0.33); no embedding model gets there.
- **Best simple model: PCA 32 + in-situ temperature, ridge.** It is at least as
  good as the embedding in every category and is the only change that helps both
  pooled and within cruises. The gain is small, on 20 KOs. To be rechecked on
  the SST variants: if the SST channel is blanked, in-situ T may matter more.

## R10. Temperature + chlorophyll index + PCs, ridge, 20 KOs (2026-11-06)

As R8/R9: the same 20 KOs, 1,281 samples / 1,093 tiles / 13 cruises (every tile
has band means), PCA fitted on the training cruises. **chl** = R7's band-ratio
index log10(max(Rrs443, 490, 510) / Rrs560); **bands** = the 16 band means.
Files `data/results/model_experiments_colour/`. Command:
`python scripts/model_experiments.py --only col:... --out-dir data/results/model_experiments_colour`
(the `col` entries of EXPERIMENTS).

Median over the KOs, pooled R2 (in brackets median within a held-out cruise):

| Features | Learner | Clear | Ahead | Similar | Worse | All 20 |
|---|---|---:|---:|---:|---:|---:|
| Embedding | ridge | 0.25 (0.24) | 0.29 (0.37) | 0.37 (0.36) | 0.06 (0.07) | 0.29 (0.29) |
| Embedding + chl | ridge | 0.25 (0.24) | 0.29 (0.37) | 0.39 (0.36) | 0.07 (0.07) | 0.29 (0.30) |
| PCA 32 + temp | ridge | 0.27 (0.24) | 0.29 (0.37) | 0.42 (0.33) | 0.15 (0.06) | 0.32 (0.31) |
| PCA 32 + chl | ridge | 0.26 (0.23) | 0.29 (0.39) | 0.40 (0.29) | 0.14 (0.00) | 0.32 (0.28) |
| PCA 32 + temp + chl | ridge | 0.27 (0.24) | 0.29 (0.39) | 0.42 (0.33) | 0.15 (0.07) | 0.32 (0.31) |
| PCA 16 + temp + chl | ridge | 0.30 (0.23) | 0.29 (0.33) | 0.43 (0.36) | 0.20 (0.13) | 0.31 (0.28) |
| Temp + PCA 8 + chl | ridge | 0.30 (0.30) | 0.27 (0.36) | 0.42 (0.37) | 0.21 (0.06) | 0.31 (**0.32**) |
| Temp + PCA 4 + chl | ridge | 0.16 (0.20) | 0.17 (0.28) | 0.40 (0.28) | 0.25 (0.17) | 0.22 (0.24) |
| Temp + PCA 8 | CatBoost | 0.30 (0.18) | 0.24 (0.26) | 0.35 (0.34) | 0.17 (0.20) | 0.25 (0.26) |
| Temp + PCA 8 + chl | CatBoost | 0.30 (0.22) | 0.24 (0.31) | 0.35 (0.35) | 0.16 (0.25) | 0.26 (0.25) |
| Temp + PCA 4 + chl | CatBoost | 0.24 (0.21) | 0.08 (0.20) | 0.32 (0.30) | 0.11 (0.23) | 0.21 (0.22) |
| Temp | spline | 0.15 (0.11) | 0.10 (0.09) | 0.35 (0.33) | 0.22 (0.24) | 0.22 (0.25) |
| Temp + chl | spline | 0.23 (0.25) | 0.11 (0.20) | 0.35 (0.29) | 0.20 (0.26) | 0.27 (0.27) |
| Temp + chl | ridge | 0.08 (0.17) | 0.15 (0.26) | 0.31 (0.34) | 0.27 (0.27) | 0.24 (0.27) |
| Chl | ridge | -0.04 (0.05) | -0.08 (0.07) | -0.02 (0.04) | -0.03 (0.02) | -0.03 (0.05) |
| Bands (16) | ridge | -13.5 (0.17) | -0.92 (0.26) | -41 (0.06) | -26 (-0.01) | -10.6 (0.11) |

Reading:
- **The chlorophyll index adds nothing once PCs are in.** PCA 32 + temp + chl is
  the same as PCA 32 + temp, and embedding + chl the same as the embedding. The
  embedding already carries what the band ratio says.
- **A small model is as good as the best.** Temperature + 8 PCs + chl with ridge
  (10 inputs) matches PCA 32 + temp (0.31 pooled, 0.32 within cruises). With 4
  PCs it drops (0.22), so 8 or so PCs are needed.
- Temperature + chl (spline) improves on temperature alone (0.22 to 0.27 pooled),
  mostly on the "clear" KOs, but stays below the embedding. This agrees with R7.
- CatBoost is again no better than ridge on the same inputs.
- The band means need robust handling (clipping, or medians instead of means)
  before any linear model can use them (see R7).

## R11. SST channel blanked, FAIRe cruises only, ridge, 20 KOs (2026-11-06) - partial

Question: how much does the embedding lose if the SST channel is blanked, so that
every tile could be used and in-situ temperature could stand in for SST? Only the
FAIRe blank-SST bundle is back so far (the original-cruise blank bundle and both
OSTIA bundles are still on Morrill). **So this is 8 cruises, not 13**, and is not
comparable to R8-R10. Embeddings: E-faire-1103 against the same 5,553 tiles
embedded with the SST channel set to the channel mean
(`cluster/embed_sst_variants.sbatch`, `data/faire/embeddings_sst_blank/`). The
tiles chosen per sample are the same in both (option A), so this compares only
the embedding: 803 samples / 639 tiles / 8 cruises. R8's 20 KOs, ridge (spline
for temperature alone). Files `data/results/sst_blank_faire/`. Command:
`python scripts/model_experiments.py --embeddings-dirs data/faire/<embeddings|embeddings_sst_blank> --only all:emb:ridge all:pca32:ridge all:pca32+temp:ridge all:temp:spline col:temp+pca8+chl:ridge --out-dir data/results/sst_blank_faire/<same>`

Median over the KOs, pooled R2 (in brackets median within a held-out cruise):

| Features | Current SST | SST blanked |
|---|---:|---:|
| Embedding | 0.26 (0.34) | -0.08 (0.03) |
| PCA 32 | 0.24 (0.36) | -0.03 (0.18) |
| PCA 32 + in-situ T | 0.26 (0.37) | 0.15 (0.23) |
| T + PCA 8 + chl | 0.26 (0.27) | 0.11 (0.26) |
| In-situ T alone (spline) | 0.22 (0.20) | 0.22 (0.20) |

Reading:
- **Blanking SST wrecks the embedding** (0.26 to -0.08 pooled). Even the
  "similar" KOs drop, from 0.37 to 0.06, so this is not just KO selection.
- Adding in-situ temperature back recovers only part of it (0.15), and stays
  below temperature alone (0.22). So in-situ T cannot simply replace SST inside
  the embedding.
- A likely reason: Granite never saw tiles with a constant SST channel, so
  blanking moves the whole embedding off its training distribution, not just the
  temperature part.
- **Leaning: keep SST in the embedding, and drop samples rather than blank.** The
  OSTIA fill is the remaining option for keeping more samples.
- To redo on 13 cruises when the original-cruise blank bundle is back.
