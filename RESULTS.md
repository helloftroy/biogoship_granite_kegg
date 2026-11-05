# Results

Every result is recorded with what it assumed and what went into it, so later
runs can be compared without re-running old ones. Each entry has: **inputs**
(data snapshot and sample set), **settings** (tile rule, SST handling, features,
model, validation), the **result tables**, **files**, and the **command** that
made it. Shared definitions are in section 0; an entry lists only how it differs.

Narrative, methods and earlier history (2026-10-20 to 10-23) are in
[METHODS_RESULTS.md](METHODS_RESULTS.md). Figures are indexed in [FIGURES.md](FIGURES.md).

Results that will change, and why:
- **SST from moved windows** (`check_sst_nearby.sbatch`, running on Morrill
  2026-11-04): may recover SST for more samples, especially I07N; R1 and R3 would
  then gain samples.
- **In-situ temperature as a baseline**: not yet run; tells whether a KO is
  predictable beyond what temperature explains.
- **CatBoost on a KO shortlist**: R3 uses ridge and nearest-neighbour models only.

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
