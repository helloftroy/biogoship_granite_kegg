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
- **SST variants** (R11, R12): blanking loses a lot; the OSTIA fill adds samples
  without loss and is proposed as the new standard. R1-R3, R5-R10 are on the
  pre-OSTIA embeddings.
- **Ocean-model fields** (GLORYS12, PISCES; `cluster/model_fields.sbatch`): pending.
- **External test sets** (R13-R15, R29): tiles and embeddings done for all three sets.
  JGI tested (R29): absolute levels do not transfer, weak skill with the offset
  removed. Tara Pacific and EMO BON need KO annotation with our pipeline.
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
| KO-1104 | KO counts `ko_counts_cov50_all.tsv` (repo root, not in git) | received 2026-11-04 | 3,081 libraries x 18,835 KOs; gene counts per KO, a gene with n KOs counting 1/n (not read-weighted; see R15) |
| PF-1104 | Pfam counts `pfam_counts_cov50_all.tsv` (repo root, not in git) | received 2026-11-04 | 3,084 libraries x 13,161 families (not yet analysed) |
| E-orig-l4-1106 | Original cruises re-embedded, SST gaps from OSTIA (`data/embeddings_sst_l4/`) | 2026-11-06 | 3,503 tiles; 556 OSTIA-only tiles; log `data/sst_l4_log.csv` |
| E-faire-l4-1106 | FAIRe cruises re-embedded, SST gaps from OSTIA (`data/faire/embeddings_sst_l4/`) | 2026-11-06 | 5,553 tiles; 1,552 OSTIA-only tiles; log `data/faire/sst_l4_log.csv` |
| E-blank-1106 | Both sets re-embedded with the SST channel blanked (`data/embeddings_sst_blank/`, `data/faire/embeddings_sst_blank/`) | 2026-11-06 | 3,503 + 5,553 tiles |
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

## R12. SST variants on 13 cruises: current, blanked, OSTIA-filled, ridge, 20 KOs (2026-11-06)

Inputs: E-orig-1102 + E-faire-1103 (current), and the same tiles re-embedded by
`cluster/embed_sst_variants.sbatch`. The variants are SST channel blanked
(`data/*embeddings_sst_blank/`) and SST gaps filled from OSTIA shifted toward
SLSTR skin SST (`data/*embeddings_sst_l4/`; logs `data/sst_l4_log.csv`,
`data/faire/sst_l4_log.csv`). OSTIA changed only tiles with **no** SLSTR SST at all
(556 original + 1,552 FAIRe tiles), because no tile had partial SST coverage. The
offset applied was the tile's own where it had SLSTR pixels, otherwise the cruise
median; the median SLSTR - OSTIA difference is -0.14 K (s.d. 0.5-0.6 K).

Two sample sets:
- **Shared:** the current tiles (option A), 1,281 samples / 1,093 tiles.
- **Full:** with OSTIA every tile counts as having SST, so option A picks the
  closest tile in time: 1,388 samples / 1,198 tiles, 253 of them OSTIA-only.
  "Blank, full" uses the same tile choice with the SST channel blanked.

R8's 20 KOs, ridge (spline for temperature alone). Files
`data/results/sst_variants/<variant>/`. Commands:
`python scripts/model_experiments.py --embeddings-dirs data/<dir> data/faire/<dir> --only all:emb:ridge all:pca32:ridge all:pca32+temp:ridge all:temp:spline col:temp+pca8+chl:ridge --out-dir data/results/sst_variants/<variant>`.
The views `embeddings_sst_l4_shared` (OSTIA embeddings, current
sst_valid_fraction) and `embeddings_sst_blank_full` (sst_valid_fraction = 1)
contain only a sample table and a link to the embeddings.

Median over the KOs, pooled R2 (in brackets median within a held-out cruise):

| Features | Current | Blank | OSTIA, shared | OSTIA, full | Blank, full |
|---|---:|---:|---:|---:|---:|
| Embedding | 0.29 (0.29) | 0.10 (0.04) | 0.29 (0.29) | **0.32** (0.32) | 0.14 (0.10) |
| PCA 32 | 0.32 (0.28) | 0.15 (0.11) | 0.32 (0.28) | 0.34 (0.29) | 0.18 (0.15) |
| PCA 32 + temp | 0.32 (0.31) | 0.25 (0.22) | 0.32 (0.31) | **0.34** (0.32) | 0.26 (0.29) |
| Temp + PCA 8 + chl | 0.31 (0.32) | 0.25 (0.28) | 0.31 (0.32) | **0.34** (0.34) | 0.28 (0.34) |
| Temp alone (spline) | 0.22 (0.25) | 0.22 (0.25) | 0.22 (0.25) | 0.23 (0.28) | 0.23 (0.28) |

"OSTIA, shared" equals current exactly, as expected: none of its tiles changed.

Do the OSTIA-only tiles predict as well? In the full set, scored apart (median
over KOs of pooled R2 within each group):

| Features | Trained on | SLSTR tiles | OSTIA-only tiles | SLSTR tiles, same cruises |
|---|---|---:|---:|---:|
| Embedding, ridge | all tiles | 0.30 | 0.40 | 0.33 |
| Embedding, ridge | SLSTR tiles only | 0.27 | 0.29 | 0.29 |
| PCA 32 + temp, ridge | all tiles | 0.33 | 0.37 | 0.33 |
| Temp + PCA 8 + chl, ridge | all tiles | 0.31 | 0.40 | 0.31 |
| Temp alone, spline | all tiles | 0.23 | 0.31 | 0.20 |

(The OSTIA-only tiles are spread over 12 cruises, most in A16N-2023, I07N and
C13.5/A13.5. They score higher even for temperature alone, so R2 on them is not
directly comparable to R2 on the SLSTR tiles; compare each row against the
temperature row.)

Reading:
- **Blanking SST is out.** It costs 0.15-0.2 R2 everywhere, and in-situ
  temperature only partly makes up for it (R11 now confirmed on 13 cruises).
- **The OSTIA fill works.** It adds 107 samples and 105 tiles. The embedding on
  OSTIA-filled tiles keeps its lead over temperature, and the gain over
  temperature on those tiles is as large as on SLSTR tiles. Training on them
  also improves the SLSTR tiles (0.27 to 0.30).
- **Suggest the OSTIA-filled embeddings as the new standard** (E-orig-l4 +
  E-faire-l4, option A = the closest tile in time). Re-baseline R1-R3 and the
  KO screen on them when convenient. Best simple models on it: PCA 32 + temp or
  temp + PCA 8 + chl, ridge (0.34).

## R13. JGI multi-environment KOs (Flinkstrom et al. 2024, mSystems): can they be a test set? (2026-11-06)

Source: figshare project 187989, `adata_raw_filtered_v2.h5ad` (raw, unnormalised KO
abundances, 6,539 metagenomes x 9,871 KOs; md5 9c4e85f1...) and the unfiltered
`adata_raw.h5ad` (10,160 x 18,290). Both are in `data/external/jgi_flinkstrom/`,
not in git.

Samples usable with Sentinel-3 (marine water column, not vents, sediment or mats;
collected from May 2016 on; with coordinates; depth <= 20 m or not given):

| Step | Filtered v2 | Unfiltered |
|---|---:|---:|
| Marine | 1,184 | 1,938 |
| Water column | 955 | |
| From 2016-05 | 156 | |
| + depth <= 20 m or unknown | **111** | **116** |

The 111 are 83 Monterey Bay (coastal, 2016-09-26 to 11-16, S3A only), 11
Marineland FL (Aug 2016), 7 South Atlantic transect, 6 Arctic Ocean (Jan-Feb 2020,
polar night, so no OLCI), 3 North Sea, 1 other. Dates have no time of day. Most
JGI marine metagenomes predate Sentinel-3 (2011-2016).

KO tables compared:
- IDs are `KO:Kxxxxx` (ours `count.ko.Kxxxxx`). The JGI KEGG release is older
  (highest K19791, ours K22512).
- 3,705 of our 3,921 target KOs are in JGI. 74% of the missing ones are KO
  numbers newer than JGI's release.
- All 10 single-copy marker KOs (K06942, K01889, K01887, K01875, K01883, K01869,
  K01873, K01409, K03106, K03110; Salazar et al. 2019) are in both. Their spread
  within a sample is similar (median CV 0.30 JGI, 0.34 ours).
- Copies per cell (KO / median of the 10 markers), per-KO median, ours against
  the 111 JGI samples: Pearson r 0.80, median log10 difference +0.04. So the
  scales agree overall. A few KOs differ by about 1000-fold, probably annotation
  differences (for example K02257, K18889/K18890, K13797 lower in JGI); these need
  screening before any transfer.
- 77% of the shared target KOs are present in at least 90% of the 111 samples.

Normalisation note: dividing by the marker median is a per-sample scale, and CLR
removes any per-sample scale. So "marker normalisation, then CLR" equals CLR,
except for where the pseudocount lands. What matters for a joint CLR is the
**same KO set** (the shared KOs) and the same pseudocount rule in both datasets.
An alternative that does not depend on the KO set is log copies per cell.

## R14. Other external sets: Tara Pacific and EMO BON (2026-11-06) - assessment

**Tara Pacific** (2016-2018). The plankton metagenomes are ENA PRJEB64422
(Genoscope, public 2025-05-19). There are 389 runs from 387 samples, all in the
0.2-3 um fraction, NovaSeq, about 40 Gbp per run. **Raw reads only**: no KO or
gene-abundance table was found. Sample titles carry the UTC timestamp and the
station; ENA has lat/lon (383 of 387) but no depth or collection_date fields.
Sampling depth from Lombard et al. 2023: about 2 m (pump), or 1.5 m (underway
inlet).

| Type (sample title code) | Samples | 2016 / 2017 / 2018 |
|---|---:|---|
| O-SRF, open-ocean surface | 36 | 21 / 8 / 7 |
| I-SRF, island surface water (within 200 nmi of an island) | 75 | 29 / 36 / 10 |
| S-SRF, reef surface water | 95 | 30 / 44 / 21 |
| C-CSW, water just above coral colonies | 181 | 55 / 88 / 38 |

Open-ocean plus island water is 111 samples. Reef and coral water (276) sits next
to reefs and land, where a 300 m tile is mostly shallow bottom or land.

**EMO BON** (from June 2021; coastal observatories in Europe and the Red Sea).
Water is filtered at 0.2-3 um and 3-200 um, with 2 replicates sequenced.
Per-sample MetaGOflow outputs are RO-Crates on GitHub (emo-bon/analysis-results-
cluster-01-crate: 151 water crates; cluster-02-crate: 170). Files sit in S3 at
s3.mesocentre.uca.fr, and each crate's ro-crate-metadata.json gives downloadUrls.
The cluster-02 bucket (mgf-data-products-2) refused connections on 2026-11-06.
- Each crate has `*.merged.summary.ko`: counts of predicted CDS on **merged reads**
  per KO (KOfam hmmsearch with gathering thresholds). In 4 samples checked there
  are 2,271-2,403 KOs and 0.1-0.9 million KO hits.
- **6 of the 10 single-copy markers are always 0** (K06942, K01883, K01869,
  K01873, K03106, K03110), probably because KOfam thresholds fail on read-length
  fragments. So the profiles are not like ours, and marker normalisation is not
  possible as is.
- Crates also hold `final.contigs.fa` (assembly) and the CDS protein files, so
  contigs could be re-annotated with our own KO pipeline.
- Sampling metadata are in the observatory crates (`logsheets/raw/water_sampling.csv`:
  collection_date, depth, size_frac, ENA accession; `water_observatory.csv`:
  position; `water_measured.csv`: in-situ variables).

---

## R15. External test sets: manifests, JGI targets, what the KO pipeline is (2026-11-07)

**Manifests** (`scripts/build_external_manifests.py`; `data/external/<set>/manifests/`,
`samples.csv` lists every candidate with the reason it was dropped):

| Set | Samples | Sites | Groups (section) | Dates | Time |
|---|---:|---:|---|---|---|
| JGI Flinkstrom 2024 | 109 | 57 | 6 studies; 83 samples Monterey Bay (one position), 11 Marineland FL, 7 S Atlantic, 4 North Sea, 3 Arctic, 1 other | 2016-04-26 to 2020-02-06 | date only |
| Tara Pacific O-SRF + I-SRF | 111 | 108 | 30 (open ocean I00 = 36; islands I03-I32) | 2016-07-12 to 2018-08-31 | UTC minute |
| EMO BON water with metaGOflow results | 307 | 86 | 15 observatories | 2021-06-18 to 2022-10-27 | date only |

- JGI: Ecosystem Type Marine, Habitat seawater/marine, depth <= 20 m, on or after
  2016-04-25. The filtered h5ad has 109 (the unfiltered one had 116).
- EMO BON: both size fractions (0.2-3 um 146, 3-200 um 154, unknown 7) and both
  replicates; they share a site per sampling. Depth 0-2 m except 8 at 18.7 m.
- Date-only times are set to local solar noon (12:00 - lon/15 h UTC).
- Tiles: fixed position only (no back-tracking), OSTIA-filled (R12 standard).
  `fill_sst_l4.py --split-gap-days 5` makes one OSTIA subset per run of tiles,
  so the Tara open-ocean group does not need a Pacific-wide 2-year download.

**JGI targets** (`scripts/build_external_targets.py`; `data/external/jgi_flinkstrom/targets/`):
- Shared KO set: our 3,921 target KOs, 3,705 in the JGI table, 2,832 present in
  >= 90% of the 109 JGI samples. Both datasets are CLR'd over those 2,832 (ours by
  re-centring the existing CLR, which is identical). Marker normalisation is not
  applied: it cancels in CLR.
- Mean CLR per KO, JGI vs Bio-GO-SHIP: r = 0.82; median difference +0.01;
  123 KOs differ by more than 2 (about 7-fold). Within-dataset spread is the same
  (median SD 0.49 vs 0.50).
- All 10 single-copy markers are in the shared set; JGI is 0.1-0.5 lower in CLR.
- 15 of the 20 screen KOs are in the shared set. Missing: K00436 (in 16% of JGI
  samples), K16924 (10%), K17710 (59%), K01530 (69%), K02999 (71%). K02999 and
  K17710 are eukaryotic, so the JGI samples hold fewer eukaryotes.

**What the KO pipeline is** (example output `A392013_008_meta`, "coverage0",
received 2026-11-07; the code is not available yet):
- Per taxonomic level (domain to species), tables of `taxon, ko|og|pfam, length, count`.
  The IDs are eggNOG-mapper style (`ko:K08332`, COG/KOG and eggNOG v5 OGs, Pfam names).
- Counts are gene counts: a gene with n KOs adds 1/n to each. Lengths are summed
  amino acids (median about 200 per gene).
- The taxonomy looks like PR2/MarFERReT for eukaryotes and NCBI/GTDB-like names
  for prokaryotes.
- Against our cov50 row for the same library: every cov50 KO is in coverage0, and
  log counts correlate r = 0.93. cov50 keeps 52% of all counts and 39-57% of the
  marker KOs.
- The kept share falls with gene length: 0.68 for KOs whose genes average under
  140 aa, 0.13 for those over 565 aa. This fits a >= 50% alignment-coverage
  filter, but it is not certain.
- To reprocess Tara Pacific and EMO BON the same way: assemble, predict genes,
  eggNOG-mapper, count genes per KO with the 1/n split, and apply the cov50 rule
  once it is known.

---

## R16. All KOs, Granite vs in-situ temperature, on the OSTIA-filled embeddings (2026-11-07)

As R6, but the embeddings are E-orig-l4-1106 + E-faire-l4-1106 and the tile is the
closest in time (all tiles have SST). That gives 1,429 samples, 1,216 tiles and
13 cruises. Learners are fixed rather than the best per KO: the embedding uses
ridge, temperature uses spline + ridge. Held-out predictions are saved. Files
are in `data/results/ko_temperature_all_ostia/`. Command:
`python scripts/test_ko_temperature.py --kos all --learners ridge spline --feature-sets temp embedding --embeddings-dirs data/embeddings_sst_l4 data/faire/embeddings_sst_l4 --save-predictions --out-dir data/results/ko_temperature_all_ostia`
(about 1 minute).

The same rule as R6 (difference 0.1, 9 of 13 cruises), with the same learners, on
the old and new embeddings:

| | Pre-OSTIA (R6 data) | OSTIA |
|---|---:|---:|
| Median pooled R2: embedding / temperature | 0.07 / -0.01 | 0.10 / 0.01 |
| KOs with embedding R2 >= 0.25 | 692 | 849 |
| Embedding clearly better | 757 | 747 |
| Temperature clearly better | 369 | 280 |
| Both >= 0.5 and within 0.05 | 45 | 59 |

Per-KO embedding R2, old vs new: r = 0.95.

The heatmap's stricter categories are in `ko_categories.csv`. They are: Granite
better by >= 0.15 in >= 10 cruises (360 KOs); both >= 0.6 and within 0.05 (37);
temperature better by >= 0.1 in >= 9 cruises (280). Since 2026-11-11 the Granite and temperature
categories also need the better input above R2 0.3: 87 and 26 KOs. Most of the 280
were KOs that neither input predicts (temperature R2 median 0.02).

Reading (figures: FIGURES.md, cruise line and KO heatmap):
- **The OSTIA fill helps the embedding across all KOs.** 157 more KOs reach
  R2 >= 0.25, and 89 fewer are clearly better predicted by temperature.
- **The best-predicted KOs are cyanobacterial**, and temperature predicts them as
  well as Granite does:
  - NDH-1 subunits ndhD/J/K
  - hoxH
  - phycoerythrin linker cpeS
  - Mn transport mntB/C
  - nrdJ
  - alkane synthesis K14330
  In the example metagenome, 90-100% of their genes are cyanobacterial.
- **Granite clearly wins on heterotroph genes** (mostly Pseudomonadota in the
  example metagenome):
  - phosphate transport pstS/A/B/C, which plausibly track P limitation
  - amino-acid synthesis (arg, lys, pro, aromatic)
  - purine synthesis
  - glutaminase glsA
  On P18, glsA follows the subtropical-gyre maximum near 25 S, which temperature
  misses.
- **Temperature wins on:**
  - eukaryotic and giant-virus genes: POLD1, TOP2, RRM2, NOP2, GPAT
  - purple-bacterial reaction-centre genes: puhA, pufL
  - the siderophore receptor fiu
  - ndhF

---

## R17. All KOs: time & place, raw satellite and the best model, OSTIA embeddings (2026-11-07)

As R16 (1,429 samples, 1,216 tiles, 13 cruises), adding three feature sets to
`scripts/test_ko_temperature.py`:
- **metadata:** time & place, as in R5.
- **satellite:** 16 OLCI band means + tile SST mean. Each is rank-normalised on the
  training cruises, because some FAIRe means are tens of s.d. out. Missing values
  are set to the training median. For OSTIA-only tiles, SST is OSTIA at the tile
  centre plus the applied offset, from `sst_l4_log.csv`.
- **pca32_temp:** the best model so far (R12); the PCA is fitted on the training cruises.

Files are in `data/results/ko_temperature_all_ostia_v2/`. Command:
`python scripts/test_ko_temperature.py --kos all --learners ridge spline knn --feature-sets temp embedding metadata satellite pca32_temp --embeddings-dirs data/embeddings_sst_l4 data/faire/embeddings_sst_l4 --save-predictions --out-dir data/results/ko_temperature_all_ostia_v2`
(about 6 minutes). Embedding and temperature results are identical to R16.

Median pooled R2 over all 3,921 KOs:

| Features | ridge | spline | kNN |
|---|---:|---:|---:|
| 32 PCs + temperature | **0.12** | - | 0.02 |
| Embedding | 0.10 | - | 0.00 |
| Raw satellite | 0.05 | 0.05 | -0.03 |
| Temperature | -0.01 | 0.01 | -0.18 |
| Time & place | -0.14 | -0.63 | -0.32 |

Reading (heatmap with four predictors, FIGURES.md):
- **Time & place does not transfer to a new cruise.** It is below 0.2 for every
  KO in the heatmap.
- **Raw satellite needs the spline.** It is nonlinear in SST; with ridge it scores
  about 0 on the cyanobacterial KOs, with the spline 0.42-0.54.
- **On the heterotroph KOs where Granite beats temperature, the raw satellite
  values do almost as well as Granite.** Examples: pstB 0.43 vs 0.46, glsA 0.48 vs
  0.46, citramalyl-CoA lyase 0.56 vs 0.44. So most of Granite's lead over
  temperature on those KOs is information in the tile's mean colour and SST. The
  embedding's own extra over the band means is small there.
- **Where both are high (cyanobacteria), Granite (about 0.7) beats raw satellite
  (about 0.5).**
- **On P18, the best model adds most where temperature fails.** glsA rises from
  -0.15 to 0.47; cpeS stays at 0.89 vs 0.92. For fiu, the best model is worse
  than temperature alone (0.30 vs 0.52): it uses temperature linearly. A
  spline-temperature + PCs model is the obvious next variant.


## R18. Model development on the figure KOs (2026-11-07)

Question: can a better model than R17's (Granite ridge; 32 PCs + linear temperature)
be found for the 54 heatmap KOs (28 Granite-better, 14 both-high, 12 temperature-better;
the three cruise-line genes are among them)? Same samples and cruise-held-out CV as R17.
Script `scripts/tune_figure_kos.py`; files `data/results/figure_ko_models/`. Models
(ridge throughout; "~" = cubic spline, 8 knots; "grp" = penalty chosen by leave-one-cruise-out
inside the training cruises instead of leave-one-sample-out; avg = mean of two predictions)
are listed in the script. "Best per KO" picks a model per KO and held-out cruise by
cruise-held-out CV on the other 12 cruises (nested).

Median pooled R2 (selected rows; "within" = median R2 within a held-out cruise):

| Model | All 54 | Granite-better | Both high | Temp-better | Within |
|---|---:|---:|---:|---:|---:|
| avg(32 PCs + temp~, sat~) | **0.47** | 0.41 | 0.65 | **0.48** | 0.42 |
| avg(emb, temp~ + sat~) | 0.45 | 0.42 | 0.70 | 0.44 | 0.42 |
| emb, grp penalty | 0.45 | 0.42 | 0.70 | 0.39 | **0.44** |
| best per KO (nested) | 0.44 | 0.39 | 0.71 | 0.47 | 0.43 |
| 32 PCs + temp~ | 0.44 | 0.39 | 0.71 | 0.46 | 0.41 |
| emb (R17 Granite) | 0.43 | 0.42 | 0.68 | 0.33 | 0.43 |
| 32 PCs + temp (R17 best) | 0.43 | 0.40 | 0.69 | 0.39 | 0.44 |
| 16 PCs + temp~ | 0.42 | 0.37 | 0.70 | 0.46 | 0.41 |
| sat~ (bands + SST) | 0.41 | 0.35 | 0.46 | 0.42 | 0.34 |
| sst~ (SST mean only) | 0.34 | 0.23 | 0.55 | 0.44 | 0.34 |
| temp~ | 0.26 | 0.16 | 0.69 | 0.45 | 0.33 |
| bands~ (16 band means only) | 0.24 | 0.31 | -0.15 | 0.02 | 0.24 |
| ~8 PCs + temp~ (a curve per PC) | 0.38 | 0.31 | 0.66 | 0.40 | 0.34 |

Reading:
- **Why raw satellite now looks strong (R7/R10 had it near or far below 0):** R17's
  "raw satellite" includes the tile SST mean, filled by OSTIA, and is rank-normalised,
  which removes the extreme FAIRe band means that broke linear models in R7/R10. Split
  up: on the Granite-better KOs the bands alone give 0.31 and SST alone 0.23, together
  0.35, against Granite 0.42. On the cyanobacterial KOs the bands alone fail (-0.15)
  and SST carries it (0.55). So the colour bands do carry part of Granite's lead on
  the heterotroph KOs, and Granite still adds about 0.07 over them (median).
- **A curved temperature term fixes the temperature-better KOs.** 32 PCs + temp~
  raises them from 0.39 to 0.46 (temperature alone 0.45), without losing elsewhere.
- **Choosing the penalty across cruises helps the embedding** (0.43 to 0.45, mostly on
  the temperature-better KOs, 0.33 to 0.39) and gives the best within-cruise score.
- **Averaging Granite with the raw satellite model** gives the best pooled median
  (0.47), mainly by reaching temperature-level scores on the temperature-better KOs.
- **Per-KO choice does not beat one good model** (0.44 nested). With 13 cruises the
  choice is noisy; a single model (32 PCs + temp~, or the average) is preferable.
- Curves on the PCs, fewer PCs (8) and 64 PCs do not help.
- P18 (cruise line, nested per-KO choice vs temperature): cpeS 0.92 vs 0.92, glsA 0.47
  vs -0.15, fiu 0.43 vs 0.52. 32 PCs + temp~ on P18: 0.91, 0.48, 0.44 (pooled over all
  cruises fiu is 0.51 vs temperature 0.51).
- Next: check the two leading models (32 PCs + temp~; the average with raw satellite)
  on all 3,921 KOs, and the grp penalty on them.

## R19. Why the Granite model misses the pst dip on P18: phosphate (2026-11-09) - exploratory

The pst genes (high-affinity phosphate uptake: pstS, pstC, pstA, pstB) drop on P18 from
the equator to about 18 S, and the Granite model (32 PCs + temp~, R18) does not follow.
Exploratory scripts (scratchpad, not in git); to be redone with PISCES phosphate.

- **That stretch is warm and phosphate-replete:** in-situ phosphate 0.5-0.9 umol/kg at
  22-25 C (eastern equatorial Pacific upwelling). In every other cruise, water above
  20 C has phosphate of at most 0.43 (median 0.00-0.13). So no training cruise
  has warm, phosphate-replete water, and neither temperature nor the embedding has
  learned what it looks like.
- pstC vs log10 phosphate: r = -0.78 over the 451 tiles with phosphate (-0.72 in water
  above 20 C); vs temperature r = 0.58.
- Cruise-held-out ridge on the 451 tiles with phosphate (12 cruises), pstC pooled / P18:
  temperature 0.09 / 0.44; phosphate (spline) 0.68 / 0.88; temp + phosphate 0.64 / 0.80;
  32 PCs + temp 0.24 / 0.43; 32 PCs + temp + phosphate in one ridge 0.46 / 0.51 (the PCs
  swamp it); **temp + phosphate first, then 32 PCs on the residual (cruise-held-out
  penalty) 0.65 / 0.80**. pstB, pstS, pstA alike. With phosphate the equatorial dip
  is predicted (somewhat too deep: -0.5 to 0.2 against observed 0.3 to 1.0).
- Phosphate does not help cpeS or NOP2 (temperature already explains them).
- **Next:** PISCES phosphate (cluster/model_fields.sbatch) as the input that exists
  everywhere, two-stage as above; check it reproduces the in-situ result.

## R20. Phosphate-first Granite models, and does cruise overlap help? (2026-11-09)

**Phosphate models** (`scripts/phosphate_models.py`, files `data/results/phosphate_models/`):
the 54 figure KOs on the 452 tiles with in-situ phosphate (12 cruises), cruise held out.
Pooled R2:

| KO | emb | 32 PCs + temp~ | P~ then 32 PCs | T~P~ then 32 PCs |
|---|---:|---:|---:|---:|
| pstC | 0.24 | 0.24 | **0.71** | 0.65 |
| pstB | 0.31 | 0.38 | 0.66 | **0.67** |
| pstS | 0.36 | 0.36 | 0.70 | **0.75** |
| pstA | 0.31 | 0.32 | 0.64 | **0.66** |
| cpeS | 0.72 | **0.72** | 0.40 | 0.54 |
| NOP2 | **0.47** | 0.46 | 0.13 | 0.46 |

Phosphate first is best for pstC (P18: 0.89, against temperature 0.43); for the other
pst genes temperature + phosphate first is as good or slightly better. For KOs not
driven by phosphate, phosphate first hurts (median over the heatmap's "both high"
KOs 0.37 against 0.73), so it is a per-KO choice. The cruise map (FIGURES.md) uses
32 PCs + temp~ for cpeS and NOP2 and phosphate first for pstC.

**Overlap between cruises** (`figures/make_skill_map.py`; per-tile table
`data/results/ko_temperature_all_ostia_v2/sample_skill.csv`). Per-sample skill = 1 -
mean over KOs of (error^2 / variance of the KO over all tiles), for the 920 KOs with
pooled R2 >= 0.25 under pca32_temp ridge; its mean over tiles (0.37) is the mean pooled R2.

- **Geographic overlap does not help.** Skill vs distance to the nearest sample of
  another cruise: Spearman 0.02 (1,216 tiles); per cruise, the share of samples within
  200 km of another cruise vs per-cruise median R2: 0.08 (13 cruises). I09N (2016) lies
  87% within 200 km of I09N-2025 and is not better predicted than average.
- **Similar conditions do.** Skill vs mean distance to the 5 nearest samples of other
  cruises in (8 embedding PCs, temperature): Spearman -0.28; median skill 0.62 in the
  most typical fifth, 0.38 in the most unusual. With temperature alone: -0.34.
- So what matters is whether the training cruises contain similar water, not the same
  place (as R19's warm, phosphate-replete equatorial Pacific).
- Per-sample skill and per-cruise R2 rank cruises differently: per-cruise R2 is
  relative to the variation within that cruise, skill to the variation over all
  cruises. HLY2202 (Arctic): per-cruise R2 0.64, median skill -0.06 (its large offset
  from the other cruises is partly but not fully predicted); P02-2022 (30 N, little
  variation along it): per-cruise R2 -0.11, skill 0.73.

## R21. Which KOs Granite or temperature predicts better: by carrier; taxonomy tables (2026-11-11) - exploratory

**Carriers.** Taken over the 664 KOs that one of the two inputs predicts with R2 >= 0.3
(Granite embedding ridge vs temperature spline, R16 OSTIA): Granite better by >= 0.1
in 210, temperature in 32, the rest within 0.1. Shares of each KO's genes per order
come from one example metagenome (A392013_008, `level_4_order_ko_counts`):

| Order | Mean share in Granite-better KOs | in temperature-better KOs | r with (Granite - temperature) |
|---|---:|---:|---:|
| SAR11 / Pelagibacterales | 0.24 | 0.02 | +0.58 |
| Pseudomonadales | 0.08 | 0.06 | +0.23 |
| Rhodobacterales | 0.10 | 0.07 | +0.18 |
| Flavobacteriales | 0.05 | 0.02 | +0.18 |
| Mamiellophyceae | 0.05 | 0.20 | -0.29 |
| Synechococcus / Prochlorococcus | 0.04 | 0.10 | -0.24 |
| Prymnesiophyceae | 0.04 | 0.06 | -0.22 |

- Granite is ahead on KOs carried by heterotrophic bacteria, above all SAR11;
  temperature on phytoplankton KOs (cyanobacteria and eukaryotic picophytoplankton),
  viral KOs (TOP2, RRM2, POLD1) and the aerobic anoxygenic phototroph reaction centre.
- By pathway, Granite leads on amino-acid synthesis and tRNA charging
  (aminoacyl-tRNA biosynthesis 12 of 12 KOs); photosynthesis, oxidative
  phosphorylation and eukaryotic cell machinery are even or lean to temperature.
- KOs spread across many taxa favour Granite (median difference 0.10 when the main
  phylum holds < 71% of genes, 0.01 above).
- Phosphate first (R20) helps only pstS, pstA, pstB, pstC (+0.27 to +0.47).
- Caveat: the shares come from one sample.

**Taxonomy tables.** `scripts/collect_taxa_tables.py` collects the pipeline's
`level_<n>_<rank>_abundance.tsv` (read counts per taxon, GTDB-Tk for MAGs, EUKulele
otherwise; counts summed over the two methods) into `data/targets/taxa/`
(`taxa_long.csv.gz`, `<rank>_counts.tsv`, `libraries.csv`). `build_ko_targets.py`
then makes CLR targets (`taxa_<rank>_clr.csv`, taxa in >= 90% of samples).
Full set (coverage_50 folder, 2026-11-11): abundance tables for 3,094 libraries ->
2,011 surface samples on 16 cruises, covering all 2,007 KO samples. Taxa in >= 90% of
samples: 75 phyla, 161 classes, 253 orders, 315 families, 429 genera. Unclassified
reads: median 17% (IQR 12-23%). The earlier abundance tables were identical to the
coverage_50 ones.
KO-by-taxon tables (`level_2` phylum and `level_4` order `ko_counts`) are stashed per
library in `data/targets/taxa/ko_by_taxon/` (gzipped; OG, Pfam and other levels dropped):
1,431 libraries so far. They are coverage_50: summed over taxa they reproduce
`ko_counts_cov50_all.tsv` exactly (6 libraries checked). Coverage per library is in
`data/targets/taxa/libraries.csv`.

## R22. Taxonomy as targets: Granite vs temperature per taxon (2026-11-11)

Same setup as R17 (OSTIA-filled embeddings, one of 13 cruises held out, pooled R2;
embedding ridge, temperature spline, pca32_temp ridge), with the coverage_50 taxon
read-count CLRs of R21 as targets (taxa in >= 90% of samples). 1,432 samples,
1,219 tiles. Command per rank:
`python scripts/test_ko_temperature.py --kos all --targets data/targets/taxa/taxa_<rank>_clr.csv --learners ridge spline --feature-sets temp embedding pca32_temp --embeddings-dirs data/embeddings_sst_l4 data/faire/embeddings_sst_l4 --save-predictions --out-dir data/results/taxa_temperature/<rank>`

Median pooled R2 over taxa:

| Rank | Taxa | Temperature | Granite | PCs + temperature |
|---|---:|---:|---:|---:|
| Phylum | 75 | 0.12 | 0.09 | 0.14 |
| Class | 161 | 0.05 | 0.05 | 0.07 |
| Order | 253 | 0.05 | 0.06 | 0.09 |
| Family | 315 | 0.12 | 0.10 | 0.13 |
| Genus | 429 | 0.27 | 0.28 | 0.29 |

Domain (CLR over Archaea, Bacteria, Eukaryota, Viruses): Bacteria 0.35 Granite vs
0.05 temperature; Archaea 0.18 vs -0.05; Eukaryota 0.48 vs 0.55; Viruses neither.

- **R21's guess (Granite = heterotrophs, temperature = phytoplankton) does not hold
  for taxon abundances.** Most taxa are ties. SAR11 genera split both ways
  (Pelagibacter 0.21 vs 0.12, MED727 0.71 vs 0.68, TMED170 0.48 vs 0.54).
  Synechococcus and Prochlorococcus are ties or slightly Granite (Prochlorococcus_A
  0.77 vs 0.69). Granite's clearest wins are eukaryotic phytoplankton: Pelagophyceae
  (class 0.36 vs 0.03; Pelagomonas 0.48 vs 0.18) and Mamiellophyceae (0.43 vs 0.14).
  Temperature's are few and small (Nitrosopumilus 0.34 vs 0.17, some SAR86 and SAR11
  sub-clades).
- **What the Granite-better KOs follow is the bacterial share of the community.**
  Across samples, the KO CLR correlates with log(Bacteria / Eukaryota) at median
  r = +0.73 for Granite-better KOs and -0.62 for temperature-better KOs (+0.55 for
  "both high", +0.19 for all KOs). The bacterial share is something Granite predicts
  (0.35) and temperature does not (0.05). So the KO split is mostly bacteria vs
  eukaryotes, not heterotrophs vs phototrophs.

## R23. Unclassified reads: what drives them, are they predictable, do they explain R22? (2026-11-11)

Unclassified = domain-level reads in "Unclassified", "Unclassified Archaea" and
"Unclassified Bacteria" (coverage_50 abundance tables), per sample (libraries averaged).
Targets in `data/targets/taxa/unclassified_targets.csv`; models as R22, results in
`data/results/taxa_temperature/unclassified/`.

- **Share per sample:** median 0.19 (IQR 0.15-0.25, max 0.89). By cruise from 0.11
  (AE1319) and 0.12 (I09N 2016) to 0.26 (A16N) and 0.31 (HLY2202, Arctic).
- **Environment:** higher in cold, high-latitude and eukaryote-rich water (Spearman
  with temperature -0.55, |latitude| +0.46, eukaryote CLR +0.60, log Bacteria/Eukaryota
  -0.63). Not related to library size (+0.02). Eukaryote-rich samples probably have
  more reads without a reference.
- **Batch:** on I09N between 2°S and 20°N (same water, about 30 °C both times) the
  share is 0.11 in 2016 (instrument OBI_2631) and 0.19 in 2025 (OBI_3310). Bacteria CLR
  is 0.4 lower in 2025 there; eukaryote CLR is the same. Median by instrument:
  OBI_2131 0.13, OBI_2631 0.17, OBI_3310 0.21. Year and instrument cannot be separated.
- **Granite does not see the batch.** Predicting the instrument (OBI_3310 vs older)
  with the cruise held out: pooled R2 -0.47 from the embedding, -0.37 from temperature
  (`data/results/taxa_temperature/instrument/`); the mean prediction per held-out cruise
  is 0.27-0.70 whatever the true instrument.
- **Predictability (pooled R2, best learner):**

| Target | Temperature | Time & place | Raw satellite | Granite | PCs + temperature |
|---|---:|---:|---:|---:|---:|
| log(unclassified / classified) | 0.26 | -0.28 | 0.20 | 0.19 | 0.22 |
| log(Bacteria / Eukaryota) | 0.52 | -0.43 | 0.59 | 0.56 | 0.57 |
| same, unclassified share regressed out | 0.17 | -0.05 | 0.30 | 0.36 | 0.28 |
| Bacteria CLR (4 domains) | 0.05 | 0.04 | 0.34 | 0.35 | 0.36 |
| Eukaryota CLR | 0.60 | -0.23 | 0.56 | 0.51 | 0.54 |

  The unclassified share explains 43% of the variance of log(Bacteria/Eukaryota).
  It is only weakly predictable, and mostly by temperature.
- **R22 holds up.** The bacteria/eukaryote balance is predicted about equally by
  temperature and Granite (0.52 vs 0.56), because its eukaryote side follows
  temperature. Once the unclassified share is regressed out, Granite leads clearly
  (0.36 vs 0.17); the Bacteria CLR, Granite 0.35 vs 0.05. The Granite-better KOs still
  track the adjusted ratio (median r +0.64; temperature-better -0.37).

## R24. Ocean-model fields instead of in-situ measurements (2026-11-12)

GLORYS12 and PISCES (Copernicus Marine reanalysis; forecast runs from 2021-11 for
fe and the carbonate variables) at each sample, from `cluster/model_fields.sbatch`
(`data/model_fields/`, not in git): 2,014 samples, every reanalysis variable present,
median 7-8 km to the cells used. Script: `scripts/model_field_tests.py`
(results in `data/results/model_fields/`).

**Model vs in-situ at the samples** (r on log10 for nutrients and chlorophyll):

| Variable | n | Cruises | r | Median r within a cruise | Median in situ | Median model |
|---|---:|---:|---:|---:|---:|---:|
| Temperature (GLORYS) | 1,827 | 14 | 0.93 | 0.98 | 25.85 | 25.76 |
| Salinity (GLORYS) | 1,734 | 14 | 0.96 | 0.92 | 35.33 | 35.35 |
| Phosphate (PISCES) | 728 | 13 | 0.70 | 0.56 | 0.120 | 0.015 |
| Nitrate (PISCES) | 429 | 6 | 0.68 | 0.33 | 0.020 | 0.010 |
| Chlorophyll (PISCES) | 178 | 1 | 0.57 | | 0.036 | 0.091 |

PISCES surface phosphate is far lower than measured in oligotrophic water (median
0.015 vs 0.12 umol/kg), but it ranks the samples reasonably.

How far off the model values are, at the same samples (quartiles 25 / 50 / 75%;
"model lower" = share of samples where the model value is below the measurement;
ratio = median of model / measured per sample):

| Variable | n | Measured q25 / median / q75 | Model q25 / median / q75 | Model lower | Median ratio | Median difference |
|---|---:|---|---|---:|---:|---:|
| Temperature (°C) | 1,827 | 20.64 / 25.85 / 28.99 | 21.06 / 25.76 / 28.97 | 58% | 1.00 | -0.07 °C |
| Salinity | 1,734 | 34.34 / 35.33 / 36.05 | 34.34 / 35.35 / 36.09 | 48% | 1.00 | +0.008 |
| Phosphate (umol/kg) | 728 | 0.040 / 0.120 / 0.171 | 0.001 / 0.015 / 0.229 | 67% | **0.13** (about 8x too low) | -0.033 |
| Nitrate + nitrite (umol/L) | 429 | 0.000 / 0.020 / 0.113 | 0.002 / 0.010 / 2.671 | 40% | 0.83 | +0.001 |
| Chlorophyll | 178 (1 cruise) | 0.035 / 0.036 / 0.041 | 0.081 / 0.091 / 0.139 | 0% | **2.47** (too high) | +0.052 |

- GLORYS physics is essentially unbiased (temperature within 0.1 °C, salinity within 0.01).
- PISCES phosphate is too low at the surface in two thirds of samples, about 8 times at
  the median. It is too low where phosphate is scarce and too high in the upper quartile.
- PISCES nitrate is about right at the median but far too high in the upper quartile
  (2.7 vs 0.11 umol/L): PISCES puts nitrate at the surface where little was measured.
- PISCES chlorophyll is 2.5 times the measured value on the one cruise with in-situ
  chlorophyll (FAIRe; the units of that column are not verified).
- Correlations (log10 for nutrients and chlorophyll): temperature 0.93 (0.98 within a
  cruise), salinity 0.96 (0.92), phosphate 0.70 (0.56), nitrate 0.68 (0.33),
  chlorophyll 0.57.

**Figure KOs, one cruise held out (median pooled R2, all 1,216 tiles):**

| Block | 32 PCs + in-situ T | 32 PCs + GLORYS T | PISCES P first, then 32 PCs |
|---|---:|---:|---:|
| Granite clearly better | 0.39 | 0.40 | 0.26 |
| Both high | 0.71 | 0.69 | 0.43 |
| Temperature clearly better | 0.46 | 0.45 | 0.15 |

- **GLORYS temperature can replace in-situ temperature** with no loss.
- **PISCES phosphate first fixes most of the pst problem, on all tiles:** pstC
  0.32 -> 0.59, pstB 0.41 -> 0.61, pstS 0.41 -> 0.64, pstA 0.38 -> 0.59. On P18,
  pstC goes from 0.55 to 0.85. On the 452 tiles with in-situ phosphate, PISCES gives
  0.49-0.60 against 0.64-0.71 with measured phosphate (baseline 0.24-0.38): about
  two thirds of the gain.
- As with measured phosphate (R20), phosphate first hurts the other KOs. Use it for
  the pst genes only.
- pstC still fails on I07N (-2.1) and I09N 2025 (-2.5) with PISCES. Both are Indian
  Ocean cruises with little pstC variation of their own.

Both inputs exist everywhere, so the pst genes can now be mapped globally, with a
caveat to state with any such map: **PISCES surface phosphate is too low** (median
0.015 vs 0.12 umol/kg measured; r 0.70, 0.56 within a cruise), and with it the model
recovers only about two thirds of the gain that measured phosphate gives.

## R25. Each GLORYS / PISCES variable with the current best Granite model (2026-11-12)

`scripts/model_field_first.py`: for each variable V, (a) V's spline first, then 32 PCs
+ temperature spline on the residuals, and (b) V's spline added to 32 PCs +
temperature spline in one ridge; against 32 PCs + temperature on the same tiles, figure
KOs, one cruise held out. Reanalysis variables on all 1,216 tiles; forecast-only ones
(iron, phytoplankton carbon, pH, DIC, alkalinity, pCO2) on the 637 tiles from 2021-11
on (fewer cruises, so noisier). Results in `data/results/model_fields/first_r2.csv`.

- **Putting a variable first hurts almost everything** (median change -0.01 to -0.81),
  because the first stage takes signal that the Granite PCs predict better. Adding it
  to the model costs little (median -0.04 to 0.00, DIC -0.14).
- **Phosphate is the only clear win:** pstC, pstS, pstA, pstB +0.20 to +0.30, first or
  added (added: pstC 0.32 -> 0.62). Nitrate does the same at half the size (+0.09 to
  +0.20), as it tracks phosphate.
- **Salinity matches the measurements best but adds almost nothing** (median 0.00
  added; pstB and pstC +0.11 to +0.12). Granite and temperature already carry it.
- Small, gene-specific gains (added, +0.05 to +0.17): chlorophyll for the pst genes
  and LACC1/udg; iron for citramalyl-CoA lyase (E4.1.3.25, +0.17 on 637 tiles);
  alkalinity for several Granite-block KOs (637 tiles). GLORYS temperature first helps
  four temperature-block KOs (TOP2, dgkA, NOP2, RRM2: +0.05 to +0.09).
- Mixed layer, current speed, sea level, silicate and oxygen help no gene reliably.

Per variable, median change in pooled R2 over the 54 figure KOs, and the KOs that gain
at least 0.05 (V first / V added):

| Variable | Tiles | First: median | First: KOs >= +0.05 | Added: median | Added: KOs >= +0.05 | Genes that gain (added unless noted) |
|---|---:|---:|---:|---:|---:|---|
| Salinity | 1,216 | -0.16 | 0 | 0.00 | 3 | pstB +0.12, pstC +0.11, pstS +0.06 |
| Mixed layer depth | 1,216 | -0.03 | 0 | -0.01 | 0 | - |
| Sea level | 1,216 | -0.23 | 0 | -0.04 | 3 | LACC1 +0.11, pntA +0.07, argG +0.06 |
| Current speed | 1,216 | -0.01 | 0 | 0.00 | 0 | - |
| GLORYS temperature | 1,216 | -0.03 | 4 | 0.00 | 0 | first: TOP2 +0.09, dgkA +0.07, NOP2 +0.06, RRM2 +0.05 |
| **Phosphate** | 1,216 | -0.19 | 4 | -0.03 | 4 | **pstC +0.30, pstS +0.26, pstA +0.24, pstB +0.22** (first: +0.20 to +0.29) |
| Nitrate | 1,216 | -0.13 | 4 | 0.00 | 4 | pstC +0.16, pstB +0.11, pstA +0.10, pstS +0.09 |
| Silicate | 1,216 | -0.11 | 0 | -0.01 | 1 | lysA +0.07 |
| Oxygen | 1,216 | -0.04 | 4 | -0.02 | 0 | first: E4.1.3.25 +0.07, dgkA +0.07, NOP2 +0.06, TOP2 +0.06 |
| Chlorophyll | 1,216 | -0.14 | 0 | 0.00 | 6 | pstA +0.08, LACC1 +0.08, pstC +0.07, udg +0.07, pstB +0.06, pstS +0.06 |
| Primary production | 1,216 | -0.10 | 0 | -0.01 | 2 | pstC +0.05, pstS +0.05 |
| Iron | 637 | -0.17 | 2 | -0.02 | 3 | E4.1.3.25 +0.17 (first +0.17), glsA +0.05 (first +0.16), ndhF +0.06 |
| Phytoplankton carbon | 637 | -0.12 | 1 | 0.00 | 5 | pstC +0.14, pstB +0.10, purK +0.08, lysA +0.08, pstS +0.06 |
| pH | 637 | -0.28 | 1 | -0.06 | 1 | NOP2 +0.10 |
| DIC | 637 | -0.81 | 0 | -0.14 | 5 | NOP2 +0.15, tldD +0.11, ndhF +0.07, glsA +0.06, dhaA +0.06 |
| Alkalinity | 637 | -0.32 | 3 | -0.02 | 8 | dhaA +0.17, pheA2 +0.15, glsA +0.13, tldD +0.09, E4.1.3.25 +0.08, udg +0.06, E3.1.2.28 +0.06, mntB +0.05 |
| pCO2 | 637 | -0.25 | 1 | -0.04 | 2 | NOP2 +0.09, E4.1.3.25 +0.07 |

The 637-tile rows (forecast runs, from 2021-11) have lower baselines and fewer cruises;
their gains are tentative. Outside the pst genes, no variable gives a consistent gain.

**Decision (2026-11-12):** R25 is a side note for the paper (one or two sentences, no
figure). These variables are not used in later models. The phosphate-first model for
the pst genes (R24) remains the only use of model fields, for global maps, with its
caveat. A caution belongs with the anecdote: the model fields are themselves off from
the measurements (PISCES phosphate about 8 times too low at the surface; within-cruise
r 0.56 for phosphate, 0.33 for nitrate; chlorophyll r 0.57 on the one cruise with
data). Part of why they add so little may be model error, not a lack of signal.

## R26. Can Granite predict the in-situ measurements? Against GLORYS / PISCES (2026-11-12)

`scripts/predict_insitu.py`: each in-situ variable as the target of
`test_ko_temperature.py` (one cruise held out, pooled R2, best of ridge / spline / knn
per input), and the GLORYS / PISCES value on the same tiles: as is, after a straight-line
calibration fitted on the training cruises (one cruise held out), and r2 (the ceiling of
any line). Nutrients on log10. In-situ chlorophyll comes from one cruise only, so it
cannot be tested this way. Results in `data/results/insitu/summary.csv`.

| Variable | Tiles | Cruises | Model as is | Model calibrated | Model r2 | Granite | 32 PCs + temperature | Raw satellite | Temperature | Time & place |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Temperature | 1,424 | 14 | 0.94 | 0.94 | 0.94 | 0.92 | - | 0.89 | - | 0.67 |
| Salinity | 1,397 | 14 | 0.94 | 0.93 | 0.94 | 0.23 | 0.10 | 0.09 | 0.21 | 0.21 |
| Phosphate | 492 | 13 | -0.47 | -0.04 | 0.40 | -0.24 | -0.09 | -0.03 | 0.03 | 0.28 |
| Nitrate | 244 | 6 | 0.30 | 0.12 | 0.51 | 0.01 | 0.02 | -0.04 | 0.09 | 0.01 |

- **Temperature:** Granite nearly matches GLORYS (0.92 vs 0.94). The OSTIA-filled tiles
  carry SST.
- **Salinity:** GLORYS by far (0.94 vs Granite 0.23).
- **Phosphate:** nothing predicts it across cruises. PISCES ranks samples (r2 0.40) but
  its offset differs between cruises, so even a calibrated line fails (-0.04). Position
  and season do best (0.28). Granite fails (-0.24).
- **Nitrate:** PISCES as is 0.30; Granite about 0 (6 cruises only).
- So Granite does not beat the ocean models on what they model, except that it nearly
  ties on temperature. Yet Granite predicts many KOs that these variables do not
  (R25). It sees biology more than water properties. PISCES phosphate helps the pst
  genes (R24) through its ranking, which the spline uses, not through its values.

**Conclusions from R24-R26 (for the manuscript, and for a physics-aware satellite model):**
1. GLORYS physics matches the ship measurements closely (temperature r 0.93, within
   0.1 °C; salinity r 0.96, within 0.01). PISCES biogeochemistry does not: surface
   phosphate is about 8 times too low at the median and only loosely right within a
   cruise (r 0.56). Nitrate is far too high where measured nitrate is low (upper
   quartile 2.7 vs 0.11 umol/L; within-cruise r 0.33). Chlorophyll is about 2.5 times too
   high on the one cruise that measured it.
2. Adding these fields to the Granite model helps almost nothing (R25). The one
   exception is phosphate for the four pst genes (+0.2 to +0.3 R2). It works through
   ranking (a spline), and recovers about two thirds of what measured phosphate gives
   (R20, R24).
3. Granite does not reproduce the water properties: salinity 0.23, phosphate -0.24,
   nitrate 0.01. It nearly matches GLORYS only on temperature (0.92 vs 0.94), which the
   OSTIA-filled tiles contain. Yet Granite predicts many KOs that the water properties
   do not.
4. So Granite's skill on the KOs does not come from standing in for temperature,
   salinity or nutrients. It comes from something the ocean models do not carry,
   consistent with R22-R23: community composition, above all the bacterial share.
5. For a satellite model made "physics-aware" with ocean-model fields: on this evidence
   the physical fields would add little for the biology targets here. The embedding
   already holds temperature, and salinity adds nothing for the KOs. The biogeochemical
   fields are too far from the measurements, except as a ranking for nutrient-specific
   genes. **Revised by R28 (all 3,921 KOs, not only the well-predicted figure KOs):**
   the nutrient fields add real, biologically coherent skill to nutrient-acquisition
   genes (the Pho regulon with phosphate; urea uptake and urease with nitrate). Physical
   fields and recent history help smaller, less consistent sets of genes. See R28.

## R27. Recent ocean history before sampling (2026-11-12)

Question: do summaries of the days *before* sampling help where the sampling-day model
values did not (R25)? Per sample (`scripts/extract_history_fields.py`,
`cluster/history_fields.sbatch`; output `data/model_fields/sample_history_fields.csv`,
`provenance_history.json`):

| Feature | Meaning | Source (dataset used) | Window | Area |
|---|---|---|---|---|
| mld_max30 | mixing: maximum mixed-layer depth | GLORYS12 mlotst (my, all 2,014) | 30 days up to sampling | +-0.125 deg mean |
| sst_change14 | cooling / upwelling: SST on the day minus the 14-day mean before | OSTIA L4 reprocessed (all 2,014) | 14 days | +-0.125 deg |
| wind_max14, wind_mean14 | storms: peak and mean wind speed | L4 hourly winds, multi-year (all 2,014) | 14 days | +-0.125 deg |
| ssh_gradient | fronts: mean sea-level gradient magnitude (m per 100 km) | GLORYS12 zos (my, all 2,014) | sampling day | +-0.25 deg |
| eke30 | eddies: mean 0.5 (u'^2 + v'^2) of geostrophic velocity anomalies (m2/s2) | DUACS L4: multi-year 1,851, NRT 163 (myint id does not exist) | 30 days | +-0.25 deg |

Every sample has every feature. Values (median, IQR from the table): 30-day max mixed
layer 37 m (now: 20 m; max 493 m); SST change +0.03 K (SD 0.47, max 2.7); 14-day peak
wind 10.9 m/s (mean 6.4); sea-level gradient 0.063 m per 100 km (max 1.17); EKE 0.009
m2/s2 (max 0.70).

**With the current best model (figure KOs, 1,216 tiles, one cruise held out;
`python scripts/model_field_first.py --history`, `data/results/model_fields/first_r2_history.csv`):**

| Feature | First: median change | First: best KO | Added: median change | Added: best KO | KOs gaining >= 0.05 |
|---|---:|---:|---:|---:|---:|
| Mixed layer, 30-day max | -0.134 | +0.027 | -0.019 | +0.047 (thiY) | 0 |
| SST change, 14 days | -0.012 | +0.030 | -0.002 | +0.015 (ispH) | 0 |
| Wind, 14-day peak | -0.078 | -0.007 | -0.007 | +0.039 (lysA) | 0 |
| Wind, 14-day mean | -0.085 | +0.010 | -0.015 | +0.007 | 0 |
| Sea-level gradient | -0.028 | +0.032 | 0.000 | +0.030 (proA) | 0 |
| Eddy kinetic energy, 30 days | -0.052 | +0.038 | -0.003 | +0.012 (enolase) | 0 |

**Are they already in Granite?** No. Pooled R2 predicting each feature (one cruise held
out, best learner; `data/results/model_fields/history_predictability/`):

| Feature | Granite | 32 PCs + temperature | Raw satellite | Temperature | Time & place |
|---|---:|---:|---:|---:|---:|
| log mixed layer, 30-day max | 0.10 | 0.05 | 0.10 | 0.20 | 0.01 |
| SST change, 14 days | -0.22 | -0.23 | -0.16 | -0.05 | -0.20 |
| Wind, 14-day peak | 0.10 | 0.11 | 0.19 | 0.23 | -0.23 |
| Wind, 14-day mean | 0.09 | 0.05 | 0.22 | 0.23 | -0.38 |
| log sea-level gradient | 0.18 | 0.20 | 0.16 | 0.19 | 0.24 |
| log EKE, 30 days | -0.07 | 0.03 | -0.01 | 0.03 | -0.03 |

- **No history feature helps any figure KO.** No gain reaches +0.05, first or added.
  The largest is +0.047 (30-day max mixed layer added, thiY). Putting one first costs
  0.01 to 0.13 at the median.
- **This is not because Granite already holds the information.** Granite predicts none
  of these features well (R2 -0.22 to 0.18). They are new information, and they still
  do not explain the KOs.
- **For the 54 figure KOs (R24-R27):** neither the sampling-day values nor the recent
  physical history add to Granite + temperature; the one exception is phosphate for the
  pst genes. **But the figure KOs were chosen for being well predicted** (the top of each
  heatmap block), which tilts the test against new inputs. On all KOs the answer
  changes: see R28.

## R28. Added variables on all KOs, not just the well-predicted figure KOs (2026-11-12)

Why: the 54 figure KOs are the top of each heatmap block (the 28 highest Granite R2
among "Granite clearly better", the 14 highest among "both high" (both >= 0.6), the 12
highest temperature R2 among "temperature clearly better"). They are well predicted by
construction, so R25/R27 were tilted against new inputs.

Test (`scripts/added_fields_all_kos.py`; `data/results/model_fields/added_all_kos.csv`,
`added_all_kos_summary.csv`, `added_useful_kos.csv`): every variable of R25 and R27 added as
a spline to 32 Granite PCs + temperature spline (one ridge, penalty by GCV per KO), all
3,921 KOs, one cruise held out, pooled R2. Controls for chance: three draws of random
noise and PISCES phosphate shuffled across tiles, on all 1,216 tiles and on the 637
forecast-period tiles (where iron, phytoplankton carbon, pH, DIC, alkalinity and pCO2
exist). Criteria per KO: gain >= 0.05; the same and improved in >= 75% of held-out
cruises; gain >= 0.05 **and** R2 >= 0.3 after ("now useful").

| Added variable | Tiles | KOs gain >= 0.05 | ... and in >= 75% of cruises | Now useful (gain >= 0.05, R2 >= 0.3) | Median gain, all KOs |
|---|---:|---:|---:|---:|---:|
| Salinity | 1,216 | 285 | 58 | 43 | -0.026 |
| Mixed layer (day) | 1,216 | 306 | 115 | 15 | -0.006 |
| Sea level | 1,216 | 450 | 45 | 39 | -0.029 |
| Current speed | 1,216 | 11 | 1 | 0 | -0.005 |
| GLORYS temperature | 1,216 | 35 | 16 | 1 | -0.002 |
| Phosphate | 1,216 | 343 | 101 | 70 | -0.028 |
| Nitrate | 1,216 | 401 | 105 | 67 | -0.010 |
| Silicate | 1,216 | 214 | 53 | 24 | -0.020 |
| Oxygen | 1,216 | 166 | 36 | 28 | -0.018 |
| Chlorophyll | 1,216 | 453 | 87 | 52 | -0.001 |
| Primary production | 1,216 | 199 | 14 | 10 | -0.017 |
| Mixed layer, 30-day max | 1,216 | 424 | 164 | 20 | -0.019 |
| SST change, 14 days | 1,216 | 20 | 7 | 0 | -0.005 |
| Wind, 14-day peak | 1,216 | 257 | 87 | 22 | -0.003 |
| Wind, 14-day mean | 1,216 | 103 | 49 | 4 | -0.009 |
| Sea-level gradient | 1,216 | 23 | 15 | 2 | -0.002 |
| Eddy kinetic energy, 30 days | 1,216 | 22 | 5 | 0 | -0.006 |
| **Control: noise (3 draws)** | 1,216 | 14-15 | 2-8 | **0** | -0.002 to 0.000 |
| **Control: shuffled phosphate** | 1,216 | 17 | 7 | **0** | -0.004 |
| Iron | 637 | 431 | 137 | 46 | -0.052 |
| Phytoplankton carbon | 637 | 600 | 233 | 29 | -0.006 |
| pH | 637 | 446 | 111 | 32 | -0.043 |
| DIC | 637 | 491 | 114 | 61 | -0.103 |
| Alkalinity | 637 | 709 | 197 | 121 | -0.050 |
| pCO2 | 637 | 411 | 69 | 24 | -0.049 |
| **Control: noise (3 draws)** | 637 | 211-259 | 48-60 | **0-3** | 0.000 to -0.002 |
| **Control: shuffled phosphate** | 637 | 210 | 40 | **1** | -0.002 |

Gains >= 0.05 by baseline R2 stratum (all tiles; KOs per stratum: below 0 1,058; 0-0.1
873; 0.1-0.3 1,404; >= 0.3 586): most gains are in KOs the baseline predicts poorly
(e.g. sea level 260 / 91 / 87 / 12; chlorophyll 178 / 118 / 133 / 24; phosphate
166 / 77 / 61 / 39; noise 12-15 / 2-3 / 0-1 / 0).

- **The figure-KO test was tilted.** On all KOs, every variable except current speed,
  GLORYS temperature, SST change, sea-level gradient and EKE helps far more KOs than
  chance (now useful: 15-70 per variable on all tiles vs 0 for the controls; 24-121
  on the 637 tiles vs 0-3).
- **Overall:** 414 distinct KOs (11% of 3,921) become useful with at least one
  variable (711 KO-variable pairs). 198 of them also improve in >= 75% of cruises (267
  pairs). They are mostly moderately predicted KOs (baseline median 0.29, IQR
  0.25-0.38). Cyanobacterial KOs are over-represented (62 of 145, 43%, vs 9-10% for
  Pseudomonadota and eukaryote KOs).
- **Median gain over all KOs is negative for every variable** (-0.001 to -0.10). Adding
  a variable to every KO hurts more than it helps; it is a per-gene choice.
- **Nutrients help nutrient-acquisition genes, coherently.**
  - Phosphate: phoU +0.38 (0.15 -> 0.52), phoR +0.30 (0.03 -> 0.33), pstC +0.30,
    pstS +0.26, pstA +0.24, adeR +0.31; phosphonate transport phnE +0.25 (with
    alkalinity).
  - Nitrate: urea uptake urtA +0.30 (0.26 -> 0.56), urtB +0.27, urtC +0.25, urtD
    +0.24 (phosphate); urease ureF +0.34, ureD +0.30, ureE +0.25, ureB +0.25.
    Cyanobacterial psaK +0.18, pcyA +0.13, menC +0.18.
  - The KO sets helped by phosphate and by nitrate overlap strongly (Jaccard 0.65); the
    PISCES nutrients co-vary, so this is "nutrient regime" more than one element. All
    other pairs overlap by 0.13 or less.
  - These gains come despite PISCES's biases (R24): phosphate about 8 times too low,
    nitrate far too high where little was measured. The spline uses the ranking.
- **Salinity helps phosphorus genes:** phnC +0.18, phnE +0.14, ppx-gppA +0.13, pstB
  +0.12, pstC +0.11 (each better in 9-11 of 13 cruises), and tartrate dehydrogenase
  ttuC +0.19. Plausibly salinity marks the saltier, phosphorus-poor North Atlantic.
- **Physical fields and recent history help small, mixed sets** (now useful: sea level
  39, 30-day max mixed layer 20, 14-day peak wind 22, mixed layer 15; median gain
  +0.06 to +0.09). Examples: 30-day max mixed layer for ynjD +0.18, mnhD +0.18, cpeD
  +0.16, phnE +0.15; sea level for aglG +0.41, tfdA +0.33. The sea-level hits are
  consistent in only 6 of 39; mixed layer 10 of 20; peak wind 16 of 22.
- **Carbonate system and iron (637 tiles, from 2021-11 only):** many hits (alkalinity
  121, DIC 61, iron 46), well above that subset's controls (0-3). Some top hits look
  like regional proxies rather than mechanism (DIC: tfdA +0.50, cycB +0.44). Tentative
  until more post-2021 cruises exist.

**What this changes for a physics-aware satellite model (supersedes the last point of
the R24-R26 conclusions):**
1. Nutrient fields (phosphate, nitrate), even biased, add real skill for
   nutrient-acquisition genes (Pho regulon, urea uptake and urease) and for some
   cyanobacterial genes, beyond Granite + temperature.
2. Physical fields and recent history (mixed layer, wind, sea level, salinity) help
   smaller, less consistent sets of genes. Salinity's clearest gains are again
   phosphorus genes.
3. Adding any field to all genes lowers the median R2. The value is in per-function
   use, so a physics-aware model would need to learn which functions use which
   field.

---

## R29. External test on the JGI marine metagenomes (Flinkstrom et al. 2024) (2026-11-16)

**Inputs.** Training: all Bio-GO-SHIP R17 units with tile SST (1,205 tiles, 13 cruises;
E-orig-l4 + E-faire-l4). Test: JGI tiles from `cluster/external_tiles.sbatch` +
`cluster/external_embed.sbatch` (jobs finished 2026-11-14; `data/external/jgi_flinkstrom/embeddings_sst_l4/`).
Same tile rule as training (fixed position, closest tile with SST, OSTIA-filled,
identical tiles merged): **90 of the 109 JGI samples got a tile, on 31 distinct tiles.**
Targets: the 2,832 shared KOs, CLR over that set in both datasets (R15).

| JGI group | Tiles | Tile SST, deg C | Dates |
|---|---:|---|---|
| Monterey Bay (one position) | 21 | 13.6-17.8 | 2016-09-26 to 11-16 |
| North Sea | 4 | 8.3-11.9 | 2016-04-26 to 05-17 |
| South Atlantic | 4 | 7.6-12.5 | 2018-06-01 to 2019-03-13 |
| Marineland FL | 1 | 27.4 | 2016-08-14 |
| Other (seawater/sediment study) | 1 | 25.8 | 2018-04-05 |

Tile SST range: Bio-GO-SHIP -1.7 to 32.9 (median 23.9), JGI 7.6 to 27.4 (median 15.3),
so JGI is inside the training range. Tile counts of the other sets (no KO tables yet):
Tara Pacific 85 samples with a tile (31 open ocean I00, the rest over 27 island groups),
EMO BON 295 samples (15 observatories).

**Settings.** No external set has in-situ temperature, so "temperature" is the tile's
own SST (SLSTR, or OSTIA where missing; Granite-normalised). Models, ridge with
per-KO GCV penalty, trained on all Bio-GO-SHIP tiles: Bio-GO-SHIP mean (no input);
granite (512-d embedding); pca32; sst (spline); pca32 + sst (R18 best model with tile
SST for in-situ temperature); raw satellite (16 OLCI band means + SST, rank-normalised,
spline); time+place (lat, sin/cos lon, sin/cos day of year). Metrics per KO on the JGI
tiles: pooled R2 against the JGI mean; **shift R2**, the R2 after removing the mean
error (one constant per KO fitted on JGI, which forgives a dataset offset and is
slightly optimistic); Pearson r. **Chance**: the expected number of KOs with shift
R2 >= 0.3 when the same predictions are given to permuted tiles (200 draws).

**Cruise-held-out reference on the shared KOs.** Re-centring the CLR over 2,832 KOs
lowers the CV skill: pca32 + in-situ temperature has 230 KOs with R2 >= 0.3 (median
0.043). The R28 baseline on the same 2,832 KOs (CLR over all 3,921) has 292 (median
0.082); per-KO r between the two = 0.92. The re-centring shift per sample has SD 0.12
(median KO SD 0.53).

| Model (CV, 1,205 tiles) | Median R2 | KOs R2 >= 0.3 |
|---|---:|---:|
| mean | -0.048 | 0 |
| granite | 0.029 | 201 |
| pca32 | 0.051 | 208 |
| sst | -0.009 | 128 |
| pca32 + sst | 0.038 | 203 |
| pca32 + in-situ temperature | 0.043 | 230 |
| raw satellite | -0.011 | 126 |
| time+place | -0.142 | 2 |

Using tile SST instead of in-situ temperature costs 27 KOs (230 to 203).

**Dataset offset.** JGI mean minus Bio-GO-SHIP mean per KO (CLR): median +0.008,
median absolute difference 0.450, more than 1 in 590 of 2,832 KOs. Within-dataset SD
per KO: JGI median 0.351, Bio-GO-SHIP 0.479. **The offset is larger than the spread
within JGI**, so plain R2 is strongly negative for every model, including the mean
(median -1.66).

**All 31 JGI tiles.** n = KOs in the stratum; "CV R2 >= 0.3" = KOs the best model
predicts well on Bio-GO-SHIP; "heatmap granite" = the strict Granite-better category
(84 of the 87 are shared).

| Model | KOs | n | CV median R2 | JGI median R2 | median shift R2 | median r | KOs shift R2 >= 0.3 | chance |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| mean | all | 2,832 | -0.048 | -1.664 | 0 | 0 | 0 | 0 |
| granite | all | 2,832 | 0.029 | -1.466 | -0.228 | 0.011 | 70 | 2.3 |
| granite | CV >= 0.3 | 230 | 0.370 | -1.922 | -0.243 | 0.125 | 10 | 0.3 |
| granite | heatmap granite | 84 | 0.244 | -1.116 | -0.055 | 0.022 | 1 | 0.02 |
| pca32 | all | 2,832 | 0.051 | -1.397 | -0.146 | 0.058 | 139 | 1.2 |
| pca32 | CV >= 0.3 | 230 | 0.362 | -2.122 | -0.179 | 0.230 | 24 | 0.2 |
| pca32 | heatmap granite | 84 | 0.238 | -0.977 | 0.068 | 0.283 | 5 | 0.01 |
| sst | all | 2,832 | -0.009 | -1.210 | -0.018 | 0.113 | 88 | 1.6 |
| sst | CV >= 0.3 | 230 | 0.298 | -2.210 | -0.143 | 0.186 | 12 | 0.5 |
| sst | heatmap granite | 84 | 0.130 | -1.138 | 0.034 | 0.243 | 4 | 0.04 |
| pca32 + sst | all | 2,832 | 0.038 | -1.442 | -0.186 | 0.044 | 123 | 1.5 |
| pca32 + sst | CV >= 0.3 | 230 | 0.367 | -1.937 | -0.207 | 0.223 | 23 | 0.4 |
| pca32 + sst | heatmap granite | 84 | 0.231 | -0.861 | 0.059 | 0.266 | 4 | 0.01 |
| raw satellite | all | 2,832 | -0.011 | -1.784 | -0.510 | 0.053 | 36 | 1.2 |
| raw satellite | CV >= 0.3 | 230 | 0.303 | -2.246 | -0.683 | 0.090 | 6 | 0.2 |
| time+place | all | 2,832 | -0.142 | -1.521 | -0.182 | -0.099 | 88 | 0.9 |
| time+place | CV >= 0.3 | 230 | -0.189 | -2.618 | -0.438 | -0.153 | 4 | 0.07 |

**Monterey Bay alone (21 tiles, one position, 7 weeks: change over time).**

| Model | KOs | median shift R2 | median r | KOs shift R2 >= 0.3 | chance |
|---|---|---:|---:|---:|---:|
| granite | all | -0.354 | 0.006 | 14 | 7.4 |
| pca32 | all | -0.221 | 0.029 | 31 | 7.6 |
| pca32 | CV >= 0.3 | -0.231 | 0.027 | 4 | 0.5 |
| sst | all | -0.056 | -0.036 | 0 | 2.8 |
| pca32 + sst | all | -0.282 | 0.018 | 33 | 6.7 |
| pca32 + sst | CV >= 0.3 | -0.364 | 0.011 | 4 | 0.4 |
| raw satellite | all | -1.365 | 0.032 | 38 | 3.0 |
| time+place | all | -0.016 | -0.030 | 63 | 1.6 |

**The other 10 JGI tiles (North Sea, South Atlantic, Florida, one other: a spatial contrast).**

| Model | KOs | median shift R2 | median r | KOs shift R2 >= 0.3 | chance |
|---|---|---:|---:|---:|---:|
| granite | all | -0.171 | 0.047 | 214 | 75.1 |
| granite | CV >= 0.3 | -0.129 | 0.252 | 29 | 7.3 |
| pca32 | all | -0.050 | 0.124 | 296 | 77.5 |
| pca32 | CV >= 0.3 | 0.041 | 0.387 | 60 | 8.7 |
| pca32 | heatmap granite | 0.065 | 0.291 | 9 | 1.0 |
| sst | all | -0.016 | 0.152 | 275 | 70.0 |
| sst | CV >= 0.3 | -0.168 | 0.289 | 38 | 9.3 |
| pca32 + sst | all | -0.083 | 0.106 | 331 | 80.2 |
| pca32 + sst | CV >= 0.3 | -0.043 | 0.383 | 59 | 8.3 |
| pca32 + sst | heatmap granite | 0.066 | 0.305 | 8 | 1.1 |
| raw satellite | all | -0.266 | 0.094 | 184 | 80.0 |
| time+place | all | -0.160 | -0.074 | 131 | 54.5 |

**KOs that transfer best** (pca32 + sst, among the 230 with CV R2 >= 0.3, ranked by
shift R2 on all 31 tiles):

| KO | Gene | Function | CV R2 | Shift R2 (all) | r (all) | Shift R2 (other 10) | r (Monterey) |
|---|---|---|---:|---:|---:|---:|---:|
| K01595 | ppc | phosphoenolpyruvate carboxylase | 0.30 | 0.59 | 0.77 | 0.33 | 0.17 |
| K03147 | thiC | phosphomethylpyrimidine synthase | 0.42 | 0.53 | 0.73 | 0.46 | 0.50 |
| K00005 | gldA | glycerol dehydrogenase | 0.40 | 0.53 | 0.80 | 0.45 | 0.15 |
| K03154 | thiS | sulfur carrier protein | 0.26 | 0.51 | 0.73 | 0.56 | 0.35 |
| K09920 | - | uncharacterized protein | 0.42 | 0.49 | 0.71 | 0.56 | 0.39 |
| K05577 | ndhF | NAD(P)H-quinone oxidoreductase subunit 5 | 0.49 | 0.48 | 0.78 | 0.34 | 0.12 |
| K00847 | E2.7.1.4 | fructokinase | 0.34 | 0.45 | 0.69 | 0.33 | 0.29 |
| K04518 | pheA2 | prephenate dehydratase | 0.39 | 0.44 | 0.71 | 0.30 | -0.06 |
| K03149 | thiG | thiazole synthase | 0.33 | 0.43 | 0.67 | 0.42 | 0.58 |
| K04097 | HPGDS | prostaglandin-H2 D-isomerase / GST | 0.40 | 0.41 | 0.73 | 0.55 | -0.12 |
| K07729 | - | putative transcriptional regulator | 0.36 | 0.39 | 0.63 | 0.62 | 0.14 |
| K02441 | glpG | rhomboid protease GlpG | 0.32 | 0.38 | 0.74 | 0.55 | 0.12 |

(The CV R2 column is the pca32 + sst model; the stratum was chosen on pca32 + in-situ
temperature.) Three thiamine-synthesis genes (thiC, thiS, thiG) transfer, and thiG
and thiC also follow the Monterey Bay time series (r 0.58, 0.50).

**Conclusions.**
1. **Absolute KO levels do not transfer.** The per-KO offset between JGI and
   Bio-GO-SHIP (median |difference| 0.45 CLR) exceeds the spread within JGI (SD 0.35),
   so every model, the mean included, has strongly negative R2. The offset is a
   pipeline difference as much as biology (JGI IMG annotation, older KEGG release,
   different assembly and gene counting; R13, R15), so a new dataset needs our KO
   pipeline, or at least a per-KO calibration.
2. **With the offset forgiven, transfer is weak but above chance.** On all 31 tiles,
   pca32 gives 139 KOs with shift R2 >= 0.3 against 1.2 by chance (pca32 + sst 123
   vs 1.5). Of the 230 KOs the model predicts well on Bio-GO-SHIP, 23-24 (10%)
   transfer this way; their median r is 0.22-0.23.
3. **No skill over time at one coastal site.** Monterey Bay (21 tiles over 7 weeks)
   is 68% of the test, and no model tracks its changes: median r 0.01-0.03,
   counts near chance. Bio-GO-SHIP trains on spatial contrast between open-ocean
   tiles; this is a different kind of variation (coastal, short time scale).
4. **Some spatial transfer across the 10 other tiles.** For KOs the model predicts
   well, median r 0.38-0.39 (pca32, pca32 + sst) and 59-60 KOs at shift R2 >= 0.3
   against 8-9 by chance. With 10 tiles in 3 regions this is thin evidence.
5. **Tile SST alone transfers about as well as Granite** on JGI (all tiles: sst 88
   KOs, pca32 + sst 123). This test cannot show Granite's advantage over
   temperature.
6. **JGI is a weak test set.** It has 31 tiles, mostly one coastal site, and a
   different annotation pipeline. Tara Pacific (85 samples with tiles, 31 open
   ocean, 27 island groups, 2016-2018) is much closer to Bio-GO-SHIP and is the
   better external test, but it needs our KO pipeline run on the raw reads (R14,
   R15). EMO BON (295 samples, coastal Europe) would test coastal transfer.

**Files.** `data/results/external_jgi/`: `external_kos.csv` (KO x model: CV R2,
R2, shift R2, r and chance probability for all / Monterey Bay / other JGI; JGI-minus-
Bio-GO-SHIP mean in `bias` on the mean rows; heatmap category), `summary.csv`,
`units.csv` (the 31 JGI tiles with tile SST in deg C).

**Command.** `python scripts/external_test_jgi.py` (about 2 minutes).

---

## R30. The Bio-GO-SHIP KO pipeline (biogoship/biogoship-metagenomics) (2026-11-17)

Source: github.com/biogoship/biogoship-metagenomics, commit 402e528 (cloned read-only
to `../biogoship-metagenomics`; not modified). It is a **read-mapping** pipeline, not
per-sample assembly. This corrects the guesses in R15.

Steps (`mapping_pipeline/Snakefile`, Snakemake 7):
1. Trimmomatic 0.39 PE: `ILLUMINACLIP:adapters.fa:2:30:7 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:50`;
   only the paired reads are mapped.
2. minimap2 2.28 `-ax sr` against a non-redundant nucleotide reference
   (`biogoship_2026_database_v2`: mmseqs cluster representatives of several source
   databases, e.g. OMDB v2 NR95; built with `clustering_pipeline/`), with a pre-built
   `.mmi` index (`-I 32g` mapping memory; the README says indexing needed a 2,000 GB
   bigmem node at WHOI).
3. SAM header filtered to the sequences with hits, sorted BAM, `samtools coverage`.
4. `count_biogoship_taxa.py`: reads per taxon at each rank, from a taxonomy table
   (sequence, method, domain ... species).
5. `count_biogoship_annotations.py`: per reference sequence, **numreads from
   samtools coverage**, split 1/n over its n eggNOG-mapper KOs (column 12); also Pfam
   (column 21) and the root OG (column 5). Lengths are the reference sequence lengths
   summed. Output per rank: `level_<i>_<rank>_{ko,pfam,og}_counts.tsv`.

So the counts are **mapped reads per KO** (not gene counts as R15 assumed), and
**cov50 is `--coverage_threshold 50`**: only reference sequences with >= 50% of their
bases covered (samtools coverage column 6) count. The Snakefile's default is 0
("coverage0"). This fits the R15 observation that the kept share falls with gene
length: a long reference needs more reads to reach 50% covered.

Inputs needed to run it on other data (all at WHOI, `/proj/omics/alexander/...`):
the reference fasta or its `.mmi` index, `biogoship_2026_taxonomy.tsv`,
`biogoship_eggnogmapper.tsv`. None are in the repository.

---

## R31. Raw reads for the external test sets; OMDB v2.0 as interim reference (2026-11-17)

**Read lists** (`scripts/build_read_manifests.py`, run locally 2026-11-17 against the ENA
portal API; `data/external/<set>/reads/runs.csv`, `unresolved.csv`). Samples: those kept
in the manifests and with a Sentinel-3 tile.

| Set | Samples with a tile | With reads | Runs | Size (fastq.gz) | Bases | Instruments | How found |
|---|---:|---:|---:|---:|---:|---|---|
| Tara Pacific | 85 | 84 | 84 | 1.74 TB | 3.57 Tbp | NovaSeq 6000 78, HiSeq 4000 6 | run accessions in samples.csv |
| JGI Flinkstrom | 90 | 90 | 90 | 0.53 TB | 0.87 Tbp | HiSeq 2500 85, NovaSeq 6000 5 | ENA sample_title = JGI sample name |
| EMO BON | 295 | 131 | 131 | 1.19 TB | 2.48 Tbp | NovaSeq 6000 131 | biosample accession |

- All runs are paired, two files per run, one run per sample.
- Tara Pacific: TARA_SAMEA6040045 (run ERR11733808, 52.6 Gbp) has no fastq on ENA
  (strategy OTHER, no fastq_ftp); left out.
- Tara Pacific runs are labelled library_strategy OTHER (source METAGENOMIC); the
  filter keeps WGS or OTHER with a metagenomic source.
- JGI: the IMG genome id is not on the SRA run, but the sample name is the ENA
  sample_title (e.g. Helgoland_Spring_Bloom_20160426_1 -> SRR7058375, 36.3 Gbp). A
  fallback via the NCBI BioSample carrying IMG_genome_ID (often a MAG record, which
  names the read experiment) was not needed.
- **EMO BON: 164 of 295 samples have no public reads.** They have no biosample
  accession in the logsheets (160 from 2022, 4 from 2021), and none of their
  observatory + date combinations appears among the 192 EMO BON water runs on ENA
  (search sample_title "EMOBON_*", 2026-11-17). The 131 with reads are all 2021.
  ENA titles look like `EMOBON_AAOT_Wa_210622_3um_1`.

**Interim reference: OMDB v2.0 NR95** (Ocean Microbiomics Database, Sunagawa lab,
CC BY 4.0; https://sunagawalab.ethz.ch/share/microbiomics/ocean/db/2.0/data/catalogs/).
- `OMDBv2.0_NT_G_NR95.fna.gz` (27 GB; 103,044,829 genes, mmseqs clusters at 95%
  identity, 90% coverage of the shorter, from 508.8 M redundant genes; the published
  MD5 41c196d8... is of the uncompressed fasta). `OMDBv2.0_NT_G_NR95.cluster.tsv.gz`
  (3.4 GB: CLUSTER, LENGTH, SIZE, REPRESENTATIVE, MEMBERS). `OMDBv2.0_data.tsv.gz`
  (per genome: sample, study and file URLs).
- **The Bio-GO-SHIP reference contains OMDB NR95 verbatim**: the first OMDB gene
  (`OMDBv2.0_NT_G_NR95_000000000001`, 189 bp) is the example in their
  `filter_sam_header.py` (`LN:189`). Their final database adds other sources
  (eukaryotes among them) and re-clusters at 90% identity / 90% coverage
  (`clustering_pipeline/`, mmseqs 16).
- **No KO annotation and no taxonomy are published for the OMDB catalog** (none in
  catalogs/, databases/ or genomes/). Representatives are named after their genome
  (`<study>_<sample>_MAG_<n>-scaffold_<k>_<gene>`), so gene taxonomy can come from a
  genome table. KOs need eggNOG-mapper run ourselves, or the Bio-GO-SHIP
  `biogoship_eggnogmapper.tsv`, which should cover these ids.
- OMDB is prokaryotic genomes only (MAGs, SAGs, isolates); our targets also hold
  eukaryotic KOs. A test on OMDB-mapped external data must use Bio-GO-SHIP libraries
  mapped the same way.

**Download job:** `cluster/download_reads.sbatch` (service partition, 48 h, resumable,
MD5 per file).

**Update 2026-11-17: the Bio-GO-SHIP reference itself is available** (team Google Drive
folder 1c4wPBOLKLz7ZIIccmnKc8Rfz08r7vJXB): `biogoship_2026_database_v2.mmi` (the index
config.yaml names), a `biogoship_2026_database_v2_32G.mmi` variant (probably the same
index built in larger parts, `-I 32G`; not needed), `biogoship_2026_taxonomy.tsv`,
`biogoship_eggnogmapper.tsv`. These go straight to Morrill with
`cluster/download_reference.sbatch` (rclone). OMDB is no longer needed; the reads job
skips it unless `OMDB=1`. The three files together are at least 74 GB (partial
browser downloads stopped at 26.3, 24.5 and 22.8 GB).

## R32. The Bio-GO-SHIP reference on Morrill; mapping jobs for the external reads (2026-11-18)

**Reference downloaded** (`cluster/download_reference.sbatch`, job 88365, 2026-11-17,
4 h 3 min, rclone from the team Drive folder, MD5 checked) to
`data/reference/biogoship_2026/` on Morrill. The earlier estimate of "at least 74 GB"
(R31) was the size of stopped browser downloads; the files are 605 GB.

| File | Bytes | Drive date | Use |
|---|---:|---|---|
| `biogoship_2026_database_v2.mmi` | 470,496,221,205 | 2026-05-01 | minimap2 index (config.yaml `input_database`) |
| `biogoship_2026_taxonomy.tsv` | 29,266,349,170 | 2026-06-24 | reference sequence -> method, domain ... species |
| `biogoship_eggnogmapper.tsv` | 105,721,874,584 | 2026-06-30 | eggNOG-mapper annotations |

Not downloaded (also in the folder): `biogoship_2026_database_v2.fasta` (212.0 GB),
`biogoship_2026_database_v2.tsv` (29.5 GB), `biogoship_2026_database_v2_32G.mmi`
(427.7 GB), and the first version (`biogoship_2026_database.{fasta,mmi,tsv}`, 2026-04-13).

Headers, checked against their scripts:
- Taxonomy: `transcript, classification_method, domain, phylum_division, class, order,
  family, genus, species` (as `count_biogoship_*.py` expect). Sequence ids look like
  `A16;unbinned;A16_177843044` (assembly; bin; gene). Methods seen: EUKulele.
- eggNOG: the standard eggNOG-mapper columns, `#query` first. `KEGG_ko` is column 12 and
  `PFAMs` column 21, the indices 11 and 20 their counting script reads; `eggNOG_OGs` is
  column 5. KOs are written `ko:K03320`.

**How the external reads are processed.** We do not run their Snakefile. Instead, we
use the same commands and their scripts unchanged (copy at 402e528, read-only), from
three jobs:
- `setup_mapping_env.sbatch`: builds the tools at their versions (minimap2 2.28,
  samtools 1.18, trimmomatic 0.39).
- `map_reads.sbatch`: one array task per run. Steps, with their parameters:
  1. trimmomatic PE with `ILLUMINACLIP:adapters.fa:2:30:7 LEADING:2 TRAILING:2
     SLIDINGWINDOW:4:2 MINLEN:50`; unpaired survivors are dropped, as in theirs.
  2. `minimap2 -I32g -ax sr v2.mmi --split-prefix`.
  3. `filter_sam_header.py`.
  4. samtools sort, then samtools coverage and flagstat.

  Where it departs from theirs:
  - Thread counts.
  - The trimmed pairs are mapped in chunks of 50 M pairs and the chunk BAMs merged
    before samtools coverage. Each chunk is a checkpoint, so runs larger than one 48 h
    job can finish. minimap2 maps each pair on its own, so the chunking does not change
    alignments.
- `count_reads.sbatch`: runs `count_biogoship_annotations.py` and `count_biogoship_taxa.py`
  at coverage threshold 0 and 50. Coverage 50 matches our Bio-GO-SHIP targets; coverage
  0 is their default.

**Why the training-table total equals their output.** The training table
(`ko_counts_cov50_all.tsv`) is per library: the KO counts summed over taxa (R21).
Their script counts only reference sequences that have a taxonomy row. The domain
level (`level_1_domain_ko_counts`) summed over taxa gives the same total.

**Size of the job** (base counts from ENA; trimming and mapping times still unknown):

| Set | Runs | Median Gbp | Largest Gbp | Total Tbp |
|---|---:|---:|---:|---:|
| Tara Pacific | 84 | 39.7 | 95.6 | 3.57 |
| JGI | 90 | 7.9 | 43.1 | 0.87 |
| EMO BON | 131 | 18.7 | 26.3 | 2.48 |

**Plan:**
1. Run a pilot on the smallest run first. It gives the minimap2 time per index part
   and the number of parts. The index is 470 GB, likely in ~8 G-base parts, and each
   chunk reads it all once.
2. Use the pilot to set the chunk size and memory, then run all 305.
3. Collect cov50 KO tables per set in the layout of `ko_counts_cov50_all.tsv`, build CLR
   targets, and redo the external test (R29) for Tara Pacific, JGI and EMO BON.


## R33. What the satellite misses: shared structure in the residuals (2026-11-19)

**Question.** After the best Granite model, is what remains unexplained structured?
If KOs share residual axes, those axes are functional signals that no surface input
captures. They could reflect subsurface drivers or batch effects.

**Method.** Script `scripts/residual_modules.py`; outputs in `data/results/residual_modules/`.
A copy of both is in `transfer/`.

To reproduce, run from the repository root:

    OMP_NUM_THREADS=4 ~/miniforge3/bin/python scripts/residual_modules.py --targets data/targets/ko_clr_incl_protocolA.csv

- Since R36, `data/targets/ko_clr.csv` leaves out the A16N Protocol A samples. R33 was
  run on the earlier targets, kept as `ko_clr_incl_protocolA.csv`; hence the `--targets`.
- The defaults are `--components 10 --null-draws 20 --splits 20 --seed 0`.
- Software: Python 3.9.10, numpy 1.21.5, pandas 1.4.3, scikit-learn 1.6.1.
- All random steps (nulls, splits) use the seed, and the PCA and ridge fits are deterministic.
- The script's docstring still describes the permutation null as the selection rule.
  That rule was dropped (see step 5). The script reports it but keeps the first 10
  components either way.

*1. Inputs.*
- **Targets.** `data/targets/ko_clr.csv`, built by `scripts/build_ko_targets.py` from
  `ko_counts_cov50_all.tsv`, with these filters:
  - libraries with at least 1e5 KO counts, from at most 10 m depth
  - library proportions averaged per sample
  - KOs present in at least 90% of samples, giving 3,921 KOs
  - CLR = log(p + 0.5 / total) minus the sample's mean log
  The pseudocount depends on the total. This is why depth reaches the CLR (axis 3).
- **Embeddings.** Granite, from `data/embeddings_sst_l4/` and
  `data/faire/embeddings_sst_l4/`.
  - One chosen tile per sample: position "fixed", best tile "prefer-sst".
  - Only samples with in-situ temperature are kept.
- **Tiles.**
  - Samples with identical embeddings (rounded to 5 decimals) form one tile.
  - Their CLRs and temperatures are averaged.
  - The result is 1,216 tiles from 13 cruises (`tune_figure_kos.load`).

*2. Baseline* (`added_fields_all_kos.cv`). One cruise is held out per fold, and
everything below is fitted on the other 12:
- embeddings standardised, then PCA, keeping 32 components (random_state 0); the
  components are standardised again
- in-situ temperature standardised, then expanded with a cubic B-spline (8 knots,
  linear extrapolation)
- RidgeCV on the concatenated features, alphas = logspace(−1, 5, 25), one alpha per
  KO
- predictions made for the held-out cruise

Residual R = Y − Ŷ (1,216 × 3,921). KO R2 is pooled over all held-out predictions;
its median is 0.103.

*3. Variance split.*
- Residual share = Σ(R − mean R)² / Σ(Y − mean Y)², with sums over tiles and KOs.
- Between-cruise share = Σ(cruise mean of R − mean R)² / Σ(R − mean R)².

*4. PCA.*
- Two input matrices, with each KO column centred and divided by its SD:
  - "full" (F): R itself
  - "within cruise" (W): R minus its per-cruise KO means, then standardised
- Eigen-decomposition of the tile Gram matrix ZZᵀ, which is exact and faster than
  3,921².
- KO-space vectors v = Zᵀu / √λ.
- Share of residual variance = λ / Σλ.
- Score s = Zv / √λ, with unit sum of squares.
- Loading = v·√(λ / (n−1)): the correlation of a KO's standardised residual with the
  score.
- Sign: chosen so that the loadings sum to positive.

*5. Which components are real.*
- **Permutation null** (reported, not used). 20 draws; each KO's residuals are
  permuted across tiles on their own. A component passes if its eigenvalue beats the
  95th percentile at its rank. It breaks the KO–KO covariance, so any shared
  structure passes.
- **Split-half reproducibility** (used).
  - Each of 20 splits divides the 13 cruises at random into 6 and 7.
  - Each half's matrix (F or W) is re-standardised, and its first 30 vectors are
    computed.
  - For each full-data component, the best match in each half is the vector with the
    largest |cosine|.
  - Congruence = |cosine| between the two matches.
  - Reported: the median over the 20 splits, and the share of splits with congruence
    ≥ 0.85 (Lorenzo-Seva and ten Berge's "fair similarity" bar).

*6. Candidate drivers.*
- **Values per tile.** Per sample, then averaged over the samples in a tile.
- **Sources:**
  - Sampling time and position: `data/model_fields/sample_model_fields.csv`.
    - Year as a decimal.
    - Season as cos/sin of (day of year − 172), multiplied by the hemisphere sign.
    - Local solar hour = UTC hour + lon/15, entered as cos and sin.
  - Depth: `ko_samples.csv`.
  - In situ (manifest): temperature, log10 phosphate (floor 0.005), log10
    nitrate+nitrite (floor 0.01).
  - CCHDO bottles (`biogoship_samples_cchdo.csv`): salinity, oxygen, log silicate,
    log nitrite, DIC, alkalinity, pH.
  - Model fields: GLORYS and PISCES day fields, plus the R27 history fields, with the
    log floors of `model_field_first.VARIABLES` / `HISTORY`.
  - Sequencing:
    - log10 total KO counts
    - number of libraries
    - instrument OBI_3310 share
    - share from 2021 on
  - Taxa: domain CLRs, the unclassified-read logit, and all phylum CLRs.

*7. Description statistics.*
- **Cruise eta2** = var(per-cruise mean of the score, broadcast to tiles) / var(score).
- **Spearman correlations.**
  - "All" uses ranks over all tiles.
  - "Within" ranks over all tiles, then demeans the ranks per cruise for both
    variables, then takes the Pearson r.
  - A driver is used only with ≥ 100 non-missing tiles in ≥ 3 cruises.
  - Sign consistency: among cruises with ≥ 15 tiles and a non-constant driver, the
    number whose own Spearman r has the sign of the pooled within r.
  - The table lists the 3 phyla with the largest |within r|.
- **Group CV R2.**
  - Leave one cruise out.
  - Missing values are filled with training medians.
  - Standardised, then RidgeCV with alphas = logspace(−2, 4, 25).
  - Pooled R2 = 1 − SSE / SST over all tiles. Negative means worse than predicting
    the mean.
  - For W components the predictors are centred per cruise first.
  - Groups:
    - sequencing: all "sequencing:" columns
    - sampling: all "sampling:" columns
    - GLORYS+history: model columns with > 90% tile coverage
    - all environment: sequencing + sampling + GLORYS+history
    - taxa: domain CLRs, unclassified logit, phylum CLRs
    - PISCES comparison: on the 637 tiles where every model column exists,
      GLORYS+history with and without the PISCES fields
- **KEGG.**
  - `data/reference/kegg_brite_ko00001.csv`, restricted to the four pathway A-classes:
    Metabolism, Genetic and Environmental Information Processing, Cellular Processes.
  - Categories at levels B and C with ≥ 10 of the target KOs.
  - z = (mean loading in the category − mean over all KOs) / (SD over all KOs / √k).

*8. Outputs.*
- `components.csv`: share, null share, split-half congruence, eta2, group R2.
- `drivers.csv`: every driver × component.
- `categories.csv`: KEGG B and C z per component.
- `top_kos.csv`: the 10 KOs at each pole.
- `scores.csv`: tile × component, with section and sample_id.
- `loadings.csv.gz`: KO × component.
- `summary.txt`: the console log.

**Residual variance.** 83.7% of the KO CLR variance is left after the baseline.
- 20.1% of it lies between cruises: per-cruise offsets, which a model with the cruise
  held out cannot learn.
- 79.9% lies within cruises.

**Components** (share = of the standardised residual variance):

| Component | Share of residual variance | Split-half congruence | Splits >= 0.85 | Cruise eta2 | CV R2: sequencing | CV R2: sampling | CV R2: GLORYS+history | CV R2: all environment | CV R2: taxa (carriers) | CV R2: GLORYS+history (PISCES tiles) | CV R2: + PISCES (PISCES tiles) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| F1 | 19.4% | 0.90 | 85% | 0.36 | -0.02 | +0.03 | -0.37 | -0.16 | +0.56 | -0.14 | -0.33 |
| F2 | 9.9% | 0.77 | 10% | 0.24 | +0.10 | -0.18 | -0.27 | -0.07 | +0.10 | -0.44 | -0.09 |
| F3 | 6.8% | 0.65 | 0% | 0.50 | +0.32 | +0.05 | -0.38 | +0.10 | +0.10 | -0.17 | -0.41 |
| F4 | 5.1% | 0.58 | 0% | 0.47 | -0.35 | -0.88 | -0.86 | -0.67 | -0.62 | -0.36 | -0.28 |
| F5 | 4.7% | 0.54 | 0% | 0.19 | +0.08 | -0.11 | -0.19 | -0.14 | +0.42 | -0.43 | -0.63 |
| F6 | 3.0% | 0.42 | 0% | 0.48 | -0.10 | -0.42 | -0.61 | -0.37 | +0.14 | +0.06 | -0.05 |
| F7 | 2.7% | 0.37 | 0% | 0.16 | +0.13 | -0.25 | -0.19 | -0.08 | -0.17 | -0.36 | -0.42 |
| F8 | 2.1% | 0.31 | 0% | 0.25 | -0.07 | -0.38 | -0.46 | -0.86 | +0.01 | -0.35 | -0.72 |
| F9 | 1.8% | 0.25 | 0% | 0.37 | -0.30 | -0.40 | -0.19 | -0.32 | -0.19 | -0.34 | -0.43 |
| F10 | 1.6% | 0.30 | 0% | 0.33 | -0.17 | -0.51 | -0.33 | -1.53 | +0.04 | -0.69 | -0.72 |
| W1 | 17.2% | 0.95 | 100% | 0.00 | +0.00 | -0.06 | -0.05 | -0.04 | +0.46 | +0.03 | +0.01 |
| W2 | 9.8% | 0.83 | 25% | 0.00 | +0.00 | -0.14 | -0.06 | -0.06 | +0.52 | -0.37 | -0.42 |
| W3 | 5.8% | 0.66 | 5% | 0.00 | +0.28 | -0.02 | +0.01 | +0.18 | +0.40 | -0.03 | -0.10 |
| W4 | 4.7% | 0.56 | 0% | 0.00 | +0.00 | -0.10 | -0.08 | -0.25 | +0.49 | -0.19 | -0.29 |
| W5 | 3.3% | 0.41 | 0% | 0.00 | +0.21 | -0.01 | -0.13 | +0.14 | +0.28 | -0.22 | -0.35 |
| W6 | 2.6% | 0.37 | 0% | 0.00 | +0.01 | -0.04 | -0.01 | -0.16 | -0.09 | -0.26 | -0.26 |
| W7 | 2.0% | 0.23 | 0% | 0.00 | -0.01 | -0.13 | -0.18 | +0.03 | -0.32 | +0.02 | +0.08 |
| W8 | 1.9% | 0.33 | 0% | 0.00 | +0.04 | -0.28 | -0.26 | -0.45 | +0.01 | -0.04 | -0.21 |
| W9 | 1.6% | 0.27 | 0% | 0.00 | +0.01 | -0.05 | -0.09 | -0.20 | +0.08 | -0.07 | -0.26 |
| W10 | 1.5% | 0.17 | 0% | 0.00 | -0.00 | -0.13 | -0.23 | -0.40 | +0.01 | -0.34 | -0.40 |

**Driver correlations** for the stable and near-stable components. Each cell shows the
all-tile r, then the within-cruise r, then in how many cruises the within-cruise sign
holds. "—" means the variable is constant within cruises. All driver correlations are
in `drivers.csv`.

| Driver (n tiles) | F1 | F2 | F3 | F5 | W1 | W2 | W3 | W4 |
|---|---|---|---|---|---|---|---|---|
| taxa: Bacteria (CLR) (1216) | +0.53 / +0.45 (13/13) | -0.10 / +0.08 (8/13) | -0.15 / -0.33 (11/13) | +0.04 / +0.20 (7/13) | +0.39 / +0.46 (13/13) | +0.00 / +0.01 (6/13) | -0.20 / -0.22 (10/13) | +0.20 / +0.28 (10/13) |
| taxa: Eukaryota (CLR) (1216) | -0.30 / -0.15 (9/13) | -0.03 / -0.15 (9/13) | -0.10 / +0.09 (10/13) | +0.48 / +0.47 (11/13) | -0.11 / -0.14 (9/13) | -0.10 / -0.15 (8/13) | +0.21 / +0.24 (10/13) | +0.32 / +0.40 (9/13) |
| taxa: Archaea (CLR) (1216) | -0.13 / -0.27 (13/13) | -0.11 / -0.03 (9/13) | +0.40 / +0.27 (10/13) | -0.54 / -0.51 (13/13) | -0.24 / -0.30 (13/13) | +0.01 / +0.02 (6/13) | +0.01 / +0.03 (8/13) | -0.45 / -0.57 (13/13) |
| taxa: unclassified (logit) (1216) | -0.23 / -0.21 (11/13) | +0.12 / +0.05 (7/13) | -0.23 / +0.05 (6/13) | +0.29 / +0.18 (8/13) | -0.15 / -0.19 (10/13) | +0.06 / +0.07 (8/13) | +0.10 / +0.11 (10/13) | +0.16 / +0.19 (9/13) |
| sequencing: log total counts (1216) | -0.15 / -0.07 (7/13) | -0.37 / -0.17 (9/13) | -0.59 / -0.45 (12/13) | -0.23 / -0.31 (13/13) | -0.06 / -0.09 (7/13) | -0.06 / -0.08 (6/13) | -0.30 / -0.47 (12/13) | -0.02 / -0.07 (9/13) |
| sequencing: instrument 3310 (1216) | -0.41 / — | +0.02 / — | -0.46 / — | +0.31 / — | +0.04 / — | +0.02 / — | +0.03 / — | +0.05 / — |
| sampling: year (1216) | -0.44 / -0.06 (7/13) | +0.08 / -0.19 (6/13) | -0.47 / -0.00 (5/13) | +0.32 / +0.11 (4/13) | +0.02 / -0.06 (7/13) | -0.00 / -0.17 (6/13) | +0.05 / +0.03 (10/13) | +0.04 / +0.13 (5/13) |
| sampling: depth (m) (1216) | +0.28 / +0.00 (4/6) | -0.03 / -0.04 (4/6) | +0.34 / -0.04 (5/6) | -0.23 / -0.06 (6/6) | -0.01 / +0.01 (4/6) | -0.03 / -0.03 (3/6) | -0.06 / -0.11 (5/6) | -0.07 / -0.03 (5/6) |
| sampling: solar hour (cos, noon +) (1216) | -0.03 / -0.04 (8/13) | -0.09 / -0.08 (11/13) | +0.02 / +0.08 (11/13) | -0.02 / -0.06 (10/13) | -0.05 / -0.06 (10/13) | -0.07 / -0.07 (11/13) | +0.02 / +0.01 (7/13) | -0.09 / -0.07 (11/13) |
| sampling: solar hour (sin, morning +) (1216) | -0.07 / -0.08 (10/13) | -0.03 / -0.00 (8/13) | -0.06 / -0.00 (6/13) | -0.05 / -0.02 (6/13) | -0.07 / -0.07 (10/13) | +0.02 / +0.02 (5/13) | -0.03 / -0.03 (9/13) | +0.01 / +0.00 (7/13) |
| sampling: season (cos, local summer +) (1216) | +0.13 / -0.06 (7/13) | -0.14 / +0.03 (6/13) | -0.23 / -0.08 (6/13) | -0.07 / +0.05 (4/13) | -0.01 / -0.02 (7/13) | +0.03 / +0.04 (6/13) | -0.03 / -0.02 (7/13) | +0.09 / +0.12 (7/13) |
| sampling: |latitude| (1216) | -0.08 / +0.09 (5/13) | -0.04 / -0.18 (10/13) | -0.00 / +0.18 (11/13) | +0.14 / -0.01 (8/13) | +0.05 / +0.08 (5/13) | -0.14 / -0.19 (10/13) | +0.12 / +0.12 (10/13) | -0.09 / -0.13 (9/13) |
| in situ: temperature (526) | -0.22 / -0.33 (4/5) | +0.17 / +0.12 (3/5) | +0.13 / -0.22 (4/5) | -0.18 / -0.19 (3/5) | -0.23 / -0.33 (4/5) | +0.15 / +0.19 (3/5) | -0.19 / -0.21 (4/5) | +0.03 / -0.01 (3/5) |
| in situ: log phosphate (266) | -0.05 / -0.06 (3/4) | -0.26 / +0.01 (2/4) | -0.15 / -0.01 (2/4) | +0.22 / -0.06 (2/4) | -0.10 / -0.05 (3/4) | +0.02 / +0.01 (2/4) | -0.05 / -0.03 (3/4) | -0.02 / -0.03 (2/4) |
| in situ: log nitrate+nitrite (176) | -0.04 / -0.05 (2/4) | -0.33 / -0.13 (2/4) | -0.13 / +0.05 (2/4) | +0.09 / -0.07 (2/4) | -0.10 / -0.06 (2/4) | -0.09 / -0.11 (2/4) | -0.20 / -0.18 (3/4) | -0.03 / -0.04 (2/4) |
| in situ (CCHDO): oxygen (296) | +0.16 / +0.20 (3/4) | -0.26 / -0.14 (2/4) | -0.05 / +0.38 (3/4) | +0.04 / +0.08 (1/4) | +0.14 / +0.19 (3/4) | -0.18 / -0.21 (2/4) | +0.25 / +0.31 (3/4) | -0.19 / -0.23 (2/4) |
| model: mixed layer (1216) | -0.13 / -0.08 (9/13) | +0.19 / +0.02 (5/13) | +0.31 / +0.16 (11/13) | -0.08 / -0.14 (11/13) | -0.08 / -0.10 (10/13) | +0.01 / +0.02 (5/13) | +0.06 / +0.06 (8/13) | -0.17 / -0.19 (11/13) |
| model: mixed layer, 30-day max (1216) | -0.22 / -0.04 (8/13) | +0.16 / -0.11 (8/13) | +0.25 / +0.23 (9/13) | +0.12 / +0.02 (4/13) | -0.03 / -0.04 (7/13) | -0.08 / -0.12 (7/13) | +0.10 / +0.13 (8/13) | -0.06 / -0.06 (8/13) |
| model: chlorophyll (1216) | -0.21 / -0.08 (8/13) | +0.02 / -0.07 (8/13) | +0.09 / +0.01 (7/13) | +0.06 / -0.03 (10/13) | -0.07 / -0.07 (9/13) | -0.02 / -0.03 (7/13) | -0.03 / -0.06 (10/13) | +0.04 / +0.05 (6/13) |
| model: primary production (1216) | -0.21 / -0.07 (7/13) | +0.01 / -0.03 (9/13) | +0.09 / -0.03 (7/13) | +0.00 / -0.07 (7/13) | -0.06 / -0.05 (9/13) | +0.01 / +0.01 (6/13) | -0.09 / -0.11 (9/13) | +0.04 / +0.04 (5/13) |
| model: nitrate (1216) | -0.17 / +0.01 (5/13) | -0.02 / -0.13 (7/13) | +0.03 / +0.06 (7/13) | +0.10 / -0.06 (8/13) | -0.01 / +0.00 (4/13) | -0.08 / -0.12 (7/13) | -0.04 / -0.05 (8/13) | -0.01 / -0.03 (7/13) |
| model: iron (637) | +0.22 / -0.05 (4/7) | +0.10 / +0.08 (6/7) | -0.32 / -0.09 (4/7) | +0.02 / +0.00 (2/7) | -0.06 / -0.06 (4/7) | +0.07 / +0.10 (6/7) | +0.05 / +0.08 (4/7) | +0.01 / +0.01 (3/7) |
| model: wind, 14-day mean (1216) | -0.21 / -0.12 (9/13) | +0.04 / -0.14 (10/13) | +0.19 / +0.26 (11/13) | -0.00 / -0.08 (9/13) | -0.12 / -0.13 (10/13) | -0.11 / -0.13 (10/13) | +0.14 / +0.15 (11/13) | -0.14 / -0.16 (10/13) |
| model: eddy kinetic energy, 30 days (1216) | -0.02 / +0.08 (8/13) | -0.18 / -0.10 (6/13) | +0.11 / -0.05 (8/13) | -0.04 / +0.03 (5/13) | +0.06 / +0.07 (8/13) | -0.09 / -0.10 (6/13) | -0.03 / -0.04 (9/13) | +0.04 / +0.04 (6/13) |

**What the components are** (BRITE C pathways by z; the KOs with the largest loadings):

- **F1** + pathways: ABC transporters (+5.7); Quorum sensing (+5.6); Glycine, serine and threonine metabolism (+4.9) | − pathways: Motor proteins (-9.3); Calcium signaling pathway (-7.1); Autophagy - animal (-7.0)
  + KOs: K07258 dacC, dacA, dacD (+0.89); K02897 RP-L25, rplY (+0.88); K07090 K07090 (+0.88); K00605 gcvT, AMT (+0.88); K00674 dapD (+0.88) | − KOs: K10357 MYO5 (-0.74); K02183 CALM (-0.74); K14572 MDN1, REA1 (-0.73); K11654 SMARCA5, SNF2H, ISWI (-0.73); K21867 AKT, KAT, GORK, SKOR (-0.72)
- **F2** + pathways: Photosynthesis (+12.5); Ribosome (+9.4); Porphyrin metabolism (+7.9) | − pathways: Flagellar assembly (-6.8); Benzoate degradation (-4.5); Butanoate metabolism (-4.4)
  + KOs: K01737 queD, ptpS, PTS (+0.91); K10206 E2.6.1.83 (+0.89); K03789 rimI (+0.89); K01739 metB (+0.88); K05895 cobK-cbiJ (+0.88) | − KOs: K03092 rpoN (-0.50); K02396 flgK (-0.49); K03572 mutL (-0.49); K03694 clpA (-0.49); K02407 fliD (-0.48)
- **F3** + pathways: Ribosome (+11.0); Oxidative phosphorylation (+7.3); Ribosome biogenesis in eukaryotes (+3.7) | − pathways: Flagellar assembly (-8.0); Porphyrin metabolism (-5.6); Homologous recombination (-5.4)
  + KOs: K02907 RP-L30, MRPL30, rpmD (+0.78); K09699 DBT, bkdB (+0.72); K08681 pdxT, pdx2 (+0.70); K00340 nuoK (+0.70); K06215 pdxS, pdx1 (+0.70) | − KOs: K03631 recN (-0.75); K03695 clpB (-0.74); K00748 lpxB (-0.68); K03555 mutS (-0.67); K06158 ABCF3 (-0.67)
- **F4** + pathways: Motor proteins (+5.0); Calcium signaling pathway (+4.6); cGMP-PKG signaling pathway (+4.2) | − pathways: Bacterial chemotaxis (-5.4); Biofilm formation - Vibrio cholerae (-5.0); Bacterial secretion system (-3.5)
  + KOs: K04566 lysK (+0.58); K04782 pchB (+0.55); K01477 alc, ALLC (+0.54); K22073 IBA57 (+0.53); K16846 suyB (+0.52) | − KOs: K07004 K07004 (-0.68); K03615 rnfC (-0.65); K03412 cheB (-0.64); K09781 K09781 (-0.63); K07787 cusA, silA (-0.63)
- **F5** + pathways: Motor proteins (+7.8); cGMP-PKG signaling pathway (+5.6); Calcium signaling pathway (+5.6) | − pathways: RNA polymerase (-3.7); Ribosome biogenesis in eukaryotes (-3.2); Valine, leucine and isoleucine degradation (-3.0)
  + KOs: K11412 SIRT2, SIR2L2 (+0.59); K08296 sixA (+0.56); K00472 P4HA (+0.56); K10901 BLM, RECQL3, SGS1 (+0.53); K00476 ASPH (+0.52) | − KOs: K02319 pol (-0.69); K02322 polC (-0.68); K17734 aprX (-0.67); K02118 ATPVB, ntpB, atpB (-0.67); K10725 cdc6A (-0.67)
- **W2** + pathways: Photosynthesis (+12.8); Porphyrin metabolism (+8.8); Ribosome (+7.5) | − pathways: Flagellar assembly (-5.2); Butanoate metabolism (-5.1); Valine, leucine and isoleucine degradation (-4.2)
  + KOs: K01739 metB (+0.91); K01737 queD, ptpS, PTS (+0.88); K10206 E2.6.1.83 (+0.86); K00507 SCD, desC (+0.85); K03087 rpoS (+0.84) | − KOs: K00856 ADK, adoK (-0.50); K00626 ACAT, atoB (-0.48); K01006 ppdK (-0.48); K00525 E1.17.4.1A, nrdA, nrdE (-0.48); K01907 AACS, acsA (-0.46)

**Conclusions**
1. **Two residual axes are robust and one more is near-robust.** Beyond them, the
   components rotate between halves of the data (congruence ≤ 0.58 from F4 and W4 on).
   They account for about 36% of the standardised residual variance; the remaining
   ~64% is not stable shared structure.
2. **Axis 1 (F1 / W1; 19% / 17%; congruence 0.90 / 0.95): bacterial versus eukaryotic
   share.**
   - The + pole is core bacterial genes: dacA, rplY, gcvT, dapD and leuA, with
     loadings around +0.88. They are enriched in ABC transporters and amino acid
     metabolism.
   - The − pole is eukaryotic genes: myosin, calmodulin, chromatin remodellers,
     autophagy.
   - It follows Bacteria CLR within every cruise (13 of 13). Taxa predict it (CV R2
     0.56 / 0.46); no environmental group does (all ≤ 0.03).
   - Part of it is cruise-level (eta2 0.36). Across cruises it lines up with the
     sequencing instrument (r −0.41) and year (−0.44), which fits a protocol batch:
     filter size or extraction.
   - Within cruise it falls with in-situ temperature (r −0.33, 4 of 5 cruises),
     although temperature is in the model: the temperature spline misses something.
   - **Correction (R34):** that −0.33 used only the 526 tiles of the 5 older cruises
     with manifest temperature. On all 1,216 tiles (the model's own temperature) the
     within-cruise r is −0.03, and the sign splits by sequencing era: negative in 4 of
     the 5 older cruises, positive (+0.03 to +0.30) in all 8 newer ones.
3. **Axis 2 (F2 / W2; 10%; congruence 0.77 / 0.83): phototroph share.**
   - The + pole is picocyanobacterial genes (metB, queD, K10206, cobK, desC, rpoS) and
     the photosynthesis and porphyrin pathways, together with Cryptophyta.
   - The − pole is flagellar and motility genes and heterotrophic degradation.
   - Within cruise, taxa explain it (R2 0.52). No surface or model field explains it:
     chlorophyll, primary production, mixed layer, nitrate and iron all have
     |r| ≤ 0.13 within cruise.
   - It is the most plausible subsurface-driver axis: the share of picocyanobacteria
     given the surface state. Candidates are nutricline or deep chlorophyll maximum
     depth and light history. None of these is among the fields tested.
4. **Axis 3 (F3 / W3; 6–7%; congruence 0.65): sequencing depth, a technical axis.**
   - The + pole is ribosome and oxidative phosphorylation genes; the − pole is
     flagellar, recombination and repair genes.
   - It tracks log total counts: r −0.59 over all tiles and −0.45 within cruise, in 12
     of 13 cruises. Sequencing variables predict it (CV R2 0.32 / 0.28).
   - This is a count-depth effect on the CLR, from the pseudocount on low counts. It
     belongs in the batch column. Log total counts could be added as a covariate, or
     the libraries rarefied, before modelling.
5. **F5 / W4 (5%; congruence 0.54, not stable): archaea versus eukaryotes.**
   - Its − pole is Thermoplasmatota, marine group II archaea (r −0.51 within cruise,
     11 of 13 cruises).
   - It weakly follows a shallower mixed layer within cruise (−0.19, 11 of 13).
   - It is a weak candidate for a mixing or depth-structure driver.
6. **Ocean-model fields capture none of the residual axes.**
   - GLORYS plus history, and PISCES on its 637 tiles, give CV R2 ≤ 0.03 for every
     component.
   - No model field reaches |r| 0.3 within cruise on a stable component.
   - This matches R25 and R28: the model fields add nothing because what is left is
     not in them.
7. **Diel timing and sampling depth explain nothing.** Local solar hour has |r| ≤ 0.09
   on every component. Sampling depth has within-cruise |r| ≤ 0.11 (it is all surface).

**For the manuscript ("what remains unexplained"):**
- About 84% of KO CLR variance is left after Granite plus temperature, a fifth of it
  as cruise offsets.
- About a third of the residual variance sits on three reproducible axes:
  - the bacterial versus eukaryotic share of the library, partly protocol batch
  - the picocyanobacterial or phototroph share, a candidate subsurface or light driver
  - sequencing depth, technical
- None is explained by ocean-model fields.

**Next checks:**
- Per-cruise protocol metadata (filter pore size, volume, extraction kit) against
  axis 1.
- Nutricline or deep-chlorophyll-maximum depth (CCHDO profiles, BGC-Argo) against
  axis 2.
- Refit with log total counts as a covariate to remove axis 3.
- When the external sets are re-mapped (R32), project their residuals on these axes.

## R34. The residual axes, followed up: protocol, sequencing depth, temperature × region, spatial scale (2026-11-19)

**Question.** Follow-ups to R33:
- Do protocol differences explain axis 1 (bacteria vs eukaryotes) or axis 3
  (sequencing depth)? Candidates: filter, volume, extraction, sequencing run.
- Does adding sequencing depth as a covariate remove axis 3?
- Does the temperature effect differ by region?
- Is axis 2 (phototrophs) small-scale and stochastic, or smooth over hundreds of km?
  If smooth, a large-scale driver is missing.

**Scripts and outputs.** All new files for this section are in `transfer/`:
- `transfer/scripts/protocol_metadata.py` builds the protocol table, written to
  `transfer/results/protocol/`:
  - `sample_protocol.csv`: one row per modelled sample
  - `cruise_coverage.csv`: per cruise, the share of samples with each field and its
    values
- `transfer/scripts/residual_axes_followup.py --parts protocol depth region spatial`
  runs the tests, written to `transfer/results/residual_followup/`. Per part:
  - protocol: `tile_protocol.csv`, `cruise_protocol.csv`, `cruise_level_tests.csv`,
    `within_cruise_labels.csv`, `within_cruise_numeric.csv`
  - depth: `depth_refit_components.csv`, `depth_refit_ko_r2.csv`
  - region: `region_interaction.csv`
  - spatial: `correlogram.csv`, `correlogram_fit.csv`, `replicates.csv`
  - all parts: the `summary_<parts>.txt` console logs
- Figures: `transfer/figures/w2_along_track.png` and
  `transfer/figures/residual_axes_correlogram.png`.

To reproduce, run from the repository root:

    ~/miniforge3/bin/python transfer/scripts/protocol_metadata.py
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/residual_axes_followup.py --parts protocol depth region spatial

- Same software as R33. Seed 0.
- The script reads the R33 outputs (`data/results/residual_modules/scores.csv`,
  `loadings.csv.gz`), so run `scripts/residual_modules.py` first.
- It recomputes the R33 baseline and caches it in `residuals.npz`.
- Run time is about 2 minutes on a laptop.

### 1. What protocol metadata exist

**Source.** The FAIRe workbook `FAIRe-NOAA_uci-ess-biogoship.xlsx`:
- `sampleMetadata`
- `experimentRunMetadata`, rows with assay `MetaG_LibPrep_Illumina`

**Linking libraries to samples.** Each library is linked to its FAIRe `samp_name` by
`build_ko_targets.library_samples`, the link the KO targets use. 2,058 libraries
link to the 2,007 modelled samples.

**Coverage**, per sample (`cruise_coverage.csv`):

| Field | Recorded for | Values |
|---|---|---|
| Sequencing run (`seq_run_id`: date, flowcell, PE150) | all 16 cruises | 10 flowcells. Cruises share flowcells: 23LT2GLT4 (sequenced 2026) holds A13.5-2024, GOMECC4, I05-2023 and 53 A16N-2023 samples; 22YJMLLT4 (2025) holds I08S, I09N-2025, P02; HMMFVDSXX (2019) holds AMT28, P18, AE1319, BVAL46, NH1418. |
| Instrument | all | OBI_2131 / OBI_2631 for the 2011–2020 cruises, OBI_3310 for the 2021–2025 cruises. Fully collinear with era. |
| Input read count | all | median 1.0e7 – 4.1e7 per cruise |
| Collection method | 14 of 16 cruises | underway / inline vs CTD / Bio-cast. Five cruises have ≥ 5 tiles of each. |
| Volume filtered | A13.5-2024, A16N-2023, GOMECC4, HLY2202, I05-2023 | about 8 L Sterivex underway; GOMECC4 about 2 L from Niskin bottles. The workbook says "L" for values of 830–10,000, which must be mL. |
| Storage solution | A13.5-2024, A16N-2023, I05-2023 (lysis buffer), GOMECC4 (none) | — |
| Extraction plate and well | A13.5-2024, A16N-2023, GOMECC4, I05-2023 | 2–5 plates per cruise |
| Size fraction, filter, extraction kit, DNA yield | GOMECC4 only | 0.22 µm Sterivex-GP PES, no pre-filter; ZymoBIOMICS 96 MagBead kit (KingFisher protocol); yield median 295 ng |

**Not recorded for the 9 older cruises:** filter pore size, pre-filter, volume,
storage, extraction kit, extraction date. For them, the cruise is the only protocol
label. Answering whether the protocol differs between the 2011–2020 and 2021–2025
cruises needs Bio-GO-SHIP: pre-filter and extraction kit in particular. The older
cruises' protocol is in their methods papers; it isn't in the workbook.

**Derived per tile** (mean over the tile's samples; mode for labels):
- storage time = sequencing date − sampling date. The run date "20189999" is taken as
  2018-07-01.
- log10 input reads
- log10 KO yield = total KO counts per input read pair. This is the share of reads
  that map to annotated references at ≥ 50% coverage.

### 2. Protocol against the axes

**a) Cruise level.**
- Per-cruise mean score, in SDs of the tile scores (`cruise_protocol.csv`):

| Cruise | Tiles | Flowcell | Sequenced | Sampled | Storage (y) | Lysis buffer | Instrument | log KO yield | log input reads | F1 | F2 | F3 |
|---|---:|---|---:|---:|---:|---|---|---:|---:|---:|---:|---:|
| A13.5-2024 | 93 | 23LT2GLT4 | 2026.26 | 2024.15 | 2.11 | yes | OBI_3310 | −0.61 | 7.33 | −0.22 | −0.10 | −0.69 |
| A16N-2023 | 175 | 22NCJYLT4 | 2025.26 | 2023.26 | 2.02 | yes | OBI_3310 | −0.65 | 7.06 | +0.23 | +0.55 | −0.69 |
| AMT28 | 46 | HMMFVDSXX | 2019.57 | 2018.78 | 0.80 | not recorded | OBI_2631 | −0.41 | 7.24 | +0.25 | −0.05 | −0.13 |
| C13.5/A13.5 | 216 | HKNJLDRXX | 2020.69 | 2020.27 | 0.42 | not recorded | OBI_2631 | −0.62 | 7.01 | +0.47 | +0.16 | +0.53 |
| GOMECC4 | 53 | 23LT2GLT4 | 2026.26 | 2021.76 | 4.50 | none | OBI_3310 | −0.55 | 7.21 | +0.10 | −0.85 | −0.24 |
| HLY2202-2022 | 32 | HFJLFDSX7 | 2023.59 | 2022.70 | 0.89 | not recorded | OBI_2631 | −0.69 | 7.37 | +0.78 | +0.18 | −1.01 |
| I05-2023 | 141 | 23LT2GLT4 | 2026.26 | 2023.64 | 2.62 | yes | OBI_3310 | −0.90 | 7.10 | −1.29 | +0.70 | +0.60 |
| I07N | 146 | H357JDSXY | 2020.21 | 2018.36 | 1.86 | not recorded | OBI_2631 | −0.59 | 7.03 | +0.27 | −0.10 | +0.71 |
| I08S-2024 | 23 | 22YJMLLT4 | 2025.58 | 2024.18 | 1.40 | not recorded | OBI_3310 | −0.58 | 7.61 | −0.16 | −0.49 | −0.40 |
| I09N | 54 | HHGJFDSXX | 2019.21 | 2016.31 | 2.90 | not recorded | OBI_2631 | −0.74 | 6.88 | +0.62 | +0.15 | +1.46 |
| I09N-2025 | 84 | 22YJMLLT4 | 2025.58 | 2025.27 | 0.31 | not recorded | OBI_3310 | −0.50 | 7.61 | −0.78 | −0.38 | −1.08 |
| P02-2022 | 89 | 22YJMLLT4 | 2025.58 | 2022.42 | 3.16 | not recorded | OBI_3310 | −0.40 | 7.62 | −0.29 | −0.91 | −0.64 |
| P18 | 64 | HMMFVDSXX | 2019.57 | 2016.96 | 2.61 | not recorded | OBI_2631 | −0.45 | 7.18 | +0.63 | −0.74 | +0.17 |

- Spearman r across the 13 cruise means:

| Cruise median | F1 | F2 | F3 |
|---|---:|---:|---:|
| sampling year | −0.71 | +0.05 | −0.72 |
| sequencing year | −0.78 | −0.15 | −0.45 |
| storage time | −0.14 | −0.30 | +0.31 |
| log KO yield | −0.16 | −0.82 | −0.24 |
| log input reads | −0.45 | −0.52 | −0.73 |
| flowcell grouping, eta2 over cruises (permutation p, 9,999) | 0.68 (0.36) | 0.49 (0.70) | 0.84 (0.08) |

- Mean by instrument:
  - OBI_2631 (6 cruises): F1 +0.50, F2 −0.07, F3 +0.29
  - OBI_3310 (7 cruises): F1 −0.34, F2 −0.21, F3 −0.45

**b) Within cruise, labels.**
- Statistic: the pooled within-cruise eta2 of the W score, by label.
- Null: the labels are rotated along the sampling order inside each cruise, 999 times.
  Contiguous blocks stay together, so a label that only marks a stretch of the track
  isn't mistaken for batch.

| Label (cruises, tiles) | W1 eta2 (p) | W2 eta2 (p) | W3 eta2 (p) |
|---|---|---|---|
| collection, flow-through vs cast (A13.5-2024, A16N-2023, I05-2023, I07N, P02-2022; 644) | 0.016 (0.12) | 0.011 (0.29) | 0.020 (0.012) |
| sequencing flowcell (A16N-2023 only; 175) | 0.007 (0.58) | 0.041 (0.14) | 0.060 (0.13) |
| extraction plate (A13.5-2024, A16N-2023, GOMECC4, I05-2023; 462) | 0.030 (0.91) | 0.095 (0.21) | 0.103 (0.13) |

- Per-cruise plate eta2 for W3: A13.5 0.06, A16N 0.07, GOMECC4 0.47, I05 0.09.
- GOMECC4's A/B/C replicates are spread across plates, and W3 is mostly within-tile
  noise (section 5). Its 0.47 is from 53 tiles and isn't significant pooled.

**c) Within cruise, numbers.** Spearman r, all tiles / within cruise; in brackets, the
number of cruises with the pooled sign.

| Variable (tiles, cruises) | F1 | F2 | F3 | W1 | W2 | W3 |
|---|---|---|---|---|---|---|
| volume filtered (494, 5) | −0.48 / +0.03 (2/5) | +0.33 / −0.08 (2/5) | +0.40 / −0.12 (4/5) | −0.02 / −0.02 (3/5) | −0.02 / −0.06 (2/5) | +0.00 / −0.02 (2/5) |
| log input reads (1,216, 13) | −0.35 / −0.22 (12/13) | −0.33 / −0.20 (11/13) | −0.55 / −0.36 (11/13) | −0.16 / −0.24 (12/13) | −0.06 / −0.10 (9/13) | −0.23 / −0.40 (11/13) |
| log KO yield (1,216, 13) | +0.18 / +0.15 (10/13) | −0.27 / −0.04 (8/13) | −0.38 / −0.33 (12/13) | +0.11 / +0.14 (11/13) | −0.01 / −0.01 (6/13) | −0.27 / −0.34 (11/13) |
| storage time (1,216, 13) | −0.11 / −0.00 (6/13) | −0.07 / +0.03 (5/13) | +0.00 / +0.10 (6/13) | +0.01 / −0.01 (6/13) | +0.02 / +0.01 (5/13) | −0.00 / +0.06 (4/13) |

**Reading.**
- **Axis 1, bacteria vs eukaryotes.**
  - Its cruise offsets follow era: sampling year r −0.71, sequencing year −0.78, and
    OBI_2631 vs OBI_3310 +0.50 vs −0.34. The newer cruises are more eukaryote-rich in
    the residuals.
  - Year, instrument and sequencing run can't be separated: they change together.
  - Not explained by:
    - flowcell (cruises sharing one are no more alike: p 0.36)
    - storage time (−0.14)
    - volume (within cruise +0.03)
    - collection method (eta2 0.016)
    - extraction plate (0.030, p 0.91)
  - So it is either a protocol change between eras that the workbook doesn't record (a
    pre-filter, the extraction kit, or library prep), or a real change. Its within-cruise
    part (W1) has no protocol correlate. Within cruise, a deeper library tilts slightly
    to eukaryotes (input reads r −0.24, 12 of 13 cruises).
- **Axis 3, sequencing depth.** Confirmed technical. It falls with input reads
  (within r −0.40) and with KO yield (−0.34). It is mostly sample-to-sample noise
  (section 5).
- **Axis 2, phototrophs.** No protocol label explains it within cruise:
  - collection 0.011
  - plate 0.095, p 0.21
  - volume −0.06

  Its cruise offsets follow KO yield (r −0.82 across cruises): cruises whose reads map
  less to the references sit higher on axis 2. That is a reference-coverage effect at
  the cruise level, and it should be checked when the external sets are mapped.

### 3. Axis 3: sequencing depth as a covariate

**Change.** The R33 baseline plus log10 total KO counts, as its own spline (8 knots,
cubic), with one cruise held out (`added_fields_all_kos.cv` with the variable).

**Prediction.**
- Median KO R2: 0.103 → 0.145. The median per-KO gain is +0.007.
- 40% of KOs gain more than 0.02; 17% lose more than 0.02.
- Residual share of KO CLR variance: 83.7% → 82.2%.

**Residual PCA after the refit.** The F3 / W3 depth axis is gone: no component has
|r| > 0.16 with log total counts. The two robust axes stay unchanged.

| After refit | Share | Split-half | Best match in R33 (cosine) | Cosine with old F3 / W3 | r with log total, all / within |
|---|---:|---:|---|---:|---|
| F1' | 20.7% | 0.93 | F1 (+1.00) | −0.03 | +0.04 / +0.07 |
| F2' | 10.2% | 0.82 | F2 (+0.94) | −0.32 | −0.02 / +0.06 |
| F3' | 5.4% | 0.52 | F4 (+0.87) | −0.31 | +0.14 / −0.14 |
| F4' | 5.4% | 0.54 | F5 (+0.76) | −0.42 | +0.03 / −0.09 |
| F5' | 3.6% | 0.38 | F7 (−0.62) | +0.51 | −0.12 / +0.07 |
| F6' | 3.1% | 0.30 | F6 (+0.98) | +0.04 | +0.16 / +0.12 |
| W1' | 18.2% | 0.95 | W1 (+1.00) | −0.02 | +0.05 / +0.05 |
| W2' | 10.1% | 0.86 | W2 (+1.00) | −0.07 | +0.05 / +0.06 |
| W3' | 5.0% | 0.65 | W4 (+0.86) | −0.47 | −0.02 / −0.06 |
| W4' | 4.5% | 0.71 | W3 (+0.67) | +0.67 | −0.00 / −0.01 |
| W5' | 2.6% | 0.42 | W6 (+0.98) | +0.05 | +0.10 / +0.15 |
| W6' | 2.2% | 0.32 | W8 (+0.75) | −0.24 | +0.04 / +0.04 |

**Remainder.** W4' keeps part of the old W3 pattern (cosine 0.67) without the depth
correlation. This is the ribosome / oxidative-phosphorylation vs repair / motility
contrast, independent of depth.

**Recommendation.**
- Use log total counts as a nuisance covariate in the models and the maps. For a map,
  predict at a fixed reference depth, such as the median.
- The cause is the CLR pseudocount 0.5 / total. Rarefying, or a depth-matched
  pseudocount, would remove it at the source.

### 4. Temperature × region

**Change.**
- The temperature spline columns are multiplied by each region's indicator and added
  to the R33 baseline, cruise held out.
- A region with no training cruise falls back to the shared curve.
- The interaction columns are scaled by shrink = 1, 0.3 or 0.1, which is a 1×, 11× or
  100× stronger ridge penalty on them.

**Regions:**
- basin:
  - Atlantic: A13.5-2024, A16N-2023, AMT28, C13.5/A13.5, GOMECC4
  - Indian: I05-2023, I07N, I08S-2024, I09N, I09N-2025
  - Pacific: P02-2022, P18
  - Arctic: HLY2202-2022
- latitude band: |lat| < 20, 20–35, > 35
- sequencing era: OBI_2131/2631 vs OBI_3310. Added because of the sign split below.

| Interaction | Shrink | Median KO R2 (baseline 0.103) | KOs gaining > 0.02 | KOs losing > 0.02 |
|---|---:|---:|---:|---:|
| ocean basin | 1 | 0.041 | 16% | 69% |
| ocean basin | 0.3 | 0.049 | 16% | 65% |
| ocean basin | 0.1 | 0.079 | 19% | 50% |
| latitude band | 1 | 0.065 | 13% | 58% |
| latitude band | 0.3 | 0.082 | 14% | 46% |
| latitude band | 0.1 | 0.093 | 13% | 28% |
| sequencing era | 1 | 0.113 | 36% | 43% |
| sequencing era | 0.3 | 0.113 | 36% | 42% |
| sequencing era | 0.1 | 0.126 | 41% | 30% |

**Within regions.** Basin, shrink 1: Atlantic −0.045, Indian −0.081, Pacific −0.078,
Arctic +0.018. Shrink 0.1: Pacific +0.008 is the only region that gains.

**Reading.**
- Separate temperature curves per region don't carry over to a held-out cruise. They
  fit that region's other cruises, and every basin has only 1–5 cruises. One shared
  curve transfers best.
- A curve per sequencing era gives a small gain at the strongest shrinkage (median KO
  R2 0.126). Era is a batch, not a region.
- **R33 correction.** The within-cruise W1–temperature r of −0.33 came from 526 tiles
  of 5 older cruises:
  - On all 1,216 tiles it is −0.03 (4 of 13 cruises with that sign).
  - Per cruise, the older cruises are negative (AMT28 −0.34, C13.5 −0.45, I07N −0.27,
    P18 −0.13), except I09N (+0.44).
  - All 8 newer cruises are positive (+0.03 to +0.30).

  Whatever W1 tracks along temperature differs by era, which again points to a
  protocol difference between eras rather than a regional temperature response.

### 5. Axis 2 along the transects: spatial scale

**Method.**
- Positions are sample positions, averaged per tile. Tiles' nearest neighbour in the
  same cruise: median 44 km (quartiles 22–63 km).
- Correlogram:
  - pairs of tiles in the same cruise
  - the score is centred per cruise
  - r = Σ xᵢxⱼ / √(Σxᵢ² Σxⱼ²) per distance bin, pooled over cruises, with ≥ 30 pairs
    per bin
- Comparators:
  - in-situ temperature, PISCES chlorophyll and GLORYS mixed layer (smooth fields)
  - the median correlogram of 500 random single-KO within-cruise residuals
- Replicate check: the 145 tiles holding two or more samples. GOMECC4 contributes 52,
  as A/B/C replicates; the rest are co-located samples, median distance 0 km.
  - Each sample's residual = its CLR − the tile prediction − the cruise mean residual.
  - It is standardised with the tile matrix's KO SDs and projected on the
    component's loadings.
  - Statistic: the correlation between the two samples of a tile.

**Correlation of pairs by distance (km):**

| Series | 0–15 | 15–40 | 40–75 | 75–150 | 150–300 | 300–600 | 600–1200 | 1200–2400 | 2400–5000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| W1 (bacteria vs eukaryotes) | 0.57 | 0.42 | 0.44 | 0.30 | 0.19 | 0.08 | −0.02 | −0.02 | −0.04 |
| **W2 (phototrophs)** | **0.55** | **0.50** | **0.49** | **0.39** | **0.32** | **0.25** | **0.15** | −0.02 | −0.09 |
| W3 (sequencing depth) | 0.35 | 0.29 | 0.30 | 0.19 | 0.10 | 0.06 | −0.00 | −0.05 | −0.03 |
| single-KO residuals (median of 500) | 0.38 | 0.34 | 0.32 | 0.25 | 0.18 | 0.11 | 0.03 | −0.03 | −0.06 |
| in-situ temperature | 1.00 | 0.99 | 0.99 | 0.97 | 0.95 | 0.90 | 0.73 | 0.34 | −0.34 |
| log PISCES chlorophyll | 1.00 | 0.97 | 0.96 | 0.92 | 0.83 | 0.63 | 0.25 | −0.16 | −0.23 |
| log GLORYS mixed layer | 0.93 | 0.89 | 0.84 | 0.74 | 0.65 | 0.47 | 0.16 | −0.10 | −0.21 |
| pairs | 102 | 394 | 811 | 1,542 | 2,898 | 5,975 | 11,502 | 16,583 | 23,937 |

**Exponential fit** r(d) = c·exp(−d/L), weighted by pairs:
- W1: c 0.53, L 201 km
- W2: c 0.51, L 557 km
- W3: c 0.34, L 206 km
- temperature: c 1.30, L 1,029 km
- chlorophyll: c 1.25, L 480 km
- mixed layer: c 1.03, L 442 km

The smooth fields fit poorly (c > 1): their correlograms are flat and then fall,
rather than exponential. For them, read the table, not the fit.

**Two samples in one tile:**

| Component | r between the two samples | Sample-to-sample variance, as a share of the tile-score variance |
|---|---:|---:|
| W1 | +0.83 | 18% |
| W2 | +0.72 | 26% |
| W3 | +0.55 | 45% |

**Reading.** The variance of W2 splits into three parts:
- **About 25–30% is sample-to-sample.** Two samples of the same water at the same
  place agree at r 0.72. This part is measurement noise or patchiness below a tile;
  no predictor can reach it.
- **About 20% changes between neighbouring samples**, from r 0.72 at 0 km to about
  0.5 at 15–75 km. That is small-scale or day-to-day variation.
- **About half is smooth over hundreds of km.** r stays at 0.39 at 75–150 km, 0.32 at
  150–300 km, 0.25 at 300–600 km and 0.15 at 600–1,200 km; the fitted decorrelation
  length is about 560 km. That is the longest of the three axes and the slowest-decaying
  residual series. It is about the scale of model chlorophyll (r 0.63 at 300–600 km)
  and mixed layer (0.47).
- The along-track plots (`w2_along_track.png`) show the same picture: long runs of
  one sign over thousands of km on top of point-to-point scatter. Examples:
  A13.5-2024 is positive from 2,000 to 7,000 km, then negative; A16N-2023, I07N and
  I09N-2025 have similar runs.

W1 and W3 behave differently:
- W1 decays over about 200 km and has high replicate agreement (0.83). It is
  measurable, but patchy at the mesoscale.
- W3 is the least structured: replicate r 0.55 and correlation 0.35 even for the
  closest pairs. That fits a technical, per-library effect.

**Conclusion for axis 2.** It isn't mainly stochastic: a large-scale driver is
missing, one that none of the tested surface or model fields captures.
- Within cruise, chlorophyll, primary production, mixed layer, nitrate and iron all
  have |r| ≤ 0.13 with W2 (R33).
- The next tests are the nutricline / deep chlorophyll maximum depth and the recent
  chlorophyll trend.
- The predictability limit: about a quarter of W2 variance (the replicate noise) is
  beyond any predictor.

**Chlorophyll trend over the previous 2–4 weeks.** We don't have it yet.
- The R27 history fields cover SST change over 14 days, the 30-day maximum mixed
  layer, wind over 14 days, the sea-level gradient and 30-day EKE. None is a
  chlorophyll trend.
- PISCES chlorophyll is a single day field, on 637 tiles.
- Options:
  - The PISCES daily analysis, 28 days before each sample. The same Copernicus Marine
    service as R25, on Morrill with the stored credentials. It covers the
    forecast-period tiles only.
  - A gap-free satellite L4 chlorophyll product (GlobColour / OC-CCI daily), for all
    tiles.
- Summarise each as the log-chlorophyll slope over 14 and 28 days, and the 28-day
  maximum relative to the day itself.

### Summary
1. **Protocol.**
   - Axis 1's cruise offsets follow era: year, instrument and sequencing run, which
     change together. Flowcell, storage time, volume, collection method and
     extraction plate don't explain it.
   - The pre-filter and kit of the older cruises aren't in the workbook. They have to
     come from Bio-GO-SHIP.
   - No protocol variable explains the within-cruise axes.
2. **Sequencing depth.** As a covariate it removes axis 3 completely and raises the
   median KO R2 from 0.103 to 0.145. Use it as a nuisance covariate.
3. **Temperature × region.** Doesn't help under cruise hold-out at any shrinkage.
   R33's W1–temperature link was an artefact of the 5-cruise subset; its sign differs
   by era.
4. **Axis 2 is spatially smooth**: decorrelation length about 560 km, half its
   variance on scales of hundreds of km. A large-scale driver is missing. About a
   quarter is replicate-level noise, the predictability floor.
5. **Next:**
   - the nutricline / deep chlorophyll maximum depth (CCHDO CTD profiles, BGC-Argo),
     and the 14- and 28-day chlorophyll trend, against W2
   - from Bio-GO-SHIP: the older cruises' pre-filter, volume and kit, against axis 1

## R35. Cruise methods table, paired protocol samples, nutricline / DCM, chlorophyll history job (2026-11-19)

**Question.**
- Do the documented methods explain the residual axes of R33 and R34? Methods here
  means pre-filter, volume, preservation and extraction, plus the paired protocol
  tests on A16N and P02.
- Does subsurface structure explain axis 2 (phototrophs)? Structure here means
  nitracline depth, DCM depth and strength, and mixed layer.

**Scripts.** All in `transfer/`:
- `scripts/cruise_methods_table.py` writes `results/protocol/cruise_methods.csv`.
- `scripts/cchdo_profiles.py` writes `results/nutricline/station_profiles.csv`.
  - It downloads the public CCHDO bottle and CTD files into `transfer/data/cchdo/`,
    75 MB.
- `scripts/residual_axes_followup.py --parts methods pairs nutricline` runs the
  tests.
  - Outputs go to `results/residual_followup/`:
    - `cruise_methods_scores.csv`, `methods_within_cruise.csv`
    - `pairs.csv`, `pairs_effects.csv`, `pairs_ko_effect.csv`
    - `tile_profiles_{100,30}km.csv`, `nutricline_drivers.csv`
    - `summary_*.txt`
  - It needs the R34 `protocol` part's `cruise_protocol.csv` and `tile_protocol.csv`.

To reproduce, run in this order, from the repository root:

    ~/miniforge3/bin/python transfer/scripts/cruise_methods_table.py
    ~/miniforge3/bin/python transfer/scripts/cchdo_profiles.py          # about 5 min, downloads
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/residual_axes_followup.py --parts protocol methods pairs nutricline

### 1. Master cruise-methods table

**Source.** Compiled from the cruise papers and reports and supplied on 2026-11-19.
- "?" = not confirmed.
- "ND" = not described in the source, which is not proof that the step was absent.
- Volumes refer to DNA samples.
- The GOMECC4 row comes from the FAIRe workbook.

| Cruise | Year | DNA pre-filter | DNA filter | Volume | Preservation / storage | Beads | DNA extraction / purification | Important differences |
|---|---|---|---|---|---|---|---|---|
| BVAL46 | 2011 | None reported | 0.22 µm Sterivex | 2 L | 1.62 mL lysis buffer; −20 °C | No, enzymatic | Lysozyme, proteinase K/SDS; Zymo Genomic DNA Clean & Concentrator | Original protocol |
| AE1319 | 2013 | None reported | 0.22 µm Sterivex | 2 L | 1.62 mL lysis buffer; −20 °C | No, enzymatic | Same original extraction | Original protocol |
| NH1418 | 2014 | None reported | 0.22 µm Sterivex | 2 L | 1.62 mL lysis buffer; −20 °C | No, enzymatic | Same original extraction | Original protocol |
| I09N | 2016 | None reported | 0.22 µm Sterivex | 10 L | 1.62 mL lysis buffer; −20 °C | No, enzymatic | Same original extraction | Underway sampling |
| P18 | 2016–17 | None reported | 0.22 µm Sterivex | 2 L | 1.62 mL lysis buffer; −20 °C | No, enzymatic | Same original extraction | CTD sampling |
| I07N | 2018 | None reported | 0.22 µm Sterivex | 2–4 L CTD; 4 L underway | 1.62 mL lysis buffer; −20 °C | No, enzymatic | Same original extraction | Mixed underway/CTD |
| AMT28 | 2018 | None reported | 0.22 µm Sterivex | 2 L | 1.62 mL lysis buffer; −20 °C | No, enzymatic | Same original extraction | CTD sampling |
| C13.5/A13.5 | 2020 | None reported for DNA | 0.22 µm Sterivex | 5–10 L | 1.62 mL lysis buffer; −20 °C | No, enzymatic | Same original extraction | Underway; POM sampling separately used 30 µm mesh |
| P02W | 2022 | 30 µm underway; CTD ? | 0.22 µm Sterivex | Up to 8 L | 1.62 mL lysis buffer; −80 °C | ND | ? | Source-dependent pre-filtration |
| P02E | 2022 | ? | 0.22 µm Sterivex | 4–8 L | Lysis buffer normally; some DNA/RNA Shield; frozen | ND | ? | 10 paired buffer-comparison samples |
| HLY2202 | 2022 | ? | 0.22 µm Sterivex | 2–8 L | ? | No, enzymatic extraction | Lysozyme, proteinase K/SDS; Zymo purification | Old extraction method on newer cruise |
| A16N, leg 1 | 2023 | ND | 0.22 µm Sterivex | ~8 L | Predominantly DNA/RNA Shield; −80 °C | Yes, usual protocol | ? | Paired old/new protocols; underway pump changed |
| A16N, leg 2 | 2023 | ND | 0.22 µm Sterivex | ? | DNA/RNA Shield or lysis buffer by protocol | Protocol-dependent | ? | Paired old/new protocols |
| I05 | 2023 | ND | 0.22 µm Sterivex | ~8 L | 1 mL DNA/RNA Shield; −80 °C | Yes | ? | Newer bead/Shield workflow |
| A13.5 | 2024 | ND | 0.22 µm Sterivex | ~8 L | 1 mL DNA/RNA Shield; −80 °C | Yes | ? | Newer bead/Shield workflow |
| I08S | 2024 | ND | Sterivex; pore size ? | ~8 L | 1 mL DNA/RNA Shield; −80 °C | Yes | ? | Newer bead/Shield workflow |
| I09N | 2025 | ND | 0.22 µm Sterivex | ~9 L | 1 mL DNA/RNA Shield; −80 °C | ? | ? | Reoccupation of 2016 line |
| GOMECC4 (workbook) | 2021 | no pre-filter | 0.22 µm Sterivex-GP PES | ~2 L, Niskin | none; −80 °C | Yes | ZymoBIOMICS 96 MagBead kit, KingFisher | AOML Niskin samples, A/B/C replicates |

**Coded columns** in the CSV, used for the tests:
- preservation: lysis buffer, Shield, mixed, none or unknown
- extraction: enzymatic, beads, mixed or unknown
- pre-filter
- volume_l, the midpoint of the stated range
- store_c
- the leg split dates: P02 at 2022-06-13, A16N at 2023-04-10

For the cruise-level test, A16N is coded Shield and P02 lysis buffer, their usual
protocols. The paired samples are tested separately in section 3.

### 2. Cruise level: the coded methods against the cruise mean scores

**Per-cruise means** (SD of the tile scores):

| Cruise | F1 | F2 | F3 | Preservation | Extraction | Pre-filter | Volume (L) | Store (°C) | Instrument | Sampled |
|---|---:|---:|---:|---|---|---|---:|---:|---|---:|
| I09N | +0.62 | +0.15 | +1.46 | lysis buffer | enzymatic | none reported | 10 | −20 | OBI_2631 | 2016.3 |
| P18 | +0.63 | −0.74 | +0.17 | lysis buffer | enzymatic | none reported | 2 | −20 | OBI_2631 | 2017.0 |
| I07N | +0.27 | −0.10 | +0.71 | lysis buffer | enzymatic | none reported | 3.5 | −20 | OBI_2631 | 2018.4 |
| AMT28 | +0.25 | −0.05 | −0.13 | lysis buffer | enzymatic | none reported | 2 | −20 | OBI_2631 | 2018.8 |
| C13.5/A13.5 | +0.47 | +0.16 | +0.53 | lysis buffer | enzymatic | none reported | 7.5 | −20 | OBI_2631 | 2020.3 |
| HLY2202-2022 | +0.78 | +0.18 | −1.01 | unknown | enzymatic | unknown | 5 | ? | OBI_2631 | 2022.7 |
| GOMECC4 | +0.10 | −0.85 | −0.24 | none | beads | none reported | 2 | −80 | OBI_3310 | 2021.8 |
| A16N-2023 | +0.23 | +0.55 | −0.69 | Shield | beads | unknown | 8 | −80 | OBI_3310 | 2023.3 |
| I05-2023 | −1.29 | +0.70 | +0.60 | Shield | beads | unknown | 8 | −80 | OBI_3310 | 2023.6 |
| A13.5-2024 | −0.22 | −0.10 | −0.69 | Shield | beads | unknown | 8 | −80 | OBI_3310 | 2024.2 |
| I08S-2024 | −0.16 | −0.49 | −0.40 | Shield | beads | unknown | 8 | −80 | OBI_3310 | 2024.2 |
| P02-2022 | −0.29 | −0.91 | −0.64 | lysis buffer | unknown | 30 µm (underway) | 8 | −80 | OBI_3310 | 2022.4 |
| I09N-2025 | −0.78 | −0.38 | −1.08 | Shield | unknown | unknown | 9 | −80 | OBI_3310 | 2025.3 |

**Group means** (number of cruises in brackets):

| Grouping | F1 | F2 | F3 |
|---|---|---|---|
| extraction: enzymatic (6) | +0.50 | −0.07 | +0.29 |
| extraction: beads (5) | −0.27 | −0.04 | −0.28 |
| extraction: unknown (2: P02, I09N-2025) | −0.53 | −0.65 | −0.86 |
| preservation: lysis buffer (6, incl. P02) | +0.33 | −0.25 | +0.35 |
| preservation: Shield (5) | −0.44 | +0.05 | −0.45 |
| instrument OBI_2631 (6) | +0.50 | −0.07 | +0.29 |
| instrument OBI_3310 (7) | −0.34 | −0.21 | −0.45 |

**Spearman r with sampling year, inside each extraction group:**
- enzymatic (6 cruises): F1 +0.09, F2 +0.66, F3 −0.71
- beads (5 cruises): F1 −0.50, F2 +0.10, F3 −0.10

**Spearman r with typical volume, all 13 cruises:** F1 −0.44, F2 +0.19, F3 −0.16.

**Within-cruise contrasts.** W scores in SD; rotation null as in R34.

| Contrast (tiles) | W1 | W2 | W3 |
|---|---|---|---|
| P02E (leg 2) − P02W (leg 1) (36 vs 53) | −0.21 (p 0.52) | −0.22 (p 0.20) | +0.06 (p 0.67) |
| A16N leg 2 − leg 1 (74 vs 101) | −0.14 (p 0.51) | −0.20 (p 0.70) | +0.19 (p 0.63) |
| P02W flow-through, 30 µm pre-filtered, − cast (38 vs 15) | −0.42 (p 0.17) | +0.18 (p 0.47) | +0.05 (p 0.74) |

**Reading.**
- **The axis-1 offset follows the extraction workflow, not the sampling year.**
  - Among the six enzymatic-extraction cruises, F1 doesn't trend with year (r +0.09).
  - HLY2202 is the key cruise: sampled in 2022 but extracted the old way, it sits with
    the old cruises (F1 +0.78).
  - Extraction (enzymatic vs beads/Shield) and instrument (OBI_2631 vs OBI_3310) are
    still fully collinear among the 11 cruises where extraction is known. HLY2202
    doesn't separate them, since it was also sequenced on OBI_2631.
  - So: extraction/preservation and sequencing platform remain confounded, but year is
    not the driver.
- **No leg or pre-filter contrast reaches significance.**
  - P02W's 30 µm pre-filtered underway samples score +0.18 SD on W2 against the casts
    (p 0.47). With 15 casts there is little power. There's no sign that the pre-filter
    removed a particle-associated phototroph signal.
  - The A16N underway pump change shows no leg offset.

### 3. Paired protocol samples, treated as technical treatments

**Pairs.** From the workbook notes, samples ≤ 10 m only:
- **A16N.** 18 underway samples are noted "Protocol A". Each has one partner MG sample
  at the same eventDate, presumably the standard Shield/bead protocol, giving 17
  pairs.
  - The notes don't say which protocol "A" is. The methods table calls these "paired
    old/new protocols", so Protocol A is probably the old lysis-buffer protocol. **To
    be confirmed.**
- **P02.** 11 samples are noted "shield", each with 1–2 partners at the same eventDate
  (lysis buffer), 21 pairs.
  - **All P02 Shield samples and their partners were sequenced on run 20230419_HWKVFDSX5,
    and none of them are in `ko_counts_cov50_all.tsv`.** The P02 libraries in the KO
    table are all from the 2025 run 22YJMLLT4. The P02 test needs those libraries
    counted: they could be run through our mapping jobs (R32), from the raw reads.

**What the R33/R34 tiles did with them.** 16 A16N pairs are among the modelled samples.
In 12 of them both samples fall in the same tile, and the tile model averaged their
CLRs. That mixes the two protocols, so the averaging should change: keep the standard
protocol sample, or carry protocol as a covariate.

**Method.**
- Per pair: treated minus partner, on the sample-level domain CLRs, the
  unclassified-read logit and log10 total KO counts.
- Also on the R33 axis scores at sample level: the sample's CLR minus its tile's
  cross-validated prediction (minus the cruise mean, for W), standardised with the tile
  matrix's KO SDs, projected on the axis loadings, in SD of the tile scores.
- Wilcoxon signed-rank test.
- KO-level effect: the mean paired difference per KO, divided by the KO's within-cruise
  residual SD. It is compared with the axis loadings by cosine, with a null of 999
  random sign flips of whole pairs.

**A16N: Protocol A minus partner** (16 pairs; 12 for the axis scores, which need both
samples in a modelled tile):

| Variable | Mean difference | SD | Same sign | Wilcoxon p |
|---|---:|---:|---:|---:|
| Archaea CLR | +0.381 | 0.351 | 14/16 | 0.002 |
| Bacteria CLR | −0.189 | 0.252 | 10/16 | 0.016 |
| Eukaryota CLR | −0.156 | 0.280 | 12/16 | 0.051 |
| bacteria vs eukaryota (logit) | −0.034 | 0.406 | 8/16 | 0.86 |
| unclassified logit | +0.023 | 0.482 | 9/16 | 0.86 |
| log10 total KO counts | −0.175 | 0.336 | 11/16 | 0.065 |
| F1 (SD) | +0.241 | 0.519 | 8/12 | 0.15 |
| F2 (SD) | −0.413 | 0.682 | 11/12 | 0.016 |
| F3 (SD) | +0.773 | 0.444 | 11/12 | 0.001 |
| W1 (SD) | +0.237 | 0.647 | 8/12 | 0.30 |
| W2 (SD) | −0.640 | 0.766 | 10/12 | 0.012 |
| W3 (SD) | +0.780 | 0.783 | 10/12 | 0.005 |

**KO-level effect.**
- RMS 0.46 residual SD per KO; 16.0% of KOs have |t| > 3.
- Cosine with the loadings: F1 +0.09, F2 −0.31, F3 +0.58, F4 +0.19, F5 −0.02,
  W1 +0.02, W2 −0.42, W3 +0.49, W4 −0.30. The sign-flip p for |cos F1| is 0.74.
- Most raised in Protocol A: K11779, K09009, K02007, K12234, K12369, K01926.
- Most lowered: K08364, K16363, K00854, K21572, K02652, K00284.

**Reading.**
- **The protocol switch is a large technical effect**: about half a residual SD per KO.
- **It loads on axis 3 and partly on axis 2, not on axis 1.**
  - Protocol A raises the Archaea share (+0.38 CLR, 14 of 16 pairs).
  - It scores +0.78 SD on W3, and its KO effect lines up with F3 (cosine 0.58).
  - Part of the W3 shift could come from its lower KO counts (−0.175 log10). That is
    only about half the shift expected from the R33 depth slope, so the extraction
    itself also contributes. The Archaea shift comes from the reads, not from the depth
    effect.
  - It lowers W2 (−0.64 SD).
- **It doesn't reproduce the era offset of axis 1.** Cosine with F1 is 0.09, and
  bacteria vs eukaryota doesn't move (−0.03).
  - If Protocol A is the old protocol, the old→new switch in preservation and
    extraction doesn't explain the axis-1 era offset. The sequencing platform and the
    library prep, which changed at the same time, become the leading candidates.
  - If Protocol A is something else, this test says nothing about the era.

**Consequences for the model.**
- Protocol A samples, and P02 Shield samples once they're counted, should be modelled
  as a technical treatment: a covariate, or kept out of the tile averages. They aren't
  replicates.
- The A16N Protocol A effect is large enough that a model trained on mixed protocols
  learns some protocol noise. Remaining checks:
  - drop the Protocol A samples
  - refit
  - check whether A16N's W2 / W3 scatter shrinks

### 4. Nutricline and deep chlorophyll maximum (CCHDO profiles)

**Data.** CCHDO bottle files (nitrate) and CTD files (2 dbar, fluorescence) for 11 of
the 13 modelled cruises. GOMECC4 and HLY2202 have no CCHDO files used here.

| Cruise | Casts | With nitracline | With DCM (source) | With CTD mixed layer | Bottles ≤ 200 m (median) |
|---|---:|---:|---|---:|---:|
| AMT28 | 63 | 60 | 62 (CTD) | 62 | 21 |
| I07N | 126 | 123 | 124 (CTD) | 124 | 5 |
| I09N | 116 | 113 | 116 (CTD) | 116 | 7 |
| P18 | 213 | 212 | 188 (CTD) | 212 | 5 |
| C13.5/A13.5 | 8 | 8 | 8 (CTD) | 8 | 6 |
| P02-2022, leg 1 | 149 | 116 | 149 (148 CTD, 1 bottle) | 148 | 7 |
| P02-2022, leg 2 | 127 | 93 | 127 (125 CTD, 2 bottle) | 125 | 7 |
| A16N-2023, leg 1 | 75 | 75 | 0 (no fluorescence) | 73 | 6 |
| A16N-2023, leg 2 | 74 | 73 | 0 | 69 | 8 |
| I05-2023 | 232 | 192 | 232 (231 CTD, 1 bottle) | 231 | 9 |
| A13.5-2024 | 112 | 112 | 0 (no CTD file) | 0 | 5 |
| I08S-2024 | 75 | 66 | 75 (CTD) | 74 | 9 |
| I09N-2025 | 114 | 97 | 114 (CTD) | 114 | 8 |

**Per-cast quantities** (`cchdo_profiles.py`):
- nitracline: the depth where NO3 first reaches 1 µmol/kg, interpolated between
  bottles; 0 if the surface is already ≥ 1
- the depth of the steepest NO3 increase (0–400 m)
- DCM depth: the fluorescence maximum in 0–300 m, after 2-m binning and a 5-bin
  running median
- DCM / surface (0–10 m) fluorescence ratio. Fluorometers aren't cross-calibrated, so
  only depths and ratios are used.
- CTD mixed layer: 0.2 °C below the 10 m temperature
- surface nitrate

**Per-cruise medians:**

| Cruise | Surface NO3 | Nitracline (1 µmol/kg) | Steepest gradient | DCM | DCM ratio | CTD mixed layer |
|---|---:|---:|---:|---:|---:|---:|
| A13.5-2024 | 0.0 | 42.5 | 92.0 | — | — | — |
| A16N-2023 | 0.0 | 44.8 | 84.8 | — | — | 31.0 |
| AMT28 | 0.0 | 87.9 | 126.7 | 73.0 | 5.7 | 48.0 |
| C13.5/A13.5 | 0.0 | 89.8 | 124.3 | 97.0 | 3.2 | 35.0 |
| I05-2023 | 0.1 | 106.1 | 152.4 | 50.0 | 1.4 | 65.0 |
| I07N | 0.0 | 53.8 | 97.9 | 68.0 | 1.4 | 43.0 |
| I08S-2024 | 22.5 | 0.0 | 86.9 | 49.0 | 1.4 | 59.0 |
| I09N | 0.0 | 63.5 | 77.8 | 75.0 | 4.2 | 27.0 |
| I09N-2025 | 0.0 | 67.5 | 89.2 | 75.0 | 2.8 | 25.0 |
| P02-2022 | 0.0 | 110.8 | 162.9 | 93.0 | 2.8 | 21.0 |
| P18 | 6.4 | 0.0 | 99.6 | 64.0 | 1.9 | 35.0 |

**Matching.** Each tile is matched to the nearest cast of its own cruise within 100 km
and 3 days: 820 of 1,216 tiles, median 0 km and 0.07 days, because most samples were
taken at stations. A stricter 30 km / 1 day match gives 792 tiles and the same results.
C13.5 (12 tiles) and I08S (9) contribute few.

**Spearman r with the axes**, all tiles / within cruise (cruises with the pooled sign);
100 km / 3 d match:

| Profile feature (tiles) | F1 | F2 | W1 | W2 | W3 |
|---|---|---|---|---|---|
| nitracline depth, 1 µmol/kg (756) | −0.06 / +0.11 (6/9) | +0.01 / −0.01 (4/9) | +0.09 / +0.10 (6/9) | −0.05 / −0.04 (5/9) | +0.02 / −0.01 (3/9) |
| nitracline, steepest gradient (756) | −0.17 / −0.00 (3/9) | −0.04 / −0.06 (6/9) | −0.02 / −0.01 (3/9) | −0.09 / −0.09 (6/9) | +0.07 / +0.04 (6/9) |
| DCM depth (618) | +0.19 / +0.04 (5/7) | −0.15 / +0.08 (4/7) | +0.05 / +0.05 (5/7) | +0.05 / +0.06 (4/7) | −0.01 / −0.01 (4/7) |
| DCM / surface fluorescence, log (618) | +0.25 / +0.19 (5/7) | −0.12 / +0.13 (5/7) | +0.14 / +0.19 (5/7) | +0.05 / +0.09 (4/7) | −0.04 / −0.06 (4/7) |
| CTD mixed layer (727) | −0.22 / −0.15 (6/8) | +0.19 / +0.03 (4/8) | −0.14 / −0.15 (6/8) | +0.03 / +0.04 (5/8) | +0.04 / +0.04 (6/8) |
| nitracline − mixed layer (667) | +0.06 / +0.20 (6/8) | −0.08 / −0.02 (4/8) | +0.18 / +0.19 (6/8) | −0.05 / −0.06 (4/8) | +0.01 / −0.02 (4/8) |
| surface nitrate, log (746) | −0.29 / −0.24 (5/9) | −0.10 / −0.15 (6/9) | −0.19 / −0.21 (5/9) | −0.10 / −0.12 (5/9) | +0.14 / +0.15 (8/9) |

**W2 from the profile, cross-validated.** Within-cruise ridge, cruise held out, on the
756 matched tiles:

| Predictors | R2 |
|---|---:|
| PISCES chlorophyll | −0.025 |
| nitracline | −0.010 |
| chlorophyll + nitracline | −0.036 |
| chlorophyll + nitracline + CTD mixed layer | −0.036 |

With the 30 km / 1 d match: −0.031, −0.013, −0.042, −0.043.

**Reading.**
- **The nutricline and DCM don't explain axis 2.** Every profile feature has
  within-cruise |r| ≤ 0.09 with W2, and adding them gives no cross-validated R2.
- Together with R34's smooth, roughly 560 km structure, this means the missing driver
  of the phototroph axis is large-scale but isn't the vertical nutrient structure at
  the station. Remaining candidates:
  - the recent chlorophyll history: bloom phase (section 5)
  - light: cloud and PAR history, mixed-layer light
  - grazing and viral lysis, which no field resolves
- **Axis 1 has weak but consistent profile links.** Its bacterial pole sits in
  stratified, oligotrophic water: a stronger DCM peak (W1 r +0.19), a nitracline well
  below the mixed layer (+0.19), lower surface nitrate (−0.21) and a shallower mixed
  layer (−0.15), each in 5–6 of 8–9 cruises. That is ecologically sensible
  (eukaryotes favoured where nitrate reaches the surface), but small next to its
  protocol / era offset.

### 5. Chlorophyll history: job for Morrill

**What it does.** `transfer/scripts/extract_chl_history.py`, run by
`transfer/cluster/chl_history.sbatch` (service partition, Copernicus Marine login):
- Data: gap-free daily GlobColour L4 chlorophyll at 4 km. It tries multi-year, then
  interim, then near-real-time, which covers 2011–2025.
- Per sample: the box mean (±0.125°) for each of the 29 days from day −28 to the
  sampling day.
- Summaries:
  - chl_day0, chl_mean28, chl_max28
  - log10 slopes over 14 and 28 days, per 10 days: > 0 growing, < 0 decaying
  - chl_peak_ratio28 = log10(28-day max / day 0); 0 = at the peak, > 0 = after a peak
  - days since the maximum
  - CV
- Resumable.
- Checked offline on synthetic series:
  - a growing bloom gives slope14 +0.18
  - a bloom past its peak gives slope14 −0.53 and peak ratio 0.70

**Run on Morrill:**

    cd /scratch/morrill/users/$USER/biogoship_granite_kegg && git pull
    sbatch --account=191001-364393 transfer/cluster/chl_history.sbatch

**Bring back:** `transfer/results/chl_history/sample_chl_history.csv` and
`provenance_chl_history.json`.

### Summary
1. **Axis-1 offsets follow extraction workflow / platform, not sampling year.**
   - HLY2202 (2022, old enzymatic extraction, OBI_2631) sits with the old cruises.
   - Within the enzymatic cruises there is no year trend (r +0.09).
   - Extraction and platform can't be separated with these cruises.
   - The A16N Protocol A pairs don't reproduce the axis-1 shift. If Protocol A is the
     old protocol, that points to sequencing platform / library prep.
2. **The A16N Protocol A vs standard pairs differ technically, and strongly.**
   - +0.38 Archaea CLR; +0.78 SD on axis 3; −0.64 SD on axis 2; 16% of KOs at
     |t| > 3.
   - These are technical treatments. They should not be averaged as replicates in
     the tiles.
   - The P02 Shield pairs can't be tested: their 2023 libraries are not in the KO
     table.
3. **Nutricline and DCM don't explain axis 2** (|r| ≤ 0.09; no CV gain). Axis 1 has
   weak stratification links (r ≈ 0.2).
4. **The chlorophyll-history job is ready to run on Morrill.**

**Open questions for the Bio-GO-SHIP team:**
- What is A16N "Protocol A"?
- Where are the P02 HWKVFDSX5 libraries (Shield pairs)? Their counts are missing from
  the KO table.
- Library prep kit per sequencing run.

## R36. Protocol A out of the real model; axis 2 is workflow-independent; C13.5 v1/v2 are duplicates; sequencer and library prep (2026-11-19)

**Decisions.**
- **A16N Protocol A samples are dropped from the real model**, and their Protocol B
  partners are kept.
  - They stay available for testing, in `data/targets/ko_clr_incl_protocolA.csv`.
  - Defined in the A16N 2023 leg 1 cruise report, section 6.13.3, p. 57. Both are
    0.22 µm Sterivex, −80 °C.
    - Protocol A: a duplicate underway DNA sample, no bashing beads, 1,600 µL lysis
      buffer, about every 3 days.
    - Protocol B: the standard underway sample, Zymo ZR BashingBeads, 1,000 µL DNA/RNA
      Shield.
- **Other multi-sample tiles stay averaged.** These are:
  - C13.5 v1/v2 (tested below)
  - GOMECC4 A/B/C
  - HLY2202 re-sequencing: 50 samples on two NovaSeq 6000 runs, HFJLFDSX7 and
    HMY7WDSX7
  - co-located samples

**Implementation.**
- `data/manifests/excluded_samples.csv` (tracked; copy in `transfer/results/protocol/`)
  lists the 17 Protocol A sample keys with the reason.
  - They are the A16N samples whose workbook note says "Protocol A": T-008 … T-121 and
    UW944, all inline.
- `scripts/build_ko_targets.py` has a new `--exclude-samples` option. It reads that
  file by default; `none` keeps all samples.
- The earlier targets are kept as `data/targets/ko_clr_incl_protocolA.csv` and
  `ko_samples_incl_protocolA.csv`.
- Check: rebuilding with `--exclude-samples none` reproduces them exactly (2,007
  samples × 3,921 KOs, max |difference| 0).
- New targets: 1,990 samples. Leaving out the 17 libraries lets 12 more KOs pass the 90%
  prevalence filter, so there are **3,933 KOs**; none are lost. The model then has
  **1,215 tiles**: one tile held only a Protocol A sample.
- **R33–R35 used the earlier targets.** `residual_axes_followup.py` now defaults to
  them, and R33's command is updated to pass them.

**Scripts.**
- `scripts/residual_modules.py --out-dir transfer/results/residual_modules_noA`: the
  refit.
- `transfer/scripts/axis2_protocol_check.py`, written to `transfer/results/axis2_protocol/`:
  refit vs R33, the protocol-group PCA and the paired tests.
- `transfer/scripts/cruise_methods_table.py`, now with library prep, sequencer and
  flowcells.

To reproduce, run from the repository root:

    ~/miniforge3/bin/python scripts/build_ko_targets.py        # real targets, Protocol A excluded
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python scripts/residual_modules.py --out-dir transfer/results/residual_modules_noA
    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/axis2_protocol_check.py
    ~/miniforge3/bin/python transfer/scripts/cruise_methods_table.py

### 1. Residual PCA without Protocol A

**Baseline.** Median KO R2 is 0.103. Residual share is 83.8%, 20.1% of it between
cruises.

**Loadings vs R33**, on the 3,921 shared KOs:

| Component | Cosine with R33 | Share (R33) | Split-half (R33) |
|---|---:|---|---|
| F1 | 1.000 | 19.4% (19.4%) | 0.90 (0.90) |
| F2 | 1.000 | 9.9% (9.9%) | 0.80 (0.77) |
| F3 | 1.000 | 6.8% (6.8%) | 0.55 (0.65) |
| F4 | 0.998 | 5.1% (5.1%) | 0.37 (0.58) |
| F5 | 0.999 | 4.7% (4.7%) | 0.56 (0.54) |
| W1 | 1.000 | 17.3% (17.2%) | 0.95 (0.95) |
| W2 | 1.000 | 9.8% (9.8%) | 0.83 (0.83) |
| W3 | 1.000 | 5.8% (5.8%) | 0.65 (0.66) |
| W4 | 1.000 | 4.6% (4.7%) | 0.59 (0.56) |
| W5 | 1.000 | 3.3% (3.3%) | 0.49 (0.41) |

The F3–F5 split-half values move because a different random split sequence follows from
the changed tile set.

**Interpretation.** Axis 2 comes back with identical loadings. But the refit changes
only the 12 tiles where A and B were averaged, plus the A-only tiles, so on its own it
cannot separate biology from protocol. The test that can is next.

### 2. Is axis 2 workflow-dependent? PCA within each workflow

**Method.**
- The within-cruise residuals of the new real model are PCA'd separately in two
  groups:
  - old workflow: enzymatic extraction, lysis buffer, Nextera, NovaSeq 6000. I09N,
    P18, I07N, AMT28, C13.5/A13.5 and HLY2202; HLY used DNA Flex prep.
  - new workflow: bead beating, DNA/RNA Shield, "LT4" flowcells. A16N, I05, A13.5,
    I08S and GOMECC4; GOMECC4 used beads and no preservative.
  - P02 and I09N-2025 are left out: their workflow is not verified.
- Reference: 200 random splits of the same 11 cruises into 6 + 5.

| Axis | Old workflow: rank (share), cos with full fit | New workflow: rank (share), cos | Old vs new congruence | Random 6+5 splits: median (5th pct) | Random splits lower |
|---|---|---|---:|---|---:|
| W1 | 1 (16.2%), 0.97 | 1 (18.5%), 0.98 | 0.91 | 0.94 (0.84) | 25% |
| **W2** | **2 (9.7%), 0.95** | **2 (11.4%), 0.94** | **0.79** | **0.81 (0.62)** | **38%** |
| W3 | 3 (6.6%), 0.95 | 4 (5.4%), 0.92 | 0.78 | 0.68 (0.47) | 82% |

**Reading.**
- **Axis 2 is the same in both workflows.** It is the second component in each, at
  9.7% and 11.4% of the variance.
- Its loadings match the full fit at cosine 0.95 and 0.94.
- The old and new versions agree as well as random cruise splits do: 0.79 against a
  median of 0.81.
- So the within-cruise phototroph axis doesn't depend on extraction, preservation,
  library prep or sequencer. **It is biology**, or at least nothing about the workflow.
- The same holds for W1 and W3.
- The R34 cruise-level offsets of axis 1 are a separate question; the next section
  takes it up.

### 3. Paired samples: A16N Protocol A vs B, and C13.5 v2 vs v1

**Method.**
- Data: the targets that keep every sample, on the 3,921 shared KOs.
- Pair difference: the difference of CLRs, standardised by the within-cruise residual
  SDs of the real model, then projected on its loadings (in SD of the tile scores).
  Both samples sit in one tile, so the tile prediction cancels.
- C13.5 pairs: the 29 v1/v2 sample pairs share time and position and come from
  underway collection on run HKNJLDRXX. They have separate BioSamples and SRA runs. 24
  pairs have both libraries.
- **The records show no protocol difference between v1 and v2.** The FAIRe workbook
  holds a single row per pair, with serial number "…v1", and has no notes.

| Variable | Protocol A − B (16 pairs): mean (SD; same sign; Wilcoxon p) | C13.5 v2 − v1 (24 pairs) |
|---|---|---|
| F1 | +0.069 (0.559; 9/16; 0.67) | −0.058 (0.368; 13/24; 0.57) |
| F2 | −0.418 (0.620; 14/16; 0.009) | −0.149 (0.601; 16/24; 0.17) |
| F3 | +0.899 (0.524; 15/16; < 0.001) | +0.205 (0.470; 17/24; 0.016) |
| W1 | +0.027 (0.691; 8/16; 0.98) | −0.101 (0.495; 14/24; 0.41) |
| W2 | −0.628 (0.719; 13/16; 0.003) | −0.180 (0.565; 16/24; 0.12) |
| W3 | +0.935 (0.824; 14/16; < 0.001) | +0.069 (0.828; 14/24; 0.49) |
| Archaea CLR | +0.381 (0.351; 14/16; 0.002) | −0.043 (0.153; 15/24; 0.25) |
| Bacteria CLR | −0.189 (0.252; 10/16; 0.016) | −0.120 (0.176; 17/24; 0.005) |
| Eukaryota CLR | −0.156 (0.280; 12/16; 0.051) | +0.011 (0.241; 11/24; 0.51) |

**KO-level systematic effect.** The RMS over KOs of the mean paired difference, against
its sign-flip null (999 flips of whole pairs).

| Contrast | RMS effect | Sign-flip null median | p | KOs with \|t\| > 3 |
|---|---:|---:|---:|---:|
| A16N Protocol A − B | 0.463 | 0.260 | 0.001 | 16.0% |
| C13.5 v2 − v1 | 0.210 | 0.174 | 0.070 | 1.5% |

- Cosine of the A − B effect with the loadings: F1 +0.09, F2 −0.32, F3 +0.58, W1 +0.02,
  W2 −0.42, W3 +0.49.
- For C13.5 v2 − v1: F1 −0.13, F2 −0.25, F3 +0.30, W1 −0.20, W2 −0.27, W3 +0.08.

**C13.5 v2 vs v1, sequencing depth.** log10 total counts differ by −0.063 on average,
12 of 24 negative, p 0.47. The Bacteria difference does not track the depth difference
(Spearman −0.12, p 0.58).

**Reading.**
- **C13.5 v1/v2 behave as duplicates, not a protocol contrast.**
  - The KO-level systematic effect isn't significant (p 0.07).
  - Strongly shifted KOs are 1.5% for C13.5 against 16% for Protocol A.
  - No axis score differs, apart from a small F3 shift.
  - One small consistent difference remains: v2 has 0.12 lower Bacteria CLR (17 of 24
    pairs). It isn't a depth effect. It is a third of the size of the Protocol A archaeal
    shift, and it could come from filter order or extraction batch.
  - **They stay averaged as replicates.**
- **Protocol A vs B (no beads + lysis buffer vs bead beating + Shield)** shifts axes 2
  and 3 and the Archaea share. It doesn't shift axis 1 (W1 +0.03; F1 +0.07).
  - The A/B pairs span both A16N LT4 runs: 8 A samples on 22NCJYLT4, 9 on 23LT2GLT4.
  - Whether each A and its B partner ran on the same flowcell isn't checked here.
  - As far as recorded, the within-pair contrast is preservation plus bead beating. The
    bench extraction after that isn't documented for either protocol.
  - That technical push along axis 2, −0.6 SD per sample, is real. But within cruises
    axis 2 can't come from it: every other sample in a cruise used one protocol, and
    section 2 shows the same axis in both workflows.

### 4. Sequencer and library preparation; what is left for the axis-1 offsets

**Library prep** (published run table, supplied 2026-11-19):
- All six runs of the original 2011–2020 collection used the same custom Nextera-style
  prep: Illumina Tagment DNA Enzyme and Buffer kit 20034197, custom 8-bp dual indexes,
  modified low-input tagmentation.
- Runs:

| Published run | Cruises | Sequencer | Libraries |
|---|---|---|---:|
| 1 | I09N | HiSeq | 24 |
| 2 | I09N, AE1319, BVAL46, NH1418 | HiSeq | 30 |
| 3 | I09N | NovaSeq | 215 |
| 4 | P18, AMT28, AE1319, BVAL46, NH1418 | NovaSeq | 203 |
| 5 | I07N | NovaSeq | 248 |
| 6 | C13.5 | NovaSeq | 251 |

- HLY2202: Nextera DNA Flex, low-volume, custom dual indexes, NovaSeq 6000.
- P02, A16N, I05, A13.5, I08S, I09N-2025: library prep and sequencer not verified.

**Matched to the FAIRe flowcells.** The coded columns are now in `cruise_methods.csv`:
library_prep, sequencer, flowcells, sequencer_coded.
- HiSeq = the "…BBXX" flowcells: HKMFTBBXX for I09N-2016, and HWKN2BBXX for I09N,
  AE1319, BVAL46 and NH1418.
- NovaSeq 6000: HHGJFDSXX, HMMFVDSXX, H357JDSXY, HKNJLDRXX, and HLY's HFJLFDSX7 /
  HMY7WDSX7.
- The 2023–2025 cruises are on "…LT4" flowcells: 22NCJYLT4, 23LT2GLT4, 22YJMLLT4. That
  is the NovaSeq X-series ID format, inferred, **not verified**.
- **No HiSeq library is in a modelled tile.** I09N-2016's 54 tiles are all NovaSeq
  (HHGJFDSXX), and AE1319, BVAL46 and NH1418 predate Sentinel-3 matching. So HiSeq
  vs NovaSeq doesn't enter the model.
- There is no direct sequencer contrast on the same DNA.
  - P02's 2023 NovaSeq 6000 run (HWKVFDSX5, 58 libraries) holds different samples from
    its 2025 LT4 run, and it isn't in the KO table.
  - HLY's two runs are both NovaSeq 6000.

**What the axis-1 cruise offsets can still be.**
- Ruled out as the main driver:
  - **Library prep**, Nextera tagmentation vs DNA Flex: HLY2202 (DNA Flex) sits with
    the tagmentation cruises (F1 +0.78).
  - **Preservation plus bead beating**: the A16N A/B pairs show no axis-1 shift.
  - **Sampling year**, within the enzymatic cruises: R35.
- Still confounded with the offset:
  - the sequencing platform: NovaSeq 6000 vs the LT4 flowcells
  - the bench extraction kit and the library prep of the 2023–2025 cruises, both not
    verified
  - a real change between the regions and years sampled
- The decisive check would be any one of these:
  - the library prep / sequencer of the new cruises
  - the same DNA re-sequenced on both platforms
  - the P02 2023 NovaSeq 6000 libraries mapped through our pipeline (R32), then
    compared with neighbouring P02 samples on 22YJMLLT4. Not paired, but the same
    cruise and the same water masses.

## R37. Bloom phase: satellite chlorophyll history against the residual axes (2026-11-20)

**Question.** Is axis 2 (phototrophs) the bloom phase: whether chlorophyll was rising
or falling over the 2–4 weeks before sampling? Does chlorophyll history help predict
KOs?

**Data.** `transfer/results/chl_history/sample_chl_history.csv` and
`provenance_chl_history.json`, from Morrill job `transfer/cluster/chl_history.sbatch`
(script `transfer/scripts/extract_chl_history.py`, R35 §5).
- Product: GlobColour gap-free L4 daily chlorophyll, 4 km, dataset
  `cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D` (1997-09-04 to
  2026-11-11). The multi-year product covered every sample, so the interim and NRT
  products weren't needed.
- 2,014 samples, including P04E-2026, which is not in the model. Each has 26–29 of
  the 29 days, median 29.
- Box mean ±0.125°, days −28 to 0.

**Sample distributions:**

| Feature | Median | IQR | Range |
|---|---:|---|---|
| chl day 0 (mg m⁻³) | 0.099 | 0.064–0.142 | 0.013–5.81 |
| slope 14 d (log10 / 10 d) | −0.003 | −0.060 to +0.059 | −0.68 to +0.82 |
| slope 28 d | −0.006 | −0.037 to +0.025 | −0.24 to +0.29 |
| peak ratio 28 d (log10 max/day 0) | 0.130 | 0.064–0.220 | 0–1.09 |
| days since 28-day max | 15 | 5–24 | 0–28 |

**Analysis.** Script `transfer/scripts/chl_history_test.py`; outputs in
`transfer/results/chl_history_test/`. It runs on the real model: Protocol A excluded,
1,215 tiles, 3,933 KOs, with the residual PCA of R36
(`transfer/results/residual_modules_noA/`). Every tile has history.

To reproduce, run from the repository root:

    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/chl_history_test.py

**1. Spearman r**, all tiles / within cruise (cruises with the pooled sign, of 13):

| Feature | F1 | F2 | W1 | W2 | W3 |
|---|---|---|---|---|---|
| chlorophyll, day 0 (log) | −0.10 / −0.05 (8) | +0.02 / −0.09 (7) | −0.06 / −0.06 (8) | −0.05 / −0.06 (6) | −0.08 / −0.09 (10) |
| chlorophyll, 28-day mean (log) | −0.06 / +0.02 (7) | +0.03 / −0.08 (6) | +0.01 / +0.01 (6) | −0.05 / −0.06 (6) | −0.07 / −0.07 (8) |
| chlorophyll, 28-day max (log) | −0.06 / +0.02 (7) | +0.04 / −0.09 (6) | +0.01 / +0.02 (5) | −0.05 / −0.07 (6) | −0.07 / −0.08 (7) |
| slope, 14 days | −0.07 / −0.12 (9) | +0.05 / +0.00 (6) | −0.12 / −0.12 (9) | +0.02 / +0.02 (6) | −0.01 / −0.02 (8) |
| slope, 28 days | −0.14 / −0.11 (7) | +0.09 / +0.02 (7) | −0.12 / −0.12 (7) | +0.04 / +0.04 (8) | +0.08 / +0.08 (7) |
| past the 28-day peak | +0.07 / +0.14 (8) | +0.02 / +0.01 (7) | +0.14 / +0.14 (8) | −0.01 / −0.02 (7) | −0.01 / −0.01 (7) |
| days since the 28-day max | +0.09 / +0.05 (9) | −0.04 / +0.03 (7) | +0.08 / +0.08 (8) | +0.02 / +0.01 (8) | −0.04 / −0.04 (8) |
| CV over 28 days | −0.07 / −0.00 (4) | +0.02 / −0.08 (6) | +0.00 / +0.00 (9) | −0.05 / −0.06 (7) | −0.01 / −0.03 (7) |
| PISCES chlorophyll, day 0 (R25) | −0.21 / −0.08 (8) | +0.02 / −0.07 (8) | −0.07 / −0.07 (9) | −0.02 / −0.03 (7) | −0.04 / −0.06 (10) |

**2. W2 from the features.** Within-cruise ridge, cruise held out, 1,215 tiles:

| Predictors | CV R2 |
|---|---:|
| satellite chlorophyll, day 0 | −0.015 |
| PISCES chlorophyll, day 0 | −0.018 |
| slopes + peak ratio + days since max | −0.016 |
| all 8 satellite history features | −0.038 |

**3. As a KO predictor.** The baseline (32 PCs + temperature spline, cruise held out)
plus one feature's spline, on 1,215 tiles. Baseline median KO R2 is 0.103.

| Added | Median KO R2 | Median gain | KOs gaining > 0.02 | KOs losing > 0.02 |
|---|---:|---:|---:|---:|
| chlorophyll, day 0 (log) | 0.089 | −0.0118 | 8% | 35% |
| chlorophyll, 28-day mean (log) | 0.080 | −0.0188 | 9% | 48% |
| slope, 14 days | 0.098 | −0.0035 | 7% | 12% |
| slope, 28 days | 0.101 | −0.0052 | 13% | 20% |
| past the 28-day peak | 0.106 | −0.0019 | 14% | 9% |

**4. Spatial scale.** Within-cruise correlation of pairs by distance (km):

| Series | 0–15 | 15–40 | 40–75 | 75–150 | 150–300 | 300–600 | 600–1200 | 1200–2400 | 2400–5000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| slope, 14 days | 0.83 | 0.69 | 0.68 | 0.41 | 0.28 | 0.18 | −0.03 | −0.12 | −0.01 |
| past the 28-day peak | 0.82 | 0.76 | 0.69 | 0.46 | 0.28 | 0.06 | 0.00 | −0.06 | −0.09 |
| satellite chlorophyll, day 0 | 0.99 | 0.94 | 0.92 | 0.81 | 0.72 | 0.50 | 0.22 | −0.17 | −0.29 |
| W2 (real model) | 0.56 | 0.49 | 0.49 | 0.39 | 0.32 | 0.24 | 0.15 | −0.02 | −0.09 |

**Reading.**
- **Bloom phase doesn't explain axis 2.**
  - Every chlorophyll-history feature has within-cruise |r| ≤ 0.07 with W2.
  - The ridge gives no cross-validated R2.
  - The scales don't match either. Bloom phase decorrelates faster than W2: the
    14-day slope falls to 0.18 by 300–600 km and the peak ratio to 0.06, while W2
    keeps 0.24 there and 0.15 at 600–1,200 km.
- **Axis 1 has a weak bloom-phase link.** Its bacterial pole lies slightly past a
  chlorophyll peak and on decaying slopes: W1 +0.14 with the peak ratio, −0.12 with
  the 14- and 28-day slopes, in 7–9 of 13 cruises. That fits the R35 stratification
  link: bacteria are favoured after the bloom. It is small.
- **Chlorophyll history doesn't improve KO prediction.**
  - Day-0 or mean chlorophyll lowers median KO R2 by 0.012–0.019. Granite already sees
    the ocean colour, so a second, noisier copy overfits.
  - The slope and peak-ratio splines are neutral to slightly negative. The best is
    the peak ratio at median R2 0.106, but the median per-KO gain is −0.002.
  - They are not worth adding.
- **Axis 2 after R34–R37.** A smooth, roughly 560 km, workflow-independent phototroph
  axis that none of these explains:
  - surface or model fields
  - the vertical nutrient structure
  - the bloom phase
  - the protocols
  - sampling time
- Candidates not yet tested:
  - light history: cloud cover / PAR over the previous weeks, and mixed-layer light
  - temperature history beyond 14 days
  - water-mass origin: Lagrangian back-trajectory age and source region; the
    backtrack job exists, cluster/backtrack.sbatch
  - top-down control: grazing and viruses, which no product resolves

## R38. Is axis 2 ridge shrinkage rather than a missing driver? No (2026-11-20)

**Hypothesis** (suggested 2026-11-20).
- Ridge predictions might be compressed. With a calibration slope b > 1, the residual
  ≈ (b − 1) × prediction + noise.
- Then a residual axis would inherit the target's smooth spatial pattern without any
  missing variable, and new features couldn't explain it.

**Signatures that would show it.** Correlation with the *observed* values doesn't
discriminate, since residual = observed − predicted for any residual axis. Under
compression we would see all of these:
- calibration slopes > 1 for the axis-2 KOs
- the axis-2 score correlating **positively** with the predicted axis-2 signal
- the axis vanishing once the predictions are rescaled

**Method.** Script `transfer/scripts/axis2_shrinkage_check.py`; outputs in
`transfer/results/axis2_shrinkage/`. Real model: Protocol A excluded, 1,215 tiles,
3,933 KOs, baseline 32 PCs + temperature spline, cruise held out.
- W2 is recomputed: within-cruise component 2, 9.8%. Its cosine with the R36 W2 is
  1.000; split-half 0.84.
- **Slopes:** OLS of observed on cross-validated predicted, per KO, both centred per
  cruise ("within") or overall ("pooled").
  - KO sets: the 50 KOs loading most positively and most negatively on W2, the KEGG
    photosynthesis KOs (BRITE C "Photosynthesis…"), and all KOs.
  - Comparison: the median slope of KOs in the same R2 decile.
- **Composite:** the mean within-cruise z of the +50 KOs minus the −50 KOs, z-scored
  with the observed SDs. It is computed from the observed and from the predicted CLRs.
  - Also the Cyanobacteriota and Cryptophyta phylum CLRs, predicted by the same
    baseline.
- **Recalibration:**
  - For each held-out cruise and KO, the within-cruise slope is estimated on the other
    12 cruises' out-of-fold predictions and applied to the held-out cruise.
  - Then the within-cruise residual PCA is redone.
  - There is a mild leak: those predictions came from models that saw the held-out
    cruise. It can only help remove the axis.

To reproduce, run from the repository root:

    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/axis2_shrinkage_check.py

**1. Calibration slopes.** 1 = calibrated, > 1 = compressed, < 1 = over-dispersed.

| KO set | KOs | Median R2 | Median slope, within cruise | Median slope, pooled | Same-R2-decile KOs | Slope > 1.5 |
|---|---:|---:|---:|---:|---:|---:|
| W2 + pole (top 50) | 50 | 0.248 | 0.78 | 0.76 | 0.69 | 0% |
| W2 − pole (top 50) | 50 | 0.083 | 0.50 | 0.66 | 0.58 | 0% |
| KEGG photosynthesis KOs | 65 | 0.106 | 0.55 | 0.65 | 0.62 | 0% |
| all KOs | 3,933 | 0.103 | 0.59 | 0.67 | 0.58 | 0% |

**2. W2 score against the observed and predicted axis-2 signal**, within cruise
(Pearson r):

| Series | r(W2, observed) | r(W2, predicted) | r(observed, predicted) | SD predicted / observed | Calibration slope |
|---|---:|---:|---:|---:|---:|
| KO composite (+50 minus −50) | +0.75 | −0.20 | +0.38 | 0.62 | 0.61 |
| Cyanobacteriota CLR | +0.38 | −0.04 | +0.79 | 0.80 | 0.99 |
| Cryptophyta CLR | +0.56 | −0.03 | +0.43 | 0.64 | 0.67 |

**3. Recalibrated predictions:**

| Residuals | W2 rank (share) | Cosine with R36 W2 | Split-half | Photosynthesis z | Median within-cruise KO R2 |
|---|---|---:|---:|---:|---:|
| baseline | 2 (9.8%) | 1.000 | 0.84 | +12.8 | 0.051 |
| recalibrated | 2 (9.8%) | 0.985 | 0.84 | +12.5 | 0.085 |

Recalibrated vs baseline W2: loading cosine 0.985, score r 0.958.

**Reading.**
- **The predictions aren't compressed; they are over-dispersed.**
  - Every set has a median slope below 1. The axis-2 KOs have 0.78 / 0.50, and all KOs
    0.59.
  - No KO set has slopes above 1.5.
  - The cross-validated predictions overshoot: they vary more than the part of the
    observed variation they capture. The axis-2 KOs are no more attenuated than other
    KOs of the same R2.
- **The W2 score correlates negatively with the predicted signal (−0.20), not
  positively.**
  - With slope b < 1, the residual carries −(1 − b) × prediction, which is the
    opposite of the shrinkage signature.
  - A negative correlation with the predicted values therefore points to
    over-dispersed predictions, not to attenuation.
- **Rescaling the predictions doesn't remove axis 2.**
  - After per-KO recalibration, W2 is still component 2, at 9.8%, with the same
    loadings (cosine 0.985) and the same Photosynthesis enrichment (z +12.5).
- **The phototroph signal the model does capture is separate from W2.**
  - The baseline predicts Cyanobacteriota well within cruise: observed–predicted r
    0.79, slope 0.99.
  - W2 is uncorrelated with that prediction (−0.04) and tracks Cryptophyta (+0.56)
    and the KO composite (+0.75).
  - So **axis 2 is the part of the phototroph community orthogonal to what Granite and
    temperature see.** It is not an attenuated copy of what the model already knows.
- **Side finding for the real model.** Recalibrating raises the median within-cruise
  KO R2 from 0.051 to 0.085. The leak above flatters this number, but the
  over-dispersion is real.
  - The ridge alpha chosen per KO (RidgeCV, internal leave-one-out across tiles
    regardless of cruise) is too weak for transfer to a new cruise.
  - Choosing alpha by leave-cruise-out (GroupKFold) inside the training cruises, or a
    post-hoc shrink factor learned the same way, should improve held-out R2. This is
    worth testing on the main model.

## R39. Per-library sequencing metadata; platform contrasts along axis 1 (2026-11-20)

### 1. The library table

**Files.**
- `transfer/scripts/library_metadata.py` writes
  `transfer/results/library_metadata/library_metadata.csv` (3,222 rows) and
  `cruise_summary.csv`.
- ENA run records are cached in `transfer/data/ena_PRJNA656268_runs.tsv`.

**Rows.**
- Every MetaG library in the FAIRe workbook.
- Plus the KO-table libraries that have no workbook row:
  - 35 ERR accessions from another project
  - SRR12480044 (I7N2018_248), which isn't in the targets
- 3,045 of the 3,081 KO-table libraries match a FAIRe run row.

**Columns, with the source of each.** A per-row `source_of_information` column lists
them.

| Column | Source |
|---|---|
| lib_id, faire_lib_id, sample_id, samp_name, cruise | KO table + FAIRe; sample key via `build_ko_targets.library_samples` |
| in_ko_table, in_real_targets, excluded_protocol_a | KO table, `data/targets/ko_samples.csv`, `excluded_samples.csv` |
| seq_run_id, sequencing_date, flowcell_id, read_length | FAIRe `seq_run_id` (date_flowcell_MetaG_PE150; "20189999" = 2018, day unknown). Read length is checked against ENA bases per spot: 302 for 937 runs, 300 for the 24 HKMFTBBXX runs. |
| lane | library / file name (`_L00n`, `-Ln-`) where encoded. Not encoded for the 2016–2020 libraries. |
| instrument_id | FAIRe instrument (OBI term) |
| sequencer_model | (1) ENA/SRA run record; (2) the published cruise table; (3) inferred from the flowcell-ID format, flagged unverified. `sequencer_source` says which. |
| library_prep_kit, extraction, documentation | the cruise table supplied 2026-11-20 (below) |
| extraction_batch | FAIRe `extract_id` plate (+ `date_ext`). Recorded only for A13.5-2024 (2 plates), A16N (4), GOMECC4 (6) and I05 (4). |
| preservation_protocol | methods table (R35), per sample for A16N Protocol A/B and the P02 "shield" samples |
| total_reads, ena_read_count | FAIRe `input_read_count`; ENA `read_count` |
| ko_assigned_reads | the library's row sum in `ko_counts_cov50_all.tsv` |
| mapped_reads | **empty**: mapping (flagstat) counts aren't available to us |
| sra_run | FAIRe `associatedSequences` / ENA |

**Shortcuts tried.**
- **ENA/SRA run records.** Only the 2016–2020 cruises are accessioned:
  PRJNA656268, 971 WGS runs, all matched.
  - instrument_model gives Illumina NovaSeq 6000 for 917 runs and HiSeq 4000 for 54.
  - library_construction_protocol is empty for every run.
  - The 2022–2025 cruises have no accessions, so nothing exists for them in SRA.
- **The flowcell-ID convention, checked against ENA on all accessioned runs:**
  - "…BBXX" flowcells are always HiSeq 4000: HKMFTBBXX (24), HWKN2BBXX (30).
  - "…DSXX", "…DSXY" and "…DRXX" are always NovaSeq 6000: H357JDSXY, HHGJFDSXX,
    HKNJLDRXX, HMMFVDSXX (906).
  - So the flowcell suffix identifies the hardware family for these runs. The "…LT4"
    flowcells of 2023–2025 (22NCJYLT4, 23LT2GLT4, 22YJMLLT4) have the NovaSeq X-series
    format, but nothing can check that against a record. **Unverified.**
- **FASTQ headers.** The ENA copies are rewritten (`@SRR12479652.1 1/1`), so instrument,
  run, flowcell and lane can't be read from public files.
  - The original files of any 2023–2025 library would settle the sequencer. A NovaSeq X
    instrument ID starts with "LH".
  - This needs the Bio-GO-SHIP team's raw files. Our Morrill downloads are the external
    sets only.

**Per cruise** (libraries; in the KO table; in the real targets; flowcells; sequencer;
SRA runs; extraction batches):

| Cruise | Libraries | In KO table | In targets | Flowcells (libraries) | Sequencer (source) | SRA runs | Extraction batches |
|---|---:|---:|---:|---|---|---:|---:|
| BVAL46-2011 | 16 | 16 | 12 | HMMFVDSXX (10); HWKN2BBXX (6) | NovaSeq 6000 (10), HiSeq 4000 (6) (ENA) | 16 | — |
| AE1319-2013 | 18 | 18 | 13 | HWKN2BBXX (11); HMMFVDSXX (7) | HiSeq 4000 (11), NovaSeq 6000 (7) (ENA) | 18 | — |
| NH1418-2014 | 30 | 29 | 23 | HMMFVDSXX (19); HWKN2BBXX (11) | NovaSeq 6000 (19), HiSeq 4000 (10) (ENA); 1 published | 29 | — |
| I09N-2016 | 242 | 242 | 239 | HHGJFDSXX (215); HKMFTBBXX (24); HWKN2BBXX (3) | NovaSeq 6000 (215), HiSeq 4000 (27) (ENA) | 242 | — |
| P18-2016 | 104 | 103 | 102 | HMMFVDSXX (104) | NovaSeq 6000 (ENA) | 104 | — |
| I07N-2018 | 237 | 235 | 235 | H357JDSXY (237) | NovaSeq 6000 (ENA) | 237 | — |
| AMT28-2018 | 63 | 62 | 62 | HMMFVDSXX (63) | NovaSeq 6000 (ENA) | 63 | — |
| C13.5/A13.5-2020 | 251 | 249 | 249 | HKNJLDRXX (251) | NovaSeq 6000 (ENA) | 251 | — |
| HLY2202-2022 | 200 | 192 | 186 | HFJLFDSX7 (142); HMY7WDSX7 (50); HWKVFDSX5 (8) | NovaSeq 6000 (published) | 0 | — |
| P02-2022 | 391 | 333 | 123 | 22YJMLLT4 (333); HWKVFDSX5 (58, not in KO table) | LT4: NovaSeq X? (inferred); HWKVFDSX5: NovaSeq 6000? (inferred) | 0 | — |
| GOMECC4 | 472 | 458 | 172 | 23LT2GLT4 (472) | NovaSeq X? (inferred) | 0 | 6 |
| A16N-2023 | 345 | 341 | 203 | 22NCJYLT4 (279); 23LT2GLT4 (66) | NovaSeq X? (inferred) | 0 | 4 |
| I05-2023 | 312 | 308 | 160 | 23LT2GLT4 (312) | NovaSeq X? (inferred) | 0 | 4 |
| A13.5-2024 | 134 | 132 | 131 | 23LT2GLT4 (134) | NovaSeq X? (inferred) | 0 | 2 |
| I08S-2024 | 44 | 44 | 28 | 22YJMLLT4 (44) | NovaSeq X? (inferred) | 0 | — |
| I09N-2025 | 190 | 190 | 103 | 22YJMLLT4 (190) | NovaSeq X? (inferred) | 0 | — |
| InterCal (not a cruise) | 30 | 0 | 0 | HWKVFDSX5 (30) | NovaSeq 6000? (inferred) | 0 | — |
| no cruise (controls, blanks, mocks, unlinked) | 143 | 129 | 0 | mixed | — | 1 | 6 |

**Documented processing per cruise** (supplied 2026-11-20; it fills `extraction`,
`library_prep_kit` and `documentation`):

| Cruises / batch | DNA extraction | Library preparation | Sequencer | Confidence |
|---|---|---|---|---|
| Original I09N, P18, I07N, AMT28, C13.5 and earlier cruises | Enzymatic lysis + Zymo purification | Modified Nextera tagmentation, Illumina kit 20034197, custom dual indexes | HiSeq 4000 or NovaSeq 6000 | Published |
| HLY2202 (2022) | Original-style enzymatic lysis + Zymo purification | Nextera DNA Flex, low-volume protocol, custom dual indexes | NovaSeq 6000, PE150 | Published |
| P02W/P02E (2022) | P02W shipped to UC Irvine; final extraction needs confirmation | ? | 2023 NovaSeq 6000 run (HWKVFDSX5) for some libraries; accession mapping missing | Partial |
| A16N (2023) | Processing assigned to NOAA AOML / Northern Gulf Institute; lab workflow not verified | ? | ? | Collection documented |
| I05 (2023), A13.5 (2024), I08S (2024), I09N (2025) | ? | ? | ? | Collection documented |

### 2. Library-level platform contrasts along axis 1

**Method** (`transfer/scripts/sequencer_contrasts.py`; outputs `sequencer_contrasts.csv`
and `sequencer_contrasts_summary.txt` in `transfer/results/library_metadata/`):
- **Rarefying.** Every library is rarefied to 250,000 KO-assigned reads: one multinomial
  draw over all 18,836 KO features, seed 0. 2,943 libraries pass; 103 below that depth
  are left out.
  - This is needed. Without it, depth differences between runs dominate. HLY2202's
    second run, HMY7WDSX7, is a shallow top-up: median 1.2 M KO reads against 4.2 M.
    NH1418's HiSeq libraries are about 3× shallower than its NovaSeq ones.
  - The unrarefied first pass gave a same-DNA run-to-run SD of 7× the era offset. That
    was a depth artefact.
- **CLR.** Over the 3,933 real-model KOs, pseudocount 0.5 / N.
- **Projection.** Centred on the contrast's reference group, divided by the KO's pooled
  within-cruise SD of the observed tile CLR, projected on the unit F1 loading of the
  real model.
- **Scale ("era").** The same projection of cruise means of modelled surface libraries:
  NovaSeq 6000 cruises minus LT4 cruises = +14.96 projection units. Cruise means run from
  I05-2023 −51.1 to BVAL46 +31.9.
- **Domain CLRs** come from `data/targets/taxa/domain_counts.tsv`, not rarefied.
- Surface libraries only (≤ 10 m).

| Contrast | Libraries / pairs | Axis-1 (F1) difference, × era offset | Bacteria CLR diff | Eukaryota CLR diff |
|---|---:|---:|---:|---:|
| **a. HiSeq 4000 − NovaSeq 6000, different samples of one cruise** | | | | |
| I09N-2016: HiSeq (26) − NovaSeq (206) | 232 | +0.86 | +0.134 | +0.140 |
| AE1319-2013: HiSeq (6) − NovaSeq (7) | 13 | +1.76 | +0.476 | +0.009 |
| BVAL46-2011: HiSeq (2) − NovaSeq (10) | 12 | +0.89 | +0.608 | +0.430 |
| NH1418-2014: HiSeq (4) − NovaSeq (19) | 23 | +2.63 | +0.476 | +0.007 |
| pooled (weighted by the smaller group; label permutation within cruise, 9,999) | 280 | **+1.19, p < 0.001** | | |
| **b. Same DNA, two NovaSeq 6000 runs** (HLY2202, HMY7WDSX7 − HFJLFDSX7) | 21 pairs | mean −0.03 (SD of pair differences 1.20) | −0.003 | +0.097 |
| **c. Two libraries of one sample, same flowcell** (I09N-2025) | 3 pairs | mean +0.12 (SD 0.58) | −0.057 | −0.008 |
| **d. Two LT4 flowcells within A16N** (23LT2GLT4 (44) − 22NCJYLT4 (158), different samples) | 202 | +0.68 | +0.170 | −0.099 |

**Reading.**
- **Sequencing the same DNA again on the same platform doesn't move axis 1.**
  - HLY2202's 21 same-DNA pairs on two NovaSeq 6000 runs differ by −0.03× the era
    offset on average.
  - The per-library scatter (SD 1.2× era) is single-library noise at 250 k reads. It
    averages out over a cruise.
- **The platform does move it.**
  - Within each of the four early cruises, the HiSeq 4000 libraries sit further toward
    the bacterial (+) pole of axis 1 than the NovaSeq 6000 libraries: +0.86 to +2.63×
    the era offset, pooled +1.19×, p < 0.001.
  - Their Bacteria CLR is also higher, by 0.13–0.61.
  - Caveat: these are different samples. The HiSeq libraries were early batches, so
    station selection could add to the difference. The permutation treats libraries
    within a cruise as exchangeable.
  - The sign is the same in all four cruises, and the effect is as large as the era
    offset itself.
  - HiSeq 4000 (four-colour SBS) and NovaSeq 6000 (two-colour SBS) differ in chemistry.
    A shift from bacterial toward eukaryotic genes between them fits chemistry-dependent
    GC or coverage bias.
- **No HiSeq library is in a modelled tile** (R36), so this contrast doesn't affect the
  model directly. It shows that **axis 1 responds to the sequencing platform.**
  - The 2023–2025 cruises' "LT4" flowcells: if they are NovaSeq X, that is another
    chemistry (XLEAP-SBS).
  - That makes the platform change the leading explanation of the axis-1 era offset,
    ahead of the unverified extraction and library prep of those cruises.
  - It also fits R36: library prep (HLY's DNA Flex) and preservation / bead beating
    (A16N A/B) don't move axis 1.
- **The A16N two-flowcell difference** (+0.68×) is between different samples, and the
  R34 rotation test within A16N found no flowcell effect on W1 (p 0.58). It isn't
  evidence of a run effect.

**Next.**
1. Confirm the 2023–2025 sequencer: one original FASTQ header from any LT4 library, or
   the sequencing-facility record.
2. A GC check: if the axis-1 loadings follow the KOs' GC content, a platform GC bias is
   the mechanism. This needs per-KO GC from the reference genes (eggNOG table on
   Morrill).
3. If platform is confirmed, model it: era/platform as a covariate, or a
   batch-correction step for the KO targets.

## R40. Sequencer confirmed from a FASTQ header; the team's run list; axis 1 follows platform, not lab; GC job (2026-11-20)

### 1. New sources (in `transfer/`)

**`maybe_useful_README.txt`** (downloaded 2026-11-20): the team's list of metagenome
sequencing runs, giving library count, lab group, `seq_run_id`, cruises and run name.

| Libraries | Lab | seq_run_id | Cruises | Run name | Notes |
|---:|---|---|---|---|---|
| 24 | whoi | 20170925_HKMFTBBXX_MetaG_PE150 | I09N 2016 | larkin_1 | |
| 31 | whoi | 20189999_HWKN2BBXX_MetaG_PE150 | I09N 2016, AE1319 2013, BVAL46 2011, NH1418 2014 | larkin_2 | |
| 215 | whoi | 20190318_HHGJFDSXX_MetaG_PE150 | I09N 2016 | larkin_3 | |
| 203 | whoi | 20190730_HMMFVDSXX_MetaG_PE150 | P18 2016–2017, AMT-28 2018, AE1319 2013, BVAL46 2011, NH1418 2014 | larkin_4 | |
| 237 | whoi | 20200318_H357JDSXY_MetaG_PE150 | I07N 2018 | larkin_5 | |
| 251 | whoi | 20200909_HKNJLDRXX_MetaG_PE150 | C13.5 2020 | larkin_6 | |
| 82 | whoi | 20230419_HWKVFDSX5_MetaG_PE150 | ARC 2022, P02 2022, Intercal | 2023_mixed | |
| 143 | uci | 20230803_HFJLFDSX7_MetaG_PE150 | HEALY2202 2022 | mlbrock23080541 | |
| 50 | uci | 20231029_HMY7WDSX7_MetaG_PE150 | HEALY2202 2022 | tchungsw23103138 | resequenced samples from mlbrock23080541? |
| 288 | whoi | 20250405_22NCJYLT4_MetaG_PE150 | A16N 2023 | thompson_2025 | |
| 634 | uci | 20250730_22YJMLLT4_MetaG_PE150 | I08S 2024, I09N 2025, P02 2022, NP | BGS1:7 | |
| 813 | noaa | 20260406_23LT2GLT4_MetaG_PE150 | GSI052023, GOMECC4 | thompson_2026 | |
| 200 | whoi | 20260406_23LT2GLT4_MetaG_PE150 | GSA1352024, GSA16N2023 | thompson_2026 | |

- HMY7WDSX7 re-sequenced HFJLFDSX7 samples. R39 confirms this: 50 HLY samples are on
  both runs.
- One flowcell, 23LT2GLT4, was shared by NOAA (I05-2023, GOMECC4) and WHOI (A13.5-2024,
  A16N-2023).

**`xR058-L2-P01-TAAGGCGA-AGGCTTAG-READ1-Sequences.txt.gz`** (2.2 GB): an original READ1
file from run 22YJMLLT4, lane 2.
- First header: `@LH00626:65:22YJMLLT4:2:1101:2991:1056 1:N:0:TNAGGCGA+AGGCTTAG`.
  - Instrument **LH00626**: the "LH" prefix is the NovaSeq X series.
  - Run 65, flowcell 22YJMLLT4, lane 2.
  - The index pair matches the file name.
- All of the first 100,000 reads carry the same instrument, run, flowcell and lane.
- Their quality scores use four bins (in Phred+33):
  - `I` (Q40): 14,097,930
  - `9` (Q24): 619,474
  - `-` (Q12): 382,445
  - `#` (Q2): 151

  This is the NovaSeq X binning. NovaSeq 6000 bins as F / : / , / #, as in the ENA
  reads of R39.
- **22YJMLLT4 (UCI; I08S-2024, I09N-2025, P02-2022) is therefore a NovaSeq X run,
  confirmed.**
  - The other two "…LT4" flowcells (22NCJYLT4: A16N; 23LT2GLT4: A13.5, A16N, I05,
    GOMECC4) share the format but are still unconfirmed. (Both were confirmed NovaSeq X
    later, in R41.)
  - One header from each would settle them.
- The file is large and not needed further. It shouldn't stay in `transfer/` unless
  wanted.

**Library table updated** (`transfer/scripts/library_metadata.py`):
- New columns `seq_lab` and `run_name`, from the README.
  - 23LT2GLT4 libraries are split by name prefix: GSI05 / GOMECC to NOAA, GSA to WHOI.
- For 22YJMLLT4, `sequencer_model` is now "Illumina NovaSeq X", with the header as its
  source. Per flowcell:

| Flowcell | Sequencer | Source |
|---|---|---|
| HKMFTBBXX, HWKN2BBXX | Illumina HiSeq 4000 | ENA |
| HHGJFDSXX, HMMFVDSXX, H357JDSXY, HKNJLDRXX | Illumina NovaSeq 6000 | ENA |
| HFJLFDSX7, HMY7WDSX7, HWKVFDSX5 | NovaSeq 6000 | published (HLY) / flowcell format |
| 22YJMLLT4 | Illumina NovaSeq X | FASTQ header |
| 22NCJYLT4, 23LT2GLT4 | NovaSeq X? | flowcell format, unconfirmed |

### 2. Axis-1 cruise offsets: platform or lab?

**Data.** Real model, R36 scores. Per cruise: the mean F score in SD of the tile scores.
Output in `transfer/results/library_metadata/cruise_axis_by_lab_platform.csv`.

| Cruise | Lab | Platform | F1 | F2 | F3 |
|---|---|---|---:|---:|---:|
| HLY2202-2022 | uci | NovaSeq 6000 | +0.79 | +0.16 | −0.95 |
| AMT28 | whoi | NovaSeq 6000 | +0.25 | −0.06 | −0.12 |
| C13.5/A13.5 | whoi | NovaSeq 6000 | +0.47 | +0.16 | +0.54 |
| I07N | whoi | NovaSeq 6000 | +0.26 | −0.11 | +0.71 |
| I09N | whoi | NovaSeq 6000 | +0.62 | +0.14 | +1.46 |
| P18 | whoi | NovaSeq 6000 | +0.64 | −0.73 | +0.17 |
| GOMECC4 | noaa | NovaSeq X (?) | +0.10 | −0.84 | −0.26 |
| I05-2023 | noaa | NovaSeq X (?) | −1.28 | +0.69 | +0.63 |
| I08S-2024 | uci | NovaSeq X | −0.16 | −0.49 | −0.40 |
| I09N-2025 | uci | NovaSeq X | −0.77 | −0.38 | −1.09 |
| P02-2022 | uci | NovaSeq X | −0.29 | −0.91 | −0.64 |
| A13.5-2024 | whoi | NovaSeq X (?) | −0.21 | −0.11 | −0.69 |
| A16N-2023 | whoi | NovaSeq X (?) | +0.22 | +0.57 | −0.74 |

**Group means.**

| Grouping | F1 | F2 | F3 |
|---|---:|---:|---:|
| NovaSeq 6000, uci (1 cruise) | +0.79 | +0.16 | −0.95 |
| NovaSeq 6000, whoi (5) | +0.45 | −0.12 | +0.55 |
| NovaSeq X, noaa (2) | −0.59 | −0.08 | +0.18 |
| NovaSeq X, uci (3) | −0.41 | −0.59 | −0.71 |
| NovaSeq X, whoi (2) | +0.01 | +0.23 | −0.72 |
| NovaSeq 6000 (6) | +0.50 | −0.07 | +0.30 |
| NovaSeq X (7) | −0.34 | −0.21 | −0.46 |

**Eta2 over the 13 cruises** (label permutation, 9,999):

| Axis | Platform | Lab |
|---|---|---|
| F1 | **0.55 (p 0.003)** | 0.34 (p 0.12) |
| F2 | 0.02 (p 0.64) | 0.13 (p 0.50) |
| F3 | 0.27 (p 0.072) | 0.37 (p 0.10) |

**Reading.**
- **The axis-1 offset follows the platform, not the lab.**
  - Within UCI, the NovaSeq 6000 cruise HLY2202 (+0.79) sits far above its three
    NovaSeq X cruises (−0.41). Same lab, library prep at least partly the same group,
    different platform: a 1.20 SD drop.
  - Within WHOI the step has the same sign: NovaSeq 6000 +0.45 against NovaSeq X +0.01,
    a 0.44 SD drop.
  - Together with R39, three findings agree:
    - HiSeq 4000 → NovaSeq 6000 within cruises moves libraries the same way
    - same-platform re-sequencing doesn't move them
    - library prep and preservation don't move them

    So **axis 1's cruise offsets are largely a sequencing-platform batch effect.**
- **The bacterial / eukaryotic composition (axis 1) still has real within-cruise
  structure.** W1 is identical in both workflows (R36) and links weakly to
  stratification and bloom phase (R35, R37). Only the cruise-level offset is the
  platform's.
- **Axis 2 doesn't follow platform** (eta2 0.02), in line with R36.

**Consequences for the model.**
- Carry the platform (NovaSeq 6000 vs NovaSeq X) as a covariate. Or correct it before
  training, for example by matching the axis-1 offsets per platform.
- When the external sets are mapped (R32), record each one's platform. A held-out set
  on a different platform from most of training will carry this offset.

### 3. GC check: job for Morrill

**Question.** If the axis-1 loadings follow the GC content of the KOs' genes, the
mechanism is a chemistry-dependent GC bias:
- HiSeq 4000 four-colour SBS
- NovaSeq 6000 two-colour
- NovaSeq X XLEAP-SBS

**Script** `transfer/scripts/ko_gc_content.py`, run by `transfer/cluster/ko_gc.sbatch`
(morrill partition, 96 GB):
- It streams the eggNOG table (gene → KOs), the taxonomy (gene → domain) and the
  reference FASTA.
- It keeps a deterministic 1-in-10 sample of genes (crc32 of the id).
- Per KO it writes the weighted mean and SD of gene GC, mean length, and mean GC per
  domain. Multi-KO genes are split 1/n, as in the counting.
- Tested on a synthetic set: GC and per-domain values are right.

**Run order.** The FASTA (212 GB) wasn't part of the R32 download.

    cd /scratch/morrill/users/$USER/biogoship_granite_kegg && git pull
    df -h /scratch/morrill/users/$USER        # about 215 GB free needed
    jid=$(sbatch --parsable --account=191001-364393 \
           --export=ALL,FILES="biogoship_2026_database_v2.fasta" cluster/download_reference.sbatch)
    sbatch --account=191001-364393 --dependency=afterok:${jid} transfer/cluster/ko_gc.sbatch

**Bring back** `transfer/results/ko_gc/ko_gc.csv` and the job's .out.

**Then, locally:** correlate GC with the F1 / W1 loadings (overall and within domain),
with the platform contrasts of R39, and with the Protocol A effect.

## R41. All flowcells identified; platforms in training and test sets; summary of the residual-axis work (R33–R41) (2026-11-20)

### 1. The last two flowcells: NovaSeq X

Two more original FASTQ files were downloaded 2026-11-20; the first 100,000 reads of
each were read, then the files were deleted.

| File | First header | Instrument | Run | Flowcell | Lane | Quality bins (first 100 k reads) |
|---|---|---|---:|---|---:|---|
| GOMECC4011_SA_L001_R2_001.fastq.gz | `@LH00316:31:23LT2GLT4:5:1101:34662:1014 2:N:0:NTGTTGGTCC+NACGGTCTAT` | LH00316 | 31 | 23LT2GLT4 | 5 | I 14,356,485; 9 499,540; * 243,598; # 377 |
| GS-A16N-2023-011_SA_L001_R1_001.fastq.gz | `@LH00162:255:22NCJYLT4:8:1101:1671:1042 1:N:0:CGCAGCAATT+ATAGCGGAAT` | LH00162 | 255 | 22NCJYLT4 | 8 | I 12,272,457; 9 424,453; * 125,922; # 1,337 |
| (R40) xR058-L2-P01 READ1 | `@LH00626:65:22YJMLLT4:2:1101:2991:1056 …` | LH00626 | 65 | 22YJMLLT4 | 2 | I / 9 / - / # |

- **All three 2023–2025 flowcells are NovaSeq X** (LH instruments).
- They ran on three different instruments.
- The two WHOI/NOAA runs bin the third quality level as `*` (Q9), while the UCI-provider
  run uses `-` (Q12). That is a different instrument software setting, a possible
  sub-batch within NovaSeq X. Nothing so far shows an effect: the within-X lab means
  differ but not significantly (R40).
- `library_metadata.py` records the header evidence for all three.

**Platform per library now** (KO-table libraries; in brackets, in the real targets):

| Platform | Libraries | In real targets | Source |
|---|---:|---:|---|
| Illumina HiSeq 4000 | 54 | 38 | ENA. None of these samples is in a modelled tile (R36). |
| Illumina NovaSeq 6000 | 901 + 192 | 897 + 186 | ENA (2016–2020); published + flowcell format (HLY2202) |
| Illumina NovaSeq X | 1,899 | 920 | FASTQ headers (all three flowcells) |
| unknown | 35 | 0 | ERR libraries of another project in the KO table |

### 2. Are other platforms missing?

**Training set:** no. Every library in the real targets is on a known platform:
- NovaSeq 6000: 2016–2022
- NovaSeq X: 2023–2025 cruises, including P02-2022 and GOMECC4
- HiSeq 4000: only samples that don't reach a modelled tile

**External test sets** (ENA `instrument_model`, in `data/external/<set>/reads/runs.csv`):

| Set | Platforms (runs) |
|---|---|
| Tara Pacific | NovaSeq 6000 (78), HiSeq 4000 (6) |
| JGI (Flinkstrom) | **HiSeq 2500 (85)**, NovaSeq 6000 (5) |
| EMO BON | NovaSeq 6000 (131) |

- **HiSeq 2500 is in no training sample.** HiSeq 2500 is four-colour SBS like the
  HiSeq 4000, which R39 found shifted toward the bacterial pole of axis 1. So JGI's
  predictions should carry a platform offset of unknown size.
- No test set is on NovaSeq X.
- When the platform covariate is added, the test sets would be predicted at
  "NovaSeq 6000", except JGI. JGI would need either a HiSeq-family value or to be
  reported separately.
- A further check is planned for when external sets are mapped: R39's
  HiSeq 4000 − NovaSeq 6000 contrast, repeated on Tara's 6 HiSeq 4000 runs.

### 3. Summary of findings, R33–R41

**What remains after the satellite model** (Granite 32 PCs + temperature spline,
cruise held out; median KO R2 0.103):
- 84% of KO CLR variance is left, a fifth of it as cruise offsets.
- About a third of the residual variance lies on three reproducible axes. Beyond them
  the residual has no stable shared structure. Ocean-model fields, sampling time and
  depth explain none of the axes.

**Axis 1: bacterial vs eukaryotic genes (17–19% of the residual).**
- **Its cruise offsets are mostly a sequencing-platform batch effect.**
  - They split NovaSeq 6000 vs NovaSeq X (eta2 0.55, p 0.003). Lab doesn't explain
    them (0.34, p 0.12).
  - Within UCI and within WHOI, moving from NovaSeq 6000 to NovaSeq X lowers the cruise
    means: by 1.20 and 0.44 SD.
  - HiSeq 4000 vs NovaSeq 6000 within the same cruises shifts libraries the same way:
    +1.19× the era offset, p < 0.001, in all four cruises.
  - Re-sequencing the same DNA on the same platform shifts nothing (HLY2202).
- Ruled out as drivers of the offset:
  - library prep: Nextera tagmentation vs DNA Flex (HLY2202)
  - preservation and bead beating: the A16N Protocol A/B pairs
  - sampling year, within one extraction workflow
  - flowcell, storage time, volume, collection method, extraction plate
- **Its within-cruise part (W1) is real and the same in every workflow.** It links
  weakly to stratification (DCM strength, nitracline below the mixed layer, low
  surface nitrate; r ≈ 0.2) and to bloom phase (after a peak; r ≈ 0.14).
- The R33 "temperature" link was an artefact of a five-cruise subset; its sign differs
  by era.

**Axis 2: phototrophs (Cyanobacteria / Cryptophyta, photosynthesis genes; 10%).**
- **It is biology.**
  - It has the same loadings in the old and new workflows: cosine 0.95 and 0.94 with
    the full fit, old vs new 0.79 against a random-split median of 0.81.
  - It isn't related to platform (eta2 0.02).
  - It survives dropping the Protocol A samples.
- **It isn't ridge shrinkage.**
  - Calibration slopes are below 1: the predictions are over-dispersed, not compressed.
  - The axis correlates negatively with the predicted signal.
  - It survives recalibration (cosine 0.985).
- **It is spatially smooth.** Decorrelation length is about 560 km. About half its
  variance varies over hundreds of km, about a quarter is replicate-level noise, which
  is the predictability floor.
- **No tested variable explains it:**
  - surface and model fields
  - nutricline and DCM, from CCHDO profiles
  - bloom phase, from 28-day satellite chlorophyll
  - protocols and sampling time

  It is the part of the phototroph community orthogonal to what Granite and
  temperature capture. The model does predict the Cyanobacteria share itself well
  within cruise (r 0.79).
- **Untested candidates:** light history (cloud / PAR, mixed-layer light), water-mass
  origin (back-trajectories), grazing and viruses.

**Axis 3: sequencing depth (6–7%).**
- Technical: it comes from the CLR pseudocount 0.5 / total.
- Log total counts as a covariate removes it completely and raises the median KO R2 to
  0.145.
- Protocol A samples also load on it.

**Protocol decisions made.**
- A16N Protocol A samples are out of the real targets: no beads + lysis buffer, a
  strong technical effect.
- C13.5 v1/v2, GOMECC4 A/B/C and HLY re-sequencing are replicates and stay averaged.
- The P02 Shield pairs can't be tested: their 2023 libraries aren't in the KO table.

**Model changes agreed, not yet tested:**
- platform covariate: NovaSeq 6000 / NovaSeq X / HiSeq family
- sequencing-depth covariate
- alpha chosen by leave-cruise-out (R38 found the predictions over-dispersed)

These wait until the remaining variables are worked out.

**Pending:**
- the GC job on Morrill (R40)
- light history and back-trajectories for axis 2
- whether the P02 HWKVFDSX5 libraries can be counted

## R42. Axis 2: viral read share (null); water-mass history and light-history job (2026-11-20)

**Order of tests** (suggested 2026-11-20, now that axis 2 isn't shrinkage):
1. water-mass origin and age, along Parcels back-trajectories
2. light history
3. grazing and viral lysis

The viral share is available locally, so it ran first. The other two need a Morrill
job, which is ready.

### 1. Viral read share: no link to axis 2

**Script** `transfer/scripts/viral_share_test.py`; outputs in `transfer/results/viral_share/`.
- Data: domain and phylum read counts per library (`data/targets/taxa/`), summed per
  sample over the real-target libraries.
- Per tile: the mean over its samples.
- Real model: 1,215 tiles, R36 scores.

To reproduce, run from the repository root:

    OMP_NUM_THREADS=4 ~/miniforge3/bin/python transfer/scripts/viral_share_test.py

**Spearman r**, all tiles / within cruise (cruises with the pooled sign, of 13):

| Feature | F1 | F2 | W1 | W2 | W3 |
|---|---|---|---|---|---|
| Viruses CLR (4 domains) | +0.15 / +0.18 (8) | +0.15 / +0.04 (8) | +0.17 / +0.19 (9) | +0.02 / +0.02 (9) | −0.17 / −0.19 (10) |
| log10 Viruses / Bacteria | −0.22 / −0.12 (8) | +0.17 / −0.01 (6) | −0.10 / −0.12 (9) | +0.01 / +0.01 (8) | +0.01 / −0.00 (7) |
| log10 Viruses / Cyanobacteriota | −0.01 / +0.01 (6) | −0.23 / −0.38 (11) | −0.03 / −0.03 (7) | −0.34 / −0.38 (11) | +0.08 / +0.09 (11) |
| viral share of classified reads | −0.13 / −0.04 (7) | +0.17 / +0.02 (9) | −0.03 / −0.03 (7) | +0.03 / +0.04 (10) | −0.05 / −0.07 (9) |

**Does the viral share itself carry the platform batch?** Eta2 of cruise means,
NovaSeq 6000 vs X:

| Feature | Eta2 | p |
|---|---:|---:|
| Viruses CLR | 0.07 | 0.40 |
| Viruses / Bacteria | 0.11 | 0.27 |
| Viruses / Cyanobacteriota | 0.05 | 0.66 |
| viral share | 0.13 | 0.23 |

None does.

**Reading.**
- **No viral signal on axis 2.** The viral share and Viruses / Bacteria have
  within-cruise r ≤ 0.04 with W2.
- **The Viruses / Cyanobacteria ratio is circular.** It correlates −0.38 with W2, but
  only through its denominator: W2 tracks Cyanobacteriota at +0.38 (R38).
- The viral CLR goes weakly with W1 (+0.19) and against W3 (−0.19).
- Caveats:
  - On 0.22 µm Sterivex filters most viral reads are prophages, large viruses inside
    cells, or virions caught on the filter. So this is a weak lysis proxy.
  - It is descriptive, like taxa.
- Grazing has no measurable proxy in these data.

### 2. Water-mass history and light: job for Morrill

**Script** `transfer/scripts/water_mass_history.py`, run by
`transfer/cluster/water_mass.sbatch` (service partition, Copernicus Marine login,
biogoship conda env with parcels).

**Steps.**
1. **Currents.** GlobCurrent at 0.25°, daily, 15 m, per section:
   - the samples' box ±15°, from 62 days before the first sample to 1 day after the
     last
   - multi-year, near-real-time after 2026-03-31
   - boxes crossing the date line (P02, P04E, HLY2202 at the Bering Strait) are
     fetched in two pieces and stored in 0–360°
   - stored in `transfer/data/currents60/`
2. **Tracking.** One particle per sample (2,014), released at the sampling position
   and time and advected **60 days backward**:
   - OceanParcels RK4, dt 30 min, no diffusion; positions every 6 h
   - the same kernel as `biogoship/advection.py`, backward only
   - output `trajectories.csv`, with a `stuck_at_end` flag for particles that reached
     land
3. **Fields** at the daily positions, days 0 to −60, nearest grid value. Output
   `along_track.csv`, in long format.
   - OSTIA SST
   - GlobColour gap-free chlorophyll
   - GLORYS mixed layer
   - GlobColour gap-free KD490, if the transparency dataset opens
   - GLORYS sea-floor depth
4. **Light probe.** It lists which Copernicus ocean-colour datasets open, with their
   variables. That finds a daily PAR product for the light test. If Copernicus has
   none, NASA OB.DAAC daily PAR (MODIS-Aqua / VIIRS) through the Earthdata login is the
   fallback (`build_new_cruise_manifest.py` already has the Earthdata download code).

**Checks run locally.**
- The section boxes from the local manifests: 14 sections, 2,014 samples.
- The date-line fix for HLY2202: its box is now 176.6–225.0°.
- **Not run locally:** the Parcels step (parcels isn't on the Mac). It uses the same
  Parcels calls as `cluster/backtrack.sbatch`, which has run on Morrill.

**Run on Morrill:**

    cd /scratch/morrill/users/$USER/biogoship_granite_kegg && git pull
    sbatch --account=191001-364393 transfer/cluster/water_mass.sbatch

**Bring back:** `transfer/results/water_mass/trajectories.csv`, `along_track.csv` and
the job's .out.

**Features to compute on the Mac** (`water_mass_test.py`, to be written when the data
are back), against W2 and the other axes:
- **Water-mass age:**
  - days since the track was last over the shelf (sea floor < 200 m; censored at 60)
  - days since the water was ≥ 1 °C colder than at sampling (time since upwelling or
    mixing)
- **Trajectory means:** SST and log chlorophyll over 30 and 60 days; the SST and
  chlorophyll change along the track; the 60-day maximum mixed layer.
- **Origin:** position at −30 and −60 days, displacement and path length, and origin
  region (basin / latitude band, or Longhurst province if a table is added).
- **Light:**
  - mean PAR over 14 and 28 days along the track
  - mixed-layer light = PAR × (1 − e^(−Kd·MLD)) / (Kd·MLD), with Kd(PAR) from KD490
    (Morel et al. 2007)

### R42 addendum: how much variance is the sequencing platform? (2026-11-20)

**Question.** Does axis 1's 17–19% mean that this much of the KO variance is platform?

**No.**
- The 17–19% is axis 1's share of the *standardised residual* variance, with every KO
  weighted equally.
- Only axis 1's *cruise-level* part follows platform. Its within-cruise part (W1) is
  the same in every workflow (R36).

**Direct estimate.**
- Real model: 1,215 tiles, 3,933 KOs, baseline 32 PCs + temperature spline, cruise
  held out.
- Statistic: the share of variance explained by the two platform group means:
  - NovaSeq 6000: I09N, P18, I07N, AMT28, C13.5/A13.5, HLY2202
  - NovaSeq X: the other 7 cruises
- Computed inline, not saved as a script.

| Quantity | Value |
|---|---:|
| residual share of KO CLR variance | 83.8% |
| platform, share of the observed KO CLR variance | 4.5% |
| platform, share of the residual variance | 4.5%, = **3.8% of the total** |
| platform, share of the standardised residual (KOs weighted equally) | 5.8% |
| all 13 cruise offsets, share of the residual | 20.1% (16.9% of the total) |
| per KO, platform share of the residual: median / 90th percentile / max | 2.8% / 15.6% / 44.1% |

**Reading.**
- Platform accounts for about **4% of KO variance overall** (3.8% of the total, 4.5% of
  the residual). That is roughly a fifth of all cruise-to-cruise offsets.
- It is concentrated: for a tenth of the KOs it is more than 15% of their residual,
  and for one KO 44%.
- It is an upper bound for the platform itself. The two platform groups are also two
  sets of years and regions, so part of it may be real between-cruise difference.
- It is part of the unexplained variance, not of the model's predictions. Under cruise
  hold-out the model can't learn a platform offset, so the offset stays in the
  residual. A platform covariate is the way to take it out.

## R43. Variance budget (first draft) and residual scree plot (2026-11-20)

**Script** `transfer/scripts/variance_budget.py`.
- Table: `transfer/results/variance_budget/budget.csv`.
- Figures: `transfer/figures/variance_budget.png` and `transfer/figures/residual_scree.png`.
  Arial, no titles.

**Data.** Real model: Protocol A excluded, 1,215 tiles × 3,933 KOs, baseline 32 PCs +
temperature spline, cruise held out.

**Partition.** Each step works on what the previous ones left, summing squares over
tiles and KOs:
1. Explained by the model: total − the cross-validated residual.
2. **Platform**: the two group means of the residual, NovaSeq 6000 cruises vs NovaSeq X
   cruises.
3. **Other cruise offsets**: cruise means minus the platform means.
4. **Sequencing depth**: the within-cruise residual regressed on a spline of log10 total
   KO counts. 8 knots, cubic, basis centred per cruise, in-sample OLS per KO.
5. **W1**, then 6. **W2**: regression on the R36 within-cruise axis scores.
7. **Replicate noise**:
   - σ² per KO = mean (a − b)² / 2 over 305 same-tile sample pairs from 133 tiles:
     GOMECC4 A/B/C, C13.5 v1/v2 and co-located samples.
   - Scaled by Σ 1/k over the tiles, where k is the number of samples averaged in the
     tile.
   - Capped at what is left.
8. The rest.

**Two scalings.** Raw CLR units, where high-variance KOs weigh more, and every KO
standardised.

| Component | Raw | KOs standardised | Status |
|---|---:|---:|---|
| Explained by the satellite model (CV) | 16.2% | 11.4% | model |
| Sequencing platform, NovaSeq 6000 vs X (upper bound) | 3.8% | 5.1% | technical, known |
| Other cruise offsets | 13.1% | 14.5% | **unknown**: region, year, other batch |
| Sequencing depth | 2.3% | 3.5% | technical, known |
| Axis 1 within cruise (bacteria vs eukaryotes) | 7.4% | 10.6% | biology; weak links to stratification and bloom phase |
| Axis 2 within cruise (phototrophs) | 3.6% | 5.9% | biology; **driver unknown** |
| Replicate noise (unpredictable floor) | 27.0% | 26.0% | noise floor |
| Remaining within cruise | 26.6% | 23.1% | **unknown** |

**Scree** (`residual_scree.png`): share of standardised residual variance per component,
for the full (F) and within-cruise (W) residual PCA of R36, against the permutation-null
95th percentile.
- Filled markers have split-half congruence ≥ 0.75: F1, F2, W1, W2. W3 (0.65) is just
  under.
- The curve flattens after component 3 to 5. Every component is far above the null
  (about 0.2%), which is why the null was too weak a criterion (R33).

**Reading.**
- **A large unpredictable floor.**
  - Two samples of the same water differ by enough to account for about 27% of all KO
    variance. That covers separate filters, extractions and libraries, plus patchiness
    below a tile.
  - So no predictor can push the cross-validated KO R2 much past about 0.7 overall.
  - The model currently explains 16%. Known technical factors (platform and depth) take
    about 6% more.
- **Known biology without a driver:** about 11%. Axis 1 within cruise is 7.4%, axis 2
  3.6%.
- **Unknown:** about 40%.
  - 13% is cruise offsets beyond platform: region, year, or a batch not yet found. Under
    cruise hold-out this part can't be learned at all.
  - 27% is within-cruise variation that is neither noise nor a reproducible shared
    axis: KO-specific or small-scale structure.
- **Caveats.**
  - The platform share is an upper bound, because platform is confounded with years and
    regions.
  - The noise floor may be slightly overestimated. Replicate libraries also differ in
    depth (overlap with step 4), and C13.5 v2 has a small Bacteria shift (R36).
  - W1 and W2 come from a PCA of standardised residuals, so their raw-unit shares are
    smaller than the R33 percentages.
- **What changes this draft.**
  - Water-mass and light history may move part of "remaining" and of axis 2 into a
    known driver.
  - The GC check may tighten the platform attribution.
  - A platform covariate and a depth covariate in the model would move their slices
    into "explained".

## R44. GC content per KO: the platform offset is not a sequencing GC bias; water-mass fields step fixed (2026-11-20)

**Question.** Is the axis-1 platform offset (R39/R40) caused by GC-dependent sequencing chemistry? This would
mean HiSeq 4000 (four-colour) vs NovaSeq 6000 / NovaSeq X (two-colour) cover low- and high-GC genes
differently. If so, the per-KO platform effect should follow the GC content of the KO's genes, and GC
should explain why that effect resembles the axis-1 loading.

**Data.**
- Per-KO GC came from `transfer/scripts/ko_gc_content.py` on Morrill (`transfer/cluster/ko_gc.sbatch`, after
  downloading `biogoship_2026_database_v2.fasta`, 212 GB, job 88833).
  - It uses a deterministic 1-in-10 sample (zlib.crc32 of the sequence id) of the reference genes with a
    KEGG KO.
  - GC of each gene is GC / (A+C+G+T). Each KO gets the gene-weighted mean, with weight 1/n for genes
    carrying n KOs.
  - The output also has means overall and per domain, plus mean gene length and the eukaryote share
    (genes_Eukaryota / genes).
  - Output: `transfer/results/ko_gc/ko_gc.csv`, 19,210 KOs.
- All 3,933 real-model KOs have a GC value. The median is 1,044 sampled genes per KO; GC median 0.517,
  IQR 0.490–0.540.

**Method.** Script: `transfer/scripts/ko_gc_test.py`. Per-KO effects are in units of the KO's pooled
within-cruise SD, as in R39. Library CLR is rarefied to 250,000 KO reads (seed 0) and cached in
`rarefied_clr.npz`. The effects:
- `hiseq`: HiSeq 4000 − NovaSeq 6000, different samples within a cruise, surface libraries.
  - Cruises: I09N-2016, AE1319-2013, BVAL46-2011, NH1418-2014.
  - Per-cruise mean difference, pooled with each cruise weighted by its smaller group.
- `era`: mean of the cruise-mean library CLR, NovaSeq 6000 cruises minus NovaSeq X ("…LT4") cruises.
  - Modelled surface libraries only.
  - This contrast mixes platform, years and regions.
- `rerun`: HLY2202, the same DNA on two NovaSeq 6000 flowcells (HMY7WDSX7 − HFJLFDSX7), 21 pairs.
  - This is the control: same chemistry.
- `protocolA`: A16N Protocol A − B, 16 pairs.
- `c13_v2`: C13.5 v2 − v1, 24 pairs (duplicates).
- The loadings F1–F3 and W1–W3 (`residual_modules_noA`).

**Statistics.**
- Spearman ρ of each effect with GC.
- Partial Spearman controlling for the eukaryote share and mean gene length (ranks residualised).
- ρ with bacterial-gene GC, for KOs with ≥ 20 sampled bacterial genes (3,507 KOs).
- Share of each effect's per-KO variance explained by a cubic spline in GC (6 knots).
- Whether controlling for GC removes the resemblance of the technical contrasts to F1/W1.
- KOs are not independent, so p values are not reported. The comparison against the rerun control is
  what carries the argument.

**Results** (`transfer/results/ko_gc/gc_correlations.csv`, `gc_deciles.csv`, `gc_effects.csv`, `summary.txt`)

| Effect | ρ with GC | Partial ρ (euk share, length) | ρ with bacterial-gene GC | ρ with eukaryote share |
|---|---:|---:|---:|---:|
| F1 | −0.268 | −0.211 | −0.198 | −0.41 |
| W1 | −0.299 | −0.240 | −0.242 | −0.41 |
| F2 | −0.213 | −0.235 | −0.259 | +0.16 |
| W2 | −0.188 | −0.224 | −0.229 | +0.19 |
| F3 | −0.122 | −0.114 | −0.093 | +0.11 |
| W3 | −0.090 | −0.131 | −0.195 | +0.32 |
| HiSeq − NovaSeq 6000 (within cruise) | −0.257 | −0.250 | −0.304 | −0.02 |
| NovaSeq 6000 − X cruises (era) | −0.207 | −0.177 | −0.133 | −0.22 |
| Same DNA, two NovaSeq 6000 runs (control) | +0.187 | +0.196 | +0.200 | −0.03 |
| A16N Protocol A − B | −0.077 | −0.050 | −0.015 | −0.01 |
| C13.5 v2 − v1 | +0.146 | +0.143 | +0.160 | +0.04 |

GC vs eukaryote share: ρ +0.223. Gene length vs the effects: |ρ| ≤ 0.10 for all of them.

| Effect | ρ with F1 | Partial ρ with F1 (GC, euk share, length) | ρ with W1 | Partial ρ with W1 | Variance explained by a GC spline | RMS (SD) | GC-decile range (SD) |
|---|---:|---:|---:|---:|---:|---:|---:|
| HiSeq − NovaSeq 6000 | +0.617 | +0.665 | +0.665 | +0.713 | 0.083 | 0.445 | 0.382 |
| Era (6000 − X) | +0.587 | +0.538 | +0.497 | +0.432 | 0.046 | 0.402 | 0.260 |
| Rerun (control) | −0.037 | −0.018 | −0.034 | −0.007 | 0.039 | 0.236 | 0.144 |

Mean effect by GC decile (SD units; loadings as they are):

| Decile | GC | F1 | W1 | F2 | W2 | F3 | W3 | HiSeq | Era | Rerun | Protocol A | C13.5 v2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.439 | 0.231 | 0.239 | 0.086 | 0.082 | 0.007 | 0.022 | 0.172 | 0.117 | −0.091 | 0.013 | −0.044 |
| 1 | 0.475 | 0.337 | 0.353 | 0.141 | 0.097 | 0.061 | 0.076 | 0.310 | 0.170 | −0.060 | 0.087 | −0.086 |
| 2 | 0.490 | 0.317 | 0.321 | 0.122 | 0.091 | 0.035 | 0.037 | 0.273 | 0.194 | −0.030 | 0.026 | −0.061 |
| 3 | 0.501 | 0.301 | 0.304 | 0.108 | 0.083 | 0.027 | 0.013 | 0.284 | 0.205 | −0.007 | −0.007 | −0.032 |
| 4 | 0.512 | 0.173 | 0.174 | 0.071 | 0.054 | 0.011 | −0.011 | 0.156 | 0.128 | 0.010 | −0.028 | −0.006 |
| 5 | 0.521 | 0.106 | 0.106 | 0.001 | −0.013 | 0.007 | −0.009 | 0.094 | 0.114 | 0.015 | −0.029 | 0.008 |
| 6 | 0.530 | 0.023 | 0.028 | 0.015 | 0.006 | 0.002 | −0.001 | 0.030 | 0.042 | 0.031 | −0.005 | −0.008 |
| 7 | 0.540 | 0.024 | 0.025 | −0.022 | −0.040 | 0.002 | 0.027 | 0.061 | 0.028 | 0.035 | 0.000 | 0.027 |
| 8 | 0.552 | −0.005 | −0.002 | −0.032 | −0.040 | −0.031 | −0.006 | −0.014 | −0.022 | 0.043 | −0.005 | 0.009 |
| 9 | 0.574 | 0.037 | 0.027 | −0.096 | −0.095 | −0.099 | −0.053 | −0.072 | −0.055 | 0.053 | −0.057 | 0.006 |

Figure: `transfer/figures/ko_gc_effects.png`. Each panel plots per-KO effect against gene GC, with the mean
in 20 GC bins. Panels, left to right: HiSeq − NovaSeq 6000; the same-DNA rerun; NovaSeq 6000 − X; the F1
loading.

**Conclusions**
1. **GC is not the mechanism of the platform offset.**
   - GC explains 8% of the per-KO HiSeq − NovaSeq effect and 5% of the era effect.
   - The resemblance of the HiSeq effect to axis 1 does not shrink when GC is controlled for:
     F1 +0.62 → +0.67, W1 +0.67 → +0.71.
   - So the platform effect lands on axis 1 through something other than GC.
2. **GC-linked run-to-run variation exists, but it is ordinary.**
   - The same DNA on two NovaSeq 6000 runs shows a GC slope of similar size and the opposite sign:
     ρ +0.19, decile range 0.14 SD, against 0.38 SD for HiSeq.
   - The C13.5 duplicates show the same kind of slope (ρ +0.15).
   - These GC slopes have no axis-1 shape (rerun vs F1 ρ −0.04).
   - So a GC tilt is a generic batch nuisance here, not what makes axis 1.
3. **Axis 1 itself leans low-GC** (F1 vs GC ρ −0.27; low- and mid-GC deciles load +0.23 to +0.34).
   - This survives controlling for the eukaryote share (−0.21) and holds among bacterial genes (−0.20).
   - The likely reason is taxonomic: KOs typical of low-GC oligotrophs (e.g. SAR11, Prochlorococcus) load
     on the bacterial side of axis 1.
   - This is biology in the axis, not evidence of a chemistry bias.
4. **The platform offset is still unexplained mechanistically.** Remaining candidates, none tested yet:
   - library-prep kit or batch that changed together with the platform;
   - read length or insert size, and therefore mapping and KO-assignment sensitivity;
   - duplicate rate or index hopping;
   - quality-filter behaviour on two-colour data (poly-G tails).
   For modelling this changes nothing: platform enters as a covariate (agreed; testing on hold).
5. **The Protocol A effect is not GC-linked** (ρ −0.08). Its F3 resemblance (ρ +0.57) is its own pattern.

**Water-mass job (R42): fields step fixed.**
- The trajectories came back complete: 2,014 samples, 60 days back at 6-h steps, none stuck or lost
  (`transfer/results/water_mass/trajectories.csv`).
- The along-track fields, however, held day 0 only for nearly every track:
  - Parcels outputs on a 6-h clock offset by a fraction of an hour from each sample's time, and the
    script kept only exact whole-day points.
  - `daily_positions` now interpolates each track linearly to days 0 … −60. Longitude is unwrapped
    before interpolating, then wrapped back.
  - Checked: 61 points for all 2,014 tracks; day 0 equals the release position; the 18 tracks crossing
    180° wrap correctly.
- The output is now `along_track_daily.csv`; `along_track.csv` from the first run is obsolete.
- The sbatch takes `STEPS=fields` to re-run that step alone. Commit 38d2e93.

## R45. Read length, mapping yield and KO-profile spread: axis 1 is largely a KO-profile spread axis, and the platform offset runs through it (2026-11-20)

From here on, new files go in the normal repo folders (`scripts/`, `data/results/`, `figures/output/`), not
`transfer/`.

**Question.** Is the platform offset on axis 1 a mapping or annotation effect rather than sequencing
chemistry? The hypothesis has three parts:
- newer runs have longer reads;
- longer reads map more eukaryotic and divergent genes;
- so the eukaryotic share rises on newer runs.

**Pipeline facts** (biogoship-metagenomics/mapping_pipeline, read-only):
- Mapping is `minimap2 -ax sr` against the reference, with no identity or MAPQ filter, so per-read hit
  identity is never logged.
- Reads per reference gene come from `samtools coverage` (numreads: both mates; secondary, QC-fail and
  duplicate reads are excluded, but duplicates are never marked).
- Genes are kept when their breadth of coverage passes `--coverage_threshold`; this is 50% for
  `ko_counts_cov50_all.tsv`.
- KO counts come from the eggNOG table.

**Method.** Script: `scripts/platform_mapping_yield.py`.
- Libraries are taken from `transfer/results/library_metadata/library_metadata.csv`: those in the KO table,
  with blanks and controls dropped and > 100k reads.
- Per-library quantities:
  - nominal read length (FAIRe seq_run_id);
  - ENA base_count / read_count;
  - mapping rate: domain-assigned reads (sum of `data/targets/taxa/domain_counts.tsv`) / total read pairs
    (FAIRe input_read_count). Both mates count, so this is reads per pair, range 0–2;
  - KO yield: KO-assigned reads / total pairs;
  - KO share of mapped reads;
  - eukaryote fraction: Eukaryota / (Archaea + Bacteria + Eukaryota);
  - log10 reads;
  - from the library KO CLR rarefied to 250,000 KO reads (seed 0; as R39/R44):
    - KO-profile spread, the SD of the CLR over the 3,933 KOs;
    - undetected KOs, the share of KOs with zero rarefied count;
    - Shannon diversity of the KO proportions;
    - the F1 projection.
- Tests:
  1. Platform and cruise medians.
  2. HiSeq 4000 − NovaSeq 6000 within the four mixed cruises (surface libraries): pooled, each cruise
     weighted by its smaller group; permutation within cruise (1,999).
  3. HLY2202, the same DNA on two NovaSeq 6000 flowcells (40 sample pairs over all depths), as a control.
  4. Within-cruise Spearman correlations (cruise-centred; 2,863 libraries).
  5. Mediation: the HiSeq − NovaSeq F1 difference after removing a within-cruise linear fit of F1 on each
     candidate.
  6. Per-KO effects (R44 `gc_effects.csv`) vs eukaryote share of genes, mean abundance (mean CLR over
     targets) and gene length.
  7. Cruise-level Spearman correlations.
- Outputs:
  - `data/results/platform_mapping/library_yield.csv` (per library), `cruise_yield.csv`, `tests.csv`,
    `summary.txt`;
  - `figures/output/platform_mapping_yield.png`, three panels:
    - mapped reads per pair by platform in the mixed cruises;
    - the libraries' axis-1 scores there;
    - within-cruise axis-1 score vs KO-profile spread, all libraries.

**1. Read length does not step with platform.**
- Nominal length is PE150 on HiSeq 4000, NovaSeq 6000 and NovaSeq X.
- ENA bases per spot: median 302 for HiSeq 4000 (n 52) and for NovaSeq 6000 (n 900). NovaSeq X runs are
  not in ENA.

**Medians by platform:**

| Platform | Libraries | Reads (pairs) | Mapped / pair | KO yield | KO / mapped | Eukaryote fraction |
|---|---:|---:|---:|---:|---:|---:|
| HiSeq 4000 | 52 | 12.4 M | 0.318 | 0.213 | 0.686 | 0.017 |
| NovaSeq 6000 | 1,082 | 11.4 M | 0.425 | 0.273 | 0.655 | 0.024 |
| NovaSeq X | 1,873 | 18.2 M | 0.419 | 0.258 | 0.607 | 0.030 |

The cruise medians are in `cruise_yield.csv`. Across single-platform cruises (n 12), cruise mean F1
correlates (Spearman) with:
- KO / mapped: +0.94
- eukaryote fraction: −0.60
- read depth: −0.42
- mapping rate: +0.19
- KO yield: +0.25

**2. HiSeq 4000 − NovaSeq 6000 within the mixed cruises** (pooled; permutation p):

| Metric | Difference | p |
|---|---:|---:|
| Mapped / pair | −0.177 | 0.0005 |
| KO yield | −0.115 | 0.0005 |
| KO / mapped | +0.013 | 0.10 |
| Eukaryote fraction | −0.005 | 0.26 |
| log10 reads | +0.110 | 0.005 |
| KO-profile spread (SD of CLR) | +0.124 | 0.0005 |
| Undetected KOs (share) | +0.041 | 0.0005 |
| Shannon | −0.059 | 0.0005 |
| F1 | +17.44 | 0.0005 |

Within-cruise SDs, for scale: spread 0.129, undetected 0.037, Shannon 0.059. So the HiSeq libraries:
- map about 18 percentage points fewer reads per pair;
- have more reads;
- have the same KO share among mapped reads and the same eukaryote fraction;
- have a less even KO profile: about +1 within-cruise SD of spread and 4% more KOs undetected.

Per cruise (HiSeq vs NovaSeq medians of mapped / pair): I09N-2016 0.312 vs 0.364; AE1319-2013 0.451 vs
0.702; BVAL46-2011 0.366 vs 0.649; NH1418-2014 0.169 vs 0.697.

**3. The same-DNA rerun (HLY2202, n 40) shows no shift in any metric.**
- Mean differences: mapped / pair +0.012, KO / mapped +0.004, eukaryote fraction +0.023, spread +0.010,
  undetected +0.007, Shannon +0.003, F1 −0.49.
- The pair SDs are large (mapped / pair 0.21, F1 18.0), but there is no shift.

**4. Within cruise (Spearman, n 2,863):**
- **F1 vs:** spread +0.729, Shannon −0.670, KO / mapped +0.637, eukaryote fraction −0.588, undetected
  +0.540, KO yield +0.263, log10 reads −0.177, mapped / pair +0.153.
- **Eukaryote fraction vs:** spread −0.561, undetected −0.428, KO / mapped −0.410, KO yield −0.350,
  mapped / pair −0.298, Shannon +0.333, log10 reads +0.044.
- **Spread vs:** undetected +0.923, Shannon −0.493, log10 reads −0.378, KO / mapped +0.336.
- **Undetected vs log10 reads:** −0.449.

**5. Mediation of the HiSeq − NovaSeq F1 difference (+17.44):**

| Removed (within-cruise linear fit) | Within-cruise R² of F1 | F1 difference after | Share explained |
|---|---:|---:|---:|
| Mapped / pair | 0.020 | +21.67 | −24% |
| KO / mapped | 0.391 | +14.12 | +19% |
| log10 reads | 0.034 | +19.48 | −12% |
| Mapped / pair + KO / mapped + log10 reads | 0.455 | +18.73 | −7% |
| KO-profile spread | 0.539 | +1.07 | +94% |
| Shannon | 0.468 | +1.40 | +92% |
| Undetected KOs | 0.247 | +4.44 | +75% |
| Eukaryote fraction | 0.419 | +16.93 | +3% |
| Spread + eukaryote fraction | 0.669 | +4.71 | +73% |

**6. Per KO (Spearman):**

| Effect | vs eukaryote share of genes | vs mean abundance (mean CLR) | vs gene length |
|---|---:|---:|---:|
| HiSeq − NovaSeq 6000 | −0.017 | +0.570 | −0.096 |
| Era (6000 − X) | −0.223 | +0.322 | +0.091 |
| Rerun (control) | −0.032 | −0.016 | +0.046 |
| F1 loading | −0.410 | +0.709 | +0.102 |

**7. Cruise level (all 16 cruises), Spearman with mean F1:** spread +0.83, Shannon −0.88, undetected
+0.65, eukaryote fraction −0.66.

**Conclusions**
1. **The read-length and annotation-yield mechanism is not supported.**
   - Read length is PE150 everywhere.
   - HiSeq libraries do map fewer reads, but the KO share among mapped reads and the eukaryote fraction
     are unchanged.
   - Mapping rate explains none of the F1 offset (removing it slightly enlarges it).
   - The eukaryote fraction explains 3% of it.
   - The per-KO HiSeq effect is unrelated to how eukaryotic a KO is (ρ −0.02).
2. **Axis 1 is largely a KO-profile spread (evenness) axis.**
   - The F1 loading tracks KO mean abundance (ρ +0.71): a high F1 score means common KOs are more
     dominant and rare KOs fall out.
   - A library's spread explains 54% of its within-cruise F1 variance.
   - At cruise level, spread vs F1 is ρ +0.83 and Shannon vs F1 is −0.88.
   - Part of this is biology. Eukaryote-rich libraries have more even KO profiles (eukaryote fraction vs
     spread ρ −0.56), and the eukaryote fraction alone explains 42% of within-cruise F1. That is the R33
     "bacteria vs eukaryotes" reading.
3. **The platform offset runs through spread, not through composition.**
   - HiSeq libraries are less even (+1 SD; 4% more KOs undetected) with the same eukaryote fraction and
     more reads.
   - Removing spread removes 94% of the HiSeq − NovaSeq F1 difference.
   - The same DNA rerun on the same chemistry shows no spread shift.
   - So the platform changes how many distinct KOs are detected, and axis 1 picks that up.
4. **Likely mechanism (untested): library complexity combined with the 50% breadth filter.**
   - Fewer unique molecules (PCR or optical duplicates, which this pipeline does not remove; or low-input
     library prep) give fewer distinct fragments per gene.
   - Then fewer rare genes reach 50% breadth and their KOs drop to zero, even at equal read counts.
   - This is consistent with: more reads not helping the HiSeq libraries; undetected KOs falling with depth
     within cruise (ρ −0.45); and the effect landing on rare KOs.
   - It is a counting-pipeline effect, not chemistry or GC (R44).
   - Direct test, to do: duplicate rate / library complexity per library from the FASTQs on Morrill (e.g.
     `fastp` duplication or `samtools markdup` on a subsample of BAMs) for the four mixed cruises and the
     LT4 vs NovaSeq 6000 era. Also: re-count with `--coverage_threshold 0` for a few libraries and check
     whether the spread offset disappears.
5. **For modelling** (testing on hold, as agreed):
   - A per-library spread / detected-KO covariate (or a count transform that does not reward detection)
     would likely absorb the platform offset better than a platform label, but it carries some biology
     (the eukaryote fraction).
   - Dropping rare KOs from the targets is another option.
   - None of these has been tested yet.

## R46. Does the CLR amplify KO dropouts? Single-copy-marker denominator ("per genome equivalent"), coverage-50% counts (2026-11-20)

**Question.**
- In the CLR, the denominator is the geometric mean over all 3,933 KOs. When rare KOs drop to the
  pseudocount, that mean falls and every detected KO's CLR rises.
- So a library that detects fewer KOs looks like "common KOs up": the axis-1 / platform pattern of R45.
- A fixed reference of universal single-copy genes would give log copies per genome equivalent (GE):
  GE_k = log(p_k + e) − mean over markers m of log(p_m + e).
- Equivalently, from the same rarefied proportions, GE_k = CLR_k − mean_M CLR_m.
- With this denominator, rare-KO dropouts cannot move common KOs.

**Markers.**
- The set: prokaryote-type ribosomal proteins in the KO set, i.e. KEGG names with an rpl* / rps* / rpm* gene
  symbol and "ribosomal protein" (rim* / prm* modifying enzymes excluded).
- All 54 are detected in every rarefied library, and all are used (`markers.csv`).
- They are universal single-copy in Bacteria and Archaea. The mitochondrial and plastid MRPL/MRPS copies
  share these KOs.
- Eukaryotic cytosolic ribosomal proteins are separate "…e" KOs. Their 27 KOs have prevalence 0.71–0.94, so
  they are not universal here; reported, not used.
- GE is therefore per prokaryote genome equivalent: eukaryotic KOs rise with the eukaryote : prokaryote
  ratio instead of being normalised by it.
- Within-library SD of the 54 markers' log proportions (median): HiSeq 4000 0.513, NovaSeq 6000 0.504,
  NovaSeq X 0.471. So the reference mean has an SE of about 0.07.

**KOs that are not 100% prevalent.**
- A zero stays at the rarefied detection limit, log(0.5 / N) minus the marker mean.
- Every library is rarefied to the same N, so the limit is comparable across libraries.
- Unlike in the CLR, a zero no longer shifts any other KO.
- A KO's own zeros still need care. Options: treat them as censored, or restrict a given analysis to KOs
  above a prevalence cut.

**Method.** Script: `scripts/marker_normalisation.py --counts ko_counts_cov50_all.tsv --label cov50`.
- 2,865 libraries (KO table, blanks dropped), rarefied to 250,000 KO reads (seed 0); proportions + 0.5/N
  over the real-model KOs.
- Per transform (CLR, GE), per-KO effects are in units of the KO's within-cruise library SD under that
  transform:
  - HiSeq 4000 − NovaSeq 6000 within the four mixed cruises (surface);
  - era: NovaSeq 6000 − NovaSeq X cruise means (modelled surface libraries);
  - HLY2202 same-DNA rerun (control).
- Summaries: RMS over KOs; share with |effect| > 0.5 SD; Spearman with KO mean abundance (the dropout
  signature); cosine with the R45 F1 loading.
- Spread coupling: per KO, the within-cruise Spearman correlation with library KO-profile spread (SD of
  CLR; R45), as the median over all KOs and over the 500 most abundant.
- Outputs: `data/results/marker_norm/cov50/` (`summary.csv`, `ko_effects.csv`, `markers.csv`,
  `summary.txt`).

| Transform | Contrast | RMS (SD) | Share \|effect\| > 0.5 SD | ρ with KO mean abundance | Cosine with F1 |
|---|---|---:|---:|---:|---:|
| CLR | HiSeq − NovaSeq 6000 | 0.415 | 0.247 | +0.544 | +0.633 |
| CLR | Era (6000 − X) | 0.362 | 0.170 | +0.356 | +0.596 |
| CLR | Rerun (control) | 0.211 | 0.026 | −0.035 | −0.045 |
| GE | HiSeq − NovaSeq 6000 | 0.343 | 0.152 | +0.246 | +0.077 |
| GE | Era (6000 − X) | 0.329 | 0.131 | +0.261 | +0.467 |
| GE | Rerun (control) | 0.213 | 0.027 | +0.077 | +0.112 |

Spread coupling (median within-cruise ρ of a KO with library spread):

| Transform | All KOs | 500 most abundant KOs |
|---|---:|---:|
| CLR | +0.151 | +0.425 |
| GE | −0.065 | −0.030 |

Remaining per-KO effects vs KO properties (Spearman):

| Effect | Eukaryote share of genes | Prevalence | GC | Mean abundance | Mean effect |
|---|---:|---:|---:|---:|---:|
| CLR HiSeq | −0.01 | +0.56 | −0.25 | +0.54 | +0.114 |
| GE HiSeq | +0.05 | +0.25 | −0.14 | +0.25 | −0.187 (markers −0.006) |
| CLR era | −0.22 | +0.43 | −0.22 | +0.36 | +0.080 |
| GE era | −0.22 | +0.33 | −0.18 | +0.26 | +0.003 |

**Conclusions**
1. **The CLR does amplify the dropout pattern.**
   - Under GE, the HiSeq − NovaSeq effect no longer resembles axis 1 (cosine 0.63 → 0.08).
   - Its common-KO signature halves (ρ 0.54 → 0.25), and its size shrinks (RMS 0.42 → 0.34 SD; share
     > 0.5 SD 25% → 15%).
   - Common KOs stop co-varying with profile spread (top-500 median ρ +0.43 → −0.03).
   - So most of the platform's axis-1 signature in R39–R45 is the CLR redistributing dropouts onto common
     KOs.
2. **The rest of the HiSeq effect is mostly a uniform shift.**
   - Non-marker KOs are about 0.19 SD lower per genome equivalent on HiSeq; the markers themselves are
     −0.006.
   - That is, the conserved, highly mappable ribosomal proteins take a slightly larger share of KO reads on
     HiSeq.
   - This fits the R45 picture (lower mapping; fewer distinct genes detected) and stays a technical offset.
   - It is still larger than the same-DNA rerun floor (RMS 0.34 vs 0.21).
3. **The era contrast changes less** (cosine with F1 0.60 → 0.47; RMS 0.36 → 0.33). It also mixes region,
   year and biology. Part of the axis-1 era offset is therefore not a CLR artefact.
4. **Next steps**
   - Re-run on the coverage-0% counts when they arrive, to separate the breadth filter from the transform:
     `scripts/marker_normalisation.py --counts <cov0 table> --label cov0`.
   - Then consider GE targets for the real model, with the residual PCA re-run to see whether axis 1
     persists. Eukaryotic KOs under GE move with the eukaryote : prokaryote ratio, so a second, eukaryotic
     reference may be needed for those.

**Note: filter size and early processing (to test on Tara).**
- Filter pore size and early sample processing may drive part of the Bacteria : Eukaryota differences.
- We have no per-sample filter or extraction methods for Bio-GO-SHIP (R35/R36), and R45 found the platform
  offset unrelated to the eukaryote fraction. So this cannot be tested here.
- Tara Oceans does record size fractions and protocols. Test there once the Tara read counting on Morrill
  finishes (still downloading, 2026-11-20). Questions:
  - Does the size fraction shift the eukaryote fraction?
  - Does it shift KO-profile spread / GE profiles?
  - Does it move samples along our axes?

## R47. Coverage-0 vs coverage-50 counts: merge and test set up; cov50 baseline (2026-11-20; cov0 pending)

**Hypothesis.**
- Under the 50% breadth-of-coverage filter, a library with fewer unique fragments covers fewer rare genes
  far enough to count them.
- Their KOs then go to zero, the KO profile becomes less even, and through the CLR axis 1 rises (R45/R46).
- If the HiSeq − NovaSeq gap in unevenness and axis 1 mostly disappears in the unfiltered (cov0) counts,
  the filter is what converts library complexity into missing KOs. Duplicate rates would then only explain
  why some libraries have lower complexity.

**Merge.** `scripts/merge_ko_counts.py`, run on Morrill by `cluster/merge_ko_counts.sbatch`:
- `--export=ALL,ROOT=<folder with the per-library folders>`.
- It sums each library's `level_<n>_<rank>_ko_counts.tsv` over taxa, using the lowest level present.
- It writes `ko_counts_cov0_all.tsv(.gz)` in the layout of `ko_counts_cov50_all.tsv`, plus
  `ko_counts_cov0_all_domain_counts.tsv` from `level_1_*_abundance.tsv`.
- Validation on the local cov50 per-taxon stash (`data/targets/taxa/ko_by_taxon`, 1,417 libraries also in
  `ko_counts_cov50_all.tsv`):
  - every library with a level-2 file reproduces its cov50 row (relative difference < 1e-5, from
    printing; the output now uses 10 significant digits);
  - the 52 libraries with only an order-level file lose about 3% of counts (ratio 0.906–0.994, median
    0.971), because genes without an order assignment drop out. The script flags such libraries.
- The cov0 folders were not on the Mac, so the merge has to run on Morrill.

**Test.** `scripts/coverage_threshold_test.py`. It shares the rarefaction helper (now gzip-aware) with
`scripts/marker_normalisation.py`.
- Per coverage table, libraries are rarefied to 250,000 KO reads (seed 0) over the 3,933 real-model KOs.
- Per-library metrics:
  - profile spread (SD of CLR), undetected KOs, Shannon;
  - axis-1 score on CLR and on genome equivalents (the cov50 F1 loading and target SDs, so both tables
    are on one scale);
  - mapping rate.
- Contrasts:
  - HiSeq 4000 − NovaSeq 6000 within the four mixed cruises: pooled, with permutation p (1,999), also in
    within-cruise SD units, and per cruise;
  - the HLY2202 same-DNA rerun (control).
- Comparison: the share of the cov50 gap left at cov0, and per-library cov0 − cov50 changes by platform.
- Outputs: `data/results/coverage_threshold/` and `figures/output/coverage_threshold_gap.png`.

**cov50 baseline** (HiSeq − NovaSeq 6000, mixed cruises, surface; rerun = HLY2202 same DNA, 21 pairs):

| Metric | Gap | Within-cruise SD units | p | I09N | AE1319 | BVAL46 | NH1418 | Rerun mean (SD) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Profile spread (SD of CLR) | +0.1250 | +0.96 | 0.0005 | +0.093 | +0.152 | +0.037 | +0.338 | +0.006 (0.153) |
| Undetected KOs | +0.0412 | +1.12 | 0.0005 | +0.029 | +0.033 | +0.017 | +0.143 | +0.005 (0.044) |
| Shannon | −0.0594 | −1.01 | 0.0005 | −0.062 | −0.055 | +0.029 | −0.090 | +0.001 (0.038) |
| Axis-1 score, CLR | +17.82 | +0.77 | 0.0005 | +12.99 | +26.08 | +12.54 | +39.45 | −0.67 (18.11) |
| Axis-1 score, genome equivalents | +1.08 | +0.05 | 0.81 | +3.39 | +8.00 | −5.59 | −21.02 | +1.77 (10.32) |
| Mapped reads / pair | −0.187 | −1.29 | 0.0005 | −0.090 | −0.304 | −0.275 | −0.549 | −0.062 (0.141) |

**Readout.**
- On the genome-equivalent scale, the platform gap on axis 1 is already gone at cov50 (+0.05 SD, p 0.81),
  consistent with R46.
- The undetected-KO gap (+1.1 SD) is in the counts themselves; the transform cannot change it. That is the
  number cov0 has to remove for the breadth-filter explanation to hold.
- Note on the mapping rate: the domain abundance counts are also coverage-filtered (count_biogoship_taxa
  uses the same threshold), so at cov0 it measures all mapped reads.

### R47 results: cov0 arrived (Google Drive zips, 175 libraries) (2026-11-20)

**Data.**
- 13 Google Drive zips (`drive-download-20261120T231310Z-1-001…013.zip`, 6.0 GB) held 175 complete
  libraries, every level present.
- The zips had no `level_1_*_abundance.tsv`, so there is no cov0 mapping rate.
- `scripts/merge_ko_counts.py` (now reads zips) copied each library's level-1 (domain) KO file gzipped to
  `data/targets/cov0_stash/<library>/` (46 MB) and merged them into `ko_counts_cov0_all.tsv.gz`
  (175 libraries × 19,340 KOs, against 18,836 at cov50). The zips were then deleted.
- Command:
  `python scripts/merge_ko_counts.py figures/output/drive-download-*.zip --stash data/targets/cov0_stash --out ko_counts_cov0_all.tsv`.
- Libraries by cruise: I09N-2016 HiSeq 4000 4 / NovaSeq 6000 95; NH1418-2014 HiSeq 5 (2 surface) /
  NovaSeq 7; AE1319-2013 HiSeq 1; HLY2202-2022 50 (all on flowcell HMY7WDSX7); AMT28-2018 7; P18-2016 6.
  BVAL46-2011 is absent.

**Checks.**
- No KO count is lower at cov0 than at cov50 in any library: the filter only removes reads.
- KO-assigned reads cov0 / cov50: median 2.19× (range 1.40–78.3).
- KOs detected per library: median 14,580 at cov0 vs 7,063 at cov50.

**Test.** `scripts/coverage_threshold_test.py`: 158 libraries in both tables and above 250,000 KO reads.
- Mixed cruises with both platforms: I09N-2016 (HiSeq 4 vs NovaSeq 91, surface) and NH1418-2014 (HiSeq 2
  vs 7, surface).
- No rerun control: cov0 has only one HLY2202 flowcell.
- Within-cruise SDs are computed on these same 158 libraries, separately for each coverage.

| Metric | cov50 gap (raw) | cov50 (within-cruise SD) | p | cov0 gap (raw) | cov0 (within-cruise SD) | p | Share of the SD-scaled gap left at cov0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Profile spread (SD of CLR) | +0.2054 | +1.68 | 0.002 | +0.0899 | +1.64 | 0.012 | 98% |
| Undetected KOs | +0.0838 | +1.98 | 0.006 | +0.0140 | +1.91 | 0.008 | 96% |
| Shannon | −0.0796 | −1.68 | 0.0005 | −0.0408 | −1.32 | 0.019 | 79% |
| Axis-1 score, CLR | +31.82 | +1.52 | 0.004 | +14.32 | +1.12 | 0.058 | 74% |
| Axis-1 score, genome equivalents | −10.16 | −0.45 | 0.36 | −4.96 | −0.29 | 0.56 | 65% |

Per cruise:
- cov50: I09N spread +0.065, undetected +0.015, F1 +10.6; NH1418 spread +0.485, undetected +0.221,
  F1 +74.4.
- cov0: I09N spread +0.025, undetected +0.003, F1 +4.2; NH1418 spread +0.220, undetected +0.035,
  F1 +34.6.

Per-library change from cov50 to cov0 (mixed cruises, surface; mean):

| Platform | Spread | Undetected | Shannon | F1 (CLR) | F1 (GE) |
|---|---:|---:|---:|---:|---:|
| HiSeq 4000 | −0.423 | −0.135 | +0.090 | −60.5 | −0.19 |
| NovaSeq 6000 | −0.346 | −0.085 | +0.056 | −47.7 | −2.01 |

Figure: `figures/output/coverage_threshold_gap.png`. Three panels: undetected KOs, profile spread and the
CLR axis-1 score, each centred within cruise, by platform, at ≥ 50% coverage and with no filter.

**Conclusions** (two cruises and 6 surface HiSeq libraries, so low power)
1. **The 50% filter amplifies the gap but does not create it.**
   - Without the filter, the raw HiSeq − NovaSeq gaps shrink by half or more (undetected KOs
     0.084 → 0.014; F1 31.8 → 14.3).
   - The HiSeq libraries also lose more to the filter (undetected −0.135 vs −0.085).
   - But library-to-library variation shrinks in step. Relative to it, the HiSeq libraries remain about
     as unusual: unevenness keeps 96–98% of its gap, CLR axis 1 keeps 74%.
   - So the coverage threshold is not what converts complexity into missing KOs. The deficit, fewer
     distinct KOs per read, is already present in the mapped reads.
2. **This points to the reads themselves: library complexity, i.e. fewer distinct fragments.**
   - The duplicate-rate / library-complexity check is now the main test, not just a confirmation:
     `fastp` duplication or `samtools markdup` on subsampled reads for the mixed cruises, on Morrill.
   - Other read-level options are still open: a different library-prep kit or input amount per platform
     batch.
3. **On the genome-equivalent scale there is no significant platform gap at either coverage**
   (−0.45 / −0.29 SD, p ≥ 0.36), consistent with R46. The transform choice matters more than the
   coverage threshold for axis 1.
4. **Limits.**
   - The gap rests mostly on NH1418 (2 surface HiSeq libraries, one extreme); I09N's gap is small at
     both coverages.
   - No same-DNA control at cov0.
   - Re-run with the remaining cov0 libraries if more arrive (e.g. BVAL46, AE1319 NovaSeq, the second
     HLY2202 flowcell).

## R48. Toward domain-split targets: denominators checked; cov0 dropped (2026-11-20)

**cov0 dropped.**
- R47 answered what cov0 could answer: the 50% filter amplifies the platform gap but does not create it.
  On the larger cruise (I09N) the gap is small at both coverages.
- The remaining question, library complexity, is in the reads and needs duplicate rates, not a different
  counting threshold.
- The coverage threshold rescales all libraries rather than forming its own variance slice, so it is not
  a predictability bucket.
- `data/targets/cov0_stash/` (46 MB) and `ko_counts_cov0_all.tsv.gz` are kept in case.

**How much of the platform problem the GE transform fixes (R46 numbers).**
- It removes the axis-1 signature of the HiSeq − NovaSeq effect (cosine with F1 0.63 → 0.08). The
  within-cruise axis-1 gap is not significant at cov50 (+0.05 SD, p 0.81, R47 baseline).
- Excess platform variance above the same-DNA rerun floor (RMS² − 0.211²) drops from 0.128 to 0.073 SD²,
  about 43% less. The rest is mostly a uniform per-genome offset.
- The era contrast keeps part of its axis-1 shape (cosine 0.60 → 0.47).
- So: the CLR-driven artefact is fixed; about half of the platform offset is not.

**Denominator checks** (cov50 counts, 3,066 non-blank libraries; Spearman).
- Domain classes come from the sampled reference genes (`ko_gc.csv`): eukaryote share of a KO's genes
  = Eukaryota / (Eukaryota + Bacteria + Archaea). Over the 3,933 targets:

  | Cut | Prokaryote-specific | Eukaryote-specific | Shared |
  |---|---:|---:|---:|
  | 5% / 95% | 1,762 | 420 | 1,751 |
  | 10% / 90% | 2,143 | 485 | 1,305 |
  | 2% / 98% | 1,358 | 335 | 2,240 |

  In the full table (18,836 KOs) at 5% / 95%: 5,534 prokaryote-specific and 8,818 eukaryote-specific.
- **Prokaryote-type ribosomal proteins** (56 KOs in the full table), pooled: ρ 0.987 with Bacteria domain
  reads. A good prokaryotic denominator.
- **Eukaryote-type cytosolic ribosomal proteins** (77 KOs), pooled:
  - median 3,579 reads per library; 634 libraries < 1,000; 1,762 < 5,000;
  - ρ 0.541 with Eukaryota reads; as shares of KO reads vs the eukaryotic read fraction, ρ 0.056;
  - log(eukaryotic RP / prokaryotic RP) vs log(Eukaryota / prokaryote reads): ρ 0.047.
  - Not usable as a denominator. A likely reason: the eukaryotic references are largely transcript
    assemblies, so intron-rich, multi-copy ribosomal-protein genes map poorly from DNA. Untested.
- **All eukaryote-specific KOs** (8,818), pooled:
  - ρ 0.961 with Eukaryota reads;
  - log(eukaryote-specific / prokaryote-specific KO reads) vs the domain read ratio: ρ 0.890;
  - median 44,788 reads per library (1% 311, 5% 4,621, 10% 7,840); < 1,000: 52 libraries;
    < 5,000: 172; < 10,000: 413.

**Proposed targets (not built yet).**
- **A. Prokaryote-specific KOs:** log copies per prokaryote genome equivalent (pooled prokaryote
  ribosomal proteins).
- **B. Eukaryote-specific KOs:** within-eukaryote composition, i.e. a log-ratio against the pooled
  eukaryote-specific KO reads (or their geometric mean).
  - This is not "per eukaryotic genome", which ploidy, genome size and ribosomal-protein copy number make
    fuzzy.
  - Libraries with < 5,000 eukaryote-specific reads are dropped from B rather than imputed.
- **C. Eukaryotic amount:** one target per sample, log(eukaryote-specific / prokaryote-specific KO reads).
  Filter size and early processing would act here (Tara test, R46 note).
- **Shared KOs:** whole-community model only, or dropped from A/B.
- **Next:** build A/B/C targets and re-run the baseline model, residual PCA and variance budget per part.
  That shows which predictability buckets remain once the CLR artefact is gone.

## R49. Plastid / organellar leak into the prokaryote genome-equivalent denominator (2026-11-20)

**Question.**
- Prokaryote-type ribosomal proteins also occur in plastids and mitochondria: plastid-encoded copies (red
  plastids nearly all, green plastids fewer), and nucleus-encoded organellar copies that annotate to the
  bacterial-type KOs.
- Plastids are high-copy, so in eukaryote-rich samples the R46 "prokaryote genome equivalents" may be
  inflated.
- Cyanobacteria carry the same genes but are genuine prokaryotes.

**Method.** Script: `scripts/plastid_marker_check.py`. Outputs: `data/results/plastid_markers/`
(`markers.csv`, `library_denominators.csv`, `deciles.csv`, `summary.txt`) and
`figures/output/plastid_marker_check.png`.
1. **Literature split of the 54 markers.** Plastid-encoded in at least one major algal lineage (union of
   green and red plastid genomes): L1–L6, L11, L7/L12, L13, L14, L16, L18–L24, L27, L29, L31–L36, S2–S20
   (45 KOs). Rarely or never plastid-encoded: L9, L10, L15, L17, L25, L28, L30, S1, S21 (9 KOs).
2. **Empirical split.** From the phylum-level per-taxon KO counts (`data/targets/taxa/ko_by_taxon`,
   1,356 libraries):
   - per marker KO, the share of its reads assigned to eukaryotic taxa and to Cyanobacteriota, over all
     libraries and over the eukaryote-rich quartile (domain eukaryote fraction ≥ 0.058);
   - "clean" = eukaryote-assigned share below 2% in that quartile.
3. **Divergence on all libraries.** cov50 table: 3,005 non-blank libraries with > 100k KO reads.
   - Per subset, the log denominator is the mean of log(count + 0.5) over its KOs.
   - Differences between subsets are tested against the eukaryote fraction and the Cyanobacteriota read
     share (Spearman, overall and within cruise), and summarised by eukaryote-fraction decile.
4. **Truth check.** Against the log mean of the same 54 markers' reads assigned to Bacteria + Archaea only,
   for the 1,337 libraries with per-taxon counts.

**Per-marker eukaryote-assigned share** (eukaryote-rich quartile; full table in `markers.csv`):
- Range 0.6% (S21) to 27.1% (L6).
- Median 4.3% for the 9 non-plastid KOs; 9.1% for the 45 plastid-encodable KOs.
- Highest: L6 0.271, L20 0.221, L2 0.181, S14 0.180, L14 0.177, S12 0.175.
- Lowest: S21 0.006, L25 0.015, L29 0.029, L30 0.031, S1 0.035, L32 0.037.
- Only 2 KOs (S21, L25) fall below 2%.
- Cyanobacteriota share of marker reads (all libraries): 0.1–19% per KO.

**Divergence** (log difference between denominators; top − bottom eukaryote-fraction decile):

| Comparison | ρ with eukaryote fraction | Within cruise | ρ with cyano share | Within cruise | Top − bottom decile |
|---|---:|---:|---:|---:|---:|
| Plastid-encodable − non-plastid (literature) | +0.240 | +0.242 | +0.302 | +0.325 | +0.052 (+5%) |
| Eukaryote-carried − clean (empirical, 52 vs 2 KOs) | −0.073 | +0.029 | +0.211 | +0.171 | −0.078 (−8%) |
| All 54 − non-plastid 9 | +0.240 | +0.242 | +0.302 | +0.325 | +0.044 (+4%) |
| All 54 − prokaryote-assigned (n 1,337) | +0.552 | +0.419 | +0.742 | +0.685 | +0.201 (+22%) |
| Clean 2 − prokaryote-assigned | +0.228 | +0.199 | +0.247 | +0.157 | +0.276 (+32%) |
| Non-plastid 9 − prokaryote-assigned | +0.294 | +0.238 | +0.450 | +0.417 | +0.159 (+17%) |

The full decile table is in `deciles.csv`. For all 54 vs prokaryote-assigned, the inflation rises to
+0.09 to +0.13 log by deciles 2–6 (eukaryote fraction about 0.017–0.036) and reaches +0.20 in the top decile
(median eukaryote fraction 0.27).

Figure: `figures/output/plastid_marker_check.png`. Left: plastid-encodable − non-plastid log denominator vs
eukaryote fraction. Right: all 54 markers − prokaryote-assigned reads only. Each panel has binned means.

**Conclusions**
1. **The leak is real.** In eukaryote-rich libraries the 54-marker genome equivalents are inflated by
   about 22% (top decile), and by about 10–13% already at modest eukaryote fractions (2–4%).
   - The inflation tracks eukaryote fraction (ρ 0.55; 0.42 within cruise).
   - It tracks the cyanobacterial share even more closely (ρ 0.74). Eukaryotic phytoplankton and
     cyanobacteria co-vary, so eukaryote-assigned marker reads rise with both.
   - Uncorrected, this would push eukaryote-rich samples' prokaryotic KOs down, along the axis-1
     (bacteria vs eukaryotes) direction.
2. **A literature split does not fix it.**
   - The 9 non-plastid KOs still carry eukaryotic reads (median 4.3%), because nucleus-encoded
     organellar (mostly mitochondrial) ribosomal proteins annotate to the same KOs.
   - The non-plastid subset still inflates by 17%, against 22% for all 54.
   - Only 2 KOs are clean, too few for a stable denominator (+32% noise-driven bias).
3. **Fix: count only prokaryote-assigned reads of the markers**, using the reference taxonomy, so
   eukaryotic plastid, mitochondrial and nuclear copies drop out by construction.
   - Cyanobacterial reads stay in, correctly: they are prokaryotes. Plastid reads mis-binned to
     Cyanobacteriota in the reference cannot be detected this way.
   - This needs per-taxon KO counts for every library. Local per-taxon counts cover 1,356 of about 3,000
     libraries (phylum level).
   - For the rest, only `level_1_domain_ko_counts.tsv` per library is needed from the Drive outputs (cov50
     run). `scripts/merge_ko_counts.py --stash` can pull exactly those from zips. No re-mapping is needed.

## R50. cov50 per-taxon KO counts for all libraries; level-1 vs level-2 sums; R47 rerun like-for-like (2026-11-20)

**Data.**
- 38 Google Drive zips (`coverage_50-20261120T235254Z-1-001…038.zip`, 42 GB) held all 3,098 cov50 library
  folders, each with KO / OG / Pfam counts per taxon at all seven levels (domain to species). No
  `level_*_abundance.tsv` files (domain read counts are already in `data/targets/taxa/domain_counts.tsv`).
- `scripts/merge_ko_counts.py` (now with `--stash-levels`, parallel `--jobs`, and `--level`) copied each
  library's `level_1_domain_ko_counts` and `level_2_phylum_division_ko_counts` gzipped to
  `data/targets/taxa/ko_by_taxon/<library>/`: 4,817 new files; the 1,379 level-2 files already there were
  kept. All 3,098 libraries now have both levels (folder 1.7 GB, gitignored). OG, Pfam and levels 3–7 were
  not kept; the zips were deleted.
- Command: `python scripts/merge_ko_counts.py coverage_50-*.zip --stash data/targets/taxa/ko_by_taxon --stash-levels 1 2 --out <table>`.

**Which level the official table sums.**
- `ko_counts_cov50_all.tsv` = sum over taxa of the **level-2 (phylum)** file: rebuilt for all 3,081
  libraries, max absolute difference 0.0003 (printing), same 18,836 KOs.
- The **level-1 (domain)** file holds more counts: per library level-1 / level-2 median 1.082 (5–95%
  1.02–1.17, max 1.44; 3,081 libraries). The extra rows are a taxon "Unclassified" (genes with no
  taxonomy) plus genes with a domain but no phylum. It also has 237 extra KOs (6e-7 of reads).
- The 17 libraries not in the official table have zero level-2 KO counts (blanks, mocks, empty runs:
  EXTRBLANKBGSP8, ZYMOMOCKBGSP5/9, 9 GOMECC, 4 GSI05, 2 GS-A16N).
- The merge script now defaults to `--level 2` (reproduces the official table); a library without that
  level falls back to the nearest and is flagged. The R47 text "lowest level present" no longer applies.
- `data/targets/ko_counts_cov50_level1.tsv` (3,098 × 19,073, gitignored): the cov50 table at level 1.

**Consequence for R47.** The cov0 table was built from level-1 files (only those were stashed), so R47
compared cov0 at level 1 with cov50 at level 2. On the 175 cov0 libraries: cov0 / cov50 KO reads median
1.948 like-for-like (range 1.25–78.3) instead of 2.19 (level-1 / level-2 at cov50 here: 1.094); KOs detected
median 7,350 (cov50 level 1) vs 7,063 (level 2) vs 14,617 (cov0). Still no cell with fewer counts at cov0.

**R47 rerun, both tables at level 1** (`scripts/coverage_threshold_test.py --tag level1 --tables
cov50=data/targets/ko_counts_cov50_level1.tsv cov0=ko_counts_cov0_all.tsv.gz`; 158 libraries; HiSeq 4000 −
NovaSeq 6000, I09N-2016 and NH1418-2014, surface):

| Metric | cov50 L1 gap (raw) | SD units | p | cov0 gap (raw) | SD units | p | Left at cov0 | R47 (L2 vs L1) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Profile spread (SD of CLR) | +0.1998 | +1.69 | 0.002 | +0.0899 | +1.64 | 0.012 | 97% | 98% |
| Undetected KOs | +0.0781 | +1.95 | 0.006 | +0.0140 | +1.91 | 0.008 | 98% | 96% |
| Shannon | −0.0770 | −1.65 | 0.0005 | −0.0408 | −1.32 | 0.019 | 80% | 79% |
| Axis-1 score, CLR | +31.48 | +1.57 | 0.004 | +14.32 | +1.12 | 0.058 | 71% | 74% |
| Axis-1 score, genome equivalents | −10.67 | −0.47 | 0.31 | −4.96 | −0.29 | 0.56 | 62% | 65% |

Per cruise, cov50 level 1: I09N spread +0.065, undetected +0.012, F1 +11.1; NH1418 spread +0.469,
undetected +0.210, F1 +72.2. Mapped reads / pair gap −0.255 (−1.62 SD, p 0.003).

Per-library change cov50 L1 → cov0 (mixed cruises, surface; mean):

| Platform | Spread | Undetected | Shannon | F1 (CLR) | F1 (GE) |
|---|---:|---:|---:|---:|---:|
| HiSeq 4000 | −0.410 | −0.127 | +0.085 | −60.7 | −0.15 |
| NovaSeq 6000 | −0.336 | −0.081 | +0.052 | −48.1 | −2.42 |

Outputs: `data/results/coverage_threshold_level1/`, `figures/output/coverage_threshold_gap_level1.png`
(same three panels as R47). **R47's conclusions stand unchanged** (filter amplifies, does not create, the
gap; within 1–4 percentage points).

**Next.** The prokaryote-assigned marker denominator (R49 fix) can now be built for every library from the
level-2 files (`Bacteria;*` + `Archaea;*` rows, consistent with the official table).

## R51. Domain-split targets A / B / C / D: which satellite input predicts each, variance budgets, two-stage model (2026-11-20)

**Targets** (`scripts/domain_split_targets.py` → `data/targets/domain_split/`, gitignored). The same 1,990
samples and 2,018 cov50 libraries as `ko_clr.csv` (rebuilt CLR matches it to 5e-6). Each of the 3,933 model
KOs is split per library into reads assigned to Bacteria + Archaea and to Eukaryota (phylum-level per-taxon
files, R50). Proportions are per library, averaged per sample, with pseudocount 0.5 / sample total.
- Classes (5% / 95% eukaryote share of reference genes, R48): 1,759 prokaryote-specific, 420
  eukaryote-specific, 1,700 shared, plus 54 ribosomal-protein markers (the denominator, not targets).
- **A** prokaryote-specific KOs (1,759): log copies per prokaryote genome equivalent. The denominator is
  the mean log of the 54 markers' **prokaryote-assigned** reads (the R49 fix). Prokaryote-assigned share
  of marker reads: median 0.878 over samples. All-54 minus prokaryote-assigned log denominator: median
  0.127 (ρ 0.10 with the KO-based eukaryote fraction).
- **B** eukaryote-specific KOs (420): log share of pooled eukaryote-specific reads (all 8,823 such KOs in
  the full table). 84 samples with < 5,000 pooled eukaryote-specific reads dropped (1,906 left).
- **C** one column: log(pooled eukaryote-specific / pooled prokaryote-specific reads; 8,823 / 5,544 KOs).
- **D** shared KOs (1,700): the current whole-community CLR, unchanged.
- **Dp / De**: the shared KOs' prokaryote-assigned part per prokaryote genome equivalent (1,561 KOs with
  ≥ 90% prevalence) and eukaryote-assigned part within eukaryotes (312 KOs; 1,906 samples).
- **The shared KOs are almost all prokaryotic reads.** Over the 3,933 model KOs 97.2% of reads are
  prokaryote-assigned, 2.7% eukaryote, 0.1% other. Per shared KO, prokaryote-assigned share median 0.978
  (IQR 0.933–0.990); 1,372 are > 90% prokaryotic, 7 are < 10%.

**Models** (`scripts/domain_split_models.py` → `data/results/domain_split/`; figures
`figures/output/domain_split_skill.png`, `domain_split_budget.png`). R43 tiles (1,215; 1,169 for B / De),
one cruise held out, ridge with per-KO GCV alpha. Feature sets: in-situ temperature spline (not a
satellite input), satellite SST spline, 16 OLCI band means (splines), OLCI + SST, 32 Granite PCs, and the
baseline (32 Granite PCs + temperature spline). Each new target is paired with the current CLR on the same
KOs and samples.

Variance explained, cruise held out (raw units, pooled over tiles and KOs; KOs standardised in brackets):

| Target | Temperature (in situ) | Satellite SST | OLCI bands | OLCI + SST | Granite PCs | Baseline |
|---|---:|---:|---:|---:|---:|---:|
| Current CLR, all 3,933 KOs | 11.6 (3.4) | 12.2 (4.7) | −2.1 (−1.3) | 12.4 (7.0) | **17.2 (12.7)** | 16.0 (11.1) |
| A prokaryote-specific per GE | 6.6 (1.8) | 6.2 (0.9) | −9.6 (−9.7) | 4.5 (−2.1) | **11.0 (3.5)** | 10.5 (2.5) |
|   same KOs, current CLR | 8.7 (1.3) | 8.3 (2.1) | −7.5 (−3.7) | 7.5 (4.3) | **13.7 (10.7)** | 13.1 (9.3) |
| B eukaryote-specific within eukaryotes | 7.7 (3.0) | 8.7 (4.5) | 0.2 (−1.4) | 9.4 (4.9) | **10.9 (6.4)** | 8.7 (3.8) |
|   same KOs, current CLR | 27.7 (25.2) | 30.0 (27.8) | 10.0 (8.5) | **31.3 (28.7)** | 30.6 (28.2) | 28.2 (25.5) |
| **C eukaryote / prokaryote amount** | 48.5 | 50.7 | 29.5 | **56.2** | 54.9 | 52.5 |
| D shared, whole-community CLR | 4.8 (0.8) | 5.4 (2.2) | −4.0 (−1.3) | 6.1 (4.9) | **12.9 (11.3)** | 11.4 (9.8) |
| Dp shared, prokaryote part per GE | 6.6 (3.0) | 5.8 (1.8) | −10.8 (−10.0) | 5.4 (−1.2) | **13.8 (4.9)** | 13.2 (3.9) |
|   same KOs, current CLR | 5.0 (0.6) | 5.0 (1.8) | −5.4 (−1.5) | 5.2 (4.6) | **12.7 (11.3)** | 11.6 (9.9) |
| De shared, eukaryote part within eukaryotes | 22.7 (13.5) | 23.2 (14.5) | 8.1 (3.9) | **23.8 (15.0)** | 23.7 (15.3) | 23.3 (14.3) |
|   same KOs, current CLR | 8.6 (0.6) | 10.7 (3.4) | 3.8 (5.0) | 14.1 (10.2) | **18.2 (15.6)** | 16.0 (14.4) |

Best feature set per KO (share of KOs; "none" = no set above R² 0.05):

| Target | Granite PCs | Baseline | Temperature | Satellite SST | OLCI + SST | OLCI | None |
|---|---:|---:|---:|---:|---:|---:|---:|
| Current CLR, all | 30% | 17% | 8% | 6% | 7% | 4% | 28% |
| A | 12% | 17% | 13% | 6% | 3% | 2% | 48% |
| B | 17% | 10% | 7% | 9% | 9% | 1% | 49% |
| B KOs, current CLR | 22% | 3% | 11% | 22% | 36% | 0% | 5% |
| D | 34% | 17% | 8% | 4% | 4% | 4% | 30% |
| Dp | 16% | 19% | 12% | 7% | 3% | 1% | 43% |
| De | 14% | 10% | 15% | 11% | 15% | 1% | 34% |

C as a covariate (spline added to the baseline). "Observed" is an oracle, since C is not seen from space;
"predicted" is the out-of-fold baseline prediction of C, a legitimate two-stage model:

| Target | Baseline | + C predicted | + C observed (oracle) |
|---|---:|---:|---:|
| Current CLR, all 3,933 KOs | 16.0 (11.1) | **19.4 (14.5)** | 26.1 (22.5) |
| A | 10.5 (2.5) | 9.1 (0.7) | 10.8 (3.0) |
| D shared CLR | 11.4 (9.8) | 12.7 (12.2) | 16.5 (19.3) |
| Dp | 13.2 (3.9) | 12.7 (3.3) | 14.0 (4.4) |

Variance budget, baseline model (raw %, sequential as R43; within-cruise PCs from a standardised PCA of each
target's own residual):

| Target | Explained | Platform | Cruise offsets | Depth | Within PC1 | Within PC2 | Replicate noise | Remaining |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Current CLR, all | 16.2 | 3.8 | 13.1 | 2.3 | 7.5 | 3.7 | 27.0 | 26.4 |
| A | 10.8 | 2.9 | 16.9 | 4.9 | 5.2 | 3.0 | 25.9 | 30.4 |
| A KOs, current CLR | 13.4 | 3.1 | 15.7 | 2.5 | 4.9 | 3.8 | 24.6 | 31.9 |
| B | 8.9 | 3.0 | 9.8 | 1.6 | 6.0 | 4.2 | 57.6 | 8.8 |
| B KOs, current CLR | 28.6 | 6.3 | 7.1 | 1.3 | 23.5 | 1.5 | 29.8 | 1.9 |
| C | 52.9 | 7.8 | 3.8 | 0.9 | – | – | 6.5 | 28.0 |
| D | 11.6 | 3.2 | 14.2 | 2.4 | 6.5 | 5.6 | 27.0 | 29.6 |
| Dp | 13.5 | 2.6 | 15.6 | 5.2 | 5.1 | 3.5 | 23.5 | 31.1 |
| De | 23.5 | 1.9 | 8.3 | 2.3 | 11.5 | 3.1 | 45.1 | 4.3 |

(The current-CLR row reproduces R43: 16.2 / 3.8 / 13.1 / 2.3 / 7.5 / 3.7 / 27.0.)

**Reading.**
1. **The largest predictable bucket is the eukaryote : prokaryote amount (C): 53–56% explained**, with low
   replicate noise (6.5%) and small cruise offsets (3.8%). Its residual has a platform part (7.8%).
2. **Most of what the current CLR predicts in eukaryote-specific KOs is C.**
   - Under the CLR those 420 KOs score 28–31%. As composition within eukaryotes (B) they drop to 9–11%.
   - B is mostly replicate noise (58%): which eukaryotic genes are present, given how many eukaryotes, is
     barely reproducible between filters of the same water at these read depths.
   - B is the worst part, as expected.
3. **Prokaryotic gene content per genome (A, Dp) is only weakly predictable.**
   - Raw 10–14%. Standardised 2.5–4% against 9–11% under the CLR.
   - So much of the per-KO CLR skill on prokaryotic KOs was also the shared eukaryote-amount signal leaking
     in through the CLR denominator.
   - A and Dp carry more sequencing depth (4.9–5.2% against 2.5%) and larger cruise offsets (16–17%).
     Genome-equivalent normalisation exposes the per-cruise offsets that the CLR partly absorbed.
4. **Shared KOs.** They are 98% prokaryotic reads, so as a whole-community target (D) they behave like
   prokaryotic KOs (11–13%). Their eukaryotic part within eukaryotes (De) scores 23–24%, against 16–18% for
   the same KOs under the CLR, but with 45% replicate noise.
5. **Which input predicts best.**
   - **Granite PCs alone are best or tied on every target except C and the within-eukaryote targets.** There
     raw OLCI + SST is level or slightly ahead (C 56.2 vs 54.9; De 23.8 vs 23.7; B KOs under the CLR 31.3 vs
     30.6).
   - The eukaryote-side signal is essentially a temperature / SST signal. Satellite SST alone gives C 50.7%
     and B-CLR 30.0%, and the OLCI bands add a little on top.
   - **OLCI bands alone fail on prokaryotic targets** (−10% on A / Dp): they overfit without SST.
   - Adding in-situ temperature to the Granite PCs (the baseline) is slightly worse than the PCs alone
     everywhere (all KOs 16.0 vs 17.2). The temperature spline overfits under cruise hold-out.
   - Per KO the skill is spread: on the full set Granite wins 30% of KOs, baseline 17%, temperature / SST /
     raw satellite 21% together, none 28%.
6. **The two-stage model raises the headline.** Predicting C first and feeding the prediction to the KO
   model: all 3,933 KOs **16.0% → 19.4%** raw (11.1 → 14.5% standardised), D 11.4 → 12.7%. The oracle with
   observed C gives 26.1%.
   - C is a strong, low-noise, one-dimensional signal that ridge on 40 features does not recover KO by KO.
   - On A and Dp, C adds nothing (9.1 / 12.7%), consistent with genome-equivalent targets having C removed.

**So the predictability buckets are:**
- eukaryote amount (C, about 55%);
- the whole-community CLR on top of C (two-stage, +3.4 points);
- per-genome prokaryotic content (weak, about 11% raw, dominated by cruise offsets, noise and an unknown
  within-cruise rest);
- within-eukaryote composition (weak, noise-dominated).

**Next options:**
- the two-stage model with the best C predictor (OLCI + SST or Granite + SST);
- splitting C's residual by platform (7.8%) and by filter size (Tara);
- whether the A / Dp cruise offsets are the per-genome platform offset of R46.

## R52. Water-mass fields step returned: along_track_daily.csv checked (2026-11-21)

`transfer/results/water_mass/along_track_daily.csv` (29 MB, not committed; from the `STEPS=fields` re-run of
`transfer/cluster/water_mass.sbatch`). Long format: sample_id, day (0 … −60), lat, lon, variable, value.
- **Complete:** 2,014 samples × 61 days × 5 variables (sst, chl, kd490, mld, sea-floor depth) = 614,270 rows.
  The interpolation fix worked: every sample has all 61 days (the first run had only day 0).
- **Missing values: 2.3%.** They grow backward in time (day −10: 9 samples per variable; day −60: about 95).
  At day 0 only chl / kd490 are missing, for 2 samples.
- **The missing days are tracks stuck at the coast.** 109 samples have days with every variable missing.
  On those days the particle moves a median 0.001° per day, against 0.145° on valid days. 85% of the runs
  are contiguous and 73% reach day −60. They are mostly coastal: OA_GOM 69, GS_HLY 21, I9N 14, GS_A16 9,
  GS_I09 8, GS_P02 7.
  - `stuck_at_end` in trajectories.csv flags none of them, because the particle is dropped or frozen
    off-grid rather than at the end.
  - Rule for the test: use the history only up to the first missing day, and flag the sample.
- **Sanity check at day 0:** OSTIA SST against in-situ temperature, original cruises r 0.996 (n 673, median
  |diff| 0.20 °C); FAIRe r 0.908 (n 1,153, median |diff| 0.20 °C).
- **First look at features:**
  - Shelf water (sea floor < 200 m) is crossed within 60 days by 289 of 2,014 samples; 120 sit on the shelf
    at day 0. For the 289, days since last on shelf: median 11 (IQR 0–28).
  - SST day 0 − day −30: median +0.30 °C (IQR −0.41 to +0.95).
  - log10 chlorophyll, day 0 over the mean of the previous 30 days: median −0.023 (IQR −0.091 to +0.041).
    So most water parcels see little bloom-scale change over a month.
- **Next:** `transfer/scripts/water_mass_test.py`.
  - Features: time since shelf, trajectory-mean and trend of SST and chlorophyll, mixed-layer light
    (MLD × KD490), origin region.
  - Tested against the R36 residual axes, the budget's "remaining" slice and the R51 targets (C, A).

## R53. Water-mass history as predictors: mostly null; 30-day mixed-layer depth helps the eukaryote amount (C) (2026-11-21)

**Script** `scripts/water_mass_test.py` → `data/results/water_mass/{features_sample,features_tile,added_value,axis_correlations}.csv`,
`summary.txt`; figure `figures/output/water_mass_test.png`.

**Features.**
- Computed per sample from `along_track_daily.csv` (R52), then averaged per tile.
- The 109 stuck tracks are truncated at the first missing day (valid days median 61, min 3).
- Tiles with every feature: 99.9%.
- The features:
  - Days since shelf (sea floor < 200 m; 61 = never). The 10th percentile is 25; the median is 61.
  - SST: change over 30 days, SD over 60 days.
  - Chlorophyll: 60-day mean log10, 30-day change, past peak above day 0.
  - MLD: 30-day mean log10, 30-day change.
  - Mixed-layer mean light fraction (1 − e^(−K·MLD)) / (K·MLD), K = KD490: at day 0 and as a 30-day mean.
  - 60-day displacement (median 376 km, 10–90%: 119–929 km) and the change in latitude from the track
    origin.
  - A control: displacement shuffled within cruise.

**Test 1: added value** (R43 tiles, cruise held out, baseline = 32 Granite PCs + temperature spline; each
feature added as a spline; gain in percentage points of explained variance, raw; KOs standardised in
brackets for the KO targets):

| Feature | Current CLR (all 3,933 KOs) | C eukaryote / prokaryote | A prokaryote KOs per GE |
|---|---:|---:|---:|
| Baseline (explained) | 15.95 (11.11) | 52.52 | 10.47 (2.55) |
| Days since shelf | −1.17 (−1.05) | −1.71 | −1.29 (−1.33) |
| SST change, 30 d | −0.75 (−0.62) | −1.39 | −0.71 (−1.02) |
| SST SD, 60 d | −0.84 (−0.76) | −1.30 | −0.69 (−0.78) |
| Chlorophyll mean, 60 d | −2.18 (−2.78) | −2.36 | −2.43 (−2.27) |
| Chlorophyll change, 30 d | +0.07 (+0.08) | −0.20 | −0.13 (−0.10) |
| Past chlorophyll peak | +0.04 (+0.05) | −0.14 | −0.18 (−0.25) |
| **MLD mean, 30 d** | −0.73 (−0.12) | **+6.04** | −2.00 (−0.94) |
| MLD change, 30 d | −0.74 (−0.68) | +0.57 | −0.94 (−0.71) |
| Mixed-layer light, day 0 | −2.46 (−2.19) | −6.68 | −0.95 (−0.82) |
| **Mixed-layer light, 30 d** | −0.10 (−0.62) | **+5.22** | −1.28 (−1.10) |
| Displacement, 60 d | −0.29 (−0.17) | −1.04 | −0.14 (−0.24) |
| Latitude of origin shift | −0.38 (−0.18) | −1.60 | −0.36 (−0.42) |
| Control (shuffled) | +0.05 (+0.05) | −0.38 | −0.30 (−0.35) |
| All features | −2.48 (−2.19) | −0.05 | −3.02 (−2.71) |
| Granite PCs only + all features (vs Granite alone: 17.24 / 54.94 / 10.98) | −2.05 (−2.26) | +0.13 | −2.38 (−2.69) |

**MLD on C, checked further.**
- **The history matters:** with day-0 MLD (GLORYS mlotst at sampling) C scores 52.26% (−0.26). The 30-day
  mean gives 58.56%. Spearman between the two at tile level is 0.66.
- **Combinations:** Granite + 30-day MLD 57.06% (Granite alone 54.94%); OLCI + SST + 30-day MLD **60.30%**
  (alone 56.20%), the best C model so far.
- **Correlations:** 30-day MLD vs C ρ +0.23 (all tiles); vs in-situ temperature ρ −0.48.
- **Per held-out cruise** (relative MSE reduction, baseline → + 30-day MLD; 9 of 13 improve):

  | Cruise | Tiles | Reduction |
  |---|---:|---:|
  | I05-2023 | 141 | 37.0% |
  | GOMECC4 | 53 | 17.7% |
  | HLY2202-2022 | 32 | 14.2% |
  | P18 | 64 | 14.0% |
  | AMT28 | 46 | 13.8% |
  | I07N | 146 | 6.9% |
  | C13.5/A13.5 | 216 | 4.7% |
  | P02-2022 | 89 | 4.6% |
  | I08S-2024 | 23 | 0.4% |
  | I09N-2025 | 84 | −4.4% |
  | A16N-2023 | 174 | −4.7% |
  | I09N | 54 | −4.8% |
  | A13.5-2024 | 93 | −5.1% |

**Test 2: residual axes** (Spearman within cruise, tiles; the shuffled control reaches |ρ| ≤ 0.055):

| Feature | Axis 1 (W1) | Axis 2 (W2) | C residual |
|---|---:|---:|---:|
| Days since shelf | +0.024 | −0.014 | −0.006 |
| SST change, 30 d | +0.102 (p 3e-4) | −0.057 | −0.032 |
| SST SD, 60 d | +0.103 (p 3e-4) | −0.086 (p 0.003) | −0.076 (p 0.008) |
| Chlorophyll mean, 60 d | −0.015 | −0.063 | −0.070 |
| Chlorophyll change, 30 d | −0.122 (p 2e-5) | −0.006 | +0.072 |
| Past chlorophyll peak | +0.072 | −0.002 | −0.055 |
| MLD mean, 30 d | −0.092 (p 0.001) | −0.028 | +0.067 |
| MLD change, 30 d | −0.033 | +0.040 | −0.013 |
| Mixed-layer light, day 0 | +0.069 | +0.000 | +0.018 |
| Mixed-layer light, 30 d | +0.095 (p 9e-4) | +0.073 | −0.047 |
| Displacement, 60 d | −0.124 (p 1e-5) | +0.020 | +0.051 |
| Latitude of origin shift | −0.041 | +0.095 (p 9e-4) | −0.021 |
| Control (shuffled) | +0.024 | +0.016 | −0.055 |

**Reading.**
1. **Water-mass history does not explain the KO residuals.**
   - No feature adds to the current CLR or to A. The best gain, +0.07 points for chlorophyll change, is the
     size of the shuffled control (+0.05).
   - Most features cost 0.3–2.5 points: an extra spline overfits under cruise hold-out.
   - Within cruise, every link to axis 1 or axis 2 is |ρ| ≤ 0.12, i.e. ≤ 1.5% of variance:
     - Axis 1 is slightly higher in water that warmed, varied more in SST, lost chlorophyll and travelled
       less.
     - Axis 2 is slightly higher in water from the pole side.
   - Axis 2's driver is not water-mass origin either. The "remaining within cruise" slice stays unknown.
2. **The exception is the eukaryotic amount (C).**
   - The 30-day mean mixed-layer depth adds 6 points (52.5 → 58.6%), and the 30-day mixed-layer light adds 5.
   - It is the history: today's MLD adds nothing.
   - The gain holds for 9 of 13 held-out cruises, largest for I05-2023 (37% lower error), GOMECC4, HLY2202,
     P18 and AMT28.
   - Deeper recent mixing goes with more eukaryotes relative to prokaryotes at fixed temperature. That is
     consistent with mixing favouring larger eukaryotic phytoplankton.
   - Best C model: OLCI + SST + 30-day MLD, 60.3%.
3. **Caveats.**
   - MLD is a GLORYS ocean-model field, not a satellite observation. It is still available anywhere without
     a ship.
   - Mixed-layer light at day 0 hurts (−6.7): KD490 × MLD at a single day is noisy at the extremes.
4. **For the two-stage model (R51):** use the improved C predictor (OLCI + SST + 30-day MLD) as the covariate
   for the whole-community and shared-KO models.

## R54. Two-stage model done properly; technical variance out of the denominator; R51's 19.4% retracted (2026-11-21)

**Framing (decided 2026-11-21).**
- **Story:** environment selects who is there, and who is there determines which genes are present. In
  model terms: satellite (+ mixing history) → community structure (domain balance, C) → gene content
  (A, B, D).
- **Reporting:**
  - Technical variance (platform, sequencing depth) is removed from the denominator, not added as model
    covariates. Fit on training folds only, subtract from the targets, and report the share of the
    remaining, biological variance that the satellite explains.
  - A scree-style plot of what makes up 100% of the variance comes later, followed by headline numbers in
    that form.

**Script** `scripts/two_stage_models.py --parts spline forms headline` →
`data/results/two_stage/{scores.csv, stage1_C.csv, summary.txt, stage2_forms.csv, plugin_headline.csv}`;
figure `figures/output/two_stage.png`.
- **Stage 1** predicts C from one of: `base` (32 Granite PCs + temperature spline), `sat` (16 OLCI band
  means + satellite SST, splines) or `sat_mld` (sat + 30-day mean MLD along the back-trajectory, R53).
- **Stage 2:** 32 Granite PCs (or + temperature spline) + a C term.
- Everything is held out by cruise, ridge with per-target GCV alpha, on the R43 tiles.
- **Technical mode:** per fold, on training tiles only, OLS per target on platform (NovaSeq 6000 vs X cruise
  groups) + a cubic spline (8 knots) of log10 KO reads. The fitted part is subtracted from training and
  test targets (test with the training coefficients), and the model is fitted and scored on the adjusted
  targets.

**1. Stage-1 accuracy for C** (outer, cruise held out): base 52.52%, sat 56.20%, sat_mld 60.30%.

**2. R51's two-stage gain was not real.**
- R51 used, for the training tiles, C predictions from leave-one-cruise-out models over all cruises,
  including the held-out cruise. Reproduced: 19.44%.
- **Nested properly** (training-tile C from inner leave-one-cruise-out within the training cruises), with C
  as a spline in stage 2:

  | Stage 2 | C from base | C from sat | C from sat_mld | No C |
  |---|---:|---:|---:|---:|
  | All KOs, baseline (+temp) | 16.42 | 14.94 | 11.63 | 15.95 |
  | All KOs, Granite PCs | 16.74 | 16.06 | 14.32 | 17.24 |
  | D, Granite PCs | 12.74 | 11.60 | 9.72 | 12.88 |
  | A, Granite PCs | 10.81 | 9.39 | 8.07 | 10.98 |

  A better stage-1 predictor makes stage 2 worse.
- **Placebo** (R51 scheme, all KOs): with the held-out cruise's C scrambled inside stage 1, skill drops to
  0.44%; with it set to a constant, 1.64%. The stage-2 spline extrapolates badly when the training-tile C
  values are distributed differently from the test-tile values.
- **Conclusion:** the R51 two-stage figure (16.0 → 19.4%) is retracted. The R51 reading on C, A, B and D
  stands.

**3. A robust stage-2 form: plug-in, linear.** Stage 2 is fitted on the observed C of the training tiles,
which is known there (a measurement, not a prediction). The held-out cruise gets its stage-1 prediction,
from a model fitted on the training cruises only. No inner loop and no test-cruise information.

| Target | Stage 1 | Nested, linear | Plug-in, spline | Plug-in, linear |
|---|---|---:|---:|---:|
| All KOs, stage 2 baseline | base | 15.80 | 15.29 | 15.86 |
| All KOs, stage 2 Granite | base | 16.74 | 15.16 | 16.39 |
| All KOs, stage 2 baseline | sat_mld | 14.85 | 17.30 | 17.72 |
| All KOs, stage 2 Granite | sat_mld | 16.08 | 17.22 | **18.27** |
| D, stage 2 baseline | base | 11.13 | 10.72 | 11.20 |
| D, stage 2 Granite | base | 12.75 | 10.67 | 12.06 |
| D, stage 2 baseline | sat_mld | 10.70 | 12.12 | 12.47 |
| D, stage 2 Granite | sat_mld | 11.89 | 12.24 | **13.40** |
| A, stage 2 baseline | base | 10.54 | 10.24 | 10.42 |
| A, stage 2 Granite | base | 10.66 | 9.66 | 10.71 |
| A, stage 2 baseline | sat_mld | 10.53 | 10.30 | 10.91 |
| A, stage 2 Granite | sat_mld | 10.40 | 9.92 | **11.16** |

**4. Headline** (stage 2 = Granite PCs + linear C; stage 1 = OLCI + SST + 30-day MLD; raw units, KOs
standardised in brackets):

| Target | Baseline (Granite + temp) | Granite PCs | Two-stage | Baseline, tech removed | Granite, tech removed | Two-stage, tech removed | Adjusted denominator / total |
|---|---:|---:|---:|---:|---:|---:|---:|
| All KOs (current CLR) | 15.95 (11.11) | 17.24 (12.66) | **18.27 (13.74)** | 17.79 (13.84) | 18.44 (14.71) | **19.15 (15.93)** | 1.018 |
| C eukaryote / prokaryote | 52.52 | 54.94 | 60.30 (stage-1 model itself) | 55.94 | 55.25 | 55.38 | 0.978 |
| D shared CLR | 11.38 (9.81) | 12.88 (11.30) | 13.40 (12.15) | 13.43 (12.56) | 14.65 (13.54) | 15.26 (14.77) | 1.033 |
| A prokaryote-specific per GE | 10.47 (2.55) | 10.98 (3.49) | 11.16 (4.03) | 13.44 (5.68) | 13.78 (6.76) | 14.01 (7.21) | 0.987 |
| Dp shared prokaryote part per GE | 13.23 (3.91) | 13.77 (4.89) | 13.92 (5.62) | 16.21 (6.86) | 16.69 (7.99) | 16.92 (8.46) | 0.993 |
| B eukaryote-specific within eukaryotes | 8.74 (3.76) | 10.94 (6.45) | 11.55 (6.77) | 9.13 (3.93) | 11.05 (6.39) | 11.46 (6.63) | 1.013 |
| De shared eukaryote part within eukaryotes | 23.34 (14.32) | 23.68 (15.32) | 24.38 (15.78) | 21.09 (13.32) | 21.05 (14.00) | 21.23 (14.11) | 1.014 |

**Reading.**
1. **The two-stage gain is real but modest.**
   - All KOs: 17.24 → 18.27% (+1.0 point over Granite PCs alone, +2.3 over the R43 baseline).
   - D +0.5, A +0.2, Dp +0.2, B +0.6, De +0.7.
   - It works only in the plug-in linear form with the best stage-1 C (OLCI + SST + 30-day MLD).
   - Story check: the domain balance carries some information the Granite features alone do not. But most
     of what the satellite says about gene content already goes through the same channel; C is largely a
     re-expression of it.
2. **Granite PCs alone beat the R43 baseline everywhere.** The in-situ temperature spline overfits under
   cruise hold-out. Stage 2 should be Granite only, which is also cleaner for the story (no ship
   measurement).
3. **Technical variance out, fold-wise:**
   - All KOs: two-stage 19.15%, Granite 18.44%, baseline 17.79%. A gains most (+2.9), C barely moves.
   - **The adjusted denominator does not shrink** (all KOs 1.018 × the original total). Platform is
     cruise-level (NovaSeq 6000 vs X groups), so for a held-out cruise the platform offset is extrapolated
     from other cruises. Subtracting it adds about as much variance as it removes.
   - The gains in this mode come from cleaner training targets, not from a smaller denominator.
   - The in-sample R43 budget puts the technical share at 6.1% (platform 3.8 + depth 2.3). Dividing by
     that instead (explained / 0.939): baseline 16.99%, Granite 18.36%, two-stage 19.46%.
   - Relative to the reproducible variance (1 − replicate noise 0.270): 21.8 / 23.6 / 25.0%.
   - Relative to the reproducible, non-technical variance (1 − 0.270 − 0.061 = 0.669): 23.8 / 25.8 /
     **27.3%**.
4. **Caveat on C and MLD.** With technical variance removed, the stage-1 model OLCI + SST + 30-day MLD drops
   from 60.3 to 55.4%, while the baseline rises 52.5 → 55.9%. Part of the 30-day MLD gain on C (R53)
   coincides with the platform grouping, since the NovaSeq X cruises are 2023–2025 and include I05-2023,
   the biggest MLD gainer. To check: the MLD gain on C within platform groups only.
5. **For the scree / 100% figure:** use the in-sample budget partition (R43 style) per target. The headline
   is then explained / (total − technical) and explained / (total − technical − replicate noise),
   reported side by side with the fold-wise number.

## R55. External read download: status after job 88813; the filter-size test has to use EMO BON, not Tara (2026-11-22)

**Job 88813** (`cluster/download_reads.sbatch`, all three sets, 3 in parallel) hit the 48 h limit on
2026-11-21. The script is resumable; resubmit the same command.
- ENA errors throughout: curl 52 (empty reply), curl 18 (transfer closed part-way), 403 throttling. One
  full-size file failed MD5 and restarted (ERR12071056_1).
- One file given up on in that run: ERR11733532_2 (Tara); the next submission retries it.

| Set | Files verified (.ok) | Files total | On disk (incl. partial) | Set size |
|---|---:|---:|---:|---:|
| Tara Pacific | 43 | 168 | 593 GB | 1.74 TB |
| JGI (Flinkstrom) | 5 | 180 | 78 GB | 0.53 TB |
| EMO BON | 0 | 262 | 23 GB | 1.19 TB |
| All | 48 (`download_reads_done.txt`) | 610 | about 0.69 TB | 3.46 TB |

About 0.7 TB in 48 h, roughly 14 GB/h, so the remaining 2.8 TB needs about 4 more 48 h jobs at this rate.

**The filter-size test cannot use Tara Pacific.** All 387 Tara Pacific samples in `samples.csv` (84 runs)
are the 0.2–3 µm fraction. **EMO BON has both:**
- 131 runs: 0.2–3 µm 66, 3–200 µm 63, unknown 2;
- **33 site-dates with both fractions**, a paired design for "does the filter size shift the eukaryote
  fraction (C)?" (R46 / R54 note).
- Median 130 M reads per run, against 267 M for Tara.

**Priority:** EMO BON first for the filter-size test; Tara and JGI for the external transfer tests.

**Option:** only the first N read pairs of each run (for example 20 M pairs, about 7–15% of the bytes) would
cut the download time about tenfold. Profiles are rarefied to 250,000 KO reads anyway. Files would then be
size-checked only, with no MD5 check. Not implemented; the user decides.

**R55 addendum: subsample mode (2026-11-22).** `cluster/download_reads.sbatch` takes `SUBSAMPLE_PAIRS=N`:
- Each file is streamed, cut after the first 4N lines (N reads, so N pairs per run) and re-compressed.
- Checks replace the ENA MD5, which cannot apply to a cut file: exactly 4N lines; for a run with fewer reads, the whole
  stream must end cleanly (curl exit 0 and gzip exit 0, so gzip's CRC and length pass).
- Files already complete on disk are cut locally and the full originals moved to `reads/fastq_full/`; earlier `.part` files
  are reused if long enough. Done files get `.ok` plus `.sub` (holds N), so `map_reads.sbatch --list` sees them.
- Tested locally on synthetic FASTQ (cut, short-whole, local-cut and truncated-stream cases).

Choice of N = 20 M pairs:
- Bio-GO-SHIP libraries (`library_metadata.csv` total_reads, 3,174 libraries): 10% 5.7 M, median 14.4 M, 90% 40.6 M read pairs.
  20 M is inside the training range, so external runs are also depth-matched, which matters for the coverage-50 threshold (R47).
- Bytes needed at 20 M pairs (from runs.csv read_count and fastq_bytes):

| Set | Runs | Full | First 20 M pairs | Runs with < 20 M pairs |
|---|---:|---:|---:|---:|
| EMO BON | 131 | 1.19 TB | 0.180 TB | 0 |
| Tara Pacific | 84 | 1.74 TB | 0.122 TB | 0 |
| JGI | 90 | 0.53 TB | 0.325 TB | 15 |
| Total | 305 | 3.46 TB | 0.627 TB | 15 |

- Caveat: the first reads of a file come from the first flowcell tiles; quality there can be slightly lower, but trimming
  handles it and taxon/KO composition does not depend on tile.

## R56. Within-platform check of the MLD gain on C; the 100% variance partition (2026-11-22)

Script: `scripts/mld_platform_partition.py` (`--parts mld halves partition figure`). Outputs: `data/results/mld_platform/`
(`mld_check.csv`, `mld_per_cruise.csv`, `mld_random_halves.csv`, `partition.csv`, `summary_*.txt`);
`figures/output/variance_partition.png`. All models: R43 tiles, one cruise held out, ridge with per-target GCV alpha.
Platform groups as R43: NovaSeq 6000 = 6 cruises, 558 tiles; other = 7 cruises, 657 tiles.

### 1. Is the 30-day MLD gain on C a platform proxy?

Question: R54 found that removing technical variance lowers C (raw_sat + MLD30) from 60.3 to 55.4%. Does MLD30
stand in for the platform difference?

**a/b. Variance of C explained (%), cruise held out, scored on each group's tiles:**

| Training | Model | All tiles | NovaSeq 6000 cruises | Other cruises |
|---|---|---:|---:|---:|
| Pooled | raw_sat | 56.20 | 45.79 | 53.27 |
| Pooled | raw_sat + MLD30 | 60.30 | 49.44 | 58.44 |
| Pooled | raw_sat + platform | 61.58 | 55.24 | 57.45 |
| Pooled | raw_sat + MLD30 + platform | 64.95 | 57.00 | 62.40 |
| Within platform group | raw_sat | 61.04 | 53.36 | 57.56 |
| Within platform group | raw_sat + MLD30 | 61.57 | 51.85 | 59.33 |

MLD30 gain (percentage points):

| Comparison | All | NovaSeq 6000 | Other |
|---|---:|---:|---:|
| Pooled, raw_sat → + MLD30 | +4.10 | +3.65 | +5.17 |
| Pooled, with platform covariate | +3.37 | +1.77 | +4.95 |
| Within-platform training | +0.53 | −1.51 | +1.77 |

**Control: random cruise splits.** Within-group training on random 6 / 7 cruise splits (20 splits, seed 0) instead of the
platform split:

| Statistic | MLD30 gain |
|---|---|
| Median | +0.77 pp |
| 10–90% | −1.59 to +8.50 pp |
| Splits with gain > 0 | 13 of 20 |
| Platform split | +0.53 pp, typical of the random splits |

So the small within-platform gain comes from halving the number of training cruises, not from MLD30 acting as a platform
proxy. MLD30 is a between-cruise feature: 58.7% of its tile variance lies between cruises.

**c. MLD30 by platform group** (log10 m):

| Group | Mean | SD |
|---|---:|---:|
| NovaSeq 6000 | 1.321 | 0.191 |
| Other | 1.399 | 0.273 |

- Point-biserial r −0.161.
- Platform group explains 2.6% of MLD30 variance.

**d. Within-cruise Spearman of MLD30** (both centred per cruise):

| Tiles | With raw_sat C residual | With C itself |
|---|---|---|
| All | +0.029 (p 0.31) | +0.149 (p 1.8e-7) |
| NovaSeq 6000 | −0.162 (p 1.2e-4) | −0.060 (p 0.16) |
| Other | +0.173 (p 8.7e-6) | +0.297 (p 6.9e-15) |

**e. Per cruise, mean squared error of C** (raw_sat → + MLD30, pooled training; bias b = observed − predicted):

| Cruise | NovaSeq 6000 | Tiles | MLD30 | Mean C | b raw_sat | b + MLD30 | MSE raw_sat | MSE + MLD30 | Change, pooled | Change, within-platform |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| I05-2023 | no | 141 | 1.793 | −1.784 | 0.599 | 0.142 | 1.143 | 0.757 | −33.8% | −12.6% |
| GOMECC4 | no | 53 | 1.187 | −4.002 | −0.008 | −0.098 | 0.321 | 0.280 | −12.9% | −7.7% |
| P02-2022 | no | 89 | 1.151 | −3.349 | 0.237 | 0.291 | 0.531 | 0.536 | +1.0% | −1.6% |
| A16N-2023 | no | 174 | 1.386 | −3.113 | 0.139 | 0.248 | 0.559 | 0.568 | +1.6% | −5.1% |
| I09N-2025 | no | 84 | 1.266 | −3.486 | 0.475 | 0.522 | 0.526 | 0.540 | +2.7% | +3.9% |
| A13.5-2024 | no | 93 | 1.281 | −3.105 | 0.098 | 0.191 | 0.572 | 0.629 | +10.0% | +0.9% |
| I08S-2024 | no | 23 | 1.492 | −2.232 | 0.106 | 0.218 | 0.202 | 0.255 | +26.2% | +47.7% |
| HLY2202-2022 | yes | 32 | 1.067 | −1.331 | −0.572 | −0.348 | 2.505 | 2.155 | −14.0% | −6.2% |
| C13.5/A13.5 | yes | 216 | 1.415 | −4.025 | −0.350 | −0.346 | 0.359 | 0.317 | −11.6% | −4.8% |
| P18 | yes | 64 | 1.278 | −3.177 | −0.121 | −0.085 | 0.569 | 0.507 | −10.9% | −10.5% |
| I07N | yes | 146 | 1.307 | −3.811 | −0.072 | −0.034 | 0.219 | 0.195 | −10.7% | −4.4% |
| AMT28 | yes | 46 | 1.419 | −3.644 | −0.443 | −0.444 | 0.529 | 0.515 | −2.6% | −0.8% |
| I09N | yes | 54 | 1.100 | −4.406 | −0.543 | −0.670 | 0.472 | 0.651 | +37.8% | +72.6% |

- MLD30 lowers error in 5 of 6 NovaSeq 6000 cruises (both schemes).
- It lowers error in 2 of 7 other cruises (pooled) and 4 of 7 (within-platform).
- Every NovaSeq 6000 cruise has a negative raw_sat bias (−0.07 to −0.57): C is lower on that platform than the satellite
  predicts. The platform covariate alone adds +5.38 pp, and the within-platform raw_sat score (61.0%) is higher than all
  20 random splits (maximum 59.2%).

**Conclusion.**
- The MLD30 gain is not a platform proxy:
  - MLD30 barely differs between platforms;
  - the gain survives a platform covariate (+3.37 pp);
  - it helps 5 of 6 cruises on one platform;
  - the platform-split loss matches random cruise splits.
- The gain is a between-cruise effect. It is carried most by I05-2023, the deepest-mixing cruise (MLD30 1.79, highest C).
  Its error drop (tiles × MSE change, −54.4) is 83% of the pooled total (−65.6).
- Separately, C carries a real platform offset (NovaSeq 6000 lower). That is why technical removal lowers C (R54).

### 2. The 100% variance partition

Each target's total tile variance is split into five parts:
- **Technical:** in-sample OLS on platform + cubic spline (8 knots) of log10 KO reads, all tiles.
- **Explained:** cruise held out, on the technically adjusted targets.
  - KO targets: the R54 plug-in linear two-stage (32 Granite PCs + standardised C; observed C in training, C predicted from
    raw_sat + MLD30 for the test cruise).
  - C: raw_sat + MLD30.
- **Cruise offsets:** per-cruise means of the held-out residual.
- **Replicate noise:** within-tile replicate variance / replicates per tile (R51 budget), capped at what is left.
- **Within-cruise:** the rest.

Parts sum to 99.86–100.00%; the shortfall is the overall mean bias of the held-out predictions.

| Target | Technical | Explained | Cruise offsets | Replicate noise | Within-cruise | Explained / non-technical | Explained / reproducible non-technical | Strict fold-wise / non-technical (R54) | Granite alone / non-technical |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| All KOs (CLR) | 9.09 | 15.46 | 12.32 | 27.04 | 35.99 | **17.01** | **24.21** | 19.15 | 16.96 |
| C | 15.95 | 43.55 | 7.32 | 6.49 | 26.68 | **51.81** | **56.15** | — (plain 60.30) | — |
| D | 7.82 | 12.35 | 12.35 | 27.00 | 40.37 | 13.40 | 18.95 | 15.26 | 13.28 |
| A | 11.69 | 11.25 | 13.86 | 25.92 | 37.14 | 12.74 | 18.04 | 14.01 | 12.78 |
| Dp | 11.27 | 13.83 | 12.85 | 23.47 | 38.45 | 15.59 | 21.20 | 16.92 | 15.63 |
| B | 6.99 | 8.87 | 9.29 | 57.60 | 17.17 | 9.54 | 25.06 | 11.46 | 9.56 |
| De | 8.56 | 16.55 | 10.70 | 45.06 | 19.02 | 18.10 | 35.67 | 21.24 | 18.44 |

**Notes.**
- **Supersedes R54's quick 6.1% technical share and its 19.46% / 27.3%.** Those took platform + depth from the
  *residual* budget, so they only counted technical variance the model had not already absorbed. Fitted on the targets
  themselves, the technical share for all KOs is 9.09%. The corrected headlines are 17.0% of non-technical and 24.2% of
  reproducible non-technical variance.
- **Two technical-removal numbers.**
  - In-sample removal (this table) is the conservative one: removing platform on all tiles also takes out whatever part
    of the cruise-level platform contrast the satellite could predict (platform is confounded with year and region).
  - The strict fold-wise number (R54) is higher for every target (19.15% for all KOs).
- **In-sample, the two-stage model adds ≈ 0 over Granite alone** (17.01 vs 16.96 for all KOs; 12.74 vs 12.78 for A).
  R54's +0.7 pp gain (fold-wise) does not survive in-sample technical removal. With the platform offset of C taken out
  of the targets, the C term has little left to add.
- **Replicate noise is the largest single part for B (57.6%) and De (45.1%).** Against reproducible non-technical
  variance, the eukaryotic targets are as predictable as the rest (B 25.1%, De 35.7%).
- **C is the only target where explained is the largest part** (43.6% of total).
- **The largest unexplained biological part for the KO targets is within-cruise** (36–40%), not cruise offsets (12–14%).

## R57. Metadata status per cruise, and the open questions for the Bio-GO-SHIP team (2026-11-23)

Script: `scripts/metadata_status.py`. Outputs in `data/results/metadata_status/`:
- `cruise_metadata_status.csv` (one row per cruise)
- `faire_fill_rate.csv` (share of samples with each FAIRe field filled)
- `biogoship_metadata_check.xlsx`, to send to the team, with two sheets:
  - **Questions**: the six asks below, with an answer column.
  - **Cruise methods**: our values. Yellow = unknown ("?"), orange = our two sources disagree. Columns for
    "All correct?" and corrections.

The table replaces the stale sequencer column of `transfer/results/protocol/cruise_methods.csv`, which still says
"NovaSeq X? (unverified)". All three 2023–2025 flowcells are NovaSeq X, confirmed by original FASTQ headers (R40, R41).

**Clarifications recorded.**
- **The `seq_run_id` field** comes from the FAIRe workbook (`FAIRe-NOAA_uci-ess-biogoship.xlsx`, sheet
  `experimentRunMetadata`). It is `date_flowcell_MetaG_PE150`.
- **NovaSeq 6000 flowcells** (…DSXX / DSXY / DRXX / DSX5 / DSX7): ENA `instrument_model` for 906 runs. HLY2202 comes from
  the publication plus the flowcell format.
- **NovaSeq X flowcells** (…LT4, e.g. 22YJMLLT4): instrument IDs LH00626, LH00316 and LH00162 in the original headers.
- **"NovaSeq X" is Illumina's product name, not a placeholder.** NovaSeq X / X Plus were released in 2023 and use
  XLEAP-SBS chemistry.

**New conflict found.** For A16N-2023, I05-2023 and A13.5-2024, the methods table (2026-11-19) says 1 mL DNA/RNA Shield,
while the FAIRe workbook's `samp_store_sol` says "lysis buffer".

**What is missing, and why it matters for this project:**

| # | Ask | Cruises | Why |
|---:|---|---|---|
| 1 | Library prep kit, input DNA, PCR cycles, size selection | all 2023–2025 (NovaSeq X) runs | prep is confounded with the sequencer change; PCR cycles drive GC bias |
| 2 | Extraction kit and lysis | P02-2022, A16N, I05, A13.5, I08S, I09N-2025 | bead beating changes eukaryote / Gram-positive capture |
| 3 | Preservation: Shield or lysis buffer; HLY2202 buffer and storage | A16N, I05, A13.5, HLY2202 | candidate cause of the older vs newer cruise shift |
| 4 | Pre-filtration; I08S pore size | HLY2202, P02 (CTD), A16N, I05, A13.5, I08S, I09N-2025 | changes the eukaryote : prokaryote ratio (C) directly |
| 5 | Same DNA on NovaSeq 6000 and X, or the raw HWKVFDSX5 files (P02 pairs, 58 libraries) | P02-2022 | the only direct test of the sequencer effect |
| 6 | Per-sample volume, extraction plate, DNA yield | 2016–2022 cruises, P02, I08S, I09N-2025 | within-cruise batch tests (lower priority) |

**Per-sample FAIRe fill** (share of samples):
- Volume: A13.5, A16N, GOMECC4, HLY2202 and I05 only.
- Extraction plate: A13.5, A16N, GOMECC4 and I05.
- DNA yield: GOMECC4 only.
- Collection method: every cruise but GOMECC4 and HLY2202.

**R57 update (2026-11-23): asks 1 and 5 dropped.**
- Dropped: the 2023–2025 library prep (ask 1) and same DNA on both sequencers (ask 5). The platform term already captures
  that prediction shift.
- The paper will state that the shift can't be attributed to the sequencer rather than the library prep or workflow
  change that came with it.
- The sheet now has four asks: extraction (2022–2025), the preservation conflict, pre-filtration, and per-sample records.
- What is tested vs still open on the non-biological side: see the chat summary of 2026-11-23, recorded here:
  - **Tested, no effect on the residual axes:**
    - library prep (Nextera vs DNA Flex, HLY2202)
    - preservation + bead beating within Protocol B
    - same-DNA re-sequencing
    - flowcell
    - storage time
    - volume
    - collection method (underway vs CTD)
    - extraction plate (the 4 cruises that record it)
    - sampling year within a workflow
    - gene GC content (R44)
  - **Removed:** platform and sequencing depth (R54, R56). Protocol A samples were excluded (R41).
  - **Open:**
    - pre-filtration (unknown for 7 cruises; acts on C directly)
    - extraction/preservation differences *within* the NovaSeq X era (GOMECC4: no preservative, 2 L Niskin, MagBead)
    - extraction plate and lane for the 2016–2022 cruises (not recorded)
    - JGI's HiSeq 2500 and EMO BON's filter fractions, for the external tests

**R57 note (2026-11-23): what "replicate noise" in the R51/R56 budgets is.**
- **Formula.** σ² = mean over pairs of (difference between two samples in the same satellite tile)² / 2, per KO. The
  budget term is Σ_tiles σ² / k_tile, where k is the number of samples per tile.
- **What it measures.** The target is the tile mean, and 1,082 of 1,215 tiles have only one sample. So the term is
  essentially the single-sample scatter, as a share of the total variance across tiles. It is not the noise of averaged
  replicates.
- **The pairs** (C target tiles): 1,416 samples; 133 tiles hold ≥ 2 (84 with 2, 42 with 3, 7 with 5–6); 305 pairs.
  - 90% of pairs are the same place and time (< 0.5 km, < 1 h). The median depth difference is 0, but 25% differ by
    ≥ 4 m (maximum 7 m).
  - Pairs by cruise: GOMECC4 222 (A/B/C replicates), I07N 27, C13.5 22, A13.5 12, HLY2202 7, A16N 4, I05 4, P02 3,
    P18 2, I09N-2025 2.
- **So it mixes field-replicate differences** (filter, extraction, library, sequencing, read sampling) with a little
  small-scale and depth biology.
- **Not driven by GOMECC4 alone.** Single-sample replicate variance as a share of total tile variance:

| Target | From GOMECC4 pairs | From the other cruises' pairs |
|---|---:|---:|
| clr_all | 27.3% (222 pairs) | 33.0% (83) |
| C | 6.5% (222) | 8.1% (83) |
| B | 63.8% (204) | 55.6% (80) |

## R58. FAIRe `instrument` field filled; "replicate noise" renamed "within-tile variability" (2026-11-23)

**FAIRe instrument.** Script: `scripts/faire_instrument_fill.py`. The original workbook is not modified. Outputs in
`data/results/metadata_status/`:
- `FAIRe-NOAA_uci-ess-biogoship_instrument.xlsx`: a copy, only `experimentRunMetadata.instrument` changed, changed
  cells yellow
- `instrument_column.txt`: the column for pasting, 4,993 lines in sheet order
- `instrument_changes.csv`: per change, the old value, new value and evidence

**OBI terms** (OLS4 / Ontobee, checked 2026-11-23):

| Instrument | OBI term |
|---|---|
| Illumina NovaSeq 6000 | OBI_0002630 |
| Illumina NovaSeq X | OBI_0003663 |
| Illumina NovaSeq X Plus | OBI_0003664 |
| Illumina HiSeq 4000 | OBI_0002049 |
| Illumina HiSeq 2500 | OBI_0002002 |

- The workbook's values were placeholders: "INSTRUMENT NAME [OBI_2131 / 2631 / 3310]".
  - Those IDs are not sequencers: OBI_0002131 is a methylation BeadChip, OBI_0002631 an scRNA-seq assay and
    OBI_0003310 a SPRITE assay.
  - But they map one-to-one onto our flowcell assignment (2131 → HiSeq 4000, 2631 → NovaSeq 6000, 3310 → NovaSeq X), so
    the team's grouping was right.

**Filled:** 3,186 MetaG rows.

| Value | Rows | Evidence |
|---|---:|---|
| Illumina HiSeq 4000 [OBI_0002049] | 55 | ENA |
| Illumina NovaSeq 6000 [OBI_0002630] | 906 | ENA |
| Illumina NovaSeq 6000 [OBI_0002630] | 193 | publication + flowcell (HFJLFDSX7, HMY7WDSX7) |
| Illumina NovaSeq 6000 [OBI_0002630] | 96 | flowcell format only (HWKVFDSX5, unverified) |
| Illumina NovaSeq X [OBI_0003663] | 1,936 | FASTQ headers |

- NovaSeq X vs X Plus is unknown: both use "LH" instrument IDs.
- Left as is: 1,807 amplicon and MetaT rows, whose `seq_run_id` has no flowcell.

**Naming.** The R51/R56 "replicate noise" term is now called **within-tile variability**: how much samples that share a
satellite tile (the same model input) differ, as a share of the total variance across tiles. It is the ceiling on
satellite skill, not purely measurement noise. The `partition.csv` column keeps the name `replicate_noise`, and the
figure is relabelled.

**Unit of modelling.** The model already uses the tile as the unit, and has since R17: `tune_figure_kos.load` averages
all samples sharing a tile (same Granite input) into one target row (1,416 samples → 1,215 tiles; 133 tiles hold ≥ 2
samples). The budget term already scales with 1/k per tile.

## R59. Main model for the paper: Granite PCs alone; C as its own result; within-tile variability by k (2026-11-23)

**Decision (user).** The chained / two-stage models (R51, R54, R56) are dropped as the headline: more complex and
not clearly better. The main model is **32 Granite PCs, ridge, one cruise held out**. It uses no ship data. C is
reported as a separate result. "Environment selects who is there" becomes two parallel findings: C (domain balance)
is highly predictable, while gene content within domains is weakly predictable.

**Setup** (`scripts/main_model.py`, outputs in `data/results/main_model/`)
- Tiles: R43 set, 1,215 (1,169 for B and De). Samples sharing a tile are averaged.
- Technical removal: in-sample, as in R56. Per target column, OLS on platform + spline of log10 KO reads, fitted on all
  tiles and removed. Models are fitted and scored on the adjusted targets.
- Ridge penalty three ways: GCV per column (as before), leave-cruise-out (inner, training cruises only) per column
  (`cv`), and leave-cruise-out with one penalty per target set (`cv_shared`).

**Explained, % of non-technical variance (pooled, cruise held out)**

| Target | Granite GCV | Granite cv | Granite cv_shared | Granite + temp GCV | Granite + temp cv | Granite + temp cv_shared | Temp alone |
|---|---:|---:|---:|---:|---:|---:|---:|
| All KOs (CLR) | 16.96 | **17.60** | 17.45 | 15.90 | 17.51 | 17.61 | 10.37 |
| D | 13.28 | **13.96** | 13.83 | 11.53 | 13.71 | 14.03 | 4.27 |
| A | 12.78 | **13.71** | 13.33 | 12.03 | 13.93 | 13.46 | 6.99 |
| Dp | 15.63 | **16.36** | 15.83 | 14.72 | 16.34 | 15.97 | 7.44 |
| B | 9.56 | **10.32** | 9.77 | 7.12 | 10.07 | 9.96 | 5.94 |
| De | 18.44 | **18.86** | 18.79 | 17.89 | 18.96 | 18.77 | 17.40 |

- **Leave-cruise-out penalty is a small free gain:** +0.4 to +0.9 pp over GCV for every target (All KOs +0.64).
  Median log10 alpha rises from 2.0 (GCV) to 2.25–2.75. The median per-column calibration slope (observed on
  predicted) moves toward 1 (All KOs 0.755 → 0.796, cv_shared 0.813; A 0.609 → 0.680; B 0.585 → 0.660), so R38's over-dispersion
  is partly a penalty-choice effect. Pooled slopes were already 0.94–0.99.
- Main KO model: Granite PCs, per-column leave-cruise-out penalty (`cv`).
- **Granite alone vs Granite + temperature:** with GCV, Granite alone wins clearly (16.96 vs 15.90; the
  temperature spline got tiny penalties, log10 alpha −0.75 to −1). With the leave-cruise-out penalty the two are
  equal within 0.3 pp (17.60 vs 17.51 / 17.61). Temperature adds nothing once the penalty is tuned properly, so the
  ship-free model is the main one.
- Per-KO median R² (All KOs): cv 0.1157, cv_shared 0.1275, GCV 0.1200, Granite + temp (cv) 0.1111, temperature alone
  0.0089. Per-column penalties overfit for low-variance KOs: they win on the pooled (variance-weighted) number and
  lose on the per-KO median. The difference is ≤ 1.2 pp either way.

**C (log euk / prok), % of non-technical variance, cruise held out**

| Model | GCV | cv |
|---|---:|---:|
| OLCI + SST (raw_sat) | 50.91 | 50.32 |
| OLCI + SST + 30-day MLD | 50.54 | 50.63 |
| Granite PCs | 53.88 | 53.80 |
| Granite PCs + temperature | 54.79 | 55.30 |

- **The 30-day MLD gain does not survive full technical removal.** Which technical term is removed matters (cv
  penalty):

| C adjusted for | OLCI + SST | + 30-day MLD | gain (pp) | Granite PCs | Granite + temp |
|---|---:|---:|---:|---:|---:|
| nothing (raw C) | 55.72 | 60.44 | +4.72 | 53.95 | 54.82 |
| platform only | 54.37 | 55.84 | +1.47 | 53.66 | 54.61 |
| read depth only | 53.86 | 56.57 | +2.71 | 52.15 | 53.67 |
| platform + depth | 50.32 | 50.63 | +0.31 | 53.80 | 55.30 |

  - MLD30 correlates with log10 KO reads (r −0.27 across tiles; −0.09 within cruises) and with platform (r −0.16).
  - Read depth is not predictable from the satellite (held-out R² −0.11 to −0.12 for OLCI + SST, + MLD and Granite),
    so treating it as technical is justified.
  - Consistent with R56: within-platform MLD gain on raw C was only +0.53 pp.
  - Conclusion: with technical variance removed, **Granite PCs predict C better than OLCI + SST + MLD30** (53.8% vs
    50.6%). The R53/R56 MLD signal on C was largely shared with platform and depth.
- C figure built with OLCI + SST + 30-day MLD as specified (`figures/output/c_prediction.png`): a observed vs
  held-out predicted (R² 0.51, 1,215 tiles), b per-cruise R² for the three inputs (* = NovaSeq 6000 cruise), c C
  against latitude (41-tile running means). Per-cruise R² (`c_per_cruise.csv`):

| Cruise | OLCI + SST | + MLD30 | Granite PCs |
|---|---:|---:|---:|
| I07N | −0.697 | −0.551 | −1.075 (drawn at −1) |
| P02-2022 | −0.226 | −0.118 | −0.040 |
| I09N-2025 | 0.004 | −0.051 | 0.251 |
| GOMECC4 | 0.387 | 0.398 | 0.460 |
| I09N | 0.484 | 0.446 | 0.820 |
| C13.5/A13.5 | 0.493 | 0.489 | 0.475 |
| A16N-2023 | 0.473 | 0.527 | 0.474 |
| I05-2023 | 0.610 | 0.572 | 0.598 |
| A13.5-2024 | 0.613 | 0.576 | 0.647 |
| AMT28 | 0.629 | 0.616 | 0.557 |
| P18 | 0.617 | 0.622 | 0.641 |
| HLY2202-2022 | 0.732 | 0.738 | 0.818 |
| I08S-2024 | 0.905 | 0.885 | 0.865 |

**Variance partition, main models** (share of total variance, %; `partition.csv`, `figures/output/variance_partition.png`)

| Target | Model | Technical | Explained | Cruise offsets | Within-tile | Within-cruise | Explained / non-technical | Explained / reproducible non-technical |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| All KOs | Granite | 9.09 | 16.00 | 12.19 | 27.04 | 35.59 | **17.60** | **25.05** |
| C | OLCI + SST + MLD30 | 15.95 | 42.56 | 7.82 | 6.49 | 27.05 | **50.63** | **54.87** |
| C (alt.) | Granite | 15.95 | 45.22 | 5.30 | 6.49 | 26.96 | 53.80 | 58.30 |
| D | Granite | 7.82 | 12.87 | 12.79 | 27.00 | 39.43 | 13.96 | 19.75 |
| A | Granite | 11.69 | 12.11 | 13.97 | 25.92 | 36.21 | 13.71 | 19.41 |
| Dp | Granite | 11.27 | 14.51 | 13.12 | 23.47 | 37.54 | 16.36 | 22.24 |
| B | Granite | 6.99 | 9.60 | 9.39 | 57.60 | 16.37 | 10.32 | 27.11 |
| De | Granite | 8.56 | 17.24 | 10.22 | 45.06 | 18.84 | 18.86 | 37.18 |

- Rows sum to 99.87–99.95% (the held-out residual mean is not zero).
- Supersedes R56's two-stage numbers (All KOs 17.01 / 24.21 → 17.60 / 25.05 with Granite alone).

**Per-KO scatter and heatmap, remade** (`figures/make_temperature_gap.py --main` → `temperature_gap_main.png`;
`figures/make_ko_heatmap.py --main` → `ko_heatmap_main.png`; per-KO R² on technically adjusted CLR, Granite PCs with
`cv` penalty vs temperature spline)
- 3,933 KOs. Median R²(Granite) − R²(temperature) +0.083. Granite ahead by ≥ 0.1 in 1,769 KOs, temperature in 58.
- "Clearly better" (≥ 0.1 and ≥ 9 of 13 cruises): Granite 1,020, temperature 38, neither 2,875.
- Heatmap blocks (unchanged thresholds): Granite clearly better 55 (28 shown), both high and equal 40 (14 shown),
  temperature clearly better 10. Right column is now Granite PCs + temperature from the same run; no footnote.
- Temperature alone loses most of its per-KO skill once platform and depth are removed (median 0.0089), because
  platform is confounded with region and hence temperature.

**Within-tile variability by samples per tile** (`within_tile_k.csv`, `within_tile_gain.csv`, `scripts/within_tile_k.py`).
Question: as k goes 2 → 3 → 5–6, does the tile mean get closer to "truth"?
- Single-sample within-tile variance σ² is about the same whatever k (All KOs: 918, 790, 857 for k = 2, 3, 4–6), so the
  variance left in the tile mean falls as σ²/k: 459 → 263 → 150 (sums over KOs).
- Inside the same tiles, the held-out error against the tile mean is lower than against a single sample by exactly
  σ²(1 − 1/k) (an identity, not an independent test): All KOs −19.9% (k = 2), −20.8% (k = 3), −18.5% (k = 4–6); B −29% to
  −34%; C −5% to −8%.
- The remaining error against the tile mean (All KOs 1,843 / 2,005 / 3,084) is far above σ²/k, so it is mostly model
  error. Beyond 2–3 samples, extra samples per tile reduce the target scatter further but cannot raise skill much.
- Error by k across tiles is confounded with cruise: k ≥ 3 tiles are almost all GOMECC4 (48 of 49).

## R60. One model for everything: C on Granite PCs; time & place and raw satellite as comparisons; heatmap of four inputs (2026-11-23)

**Decision (user).** C is predicted with the same model as every KO target: 32 Granite PCs, ridge, penalty by
leave-cruise-out inside the training cruises (`cv`), technical variance (platform + read-depth spline) removed
in-sample. OLCI + SST + 30-day MLD is dropped; it is reported in a few sentences only (R59: 50.6% vs Granite 53.8%; the
MLD gain on C is +4.7 pp on raw C and +0.3 pp after technical removal). Technical removal applies to every model and
every figure.

**Comparison inputs, all on technically adjusted targets, all with the `cv` penalty**
- Temperature alone: in-situ temperature spline. It was GCV in R59. With `cv`, All KOs is 11.12% (GCV 10.37%).
- Raw satellite: 16 OLCI band means + satellite SST, quantile-normalised, splines.
- Time & place: lat / 90, **|lat| / 90**, sin and cos of longitude, sin and cos of day of year; standardised, linear.
  - R5's version without |lat| scores below zero once cruises are held out: C −14.96%, All KOs −4.67%; for C the
    slope of observed on predicted is −2.9. It cannot follow the increase toward both poles.
  - Variants for C / All KOs: linear 5 features −14.96 / −4.67; + |lat| 17.17 / 6.51; latitude spline + |lat|
    26.02 / −1.68. The + |lat| version is used throughout, as the simplest one that is positive for both.

**Explained, % of non-technical variance (pooled, cruise held out)**

| Target | Granite PCs | Granite + temperature | Raw satellite | Temperature | Time & place |
|---|---:|---:|---:|---:|---:|
| All KOs | **17.60** | 17.51 | 13.33 | 11.12 | 6.51 |
| Eukaryote : prokaryote reads (C) | **53.80** | 55.30 | 50.32 | 45.86 | 17.17 |
| Shared KOs (D) | **13.96** | 13.71 | 8.64 | 5.99 | 6.19 |
| Prokaryote-specific (A) | **13.71** | 13.93 | 9.47 | 7.19 | 6.56 |
| Shared, prokaryote part (Dp) | **16.36** | 16.34 | 9.96 | 7.06 | 9.47 |
| Eukaryote-specific (B) | **10.32** | 10.07 | 8.82 | 7.18 | 2.20 |
| Shared, eukaryote part (De) | **18.86** | 18.96 | 18.56 | 17.99 | 9.75 |

- Granite PCs beat every other single input on every target. Adding temperature changes nothing beyond ±0.3 pp, except
  for C (+1.5 pp).

**Per-KO median R² (All KOs, 3,933 KOs)**
- Granite PCs 0.1157; Granite + temperature 0.1111; raw satellite 0.0682; time & place 0.0387; temperature 0.0069.
- Temperature reaches 11.12% pooled but a median of 0.007 per KO. Its skill sits in a few high-variance KOs (the
  cyanobacterial photosynthesis / clock block: ndhD 0.76, hoxH 0.75, kaiB 0.71, cpeS 0.69). For the typical KO it
  carries almost nothing once platform and depth are removed.

**C variance partition** (`partition.csv`): technical 15.95, explained 45.22, cruise offsets 5.30, within-tile 6.49,
within-cruise 26.96 (% of total) → **53.80% of non-technical, 58.30% of reproducible non-technical**. The KO rows are
unchanged from R59.

**C per cruise** (R² within the held-out cruise; `c_per_cruise.csv`)

| Cruise | Granite PCs | OLCI + SST | Time & place |
|---|---:|---:|---:|
| I07N | −1.075 | −0.697 | −1.558 |
| P02-2022 | −0.040 | −0.226 | −0.566 |
| I09N-2025 | 0.251 | 0.004 | 0.158 |
| GOMECC4 | 0.460 | 0.387 | −0.775 |
| A16N-2023 | 0.474 | 0.473 | 0.334 |
| C13.5/A13.5 | 0.475 | 0.493 | 0.208 |
| AMT28 | 0.557 | 0.629 | 0.438 |
| I05-2023 | 0.598 | 0.610 | 0.450 |
| P18 | 0.641 | 0.617 | −0.068 |
| A13.5-2024 | 0.647 | 0.613 | 0.405 |
| HLY2202-2022 | 0.818 | 0.732 | 0.395 |
| I09N | 0.820 | 0.484 | 0.435 |
| I08S-2024 | 0.865 | 0.905 | 0.731 |

**Figures**
- `figures/output/variance_partition.png`: the main figure. Tick labels now use names, not letters.
- `figures/output/c_prediction.png`: overwritten.
  - a: observed vs held-out predicted, Granite PCs, R² 0.54, 1,215 tiles.
  - b: per-cruise R² for Granite PCs, OLCI + SST, and time & place (* = NovaSeq 6000 cruise; values below −1 drawn
    at −1).
  - c: against latitude, 41-tile running means.
  - Axes read "log(eukaryote / prokaryote reads)".
- `figures/output/ko_heatmap_inputs_main.png` (`figures/make_ko_heatmap.py --main --columns inputs`): four inputs
  (Granite PCs, in-situ temperature, time & place, raw satellite) plus Granite + temperature at the right.
  - Technical variance is removed throughout; no field block, no footnote.
  - Blocks: Granite clearly better 54 (28 shown), both high and about equal 40 (14 shown), temperature clearly
    better 11.
  - The R59 two-column heatmap (`ko_heatmap_main`) is deleted.
- `figures/output/temperature_gap_main.png`: redone with the `cv` temperature model.
  - Median difference +0.066.
  - Granite ahead by ≥ 0.1 in 1,573 KOs, temperature in 51.
  - "Clearly better" (also ≥ 9 of 13 cruises): Granite 835, temperature 38, neither 3,060.

**Candidate marker genes** (per-KO R², technical removed; Granite / temperature / raw satellite)

| Group | Gene | KO | Granite | Temperature | Raw satellite |
|---|---|---|---:|---:|---:|
| Phosphorus | pstS | K02040 | 0.40 | 0.28 | 0.46 |
| Phosphorus | pstB | K02036 | 0.41 | 0.23 | 0.45 |
| Phosphorus | phnD | K02044 | 0.38 | 0.41 | 0.33 |
| Phosphorus | phoB | K07657 | 0.17 | −0.27 | 0.24 |
| Phosphorus | phoD | K01113 | −0.02 | −0.07 | −0.06 |
| Phosphorus | phoA | K01077 | 0.01 | −0.18 | 0.05 |
| Nitrogen | narB | K00367 | 0.44 | 0.38 | 0.27 |
| Nitrogen | ureC | K01428 | 0.28 | 0.36 | 0.26 |
| Nitrogen | amt | K03320 | 0.26 | 0.12 | 0.20 |
| Nitrogen | urtA | K11959 | 0.19 | 0.19 | 0.16 |
| Iron | fiu | K16090 | 0.36 | 0.51 | 0.29 |
| Iron | petF | K02639 | 0.27 | 0.25 | 0.17 |
| Iron | isiB / fldA | K03839 | 0.05 | −0.07 | −0.02 |
| Light | cpeS | K05382 | 0.67 | 0.69 | 0.55 |
| Light | pufL | K08928 | 0.25 | 0.46 | 0.44 |
| Light | puhA | K13991 | 0.30 | 0.45 | 0.42 |
| Light | rbcL | K01601 | 0.17 | 0.12 | 0.05 |
| Light | psbA | K02703 | 0.03 | 0.04 | −0.03 |

- Not among the 3,933 modelled KOs: nifH, nifD, amoA / pmoA, nasA, narG, phnJ, bacteriorhodopsin.
- 620 KOs have Granite R² > 0.3.

## R61. Poorly predicted cruises for eukaryote : prokaryote reads; temperature on top of Granite; pathway-level predictability (2026-11-23)

**Why I07N, P02-2022 and I09N-2025 score badly on eukaryote : prokaryote reads (C)** (Granite PCs, technical removed)

| Cruise | Where | When | In-situ T (°C) | Tiles | C mean | C SD within cruise | Predicted SD | r (obs, pred) | Bias (pred − obs) | R² in cruise |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| I07N | W Indian Ocean, 30 °S–18 °N, 40–69 °E (incl. Arabian Sea) | 2018-04-25 to 06-03 | 22.4–31.0 | 146 | −3.41 | 0.44 | 0.35 | −0.10 | −0.28 | −1.075 |
| P02-2022 | N Pacific zonal line, 27–33 °N, dateline crossing | 2022-05-03 to 07-12 | 17.8–25.5 | 89 | −3.36 | 0.64 | 0.30 | 0.26 | −0.17 | −0.040 |
| I09N-2025 | E Indian Ocean, 29 °S–20 °N, 87–100 °E | 2025-03-23 to 04-23 | 22.1–30.9 | 84 | −3.61 | 0.58 | 0.29 | 0.36 | −0.19 | 0.251 |
| (all tiles) | | | | 1,215 | | 1.05 | | | | |

- All three are warm, low-latitude or subtropical lines with little range in C: within-cruise SD 0.44–0.64, against
  1.05 over all tiles and 1.0–1.4 for the cruises that reach high latitudes (P18, A16N-2023, A13.5-2024, I08S-2024,
  HLY2202-2022).
- The model gets the cruise level roughly right (bias −0.17 to −0.28) but hardly tracks the variation along the line
  (r −0.10 to 0.36; predictions compressed, SD 0.29–0.35).
- The per-cruise R² is scored against the training-cruise mean. A cruise whose mean C sits near the global mean and
  varies little has a small denominator, so small errors give R² ≤ 0. GOMECC4 shows the reverse: r 0.05, but its mean
  is far from the rest (−4.40), so it gets R² 0.46 from the level alone. Per-cruise R² mostly measures whether a
  cruise has a gradient the model can follow.
- The variation of C along these warm lines (e.g. Arabian Sea at monsoon onset) is not something the surface image
  resolves at present.

**Granite PCs + in-situ temperature vs Granite alone** (pp of non-technical variance, cv penalty, from R59 / R60):
All KOs −0.09, eukaryote : prokaryote reads +1.50, shared KOs −0.25, prokaryote-specific +0.22, shared prokaryote part
−0.02, eukaryote-specific −0.25, shared eukaryote part +0.10. Only C gains more than 0.3 pp.

**Pathway-level predictability** (`scripts/module_predictability.py`, `figures/output/module_predictability.png`;
`module_predictability.csv`, `category_predictability.csv`, `ko_ceiling.csv`)
- Per KO: ceiling = 1 − within-tile variance / non-technical variance across tiles.
  - Within-tile variance = σ² from the sample pairs × Σ 1/k over tiles.
  - Median ceiling 0.740 (10–90%: 0.498–0.864). 25 KOs fall below 0.1 and are floored at 0.1.
- Fraction = held-out R² (Granite PCs) / ceiling. KO median: 0.168 for Granite PCs, 0.010 for temperature.
- KEGG modules (REST `link/module/ko` and `br:ko00002`, downloaded 2026-11-23 to `data/reference/kegg_modules/`):
  172 modules with ≥ 4 modelled KOs. The median over modules is 0.186 for Granite PCs and −0.007 for temperature.
  - Top modules:
    - NAD(P)H:quinone oxidoreductase, chloroplasts / cyanobacteria (M00145): 0.684; R² 0.605, ceiling 0.920.
    - F-type ATPase, eukaryotes: 0.478.
    - Thiamine salvage: 0.456; thiamine biosynthesis (archaea): 0.456.
    - Pyrimidine ribonucleotide biosynthesis: 0.427.
    - Isoleucine biosynthesis: 0.418; leucine biosynthesis: 0.395.
    - Prokaryotic thiamine biosynthesis: 0.378; methionine degradation: 0.378; C5 isoprenoid (MEP): 0.377;
      serine-pathway formaldehyde assimilation: 0.377.
  - Bottom: methanogenesis, catechol ortho-cleavage, glucuronate degradation, eukaryotic ubiquinone, CDP-sugar
    biosynthesis (ceiling ≈ 0) — all ≤ 0.
  - Class medians (module count in brackets):
    - Highest: branched-chain amino acids 0.355 (5), cysteine / methionine 0.337 (3), arginine / proline 0.330 (6),
      fatty acids 0.266 (4), pyrimidines 0.253 (5), lysine 0.239 (7), cofactors and vitamins 0.224 (34),
      ATP synthesis 0.210 (11).
    - Lowest: central carbohydrate 0.119 (12), nucleotide sugars 0.111 (15), photosynthesis 0.087 (2),
      sulfur 0.086 (2), aromatics degradation 0.061 (4), glycosaminoglycans −0.026 (2).
  - Temperature reaches a higher module fraction than Granite in only a few modules: cyanobacterial NDH 0.718 vs 0.684,
    and eukaryotic F-type ATPase 0.512 vs 0.478.
- KEGG modules cover metabolism only, so panel b uses BRITE level-B categories (all KOs, ≥ 15 per category).
  Median fraction, Granite PCs [IQR], with temperature after the semicolon:
  - Highest:
    - Cellular community – eukaryotes 0.502 [0.372–0.556]; 0.446.
    - Signalling molecules 0.444; 0.346.
    - Transcription 0.364; 0.183.
    - Transport and catabolism 0.355; 0.198.
    - Cell growth and death 0.280; 0.233.
    - Folding / sorting / degradation 0.278; 0.044.
    - Signal transduction 0.248; 0.146.
    - Replication and repair 0.213; 0.016.
    - Translation 0.210; −0.042.
  - Middle:
    - Cofactors and vitamins 0.187; nucleotides 0.185; membrane transport 0.173 (temperature 0.005); amino acids
      0.159; energy 0.146.
  - Lowest:
    - Carbohydrate 0.117; prokaryotic cellular community (biofilm, quorum sensing) 0.113; xenobiotics 0.084;
      cell motility 0.038.
- **Caveat on interpretation:** the top BRITE categories are mostly eukaryote-only KOs (cell-cycle, signalling,
  vesicle transport). In the all-KO CLR, their skill largely rides on the eukaryote : prokaryote balance, which is the
  most predictable quantity (54%). The module ranking within prokaryote metabolism (amino-acid and cofactor
  biosynthesis high, carbohydrate and degradation pathways low) is less exposed to this. A version on the domain-split
  targets (prokaryote-specific per genome equivalent, eukaryote-specific within eukaryotes) would separate the two.

**Ocean-model fields added to the main model** (`scripts/main_fields.py`, run in 3 shards; `added_fields.csv`,
`added_fields_summary.csv`).
- Redo of R28 with Granite PCs as the base (no temperature), the cv penalty and technically adjusted targets.
- "Hit" = the heatmap rule: gain ≥ 0.05, R² after ≥ 0.3, and better in ≥ 10 of 13 held-out cruises. The R28
  hits (on the old baseline) are given in brackets.

| Field | KOs gaining ≥ 0.05 | Hits (R28) | Median gain |
|---|---:|---:|---:|
| sea level | 357 | 57 (6) | +0.0008 |
| phosphate | 256 | 49 (36) | −0.0016 |
| nitrate | 238 | 22 (22) | −0.0001 |
| GLORYS temperature | 143 | 14 (1) | −0.0015 |
| oxygen | 129 | 8 (2) | −0.0008 |
| silicate | 63 | 7 (11) | −0.0011 |
| chlorophyll | 121 | 6 (14) | +0.0026 |
| salinity | 37 | 5 (7) | −0.0010 |
| mixed layer | 67 | 2 (8) | −0.0009 |
| mixed layer, 30-day max | 138 | 2 (10) | −0.0011 |
| wind, 14-day mean | 68 | 2 (3) | +0.0009 |
| wind, 14-day peak | 73 | 1 (16) | −0.0003 |
| primary production | 30 | 0 (1) | +0.0005 |
| SST change 14 d / current speed / sea-level gradient / EKE 30 d | 0–2 | 0 | ~0 |

- Forecast-period fields (iron, DIC, alkalinity, pH, pCO2, phytoplankton carbon; 636 tiles, 6 cruises) cannot meet
  the ≥ 10-cruise rule. Their KO counts gaining ≥ 0.05 are 259 / 100 / 288 / 86 / 104 / 35.
- Controls: three noise columns and shuffled phosphate give 0 KOs gaining ≥ 0.05 on all tiles. Noise on the forecast
  tiles gives 0, shuffled phosphate on the forecast tiles 4. All controls give 0 hits.
- 154 distinct KOs are hits (best field per KO): sea level 57, phosphate 39, nitrate 18, GLORYS temperature 14,
  silicate 6, oxygen 5, chlorophyll 5, salinity 4, mixed layer 2, 30-day max mixed layer 2, wind mean 2.
- **The field gains survive technical removal.** Phosphate and nitrate keep their hits. Sea level and model
  temperature gain more than before, because temperature is no longer in the base.
- Heatmap field block (top 16 by gain, not already shown above; 144 eligible):
  - + nitrate: the whole urease / urea-transport set. ureD 0.10 → 0.55, ureF 0.15 → 0.57, ureC 0.28 → 0.64,
    ureB, ureE, ureA, urtA–D (+0.30 to +0.45); also mqo, menC.
  - + phosphate: the pho regulon. phoU 0.17 → 0.53, pstC 0.33 → 0.63, pstS 0.40 → 0.67, adeR 0.26 → 0.52.
- `figures/output/ko_heatmap_inputs_main.png` now has the four input blocks plus this field block.

## R62. Per-cruise error as level + pattern instead of per-cruise R²; module figure fixes (2026-11-23)

**Decisions (user)**
- Granite PCs alone is the main model. The C gain from adding temperature (+1.5 pp) gets one sentence.
- Pooled R² stays the headline. Per-cruise R² is dropped: it is scored against the training-cruise mean, so it rewards
  cruises whose mean is far from the rest (GOMECC4) and penalises flat ones (R61).

**Per-cruise reporting** (`main_model.py --parts cruise`; `cruise_decomposition.csv` for every target,
`c_per_cruise.csv` for C). For each cruise and target column:

MSE = bias² + (SD_pred − SD_obs)² + 2 · SD_pred · SD_obs · (1 − r)

- bias = mean(observed − predicted); SDs and r are within the cruise (population moments). The identity is exact:
  max |sum / MSE − 1| = 8.9e−16.
- **Level** = bias. **Pattern** = within-cruise r. Each cruise's within-cruise SD (observed and predicted) is shown
  beside them. Shares are given as parts of each cruise's MSE.

**Eukaryote : prokaryote reads (C), Granite PCs, sorted by within-cruise SD of the observed values** (* NovaSeq 6000
cruise; the last two columns are bias / r for OLCI + SST and for time & place)

| Cruise | Tiles | SD obs | SD pred | Level (bias) | Pattern (r) | Level share | Spread share | Pattern share | MSE | OLCI + SST | Time & place |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| * I09N | 54 | 0.364 | 0.131 | −0.106 | 0.124 | 0.08 | 0.36 | 0.56 | 0.149 | −0.410 / 0.140 | −0.583 / 0.170 |
| * I07N | 146 | 0.434 | 0.349 | +0.283 | −0.104 | 0.19 | 0.02 | 0.79 | 0.422 | +0.339 / 0.108 | +0.204 / −0.315 |
| GOMECC4 | 53 | 0.486 | 0.414 | −0.668 | 0.052 | 0.54 | 0.01 | 0.46 | 0.833 | −0.740 / 0.053 | −1.582 / 0.124 |
| * C13.5/A13.5 | 216 | 0.487 | 0.485 | +0.043 | 0.531 | 0.01 | 0.00 | 0.99 | 0.224 | +0.004 / 0.470 | +0.050 / 0.200 |
| I09N-2025 | 84 | 0.576 | 0.285 | +0.195 | 0.362 | 0.11 | 0.25 | 0.63 | 0.332 | +0.234 / 0.391 | −0.101 / 0.130 |
| P02-2022 | 89 | 0.640 | 0.300 | +0.168 | 0.256 | 0.07 | 0.27 | 0.67 | 0.430 | +0.279 / 0.148 | −0.482 / −0.125 |
| * AMT28 | 46 | 0.825 | 0.918 | +0.264 | 0.827 | 0.21 | 0.03 | 0.77 | 0.340 | +0.237 / 0.855 | +0.141 / 0.649 |
| I05-2023 | 141 | 0.939 | 0.465 | +0.290 | 0.449 | 0.11 | 0.28 | 0.61 | 0.789 | −0.037 / 0.369 | +0.473 / 0.302 |
| I08S-2024 | 23 | 0.996 | 1.182 | −0.036 | 0.909 | 0.01 | 0.14 | 0.86 | 0.249 | +0.013 / 0.935 | +0.195 / 0.901 |
| A13.5-2024 | 93 | 1.115 | 0.853 | −0.227 | 0.831 | 0.12 | 0.16 | 0.73 | 0.442 | −0.197 / 0.804 | +0.079 / 0.682 |
| A16N-2023 | 174 | 1.121 | 0.574 | −0.317 | 0.746 | 0.14 | 0.41 | 0.45 | 0.727 | −0.469 / 0.830 | +0.006 / 0.633 |
| * P18 | 64 | 1.190 | 1.128 | +0.321 | 0.782 | 0.15 | 0.01 | 0.85 | 0.691 | +0.425 / 0.786 | +1.069 / 0.602 |
| * HLY2202-2022 | 32 | 1.353 | 0.366 | +0.060 | 0.457 | 0.00 | 0.64 | 0.36 | 1.515 | +0.310 / 0.005 | +1.809 / 0.368 |

- **The model places cruises well and follows the gradient where there is one.**
  - Level is within ±0.33 for 12 of 13 cruises, against an overall SD of 1.05.
  - Every cruise with within-cruise SD ≥ 0.8 has r 0.45–0.91; seven of those eight have r 0.75–0.91.
- **On flat, warm lines it misses the pattern:** I09N (2016) r 0.12, I07N −0.10, GOMECC4 0.05, I09N-2025 0.36,
  P02-2022 0.26. C13.5/A13.5, also flat, is the exception (r 0.53).
- The 2016 I09N line is flat too. Its R² of 0.82 in R60 came from its level alone.
- **GOMECC4 is the one cruise placed wrong:** bias −0.67 (observed below predicted). Level is 54% of its MSE. The Gulf
  of Mexico in Sep–Oct 2021 had the lowest eukaryote share of any cruise (mean −4.40). OLCI + SST and time & place
  misplace it too (−0.74, −1.58).
- **HLY2202-2022:** right level (+0.06), but the predicted spread is 0.37 against 1.35 observed. Its error is mostly
  under-dispersion (spread share 0.64).

**All KOs (CLR), Granite PCs: medians over KOs per cruise; shares of the cruise's summed MSE**

| Cruise | SD obs | SD pred | abs(bias) | r | Level | Spread | Pattern |
|---|---:|---:|---:|---:|---:|---:|---:|
| I09N-2025 | 0.250 | 0.147 | 0.084 | 0.261 | 0.14 | 0.25 | 0.61 |
| P02-2022 | 0.314 | 0.143 | 0.076 | 0.205 | 0.10 | 0.31 | 0.59 |
| GOMECC4 | 0.407 | 0.147 | 0.274 | 0.158 | 0.34 | 0.28 | 0.38 |
| I08S-2024 | 0.432 | 0.233 | 0.099 | 0.635 | 0.13 | 0.28 | 0.58 |
| A13.5-2024 | 0.443 | 0.218 | 0.072 | 0.521 | 0.08 | 0.35 | 0.57 |
| AMT28 | 0.457 | 0.232 | 0.131 | 0.369 | 0.20 | 0.26 | 0.54 |
| I07N | 0.479 | 0.126 | 0.118 | 0.099 | 0.13 | 0.44 | 0.43 |
| I05-2023 | 0.545 | 0.170 | 0.253 | 0.214 | 0.25 | 0.36 | 0.39 |
| C13.5/A13.5 | 0.561 | 0.178 | 0.082 | 0.392 | 0.07 | 0.53 | 0.40 |
| I09N | 0.580 | 0.087 | 0.165 | 0.035 | 0.20 | 0.57 | 0.22 |
| HLY2202-2022 | 0.585 | 0.120 | 0.308 | 0.184 | 0.40 | 0.38 | 0.22 |
| A16N-2023 | 0.665 | 0.211 | 0.181 | 0.336 | 0.11 | 0.44 | 0.45 |
| P18 | 0.690 | 0.251 | 0.175 | 0.489 | 0.17 | 0.42 | 0.41 |

- For single KOs, the level error is a small part of the MSE (7–40%). Most of the error is pattern and under-dispersion:
  predicted SD is 0.09–0.25 against 0.25–0.69 observed. Ridge shrinkage plus weak per-KO signal.
- The highest-r cruises for KOs (I08S-2024 0.635, A13.5-2024 0.521, P18 0.489) are the high-gradient ones, as for C.
- Tables for D, A, Dp, B, De are in `cruise_decomposition.csv` and `summary_cruise.txt`.

**C figure** (`figures/output/c_prediction.png`, overwritten). Cruises in c–e are sorted by observed within-cruise SD;
c and d also show OLCI + SST and time & place as open markers.
- a: observed vs held-out predicted (R² 0.54).
- b: against latitude.
- c: level, mean(observed − predicted).
- d: pattern, within-cruise r.
- e: within-cruise SD, observed and predicted.

**Module figure fixes** (`figures/output/module_predictability.png`)
- Panel a has x-tick labels and a dot-size legend (KOs in module: 5, 10, 20).
- "Metabolic capacity" is a KEGG class from the Signature modules branch ("Signature modules > Module set > Metabolic
  capacity"). Signature modules are marker sets, not pathways. Most are defined by other modules (e.g. M00611
  oxygenic photosynthesis), so they have no direct KO links. Only M00615 nitrate assimilation links to KOs here (4 KOs:
  fraction 0.253, R² 0.016, ceiling 0.462). Its row is now labelled "Signature module: nitrate assimilation".

## R63. Variance budget on the main model; C figure with per-cruise R² restored; figure titles and legends (2026-11-23)

**Variance budget** (`scripts/variance_budget_main.py`; `data/results/main_model/variance_budget.csv`;
`figures/output/variance_budget.png`). The R43 waterfall redone on the main model: all KOs, Granite PCs, leave-cruise-out
penalty, technical variance removed in-sample. 1,215 tiles × 3,933 KOs; parts sum to 100.00%.

| Component | Raw CLR | KOs standardised | R43 (raw) |
|---|---:|---:|---:|
| Explained by Granite PCs (cruise held out) | 16.00 | 12.41 | 16.2 |
| Sequencing platform (NovaSeq 6000 vs X; upper bound) | 4.54 | 5.38 | 3.8 |
| Sequencing depth | 4.55 | 5.90 | 2.3 |
| Other cruise offsets (unknown) | 12.29 | 13.41 | 13.1 |
| Axis 1 within cruise (eukaryotes vs prokaryotes) | 6.77 | 9.68 | 7.4 |
| Axis 2 within cruise (phototrophs; driver unknown) | 3.73 | 6.03 | 3.6 |
| Within-tile variability (unpredictable floor) | 27.04 | 26.01 | 27.0 |
| Remaining within cruise (unknown) | 25.08 | 21.17 | 26.6 |

- Platform and depth are now the in-sample technical fit, which is removed before modelling:
  - Platform = OLS on the platform group.
  - Depth = the rest of the joint platform + depth-spline fit.
  - The two sum to the 9.09% technical share of R59 / R60.
  - In R43 they were taken from the residual after the model instead.
- The axes are the two leading principal components of the new within-cruise residual (KOs standardised), fitted
  in-sample.
  - Axis 1: r 0.941 with R36 W1 (0.042 with W2); r −0.655 with the within-cruise log(eukaryote / prokaryote reads).
    Its loadings run from bacterial cell-wall and replication genes (bacA, murB, murC, holB, ubiE) to eukaryotic
    nucleolar and RNA-processing genes (NOG1, DKC1, UPF1, STT3, TSR3).
  - Axis 2: r −0.914 with R36 W2 (0.037 with W1); r 0.165 with eukaryote : prokaryote reads. Its negative end
    includes petB (cytochrome b6), desC and SulP.
- Cruise offsets include the overall mean of the held-out residual (0.09% of total). With that, the parts sum exactly.

**C figure.** `figures/output/c_prediction.png` is again the per-cruise R² version (a scatter, b per-cruise R² for
three inputs, c latitude), now with Granite PCs. The level / pattern / SD version (R62) is
`figures/output/c_prediction_split.png`, for the supplement. Per-cruise R² values are in `c_per_cruise_r2.csv` and the
split in `c_per_cruise.csv`.

### Figure titles and legends (drafts)

**Figure — Variance budget.** *Where the variance in gene content goes.*
Variance across 1,215 satellite tiles of the centred log-ratio abundances of 3,933 KEGG orthologues (KOs), summed over
KOs and partitioned sequentially (bars run left to right and sum to 100%). Explained: held-out skill of a ridge regression on 32
principal components of the Granite embedding, with each cruise held out in turn (13 cruises) and the ridge penalty chosen
by leave-one-cruise-out within the training cruises; fitted after removing technical variance. Sequencing platform:
NovaSeq 6000 vs NovaSeq X cruise groups. This is an upper bound, because the platform groups also differ in year and
region. Sequencing depth: cubic spline of log10 KO-assigned reads per tile. Platform and depth are fitted in-sample.
Other cruise offsets: per-cruise means of the held-out residual. Axes 1 and 2: the two leading principal components of
the within-cruise residual. Axis 1 separates eukaryote- from prokaryote-dominated samples; axis 2 loads on phototroph
genes and has no identified driver. Within-tile variability: variance between samples sharing a tile, estimated from
305 same-tile sample pairs and scaled by the number of samples averaged per tile. It cannot be predicted from a
tile-level input. Remaining: within-cruise variation not accounted for above.

**Figure — The eukaryote : prokaryote balance is predictable from space.**
*(Updated R68: C is now the direct read ratio.)* **a**, Observed vs predicted log ratio of eukaryote- to
prokaryote-assigned KO reads (all KOs; each read takes the domain of the reference gene it maps to) for 1,215 tiles. Predictions come
from 32 Granite principal components (ridge), with each cruise held out in turn. Technical variance (sequencing
platform, depth) has been removed from the observations. Colour shows sequencing platform; the dashed line is 1:1.
R² = 0.54 pooled over all held-out predictions (53.8% of the non-technical variance). **b**, R² within each held-out
cruise for Granite PCs (circles), OLCI band means + SST (triangles) and time & place (squares: latitude, |latitude|,
longitude, day of year). Per-cruise R² is relative to the mean of the training cruises, so flat cruises near the
global mean score low even when placed correctly (Supplementary Fig. X). Values below −1 are drawn at −1 (time &
place on I07N, −1.32). * NovaSeq 6000 cruises. **c**, Observed (grey points; black line) and predicted
(blue line) log ratio against latitude; lines are 41-tile running means.

**Supplementary Figure X — Per-cruise error as level and pattern.**
Each held-out cruise's mean squared error splits exactly into bias², a spread mismatch (SD_pred − SD_obs)², and a
pattern term 2·SD_pred·SD_obs·(1 − r). **a**, **b**, As in the main figure (a: observed vs predicted; b: against
latitude). **c**, Level: mean of observed − predicted per cruise. **d**, Pattern: within-cruise correlation between
observed and predicted. **e**, Within-cruise SD of the observed (grey) and predicted (blue) values. Cruises are
sorted by observed within-cruise SD. Filled circles: Granite PCs; open triangles: OLCI + SST; open squares: time &
place. The model places 12 of 13 cruises within ±0.33 (overall SD 1.05). It follows the gradient where there is one:
r 0.45–0.91 for cruises with SD ≥ 0.8, and 0.75–0.91 for seven of those eight. It misses the within-line pattern on
flat, warm cruises (I09N 2016, I07N, GOMECC4, I09N-2025, P02-2022). GOMECC4 is the one cruise placed wrong (−0.67).
HLY2202-2022 is under-dispersed (predicted SD 0.37 vs 1.35). * NovaSeq 6000 cruises.

**Figure — What the satellite predicts, gene by gene.**
Pooled held-out R² per KO (rows) from four inputs (columns): 32 Granite principal components (ridge); in-situ
temperature (cubic spline); time & place (latitude, |latitude|, longitude, day of year; linear); and raw satellite
values (16 OLCI band means + SST, rank-normalised splines). Each cruise is held out in turn (13 cruises, 1,215 tiles),
the ridge penalty is chosen by leave-one-cruise-out within the training cruises, and technical variance (platform,
depth) is removed from every KO. Right-hand numbers: R² of Granite PCs + temperature. Blocks:
- *Granite clearly better*: Granite R² > 0.3, at least 0.15 above temperature, and ahead in ≥ 10 of 13 cruises
  (28 shown of 54, highest Granite R²).
- *Both high, about equal*: both R² ≥ 0.6 and differing by < 0.05 (14 of 40).
- *Temperature clearly better*: temperature R² > 0.3, at least 0.1 above Granite, and ahead in ≥ 9 cruises
  (all 11).
- *Improved by adding an ocean-model field*: adding one Copernicus GLORYS / PISCES field for the sampling day (spline)
  to the Granite PCs gains ≥ 0.05, reaches R² ≥ 0.3 and lowers error in ≥ 10 of 13 cruises. The 16 largest gains of
  144 KOs are shown (not repeated from the blocks above). The second column is then Granite PCs + the field, and the
  right-hand text gives the change. Random-noise and shuffled-field controls produce no such KOs.

Within blocks, rows are grouped by KEGG BRITE functional category. Colour: R² from ≤ 0 (pale) to 1 (dark).

**Figure — Which kinds of function are visible from space.**
Skill as a share of what is predictable: for each KO, held-out R² (Granite PCs, as above) divided by its ceiling, the
share of its non-technical variance that is not within-tile variability. Median ceiling 0.74; 25 KOs with ceilings
below 0.1 are set to 0.1. **a**, KEGG modules with ≥ 4 modelled KOs (172 modules). Each dot is a module's median
over its KOs; dot size is the number of KOs. Rows are KEGG module classes, ordered by the median over their modules
(black tick); the class's best module is named at the right. **b**, KEGG BRITE level-B categories with ≥ 15 KOs,
covering functions outside modules (transport, genetic information processing, signalling). Median (filled) and
interquartile range over KOs for Granite PCs, and median for in-situ temperature (open). The median module scores 0.19 with
Granite PCs and −0.01 with temperature. Categories at the top of b consist largely of eukaryote-specific genes, whose
abundance follows the eukaryote : prokaryote balance.

## R64. Benchmark bundle for the new foundation model (2026-11-23)

- **Instructions:** `~/multimodal_seusmbol/instructions_agent_biogoship.md`, written for the agent working on the
  foundation model.
- **Bundle:** `scripts/export_benchmark.py` → `data/results/benchmark_export/` (zipped as
  `data/results/benchmark_export.zip`, 189 MB). It contains:
  - `samples.csv`: 1,990 samples, 1,416 in the benchmark, with the Granite match per sample.
  - `tiles.csv`: 1,215 tiles; `in_B_De` 1,169.
  - Raw tile-level and sample-level targets.
  - The 305 within-tile pairs.
  - Granite embeddings (1,215 × 512), raw satellite features, held-out predictions and all score tables.
  - KEGG references.
- **Tile keys:** the loader renumbers `tile_key` for each target set. The B and De tiles are mapped to the all-KO
  keys through the embedding row.
- **Standalone evaluation (`scripts/benchmark_eval.py`, copied into the bundle):** run from the bundle alone with
  `python -I`, it reproduces the main model on Granite embeddings:
  - Exact: all KOs 17.60, C 53.80, D 13.96, A 13.71, Dp 16.36.
  - Close: B 10.33 (main 10.32), De 18.89 (18.86). For those two targets the sample composition of some tiles
    differs, which shifts technical and within-tile terms by ≤ 0.0008.
- **Granite tile match:** day gap median 0.59 (75th percentile 1.31, max 9.5); window offset median 2.83 km (max
  19.48); valid fraction ≥ 0.80.
- **Samples outside the benchmark (574):**
  - I09N-2016 185, HLY2202-2022 100, I07N-2018 63, P18-2016 36, P02-2022 31, A13.5-2024 26, A16N-2023 25, I05-2023
    15, I09N-2025 13, C13.5/A13.5 11, I08S-2024 5, AMT28 16.
  - Pre-2016 cruises: NH1418 23, AE1319 13, BVAL46 12.
- **P04E-2026:** 279 MG samples in FAIRe, no KO counts yet.

## R65. Hurdle models for patchy KOs (nifH, amoA, nasA, phnJ): rarefied presence + abundance where present (2026-11-24; full run done)

**Why.** The 90% prevalence filter of the main targets drops nifH, amoA, nasA and phnJ. They are modelled separately as
a hurdle: presence / absence (logistic) plus abundance where present, at a relaxed prevalence cut.

**Absence is partly technical** (R45; the coverage-50 filter turns low coverage into zeros). Raw presence (any gene
> 0), by sample-depth quartile, lowest to highest:
- nifH 0.17 / 0.29 / 0.40 / 0.62; amoC 0.16 / 0.15 / 0.18 / 0.48; nasA 0.57 / 0.73 / 0.87 / 0.95; phnJ 0.61 / 0.84 /
  0.92 / 0.97.
- By platform, nifH is 0.08 on HiSeq 4000, 0.33 on NovaSeq 6000 and 0.42 on NovaSeq X.

**Method** (`scripts/hurdle_kos.py`; `data/results/hurdle/`)
- **Counts.** Per sample, KO gene counts summed over libraries (`sample_counts.pkl`). These come from
  `ko_counts_cov50_all.tsv` with the same library → sample link and filters as `build_ko_targets.py`.
- **Rarefied presence.** A KO is present if count / total ≥ 1 / G, with G = 1,000,000 genes; that is, if it would be
  expected at least once in a library of G genes. Samples with fewer than G genes are dropped: 180 of 1,990, leaving
  1,111 of the 1,215 benchmark tiles. Tile presence = majority of the tile's kept samples.
- **Abundance where present.** log(share) minus the mean log share of the 3,933 core KOs, which is the CLR scale of the
  main targets; averaged over present samples.
- **Hurdle set.** KOs present in 20–90% of tiles, plus the focus KOs: 5,808 KOs, 212 of them also in the core set.
- **Presence models.** L2 logistic, each cruise held out, penalty chosen by inner leave-cruise-out.
  - Technical only: platform group + spline of log10 total genes.
  - Granite PCs: the number of PCs (4 / 8 / 16 / 32) is chosen with the penalty inside the training cruises. With all
    32, sparse presence overfits: nifH pooled AUC 0.541, against 0.654 with 4 PCs.
  - Granite + technical.
- **Presence scores.**
  - Pooled AUC: dominated by cruise-level differences.
  - **Within-cruise AUC:** n-weighted mean over cruises, the "pattern" measure, used as the primary score.
  - Deviance explained.
- **The "presence is technical" check:** technical-only AUC within cruise ≥ Granite − 0.02.
- **Abundance models:** technical removal in-sample, ridge on Granite PCs, leave-cruise-out penalty (as the main
  model). Fitted for KOs with ≥ 50 present tiles.

**Focus genes** (`focus_models.csv`, `focus_predictions.npz`). Presence AUC. Columns: tile prevalence; within cruise
for Granite, technical only, in-situ temperature (spline) and raw satellite (OLCI + SST); then pooled for Granite and
technical only.

| Gene | Prevalence | Within: Granite | Within: technical | Within: temperature | Within: OLCI + SST | Pooled: Granite | Pooled: technical |
|---|---:|---:|---:|---:|---:|---:|---:|
| nifH | 0.28 | 0.629 | 0.534 | 0.625 | 0.613 | 0.655 | 0.409 |
| nifD | 0.15 | 0.696 | 0.588 | 0.705 | 0.704 | 0.721 | 0.339 |
| nifK | 0.32 | 0.484 | 0.565 | 0.513 | 0.493 | 0.376 | 0.622 |
| amoA / pmoA | 0.16 | 0.665 | 0.484 | 0.606 | 0.736 | 0.590 | 0.224 |
| amoB / pmoB | 0.15 | 0.691 | 0.548 | 0.632 | 0.669 | 0.677 | 0.251 |
| amoC / pmoC | 0.17 | 0.653 | 0.476 | 0.592 | 0.681 | 0.555 | 0.252 |
| nasA | 0.77 | 0.584 | 0.576 | 0.659 | 0.583 | 0.610 | 0.535 |
| phnJ | 0.88 | 0.683 | 0.759 | 0.641 | 0.599 | 0.647 | 0.725 |
| phnI | 0.85 | 0.547 | 0.740 | 0.586 | 0.598 | 0.584 | 0.703 |
| phnH | 0.86 | 0.643 | 0.646 | 0.689 | 0.688 | 0.626 | 0.582 |
| phnG | 0.86 | 0.664 | 0.647 | 0.677 | 0.661 | 0.694 | 0.513 |
| phnL | 0.85 | 0.607 | 0.676 | 0.615 | 0.560 | 0.593 | 0.670 |
| phnK | 0.85 | 0.580 | 0.664 | 0.679 | 0.632 | 0.563 | 0.669 |
| phnM | 0.95 | 0.664 | 0.747 | 0.512 | 0.664 | 0.598 | 0.713 |

- **Nitrogen fixation (nifH, nifD): a real but modest signal.**
  - Within-cruise AUC is 0.63 / 0.70 with Granite, against 0.53 / 0.59 for technical only, and it matches temperature
    (0.63 / 0.71).
  - nifH is present in 41% of tiles above 24 °C and 14% below; the within-cruise correlation with temperature is 0.24.
  - Cruise-level differences dominate: I09N 2016 at 30.6 °C has 0.14, I09N-2025 has 0.52.
  - nifK does not behave like nifH / nifD. Its presence follows sequencing depth: highest on the two
    deepest-sequenced cruises (P02-2022 0.78, I09N-2025 0.64), and the technical model beats Granite.
- **amoA / B / C (K10944–6; KEGG does not separate pmo and amo): a satellite signal** (Granite 0.65–0.69 against
  technical 0.48–0.55). Presence is highest on I05-2023 (0.70; cold, southern Indian Ocean in austral winter) and
  falls with temperature.
- **phn (C-P lyase) presence is technical:** the technical model alone matches or beats Granite for phnJ, phnI, phnL,
  phnK and phnM.
- **nasA:** temperature 0.66, Granite 0.58, technical 0.58.
- **Abundance where present: not predicted** for any focus gene (R² −0.13 to 0.02).
- **Technical AUC below 0.5 pooled** (nifH 0.41, amo 0.22–0.25): the depth / platform relation reverses between
  cruises, so the technical covariates do not transfer across cruises. Adding them to Granite lowers some pooled AUCs
  (nifD 0.72 → 0.60).

**Full run (5,808 KOs).**
- A first run was lost when the machine rebooted (about 86% done; results were only written at the end).
- `fit` now appends each KO to `fit_shard<i>.csv.partial` and resumes from it.
- To restart or resume: `bash cluster/run_hurdle_local.sh` (4 shards, `caffeinate`). Then
  `scripts/hurdle_kos.py --parts summary figure` for `hurdle_all.csv`, `summary.txt` and `figures/output/hurdle_kos.png`
  (a: Granite vs technical-only within-cruise AUC for all hurdle KOs; b: focus genes, four inputs; c, d: nifH and amoB
  presence against latitude).
- Restarted 2026-11-24 09:25, about 5.5 h.
- **Moved to Morrill (2026-11-24).** The local restart was stopped after a few minutes; its partial files are in
  `data/results/hurdle/local_4shard_aborted/`. Steps:
  1. Copy `data/results/hurdle_cluster_inputs.tar.gz` (12 MB) to the Morrill repo, `git pull`, and unpack it in the
     repo root.
  2. `sbatch --account=191001-364393 cluster/hurdle_fit.sbatch` (24-shard array, 2 CPUs each, resumable).
  3. Copy `data/results/hurdle/fit_shard*.csv` back, then run `scripts/hurdle_kos.py --parts summary figure` locally.
  4. **Bug found while waiting (2026-11-24 09:40).** A smoke test of `summary` on the aborted local partials failed:
     `fit` wrote the `.partial` header from the first KO, then appended each later row with only the blocks that ran for
     it. Unrarefied presence columns are written only for focus KOs and abundance columns only with ≥ 50 present tiles,
     so rows have 13 / 16 / 17 / 20 fields. The final read-back in `fit` therefore raises a pandas `ParserError`, so the
     Morrill shards (submitted with the old code) should finish all their KOs and then exit FAILED without writing
     `fit_shard<i>.csv`. No results are lost: each row's column order is fixed and the field counts are distinct.
     - Fix (`scripts/hurdle_kos.py`): `read_partial()` maps each row by its field count; new rows are always written with
       all 20 columns (`FIT_COLS`); `summary` reads `fit_shard<i>.csv` or, failing that, `fit_shard<i>.csv.partial`, and
       reports how many of the 5,808 KOs it found (written to `summary.txt`).
     - Checked: the summary and figure run on the 150 KOs of the aborted local partials (one 20-field row among the
       16-field rows parsed correctly). A 3-KO fit in scratch wrote 20-field rows and resumed correctly; nifH reproduces
       the focus table (within-cruise AUC Granite 0.629, technical 0.534).
     - **Step 3 changes to:** copy `data/results/hurdle/fit_shard*.csv.partial` back (24 files), then
       `scripts/hurdle_kos.py --parts summary figure`.

### R65 full-run results (5,808 KOs; 2026-11-24)

**Run.** Morrill array `bgs_hurdle`, 24 shards. As expected (step 4 above), the shards ended on the old final read-back,
so the 24 `fit_shard<i>.csv.partial` files were copied back. Rows: 14 with 20 fields (the 14 focus KOs) and 5,794 with 16.
`summary` read 5,808 of the 5,808 expected KOs: none missing, no duplicates. Outputs: `data/results/hurdle/hurdle_all.csv`,
`summary.txt`, `figures/output/hurdle_kos.{png,pdf}`. The full-run focus-gene values match the focus table above exactly
(e.g. nifH within-cruise AUC 0.629 Granite / 0.534 technical; amoB 0.690 / 0.548).

**Setup recap.** 1,111 tiles in 13 cruises; each cruise held out in turn. Presence classes (from `summary`):
- **Satellite:** within-cruise Granite AUC ≥ 0.60 *and* ≥ 0.05 above technical only.
- **Technical:** technical-only within-cruise AUC ≥ Granite − 0.02.
- **Weak:** neither.

**Presence: score distributions** (median, 10–90% range over 5,808 KOs)

| Score | Granite PCs | Technical only | Granite + technical |
|---|---|---|---|
| Within-cruise AUC | 0.626 (0.543–0.704) | 0.550 (0.474–0.663) | 0.639 (0.562–0.730) |
| Pooled AUC | 0.660 (0.465–0.780) | 0.512 (0.320–0.659) | 0.709 (0.552–0.808) |

- **Classes:** satellite 3,047 (52.5%), technical 1,617 (27.8%), weak 1,144 (19.7%).
- Within-cruise Granite AUC is below 0.5 for 182 KOs; technical-only AUC is below 0.5 for 1,183.
- **Mirror-tail check against chance.** If presence had no relation to Granite, within-cruise AUC would scatter
  roughly symmetrically around 0.5. Only 1 KO has Granite AUC < 0.40, against 3,853 above 0.60 (35 below 0.45, against
  5,109 above 0.55). So the satellite class is not a multiple-testing artefact. It is still a descriptive threshold, not
  a per-KO test.
- **Adding the technical covariates to Granite** adds little within cruise: median change +0.004; +0.02 or more for
  1,853 KOs.
- **Pooled AUC exceeds within-cruise AUC** for Granite by a median 0.031, because cruise-level differences help the
  pooled score.
- **Granite PCs chosen** (median over folds): 4 for 905 KOs, 8 for 2,047, 16 for 2,110 and 32 for 746. The PC count
  matters for sparse targets (R65 method).
- **Near-duplicate KOs.** There are 5,799 distinct tile-presence profiles among the 5,808 KOs (3,040 of the 3,047
  satellite KOs). A few KOs are counted from the same genes, e.g. ETV1 / ETV4 / ETV5 have identical scores.

**Class by tile prevalence**

| Prevalence | Satellite | Technical | Weak | Total |
|---|---:|---:|---:|---:|
| 0–0.3 | 601 | 283 | 259 | 1,143 |
| 0.3–0.5 | 1,160 | 337 | 320 | 1,817 |
| 0.5–0.7 | 816 | 371 | 271 | 1,458 |
| 0.7–0.9 | 470 | 625 | 294 | 1,389 |
| > 0.9 | 0 | 1 | 0 | 1 |

- **Near-ubiquitous KOs (prevalence 0.7–0.9) are mostly technical (45%).** When a KO is almost always present, its
  occasional absences are dropouts in shallow or HiSeq libraries. This is the phn pattern.
- **Mid-prevalence KOs (0.3–0.5) are mostly satellite (64%).**
- **The 212 KOs also in the core set:** satellite 93, technical 46, weak 73.

**Abundance where present** (ridge on Granite PCs after in-sample technical removal; all 5,808 KOs had ≥ 50 present
tiles)
- Overall: median R² 0.087 (10–90% −0.058 to 0.402). R² > 0.1 for 2,796 KOs; > 0.3 for 1,358; > 0.5 for 38.
- **By presence class** (median R², 10–90% range):
  - Satellite: 0.206 (−0.037 to 0.427).
  - Weak: 0.070 (−0.061 to 0.396).
  - Technical: −0.012 (−0.076 to 0.171).
- Abundance R² tracks presence AUC within cruise (Spearman 0.23).
- Median share of abundance variance that is technical (in-sample): 0.093.
- **Focus genes:** abundance where present is not predicted (R² −0.134 to 0.024). This matches the focus run.

**Focus genes in the full-run classes**
- **Satellite:** nifH, nifD, amoA, amoB, amoC.
- **Technical:** nifK; phnG, phnH, phnI, phnJ, phnK, phnL, phnM; and nasA.
  - **nasA** is technical only by the class rule (Granite 0.584 is below 0.60; technical 0.576). In-situ temperature
    does best for nasA (0.659, focus table).
  - **phnG / phnH:** Granite and technical are tied (0.664 / 0.647 and 0.643 / 0.646), so these are "technical" by the
    −0.02 rule rather than because technical wins.

**Rarefying removes part of the depth signal** (focus KOs; technical-only within-cruise AUC, unrarefied vs rarefied
presence):
- nifH 0.651 vs 0.534; nifD 0.673 vs 0.588.
- amoA 0.570 vs 0.484; amoB 0.575 vs 0.548; amoC 0.570 vs 0.476.
- nasA 0.650 vs 0.576; phnH 0.768 vs 0.646; phnJ 0.788 vs 0.759.
- Unrarefied prevalence is also higher (nifH 0.44 vs 0.28).
- Rarefying to G = 10⁶ genes therefore does what it was meant to, but not completely: technical AUC stays above 0.5
  within cruise.

**Top satellite KOs** (`summary.txt`; within-cruise Granite AUC 0.82–0.94).
- Prokaryotic examples:
  - ykkD multidrug pump 0.941 (prevalence 0.27).
  - **mcrA** (methyl-coenzyme M reductase α) 0.905 (prevalence 0.24).
  - nrsD 0.890; malF 0.885; aaeB 0.867.
  - The rifamycin PKS modules rifA / rifB / rifE 0.86–0.88, and asm21 0.872.
  - abrB 0.854; pcpB 0.846; csy2 0.838.
- Eukaryote-annotated examples:
  - Kv1 potassium channels KCNA2 / 4 / 5, 0.83–0.84.
  - SLC7A14 0.869; PCSK4 0.851.
  - ETV1 / 4 / 5 0.828.
- **Caveat on the top list.** KOs like the Kv channels, rif PKS modules or mcrA in surface seawater are probably
  homolog hits that mark a lineage rather than the named function. Their presence would then follow the
  eukaryote : prokaryote balance and community-composition axes that the main model already predicts (R59–R61), so
  they should not be read as evidence for that pathway in situ. Not checked yet: correlation of each satellite KO's
  presence with the C target (eukaryote : prokaryote log ratio).

**Interpretation for the paper**
- **Presence / absence of patchy KOs is predictable from space within cruise for about half of them.** Granite is above
  technical only by ≥ 0.05 AUC, at a median within-cruise AUC of 0.63, with essentially no chance-level tail below 0.4.
- **Near-ubiquitous KOs are mostly technical:** their absences are sequencing-depth dropouts. The phn genes are the
  clearest case.
- **Where a KO is present, its abundance is modestly predictable for satellite-class KOs** (median R² 0.21) and not at
  all for technical-class KOs.
- **Focus genes.** nifH / nifD and amoA–C are satellite-class. Their signal matches what in-situ temperature gives
  (nif) or raw OLCI + SST (amo). Their abundance where present is not predicted.

**Figure** (`figures/output/hurdle_kos.{png,pdf}`; `scripts/hurdle_kos.py --parts figure`). Layout fixes made today:
the amoB / nifD labels in panel a no longer overlap, and the panel b legend now sits above the axis instead of over the
phn points.

Draft legend: **Figure — Presence of patchy genes.** **a**, Presence AUC within cruise (n-weighted mean over 13
held-out cruises) for 5,808 KEGG orthologues present in 20–90% of 1,111 tiles. Presence is rarefied to 10⁶ genes, and a
tile is present if most of its samples are. x: logistic model on sequencing platform + spline of log10 sequencing depth.
y: logistic model on Granite principal components (4–32, chosen within the training cruises). Blue: satellite signal
(Granite ≥ 0.60 and ≥ 0.05 above depth + platform; 3,047). Orange: technical (depth + platform within 0.02 of Granite
or better; 1,617). Grey: neither (1,144). Circles: focus genes. **b**, Within-cruise presence AUC of focus genes for four
inputs: depth + platform, in-situ temperature (spline), OLCI band means + SST, and Granite PCs. **c**, **d**, nifH and
amoB: observed share of tiles with the gene present (black) and mean held-out predicted probability, in 10° latitude
bins with ≥ 8 tiles.

---

## R66. Why hurdle KOs are "satellite" or "technical": the split is mostly eukaryotic vs prokaryotic KOs, not detection limit (2026-11-24)

**Questions.**
1. Does adding Granite to the technical covariates (depth + platform) improve presence AUC?
2. **Hypothesis to test:** technical-class KOs are low-abundance genes near the detection limit. Their absences are
   sampling dropouts, and when present their counts are small and dominated by counting noise, so abundance has nothing
   to predict. Satellite-class KOs are mid-prevalence genes whose presence and abundance follow real biogeography. Long
   genes should be more depth-sensitive; this would fit nifK being technical while nifH / nifD are satellite.

**Method** (`scripts/hurdle_mechanism.py`; outputs `data/results/hurdle/mechanism.csv`, `mechanism_summary.txt`,
`figures/output/hurdle_mechanism.{png,pdf}`). Same 5,808 KOs, 1,111 tiles, 13 cruises and samples as R65 (samples with
≥ 10⁶ genes). Per KO:
- **Expected count at G:** mean share over kept samples × 10⁶.
- **Counting-noise share (`noise_frac`)** of the abundance-where-present variance.
  - For a present sample with count c, var(log c) ≈ 1/c (Poisson).
  - The tile value is a mean over its n present samples, so its noise variance is Σ(1/cᵢ) / n².
  - `noise_frac` = mean noise variance over present tiles / observed variance of the tile values.
  - 1 − `noise_frac` is an approximate ceiling on abundance R².
  - Counts are gene counts with the 1/n multi-KO split (not reads), so Poisson is an approximation.
- **Within-cruise correlations of tile presence** (cruise means removed from both) with log10 depth (`r_depth_within`),
  in-situ temperature (`r_temp_within`) and observed target C, the eukaryote : prokaryote log ratio (`r_C_within`).
- **Predicted-C presence AUC** (`auc_within_Cpred`): within-cruise AUC of the Granite prediction of C for that tile
  (main model R60, cruise held out). The sign is taken from the training cruises, so this is a one-number model
  comparable with Granite's own within-cruise AUC.
- **Gene length and eukaryote share:** mean reference gene length (nt) and eukaryote share of the KO's reference genes
  (genes_Eukaryota / genes), both from R44 `transfer/results/ko_gc/ko_gc.csv` (a 1-in-10 sample of the 2026 reference
  database). No KO is missing.

### 1. Granite added to technical: a large gain except for technical-class KOs

Within-cruise presence AUC, Granite + technical minus technical only:
- **All KOs:** median +0.081 (10–90% −0.002 to +0.169). At least +0.02 for 4,716 of 5,808 KOs; −0.02 or worse for 280.
- **By class:**
  - Satellite: +0.121 (0.532 → 0.653); ≥ +0.02 for 100%.
  - Weak: +0.062 (0.522 → 0.587); ≥ +0.02 for 93%.
  - Technical: +0.010 (0.625 → 0.640); ≥ +0.02 for 38%.
- **Pooled AUC:** median +0.163 (10–90% +0.023 to +0.400). By class: satellite +0.242 (0.483 → 0.743), weak +0.125
  (0.508 → 0.664), technical +0.050 (0.583 → 0.666).
- **The reverse direction (technical added to Granite) depends on class.** The R65 figure of +0.004 is the median
  over *all* KOs. By class, Granite + technical minus Granite:
  - Satellite: −0.003 (Granite 0.656, both 0.653).
  - Weak: +0.002 (0.588 vs 0.587).
  - **Technical: +0.060** (Granite 0.571, technical 0.625, both 0.640).
- Note on medians: the class medians of the three AUCs (technical / Granite / both) are 0.532 / 0.656 / 0.653 for
  satellite, 0.625 / 0.571 / 0.640 for technical and 0.522 / 0.588 / 0.587 for weak. The "both − X" figures above are
  medians of per-KO differences, so they need not equal differences of medians.
- So each input adds to the other where the other is weak: Granite adds to technical for satellite KOs, and technical
  adds to Granite for technical KOs. Neither carries everything.

### 2. Per-class medians

| | Satellite | Weak | Technical |
|---|---:|---:|---:|
| Prevalence | 0.455 | 0.495 | 0.612 |
| Expected count in 10⁶ genes | 4.69 | 3.45 | 4.18 |
| Median count where present | 20.0 | 16.0 | 21.5 |
| Counting-noise share | 0.085 | 0.125 | 0.117 |
| Abundance R² | 0.206 | 0.070 | −0.012 |
| r(presence, depth) within cruise | −0.071 | −0.027 | +0.171 |
| r(presence, temperature) within cruise | −0.228 | −0.142 | +0.045 |
| Between-cruise share of presence variance | 0.208 | 0.162 | 0.153 |
| Mean gene length (nt) | 741 | 756 | 874 |
| **Eukaryote share of genes** | **0.907** | 0.895 | **0.041** |

Within prevalence bins (0.2–0.5 / 0.5–0.7 / 0.7–0.9; n satellite 1,757 / 816 / 470, weak 579 / 271 / 294, technical
620 / 371 / 625):
- **Expected count:** satellite 3.00 / 6.90 / 12.18; weak 1.85 / 4.79 / 10.92; technical 1.66 / 4.03 / 8.82.
- **Counting-noise share:** satellite 0.107 / 0.074 / 0.059; weak 0.175 / 0.116 / 0.066; technical 0.176 / 0.125 /
  0.083.
- **Abundance R²:** satellite 0.118 / 0.336 / 0.315; weak −0.015 / 0.164 / 0.274; technical −0.023 / −0.011 / 0.010.
- **r(presence, depth):** satellite −0.091 / −0.065 / +0.009; technical +0.079 / +0.197 / +0.214.
- **Gene length (nt):** satellite 750 / 739 / 722; technical 862 / 904 / 862.
- **Eukaryote share:** satellite 0.892 / 0.926 / 0.918; weak 0.863 / 0.912 / 0.922; technical 0.122 / 0.028 / 0.032.

### 3. The detection-limit hypothesis: half right for presence, wrong for abundance

**Presence: depth-driven for technical KOs, and long genes are more depth-sensitive (supported, for prokaryotic KOs).**
- Technical-class presence rises with sequencing depth within cruise (r +0.17); satellite-class presence does not
  (−0.07).
- Within prevalence bins, technical KOs have somewhat lower expected counts than satellite KOs (1.7 vs 3.0 at 0.2–0.5).
- Gene length vs depth sensitivity, across all KOs: Spearman ρ +0.336; vs technical-only within-cruise AUC +0.270.
  Within prevalence bins: +0.339 / +0.356 / +0.434.
- **The length effect is a prokaryotic-KO effect.**
  - Prokaryotic KOs (euk share ≤ 0.1; 1,716): ρ(length, r_depth) = +0.404, and +0.432 / +0.457 / +0.480 within the
    prevalence bins; ρ(length, technical AUC) +0.366.
  - Eukaryotic KOs (> 0.9; 2,408): +0.093 and −0.017.
  - This fits R15: the coverage-50 filter keeps 68% of counts for genes under 140 aa but 13% over 565 aa. A long
    prokaryotic gene needs a better assembly to pass, and assemblies are worse in shallow libraries.
  - Across all KOs, the length–depth curve is not monotonic (it dips at 600–900 nt). That is because mid-length KOs are
    mostly eukaryotic, and eukaryotic presence correlates *negatively* with depth (median −0.087). Panel d therefore
    splits by domain.
- **Why eukaryotic presence falls with depth:** eukaryote-rich tiles were sequenced less deeply within a cruise. The
  within-cruise r between C (eukaryote : prokaryote log ratio) and log10 depth is −0.231. This is biology carried by
  depth, not a depth artefact.
- **But "near the detection limit" in the Poisson sense does not describe the technical class.**
  - Prokaryotic technical KOs have *higher* expected counts than prokaryotic satellite KOs (4.74 vs 2.89) and lower
    counting noise (0.103 vs 0.119).
  - The phn genes (all technical) have expected counts of 16–62 per 10⁶ genes, yet are absent in 5–15% of tiles. A
    Poisson dropout at 16–62 expected counts is essentially impossible.
  - Their absences are therefore assembly / coverage-filter failures in shallow libraries, which hit long genes most.
    This is not counting noise.

**Abundance: counting noise caps R² in general, but does not explain the technical class.**
- Across KOs, noise share vs abundance R²: Spearman −0.465. Noise share quantiles: 10% 0.048, 25% 0.066, 50% 0.100,
  75% 0.165, 90% 0.277. Only 9 KOs exceed their noise ceiling (1 − noise share).
- Median abundance R² by noise-share quintile:

| Noise-share quintile | Satellite | Weak | Technical |
|---|---:|---:|---:|
| 0.008–0.060 | 0.350 | 0.331 | −0.014 |
| 0.060–0.083 | 0.316 | 0.280 | 0.003 |
| 0.083–0.119 | 0.221 | 0.173 | −0.006 |
| 0.119–0.189 | 0.080 | 0.045 | −0.009 |
| 0.189–1.894 | −0.017 | −0.019 | −0.015 |

  (n per quintile: satellite 721 / 758 / 641 / 497 / 430; weak 171 / 165 / 203 / 291 / 314; technical 270 / 238 /
  318 / 373 / 418.)
- **Technical KOs with low counting noise are still unpredictable** (R² −0.014 in the lowest-noise quintile, against
  0.350 for satellite). Their abundance has real variance, about 90% of it beyond counting noise, that Granite does not
  predict.

### 4. What does separate the classes: domain

- **The best single separator of satellite vs technical is the eukaryote share** of the KO's genes.
  - Medians: 0.907 satellite, 0.041 technical.
  - Satellite vs technical KOs (n 4,664, of which 3,047 satellite), 10-fold cross-validated logistic on ranked
    predictors (AUC; standardised coefficient):
    - Prevalence: 0.626 (−0.45).
    - Expected count: 0.536 (+0.13).
    - Gene length: 0.657 (−0.57).
    - **Eukaryote share: 0.773 (+1.12).**
    - Expected count + prevalence: 0.745.
    - All four: 0.835 (expected count +1.22, length −0.49, eukaryote share +0.77, prevalence −1.38).
- **Class by eukaryote share:**

| Eukaryote share | Satellite | Weak | Technical | n |
|---|---:|---:|---:|---:|
| ≤ 0.1 | 498 | 263 | 955 | 1,716 |
| 0.1–0.5 | 100 | 50 | 177 | 327 |
| 0.5–0.9 | 856 | 278 | 223 | 1,357 |
| > 0.9 | 1,593 | 553 | 262 | 2,408 |

- The 3,933 core KOs are mostly prokaryotic: median eukaryote share 0.064; 56.7% at ≤ 0.1 and 9.9% at > 0.9. The
  patchy set is the reverse: 41% eukaryotic (> 0.9) and 30% prokaryotic (≤ 0.1).
- **Eukaryotic KOs: the satellite signal is largely the eukaryote : prokaryote balance (target C).**
  - Median within-cruise presence AUC: Granite 0.630, *predicted C alone* 0.617, technical 0.529.
  - Presence correlates with observed C within cruise (median r +0.380).
  - Abundance where present: median R² 0.263. For the 1,593 eukaryotic satellite KOs it is +0.306; weak +0.204;
    technical +0.066.
  - Abundance is on the CLR-like scale relative to the (mostly prokaryotic) core KOs, so a eukaryotic KO's abundance
    also moves with C.
  - Not yet tested: the abundance R² of eukaryotic KOs from predicted C alone. We expect it to be most of the 0.31.
- **Prokaryotic KOs: presence mostly technical, abundance unpredictable.**
  - Median within-cruise AUC: Granite 0.609, technical 0.604, predicted C 0.532. r with C −0.052.
  - Abundance R²: median −0.017. By class: satellite −0.002, weak −0.027, technical −0.020.
  - Even at the lowest counting noise, prokaryotic KOs reach a median abundance R² of only 0.04 (satellite) and −0.02
    (technical); see figure panel c.
  - Prokaryotic KOs by prevalence bin (0.2–0.5 / 0.5–0.7 / 0.7–0.9): satellite 278 / 108 / 108; weak 153 / 60 / 50;
    technical 305 / 252 / 397.
- **Across classes:**
  - Granite minus predicted-C within-cruise AUC, median: satellite +0.022, weak +0.010, technical +0.020. Granite is
    ahead by ≥ 0.05 for 29% / 21% / 32%.
  - Median r_C_within: satellite 0.373, weak 0.263, technical −0.017.
  - Median auc_within_Cpred: satellite 0.625, weak 0.574, technical 0.545.
- **The 498 prokaryotic satellite KOs are the set with a signal beyond C.**
  - Top of the list (within-cruise AUC: Granite / predicted C / technical):
    - ykkD 0.940 / 0.725 / 0.507; mcrA 0.905 / 0.641 / 0.554.
    - nrsD 0.890 / 0.474; malF 0.885 / 0.873; aaeB 0.867 / 0.537.
    - abrB 0.854 / 0.671; pcpB 0.846 / 0.463; bar 0.840; csy2 0.838; aaeR 0.837.
    - K07721 (ArsR family) 0.832 / 0.332; tgtA2 0.815 / 0.229; adiC 0.807 (r_temp −0.503).
  - nifH, nifD and amoA–C are in this set (eukaryote share 0.000–0.010).
  - Several top ones (mcrA, sspH spore protein, ykkD, pcpB) are odd for surface seawater. They could be homolog hits
    marking a lineage, or a technical axis that depth + platform does not capture (e.g. the axis-1 KO-profile spread,
    R45).
  - Their abundance where present is still not predicted (median −0.002).
- **Other correlations** (Spearman, all KOs):
  - Granite within-cruise AUC vs |r(presence, temperature)|: +0.584; vs signed r: −0.310. Most satellite KOs are more
    often present in colder water within a cruise (median −0.228).
  - Eukaryote share vs Granite AUC: +0.105.
  - Expected count vs technical AUC: −0.094; vs r_depth +0.164. Noise share vs expected count: −0.856.

### 5. Gene length and the nif genes: does not explain the split

| Gene | Class | Prevalence | Expected count | Noise share | Abundance R² | r(depth) | r(temperature) | Length (nt) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| nifH | satellite | 0.275 | 1.95 | 0.064 | −0.134 | +0.071 | +0.236 | 836 |
| nifD | satellite | 0.151 | 1.10 | 0.055 | −0.103 | +0.118 | +0.270 | 1,334 |
| nifK | technical | 0.325 | 2.22 | 0.063 | −0.029 | +0.127 | +0.188 | 1,225 |
| amoA | satellite | 0.164 | 2.17 | 0.055 | −0.038 | +0.038 | −0.182 | 698 |
| amoB | satellite | 0.146 | 1.91 | 0.078 | −0.034 | +0.030 | −0.182 | 1,000 |
| amoC | satellite | 0.170 | 1.82 | 0.071 | −0.001 | +0.044 | −0.132 | 669 |
| nasA | technical | 0.767 | 17.88 | 0.034 | +0.014 | +0.105 | −0.159 | 1,818 |
| phnJ | technical | 0.885 | 40.01 | 0.016 | +0.024 | +0.273 | +0.183 | 843 |
| phnI | technical | 0.845 | 41.87 | 0.016 | −0.037 | +0.292 | +0.133 | 1,006 |
| phnH | technical | 0.858 | 16.07 | 0.039 | −0.027 | +0.199 | +0.144 | 612 |
| phnG | technical | 0.859 | 21.87 | 0.031 | +0.022 | +0.180 | +0.183 | 467 |
| phnL | technical | 0.853 | 25.34 | 0.025 | −0.058 | +0.225 | +0.112 | 699 |
| phnK | technical | 0.851 | 28.48 | 0.022 | −0.014 | +0.236 | +0.131 | 748 |
| phnM | technical | 0.947 | 62.18 | 0.015 | −0.117 | +0.202 | +0.111 | 1,033 |

- **nifK (1,225 nt) is shorter than nifD (1,334 nt)**, and their depth correlations are about the same (+0.127 vs
  +0.118). So length does not explain why nifK is technical.
- nifK has twice nifD's prevalence (0.33 vs 0.15), although the two genes sit in one operon. A likely reason
  (unverified) is that K02591 also collects homolog hits (nifE / nifN-like and other nitrogenase-like genes), which dilute
  the diazotroph signal.
- The phn genes are not near the detection limit (expected count 16–62; noise share 0.015–0.039), yet presence follows
  depth (r +0.18 to +0.29). This is the assembly / filter dropout described in section 3.

### Interpretation

- **Presence:** "technical" KOs are mostly *prokaryotic* KOs whose occasional absences track sequencing depth. Long
  genes are affected most, consistent with assembly failing the coverage-50 filter in shallow libraries. This is not
  Poisson dropout at the detection limit: technical KOs are not rarer than satellite KOs, and phn has 16–62 expected
  counts.
- **"Satellite" is two different things** (corrected in section 6; class is not just domain):
  - Eukaryotic satellite KOs: the signal is largely the eukaryote : prokaryote balance that the main model already
    predicts (C, R² 0.54). Granite 0.648 against predicted C 0.633.
  - Prokaryotic satellite KOs (498): the strongest Granite signal of any group, and mostly *not* C. Granite 0.685
    against predicted C 0.532.
  - The 3,047 satellite KOs should therefore not be presented as 3,047 independent results. The eukaryotic part is
    largely one axis; the prokaryotic part is where gene-level biogeography is.
- **Abundance where present:** predictable for eukaryotic KOs (median R² 0.31 in the satellite class) and not for
  prokaryotic KOs of any class (≈ 0), even at low counting noise.
  - This answers why technical-class abundance is unpredictable: technical KOs are prokaryotic. Patchy prokaryotic
    genes' abundance where present is not predicted from space.
  - That contrasts with the core (mostly prokaryotic, always present) KOs at about 16% (R59). The difference is
    probably that core KOs' variance is dominated by shared community axes, while a patchy prokaryotic KO, where
    present, varies with its particular carriers.
- **The genuinely gene-specific satellite signals** are among the 498 prokaryotic satellite KOs, which include nif and
  amo. These need the next checks:
  - (a) Granite beyond predicted C *and* temperature for each.
  - (b) Whether the oddities at the top (mcrA, ykkD, sspH) follow a residual technical axis, such as axis 1 of R45.
  - (c) The abundance R² of eukaryotic KOs from predicted C alone.

**Figure** `figures/output/hurdle_mechanism.{png,pdf}` (`python scripts/hurdle_mechanism.py`, about 2 min;
`--figure-only` redraws).
- **a:** Share of KOs in each presence class by eukaryote share of the KO's reference genes; n at right.
- **b:** Within-cruise presence AUC, Granite PCs vs the Granite-predicted eukaryote : prokaryote ratio alone (sign from
  the training cruises), coloured by eukaryote share. Eukaryotic KOs sit on the 1:1 line; prokaryotic KOs sit above it
  or near 0.5.
- **c:** Median abundance-where-present R² by quintile of counting-noise share, for eukaryotic (> 0.9) and prokaryotic
  (≤ 0.1) KOs in the satellite and technical classes.
- **d:** Within-cruise correlation of presence with log10 sequencing depth against mean reference gene length, by
  domain, with decile medians.

### 6. Correction / clarification: class is not just domain (2026-11-24)

The summary in chat said satellite vs technical is "mostly eukaryotic vs prokaryotic KOs". Eukaryote share is the
strongest single separator (AUC 0.77), but the class is not the domain. Domain: eukaryotic = eukaryote share > 0.9;
prokaryotic ≤ 0.1; mixed in between.

**Class shares by domain:**
- Eukaryotic: satellite 0.662, technical 0.109, weak 0.230.
- Mixed: satellite 0.568, technical 0.238, weak 0.195.
- Prokaryotic: satellite 0.290, technical 0.557, weak 0.153.

**Domain × class** (medians). Within-cruise AUC: Granite and predicted C, then Granite − predicted C and the share of
KOs where Granite is ahead by ≥ 0.05.

| Domain | Class | n | Granite | Predicted C | Granite − C | Ahead by ≥ 0.05 | Abundance R² | r(depth) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Eukaryotic | satellite | 1,593 | 0.648 | 0.633 | +0.014 | 11% | +0.306 | −0.091 |
| Eukaryotic | weak | 553 | 0.588 | 0.584 | +0.002 | 6% | +0.204 | −0.069 |
| Eukaryotic | technical | 262 | 0.574 | 0.575 | −0.003 | 6% | +0.066 | −0.109 |
| Mixed | satellite | 956 | 0.662 | 0.631 | +0.023 | 29% | +0.139 | −0.071 |
| Mixed | weak | 328 | 0.585 | 0.574 | +0.007 | 16% | +0.049 | −0.015 |
| Mixed | technical | 400 | 0.565 | 0.545 | +0.020 | 30% | −0.005 | +0.148 |
| Prokaryotic | satellite | 498 | **0.685** | 0.532 | **+0.152** | **87%** | −0.002 | +0.092 |
| Prokaryotic | weak | 263 | 0.597 | 0.526 | +0.062 | 58% | −0.027 | +0.102 |
| Prokaryotic | technical | 955 | 0.572 | 0.533 | +0.033 | 40% | −0.020 | +0.207 |

- **Eukaryotic KOs:** whatever their class, Granite ≈ predicted C (Granite ahead by ≥ 0.05 for only 6–11%). Their
  presence signal is the eukaryote : prokaryote axis.
  - Abundance R² still differs by class within this domain: satellite 0.31, weak 0.20, technical 0.07. So for these
    KOs, class does carry information on abundance beyond domain.
- **Prokaryotic satellite KOs:** the clearest gene-level satellite signal.
  - Granite has the highest median presence AUC of any group, and 87% beat predicted C by ≥ 0.05.
  - Their abundance where present is still not predicted.
  - nifH, nifD and amoA–C are in this group.
- **Technical KOs** are depth-driven in the prokaryotic and mixed domains (r +0.21, +0.15), but not in the eukaryotic
  domain (−0.11). The 262 eukaryotic technical KOs are simply KOs where C, and hence Granite, is weak (0.57) and depth
  + platform happens to do as well.
- **Caveat on the domain variable.** Eukaryote share comes from who carries the KO in the reference database (R44,
  1-in-10 gene sample). It is not taken from which organisms carry the KO in these samples. A KO that is mostly
  eukaryotic in the reference can be found here mainly in prokaryotes, and vice versa. A per-sample domain split of each
  KO's counts would need the per-taxon KO tables (R50).

---

## R67. Which domain carries each KO: read-based domain share replaces the reference-gene share (2026-11-24)

**Why.** R44's "eukaryote share" of a KO (`transfer/results/ko_gc/ko_gc.csv`) was used as the domain label in R44
(a control variable), in R48 / R51 (the prokaryote-specific / eukaryote-specific / shared classes behind targets A, B, C,
D) and in R66.
- **What it is.** The eukaryotic fraction of a deterministic 1-in-10 sample of the KO-annotated genes in the Bio-GO-SHIP
  2026 reference catalogue (`biogoship_2026_database_v2`).
- **What the catalogue is.** It is not one metagenome. It is the mmseqs-clustered, non-redundant catalogue that our reads
  are mapped to (R30–R32): OMDB v2 NR95 prokaryotic genomes verbatim, Bio-GO-SHIP assemblies (e.g. A16 unbinned contigs)
  and eukaryotic sources, largely transcript assemblies (R48), re-clustered at 90% identity / 90% coverage. Each
  sequence's domain comes from `biogoship_2026_taxonomy.tsv` (EUKulele for eukaryotes).
- **So it counts distinct catalogue genes, not what is in our water.** Eukaryotic catalogue entries are many per gene
  (diverse transcripts) while seawater reads are about 95% prokaryotic. A KO with many eukaryotic catalogue entries can
  still be carried mostly by prokaryotes in our samples.

**Method** (`scripts/ko_domain_reads.py`; `data/results/ko_domain/ko_domain_reads.csv`, `library_domain.npz`,
`summary.txt`; `figures/output/ko_domain_reads.{png,pdf}`).
- Each library's reads per KO, split by the domain of the reference sequence they map to.
  - Source: the R50 phylum-level per-taxon files `ko_by_taxon/<lib>/level_2_phylum_division_ko_counts.tsv.gz`, whose
    sum over taxa is `ko_counts_cov50_all.tsv`.
  - Domains: prokaryotic = Bacteria + Archaea; eukaryotic = Eukaryota; viruses separately.
- Libraries: 2,018 cov50 libraries, the same selection as `build_ko_targets.py` / `domain_split_targets.py`, covering
  all 1,990 hurdle samples.
- Check: the domain split sums back to the official table (max |diff| 0.0074 reads, of 10,537,343,456 KO reads).
  - Over all KO reads: prokaryotic 0.9545, eukaryotic 0.0447, viruses 0.0008, other 0.0000.
- Per KO:
  - Pooled eukaryotic read share.
  - Per-sample share (samples with ≥ 10 reads of the KO): median and IQR.
  - The R48 reference share, Eukaryota / (Eukaryota + Bacteria + Archaea), for comparison.

**Agreement.** Spearman correlation between the reference share and the read share:
- All KOs with reads: +0.894 (n 18,240).
- Core 3,933: +0.839.
- Hurdle 5,808: +0.862.

At the prokaryotic end the reference share is reliable:
- Reference ≤ 0.1 → read share ≤ 0.1 for 97.0% (all KOs), 98.9% (core) and 97.7% (hurdle). Never > 0.9 for core or
  hurdle KOs; 0.5% over all KOs.

At the eukaryotic end it overstates:
- Reference > 0.9 → read share > 0.9 for 81.6% of all KOs, 54.4% of core and 78.1% of hurdle KOs. Read share ≤ 0.1 for
  1.0% / 2.3% / 0.9%.
- Median read share for these KOs: 0.990 (all), 0.923 (core), 0.972 (hurdle).

Medians:
- Reference share: 0.918 over all KOs, 0.072 core, 0.975 hurdle.
- Read share: 0.619 over all KOs, 0.009 core, 0.864 hurdle.
- Pooled over these KOs' reads, the eukaryotic share is 0.045 (all), 0.027 (core) and 0.551 (hurdle).

Cross-tabulation, hurdle KOs (rows: reference share; columns: read share):

| Reference \ reads | ≤ 0.05 | 0.05–0.1 | 0.1–0.5 | 0.5–0.9 | 0.9–0.95 | > 0.95 | All |
|---|---:|---:|---:|---:|---:|---:|---:|
| ≤ 0.05 | 1,410 | 24 | 16 | 0 | 0 | 0 | 1,450 |
| 0.05–0.1 | 181 | 18 | 21 | 2 | 0 | 0 | 222 |
| 0.1–0.5 | 227 | 39 | 56 | 10 | 0 | 2 | 334 |
| 0.5–0.9 | 45 | 23 | 112 | 104 | 16 | 14 | 314 |
| 0.9–0.95 | 10 | 10 | 65 | 85 | 41 | 41 | 252 |
| > 0.95 | 2 | 8 | 90 | 494 | 394 | 2,248 | 3,236 |
| All | 1,875 | 122 | 360 | 695 | 451 | 2,305 | 5,808 |

Core KOs (rows: reference share; columns: read share):

| Reference \ reads | ≤ 0.05 | 0.05–0.1 | 0.1–0.5 | 0.5–0.9 | 0.9–0.95 | > 0.95 | All |
|---|---:|---:|---:|---:|---:|---:|---:|
| ≤ 0.05 | 1,729 | 23 | 10 | 0 | 0 | 0 | 1,762 |
| 0.05–0.1 | 349 | 19 | 13 | 0 | 0 | 0 | 381 |
| 0.1–0.5 | 758 | 103 | 75 | 0 | 0 | 0 | 936 |
| 0.5–0.9 | 117 | 70 | 167 | 14 | 1 | 0 | 369 |
| 0.9–0.95 | 5 | 4 | 36 | 16 | 2 | 2 | 65 |
| > 0.95 | 2 | 0 | 35 | 123 | 63 | 197 | 420 |

**Stability across samples.** For 11,458 KOs with ≥ 100 samples of ≥ 10 reads, the IQR of the per-sample eukaryotic
read share has median 0.000 (90th percentile 0.267). Most KOs are carried by one domain everywhere. The 754 KOs with a
mixed median share (0.1–0.9) vary a lot between samples (IQR median 0.612).

Shares are clipped to [0, 1]: float32 sums put some all-eukaryotic KOs a hair above 1, which first dropped 17 hurdle
KOs (592 over all KOs) from the cross-tabulation's top bin. Medians and percentages were unaffected.

**Spot checks** (reference share → pooled read share; per-sample median):

| KO | Reference | Reads | Note |
|---|---:|---:|---|
| nifH, nifD, nifK, amoA, amoB, mcrA, phnJ, sspH, nrsD, ykkD | 0.00–0.05 | 0.000 | prokaryotic, as expected |
| nasA | 0.038 | 0.010 | |
| TUBA, TUBB, actin, histone H3, H4 | 0.996–1.000 | 0.963–0.989 (median 1.000) | eukaryotic, as expected |
| psbA, psaA | 0.844, 0.824 | 0.087, 0.038 | catalogue has many plastid copies; reads are cyanobacterial |
| rbcL | 0.786 | 0.209 | |
| rpsL, rplB (prokaryote-type ribosomal proteins) | 0.171, 0.184 | 0.239, 0.110 | organellar ribosomes; the plastid leak of R49 |
| **KCNA2, KCNA4** (Kv channels) | 0.989, 0.987 | 0.298, 0.288 (median 0.030, 0.027) | **prokaryotic reads in most samples** |
| **SLC7A14** | 0.982 | 0.116 (median 0.000) | prokaryotic reads |
| **rifA** | 0.809 | 0.013 | prokaryotic reads |
| **PCSK4** | 0.583 | 0.006 | prokaryotic reads |
| **ETV1 / 4 / 5** | 0.707 | 0.210 (median 0.168) | mostly prokaryotic reads |

The bold KOs were described in R65 / R66 as "eukaryote-annotated" top satellite KOs. **By reads they are carried by
prokaryotes in our samples.** Their KO names are eggNOG / KEGG annotations of prokaryotic homologs (e.g. bacterial
voltage-gated K⁺ channels assigned to KCNA KOs). That correction to R65 / R66 is below.

**Consequences for earlier results**
- **R44** (GC vs platform effects, partial correlations controlling for the eukaryote share): the share was only a
  control variable, and the two measures agree (ρ 0.84 on the core KOs). Conclusions unchanged. Not rerun.
- **R48 / R51 classes and targets**, by read share:
  - Prokaryote-specific (reference ≤ 0.05; 1,759 KOs): read share median 0.002, 90th percentile 0.013, max 0.259.
    **Target A is fine.**
  - Eukaryote-specific (reference ≥ 0.95; 420): read share median 0.943, 10th percentile 0.526, minimum 0.038. 260 of
    420 are > 0.9 by reads; 37 are ≤ 0.5. **Target B contains some KOs that are partly prokaryotic in our samples.**
    Rebuilding B with read-based classes is open.
  - Shared (1,700): read share median 0.019 (90th percentile 0.173). Already noted in R51: the shared KOs are mostly
    prokaryotic reads. Dp / De split these by reads and are unaffected.
  - The 54 ribosomal-protein markers: read share median 0.099 (90th percentile 0.280). R49's prokaryote-assigned
    denominator already handles this.
  - **Target C** (log pooled eukaryote-specific / prokaryote-specific KO reads, classes from the reference share) vs
    the direct read-based ratio, log(eukaryotic / prokaryotic KO reads) per sample, averaged per tile (1,215 tiles):
    Pearson r 0.979, Spearman 0.938, within cruise 0.966. **C stands.** The direct ratio is the simpler definition for
    the paper's methods.
- **R66** (domain × class), redone with the read-based domain (eukaryotic > 0.9, prokaryotic ≤ 0.1 pooled read share;
  class labels dropped as agreed).
  - Domain moves, reference → reads:
    - Of 2,408 reference-eukaryotic KOs, 1,992 stay eukaryotic, 405 become mixed and 11 prokaryotic.
    - Of 1,684 reference-mixed KOs, 764 become eukaryotic, 611 stay mixed and 309 become prokaryotic.
    - Of 1,716 reference-prokaryotic KOs, 1,677 stay prokaryotic and 39 become mixed.
    - By reads: 2,756 eukaryotic, 1,055 mixed, 1,997 prokaryotic.
  - Medians by read domain:

| Read domain | n | Granite | Granite + technical | Technical | Predicted C alone | Granite − C | Granite ahead of C by ≥ 0.05 | Abundance R² | r(depth) | r(C) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Eukaryotic | 2,756 | 0.635 | 0.630 | 0.534 | 0.622 | +0.009 | 7% | +0.253 | −0.094 | +0.402 |
| Mixed | 1,055 | 0.618 | 0.626 | 0.527 | 0.586 | +0.022 | 30% | +0.083 | +0.011 | +0.251 |
| Prokaryotic | 1,997 | 0.610 | 0.670 | 0.607 | 0.533 | +0.067 | 56% | −0.017 | +0.158 | −0.049 |

  (Within-cruise presence AUC unless stated.)
  - **The pattern of R66 holds with read-based domain:**
    - Eukaryotic-KO presence is what predicted C gives (Granite ahead by ≥ 0.05 for only 7%).
    - Prokaryotic-KO presence depends on depth (r +0.16) but is predicted beyond C (Granite + technical 0.670 vs
      technical 0.607; Granite ahead of C by ≥ 0.05 for 56%).
    - Abundance where present: R² 0.25 for eukaryotic KOs, −0.02 for prokaryotic ones.
  - **But the R65 / R66 lists of "eukaryote-annotated" top KOs were wrong.** Of the 30 KOs with the highest Granite
    within-cruise AUC, 22 are prokaryotic and 8 mixed by reads; none is eukaryotic. The R65 caveat that these KOs "follow
    the eukaryote : prokaryote balance" does not apply to them.

**Wording correction (R65 / R66).** The counts are **mapped reads per KO** (R30: samtools coverage numreads,
coverage-50 filter, split 1/n over multi-KO references), not gene counts as R15 assumed. Where R65 / R66 say
"genes":
- "10⁶ genes" (the rarefaction depth G) means 10⁶ KO-assigned reads.
- "expected count in 10⁶ genes" means expected reads at that depth.
- "Samples with fewer than G genes" means fewer than 10⁶ KO reads.

The Poisson counting-noise approximation of R66 is, if anything, more apt for reads. The 1/n split makes counts
fractional. The coverage-50 filter makes dropout depend on reference length (R30), which is the mechanism R66 section 3
inferred.

**Going forward.**
- Use the read-based domain share (`ko_domain_reads.csv`, `euk_read_share`) wherever a KO's domain is needed. Report
  the reference share only as the R48 / R51 history.
- Open:
  - Rebuild target B with read-based classes, and possibly define C directly as the read ratio.
  - Then the patchy-gene model restructure (base model Granite + technical ± temperature; panel b with read-based
    colouring), on hold until this is agreed.

**Figure** (supplementary) `figures/output/ko_domain_reads.{png,pdf}`.
- Each point is a KO (18,240 with reads and a reference share; blue: the 5,808 patchy KOs).
- x: eukaryotic share of its reference-catalogue genes. y: eukaryotic share of its reads in our 2,018 libraries.
- Circled: nifH, rpsL, rbcL, psbA, rifA, SLC7A14, KCNA2 and TUBA.

---

## R68. Target C redefined as the direct read ratio; what depends on C, and on B (2026-11-24)

**Decision (user, 2026-11-24).**
- **C is now log(eukaryote-assigned / prokaryote-assigned KO reads)** over all KOs. Each read takes the domain of the
  reference sequence it maps to (R67).
- The R51 definition is kept as `data/targets/domain_split/C_refclass.csv` (identical to the old `C.csv`). That
  definition was log(pooled eukaryote-specific / pooled prokaryote-specific KO reads), with the classes taken from the
  R44 reference-gene share.
- **Target B is not used in the paper.** It is not rebuilt with read-based classes.

**Paper figures** (`figures/output/`, not `old/`) and their targets:
- `c_prediction` (+ `_split`): C. Remade.
- `variance_budget`: the all-KO CLR. C is used only to describe axis 1. Rerun.
- `ko_heatmap_inputs_main`, `module_predictability`: the all-KO CLR only.
- `hurdle_kos`, `hurdle_mechanism`, `ko_domain_reads`: R65–R67. `hurdle_mechanism` uses predicted C. Rerun; that
  figure is being redesigned anyway.
- **No paper figure uses target B.** It is only fitted by `main_model.py` and exported in the benchmark bundle.

**Construction** (`scripts/domain_split_targets.py`).
- Per library: eukaryote-assigned and prokaryote-assigned (Bacteria + Archaea) reads over all KOs, from the
  per-library domain split of `scripts/ko_domain_reads.py` (`data/results/ko_domain/library_domain.npz`). Each is taken
  as a share of the library's KO reads.
- Libraries of one sample are averaged, then C = log(eukaryote share / prokaryote share).
- Eukaryote-assigned share of KO reads: median over samples 0.0196 (min 0.00334), so no pseudocount is needed.
- r with the R51 version: 0.976 per sample; 0.979 per tile (R67).
- Rebuilding all domain-split targets left A, B, D, Dp, De, `ko_classes.csv` and `shared_split.csv` byte-identical.
  Only C changed.

**Main model refit for C only.**
- `scripts/main_model.py` has a new `--targets` option: fit / partition / cruise for the listed targets only. Their rows
  replace those in `scores.csv`, `partition.csv`, `within_tile_k.csv` and `cruise_decomposition.csv`.
- Ran `--parts fit partition cruise figure --targets C`. The other targets' rows are identical to before.

Explained share of the non-technical variance in C (cruise held out):

| Model | Old C | New C |
|---|---:|---:|
| **Granite PCs (main; leave-cruise-out penalty)** | **53.80%** | **53.79%** |
| Granite PCs, GCV penalty | 53.88% | 54.01% |
| OLCI + SST | 50.32% | 50.27% |
| OLCI + SST + 30-day mixed-layer depth | 50.63% | 50.84% |
| In-situ temperature | 45.86% | 46.83% |
| Granite + temperature | 55.30% | 55.20% |
| Time & place | 17.17% | 12.94% |

- Calibration slope (Granite, pooled): 0.955 → 0.955.

Variance partition of C (share of total variance):

| | Technical | Explained | Cruise offsets | Within-tile | Within-cruise |
|---|---:|---:|---:|---:|---:|
| Old C | 15.95% | 45.22% | 5.30% | 6.49% | 26.96% |
| New C | 15.29% | 45.56% | 4.57% | 6.40% | 28.09% |

- Of the non-technical variance, 53.80% → 53.79%. Of the reproducible non-technical variance, 58.30% → 58.18%.

Per-cruise R² of Granite (`c_per_cruise_r2.csv`; old → new):

| Cruise | Old | New |
|---|---:|---:|
| I07N | −1.074 | −0.573 |
| P02-2022 | −0.040 | −0.138 |
| I09N-2025 | 0.251 | 0.254 |
| GOMECC4 | 0.460 | 0.380 |
| A16N-2023 | 0.474 | 0.450 |
| C13.5/A13.5 | 0.474 | 0.356 |
| AMT28 | 0.557 | 0.550 |
| I05-2023 | 0.598 | 0.539 |
| P18 | 0.641 | 0.708 |
| A13.5-2024 | 0.647 | 0.720 |
| HLY2202-2022 | 0.818 | 0.783 |
| I09N (2016) | 0.820 | 0.712 |
| I08S-2024 | 0.865 | 0.871 |

- OLCI + SST, new: I07N −0.474, P02-2022 −0.400, I09N-2025 −0.128, C13.5 0.343, GOMECC4 0.335, A16N 0.487, I05 0.554,
  AMT28 0.613, P18 0.702, I09N 0.340, A13.5 0.667, HLY2202 0.704, I08S 0.899.
- Time & place, new: I07N −1.315 (the only value below −1), P02 −0.795, GOMECC4 −0.617, C13.5 −0.200, I09N-2025
  0.093, P18 0.169, A16N 0.246, HLY2202 0.304, A13.5 0.331, AMT28 0.352, I05 0.368, I09N 0.439, I08S 0.674.
- **The pooled result is unchanged; individual cruises move by up to ±0.1, and I07N by +0.5.**
  - I07N's within-cruise SD is small (0.339). Its R² is sensitive to how much of C's cruise-level offset the definition
    carries.
  - The same cruises remain the poorly predicted ones (I07N, P02-2022, I09N-2025), as in R61 / R62.
- Per-cruise level / spread / pattern for C (`cruise_decomposition.csv`; median within-cruise r, Granite):
  - I09N −0.078, I07N −0.053, C13.5 0.499, GOMECC4 −0.272 (level share 0.509), I09N-2025 0.368, P02 0.206.
  - AMT28 0.825, I05 0.475, I08S 0.917, A16N 0.712, A13.5 0.871, P18 0.816, HLY2202 0.424 (spread share 0.602).
  - GOMECC4 is still the cruise placed wrong in level (bias −0.633).

**Variance budget** (`scripts/variance_budget_main.py`).
- The budget numbers are unchanged; it is built on the all-KO CLR.
- Axis 1's within-cruise correlation with C: −0.602, against −0.655 with the R51 C (R63 text).
- The old C, built from the reference-class KOs, shared part of the KO-profile signal with axis 1. The direct read
  ratio is a cleaner eukaryote : prokaryote measure.

**Hurdle C comparison** (`scripts/hurdle_mechanism.py` rerun; `mechanism.csv` now uses the new C). By read-based
domain (R67 table, new C; within-cruise presence AUC unless stated):

| Read domain | n | Granite | Granite + technical | Technical | Predicted C alone | Granite − C | Granite ahead by ≥ 0.05 | Abundance R² | r(depth) | r(C) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Eukaryotic | 2,756 | 0.635 | 0.630 | 0.534 | 0.614 | +0.017 | 13% | +0.253 | −0.094 | +0.363 |
| Mixed | 1,055 | 0.618 | 0.626 | 0.527 | 0.583 | +0.025 | 32% | +0.083 | +0.011 | +0.227 |
| Prokaryotic | 1,997 | 0.610 | 0.670 | 0.607 | 0.532 | +0.067 | 58% | −0.017 | +0.158 | −0.050 |

- The R67 conclusion stands. Eukaryotic-KO presence is mostly what predicted C gives (Granite ahead by ≥ 0.05 for 13%;
  7% with the old C). Prokaryotic-KO presence goes beyond it.

**Benchmark bundle** (`scripts/export_benchmark.py` rerun; `benchmark_export.zip` rezipped, 188.8 MB).
- Compared with the previous zip, decompressed contents: only `targets/C.csv.gz`, `targets_sample/C.csv.gz`,
  `granite/predictions_C.npz` and the C rows of `granite/scores.csv`, `partition.csv`, `cruise_decomposition.csv` and
  `c_per_cruise*.csv` differ.
- The standalone `benchmark_eval.py` run from the bundle (`python -I`) reproduces C exactly: 53.79% of non-technical,
  58.18% of reproducible.
- `~/multimodal_seusmbol/instructions_agent_biogoship.md` updated:
  - The C definition.
  - The Granite reference row for C: 53.79 (58.18) / 50.27 / 46.83 / 12.94.
  - Within-tile share for C: 6.4%.
  - A note listing the files to replace if an older copy is already on the other cluster.

**Legend draft** (R63, C figure): updated in place. The definition is now "eukaryote- to prokaryote-assigned KO reads",
and the clipping note now names time & place on I07N (−1.32), since Granite on I07N is no longer below −1.

### R68 addendum: paper figure legends checked against the new C and the read-based domain (2026-11-24)

The user's legend document (Figures 1–4, Supplementary S1) was checked. **Figures 1 and 3: no change.** They are built
on the all-KO CLR. Axis 1 still separates eukaryote- from prokaryote-dominated samples (within-cruise r −0.60 with the
new C).

**Figure 2** (C):
- In a, "eukaryote-specific to prokaryote-specific KO reads" becomes "eukaryote-assigned to prokaryote-assigned KO reads
  (each read takes the domain of the reference gene it maps to)".
- R² 0.54 / 53.8% unchanged.
- In b, the clipping note becomes "Values below −1 are drawn at −1 (time & place on I07N, −1.32)". Granite on I07N is now
  −0.57.

**Supplementary S1** (new C, Granite; per cruise from `predictions_C_granite_cv.npz`):
- 12 of 13 cruises within ±0.31 (was ±0.33). Largest |bias| apart from GOMECC4: I05-2023 0.303.
- The overall SD of the adjusted C is 0.91 (was 1.05).
- Cruises with observed SD ≥ 0.8: six (was seven; the legend said eight).
  - I05-2023 0.875 (r 0.475), I08S-2024 0.972 (0.917), A16N-2023 0.980 (0.712), A13.5-2024 1.002 (0.871), P18 1.156
    (0.816), HLY2202-2022 1.303 (0.424).
  - So r is 0.42–0.92, and 0.71–0.92 for four of the six.
  - Old values: seven cruises, r 0.45–0.91, with five at 0.75–0.91. The legend's "seven of those eight" did not match
    the old numbers either.
- Flat, warm cruises with low within-cruise r: I09N −0.08, I07N −0.05, GOMECC4 −0.27, I09N-2025 0.37, P02-2022 0.21. The
  list is unchanged.
- GOMECC4 bias −0.63 (was −0.67).
- HLY2202-2022 predicted SD 0.38 vs 1.30 observed (was 0.37 vs 1.35).

**Figure 4:** the module and category numbers are unchanged (all-KO CLR). The last sentence ("the categories at the
top of b consist largely of eukaryote-specific genes") was checked with the read-based domain (R67). Share of each
category's KOs with eukaryotic read share > 0.5:
- Cellular community – eukaryotes 95% (median read share 0.98).
- Signaling molecules and interaction 89% (0.98).
- Transcription 47%; transport and catabolism 52%; cell growth and death 45%.
- Folding, sorting and degradation 30%; signal transduction 34%; replication and repair 18%.
- At the bottom: cellular community – prokaryotes 0%, xenobiotics 4%, cell motility 43%.

The sentence is reworded to: "The two categories at the top of b consist mostly of eukaryotic genes (89–95% of their
KOs are mostly eukaryotic by reads), whose abundance follows the eukaryote : prokaryote balance; the next three are
about half eukaryotic."

## R69. Patchy KOs with one base model (Granite + depth + platform): no classes; temperature adds nothing (2026-11-24)

Replaces the satellite / technical / weak classes of R65 and the R66 mechanism figure (`hurdle_kos.png`,
`hurdle_mechanism.png` are superseded by `figures/output/hurdle_base.png`). Domain is the read-based eukaryotic share
of each KO's reads (R67): eukaryotic > 0.9 (2,756 KOs), prokaryotic ≤ 0.1 (1,997), mixed (1,055). C is the R68 read
ratio.

### Methods

- Same tiles, presence and folds as R65: 1,111 tiles, 13 cruises, 5,808 KOs present in 20–90% of tiles; presence =
  rarefied share ≥ 1/10^6 KO reads, tile present if the majority of its kept samples are. Each cruise held out; L2
  logistic, penalty (and, for Granite, 4 / 8 / 16 / 32 PCs) chosen by leave-cruise-out inside the training cruises.
  Score: within-cruise AUC (n-weighted mean over held-out cruises).
- Models per KO (`scripts/hurdle_base.py`; new fits run locally in 10 shards, `cluster/run_hurdle_base_local.sh`,
  ~50 min; tech / Granite / Granite + tech reused from R65 `hurdle_all.csv`, identical folds):
  - **tech**: platform group (NovaSeq 6000 group) + spline of log10 depth.
  - **temp_tech**: in-situ temperature spline (5 knots, cubic) + tech.
  - **c_tech**: Granite-predicted C + tech. Predicted C is nested: for held-out cruise c, a training tile in cruise t
    gets C from a ridge on 32 Granite PCs fitted without c and t; tiles of c get C from all cruises but c (ridge on
    all 1,215 benchmark tiles, penalty by leave-cruise-out). Held-out predicted C on the hurdle tiles: R² 0.539.
  - **base**: Granite PCs + tech (= R65 "both").
  - **base_temp**: Granite PCs + tech + temperature spline.
- Abundance where present (R65, unchanged): tiles where present (≥ 50 tiles, ≥ 4 cruises), technical variance removed
  in-sample, ridge on 32 Granite PCs, cruise held out; R² against the technically adjusted values.
- Depth and detection (panel d): within each cruise, tiles split into fifths by log10 depth; per KO, presence rate in
  each fifth divided by the KO's overall presence rate; median over KOs in each domain × gene-length tercile (terciles
  of mean reference gene length within each domain, R44 `ko_gc.csv`).

### 1. Temperature: no gain over the base model

Within-cruise AUC, base_temp − base: median +0.001 (10–90% −0.004 to +0.007); ≥ +0.02 for 2% of KOs, ≥ +0.05 for 0%.
Same in every domain (eukaryotic +0.002, mixed +0.001, prokaryotic +0.001). Base model = Granite + depth + platform.

### 2. Within-cruise presence AUC, median (10–90%), by read domain

| Domain | n | Depth + platform | Temp + tech | Pred. C + tech | Granite | Granite + tech (base) | Base + temp |
|---|---|---|---|---|---|---|---|
| all | 5,808 | 0.550 (0.47–0.66) | 0.620 (0.54–0.71) | 0.612 (0.53–0.70) | 0.626 (0.54–0.70) | 0.639 (0.56–0.73) | 0.641 (0.56–0.73) |
| eukaryotic | 2,756 | 0.534 (0.46–0.60) | 0.616 (0.56–0.67) | 0.617 (0.56–0.69) | 0.635 (0.57–0.70) | 0.630 (0.56–0.70) | 0.632 (0.57–0.70) |
| mixed | 1,055 | 0.527 (0.47–0.61) | 0.597 (0.51–0.68) | 0.594 (0.52–0.67) | 0.618 (0.54–0.70) | 0.626 (0.55–0.70) | 0.625 (0.55–0.71) |
| prokaryotic | 1,997 | 0.607 (0.51–0.72) | 0.650 (0.53–0.75) | 0.615 (0.51–0.73) | 0.610 (0.52–0.72) | 0.670 (0.56–0.77) | 0.672 (0.56–0.77) |

Pooled (not within-cruise) AUC, median: tech 0.512, temp + tech 0.694, C + tech 0.684, Granite 0.660, base 0.709,
base + temp 0.710 (eukaryotic base 0.733, mixed 0.685, prokaryotic 0.677).

Paired differences in within-cruise AUC (median, 10–90%; share of KOs ≥ +0.02 / ≥ +0.05 / ≤ −0.02):

| Comparison | all | eukaryotic | mixed | prokaryotic |
|---|---|---|---|---|
| base − tech | +0.081 (−0.002 to +0.169); 81 / 67 / 5% | +0.097 (+0.027 to +0.166); 91 / 80 / 2% | +0.088 (+0.008 to +0.190); 86 / 73 / 3% | +0.043 (−0.020 to +0.156); 65 / 45 / 10% |
| base − C + tech | +0.016 (−0.024 to +0.097); 45 / 22 / 12% | +0.010 (−0.023 to +0.046); 33 / 8 / 12% | +0.023 (−0.020 to +0.109); 52 / 27 / 10% | +0.030 (−0.027 to +0.152); 57 / 39 / 13% |
| base − temp + tech | +0.014 (−0.029 to +0.079); 44 / 21 / 17% | +0.012 (−0.029 to +0.062); 41 / 17 / 17% | +0.023 (−0.024 to +0.097); 53 / 30 / 13% | +0.011 (−0.034 to +0.092); 43 / 24 / 19% |
| C + tech − tech | +0.051 (−0.010 to +0.132); 66 / 50 / 5% | +0.085 (+0.027 to +0.144); 92 / 78 / 0% | +0.056 (−0.006 to +0.133); 72 / 54 / 2% | +0.002 (−0.025 to +0.051); 27 / 11 / 13% |
| base_temp − base | +0.001 (−0.004 to +0.007); 2 / 0 / 1% | +0.002; 1 / 0 / 0% | +0.001; 2 / 0 / 1% | +0.001; 4 / 0 / 1% |

Base − Granite alone (median): all +0.004, eukaryotic −0.004, mixed +0.002, prokaryotic +0.039. Median number of
Granite PCs chosen (base_temp): 16 overall and for eukaryotic, 8 for prokaryotic KOs. Median prevalence: eukaryotic
0.45, mixed 0.54, prokaryotic 0.56.

Spearman across KOs: base vs C + tech 0.637; base vs eukaryotic read share −0.205; tech vs eukaryotic read share
−0.338; abundance R² vs eukaryotic read share 0.504.

### 3. Abundance where present, by eukaryotic read share

| Eukaryotic read share | n | median R² | IQR |
|---|---|---|---|
| 0.0–0.1 | 1,997 | −0.017 | −0.047 to +0.044 |
| 0.1–0.2 | 124 | −0.009 | −0.039 to +0.057 |
| 0.2–0.3 | 99 | −0.004 | −0.043 to +0.077 |
| 0.3–0.4 | 74 | −0.006 | −0.030 to +0.046 |
| 0.4–0.5 | 63 | +0.007 | −0.034 to +0.098 |
| 0.5–0.6 | 70 | +0.086 | +0.021 to +0.179 |
| 0.6–0.7 | 117 | +0.132 | −0.002 to +0.222 |
| 0.7–0.8 | 173 | +0.119 | +0.006 to +0.276 |
| 0.8–0.9 | 335 | +0.192 | +0.067 to +0.341 |
| 0.9–1.0 | 2,756 | +0.253 | +0.081 to +0.377 |

By domain: prokaryotic median −0.017 (R² > 0.1 for 15% of KOs), mixed +0.083 (47%), eukaryotic +0.253 (73%).

### 4. Sequencing depth and detection, by gene length

Relative detection (presence rate in the depth fifth / the KO's mean rate; median over KOs), shallowest → deepest
fifth within cruise. Length tercile cuts (nt): prokaryotic 763 / 1,082; mixed 694 / 806; eukaryotic 690 / 784.

| Domain | Length | n | Q1 | Q2 | Q3 | Q4 | Q5 | absolute rate Q1 → Q5 |
|---|---|---|---|---|---|---|---|---|
| prokaryotic | short | 666 | 0.89 | 0.97 | 1.00 | 1.04 | 1.09 | 0.48 → 0.63 |
| prokaryotic | medium | 665 | 0.82 | 0.95 | 1.01 | 1.06 | 1.14 | 0.47 → 0.69 |
| prokaryotic | long | 666 | 0.81 | 0.94 | 1.00 | 1.07 | 1.18 | 0.44 → 0.68 |
| mixed | short | 352 | 1.02 | 1.03 | 0.99 | 0.97 | 1.00 | 0.57 → 0.56 |
| mixed | medium | 351 | 1.00 | 1.02 | 0.99 | 0.98 | 1.02 | 0.61 → 0.60 |
| mixed | long | 352 | 0.99 | 1.00 | 0.98 | 0.98 | 1.02 | 0.46 → 0.46 |
| eukaryotic | short | 919 | 1.18 | 1.08 | 0.98 | 0.90 | 0.86 | 0.53 → 0.38 |
| eukaryotic | medium | 918 | 1.14 | 1.07 | 0.98 | 0.92 | 0.89 | 0.57 → 0.44 |
| eukaryotic | long | 919 | 1.16 | 1.08 | 0.99 | 0.90 | 0.86 | 0.48 → 0.34 |

Within cruise, log10 depth and observed C correlate at r = −0.210 (eukaryote-rich tiles were sequenced less deeply),
so eukaryotic KOs are detected less often in the deeper tiles despite the rarefaction.

### 5. Focus genes (all prokaryotic by reads)

| Gene | KO | euk read share | tech | temp + tech | C + tech | base | base + temp | abundance R² |
|---|---|---|---|---|---|---|---|---|
| nifH | K02588 | 0.000 | 0.534 | 0.621 | 0.548 | 0.619 | 0.620 | −0.134 |
| nifD | K02586 | 0.000 | 0.588 | 0.709 | 0.595 | 0.694 | 0.692 | −0.103 |
| nifK | K02591 | 0.000 | 0.565 | 0.574 | 0.568 | 0.545 | 0.543 | −0.029 |
| amoB / pmoB | K10945 | 0.000 | 0.548 | 0.548 | 0.629 | 0.723 | 0.727 | −0.034 |
| nasA | K00372 | 0.010 | 0.576 | 0.710 | 0.560 | 0.632 | 0.696 | 0.014 |
| phnJ | K06163 | 0.000 | 0.759 | 0.797 | 0.792 | 0.785 | 0.785 | 0.024 |

(nasA is the one focus gene where temperature adds to the base model: +0.064.)

### Interpretation

- Granite + depth + platform is the base model: it beats depth + platform alone for 81% of patchy KOs by ≥ 0.02
  (median +0.081) and adding temperature changes nothing (+0.001). Temperature + technical alone gets most of the way
  (base ahead by a median 0.014).
- **Eukaryotic KOs**: presence follows the eukaryote : prokaryote balance. The predicted ratio + technical reaches
  0.617 against 0.630 for the base model; Granite is ahead by ≥ 0.05 for only 8%. Depth + platform alone barely
  predicts them (0.534). Where present, their abundance is predictable (median R² 0.253).
- **Prokaryotic KOs**: the predicted ratio adds nothing over depth + platform (+0.002), but Granite does (+0.043; ahead
  of C + tech by ≥ 0.05 for 39%). Depth matters: depth + platform alone reaches 0.607, and detection rises with depth
  within cruise, most for long genes (0.81 → 1.18 relative, vs 0.89 → 1.09 for short genes) — consistent with the
  coverage filter on mapped reads penalising long genes in shallow libraries. Where present, their abundance is not
  predictable (median R² −0.017).
- So the old "satellite vs technical" split was mostly this: eukaryotic patchy KOs track community composition, which
  Granite sees; prokaryotic patchy KOs are partly depth-limited detection, with a Granite-visible signal beyond
  composition for a subset (e.g. amoB 0.723 vs 0.629 with C).

### Figure legend draft (hurdle_base)

**Patchy genes: presence predicted from Granite + sequencing depth and platform.** 5,808 KOs present in 20–90% of 1,111
tiles (presence: ≥ 1 read per 10^6 KO reads after rarefaction). Each cruise held out in turn; AUC is computed within
each held-out cruise and averaged, weighted by tiles. KOs are grouped by the share of their reads assigned to
eukaryotic reference genes (eukaryotic > 0.9, prokaryotic ≤ 0.1, otherwise mixed). **a**, Presence AUC (median and
interquartile range) for four logistic models: sequencing depth + platform; in-situ temperature + depth + platform;
the Granite-predicted eukaryote : prokaryote read ratio + depth + platform; Granite PCs + depth + platform. Adding
temperature to the last model changes the median AUC by +0.001. **b**, Per-KO presence AUC, Granite + depth + platform
against predicted eukaryote : prokaryote ratio + depth + platform, coloured by the eukaryotic share of the KO's reads;
dashed line, equal AUC; circled, nitrogen-fixation (nifH, nifD, nifK), ammonia-oxidation (amoB), nitrate-assimilation
(nasA) and phosphonate (phnJ) genes. Eukaryotic KOs lie near the line (median gain 0.010); prokaryotic KOs lie above
it (median gain 0.030; ≥ 0.05 for 39%). **c**, Abundance where present: R² of a Granite ridge model (cruise held out,
depth and platform removed) by eukaryotic read share, in tenths; line, median; band, interquartile range; numbers,
KOs per bin. **d**, Detection rate relative to each KO's average across fifths of the tiles ranked by sequencing
depth within each cruise, for prokaryotic and eukaryotic KOs split into thirds by gene length (median over KOs).
Prokaryotic genes are detected more often in deeper libraries, long genes most (0.81 to 1.18); eukaryotic genes are
detected less often (1.16 to 0.86), because eukaryote-rich tiles were sequenced less deeply (within-cruise r = −0.21
between depth and the eukaryote : prokaryote ratio).

### R69 addendum: core KOs alongside patchy KOs; temperature dropped; figure redone as a–f (2026-11-24)

Decisions: temperature is not used (base_temp − base +0.001; its numbers stay above for the record); the figure no
longer shows temperature. Core and patchy KOs are shown side by side but scored separately: the patchy abundance R²
does not enter any main-model number.

Core KOs = the 3,933 KOs of the main (full-KO, CLR) model; per-KO held-out Granite R² from
`data/results/main_model/ko_ceiling.csv` (technical variance removed in-sample, cruise held out). 212 of them are also
in the patchy set (they pass both filters). Domain by read share as above. Depth sensitivity of core KOs = within-cruise
Pearson r between the CLR abundance (before technical adjustment, `targets/clr_all.csv.gz`) and log10 KO reads per
tile, over all 1,215 tiles; length terciles within domain. Outputs: `data/results/hurdle/core_domain_depth.csv`,
`depth_detection.csv`; numbers appended to `base_summary.txt`.

**Abundance R² by read domain (median, IQR, share of KOs with R² > 0.1)**

| | Prokaryotic | Mixed | Eukaryotic |
|---|---|---|---|
| Core KOs (main model) | n 3,179: +0.102 (+0.010 to +0.213), 50% | n 489: +0.132 (+0.035 to +0.223), 59% | n 265: +0.320 (+0.236 to +0.376), 92% |
| Patchy KOs, where present | n 1,997: −0.017 (−0.047 to +0.044), 15% | n 1,055: +0.083 (−0.015 to +0.235), 47% | n 2,756: +0.253 (+0.081 to +0.377), 73% |

Core R² by eukaryotic read share (tenths; n, median): 0.0–0.1 3,179 +0.102; 0.1–0.2 190 +0.104; 0.2–0.3 74 +0.133;
0.3–0.4 46 +0.066; 0.4–0.5 26 +0.058; 0.5–0.6 32 +0.138; 0.6–0.7 30 +0.212; 0.7–0.8 40 +0.228; 0.8–0.9 51 +0.259;
0.9–1.0 265 +0.320. Spearman (R², eukaryotic read share): core 0.116, patchy 0.504.

**Core KOs: within-cruise r (abundance, depth), median (IQR), by length tercile**

| Domain | n | Tercile cuts (nt) | Short | Medium | Long | Spearman (r, length) |
|---|---|---|---|---|---|---|
| Prokaryotic | 3,179 | 746 / 1,010 | −0.112 (−0.246 to −0.010) | +0.033 (−0.070 to +0.114) | +0.133 (+0.056 to +0.207) | 0.657 |
| Mixed | 489 | 665 / 765 | −0.126 (−0.218 to −0.043) | −0.048 (−0.127 to +0.014) | +0.010 (−0.079 to +0.077) | 0.489 |
| Eukaryotic | 265 | 687 / 754 | −0.228 (−0.250 to −0.201) | −0.221 (−0.247 to −0.188) | −0.202 (−0.248 to −0.174) | 0.162 |

Patchy detection in % vs the KO's average, depth fifths 1 → 5: prokaryotic short −11 −3 +0 +4 +9, medium −18 −5 +1 +6
+14, long −19 −6 +0 +7 +18; mixed short +2 +3 −1 −3 −0, medium +0 +2 −1 −2 +2, long −1 +0 −2 −2 +2; eukaryotic
short +18 +8 −2 −10 −14, medium +14 +7 −2 −8 −11, long +16 +8 −1 −10 −14.

**Interpretation**
- The domain gradient in abundance predictability holds for core KOs: eukaryotic 0.320, mixed 0.132, prokaryotic
  0.102. Unlike patchy prokaryotic KOs (−0.017), core prokaryotic KOs are predictable (half have R² > 0.1); the
  patchy result does not carry over to them.
- The depth × length effect is the same in core KOs: within a cruise, abundance of long prokaryotic genes rises with
  depth (+0.133) and of short ones falls (−0.112; CLR, so relative). This is the coverage filter on mapped reads; the
  main model removes depth per KO, so it does not enter the R². Eukaryotic KOs fall with depth at every length
  (≈ −0.2), i.e. composition: eukaryote-rich tiles were sequenced less deeply.
- Mixed KOs sit in between in every respect: presence AUC as eukaryotic (base 0.626), detection flat with depth
  (within ±3%), abundance R² between the two.

**Legend draft (replaces the R69 draft)**

**Patchy and core genes: what Granite predicts depends on whether the gene is carried by eukaryotes or prokaryotes.**
KOs are grouped by the share of their reads mapped to eukaryotic reference genes (eukaryotic > 0.9, prokaryotic ≤ 0.1,
otherwise mixed). Patchy KOs: 5,808 KOs present in 20–90% of 1,111 tiles (≥ 1 read per 10^6 KO reads after
rarefaction). Core KOs: the 3,933 KOs of the main model. Every model holds out each cruise in turn. **a**, Presence of
patchy KOs: AUC within the held-out cruise (median, interquartile range) for logistic models on sequencing depth +
platform, the Granite-predicted eukaryote : prokaryote read ratio + depth + platform, and Granite PCs + depth +
platform. **b**, Abundance R² (Granite ridge, depth and platform removed) for core KOs (main model) and for patchy KOs
where present, scored separately; numbers, KOs. **c**, Per-KO presence AUC of patchy KOs, Granite + depth + platform
against predicted ratio + depth + platform, coloured by eukaryotic read share; dashed line, equal AUC; labelled,
nitrogen-fixation (nifH, nifD, nifK), ammonia-oxidation (amoB), nitrate-assimilation (nasA) and phosphonate (phnJ)
genes. **d**, Abundance R² against eukaryotic read share (tenths; line, median; band, interquartile range). **e**,
Patchy KOs: how often a KO is detected in each fifth of the tiles ranked by sequencing depth within their cruise,
as % above or below its average over all tiles, for short, medium and long genes (thirds by mean reference gene
length within each domain; median over KOs). **f**, Core KOs: within-cruise correlation between abundance (before
depth and platform are removed) and sequencing depth, by gene length (median, interquartile range). Longer
prokaryotic genes are detected more often and counted more in deeper libraries; eukaryotic genes go the other way
because eukaryote-rich tiles were sequenced less deeply (within-cruise r = −0.21 between depth and the eukaryote :
prokaryote ratio).
