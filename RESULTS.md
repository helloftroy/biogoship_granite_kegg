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

    OMP_NUM_THREADS=4 ~/miniforge3/bin/python scripts/residual_modules.py

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
