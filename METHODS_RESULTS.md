# Methods and results log

Working notes for the Bio-GO-SHIP x Sentinel-3 x Granite project: what was done,
why, and what came out. Running instructions are in [README.md](README.md).
Newest results are added at the end of each section, with the date.

## 1. Aim

Predict surface-ocean microbial community structure (NMF communities from eDNA)
and metabolic potential (KEGG profiles from metagenomes) from satellite data:
Sentinel-3 tiles -> Granite Geospatial Ocean embeddings -> CatBoost -> NMF / KEGG.

New relative to the Florida project (`multimodal_eDNA`): samples are global
transects rather than one region, and each sample gets two tiles per overpass:
one at the sampling position and one at the position its water occupied when the
satellite passed (section 5).

## 2. Samples

**Source.** Bio-GO-SHIP (Larkin et al. 2021), BioProject PRJNA656268.

- Sample list, positions and times: ENA run report for PRJNA656268, fetched live
  by `scripts/build_sample_manifest.py`. 2,407 runs: 971 metagenomes (WGS) and
  1,436 amplicon runs.
- Extras for the metagenomes (section, depth, sampling method, temperature,
  phosphate, nitrate+nitrite): `BGS_patch2.csv` from
  merenlab/public-marine-omics-metadata, joined on BioSample. That table only
  holds metagenomes from 100 m or shallower with a temperature value, so 27
  samples have no depth or method.

**Units.**

- *Sample* = one BioSample (one filter). 996 in total: 950 have both a
  metagenome and amplicons, 21 a metagenome only, 25 amplicons only.
- *Site* = a unique position and collection time. 954 in total. Replicates share
  a site, and the site is what gets matched to satellite data.

**Result (2026-10-12): samples by section.**

| Section | Ocean | Dates | Samples | Metagenomes | Sites |
|---|---|---|---:|---:|---:|
| C13.5/A13.5 | Atlantic | 2020-03-21 to 2020-04-17 | 263 | 251 | 234 |
| I07N | Indian | 2018-04-25 to 2018-06-04 | 261 | 248 | 259 |
| I09N | Indian | 2016-03-21 to 2016-04-24 | 242 | 242 | 242 |
| P18 | Pacific | 2016-11-20 to 2017-01-28 | 104 | 104 | 104 |
| AMT28 | Atlantic | 2018-09-25 to 2018-10-27 | 63 | 63 | 63 |
| NH1418 | Pacific | 2014-09-20 to 2014-10-06 | 29 | 29 | 27 |
| AE1319 | Atlantic | 2013-08-26 to 2013-09-09 | 18 | 18 | 13 |
| BVAL46 | Atlantic | 2011-09-29 to 2011-10-13 | 16 | 16 | 12 |
| **Total** | | | **996** | **971** | **954** |

Most samples are near-surface: median depth 7 m (underway intake); the deepest
is 80 m. Samples below the surface layer are kept in the manifest with their
depth so they can be filtered later.

## 3. Matching sites to Sentinel-3

**Products.**

- OLCI level-2 water full resolution (`OL_2_WFR`, 300 m): water-leaving
  reflectances and quality flags.
- SLSTR level-2 sea surface temperature (`SL_2_WST`): one GHRSST L2P file per orbit.
- Non-time-critical (NT) processing only; where a scene was reprocessed, the
  newest version is used.

**Search.** `scripts/search_products.py` queries the Copernicus catalogue once
per site for OLCI scenes whose footprint contains the site and whose sensing
start is within the time window of sampling. All such scenes are kept as
candidates, not just the closest, because any one overpass may be cloudy.

**Time window.** +/-3 days, the same as the Florida project. For I09N only, +/-10
days (see below).

**SST pairing.** Each OLCI scene is paired with the SST orbit file whose time
span contains the OLCI sensing start, from the same satellite. If the same
satellite has none within 180 minutes, the other satellite's closest orbit file
is used and flagged `cross_satellite`.

**Result (2026-10-12): candidate overpasses.**

| Section | Sites | Sites with a candidate scene | Candidate scenes per site (median, range) | Gap to closest scene, hours (median, max) | Window |
|---|---:|---:|---|---|---|
| C13.5/A13.5 | 234 | 234 | 6 (5-8) | 4.3, 24.0 | +/-3 d |
| I07N | 259 | 259 | 4 (2-8) | 10.2, 35.2 | +/-3 d |
| P18 | 104 | 104 | 4 (2-12) | 8.9, 35.8 | +/-3 d |
| AMT28 | 63 | 63 | 6 (4-10) | 5.2, 29.6 | +/-3 d |
| I09N | 242 | 61 | 3 (1-4) | 140.4, 237.4 | +/-10 d |
| NH1418 | 27 | 0 | | | before Sentinel-3 |
| AE1319 | 13 | 0 | | | before Sentinel-3 |
| BVAL46 | 12 | 0 | | | before Sentinel-3 |
| **Total** | **954** | **721** | | | |

- 721 sites (752 samples, 727 of the 971 metagenomes) have at least one candidate scene.
- 3,727 site x scene pairs; 414 unique OLCI scenes (160 GB) and 218 SST orbit
  files (141 GB).
- SST pairing: 333 scenes same-satellite, 81 cross-satellite. The cross-satellite
  ones are Sentinel-3B scenes from spring 2018 (I07N), before Sentinel-3B SST
  products were available.
- "Candidate" means the site is inside the scene footprint. Whether the site is
  cloud-free is only known once the scene is downloaded.

**Sections that cannot be matched, and why.**

| Section | Samples | Reason |
|---|---:|---|
| NH1418 (2014), AE1319 (2013), BVAL46 (2011) | 63 | Sampled before Sentinel-3A launched (February 2016). Dropped. |
| I09N (2016) | 181 of 242 | The cruise ended 2016-04-24. The earliest OLCI level-2 scene over the transect is 2016-04-27, so only the last part of the cruise is within reach. |

I09N by window (all matches are scenes taken *after* sampling):

| Window | I09N sites with a candidate scene (of 242) | Scenes |
|---|---:|---:|
| +/-3 days | 7 | 1 |
| +/-5 days | 24 | 2 |
| +/-7 days | 40 | 5 |
| +/-10 days | 61 | 8 |

Reaching all 242 would need about +/-37 days, which is not a meaningful match.
The 61 sites at +/-10 days are included in the download (8 scenes, about 3 GB)
with their time gap recorded (57 to 237 hours), so they can be kept or dropped
at modelling time. Note for back-trajectories: these scenes post-date sampling,
so they show where the water was after it was sampled, not before.

## 4. Download

`scripts/download_products.py` downloads every candidate OLCI scene and its
paired SST orbit file as complete product archives (about 300 GB) to
`data/products/` on cluster scratch. Whole scenes are kept, rather than tiles
around the sites, so that later steps can read any location in a scene; the
back-trajectory analysis needs positions upstream of the sampling site.

**Result (2026-10-13):** download completed on the cluster (reported by the
project owner; archive counts and total size still to be recorded here from
`python scripts/download_products.py --dry-run`).

## 5. Back-tracking the sampled water

**Why.** The satellite sees a site hours to days before or after it was sampled
(median gap 4 to 10 hours for the closest scene, up to 3 days across all
candidates, up to 10 for I09N). Surface water moves in that time. Two versions
of each tile are therefore made, as a test of whether accounting for the
movement helps:

- **fixed**: the tile at the sampling position, exactly as in the Florida project.
- **advected**: the tile at the position the sampled water occupied at the
  scene's sensing time.

Both go through Granite later; the comparison is between embeddings from the
two positions.

**Currents.** Copernicus Marine MULTIOBS_GLO_PHY_MYNRT_015_003 (GlobCurrent):
total velocity `uo`, `vo` combining altimeter-derived geostrophic currents and
modelled wind-driven Ekman currents, 0.25 degree, global. Dataset
`cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m` (multi-year reprocessing, daily
means). Levels: 0 m and 15 m; **15 m is used by default** (see choices below).
One subset per section: site extent plus 5 degrees, from 8 days before the
first to 8 days after the last sampling or overpass time (to cover 7-day tracks).

**Tracking** (`biogoship/advection.py`, `scripts/backtrack_sites.py`).

- OceanParcels 3.x, one particle per site, released at the sampling position
  and time.
- Fourth-order Runge-Kutta, 20-minute time step, positions recorded hourly.
- Run once backward and once forward, each for 7 days (or the largest
  sampling-to-overpass gap in the section if longer: 237 h for I09N).
- Deterministic: no diffusion, no ensemble, horizontal motion only.
- For each (site, scene) pair the track is linearly interpolated to the scene's
  sensing start. A scene earlier than sampling uses the backward track ("where
  the water came from"); a later scene uses the forward track ("where it went").
- Missing currents (land, coast, ice) are treated as zero velocity. A particle
  that never moves is flagged `stalled`; one that leaves the currents subset
  before the overpass time is flagged `left_domain`. Neither gets an advected tile.

**Choices to revisit.**

- *15 m vs 0 m.* The 0 m level includes the strongest direct wind drift; 15 m is
  the level of drogued drifters and is closer to the bulk flow of the mixed layer
  that an intake at ~7 m samples. 15 m is the default; `--depth-m 0` gives the other.
- *Daily means vs hourly.* The hourly dataset
  (`cmems_obs-mob_glo_phy-cur_my_0.25deg_PT1H-i`) resolves inertial and tidal
  motion; daily means were chosen for a first, smoother estimate. Set with
  `--dataset-id` in `download_currents.py`.
- *Resolution.* 0.25 degree currents do not resolve submesoscale motion, and one
  tile is 12.6 km across. Small displacements are within the uncertainty of the
  current field; the advected tile is an estimate, not a measurement.
- *Coverage.* Scenes were searched by whether they contain the sampling
  position. If the advected position falls outside a scene, that pair is logged
  `out_of_swath` and has no advected tile.

**Verification.** In a uniform 0.5 m/s eastward test current, 24 h of tracking
gives 43.2 km in both directions (expected 43.2 km), and a full synthetic run
over all 3,727 site x scene pairs matched the expected drift to within 0.2%.
This checks the tracking code, not the realism of the currents.

**Result (2026-10-13): first run on real currents** (15 m level, daily means).
All 3,727 site x scene pairs resolved (`ok`).

Straight-line distance between the sampling position and the advected position
at overpass time:

| Section | Pairs | Median | Maximum |
|---|---:|---:|---:|
| AMT28 | 431 | 9.0 km | 79.7 km |
| P18 | 472 | 9.5 km | 273.5 km |
| C13.5/A13.5 | 1,495 | 14.0 km | 308.0 km |
| I07N | 1,159 | 23.3 km | 153.7 km |
| I09N (overpasses 2 to 10 days after sampling) | 170 | 121.6 km | 223.3 km |

For comparison, a tile is 12.6 km across: the typical displacement within +/-3
days is about one tile width, and for a sizeable share of pairs it is several.
For I09N the advected position is, in the median, about ten tile widths from the
sampling position.

The first run reused currents files downloaded with only 2 days of padding, so
some 7-day tracks stopped early. `download_currents.py` now checks coverage and
re-downloads; on the re-run (same day) every track in AMT28, P18, I07N and
C13.5/A13.5 ran the full 7 days both ways, and the displacements above were
unchanged. For I09N, 20 of 61 backward tracks stop between 7.3 and 9.9 days;
this does not matter, since all I09N overpasses come after sampling.

## 6. Tiles and cloud screening (written, not run)

Carried over unchanged from the Florida project: band selection, quality flags,
normalisation, SST resampling. Changed: the tile is centred on the site instead
of being cut from a fixed grid over a region.

For each site and candidate scene (`scripts/download_tiles.py`):

1. Find the OLCI pixel nearest the site. If it is more than 1 km away, the site
   is outside the swath (`out_of_swath`).
2. Take the 42x42 pixel window (12.6 km) centred on that pixel. If it runs off
   the scene, skip (`edge`). (With `--max-offset-km`, a window that fails here or
   at step 3 is replaced by the nearest clear window within that distance.)
3. Compute the fraction of clear open-ocean pixels: WATER set, and none of
   INVALID, LAND, CLOUD, CLOUD_AMBIGUOUS, CLOUD_MARGIN, SNOW_ICE, INLAND_WATER,
   COASTLINE, SUSPECT, SATURATED, AC_FAIL. Below 0.8, skip (`cloudy`).
4. Read 16 reflectance bands (Oa01-Oa12, Oa16, Oa17, Oa18, Oa21), convert to
   physical units and z-score with the Granite pretraining mean and standard
   deviation per band. Fill pixels are set to 0 (the band mean).
5. Resample SST to the tile grid by nearest neighbour, using only pixels with
   GHRSST quality level 4 or 5, and z-score (mean 293.91 K, s.d. 2.68 K). Where
   there is no acceptable SST pixel the channel is NaN.
6. Save `(17, 42, 42)` float32 plus the valid-pixel mask and metadata.

Steps 1 to 6 are applied twice per site and scene: once at the sampling position
(`fixed`) and once at the back-tracked position (`advected`). The two tiles can
differ in cloud cover, so a site may have one without the other.

**Cloud handling.** Cloud cover is not known from the catalogue; it is measured
here, per tile, from the pixel flags. All candidate scenes within the time window
are tested, so if the closest overpass is cloudy at a position, the next closest
clear one is used instead. Nothing cloudy or out of scene reaches the embedding
step: such cases never produce a tile. They are recorded in `data/tile_log.csv`
(one line per site, scene and position, with the clear fraction and the dominant
rejection flags) and summarised per site and position in
`data/manifests/site_tile_status.csv` with a `flag`:

| Flag | Meaning | Possible response |
|---|---|---|
| `ok` | at least one clear tile | none |
| `all_cloudy` | every scene covers the position; none is 80% clear | widen the time window, or lower the threshold (`best_cloudy_valid_fraction` shows how close it came) |
| `outside_all_scenes` | the position is off the edge of, or outside, every scene | search for scenes containing the advected position |
| `cloudy_or_outside` | a mix of the two | either of the above |
| `no_position` | back-tracking gave no advected position | check `advected_positions.csv` status |

**Result (2026-10-13): advected case, scenes within the original window** (80%
clear, nearest clear tile within 20 km). 1,635 advected tiles. Over both
positions: 3,327 clear, 3,962 cloudy, 77 at a scene edge, and 88 where the
advected position fell outside the scene (searched by the sampling position).

| Section | Samples | Fixed tile | Advected tile | Both |
|---|---:|---:|---:|---:|
| C13.5/A13.5 | 263 | 251 | 246 | 243 |
| I07N | 261 | 194 | 195 | 185 |
| AMT28 | 63 | 46 | 46 | 45 |
| P18 | 104 | 67 | 65 | 62 |
| I09N | 242 | 56 | 59 | 54 |
| **Total** | **996** | **614** | **611** | **589** |

By site, 586 of 721 have an advected tile (the same count as fixed, by
coincidence; the sets differ). 135 sites lack one: 112 cloudy in every scene,
23 cloudy or outside. These are the sites for the wider search along the track.
636 samples have at least one of the two.

Every clear overpass yields a tile, so a site can have several.
`scripts/build_tile_manifest.py` joins samples to tiles. For each sample and
position it marks the tile closest in time as `is_best`, and it marks the advected
tile from the same overpass as the best fixed tile (`same_scene_as_best_fixed`),
which is the like-for-like pair for comparing the two approaches.

**Result (2026-10-13): fixed position, centred window, 80% clear, +/-3 days (I09N +/-10).**

All 3,727 site x scene combinations were tested: 674 clear (18%), 2,992 cloudy,
61 at a scene edge.

| Section | Samples | Samples with a clear tile | Share |
|---|---:|---:|---:|
| C13.5/A13.5 | 263 | 189 | 72% |
| I07N | 261 | 113 | 43% |
| I09N | 242 (61 matchable sites) | 47 | 19% of all samples |
| AMT28 | 63 | 18 | 29% |
| P18 | 104 | 29 | 28% |
| Pre-Sentinel-3 cruises | 63 | 0 | |
| **Total** | **996** | **396** | **40%** |

By site: 376 of the 721 matchable sites have a clear tile (52%); 314 are cloudy in
every scene; 31 are cloudy or off the edge in every scene. Advected position not
run yet.

So about one overpass in five is clear at a given site, and having four to six
overpasses per site lifts coverage to about half. Every candidate scene was
already tried; the misses are sites cloudy in all of them.

**Result (2026-10-13): fixed case, final rule** (80% clear, nearest clear tile
within 20 km, +/-3 days; I09N +/-10). The 3,053 cloudy or edge combinations were
re-tested with the 20 km allowance and 1,018 of them gave a tile, so 1,692 of
3,727 site x scene combinations are now clear (45%, up from 18%).

| Section | Samples | With a tile: centred only | With a tile: within 20 km | Share |
|---|---:|---:|---:|---:|
| C13.5/A13.5 | 263 | 189 | 251 | 95% |
| I07N | 261 | 113 | 194 | 74% |
| AMT28 | 63 | 18 | 46 | 73% |
| P18 | 104 | 29 | 67 | 64% |
| I09N | 242 | 47 | 56 | 23% (92% of the 61 reachable sites) |
| Pre-Sentinel-3 cruises | 63 | 0 | 0 | |
| **Total** | **996** | **396** | **614** | **62%** |

By site: 586 of the 721 reachable sites have a tile (81%, up from 52%); 123 are
cloudy in every scene even within 20 km; 12 are cloudy or off the edge in every
scene. 1,775 sample x tile rows in the manifest.

With the 20 km rule, sites with a tile by largest time gap allowed (of 660
outside I09N): 150 within 6 h, 236 within 12 h, 335 within 24 h, 477 within 48 h,
530 within 72 h. The remaining rejections are almost all true cloud (CLOUD on 72%
of the pixels of a rejected window).

Tile offsets: the 674 tiles from the first run are centred on the sampling
position (offset under one pixel). Of the 1,018 tiles added by the 20 km rule,
994 are shifted by more than 0.5 km: median 8.3 km, maximum 19.9 km. So a third
of fixed tiles are not centred on the sample, typically by less than a tile
width. That is comparable to the median water displacement at overpass time
(9 to 23 km), which matters when comparing the fixed and advected cases.

**Yield report (2026-10-13, `scripts/cloud_report.py`, fixed position, before the 20 km rule).**

Sites with a tile, of 721, as the required clear fraction is lowered:

| Clear fraction | 0.9 | 0.8 | 0.7 | 0.6 | 0.5 | 0.3 |
|---|---:|---:|---:|---:|---:|---:|
| Sites | 319 | 376 | 416 | 442 | 476 | 539 |

Sites with a tile (80% clear), of the 660 outside I09N, by largest time gap allowed:

| Gap | 6 h | 12 h | 24 h | 48 h | 72 h |
|---|---:|---:|---:|---:|---:|
| Sites | 74 | 117 | 164 | 276 | 329 |

- Only 74 sites have a clear tile within 6 hours of sampling and 164 within a
  day; for half of the sites with a tile, the best one is more than a day from
  sampling. This is the case for the advected position.
- The gain per extra day is shrinking (+112 sites from day 1 to 2, +53 from day
  2 to 3) but has not levelled off, so a wider window would add more.
- Rejections are real cloud: CLOUD covers 55% of the pixels of a rejected
  window on average, CLOUD_MARGIN 17%, CLOUD_AMBIGUOUS 7%. SUSPECT and AC_FAIL
  are about 1% each, so loosening those would gain almost nothing. (HIGHGLINT,
  15%, and OCNN_FAIL, 9%, appear in the report but are not rejection criteria,
  as in the Florida project.)
- 674 clear tiles over 376 sites is 1.8 per site with a tile, and 0.9 per
  matchable site. Averaging 3 to 5 clear views per sample would therefore need
  roughly four times as many overpasses per site, i.e. a window of about +/-10
  to 15 days along the trajectory.

**Decision (2026-10-13).** The two cases to compare are now defined as:

- **fixed**: sampling position, overpasses within +/-3 days of sampling (I09N
  +/-10), tile 80% clear, nearest clear tile within 20 km accepted. This is the
  Florida project's rule.
- **advected**: position of the sampled water at overpass time, same 80% and
  20 km rule, and the time window widened along the track for sites that have no
  clear tile within +/-3 days (out to 6 days to start with, in stages).

The aim for now is one good tile per sample in each case; averaging several
views along the track is set aside. Both cases then go through Granite.

**Widening along the track** (`scripts/search_trajectory_products.py`). Each site
is tracked 7 days each way. For every site still without a clear advected tile,
and for every day between the window already searched and the new limit, the
catalogue is queried for OLCI scenes covering the water's position at the middle
of that day. New (site, scene) pairs are marked `source=trajectory` in
`site_products.csv` and are tiled at the advected position only. The exact
advected position at each scene's sensing time comes from the hourly track. In a
test on AMT28 (all 63 sites, days 3 to 5), 53 sites gained candidate scenes: 282
new pairs from 64 new scenes.

**Result (2026-10-19): search out to +/-6 days** for the 135 sites without an
advected tile (798 catalogue queries). 133 of the 135 gained candidate scenes:
750 new site x scene pairs from 113 new OLCI scenes (32 GB) and up to 79 SST
files (52 GB). With what was already on disk, 146 new archives (about 54 GB), for
a total of 527 OLCI scenes and 251 SST files (355 GB).

**Result (2026-10-20): advected tiles after the +/-6 day search.** 176 more
advected tiles. Over both positions: 3,503 clear, 4,527 cloudy, 85 edge, 89
outside the scene.

| Section | Samples | Fixed tile | Advected tile | Both |
|---|---:|---:|---:|---:|
| C13.5/A13.5 | 263 | 251 | 261 | 251 |
| I07N | 261 | 194 | 224 | 193 |
| P18 | 104 | 67 | 91 | 66 |
| AMT28 | 63 | 46 | 59 | 46 |
| I09N | 242 | 56 | 59 | 54 |
| **Total** | **996** | **614** | **694** | **610** |

By site, 663 of 721 now have an advected tile (92%, up from 81%); 58 do not
(47 cloudy in every scene, 11 cloudy or outside). 698 samples have at least one
of the two cases. Note that the advected tiles found by the wider search are 3
to 6 days from sampling, while fixed tiles are within 3 days, so the two cases
now also differ in time gap; for a like-for-like comparison use the samples
with both, and `same_scene_as_best_fixed`.

**Ways to raise the yield** (`scripts/cloud_report.py` quantifies the second and
third from the existing log):

1. *Nearest clear window* (`--max-offset-km`): if the window centred on the
   position is not clear, take the closest clear window within a set distance.
   The Florida project accepted tiles up to 20 km from the sample, so 20 km is
   consistent with it. The offset is recorded per tile (`offset_km`).
2. *Lower the clear fraction* below 0.8. Masked pixels then make up more of the
   tile; how Granite handles them is untested.
3. *Widen the time window* beyond +/-3 days: more scenes to download, and a
   looser match in time (the advected position matters more as the gap grows).

**Possible extension.** Tiles at several points along each trajectory (for
example every 24 h back), giving a short history of the water mass. The hourly
tracks are already saved in `data/advection/trajectories.csv`; this would need a
further scene search along the tracks.

## 7. Granite embeddings (written, not run)

`scripts/embed_tiles.py`, model loading carried over from the Florida project
(`biogoship/granite.py`).

- Model: `ibm-granite/granite-geospatial-ocean`, a masked autoencoder; only the
  encoder is used. It is built as a TerraTorch `prithvi_eo_tiny` backbone
  (17 channels, 42x42 input, 2x2 patches, width 512, 12 layers, 8 heads) and the
  encoder weights are loaded from `checkpoint.pt`; decoder weights are dropped.
- Input: the `(17, 42, 42)` tile, with NaN (missing SST) set to 0, the channel mean.
- Embedding: mean over the encoder's output tokens, 512 numbers per tile.
- Every clear tile is embedded, for both positions and all clear overpasses, not
  only the best one. The per-sample table carries `position`, `is_best` and
  `same_scene_as_best_fixed` to select from them.

**Model check (2026-10-13, `cluster/cache_granite.sbatch`):** the checkpoint
loaded into the TerraTorch encoder with 150 tensors matched, 56 decoder tensors
skipped, no unmatched checkpoint tensors and no encoder tensors left without
weights; a test batch gives embeddings of shape (2, 512). torch 2.14.1.

**Result (2026-10-20):** 3,503 tiles embedded (1,692 fixed, 1,811 advected) in
101 s on one A100. Checks on the returned bundle:

- All values finite; no constant dimensions; one row per tile, in order.
- 227 tiles have an identical embedding to another tile because the tiles are
  identical: mostly the fixed and advected positions of a site snapping to the
  same nearest clear window in the same scene (95 cases), or two neighbouring
  sites doing so (70 fixed, 43 advected). Expected under the 20 km rule.
- Embeddings are strongly aligned: cosine similarity between random tiles has a
  median of 0.78, and the first three principal components hold 84% of the
  variance. Fixed and advected tiles from the same site and scene have a median
  cosine similarity of 0.999.

### Input audit against the model (2026-10-20)

Checked against the model's own `config.yaml` and checkpoint (Hugging Face
snapshot `319ac64`), IBM's backbone code (`custom_modules/prithvi_vit_S3.py` and
`configs/config-fine-tuning.yaml` in github.com/ibm-granite/geospatial) and the
paper (Dawson et al. 2025, arXiv:2509.21273). Items marked "verified" were
checked by `scripts/check_granite_inputs.py` against the model files.

| Item | Granite pretraining | Our tiles | Status |
|---|---|---|---|
| Tile size | 42 x 42 pixels (`img_size [1, 42, 42]`); the 224 x 224 of Prithvi-EO was reduced for Sentinel-3's 300 m pixels | 42 x 42 | match, verified |
| Patches | 2 x 2 pixels, 21 x 21 = 441 tokens, attention within the tile | same | match, verified |
| Pixel grid | native OLCI full-resolution (300 m) swath grid, cut into non-overlapping tiles; SST resampled onto the OLCI grid | native OLCI swath grid, no reprojection; SST nearest-neighbour onto the OLCI grid | match (their SST resampling method is not stated) |
| Padding | none needed at 42 x 42 | none; tiles that would run off the swath are skipped | match |
| Sensor and product | Sentinel-3A/B OLCI L2 WFR water reflectances + SLSTR L2 SST, 2017-2021 | the same products, 2016-2020 | match |
| Bands and order | Oa01-Oa12, Oa16, Oa17, Oa18, Oa21, SST | the same | match, verified |
| Units | physical reflectance (band means ~0.02) and SST in kelvin (mean 293.9) | scale factor and offset applied; kelvin | match |
| Normalisation | per-band mean and s.d. in config.yaml | identical (largest relative difference 4e-8, float rounding) | match, verified |
| Other inputs | IBM's backbone lists coordinate encodings (day length, distance to coast, region, biome), but TerraTorch only acts on "time" and "location", so they are inert; the paper describes only 2D sin/cos positional embeddings | none | match |
| Positional embedding | fixed 2D sin/cos | checkpoint values equal the sin/cos formula to 8e-7 | match, verified |
| Weights | encoder of the MAE checkpoint | 150 encoder tensors loaded, none missing; 56 decoder tensors dropped | match |
| Quality screening | product-notice recommended flags masked "with the exception of the cloud flag"; tiles kept if at least 80% cloud-free | tiles kept if at least 80% of pixels are clear of cloud, cloud margin, ambiguous cloud and the other quality flags | similar; ours slightly stricter |
| Masked and missing pixels | not documented for pretraining. IBM's fine-tuning config sets missing values to 0 in physical units before normalisation (TerraTorch `no_data_replace: 0`) | flagged pixels keep their measured values; missing SST and OLCI fill values set to the channel mean (0 after normalisation) | **open**: tested below |

The one convention that cannot be read from the documentation is how pixels
that are flagged or missing were filled during pretraining. Cloud pixels were
apparently kept with their measured values, as in ours. `check_granite_inputs.py`
therefore runs real tiles through Granite's own masked-autoencoder decoder under
three conventions (ours; flagged and missing pixels set to the channel mean;
set to physical zero as in IBM's fine-tuning config) and measures reconstruction
error on hidden clear pixels: the convention closest to pretraining should give
the lowest error, and all should beat a constant guess by a wide margin if the
inputs are in the form the model expects. It also reports how much the
embeddings change between conventions.

**Result (2026-10-20, 300 real tiles, cluster job 80958).**

- Configuration, band order, normalisation and positional embedding: all match
  the model's files.
- Our tiles are bluer and warmer than Granite's average. On clear pixels, after
  normalisation, Oa01-Oa03 average +1.4 to +1.8 s.d. and SST +1.6 s.d.; the red
  and near-infrared bands sit slightly below average with little spread. That is
  what clear subtropical open ocean looks like (water-leaving reflectance near
  0.04 at 400 nm, SST near 298 K), against a pretraining set that over-sampled
  coastal water. The size of the offsets rules out a unit error: a wrong scale
  (for example reflectance divided by pi, or SST in Celsius) would put them tens
  of s.d. away.
- SST is missing on 24% of pixels (71 of 300 tiles have some gaps); OLCI fill
  values on 1.4%.
- Reconstruction of hidden patches (clear pixels; normalised units): our
  convention 1.204, flagged pixels to the channel mean 1.207, IBM's physical-zero
  fill 1.217, against 1.330 for a constant guess. Our convention is the best of
  the three, by a small margin.
- Embedding sensitivity: setting flagged pixels to the channel mean barely moves
  the embeddings (median cosine similarity 0.995). The physical-zero convention
  moves some of them drastically (minimum 0.12), because a missing SST pixel then
  becomes -110 s.d.; it should not be used with these tiles.

Conclusion on gap filling: keep the current convention.

**How "bluer and warmer" is measured.** Every tile is stored normalised with
Granite's own statistics: value = (physical value - Granite mean) / Granite s.d.
So a normalised value of 0 is Granite's pretraining average and +1 is one
pretraining s.d. above it. Averaging the normalised clear pixels of our tiles
and converting back to physical units (`tile_raw_features.csv`, all 3,503 tiles):

| Quantity | Granite pretraining average | Our tiles, median |
|---|---:|---:|
| Water-leaving reflectance at 400 nm (Oa01) | 0.0235 | 0.0365 |
| Water-leaving reflectance at 560 nm (Oa06) | 0.0077 | 0.0041 |
| Blue-to-green ratio (Oa01 / Oa06) | 3.1 | 9.0 |
| SST | 20.8 C | 25.6 C |

83% of our tiles are bluer than Granite's average at 400 nm and 67% are warmer.
The blue-to-green ratio, the basis of standard chlorophyll algorithms, is three
times Granite's average: our water is much clearer, i.e. lower in chlorophyll.

| Section | Tiles | Median latitude (absolute) | Oa01, normalised | SST, normalised | SST |
|---|---:|---:|---:|---:|---:|
| C13.5/A13.5 | 1,549 | 21 | +2.85 | +1.75 | 25.5 C |
| P18 | 280 | 21 | +2.03 | +1.13 | 23.8 C |
| AMT28 | 372 | 26 | +1.48 | +0.94 | 23.3 C |
| I07N | 1,086 | 9 | +1.00 | +2.38 | 27.1 C |
| I09N | 216 | 12 | +0.88 | +3.71 | 30.7 C |

The Atlantic and Pacific transects cross the subtropical gyres (clearest water);
the Indian Ocean transects are the warmest (I09N: Bay of Bengal and eastern
Indian Ocean in April-May, near 31 C).

Our tiles also vary more than Granite's training data: between tiles, the 400 nm
band spans 1.8 Granite s.d., while within a tile it varies by only 0.4.

**Why Granite's average is not open-ocean water.** Most of the ocean by area is
indeed clear blue open water, but Granite's pretraining set was not sampled by
area (paper, section 2.1):

- An equal number of tiles (6,400) was taken from every Longhurst
  biogeographic province. Many provinces are coastal or shelf provinces, each
  small in area, so equal sampling per province gives coastal water far more
  weight than its share of the ocean; the paper notes that "by design the
  sampling strategy results in more data being used near the coast". Coastal
  water is greener and browner (chlorophyll, sediment, dissolved organic matter
  absorb blue light), which lowers the average 400 nm reflectance and raises
  green and red.
- Every month of 2017-2021 was sampled equally in every province, including
  winter in temperate and subpolar provinces (colder, more productive water).
  Only the polar provinces were excluded. Hence an average SST of 20.8 C.
- Our samples are the opposite selection: Bio-GO-SHIP transects through
  tropical and subtropical open ocean, and only tiles at least 80% cloud-free,
  which favours the clear skies over the subtropical gyres.

So relative to what Granite saw most of, our samples are an under-represented
corner: clear, warm, oligotrophic water.

**How the reconstruction numbers are measured.** Granite was pretrained as a
masked autoencoder: 75% of a tile's 2x2-pixel patches are hidden and the
decoder predicts them from the rest. `check_granite_inputs.py` does the same on
real tiles and scores only hidden pixels that are clear water with SST present,
as a root-mean-square error in normalised units (1 = one Granite s.d.):

- Granite's decoder: 1.204.
- Constant guess, predicting Granite's average (0) for every hidden pixel: 1.330.
- (1.330 - 1.204) / 1.330 = 9.5%, hence "about 10% better than a constant guess".

The constant guess is a low bar for our tiles, because they sit far from
Granite's average. The fairer yardstick is predicting each hidden pixel from the
tile's own visible pixels (their mean, per channel). Since a tile varies by only
about 0.4 s.d. internally, that simple guess should do far better than 1.2.
`check_granite_inputs.py` now reports it, with the per-channel bias, but it has
not been run on real tiles yet. On synthetic tiles with our open-ocean spectrum,
the decoder returned blue bands 1.5 s.d. too low and SST 0.75 s.d. too low: it
pulls water like ours back toward its training average instead of copying the
level it can see.

This concerns the decoder, which is not used for embeddings. (An earlier version
of this note concluded that the encoder therefore represents open-ocean water
coarsely. The follow-up below shows that conclusion was wrong.)

### Follow-up (2026-10-21): is the decoder result a preprocessing problem?

A masked autoencoder that sees 25% of a tile should copy the visible level
almost trivially, so a reconstruction error of 1.2 when tiles vary by 0.4
internally pointed at the inputs. Checked in order (local tests with the real
checkpoint; pseudo-tiles built from the band means of 96 real tiles, with mild
texture added):

1. **Fill values.** Flagged pixels (cloud, glint and the other flags) are not
   filled: they keep their measured values. Only OLCI fill pixels (1.4% of
   pixels) and missing SST (24% of SST pixels, in 71 of 300 tiles) are set to 0,
   Granite's average. Reconstruction is scored only on clear pixels with SST.
   Setting flagged pixels to the average as well changed the error from 1.204 to
   1.207 and the embeddings hardly at all (median cosine 0.995). Not the cause.
2. **Product, units, band order.** Our input is the `OaNN_reflectance` variables
   of OLCI L2 WFR, which are water-leaving reflectance (dimensionless), the
   product the paper names ("L2 OLCI full resolution water reflectances"), with
   Granite's band names being these exact variable names. Granite's mean at
   400 nm (0.0235) is plausible as water-leaving reflectance and implausibly
   bright as remote-sensing reflectance (which is 1/pi of it, typically 0.01 to
   0.015 even in clear water). Band order matches config.yaml exactly. Feeding
   the pseudo-tiles as Rrs, as pi x reflectance, with SST in Celsius, with the
   OLCI bands reversed or with SST first: none reconstructs better than our
   convention (errors 1.0 to 24 against our 1.01).
3. **The decoder itself.** Copying each pseudo-tile's visible pixels scores 0.05;
   the decoder scores 1.0. It scores just as badly on random spectra of ordinary
   size (each channel drawn from N(0, 1): 0.95), so this is not about our water
   being unusual. It reconstructs well only when all channels share one level:
   Granite's average everywhere 0.046, every channel at +1 0.040, but eight
   channels at +1 and nine at -1 gives 0.94. At a hidden pixel, the spread across
   the 17 channels is 0.97 in the input and 0.12 in the decoder output. The
   released decoder reproduces a shared level and spatial texture but not the
   differences between channels, whatever the input convention. Its
   reconstruction error therefore cannot diagnose our preprocessing.
   (Checked and ruled out: a different patch-unpacking order, and a decoder
   trained to predict per-patch standardised pixels.)
4. **What matters for us: does the encoder keep the absolute spectrum?** A
   cross-validated linear readout (ridge regression, folds grouped by cruise)
   from each tile's 512-number embedding to its own mean value per channel,
   over all 3,503 tiles: R2 = 0.99 for Oa01-Oa04, 0.96-0.97 for Oa05-Oa06,
   0.83-0.92 for Oa07-Oa11, 0.95 for Oa21 and 0.97 for SST. Only Oa12-Oa18, the
   near-infrared bands where tiles hardly differ, are poorly recovered (0.07 to
   0.53). The embeddings retain the absolute colour and temperature of the tile,
   including for our bright-blue, warm water.

Conclusion: the inputs are what Granite expects, and the embeddings carry the
spectral levels. The poor reconstruction is a property of the released decoder,
which only matters if the decoder is used (for example to measure
reconstruction, or to continue pretraining). Worth raising with IBM on the model
page. `check_granite_inputs.py` now also reports the channel spread of input and
decoder output, to confirm item 3 on real tiles.

**A finding about the model rather than the inputs.** Granite reconstructs our
tiles only about 10% better than a constant guess. Probing the decoder with
synthetic tiles (local test, real checkpoint): a uniform tile at +0.5 or +1.0 s.d.
in every channel is reproduced exactly, but a uniform tile with the open-ocean
spectrum above is pulled toward the training average: hidden pixels come back
1.5 s.d. too low in Oa01-Oa02 and 0.75 s.d. too low in SST. So Granite represents
bright-blue, warm, oligotrophic water poorly: the kind of water most Bio-GO-SHIP
samples come from. Its embeddings may compress the open-ocean variability that
matters here. This makes the raw-channel baseline comparison (section 8)
important, and makes fine-tuning Granite on this data worth considering later.
`check_granite_inputs.py` now also reports a "tile's own visible-pixel mean"
baseline and the per-channel bias, to quantify this on real tiles at the next run.

## 8. First check: do the features carry signal?

`scripts/compare_feature_sets.py`, 2026-10-20. Targets: in-situ temperature,
phosphate and nitrate + nitrite from the sample metadata (NMF and KEGG not yet
available). 610 samples with both a fixed and an advected tile; one row per
sample, using each position's best tile. CatBoost (600 trees, depth 6, learning
rate 0.05). R2 against the training-fold mean.

Spatial-block cross-validation (5 folds; blocks = section x 5-degree latitude band):

| Target | Samples | Fixed | Advected | Metadata (lat, lon, day of year) |
|---|---:|---:|---:|---:|
| Temperature | 587 | 0.77 +/- 0.08 | 0.77 +/- 0.07 | 0.93 +/- 0.02 |
| Phosphate | 302 | 0.68 +/- 0.12 | 0.65 +/- 0.14 | 0.92 +/- 0.06 |
| Nitrate + nitrite | 180 | 0.71 +/- 0.26 | 0.73 +/- 0.20 | 0.92 +/- 0.05 |

Holding out one cruise at a time (temperature, 5 folds):

| Held out | Fixed | Advected | Metadata |
|---|---:|---:|---:|
| AMT28 | 0.71 | 0.71 | 0.67 |
| C13.5/A13.5 | 0.65 | 0.68 | 0.51 |
| I07N | 0.69 | 0.65 | 0.84 |
| I09N | 0.82 | 0.82 | 0.81 |
| P18 | 0.72 | 0.68 | 0.79 |
| **Mean** | **0.72 +/- 0.06** | **0.71 +/- 0.06** | **0.73 +/- 0.12** |

Reading:

- The pipeline works end to end: Granite embeddings of the tiles explain about
  three quarters of the variance in in-situ temperature, and two thirds or more
  for the nutrients, on samples from parts of the ocean the model was not
  trained on.
- Fixed and advected are indistinguishable on these targets; the differences
  are well inside the fold-to-fold spread. Expected: the water moved about one
  tile width, the two tiles are often the same window, and the embeddings of the
  two positions are nearly identical (cosine 0.999).
- Position and season alone predict these targets very well when neighbouring
  stretches of the same transect are in the training data (0.92 to 0.93),
  because temperature and nutrients vary smoothly along a transect. When a whole
  cruise is held out, metadata drops to the level of the embeddings and becomes
  less stable from cruise to cruise, while the embeddings hold up. Held-out
  cruises are the fairer test for transfer to new regions.
- Raw-channel baseline added 2026-10-20 (below): Granite's embeddings beat the
  tile's own band averages clearly.
- These are physical targets with a direct satellite signature (temperature has
  its own SST channel). Microbial community targets are the real test.

**With the raw-channel baseline (2026-10-20).** `raw_fixed` and `raw_advected` are
the mean and s.d. of each of the 17 normalised channels over the tile's clear
pixels (34 numbers; `scripts/pool_tile_features.py`). Same 610 samples; tiles
without SST keep a gap that CatBoost handles. R2:

| Target | CV | Fixed | Raw fixed | Advected | Raw advected | Metadata |
|---|---|---:|---:|---:|---:|---:|
| Temperature | blocks | 0.77 | 0.70 | 0.77 | 0.72 | 0.93 |
| Phosphate | blocks | 0.68 | 0.36 | 0.65 | 0.39 | 0.92 |
| Nitrate + nitrite | blocks | 0.71 | 0.04 | 0.73 | -0.02 | 0.92 |
| Temperature | cruise held out | 0.72 +/- 0.06 | 0.61 +/- 0.18 | 0.71 +/- 0.06 | 0.63 +/- 0.14 | 0.73 +/- 0.12 |
| Phosphate | cruise held out | 0.40 +/- 0.24 | -0.04 +/- 0.83 | 0.34 +/- 0.40 | -0.08 +/- 0.76 | 0.51 +/- 0.41 |

- Granite's embeddings carry clearly more usable information than the band
  averages of the same tiles, most of all for the nutrients, and they transfer
  to unseen cruises more steadily. So despite reconstructing open-ocean spectra
  imperfectly (section 7), the encoder adds something: likely spatial texture
  within the tile that averages discard.
- One confound to rule out: the embeddings also see flagged pixels (up to 20% of
  a tile) that the raw averages exclude, so cloud and glint patterns, which vary
  with region, could carry some of the signal. A raw baseline that includes the
  clear fraction and flag fractions would test this.
- Fixed and advected remain indistinguishable on these targets.

### Missing SST (2026-10-21)

**Where it is missing.** SST is all or nothing per tile: 861 of 3,503 tiles (25%)
have no SST at all, the rest have it everywhere. (The nearest-neighbour
resampling fills every pixel from the nearest acceptable SLSTR pixel within the
search window, so a tile is empty only if the window holds no acceptable SST.)
The gaps are concentrated:

| Section | Tiles | Without SST |
|---|---:|---:|
| I07N | 1,086 | 493 (45%) |
| C13.5/A13.5 | 1,549 | 252 (16%) |
| AMT28 | 372 | 74 (20%) |
| P18 | 280 | 26 (9%) |
| I09N | 216 | 16 (7%) |

The main cause is the SST pairing. Where a scene's own satellite had no SST
product (Sentinel-3B in spring 2018, during I07N, and in autumn 2018, during
AMT28), the other satellite's orbit file was used; its swath and timing differ,
and 52% of those tiles end up without SST, against 20% with a same-satellite
pair. I07N Sentinel-3B scenes: 75% without SST. Fixed and advected tiles are
affected equally (24% and 25%).

**The embeddings record it strongly.** A linear classifier tells no-SST tiles
from the rest by their embedding with AUC 0.99 (cruises held out). Along that
direction, no-SST tiles sit 1.9 s.d. from the others (3.6 s.d. within I07N
alone). The rest of the content survives: blue and green bands are recovered from
the embedding equally well with or without SST (R2 0.99 and 0.96-0.97 in both).
So "SST missing" is a large, easily detected signal in the embeddings, and it is
regional, which makes it a possible confound for any target that differs by
region.

**Effect on predictions** (CatBoost, spatial-block CV, folds now balanced so each
test fold holds 23-25% of samples without SST for temperature, 22-38% for
phosphate, the closest whole blocks allow; `compare_feature_sets.py`):

| Target | Features | RMSE with SST | RMSE without SST |
|---|---|---:|---:|
| Temperature | Fixed embedding | 1.8 C | 3.7 C |
| Temperature | Advected embedding | 2.1 C | 3.5 C |
| Temperature | Metadata | 1.7 C | 1.3 C |
| Phosphate | Fixed embedding | 0.17 | 0.19 umol/kg |

Missing SST roughly doubles the temperature error from embeddings; phosphate is
much less affected. Holding out whole cruises cannot be balanced (I07N carries
57% of the no-SST tiles).

**Filling.** Missing SST is currently set to Granite's average, 20.8 C, while our
tiles with SST have a median of 25.6 C (23 to 31 C by cruise). For a whole tile
this is not a small gap to patch: it replaces the tile's temperature with a value
that is wrong by about 5 C on average, and it is the same value everywhere, which
is how the embedding recognises it. Options, in order of preference:

1. Fill from a gap-free SST analysis at the scene's date and the tile's position,
   e.g. the Copernicus Marine OSTIA L4 product (0.05 degree, daily; the same
   account as the currents). It blends SLSTR with other sensors, so it is close to
   what Granite was trained on, and differs from the tile's true SST by a fraction
   of a degree in clear-sky conditions.
2. Fill from the nearest tile of the same cruise and day that has SST.
3. Keep the current fill but give downstream models an explicit SST-missing
   indicator, so the confound is visible rather than hidden.

In-situ temperature must not be used for filling: it is a target.

### Ways to improve Granite for this data

In rough order of cost:

1. **Align the quality mask with pretraining.** The paper masked the flags
   recommended by the Sentinel-3 product notice. Our mask, inherited from the
   Florida project, does not reject HIGHGLINT, which is set on many of our
   pixels (on average 18% of a rejected window). Glint brightens all bands; if
   pretraining masked it and we do not, that adds a mismatch. Check the
   product-notice list, record per-pixel flag fractions per tile, and compare
   embeddings with and without glint pixels.
2. **Better readout of the encoder.** The embedding is the mean of the final
   layer's tokens, which discards spatial structure and earlier layers. Try the
   mean of several layers, the CLS token, the tokens nearest the tile centre, and
   embeddings concatenated with the raw band averages (which carry the absolute
   levels the model compresses).
3. **Continued self-supervised pretraining on open-ocean tiles.** The 527 OLCI
   scenes on disk hold many more clear tiles than the sample sites, almost all
   tropical and subtropical open ocean. Continuing Granite's masked-autoencoder
   training on them (no labels needed) adapts it to this water. Because the
   released decoder does not reconstruct spectra (follow-up above), this would
   in effect retrain the decoder too; judge it by the downstream comparison, not
   by reconstruction error alone. Lower priority now that the encoder is known
   to keep the spectral levels.
4. **Fine-tuning on the targets**, once NMF and KEGG are available: small sample
   sizes, so only with strong regularisation and cruise-held-out validation.

## 9. Planned next

- **Fixed vs advected.** Compare models trained on the two sets of embeddings,
  using the same overpass for both (`same_scene_as_best_fixed`).
- **Short time series per sample.** Average the embeddings of roughly 3 to 5
  clear scenes along the drift track, so each sample is represented by more than
  one view of its water. The embeddings for all clear overpasses within +/-3 days
  are already produced; scenes at earlier points on the track would need a
  further search.
- **Advected position outside the scene** (expected near coasts). Two responses,
  to be built once `site_tile_status.csv` shows how often it happens: (a) search
  for and download the scene that does contain the advected position; (b) switch
  to another overpass, for example one before sampling if the one after failed.
- **Models.** Embeddings -> CatBoost predicting NMF community weights and KEGG
  profiles.

## Open questions

- Keep the 61 I09N sites matched at 3 to 10 days after sampling, or drop I09N?
- Depth cut-off for "surface" samples (a few are from 20 to 80 m).
- Back-tracking: 15 m or 0 m currents; daily or hourly; add diffusion or an ensemble?
- Tiles along the trajectory (water-mass history), and how far back to go.
