# Methods and results log

Working notes for the Bio-GO-SHIP x Sentinel-3 x Granite project: what was done,
why, and what came out. Running instructions are in [README.md](README.md).
Newest results are added at the end of each section, with the date.

## 1. Aim

Predict surface-ocean microbial community structure (NMF communities from eDNA)
and metabolic potential (KEGG profiles from metagenomes) from satellite data:
Sentinel-3 tiles -> Granite Geospatial Ocean embeddings -> CatBoost -> NMF / KEGG.

New relative to the Florida project (`multimodal_eDNA`): samples are global
transects rather than one region, and each sample will also get satellite views
along its water mass's back-trajectory, not just a snapshot at the sampling site.

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

**Result:** not run yet.

## 5. Tiles (written, not run)

Carried over unchanged from the Florida project: band selection, quality flags,
normalisation, SST resampling. Changed: the tile is centred on the site instead
of being cut from a fixed grid over a region.

For each site and candidate scene (`scripts/download_tiles.py`):

1. Find the OLCI pixel nearest the site. If it is more than 1 km away, the site
   is outside the swath (`out_of_swath`).
2. Take the 42x42 pixel window (12.6 km) centred on that pixel. If it runs off
   the scene, skip (`edge`).
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

Every clear overpass yields a tile, so a site can have several.
`scripts/build_tile_manifest.py` joins samples to tiles and marks the one closest
in time as `is_best`.

**Result:** not run yet.

## 6. Water-mass back-trajectories (planned)

Idea: track the sampled water backwards from each site and collect satellite
views along the trajectory, giving a short history of the water mass in addition
to the snapshot at the sampling site. Method to be decided.

## 7. Embeddings and models (planned)

Tiles -> Granite Geospatial Ocean encoder -> embeddings -> CatBoost predicting
NMF community weights and KEGG profiles.

## Open questions

- Keep the 61 I09N sites matched at 3 to 10 days after sampling, or drop I09N?
- Depth cut-off for "surface" samples (a few are from 20 to 80 m).
- Back-trajectory method, and how far back in time to go.
