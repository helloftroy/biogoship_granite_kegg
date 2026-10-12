# biogoship_granite_kegg

Predicting marine microbial community structure (NMF communities) and metabolic
potential (KEGG) from space, using Bio-GO-SHIP eDNA samples, Sentinel-3 ocean
colour and sea surface temperature, and the
[IBM Granite Geospatial Ocean](https://huggingface.co/ibm-granite/granite-geospatial-ocean)
foundation model.

Planned pipeline:

1. **Samples to satellite tiles** (this repo, working): match every Bio-GO-SHIP sample
   to Sentinel-3 overpasses and cut a Granite-ready tile around it.
2. Tiles to Granite embeddings (to come).
3. Embeddings to NMF communities and KEGG profiles with CatBoost (to come).

## Data

**Bio-GO-SHIP** ([Larkin et al. 2021](https://doi.org/10.1038/s41597-021-00889-9),
BioProject [PRJNA656268](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA656268)):
996 surface-ocean samples from eight cruises, 971 with a metagenome and 975 with
amplicon sequencing. Sample positions and times come from the ENA run report;
depth, sampling method, temperature and nutrients come from the curated table in
[merenlab/public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata).

**Sentinel-3** (Copernicus Data Space): OLCI level-2 water reflectance
(`OL_2_WFR`, 300 m) and SLSTR sea surface temperature (`SL_2_WST`).

Only cruises after Sentinel-3A began delivering data (late April 2016) can be matched:

| Section | Dates | Samples | Sites | Sites with a candidate overpass (±3 days) |
|---|---|---:|---:|---:|
| C13.5/A13.5 (Atlantic) | Mar-Apr 2020 | 263 | 234 | 234 |
| I07N (Indian) | Apr-Jun 2018 | 261 | 259 | 259 |
| P18 (Pacific) | Nov 2016-Jan 2017 | 104 | 104 | 104 |
| AMT28 (Atlantic) | Sep-Oct 2018 | 63 | 63 | 63 |
| I09N (Indian) | Mar-Apr 2016 | 242 | 242 | 7 |
| NH1418, AE1319, BVAL46 | 2011-2014 | 63 | 52 | 0 |

A *site* is a unique position and time; replicate samples share a site. The 667
matchable sites are covered by 407 OLCI scenes. How many sites end up with a
*cloud-free* tile is only known after the download.

## Layout

```
biogoship/            library code
  cdse.py             Copernicus login, catalogue search, per-file downloads
  tiles.py            Sentinel-3 readers, quality masks, Granite normalisation
  tables.py           CSV, time and distance helpers
scripts/
  build_sample_manifest.py   ENA + curated metadata -> samples and sites
  search_products.py         sites -> candidate OLCI scenes + paired SST
  download_tiles.py          scenes -> one tile per (site, clear overpass)
  build_tile_manifest.py     samples x tiles table for modelling
cluster/              SLURM jobs for the morrill cluster
data/manifests/       small tables, committed
tests/                offline tests (python -m unittest discover -s tests)
```

Committed manifests:

| File | Rows | Contents |
|---|---|---|
| `biogoship_samples.csv` | one per BioSample | position, time, section, depth, run accessions, `site_id` |
| `biogoship_sites.csv` | one per site | the unit matched to satellite data |
| `site_products.csv` | one per (site, OLCI scene) | candidate overpasses and time gap in hours |
| `products.csv` | one per OLCI scene | scene and its paired SST orbit file |

## How a tile is made

For each OLCI scene that covers at least one site within ±3 days of sampling:

1. Download only `geo_coordinates.nc` and `wqsf.nc` (about 60 MB) and find each
   site's nearest pixel.
2. Take the 42x42 pixel window (12.6 km) centred on that pixel. Keep it if at least
   80% of its pixels are clear open ocean (water, and not land, cloud, cloud
   margin, ice, coastline, suspect, saturated or failed atmospheric correction).
3. Only if some site passes, download the 16 reflectance bands (about 280 MB) and
   the paired SST orbit file (about 600 MB, cached while neighbouring scenes share it).
4. Save `data/tiles/YYYY/MM/tile_<satellite>_<sensing time>_<site_id>.npz` with
   `data` of shape `(17, 42, 42)`: 16 OLCI bands (Oa01-Oa12, Oa16, Oa17, Oa18, Oa21)
   and SST, z-scored with the Granite pretraining statistics. SST is NaN where
   SLSTR has no acceptable-quality pixel.

Every clear overpass is kept, so a site can have several tiles;
`build_tile_manifest.py` flags the one closest in time as `is_best`.

Checking for cloud before fetching the bands means cloudy scenes cost 60 MB
instead of 300-600 MB each. If per-file downloads misbehave, `--download-mode zip`
falls back to whole product archives.

## Running on the cluster

Jobs that need internet run on the `service` partition. Always `cd` into the repo
before `sbatch`: the job scripts locate the repo through `$SLURM_SUBMIT_DIR`.

```bash
# once, on the login node
cd /scratch/morrill/users/$USER
git clone https://github.com/helloftroy/biogoship_granite_kegg.git
cd biogoship_granite_kegg
cp cdse_credentials.sh.example cdse_credentials.sh    # then edit: Copernicus login
chmod 600 cdse_credentials.sh

# smoke test: three scenes
sbatch --account=<allocation> --export=ALL,MAX_PRODUCTS=3 cluster/download_tiles.sbatch
tail -f logs/bgs_tiles_<jobid>.out

# full run (resubmit the same command to resume or to retry failed scenes)
sbatch --account=<allocation> cluster/download_tiles.sbatch
```

Outputs, all under the repo on scratch:

- `data/tiles/` the tiles
- `data/tile_log.csv` outcome for every (site, scene): `valid`, `cloudy`, `edge` or `out_of_swath`
- `data/product_log.csv` one line per processed scene
- `data/manifests/biogoship_tile_manifest.csv` samples joined to tiles, with a coverage summary printed at the end of the job log

`cluster/env.sh` activates the conda environment by full path
(`/scratch/morrill/users/$USER/conda_envs/faire-agent` by default) and installs
`numpy`/`h5py` into it if they are missing. Point it elsewhere with
`--export=ALL,BGS_CONDA_ENV=/path/to/env`.

Settings, passed the same way: `MIN_VALID` (default 0.8), `DOWNLOAD_MODE`
(`files` or `zip`), `NO_SST=1`, `CONNECTIONS` (default 4, the Copernicus limit).

## Rebuilding the manifests

The search needs no credentials and runs anywhere with internet:

```bash
python scripts/build_sample_manifest.py
python scripts/search_products.py --max-days 3
```

## Provenance

The Sentinel-3 handling (band selection, quality flags, normalisation, SST
pairing and resampling) is carried over from the Florida eDNA/fish project,
[helloftroy/multimodal_eDNA](https://github.com/helloftroy/multimodal_eDNA). That
project tiled whole scenes over one region; here each site gets its own centred
tile, since the samples are spread along global transects.
