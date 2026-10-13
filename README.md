# biogoship_granite_kegg

Predicting marine microbial community structure (NMF communities) and metabolic
potential (KEGG) from space, using Bio-GO-SHIP eDNA samples, Sentinel-3 ocean
colour and sea surface temperature, and the
[IBM Granite Geospatial Ocean](https://huggingface.co/ibm-granite/granite-geospatial-ocean)
foundation model.

This file is the how-to: what each step does and how to run it. The reasoning
behind the choices, and the results so far, are in
[METHODS_RESULTS.md](METHODS_RESULTS.md).

## Pipeline

| Step | What it does | Script | Status |
|---|---|---|---|
| 1 | Build the sample and site manifests from public metadata | `scripts/build_sample_manifest.py` | done, output committed |
| 2 | Find the Sentinel-3 scenes covering each site | `scripts/search_products.py` | done, output committed |
| 3 | Download those scenes (OLCI + SST) | `scripts/download_products.py` | ready to run on the cluster |
| 4 | Back-track each site's water to the time of every overpass | `scripts/download_currents.py`, `scripts/backtrack_sites.py` | ready to run on the cluster |
| 5 | Cut Granite-ready tiles at the sampling position and at the back-tracked position | `scripts/download_tiles.py`, `scripts/build_tile_manifest.py` | written, not run yet |
| 6 | Tiles to Granite embeddings | | planned |
| 7 | Embeddings to NMF communities and KEGG profiles (CatBoost) | | planned |

## Repository layout

```
biogoship/            library code
  cdse.py             Copernicus login, catalogue search, downloads
  tiles.py            Sentinel-3 readers, quality masks, Granite normalisation
  advection.py        back-tracking with OceanParcels
  tables.py           CSV, time and distance helpers
scripts/              one script per pipeline step (see table above)
cluster/              SLURM job scripts and shared environment setup
data/manifests/       small tables, committed
tests/                offline tests
METHODS_RESULTS.md    methods notes and results log
```

Committed manifests in `data/manifests/`:

| File | One row per | Contents |
|---|---|---|
| `biogoship_samples.csv` | BioSample | position, time, section, depth, run accessions, `site_id` |
| `biogoship_sites.csv` | site (unique position and time) | the unit matched to satellite data |
| `site_products.csv` | (site, OLCI scene) | every candidate overpass and its time gap in hours |
| `products.csv` | OLCI scene | the scene and its paired SST orbit file |
| `advected_positions.csv` (after step 4) | (site, OLCI scene) | where the site's water was at the scene's sensing time |

## Requirements

Python 3.10 or newer. Steps 1-3 and 5 need only `numpy`, `h5py`, `certifi` and
`pandas`; step 4 also needs `xarray`, `netCDF4`, OceanParcels 3.x and the
Copernicus Marine toolbox.

```bash
pip install -r requirements.txt
python -m unittest discover -s tests     # offline tests, no network needed
```

Two free accounts are needed, both entered in `cdse_credentials.sh` (copy it from
`cdse_credentials.sh.example`):

- [Copernicus Data Space](https://dataspace.copernicus.eu) for Sentinel-3 (steps 3 and 5)
- [Copernicus Marine](https://marine.copernicus.eu) for surface currents (step 4)

Searching (step 2) needs neither.

## Steps 1 and 2: manifests and scene search

Already run; the outputs are committed. To rebuild them (any machine with
internet, no credentials, about five minutes):

```bash
python scripts/build_sample_manifest.py
python scripts/search_products.py --max-days 3 --section-max-days I09N=10
```

`--max-days` is the largest allowed gap between sampling and a satellite overpass.
`--section-max-days` widens it for one section.

## Step 3: download the Sentinel-3 scenes

Downloads every candidate scene, not just the closest in time, so that a cloudy
overpass can be replaced by a clear one: 414 OLCI scenes and 218 SST orbit files,
about 300 GB, as complete product archives under `data/products/`.

On the cluster (SLURM). Jobs that need internet use the `service` partition.
Always `cd` into the repo before `sbatch`; the job scripts find the repo through
`$SLURM_SUBMIT_DIR`.

```bash
# once, on the login node
cd /scratch/morrill/users/$USER
git clone https://github.com/helloftroy/biogoship_granite_kegg.git
cd biogoship_granite_kegg
cp cdse_credentials.sh.example cdse_credentials.sh    # then edit: Copernicus login
chmod 600 cdse_credentials.sh
df -h .                                               # needs about 300 GB free

# smoke test: the first two scenes and their SST files (about 2 GB)
sbatch --account=<allocation> --export=ALL,MAX_PRODUCTS=2 cluster/download_products.sbatch
tail -f logs/bgs_products_<jobid>.out

# everything; resubmit the same command to resume or to retry failures
sbatch --account=<allocation> cluster/download_products.sbatch
```

Outputs:

- `data/products/olci/<product>.zip` and `data/products/sst/<product>.zip`
- `data/download_log.csv`, one line per download attempt

Options, passed with `--export=ALL,NAME=value`:

| Name | Default | Meaning |
|---|---|---|
| `MAX_PRODUCTS` | 0 (all) | only the first N OLCI scenes |
| `NO_SST` | 0 | set to 1 to skip the SST files |
| `CONNECTIONS` | 4 | parallel downloads; 4 is the Copernicus per-account limit |
| `BGS_CONDA_ENV` | `/scratch/morrill/users/$USER/conda_envs/faire-agent` | conda environment, by full path (step 4 defaults to `.../conda_envs/biogoship`) |

Anywhere else:

```bash
source cdse_credentials.sh
python scripts/download_products.py --dry-run     # counts and sizes only
python scripts/download_products.py
```

`cluster/env.sh` is shared by all jobs. It keeps caches off `$HOME`, activates the
conda environment by full path and checks that the activation took effect, and
installs `numpy`/`h5py` into the environment if they are missing.

## Step 4: back-track the sampled water

A satellite overpass is hours to days away from sampling, and the water moves in
between. This step releases one particle per site at the sampling position and
time and advects it with surface currents, backward and forward, using
[OceanParcels](https://docs.oceanparcels.org). Each (site, scene) pair then gets
the position of that water at the scene's sensing time.

Currents: Copernicus Marine
[MULTIOBS_GLO_PHY_MYNRT_015_003](https://data.marine.copernicus.eu/product/MULTIOBS_GLO_PHY_MYNRT_015_003/description)
(GlobCurrent total velocity, 0.25 degree, daily), about 100 MB for all sections.

```bash
# once: a project environment with OceanParcels and the Copernicus Marine toolbox
sbatch --account=<allocation> cluster/setup_env.sbatch

# add your Copernicus Marine login to cdse_credentials.sh, then
sbatch --account=<allocation> cluster/backtrack.sbatch
```

Outputs:

- `data/currents/currents_<section>.nc` the current fields
- `data/advection/trajectories.csv` hourly track for every site
- `data/manifests/advected_positions.csv` one row per (site, scene): advected
  position, straight-line displacement and path length in km, and a `status`

Options (`--export=ALL,NAME=value`): `DEPTH_M` (15, or 0 for the surface level),
`MARGIN_DEG` (5; room around the sites for water to drift into).

Anywhere else:

```bash
source cdse_credentials.sh
python scripts/download_currents.py
python scripts/backtrack_sites.py --depth-m 15
```

## Step 5: cut tiles (not run yet)

For each site and each clear overpass, writes a `(17, 42, 42)` array: 16 OLCI
bands plus SST, normalised for Granite. Two versions are cut from the same scene:

- `fixed`: centred on the sampling position, as in the Florida project
- `advected`: centred on where the sampled water was at the scene's sensing time (step 4)

It reads the archives from step 3 when they are on disk, and otherwise downloads
just the files it needs.

```bash
sbatch --account=<allocation> cluster/download_tiles.sbatch
sbatch --account=<allocation> --export=ALL,POSITIONS=fixed cluster/download_tiles.sbatch   # without step 4
```

Outputs:

- `data/tiles/YYYY/MM/tile_<satellite>_<time>_<site>_<position>.npz`
- `data/tile_log.csv` outcome for every site, scene and position
- `data/manifests/biogoship_tile_manifest.csv` samples joined to tiles. `position`
  is `fixed` or `advected`; `is_best` marks the tile closest in time for each
  sample and position; `same_scene_as_best_fixed` marks the advected tile from the
  same overpass as the best fixed tile.

## Data sources and citation

- Bio-GO-SHIP metagenomes and amplicons: Larkin et al. 2021, *Scientific Data* 8, 107,
  [doi:10.1038/s41597-021-00889-9](https://doi.org/10.1038/s41597-021-00889-9);
  BioProject [PRJNA656268](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA656268).
- Curated sample metadata:
  [merenlab/public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata).
- Sentinel-3 OLCI and SLSTR level-2 products: Copernicus Data Space Ecosystem.
- Surface currents: Copernicus Marine Service, MULTIOBS_GLO_PHY_MYNRT_015_003,
  [doi:10.48670/mds-00327](https://doi.org/10.48670/mds-00327).
- Particle tracking: [OceanParcels](https://oceanparcels.org) v3.
- The Sentinel-3 handling is carried over from
  [helloftroy/multimodal_eDNA](https://github.com/helloftroy/multimodal_eDNA).
