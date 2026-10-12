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
| 4 | Cut Granite-ready tiles around each site | `scripts/download_tiles.py`, `scripts/build_tile_manifest.py` | written, not run yet |
| 5 | Water-mass back-trajectories from each site | | planned |
| 6 | Tiles to Granite embeddings | | planned |
| 7 | Embeddings to NMF communities and KEGG profiles (CatBoost) | | planned |

## Repository layout

```
biogoship/            library code
  cdse.py             Copernicus login, catalogue search, downloads
  tiles.py            Sentinel-3 readers, quality masks, Granite normalisation
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

## Requirements

Python 3.9 or newer with `numpy`, `h5py`, `certifi` and (for step 1 only) `pandas`:

```bash
pip install -r requirements.txt
python -m unittest discover -s tests     # offline tests, no network needed
```

Downloading needs a free [Copernicus Data Space](https://dataspace.copernicus.eu)
account. Searching does not.

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
| `BGS_CONDA_ENV` | `/scratch/morrill/users/$USER/conda_envs/faire-agent` | conda environment, by full path |

Anywhere else:

```bash
source cdse_credentials.sh
python scripts/download_products.py --dry-run     # counts and sizes only
python scripts/download_products.py
```

`cluster/env.sh` is shared by all jobs. It keeps caches off `$HOME`, activates the
conda environment by full path and checks that the activation took effect, and
installs `numpy`/`h5py` into the environment if they are missing.

## Step 4: cut tiles (not run yet)

For each site and each clear overpass, writes a `(17, 42, 42)` array centred on
the site: 16 OLCI bands plus SST, normalised for Granite. It reads the archives
from step 3 when they are on disk, and otherwise downloads just the files it needs.

```bash
sbatch --account=<allocation> cluster/download_tiles.sbatch
```

Outputs: `data/tiles/YYYY/MM/*.npz`, `data/tile_log.csv` (outcome for every site
and scene) and `data/manifests/biogoship_tile_manifest.csv` (samples joined to
tiles, with `is_best` marking the tile closest in time).

## Data sources and citation

- Bio-GO-SHIP metagenomes and amplicons: Larkin et al. 2021, *Scientific Data* 8, 107,
  [doi:10.1038/s41597-021-00889-9](https://doi.org/10.1038/s41597-021-00889-9);
  BioProject [PRJNA656268](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA656268).
- Curated sample metadata:
  [merenlab/public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata).
- Sentinel-3 OLCI and SLSTR level-2 products: Copernicus Data Space Ecosystem.
- The Sentinel-3 handling is carried over from
  [helloftroy/multimodal_eDNA](https://github.com/helloftroy/multimodal_eDNA).
