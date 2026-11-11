#!/bin/bash
#SBATCH --job-name=bgs_sst_l4
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Inputs for the SST comparison (RESULTS.md, test 3), CPU only, needs internet:
#  1. OSTIA L4 SST per cruise section (copernicusmarine, same login as the currents),
#     then copies of every tile with SST gaps filled from OSTIA shifted toward the
#     daytime SLSTR skin SST (scripts/fill_sst_l4.py). Originals are not touched.
#     Original cruises: data/tiles_sst_l4/, manifest data/manifests/biogoship_tile_manifest_sst_l4.csv
#     FAIRe cruises:    data/faire/tiles_sst_l4/, manifest data/faire/manifests/biogoship_tile_manifest_sst_l4.csv
#  2. Band means for the FAIRe tiles (data/faire/embeddings/tile_raw_features.csv), for the
#     temperature + ocean-colour baseline on all 13 cruises.
# Re-running skips OSTIA files and tiles already written. Then: cluster/embed_sst_variants.sbatch.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/sst_variants.sbatch

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"

if [ -f cdse_credentials.sh ]; then
  source cdse_credentials.sh
fi
if [ -z "${COPERNICUSMARINE_SERVICE_USERNAME:-}" ] || [ -z "${COPERNICUSMARINE_SERVICE_PASSWORD:-}" ]; then
  echo "No Copernicus Marine credentials: add them to cdse_credentials.sh (see cdse_credentials.sh.example)." >&2
  exit 1
fi
# Keep the toolbox's own files off $HOME.
export COPERNICUSMARINE_CREDENTIALS_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
export COPERNICUSMARINE_CACHE_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
mkdir -p "${COPERNICUSMARINE_CREDENTIALS_DIRECTORY}"

echo "Job started: $(date)"
if [ ! -s data/faire/embeddings/tile_raw_features.csv ]; then
  python scripts/pool_tile_features.py --manifest-csv data/faire/manifests/biogoship_tile_manifest.csv \
    --out-csv data/faire/embeddings/tile_raw_features.csv
fi
python scripts/fill_sst_l4.py
python scripts/fill_sst_l4.py --manifest-csv data/faire/manifests/biogoship_tile_manifest.csv \
  --tiles-root data/faire/tiles --out-root data/faire/tiles_sst_l4 \
  --out-manifest data/faire/manifests/biogoship_tile_manifest_sst_l4.csv --log-csv data/faire/sst_l4_log.csv
echo "Job finished: $(date)"
