#!/bin/bash
#SBATCH --job-name=bgs_ext_tiles
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Sentinel-3 tiles for the external test sets (RESULTS.md R13, R14), one after another:
#   data/external/jgi_flinkstrom  data/external/tara_pacific  data/external/emo_bon
# Manifests come from scripts/build_external_manifests.py (in git). For each set:
# scene search -> tiles at the sampling position (fixed only) -> tile manifest ->
# OSTIA-filled copies (tiles_sst_l4, the standard since R12) -> band means of those.
# Then cluster/external_embed.sbatch (GPU). Resumable: re-submit after a stop.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/external_tiles.sbatch
#   sbatch --account=191001-364393 --export=ALL,SETS="tara_pacific" cluster/external_tiles.sbatch
#
# Overrides: SETS, MAX_DAYS=3, PAUSE=2, MIN_VALID=0.8, MAX_OFFSET_KM=20.

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"

if [ -f cdse_credentials.sh ]; then
  source cdse_credentials.sh
fi
if [ -z "${CDSE_ACCESS_TOKEN:-}" ] && { [ -z "${CDSE_USERNAME:-}" ] || [ -z "${CDSE_PASSWORD:-}" ]; }; then
  echo "No CDSE credentials: create cdse_credentials.sh from cdse_credentials.sh.example." >&2
  exit 1
fi
if [ -z "${COPERNICUSMARINE_SERVICE_USERNAME:-}" ] || [ -z "${COPERNICUSMARINE_SERVICE_PASSWORD:-}" ]; then
  echo "No Copernicus Marine credentials (needed for OSTIA): add them to cdse_credentials.sh." >&2
  exit 1
fi
export CDSE_TOKEN_METHOD="${CDSE_TOKEN_METHOD:-curl}"
export CDSE_TOKEN_TIMEOUT="${CDSE_TOKEN_TIMEOUT:-600}"
export CDSE_CONNECT_TIMEOUT="${CDSE_CONNECT_TIMEOUT:-300}"
export COPERNICUSMARINE_CREDENTIALS_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
export COPERNICUSMARINE_CACHE_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
mkdir -p "${COPERNICUSMARINE_CREDENTIALS_DIRECTORY}"

echo "Job started: $(date)"
for name in ${SETS:-jgi_flinkstrom tara_pacific emo_bon}; do
  root="data/external/${name}"
  dir="${root}/manifests"
  if [ ! -f "${dir}/biogoship_sites.csv" ]; then
    echo "Missing ${dir}/biogoship_sites.csv: git pull." >&2
    exit 1
  fi
  echo; echo "===== ${name}: $(date)"
  python scripts/search_products.py --sites-csv "${dir}/biogoship_sites.csv" --out-dir "${dir}" \
    --max-days "${MAX_DAYS:-3}" --workers 1 --pause "${PAUSE:-2}"
  python scripts/download_tiles.py --manifest-dir "${dir}" --data-dir "${root}" --products-dir data/products \
    --positions fixed --min-valid "${MIN_VALID:-0.8}" --max-offset-km "${MAX_OFFSET_KM:-20}" \
    --download-mode files --connections 4
  python scripts/build_tile_manifest.py --manifest-dir "${dir}" --tile-log "${root}/tile_log.csv" \
    --out-csv "${dir}/biogoship_tile_manifest.csv"
  python scripts/fill_sst_l4.py --manifest-csv "${dir}/biogoship_tile_manifest.csv" \
    --tiles-root "${root}/tiles" --out-root "${root}/tiles_sst_l4" \
    --out-manifest "${dir}/biogoship_tile_manifest_sst_l4.csv" --log-csv "${root}/sst_l4_log.csv" \
    --l4-dir "${root}/sst_l4" --split-gap-days 5
  mkdir -p "${root}/embeddings_sst_l4"
  python scripts/pool_tile_features.py --manifest-csv "${dir}/biogoship_tile_manifest_sst_l4.csv" \
    --out-csv "${root}/embeddings_sst_l4/tile_raw_features.csv"
done
df -h "${REPO_ROOT}" | tail -1
echo "Job finished: $(date)"
