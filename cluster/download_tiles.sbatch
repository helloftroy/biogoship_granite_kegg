#!/bin/bash
#SBATCH --job-name=bgs_tiles
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Cut Granite tiles for every Bio-GO-SHIP site, at the sampling position and at the
# back-tracked position. Reads the archives from download_products.sbatch if present,
# otherwise downloads what it needs. Safe to resubmit: it resumes.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   cp cdse_credentials.sh.example cdse_credentials.sh   # once; then fill it in
#
#   # 1. smoke test on three scenes
#   sbatch --account=191001-364393 --export=ALL,MAX_PRODUCTS=3 cluster/download_tiles.sbatch
#   # 2. everything
#   sbatch --account=191001-364393 cluster/download_tiles.sbatch
#
# A tile must be 80% clear; if the one centred on the position is not, the nearest
# clear one within 20 km is used (MAX_OFFSET_KM). Combinations tested earlier under a
# different rule are only re-tested if asked:
#   sbatch --account=191001-364393 --export=ALL,POSITIONS=fixed,RECHECK="cloudy edge" cluster/download_tiles.sbatch
#
# Other overrides: POSITIONS="fixed advected" (or just "fixed"), MIN_VALID=0.8,
# DOWNLOAD_MODE=files|zip, NO_SST=1, CONNECTIONS=4.
#
# Another sample set, with its own tiles, logs and tile manifest (MANIFEST_DIR and DATA_DIR):
#   sbatch --account=191001-364393 --export=ALL,POSITIONS=fixed,MANIFEST_DIR=data/faire/manifests,DATA_DIR=data/faire cluster/download_tiles.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"

if [ -f cdse_credentials.sh ]; then
  source cdse_credentials.sh
fi
if [ -z "${CDSE_ACCESS_TOKEN:-}" ] && { [ -z "${CDSE_USERNAME:-}" ] || [ -z "${CDSE_PASSWORD:-}" ]; }; then
  echo "No Copernicus credentials: create cdse_credentials.sh from cdse_credentials.sh.example." >&2
  exit 1
fi
# curl with generous timeouts was the reliable way to reach the CDSE login server from this cluster.
export CDSE_TOKEN_METHOD="${CDSE_TOKEN_METHOD:-curl}"
export CDSE_TOKEN_TIMEOUT="${CDSE_TOKEN_TIMEOUT:-600}"
export CDSE_CONNECT_TIMEOUT="${CDSE_CONNECT_TIMEOUT:-300}"

EXTRA_ARGS=()
if [ -n "${RECHECK:-}" ]; then
  # shellcheck disable=SC2206  # RECHECK is a space-separated list
  EXTRA_ARGS+=(--recheck ${RECHECK})
fi
if [ "${NO_SST:-0}" = "1" ]; then
  EXTRA_ARGS+=(--no-sst)
fi

echo "Job started: $(date)"
df -h "${REPO_ROOT}" | tail -1
# shellcheck disable=SC2086  # POSITIONS is a space-separated list
MANIFEST_DIR="${MANIFEST_DIR:-data/manifests}"
DATA_DIR="${DATA_DIR:-data}"
python scripts/download_tiles.py \
  --manifest-dir "${MANIFEST_DIR}" \
  --data-dir "${DATA_DIR}" \
  --products-dir data/products \
  --positions ${POSITIONS:-fixed advected} \
  --min-valid "${MIN_VALID:-0.8}" \
  --max-offset-km "${MAX_OFFSET_KM:-20}" \
  --max-products "${MAX_PRODUCTS:-0}" \
  --download-mode "${DOWNLOAD_MODE:-files}" \
  --connections "${CONNECTIONS:-4}" \
  ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}

python scripts/build_tile_manifest.py --manifest-dir "${MANIFEST_DIR}" --tile-log "${DATA_DIR}/tile_log.csv" \
  --out-csv "${MANIFEST_DIR}/biogoship_tile_manifest.csv"
echo "Job finished: $(date)"
