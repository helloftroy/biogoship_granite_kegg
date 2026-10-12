#!/bin/bash
#SBATCH --job-name=bgs_tiles
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Download Sentinel-3 OLCI + SST and cut Granite tiles for every Bio-GO-SHIP site.
# Needs internet, hence the 'service' partition. Safe to resubmit: it resumes.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   cp cdse_credentials.sh.example cdse_credentials.sh   # once; then fill it in
#
#   # 1. smoke test on three scenes
#   sbatch --account=<allocation> --export=ALL,MAX_PRODUCTS=3 cluster/download_tiles.sbatch
#   # 2. everything
#   sbatch --account=<allocation> cluster/download_tiles.sbatch
#
# Other overrides: MIN_VALID=0.8, DOWNLOAD_MODE=files|zip, NO_SST=1, CONNECTIONS=4.

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
if [ "${NO_SST:-0}" = "1" ]; then
  EXTRA_ARGS+=(--no-sst)
fi

echo "Job started: $(date)"
df -h "${REPO_ROOT}" | tail -1
python scripts/download_tiles.py \
  --min-valid "${MIN_VALID:-0.8}" \
  --max-products "${MAX_PRODUCTS:-0}" \
  --download-mode "${DOWNLOAD_MODE:-files}" \
  --connections "${CONNECTIONS:-4}" \
  ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}

python scripts/build_tile_manifest.py
echo "Job finished: $(date)"
