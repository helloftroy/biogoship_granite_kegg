#!/bin/bash
#SBATCH --job-name=bgs_products
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Download every candidate Sentinel-3 OLCI scene and its SST orbit file (about 300 GB)
# into data/products/. Needs internet, hence the 'service' partition. Resumes on resubmit.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   cp cdse_credentials.sh.example cdse_credentials.sh   # once; then fill it in
#
#   sbatch --account=<allocation> --export=ALL,MAX_PRODUCTS=2 cluster/download_products.sbatch   # smoke test
#   sbatch --account=<allocation> cluster/download_products.sbatch                               # everything
#
# Other overrides: NO_SST=1, CONNECTIONS=4.

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
echo "Free space:"; df -h "${REPO_ROOT}" | tail -1
python scripts/download_products.py \
  --max-products "${MAX_PRODUCTS:-0}" \
  --connections "${CONNECTIONS:-4}" \
  ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}
echo "Disk used by products: $(du -sh data/products | cut -f1)"
echo "Job finished: $(date)"
