#!/bin/bash
#SBATCH --job-name=bgs_sst_fill
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Fill tiles that have no SST from the closest-in-time SLSTR overpass of the same spot
# (either satellite; daytime first). Needs internet and Copernicus Data Space login.
# Then refreshes the manifest and the raw tile features; re-embed afterwards on a GPU.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 --export=ALL,DRY_RUN=1 cluster/fill_missing_sst.sbatch   # counts and sizes only
#   sbatch --account=191001-364393 cluster/fill_missing_sst.sbatch
#   sbatch --account=191001-364393 cluster/embed_tiles.sbatch
#
# Overrides: MAX_HOURS=48 (largest gap to the OLCI scene), MAX_TRIES=4, FIRST_TRY=0.
# Second pass for tiles still without SST, trying the next 12 overpasses:
#   sbatch --account=191001-364393 --export=ALL,FIRST_TRY=4,MAX_TRIES=12 cluster/fill_missing_sst.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
if [ -f cdse_credentials.sh ]; then
  source cdse_credentials.sh
fi
export CDSE_TOKEN_METHOD="${CDSE_TOKEN_METHOD:-curl}"
export CDSE_TOKEN_TIMEOUT="${CDSE_TOKEN_TIMEOUT:-600}"
export CDSE_CONNECT_TIMEOUT="${CDSE_CONNECT_TIMEOUT:-300}"

EXTRA_ARGS=()
if [ "${DRY_RUN:-0}" = "1" ]; then
  EXTRA_ARGS+=(--dry-run)
fi

echo "Job started: $(date)"
python scripts/fill_missing_sst.py --max-hours "${MAX_HOURS:-48}" --max-tries "${MAX_TRIES:-4}" --first-try "${FIRST_TRY:-0}" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}
if [ "${DRY_RUN:-0}" != "1" ]; then
  python scripts/build_tile_manifest.py
  python scripts/pool_tile_features.py
fi
echo "Job finished: $(date)"
