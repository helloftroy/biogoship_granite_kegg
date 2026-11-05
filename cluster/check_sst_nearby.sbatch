#!/bin/bash
#SBATCH --job-name=bgs_sst_nearby
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Check only, nothing changed: for samples whose best tile has no SST, how far would the
# tile have to move (within 40 km, same scene) to be clear in OLCI and have SST?
# Reads the OLCI and SST archives on disk; fetches up to FILL_CANDIDATES daytime SST
# overpasses per scene (single L2P files, needs cdse_credentials.sh). Writes
# data/sst_nearby_check.csv and prints recoveries per cruise at 10, 20 and 40 km.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/check_sst_nearby.sbatch
#   sbatch --account=191001-364393 --export=ALL,FILL_CANDIDATES=0 cluster/check_sst_nearby.sbatch   # paired SST only

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
if [ -f cdse_credentials.sh ]; then
  source cdse_credentials.sh
fi
export CDSE_TOKEN_METHOD="${CDSE_TOKEN_METHOD:-curl}"
export CDSE_TOKEN_TIMEOUT="${CDSE_TOKEN_TIMEOUT:-600}"
export CDSE_CONNECT_TIMEOUT="${CDSE_CONNECT_TIMEOUT:-300}"

echo "Job started: $(date)"
python scripts/check_sst_nearby.py --fill-candidates "${FILL_CANDIDATES:-2}"
echo "Job finished: $(date)"
