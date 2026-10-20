#!/bin/bash
#SBATCH --job-name=bgs_granite_check
#SBATCH --partition=morrill
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Due diligence on the Granite inputs: compares our tile format with the model's own
# config, and tests real tiles with Granite's masked-autoencoder reconstruction under
# different gap-filling conventions. CPU only, offline (needs cache_granite.sbatch first).
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/check_granite_inputs.sbatch

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
export HF_HUB_OFFLINE=1

echo "Job started: $(date)"
python scripts/check_granite_inputs.py --local-files-only --n-tiles "${N_TILES:-300}"
echo "Job finished: $(date)"
