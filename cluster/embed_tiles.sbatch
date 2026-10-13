#!/bin/bash
#SBATCH --job-name=bgs_embed
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Granite embeddings for every clear tile, fixed and advected positions alike.
# No internet here: run cluster/cache_granite.sbatch once beforehand.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/embed_tiles.sbatch
#   sbatch --account=191001-364393 --export=ALL,POSITIONS=fixed cluster/embed_tiles.sbatch

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
export HF_HUB_OFFLINE=1

EXTRA_ARGS=()
if [ -n "${POSITIONS:-}" ]; then
  # shellcheck disable=SC2206  # POSITIONS is a space-separated list
  EXTRA_ARGS+=(--positions ${POSITIONS})
fi

echo "Job started: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
python scripts/build_tile_manifest.py
python scripts/embed_tiles.py --local-files-only ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}
echo "Job finished: $(date)"
