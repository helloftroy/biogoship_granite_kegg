#!/bin/bash
#SBATCH --job-name=bgs_hurdle
#SBATCH --partition=morrill
#SBATCH --array=0-23
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Hurdle presence / abundance fit for the 5,808 patchy KOs (RESULTS.md R65), 24 shards as an array job.
# Resumable: each shard appends to data/results/hurdle/fit_shard<i>.csv.partial; resubmitting skips finished KOs.
#
# Inputs (not in git; unpack data/results/hurdle_cluster_inputs.tar.gz in the repo root first):
#   data/results/hurdle/{presence_tile.csv.gz, presence_unrarefied_focus.csv.gz, abundance_tile.csv.gz, tile_depth.csv}
#   data/results/benchmark_export/{tiles.csv, granite/tile_embeddings.npy}
# Needs scikit-learn in the environment (see the check below).
#
#   sbatch --account=191001-364393 cluster/hurdle_fit.sbatch
#
# Afterwards copy data/results/hurdle/fit_shard*.csv back and run locally:
#   scripts/hurdle_kos.py --parts summary figure

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
if ! python -c 'import sklearn, pandas' 2>/dev/null; then
  echo "scikit-learn / pandas missing in ${BGS_CONDA_ENV}. Install once from a node with internet, e.g.:" >&2
  echo "  srun --account=191001-364393 --partition=service --time=00:30:00 bash -c 'source cluster/env.sh && python -m pip install scikit-learn pandas'" >&2
  exit 3
fi
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
echo "Job started: $(date); shard ${SLURM_ARRAY_TASK_ID} of 24"
python scripts/hurdle_kos.py --parts fit --shard "${SLURM_ARRAY_TASK_ID}" 24
echo "Job finished: $(date)"
