#!/bin/bash
#SBATCH --job-name=bgs_combo
#SBATCH --partition=morrill
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# R73: nested search over input-group combinations (scripts/combo_search.py). Two arrays, one per target:
#
#   cd /scratch/morrill/users/hmp278/biogoship_granite_kegg && git pull
#   mkdir -p data/results/combo      # then scp data/results/combo/inputs.npz from the Mac into it
#   sbatch --account=191001-364393 --array=0-63 cluster/combo_search.sbatch C 64     # 32,767 subsets, ~15 min per task
#   sbatch --account=191001-364393 --array=0-31 cluster/combo_search.sbatch ko 32    # 127 subsets x 3,933 KOs
#
# Each task writes data/results/combo/<target>_shard<i>.npz when it finishes. Copy them back and run locally:
#   scripts/combo_search.py --parts summary --target C   (and --target ko)
# Input: data/results/combo/inputs.npz (built locally with --parts bundle; not in git).

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
TARGET="${1:?target C or ko}"
NSHARDS="${2:?number of shards (must equal the array size)}"
if ! python -c 'import sklearn, numpy, pandas' 2>/dev/null; then
  echo "scikit-learn / numpy / pandas missing in ${BGS_CONDA_ENV}" >&2
  exit 3
fi
test -s data/results/combo/inputs.npz || { echo "data/results/combo/inputs.npz missing: scp it from the Mac" >&2; exit 4; }
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}" OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
echo "Job started: $(date); ${TARGET} shard ${SLURM_ARRAY_TASK_ID} of ${NSHARDS}"
python scripts/combo_search.py --parts fit --target "${TARGET}" --shard "${SLURM_ARRAY_TASK_ID}" "${NSHARDS}"
echo "Job finished: $(date)"
