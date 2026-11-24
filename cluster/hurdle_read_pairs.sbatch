#!/bin/bash
#SBATCH --job-name=bgs_hurdle_rp
#SBATCH --partition=morrill
#SBATCH --array=0-23
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# R70 rerun of both hurdle fits with depth = log10 read pairs sequenced (not KO-assigned reads), 24 shards.
# Per shard, in order (both resumable, fixed-column .partial files):
#   1. scripts/hurdle_kos.py  --parts fit   -> data/results/hurdle/fit_shard<i>.csv(.partial)   tech / Granite / Granite + tech,
#                                                                                                 abundance where present
#   2. scripts/hurdle_base.py --parts fit   -> data/results/hurdle/base_shard<i>.csv(.partial)  predicted C + tech
#
# Inputs (not in git): unpack data/results/hurdle_rp_inputs.tar.gz in the repo root first. It holds
#   data/results/hurdle/{presence_tile.csv.gz, presence_unrarefied_focus.csv.gz, abundance_tile.csv.gz,
#                        tile_depth.csv (now with log10_read_pairs), c_nested.npz (from the read-pair main model)}
#   data/results/benchmark_export/{tiles.csv, granite/tile_embeddings.npy}
#
#   cd /scratch/morrill/users/hmp278/biogoship_granite_kegg && git pull
#   tar -xzf data/results/hurdle_rp_inputs.tar.gz
#   rm -f data/results/hurdle/fit_shard*.csv* data/results/hurdle/base_shard*.csv*     # old (KO-read depth) shards
#   sbatch --account=191001-364393 cluster/hurdle_read_pairs.sbatch
#
# Afterwards bring back data/results/hurdle/fit_shard*.csv* and base_shard*.csv* and run locally:
#   scripts/hurdle_kos.py --parts focus summary ; scripts/hurdle_base.py --parts summary figure

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
if ! python -c 'import sklearn, pandas' 2>/dev/null; then
  echo "scikit-learn / pandas missing in ${BGS_CONDA_ENV}" >&2
  exit 3
fi
python - <<'EOF'
import pandas as pd
d = pd.read_csv("data/results/hurdle/tile_depth.csv", index_col=0)
assert "log10_read_pairs" in d and d["log10_read_pairs"].notna().all(), "tile_depth.csv lacks read pairs: unpack the R70 inputs"
EOF
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
i="${SLURM_ARRAY_TASK_ID}"
echo "Job started: $(date); shard ${i} of 24"
python scripts/hurdle_kos.py --parts fit --shard "${i}" 24
echo "hurdle_kos fit done: $(date)"
python scripts/hurdle_base.py --parts fit --shard "${i}" 24
echo "Job finished: $(date)"
