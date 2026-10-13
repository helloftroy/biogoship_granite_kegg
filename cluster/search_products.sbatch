#!/bin/bash
#SBATCH --job-name=bgs_search
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Optional: the search results are already committed in data/manifests/. Re-run this
# only after changing the sample list or the time window.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/search_products.sbatch
#   sbatch --account=191001-364393 --export=ALL,MAX_DAYS=5,I09N_MAX_DAYS=10 cluster/search_products.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"

echo "Job started: $(date)"
python scripts/build_sample_manifest.py
python scripts/search_products.py --max-days "${MAX_DAYS:-3}" --section-max-days "I09N=${I09N_MAX_DAYS:-10}"
echo "Job finished: $(date)"
