#!/bin/bash
#SBATCH --job-name=bgs_extend
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# For sites still without a clear advected tile, look for Sentinel-3 scenes further
# out in time along each site's water track (beyond the +/-3 days already searched).
# Run after backtrack.sbatch and after download_tiles.sbatch has tried the advected
# position on the scenes already downloaded. Needs internet (catalogue queries).
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/extend_search.sbatch                         # out to +/-6 days
#   sbatch --account=191001-364393 --export=ALL,MAX_DAYS=7 cluster/extend_search.sbatch
#
# Then: download_products.sbatch (fetches only the new scenes), download_tiles.sbatch.
# Can be repeated with a larger MAX_DAYS; each round only searches sites still missing.
# MAX_DAYS cannot exceed the track length (TRACK_DAYS in backtrack.sbatch, default 7).

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"

echo "Job started: $(date)"
python scripts/search_trajectory_products.py --max-days "${MAX_DAYS:-6}"
echo
echo "What download_products.sbatch would now fetch:"
python scripts/download_products.py --dry-run
df -h "${REPO_ROOT}" | tail -1
echo "Job finished: $(date)"
