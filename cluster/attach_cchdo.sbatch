#!/bin/bash
#SBATCH --job-name=bgs_cchdo
#SBATCH --partition=service
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# CCHDO bottle chemistry (salinity, oxygen, nutrients, DIC, alkalinity, pH) for the
# sequenced samples of AMT28, I07N, I09N 2016, P18 and C13.5/A13.5: CTD samples by
# station number, underway samples from the nearest cast (<= 30 km, <= 12 h).
# Writes data/manifests/biogoship_samples_cchdo.csv. Needs internet; a few minutes.
#
#   sbatch --account=191001-364393 cluster/attach_cchdo.sbatch
#
# Overrides: MAX_KM=30, MAX_HOURS=12.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
echo "Job started: $(date)"
python scripts/attach_cchdo_bottles.py --max-km "${MAX_KM:-30}" --max-hours "${MAX_HOURS:-12}"
echo "Job finished: $(date)"
