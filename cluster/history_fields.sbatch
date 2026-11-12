#!/bin/bash
#SBATCH --job-name=bgs_history_fields
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Recent ocean history before each sample (scripts/extract_history_fields.py, RESULTS.md R27):
# 30-day max mixed layer (GLORYS), 14-day SST change (OSTIA), 14-day peak and mean wind
# (L4 hourly winds), sea-level gradient (GLORYS), 30-day eddy kinetic energy (DUACS).
# CPU only, needs internet and the Copernicus Marine login (cdse_credentials.sh).
# The first step lists which dataset ids open; ids that do not are skipped and reported.
# Resumable: values are saved per sample in data/model_fields/partial_history/; after a
# stop, re-submit the same command and it continues.
# Bring back: data/model_fields/sample_history_fields.csv and provenance_history.json.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/history_fields.sbatch

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"

if [ -f cdse_credentials.sh ]; then
  source cdse_credentials.sh
fi
if [ -z "${COPERNICUSMARINE_SERVICE_USERNAME:-}" ] || [ -z "${COPERNICUSMARINE_SERVICE_PASSWORD:-}" ]; then
  echo "No Copernicus Marine credentials: add them to cdse_credentials.sh (see cdse_credentials.sh.example)." >&2
  exit 1
fi
export COPERNICUSMARINE_CREDENTIALS_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
export COPERNICUSMARINE_CACHE_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
mkdir -p "${COPERNICUSMARINE_CREDENTIALS_DIRECTORY}"

echo "Job started: $(date)"
python scripts/extract_history_fields.py --list-datasets
python scripts/extract_history_fields.py --workers "${WORKERS:-4}"
echo "Job finished: $(date)"
echo "Bring back: data/model_fields/sample_history_fields.csv and provenance_history.json"
