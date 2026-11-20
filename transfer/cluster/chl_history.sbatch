#!/bin/bash
#SBATCH --job-name=bgs_chl_history
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Satellite chlorophyll over the 28 days before each sample (transfer/scripts/extract_chl_history.py,
# RESULTS.md R35): gap-free daily GlobColour L4 (4 km), box mean +-0.125 deg; growth / decay slopes
# over 14 and 28 days, distance from the 28-day peak. CPU only, needs internet and the Copernicus
# Marine login (cdse_credentials.sh), like cluster/history_fields.sbatch.
# The first step lists which dataset ids open; ids that do not are skipped and reported.
# Resumable: values are appended per sample to transfer/results/chl_history/partial_chl.csv;
# after a stop, re-submit the same command and it continues.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   git pull
#   sbatch --account=191001-364393 transfer/cluster/chl_history.sbatch
#
# Bring back (scp to the Mac, same path): transfer/results/chl_history/sample_chl_history.csv
# and provenance_chl_history.json; paste back this job's .out file.

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
python transfer/scripts/extract_chl_history.py --list-datasets
python transfer/scripts/extract_chl_history.py --workers "${WORKERS:-4}"
echo "Job finished: $(date)"
echo "Bring back: transfer/results/chl_history/sample_chl_history.csv and provenance_chl_history.json"
