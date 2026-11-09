#!/bin/bash
#SBATCH --job-name=bgs_model_fields
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# GLORYS12 and PISCES values at each sample's own position, UTC day and depth
# (scripts/extract_model_fields.py). CPU only, needs internet and the Copernicus
# Marine login (cdse_credentials.sh, as for the currents and OSTIA jobs).
# Reanalysis first, then its interim continuation; the analysis/forecast system
# only for dates neither reaches (recorded per value in *_dataset columns).
# Bring back: data/model_fields/ (3 files).
# Resumable: values are saved per model day in data/model_fields/partial/; after a
# time-limit stop, re-submit the same command and it continues where it stopped.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/model_fields.sbatch

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
python scripts/extract_model_fields.py --list-datasets
python scripts/extract_model_fields.py --allow-forecast --workers "${WORKERS:-4}"
echo "Job finished: $(date)"
echo "Bring back: data/model_fields/sample_model_fields.csv, sample_model_fields_long.csv, provenance.json"
