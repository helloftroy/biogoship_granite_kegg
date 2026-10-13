#!/bin/bash
#SBATCH --job-name=bgs_backtrack
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Download surface currents and back-track every site to the time of each overpass.
# Needs internet for the currents, hence the 'service' partition; the tracking
# itself takes a few minutes. Run cluster/setup_env.sbatch once first.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=<allocation> cluster/backtrack.sbatch
#
# Overrides: DEPTH_M=15 (or 0), MARGIN_DEG=5.

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
# Keep the toolbox's own files off $HOME.
export COPERNICUSMARINE_CREDENTIALS_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
export COPERNICUSMARINE_CACHE_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
mkdir -p "${COPERNICUSMARINE_CREDENTIALS_DIRECTORY}"

echo "Job started: $(date)"
python scripts/download_currents.py --margin-deg "${MARGIN_DEG:-5}"
python scripts/backtrack_sites.py --depth-m "${DEPTH_M:-15}"
echo "Job finished: $(date)"
