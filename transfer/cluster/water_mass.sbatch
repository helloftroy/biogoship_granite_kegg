#!/bin/bash
#SBATCH --job-name=bgs_water_mass
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# 60-day back-trajectories of every sample and the ocean along them (transfer/scripts/water_mass_history.py,
# RESULTS.md R42): currents download, OceanParcels back-tracking, then SST / chlorophyll / mixed layer /
# KD490 / sea-floor depth at the daily positions. Also probes which Copernicus ocean-colour datasets
# (PAR, KD490) open, for the light-history test. Needs internet + the Copernicus Marine login, and the
# biogoship conda env (parcels), as cluster/backtrack.sbatch. Resumable: re-submit after a stop.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg && git pull
#   sbatch --account=191001-364393 transfer/cluster/water_mass.sbatch
#
# Bring back (scp to the same paths on the Mac):
#   transfer/results/water_mass/trajectories.csv and along_track.csv, and this job's .out

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
if [ -f cdse_credentials.sh ]; then
  source cdse_credentials.sh
fi
if [ -z "${COPERNICUSMARINE_SERVICE_USERNAME:-}" ] || [ -z "${COPERNICUSMARINE_SERVICE_PASSWORD:-}" ]; then
  echo "No Copernicus Marine credentials: add them to cdse_credentials.sh." >&2
  exit 1
fi
export COPERNICUSMARINE_CREDENTIALS_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
export COPERNICUSMARINE_CACHE_DIRECTORY="${REPO_ROOT}/.cache/copernicusmarine"
mkdir -p "${COPERNICUSMARINE_CREDENTIALS_DIRECTORY}"

echo "Job started: $(date)"
python transfer/scripts/water_mass_history.py --probe-light || true
python transfer/scripts/water_mass_history.py --steps currents track fields --workers "${WORKERS:-4}"
echo "Job finished: $(date)"
echo "Bring back: transfer/results/water_mass/trajectories.csv, along_track.csv and this .out"
