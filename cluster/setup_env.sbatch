#!/bin/bash
#SBATCH --job-name=bgs_env
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# One-off: create the project's conda environment on scratch, with OceanParcels and
# the Copernicus Marine toolbox (needed by cluster/backtrack.sbatch).
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/setup_env.sbatch

set -euo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
SCRATCH_USER_DIR="${SCRATCH_USER_DIR:-/scratch/morrill/users/${USER}}"
CONDA_SH="${CONDA_SH:-${SCRATCH_USER_DIR}/miniforge3/etc/profile.d/conda.sh}"
ENV_PATH="${BGS_CONDA_ENV:-${SCRATCH_USER_DIR}/conda_envs/biogoship}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-${SCRATCH_USER_DIR}/conda_envs/pip_cache}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${SCRATCH_USER_DIR}/conda_envs/conda_pkgs}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/.cache/xdg}"

set +u
source "${CONDA_SH}"
set -u
if [ -d "${ENV_PATH}" ]; then
  echo "Environment already exists at ${ENV_PATH}; updating it."
  conda install -y -p "${ENV_PATH}" -c conda-forge --file "${REPO_ROOT}/cluster/conda_packages.txt"
else
  conda create -y -p "${ENV_PATH}" -c conda-forge --file "${REPO_ROOT}/cluster/conda_packages.txt"
fi
"${ENV_PATH}/bin/python" -c "
import numpy, h5py, pandas, xarray, parcels, copernicusmarine
print('numpy', numpy.__version__, '| xarray', xarray.__version__, '| parcels', parcels.__version__, '| copernicusmarine', copernicusmarine.__version__)"
echo "Environment ready: ${ENV_PATH}"
