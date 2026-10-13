#!/usr/bin/env bash
# Shared setup for every job. Source it from the repo root:  source cluster/env.sh
#
# Override any of these at submit time, e.g.
#   sbatch --export=ALL,BGS_CONDA_ENV=/path/to/env cluster/download_tiles.sbatch

# SLURM runs a spooled copy of the job script, so $0 is useless for finding the repo.
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [ ! -d "${REPO_ROOT}/biogoship" ] || [ ! -d "${REPO_ROOT}/scripts" ]; then
  echo "Run this from the biogoship_granite_kegg repo root (cd there before sbatch)." >&2
  exit 1
fi
cd "${REPO_ROOT}"

# Everything sizeable stays next to the repo on scratch; $HOME is at quota.
SCRATCH_ROOT="${SCRATCH_ROOT:-${REPO_ROOT}}"
SCRATCH_USER_DIR="${SCRATCH_USER_DIR:-/scratch/morrill/users/${USER}}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${SCRATCH_ROOT}/.cache/xdg}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${SCRATCH_ROOT}/.cache/mpl}"
export HF_HOME="${HF_HOME:-${SCRATCH_ROOT}/.cache/huggingface}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-${SCRATCH_USER_DIR}/conda_envs/pip_cache}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${SCRATCH_USER_DIR}/conda_envs/conda_pkgs}"
export PYTHONUNBUFFERED=1
export TMPDIR="${TMPDIR:-${SCRATCH_ROOT}/.cache/tmp}"
export PIP_DISABLE_PIP_VERSION_CHECK=1
mkdir -p "${XDG_CACHE_HOME}" "${MPLCONFIGDIR}" "${HF_HOME}" "${TMPDIR}" "${REPO_ROOT}/logs"

# The environment is activated by full path, never by name. It is either a conda
# environment (conda is not on PATH by default, so conda.sh is sourced first) or a
# plain virtual environment such as the one cluster/setup_env.sbatch builds.
CONDA_SH="${CONDA_SH:-${SCRATCH_USER_DIR}/miniforge3/etc/profile.d/conda.sh}"
BGS_CONDA_ENV="${BGS_CONDA_ENV:-${SCRATCH_USER_DIR}/conda_envs/faire-agent}"
if [ ! -d "${BGS_CONDA_ENV}" ]; then
  echo "Python environment not found at ${BGS_CONDA_ENV}; set BGS_CONDA_ENV, or run cluster/setup_env.sbatch." >&2
  exit 2
fi
export PYTHONNOUSERSITE=1  # ignore ~/.local
set +u  # activation scripts read unset variables
if [ -f "${BGS_CONDA_ENV}/pyvenv.cfg" ]; then
  source "${BGS_CONDA_ENV}/bin/activate"
else
  if [ ! -f "${CONDA_SH}" ]; then
    echo "conda.sh not found at ${CONDA_SH}; set CONDA_SH." >&2
    exit 2
  fi
  source "${CONDA_SH}"
  conda activate "${BGS_CONDA_ENV}"
fi
set -u
_actual_prefix="$(python -c 'import sys; print(sys.prefix)')"
if [ "$(cd "${_actual_prefix}" && pwd -P)" != "$(cd "${BGS_CONDA_ENV}" && pwd -P)" ]; then
  echo "Environment activation did not take effect: python is ${_actual_prefix}, expected ${BGS_CONDA_ENV}." >&2
  exit 2
fi

# The pipeline needs only numpy and h5py; install them if this environment lacks them.
# (pip needs internet, which the 'service' partition has.)
if ! python -c 'import numpy, h5py' 2>/dev/null; then
  echo "Installing missing Python requirements into ${BGS_CONDA_ENV}..."
  python -m pip install --quiet -r "${REPO_ROOT}/requirements.txt"
fi

echo "Host: $(hostname) | repo: ${REPO_ROOT} | python: $(command -v python) ($(python --version 2>&1))"
