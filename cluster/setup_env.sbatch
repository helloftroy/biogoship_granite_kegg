#!/bin/bash
#SBATCH --job-name=bgs_env
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# One-off: create the project's Python environment on scratch, with OceanParcels and
# the Copernicus Marine toolbox (needed by cluster/backtrack.sbatch).
#
# It is a plain virtual environment built with pip, not a conda environment: conda
# records every environment it creates under ~/.conda, and $HOME is at quota.
# Nothing here writes to $HOME: the environment, pip's cache and temporary build
# files all go to scratch.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/setup_env.sbatch

set -euo pipefail
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
SCRATCH_USER_DIR="${SCRATCH_USER_DIR:-/scratch/morrill/users/${USER}}"
ENV_PATH="${BGS_CONDA_ENV:-${SCRATCH_USER_DIR}/conda_envs/biogoship}"
# Python 3.11 from the existing faire-agent environment is the base interpreter.
BASE_PYTHON="${BGS_BASE_PYTHON:-${SCRATCH_USER_DIR}/conda_envs/faire-agent/bin/python}"

export PIP_CACHE_DIR="${PIP_CACHE_DIR:-${SCRATCH_USER_DIR}/conda_envs/pip_cache}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${REPO_ROOT}/.cache/xdg}"
export TMPDIR="${REPO_ROOT}/.cache/tmp"
export PYTHONNOUSERSITE=1
export PIP_DISABLE_PIP_VERSION_CHECK=1
mkdir -p "${PIP_CACHE_DIR}" "${XDG_CACHE_HOME}" "${TMPDIR}" "$(dirname "${ENV_PATH}")"

case "${ENV_PATH}" in
  /scratch/*) ;;
  *) echo "Refusing to build an environment outside /scratch: ${ENV_PATH}" >&2; exit 2 ;;
esac
if [ ! -x "${BASE_PYTHON}" ]; then
  echo "Base Python not found at ${BASE_PYTHON}; set BGS_BASE_PYTHON to a Python 3.10-3.12 on scratch." >&2
  exit 2
fi

echo "Job started: $(date)"
echo "Environment: ${ENV_PATH}"
echo "Base Python: ${BASE_PYTHON} ($("${BASE_PYTHON}" --version 2>&1))"
df -h "${SCRATCH_USER_DIR}" | tail -1

if [ -d "${ENV_PATH}" ] && [ ! -f "${ENV_PATH}/pyvenv.cfg" ]; then
  echo "Removing the incomplete environment left by an earlier attempt at ${ENV_PATH}"
  rm -rf "${ENV_PATH}"
fi
if [ ! -f "${ENV_PATH}/pyvenv.cfg" ]; then
  "${BASE_PYTHON}" -m venv "${ENV_PATH}"
fi
"${ENV_PATH}/bin/python" -m pip install --upgrade pip
"${ENV_PATH}/bin/python" -m pip install -r "${REPO_ROOT}/cluster/env_packages.txt"

# Check the pieces work together, not just that they import: the Copernicus Marine
# toolbox reads its data through xarray + zarr, which breaks if their versions drift apart.
"${ENV_PATH}/bin/python" - <<'PY'
import sys, tempfile
import numpy, h5py, pandas, xarray, zarr, parcels, copernicusmarine

print("python", sys.version.split()[0], "| prefix", sys.prefix)
print("numpy", numpy.__version__, "| xarray", xarray.__version__, "| zarr", zarr.__version__,
      "| parcels", parcels.__version__, "| copernicusmarine", copernicusmarine.__version__)
with tempfile.TemporaryDirectory() as tmp:
    data = xarray.Dataset({"uo": (("time", "x"), numpy.ones((2, 3), "float32"))})
    data.to_zarr(f"{tmp}/check.zarr", consolidated=True)
    assert float(xarray.open_zarr(f"{tmp}/check.zarr", consolidated=True)["uo"].sum()) == 6.0
print("xarray + zarr read/write check passed")
PY
echo "Environment size: $(du -sh "${ENV_PATH}" | cut -f1)"
echo "Environment ready: ${ENV_PATH}"
echo "Job finished: $(date)"
