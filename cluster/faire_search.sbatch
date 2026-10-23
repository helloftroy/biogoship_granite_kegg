#!/bin/bash
#SBATCH --job-name=bgs_faire_search
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Sentinel-3 scene search for the sequenced samples of the newer cruises, from the
# FAIRe-NOAA sample metadata (A13.5 2024, A16N 2023, GOMECC4, HLY2202, I05 2023,
# I08S 2024, I09N 2025, P02 2022, P04E 2026). Ends with what download_products.py
# would fetch; nothing is downloaded here.
#
# data/faire/ is built off the cluster (the workbook needs openpyxl) and copied in:
#   python scripts/build_faire_manifest.py FAIRe-NOAA_uci-ess-biogoship.xlsx      # locally
#   scp -r data/faire hpc2:/scratch/morrill/users/hmp278/biogoship_granite_kegg/data/
#
#   sbatch --account=191001-364393 cluster/faire_search.sbatch
#
# The search is cached and resumable: if the catalogue's rate limit stops it, re-run.
# Overrides: MAX_DAYS=3, PAUSE=2.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
DIR="data/faire/manifests"
if [ ! -f "${DIR}/biogoship_sites.csv" ]; then
  echo "Missing ${DIR}/biogoship_sites.csv: build it locally and copy data/faire/ here (see the header)." >&2
  exit 1
fi

echo "Job started: $(date)"
python scripts/search_products.py --sites-csv "${DIR}/biogoship_sites.csv" --out-dir "${DIR}" \
  --max-days "${MAX_DAYS:-3}" --workers 1 --pause "${PAUSE:-2}"
echo
echo "Sentinel-3 products these sites need:"
python scripts/download_products.py --manifest-dir "${DIR}" --dry-run
df -h "${REPO_ROOT}" | tail -1
echo "Job finished: $(date)"
