#!/bin/bash
#SBATCH --job-name=bgs_new_cruises
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Newer Bio-GO-SHIP cruises (P02 2022, A16N 2023, I05 2023, A13.5 2024, I08S 2024,
# I09N 2025, A16S 2026): candidate surface sampling events from CCHDO bottle files
# and SeaBASS biological sampling records, then a Sentinel-3 scene search for them.
# A16S is catalogued but not matched to satellite data. Nothing is downloaded from
# Sentinel-3 here; the job ends with what download_products.py would fetch.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/new_cruises.sbatch
#
# SeaBASS needs EARTHDATA_USERNAME / EARTHDATA_PASSWORD in cdse_credentials.sh;
# without them SeaBASS is listed but not downloaded. Overrides: MAX_DEPTH=10, MAX_DAYS=3,
# PAUSE=2 (seconds between catalogue requests). The scene search is cached: if the
# catalogue's rate limit stops it, re-run; SEARCH_ONLY=1 skips rebuilding the sites:
#   sbatch --account=191001-364393 --export=ALL,SEARCH_ONLY=1 cluster/new_cruises.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
if [ -f cdse_credentials.sh ]; then
  source cdse_credentials.sh
fi
NEW="data/new_cruises"

echo "Job started: $(date)"
if [ "${SEARCH_ONLY:-0}" != "1" ]; then
  python scripts/build_new_cruise_manifest.py --out-dir "${NEW}" --max-depth "${MAX_DEPTH:-10}"
fi
python scripts/search_products.py --sites-csv "${NEW}/manifests/biogoship_sites.csv" --out-dir "${NEW}/manifests" \
  --max-days "${MAX_DAYS:-3}" --workers 1 --pause "${PAUSE:-2}"
echo
echo "Sentinel-3 products the new cruises would need:"
python scripts/download_products.py --manifest-dir "${NEW}/manifests" --dry-run
df -h "${REPO_ROOT}" | tail -1
echo "Job finished: $(date)"
