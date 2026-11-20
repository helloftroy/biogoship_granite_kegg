#!/bin/bash
#SBATCH --job-name=bgs_merge_cov0
#SBATCH --partition=morrill
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Merge the per-library coverage-0 KO counts (one folder per library, level_*_ko_counts.tsv) into one wide
# table in the layout of ko_counts_cov50_all.tsv (scripts/merge_ko_counts.py; RESULTS.md R47).
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg && git pull
#   sbatch --account=191001-364393 --export=ALL,ROOT=/path/to/cov0/folders cluster/merge_ko_counts.sbatch
#
# ROOT: the directory holding the per-library folders (searched recursively; several paths may be given,
# space-separated). Bring back: ko_counts_cov0_all.tsv.gz and ko_counts_cov0_all_domain_counts.tsv
# (put them in the repo root on the Mac) and this job's .out. Delete the per-library folders only after
# the .out shows the expected number of libraries (about 3,000) and the files are safely on the Mac.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
if [ -z "${ROOT:-}" ]; then
  echo "Set ROOT to the folder holding the per-library cov0 folders." >&2
  exit 1
fi
echo "Job started: $(date)"
du -sh ${ROOT} || true
echo "library folders: $(find ${ROOT} -name 'level_1_*_ko_counts.tsv*' | wc -l)"
python scripts/merge_ko_counts.py ${ROOT} --out ko_counts_cov0_all.tsv
gzip -f ko_counts_cov0_all.tsv
ls -la ko_counts_cov0_all.tsv.gz ko_counts_cov0_all_domain_counts.tsv
echo "Job finished: $(date)"
echo "Bring back: ko_counts_cov0_all.tsv.gz, ko_counts_cov0_all_domain_counts.tsv and this .out"
