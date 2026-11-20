#!/bin/bash
#SBATCH --job-name=bgs_ko_gc
#SBATCH --partition=morrill
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# GC content of the reference genes per KO (transfer/scripts/ko_gc_content.py, RESULTS.md R40),
# for checking whether the axis-1 platform offset is a sequencing GC bias.
# Needs data/reference/biogoship_2026/biogoship_2026_database_v2.fasta (212 GB), which the
# R32 download skipped. Get it first (service partition, internet; same rclone setup as R32):
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg && git pull
#   df -h /scratch/morrill/users/$USER        # needs about 215 GB free
#   jid=$(sbatch --parsable --account=191001-364393 \
#          --export=ALL,FILES="biogoship_2026_database_v2.fasta" cluster/download_reference.sbatch)
#   sbatch --account=191001-364393 --dependency=afterok:${jid} transfer/cluster/ko_gc.sbatch
#
# Reads the eggNOG table (105 GB), the taxonomy (29 GB) and the FASTA (212 GB) once each,
# streaming; keeps 1 gene in SAMPLE (default 10). Expect several hours.
# Bring back: transfer/results/ko_gc/ko_gc.csv (a few hundred kB) and this job's .out.
# The FASTA can be deleted afterwards if space is needed (it is not used for mapping).

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
echo "Job started: $(date)"
ls -la data/reference/biogoship_2026/
python transfer/scripts/ko_gc_content.py --sample "${SAMPLE:-10}"
echo "Job finished: $(date)"
echo "Bring back: transfer/results/ko_gc/ko_gc.csv"
