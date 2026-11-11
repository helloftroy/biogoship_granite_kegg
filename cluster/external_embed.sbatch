#!/bin/bash
#SBATCH --job-name=bgs_ext_embed
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Granite embeddings of the OSTIA-filled tiles of the external test sets
# (after cluster/external_tiles.sbatch):
#
#   sbatch --account=191001-364393 cluster/external_embed.sbatch
#   sbatch --account=191001-364393 --dependency=afterok:<external_tiles job id> cluster/external_embed.sbatch

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
export HF_HUB_OFFLINE=1

echo "Job started: $(date)"
for name in ${SETS:-jgi_flinkstrom tara_pacific emo_bon}; do
  root="data/external/${name}"
  if [ ! -s "${root}/manifests/biogoship_tile_manifest_sst_l4.csv" ]; then
    echo "Skipping ${name}: no OSTIA tile manifest yet." >&2
    continue
  fi
  python scripts/embed_tiles.py --local-files-only --positions fixed \
    --manifest-csv "${root}/manifests/biogoship_tile_manifest_sst_l4.csv" --out-dir "${root}/embeddings_sst_l4"
done
echo "Job finished: $(date)"
echo "Bring back, per set: data/external/<set>/embeddings_sst_l4/ (bundle zip and tile_raw_features.csv),"
echo "  data/external/<set>/manifests/biogoship_tile_manifest_sst_l4.csv, site_tile_status.csv, data/external/<set>/sst_l4_log.csv"
