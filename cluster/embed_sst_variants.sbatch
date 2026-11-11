#!/bin/bash
#SBATCH --job-name=bgs_embed_sst
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Granite embeddings of the two SST variants, for both sample sets (after cluster/sst_variants.sbatch):
#   blank: every tile's SST channel set to missing (no temperature at all)
#          -> data/embeddings_sst_blank/, data/faire/embeddings_sst_blank/
#   l4:    SST gaps filled from OSTIA (tiles_sst_l4)
#          -> data/embeddings_sst_l4/, data/faire/embeddings_sst_l4/
# The current embeddings (data/embeddings, data/faire/embeddings) are not touched.
#
#   sbatch --account=191001-364393 cluster/embed_sst_variants.sbatch
#   sbatch --account=191001-364393 --dependency=afterok:<sst_variants job id> cluster/embed_sst_variants.sbatch

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
export HF_HUB_OFFLINE=1

echo "Job started: $(date)"
for set in "data/manifests data" "data/faire/manifests data/faire"; do
  read -r manifests root <<< "${set}"
  python scripts/embed_tiles.py --local-files-only --blank-sst \
    --manifest-csv "${manifests}/biogoship_tile_manifest.csv" --out-dir "${root}/embeddings_sst_blank"
  python scripts/embed_tiles.py --local-files-only \
    --manifest-csv "${manifests}/biogoship_tile_manifest_sst_l4.csv" --out-dir "${root}/embeddings_sst_l4"
done
echo "Job finished: $(date)"
echo "Bring back: data/embeddings_sst_blank/granite_embeddings_bundle.zip, data/embeddings_sst_l4/granite_embeddings_bundle.zip,"
echo "  the same two under data/faire/, data/sst_l4_log.csv, data/faire/sst_l4_log.csv, data/faire/embeddings/tile_raw_features.csv"
