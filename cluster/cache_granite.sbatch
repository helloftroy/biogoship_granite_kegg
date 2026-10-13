#!/bin/bash
#SBATCH --job-name=bgs_granite_cache
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# One-off, needs internet: install TerraTorch (and with it PyTorch) into the project
# environment, download the Granite checkpoint into the repo's Hugging Face cache,
# and check that the encoder loads. GPU nodes have no internet, so run this first.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/cache_granite.sbatch

set -euo pipefail
export BGS_CONDA_ENV="${BGS_CONDA_ENV:-/scratch/morrill/users/${USER}/conda_envs/biogoship}"
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"

echo "Job started: $(date)"
if ! python -c 'import torch, terratorch, huggingface_hub' 2>/dev/null; then
  echo "Installing terratorch and huggingface_hub into ${BGS_CONDA_ENV}..."
  python -m pip install --quiet terratorch huggingface_hub
fi
python - <<'PY'
import sys
sys.path.insert(0, ".")
import torch
from biogoship import granite

encoder, report = granite.build_encoder(local_files_only=False)
print(report)
with torch.no_grad():
    out = granite.encode(encoder.eval(), torch.zeros(2, 17, 42, 42))
print("torch", torch.__version__, "| test embedding shape", tuple(out.shape))
assert tuple(out.shape) == (2, granite.EMBED_DIM), "unexpected embedding shape"
PY
echo "Model cached under ${HF_HOME}"
echo "Job finished: $(date)"
