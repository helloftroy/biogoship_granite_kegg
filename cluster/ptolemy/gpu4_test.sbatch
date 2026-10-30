#!/bin/bash
#SBATCH --job-name=gpu4_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Can this account get 4 GPUs on one node at once? Needs no Python.
# Ptolemy: partition gpu-a100 (4 nodes x 8 whole A100s); QOS "normal" allows at most
# 4 A100s per user, 128 CPUs, 1 TB memory and 2 days per job. Account 191001-366746
# (not Morrill's 191001-364393).
#
#   cd /scratch/ptolemy/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-366746 cluster/ptolemy/gpu4_test.sbatch
#
# With a PyTorch environment, set TORCH_ENV=/path/to/env to also run a 4-GPU NCCL check.

echo "Job $SLURM_JOB_ID on $(hostname), partition $SLURM_JOB_PARTITION, account $SLURM_JOB_ACCOUNT"
echo "Started: $(date)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES  SLURM_GPUS_ON_NODE=${SLURM_GPUS_ON_NODE:-?}"
echo
nvidia-smi --query-gpu=index,name,memory.total,pci.bus_id --format=csv
n=$(nvidia-smi -L | wc -l)
echo
echo "GPUs visible to this job: $n"
echo
echo "Interconnect between them (NV# = NVLink, SYS/PHB = PCIe):"
nvidia-smi topo -m
echo
scontrol show job "$SLURM_JOB_ID" | grep -E "TRES|TimeLimit|Partition"

if [ -n "${TORCH_ENV:-}" ]; then
  echo; echo "NCCL check across $n GPUs with $TORCH_ENV"
  # shellcheck disable=SC1091
  if [ -f "$TORCH_ENV/bin/activate" ]; then source "$TORCH_ENV/bin/activate"; fi
  torchrun --standalone --nproc_per_node="$n" "${SLURM_SUBMIT_DIR}/cluster/ptolemy/ddp_check.py"
fi
echo "Finished: $(date)"
[ "$n" -ge 4 ] && echo "RESULT: 4 GPUs OK" || echo "RESULT: only $n GPU(s)"
