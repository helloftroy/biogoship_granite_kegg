#!/bin/bash
# What this cluster offers this account: storage, quotas, partitions, GPUs, accounts.
# Runs on a login node in a few seconds; needs no Python. Paste the output back.
#
#   bash cluster/ptolemy/check_cluster.sh 2>&1 | tee logs/ptolemy_check.txt

SCRATCH="${SCRATCH_DIR:-/scratch/ptolemy/users/$USER}"
section() { printf '\n===== %s =====\n' "$1"; }
have() { command -v "$1" >/dev/null 2>&1; }

section "Host and user"
hostname; whoami; date

section "Scratch: $SCRATCH"
if [ -d "$SCRATCH" ]; then
  df -h "$SCRATCH"
  echo "Filesystem type: $(stat -f -c %T "$SCRATCH" 2>/dev/null || echo unknown)"
  echo "Used by me (may take a minute on a big tree):"
  timeout 300 du -sh "$SCRATCH" 2>/dev/null || echo "  du timed out or failed"
else
  echo "Not found. Scratch directories here:"; ls -d /scratch/*/ 2>/dev/null
fi

section "Home"
df -h "$HOME"

section "Quotas (whichever tool this filesystem uses)"
have quota && { echo "-- quota -s"; quota -s 2>&1; }
have lfs && { echo "-- lfs quota (Lustre)"; lfs quota -h -u "$USER" "$SCRATCH" 2>&1; }
have mmlsquota && { echo "-- mmlsquota (GPFS)"; mmlsquota -u "$USER" --block-size auto 2>&1; }
for tool in showquota myquota quota_check checkquota; do
  have "$tool" && { echo "-- $tool"; "$tool" 2>&1; }
done
echo "(If none printed a limit, ask HPC2 support for the scratch quota and purge policy.)"

section "Purge policy hints"
ls /etc/motd* 2>/dev/null && grep -i -E "purge|scratch|quota" /etc/motd* 2>/dev/null | head -20

section "SLURM accounts I can charge"
have sacctmgr && sacctmgr -n -P show assoc user="$USER" format=account,partition,qos,grptres,maxtres,maxwall 2>&1 | head -40

section "Partitions: name, time limit, nodes, CPUs, memory (MB), GPUs (gres), state"
have sinfo && sinfo -o "%15P %11l %6D %5c %8m %30G %10a" 2>&1

section "GPU nodes in detail"
have sinfo && sinfo -N -o "%15N %15P %30G %8m %10T" 2>&1 | grep -i -E "gpu|NODELIST" | sort -u | head -60

section "QOS limits"
have sacctmgr && sacctmgr -n -P show qos format=name,maxwall,maxtrespu,maxjobspu,maxsubmitpu 2>&1 | head -30

section "Modules that may matter (cuda, python, conda, gcc)"
if have module || type module >/dev/null 2>&1; then
  module -t avail 2>&1 | grep -i -E "^(cuda|python|anaconda|miniconda|miniforge|gcc|nccl|cudnn)" | head -40
else
  echo "no module command"
fi

section "Internet from the login node"
curl -s -o /dev/null -w "github %{http_code}\n" --max-time 15 https://github.com || echo "github: no"
curl -s -o /dev/null -w "nasa earthdata %{http_code}\n" --max-time 15 https://urs.earthdata.nasa.gov || echo "earthdata: no"
