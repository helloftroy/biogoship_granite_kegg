#!/bin/bash
#SBATCH --job-name=bgs_setup_mapping
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Tools for the Bio-GO-SHIP read-mapping pipeline (RESULTS.md R30, R32), and a check
# of the copy of their repository that cluster/map_reads.sbatch calls.
#
# 1. Environment /scratch/morrill/users/$USER/conda_envs/bgs-mapping with the tool
#    versions of their envs/mapping.yaml: minimap2 2.28, samtools 1.18, trimmomatic 0.39.
#    Their yaml pins python 3.6 only for their small scripts, which run unchanged on 3.11.
#    conda writes bookkeeping under $HOME, which is full, so this job points HOME at
#    scratch while it runs conda (the reason CLUSTER.md warns against conda create).
# 2. Their repository, pinned at commit 402e528 and used read-only. It is private, so
#    copy it from the Mac first (about 200 KB):
#      scp -r /Users/hellpark/multimodal_seusmbol/biogoship-metagenomics \
#          hmp278@morrill-dtn-1.hpc.msstate.edu:/scratch/morrill/users/hmp278/
#    The job checks the MD5s of the four files the mapping uses.
# 3. Prints the morrill partition's node sizes (cores, memory) for sizing the mapping.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/setup_mapping_env.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
MAP_ENV="${MAP_ENV:-${SCRATCH_USER_DIR}/conda_envs/bgs-mapping}"
PIPE="${BGS_PIPELINE:-${SCRATCH_USER_DIR}/biogoship-metagenomics/mapping_pipeline}"

echo "Job started: $(date)"
set +u; source "${CONDA_SH}"; set -u
if [ ! -x "${MAP_ENV}/bin/minimap2" ]; then
  FAKE_HOME="${SCRATCH_USER_DIR}/.home_for_conda"
  mkdir -p "${FAKE_HOME}"
  set +u
  HOME="${FAKE_HOME}" conda create -y -p "${MAP_ENV}" --override-channels -c conda-forge -c bioconda \
    minimap2=2.28 samtools=1.18 trimmomatic=0.39 python=3.11 pigz
  set -u
fi
set +u; conda activate "${MAP_ENV}"; set -u
echo "== tools"
minimap2 --version
samtools --version | head -1
trimmomatic -version
python --version

echo "== Bio-GO-SHIP pipeline copy at ${PIPE}"
git -C "${PIPE}/.." log --oneline -1 2>/dev/null || echo "(no git history; checking files only)"
cd "${PIPE}"
cat > "${TMPDIR}/bgs_pipeline.md5" <<'EOF'
eebc3c690999a819d0dcb22792da9d29  resources/adapters.fa
a9e3a1b5989c4c116307e2cf0930976c  scripts/filter_sam_header.py
731434f82ce78cc035d73345bb2906d4  scripts/count_biogoship_annotations.py
aba63cd57f91a1b84a28497b3266e464  scripts/count_biogoship_taxa.py
EOF
md5sum -c "${TMPDIR}/bgs_pipeline.md5"
cd "${REPO_ROOT}"

echo "== reference"
ls -la data/reference/biogoship_2026/
echo "== morrill nodes (name, cores, memory MB, state)"
sinfo -p morrill -N -h -o "%N %c %m %t" | sort -k2,3 -u | head -20 || true
echo "Job finished: $(date)"
echo "Paste back: this .out file."
