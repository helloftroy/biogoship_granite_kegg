#!/bin/bash
#SBATCH --job-name=bgs_count_reads
#SBATCH --partition=morrill
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# KO, Pfam, OG and taxon read counts for one mapped run (array task = line of the same
# RUNLIST as cluster/map_reads.sbatch), with the Bio-GO-SHIP pipeline's own scripts,
# unchanged (biogoship-metagenomics @ 402e528; RESULTS.md R30, R32):
#   count_biogoship_annotations.py <run>.coverage eggnogmapper.tsv taxonomy.tsv OUT --coverage_threshold T
#   count_biogoship_taxa.py        <run>.coverage taxonomy.tsv OUT --coverage_threshold T
# at T = 0 (their Snakefile default, "coverage0") and T = 50 (the "cov50" tables our
# Bio-GO-SHIP targets come from). The four run side by side; each reads the 105 GB
# eggNOG table and/or the 29 GB taxonomy once, so a task takes hours, not minutes.
# Outputs, gzipped, in the layout scripts/collect_taxa_tables.py reads:
#   data/external/<set>/ko_pipeline/cov<T>/<run>/level_<n>_<rank>_{ko,pfam,og}_counts.tsv.gz
#   data/external/<set>/ko_pipeline/cov<T>/<run>/level_<n>_<rank>_abundance.tsv.gz
#
# Runs that are not mapped yet are skipped (re-submit later). To start each run as
# soon as its mapping task finishes:
#   sbatch --account=191001-364393 --array=0-<n-1>%20 --dependency=aftercorr:<map job id> \
#       --export=ALL,RUNLIST=logs/map_runs_<stamp>.tsv cluster/count_reads.sbatch
# or, after mapping:
#   sbatch --account=191001-364393 --array=0-<n-1>%20 --export=ALL,RUNLIST=... cluster/count_reads.sbatch

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
RUNLIST="${RUNLIST:?Set RUNLIST to the list from: bash cluster/map_reads.sbatch --list}"
MAP_ENV="${MAP_ENV:-${SCRATCH_USER_DIR}/conda_envs/bgs-mapping}"
PIPE="${BGS_PIPELINE:-${SCRATCH_USER_DIR}/biogoship-metagenomics/mapping_pipeline}"
REF="${REPO_ROOT}/data/reference/biogoship_2026"
THRESHOLDS="${THRESHOLDS:-0 50}"
set +u; source "${CONDA_SH}"; conda activate "${MAP_ENV}"; set -u

line="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${RUNLIST}")"
[ -n "${line}" ] || { echo "No line $((SLURM_ARRAY_TASK_ID + 1)) in ${RUNLIST}" >&2; exit 2; }
IFS=$'\t' read -r SET SAMPLE RUN R1 R2 BASES <<< "${line}"
MAPPED="${REPO_ROOT}/data/external/${SET}/mapping/${RUN}"
COV="${MAPPED}/${RUN}.coverage"
step() { echo "$(date '+%F %T') ${RUN}: $*"; }
if [ ! -f "${MAPPED}/done" ]; then step "not mapped yet; skipped"; exit 0; fi

declare -A pids
for t in ${THRESHOLDS}; do
  base="${REPO_ROOT}/data/external/${SET}/ko_pipeline/cov${t}"
  mkdir -p "${base}"
  out="${base}/${RUN}"
  [ -f "${out}/annotations.ok" ] || {
    step "count_biogoship_annotations.py, coverage >= ${t}"
    mkdir -p "${out}"
    python "${PIPE}/scripts/count_biogoship_annotations.py" "${COV}" "${REF}/biogoship_eggnogmapper.tsv" \
      "${REF}/biogoship_2026_taxonomy.tsv" "${out}" --coverage_threshold "${t}" &
    pids["annotations ${t}"]=$!
  }
  [ -f "${out}/taxa.ok" ] || {
    step "count_biogoship_taxa.py, coverage >= ${t}"
    mkdir -p "${out}"
    python "${PIPE}/scripts/count_biogoship_taxa.py" "${COV}" "${REF}/biogoship_2026_taxonomy.tsv" \
      "${out}" --coverage_threshold "${t}" &
    pids["taxa ${t}"]=$!
  }
done
failed=0
for key in "${!pids[@]}"; do
  set -- ${key}
  if wait "${pids[${key}]}"; then
    out="${REPO_ROOT}/data/external/${SET}/ko_pipeline/cov$2/${RUN}"
    if [ "$1" = "annotations" ]; then
      pigz -f -p 4 "${out}"/level_*_{ko,pfam,og}_counts.tsv
    else
      pigz -f -p 4 "${out}"/level_*_abundance.tsv
    fi
    touch "${out}/$1.ok"
    step "${key}: done"
  else
    step "${key}: FAILED" >&2
    failed=1
  fi
done
[ "${failed}" -eq 0 ] || exit 1
for t in ${THRESHOLDS}; do
  out="${REPO_ROOT}/data/external/${SET}/ko_pipeline/cov${t}/${RUN}"
  step "coverage >= ${t}: $(gzip -dc "${out}"/level_1_*_ko_counts.tsv.gz | awk -F'\t' 'NR > 1 { n++; s += $4 } END { printf "%d domain x KO rows, %.0f reads on KOs", n, s }')"
done
