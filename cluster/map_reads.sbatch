#!/bin/bash
#SBATCH --job-name=bgs_map_reads
#SBATCH --partition=morrill
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# Maps the external test sets' reads to the Bio-GO-SHIP reference with the commands of
# their mapping pipeline (biogoship/biogoship-metagenomics @ 402e528, mapping_pipeline/;
# RESULTS.md R30, R32). One array task per sequencing run.
#
#   their rule      here
#   trimmomatic     trimmomatic PE ... ILLUMINACLIP:adapters.fa:2:30:7 LEADING:2 TRAILING:2
#                   SLIDINGWINDOW:4:2 MINLEN:50 (unpaired survivors dropped, as in theirs)
#   minimap2        minimap2 -I32g -ax sr <v2.mmi> paired_1 paired_2 --split-prefix ...
#   filter_sam      their scripts/filter_sam_header.py, unchanged
#   sort_bam        samtools view -bh | samtools sort
#   coverage        samtools coverage;  flagstat: samtools flagstat
# Differences, none of which changes a read's alignment:
#   - Thread counts (theirs: -t 16 and single-threaded sort).
#   - The trimmed pairs are mapped in chunks of CHUNK_PAIRS pairs (default 50 M), and the
#     chunk BAMs are merged before samtools coverage. Each chunk is a checkpoint, so a
#     large run that hits the 48 h limit resumes at the next chunk when re-submitted.
#     minimap2 maps each pair on its own; the chunk boundary is one more mini-batch
#     boundary among many.
#   - Plain bash steps instead of Snakemake.
# The counting (their count_biogoship_annotations.py and count_biogoship_taxa.py, at
# coverage 0 and 50) is cluster/count_reads.sbatch, on few cores.
#
# Needs: cluster/setup_mapping_env.sbatch run once (tools + their repo copy), the
# reference (cluster/download_reference.sbatch) and the reads (cluster/download_reads.sbatch).
#
# 1. List the downloaded runs (login node, seconds). Smallest first; writes
#    logs/map_runs_<stamp>.tsv and prints the sbatch commands:
#      cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#      bash cluster/map_reads.sbatch --list
# 2. Pilot on the smallest run, then all (at most 20 at a time):
#      sbatch --account=191001-364393 --array=0 --export=ALL,RUNLIST=logs/map_runs_<stamp>.tsv cluster/map_reads.sbatch
#      sbatch --account=191001-364393 --array=0-<n-1>%20 --export=ALL,RUNLIST=... cluster/map_reads.sbatch
# Resumable: re-submit the same command; finished steps and runs are skipped.
# Outputs per run, data/external/<set>/mapping/<run>/:
#   <run>.sorted.bam (KEEP_BAM=0 deletes it), <run>.coverage, <run>.flagstat,
#   summary.tsv, trimmomatic.log, minimap2_chunk_*.log; done = finished.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"

if [ "${1:-}" = "--list" ]; then
  SETS="${SETS:-tara_pacific jgi_flinkstrom emo_bon}"
  stamp="$(date +%Y%m%d_%H%M%S)"
  list="logs/map_runs_${stamp}.tsv"
  python - "${list}" ${SETS//,/ } <<'EOF'
import os
import sys
import pandas as pd
out, sets = sys.argv[1], sys.argv[2:]
rows, waiting = [], 0
for s in sets:
    runs = pd.read_csv(f"data/external/{s}/reads/runs.csv", dtype=str)
    for _, r in runs.iterrows():
        files = [f"data/external/{s}/reads/fastq/" + u.rsplit("/", 1)[1] for u in r["fastq_ftp"].split(";")]
        r1 = [f for f in files if f.endswith("_1.fastq.gz")]
        r2 = [f for f in files if f.endswith("_2.fastq.gz")]
        if len(r1) != 1 or len(r2) != 1:
            print(f"  {s} {r['run_accession']}: not a _1/_2 pair, skipped", flush=True)
            continue
        if not (os.path.exists(r1[0] + ".ok") and os.path.exists(r2[0] + ".ok")):
            waiting += 1
            continue
        rows.append((s, r["sample_id"], r["run_accession"], os.path.abspath(r1[0]), os.path.abspath(r2[0]), int(r["base_count"])))
rows.sort(key=lambda x: x[5])
pd.DataFrame(rows).to_csv(out, sep="\t", header=False, index=False)
done = sum(os.path.exists(f"data/external/{r[0]}/mapping/{r[2]}/done") for r in rows)
print(f"{len(rows)} runs downloaded and listed (already mapped: {done}); {waiting} still downloading")
if rows:
    print(f"smallest: {rows[0][2]} ({rows[0][0]}, {rows[0][5] / 1e9:.1f} Gbp); largest: {rows[-1][2]} ({rows[-1][5] / 1e9:.1f} Gbp)")
EOF
  n="$(wc -l < "${list}" | tr -d ' ')"
  echo "List: ${list}"
  echo "Pilot:  sbatch --account=191001-364393 --array=0 --export=ALL,RUNLIST=${list} cluster/map_reads.sbatch"
  echo "All:    sbatch --account=191001-364393 --array=0-$((n - 1))%20 --export=ALL,RUNLIST=${list} cluster/map_reads.sbatch"
  exit 0
fi

RUNLIST="${RUNLIST:?Set RUNLIST to the list from: bash cluster/map_reads.sbatch --list}"
MAP_ENV="${MAP_ENV:-${SCRATCH_USER_DIR}/conda_envs/bgs-mapping}"
PIPE="${BGS_PIPELINE:-${SCRATCH_USER_DIR}/biogoship-metagenomics/mapping_pipeline}"
INDEX="${INDEX:-${REPO_ROOT}/data/reference/biogoship_2026/biogoship_2026_database_v2.mmi}"
MAX_MAPPING_MEMORY="${MAX_MAPPING_MEMORY:-32g}"  # their config.yaml max_mapping_memory
CHUNK_PAIRS="${CHUNK_PAIRS:-50000000}"
KEEP_BAM="${KEEP_BAM:-1}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"
set +u; source "${CONDA_SH}"; conda activate "${MAP_ENV}"; set -u
for f in "${INDEX}" "${PIPE}/resources/adapters.fa" "${PIPE}/scripts/filter_sam_header.py"; do
  [ -f "${f}" ] || { echo "Missing ${f}" >&2; exit 2; }
done

line="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${RUNLIST}")"
[ -n "${line}" ] || { echo "No line $((SLURM_ARRAY_TASK_ID + 1)) in ${RUNLIST}" >&2; exit 2; }
IFS=$'\t' read -r SET SAMPLE RUN R1 R2 BASES <<< "${line}"
OUT="${REPO_ROOT}/data/external/${SET}/mapping/${RUN}"
WORK="${OUT}/work"
mkdir -p "${WORK}"
step() { echo "$(date '+%F %T') ${RUN}: $*"; }
step "set ${SET}, sample ${SAMPLE}, $((BASES / 1000000000)) Gbp, ${THREADS} threads, node $(hostname)"
if [ -f "${OUT}/done" ]; then step "already done"; exit 0; fi
t0=$(date +%s)

# 1. Trim (their trimmomatic rule)
if [ ! -f "${WORK}/trim.ok" ]; then
  step "trimmomatic"
  trimmomatic PE -threads "${THREADS}" "${R1}" "${R2}" \
    "${WORK}/paired_1.fastq.gz" "${WORK}/unpaired_1.fastq.gz" "${WORK}/paired_2.fastq.gz" "${WORK}/unpaired_2.fastq.gz" \
    ILLUMINACLIP:"${PIPE}/resources/adapters.fa":2:30:7 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:50 \
    2> "${OUT}/trimmomatic.log"
  rm -f "${WORK}"/unpaired_*.fastq.gz
  touch "${WORK}/trim.ok"
fi
grep -h "Input Read Pairs" "${OUT}/trimmomatic.log" || true

# 2. Chunks of CHUNK_PAIRS pairs (both mates split at the same records)
if [ ! -f "${WORK}/split.ok" ]; then
  step "split into chunks of ${CHUNK_PAIRS} pairs"
  rm -f "${WORK}"/chunk_*
  pids=()
  for m in 1 2; do
    pigz -dc -p 4 "${WORK}/paired_${m}.fastq.gz" \
      | split -l $((CHUNK_PAIRS * 4)) -d -a 3 --additional-suffix="_${m}.fastq.gz" \
          --filter='pigz -1 -p 6 > "$FILE"' - "${WORK}/chunk_" &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait "${p}"; done
  ls "${WORK}"/chunk_*_1.fastq.gz | sed 's/_1\.fastq\.gz$//' > "${WORK}/chunks.txt"
  n1=$(ls "${WORK}"/chunk_*_1.fastq.gz | wc -l); n2=$(ls "${WORK}"/chunk_*_2.fastq.gz | wc -l)
  [ "${n1}" -eq "${n2}" ] || { echo "Chunk count differs: ${n1} vs ${n2}" >&2; exit 3; }
  touch "${WORK}/split.ok"
  rm -f "${WORK}"/paired_*.fastq.gz
fi
nchunks=$(wc -l < "${WORK}/chunks.txt" | tr -d ' ')

# 3. Map each chunk (their minimap2, filter_sam and sort_bam rules)
k=0
while read -r c; do
  k=$((k + 1))
  name="$(basename "${c}")"
  [ -f "${c}.bam.ok" ] && continue
  step "minimap2 ${name} (${k}/${nchunks})"
  minimap2 -I"${MAX_MAPPING_MEMORY}" -ax sr "${INDEX}" "${c}_1.fastq.gz" "${c}_2.fastq.gz" -t "${THREADS}" \
    --split-prefix "${c}.tmpindex" > "${c}.sam" 2> "${OUT}/minimap2_${name}.log"
  step "filter_sam_header ${name}"
  python "${PIPE}/scripts/filter_sam_header.py" "${c}.sam" "${c}.sam.filtered"
  rm -f "${c}.sam"
  nsq=$(awk '/^@SQ/{n++; next} !/^@/{exit} END{print n+0}' "${c}.sam.filtered")
  step "sort ${name} (${nsq} reference sequences with reads)"
  samtools view -bh "${c}.sam.filtered" | samtools sort -@ "${THREADS}" -m 2G -T "${c}.sorttmp" -O BAM -o "${c}.sorted.bam"
  rm -f "${c}.sam.filtered"
  touch "${c}.bam.ok"
  rm -f "${c}_1.fastq.gz" "${c}_2.fastq.gz"
done < "${WORK}/chunks.txt"

# 4. Merge, then their coverage and flagstat rules
BAM="${OUT}/${RUN}.sorted.bam"
if [ ! -f "${WORK}/merge.ok" ]; then
  step "merge ${nchunks} chunk BAMs"
  if [ "${nchunks}" -eq 1 ]; then
    mv "$(head -1 "${WORK}/chunks.txt").sorted.bam" "${BAM}"
  else
    samtools merge -f -@ "${THREADS}" -o "${BAM}" $(sed 's/$/.sorted.bam/' "${WORK}/chunks.txt")
  fi
  touch "${WORK}/merge.ok"
  rm -f "${WORK}"/chunk_*.sorted.bam
fi
step "samtools coverage, flagstat"
samtools coverage "${BAM}" > "${OUT}/${RUN}.coverage.part" && mv "${OUT}/${RUN}.coverage.part" "${OUT}/${RUN}.coverage"
samtools flagstat -@ "${THREADS}" "${BAM}" > "${OUT}/${RUN}.flagstat"

pairs_in=$(grep -o "Input Read Pairs: [0-9]*" "${OUT}/trimmomatic.log" | grep -o "[0-9]*$" || echo NA)
pairs_kept=$(grep -o "Both Surviving: [0-9]*" "${OUT}/trimmomatic.log" | grep -o "[0-9]*$" || echo NA)
primary=$(grep -m1 "primary$" "${OUT}/${RUN}.flagstat" | cut -d' ' -f1 || echo NA)
primary_mapped=$(grep -m1 "primary mapped" "${OUT}/${RUN}.flagstat" | cut -d' ' -f1 || echo NA)
refs=$(awk 'NR > 1 && $4 > 0' "${OUT}/${RUN}.coverage" | wc -l | tr -d ' ')
refs50=$(awk 'NR > 1 && $4 > 0 && $6 >= 50' "${OUT}/${RUN}.coverage" | wc -l | tr -d ' ')
printf "set\tsample_id\trun\tbases\tpairs_in\tpairs_after_trim\tprimary_reads\tprimary_mapped\trefs_with_reads\trefs_cov50\tchunks\thours_last_job\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n" \
  "${SET}" "${SAMPLE}" "${RUN}" "${BASES}" "${pairs_in}" "${pairs_kept}" "${primary}" "${primary_mapped}" "${refs}" "${refs50}" \
  "${nchunks}" "$(awk -v a="${t0}" -v b="$(date +%s)" 'BEGIN { print (b - a) / 3600 }')" > "${OUT}/summary.tsv"
cat "${OUT}/summary.tsv"
[ "${KEEP_BAM}" = "1" ] || rm -f "${BAM}"
rm -rf "${WORK}"
touch "${OUT}/done"
step "done"
