#!/bin/bash
#SBATCH --job-name=bgs_download_reads
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Raw metagenome reads of the external test sets (RESULTS.md R31), from the lists that
# scripts/build_read_manifests.py wrote (data/external/<set>/reads/runs.csv; ENA URLs,
# sizes, MD5s), plus the OMDB v2.0 NR95 gene catalog as the interim mapping reference.
#
#   set             samples  runs   size
#   tara_pacific         84    84   1.74 TB
#   jgi_flinkstrom       90    90   0.53 TB
#   emo_bon             131   131   1.19 TB
#   OMDB v2.0 NR95                  27 GB fasta + 3.4 GB cluster table
#   total                           about 3.5 TB: check the scratch quota first
#
# Files go to data/external/<set>/reads/fastq/<run>_{1,2}.fastq.gz and
# data/reference/omdb_v2/. Every file is checked against the ENA MD5; a file that
# passes gets a .ok marker and is never fetched again. Partial files (.part) resume.
# Resumable: after a stop or the 48 h limit, re-submit the same command.
#
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/download_reads.sbatch
#   # one set only, or more parallel downloads:
#   sbatch --account=191001-364393 --export=ALL,SETS=tara_pacific,PARALLEL=6 cluster/download_reads.sbatch
# Progress: grep -c . logs/download_reads_done.txt; failures in logs/download_reads_failed.txt

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
SETS="${SETS:-tara_pacific jgi_flinkstrom emo_bon}"
SETS="${SETS//,/ }"
PARALLEL="${PARALLEL:-4}"
OMDB="${OMDB:-1}"
DONE_LOG="${REPO_ROOT}/logs/download_reads_done.txt"
FAIL_LOG="${REPO_ROOT}/logs/download_reads_failed.txt"
touch "${DONE_LOG}"
: > "${FAIL_LOG}"

echo "Job started: $(date)"
df -h "${REPO_ROOT}" | tail -1

# One line per file: url, md5, bytes, destination
LIST="${TMPDIR}/download_reads_${SLURM_JOB_ID:-local}.tsv"
python - "${LIST}" ${SETS} <<'EOF'
import sys
import pandas as pd
out, sets = sys.argv[1], sys.argv[2:]
rows = []
for s in sets:
    runs = pd.read_csv(f"data/external/{s}/reads/runs.csv", dtype=str)
    for _, r in runs.iterrows():
        for url, md5, size in zip(r["fastq_ftp"].split(";"), r["fastq_md5"].split(";"), r["fastq_bytes"].split(";")):
            rows.append(("https://" + url, md5, size, f"data/external/{s}/reads/fastq/" + url.rsplit("/", 1)[1]))
    total = sum(int(x) for v in runs["fastq_bytes"] for x in v.split(";"))
    print(f"{s}: {len(runs)} runs, {total / 1e12:.2f} TB", flush=True)
pd.DataFrame(rows).to_csv(out, sep="\t", header=False, index=False)
EOF

fetch() {
  local url="$1" md5="$2" size="$3" dest="$4"
  [ -f "${dest}.ok" ] && return 0
  mkdir -p "$(dirname "${dest}")"
  local k
  for k in 1 2 3; do
    if curl -fsSL --retry 10 --retry-delay 30 --retry-connrefused -C - -o "${dest}.part" "${url}"; then
      if [ "$(md5sum "${dest}.part" | cut -d' ' -f1)" = "${md5}" ]; then
        mv "${dest}.part" "${dest}" && touch "${dest}.ok"
        echo "${dest}" >> "${DONE_LOG}"
        return 0
      fi
      echo "MD5 mismatch (attempt ${k}): ${dest}" >&2
      rm -f "${dest}.part"
    else
      echo "curl failed (attempt ${k}): ${url}" >&2
    fi
  done
  echo -e "${url}\t${dest}" >> "${FAIL_LOG}"
  return 0
}
export -f fetch
export DONE_LOG FAIL_LOG

# The reference first (small, and needed before any mapping)
if [ "${OMDB}" = "1" ]; then
  REF="${REPO_ROOT}/data/reference/omdb_v2"
  B="https://sunagawalab.ethz.ch/share/microbiomics/ocean/db/2.0/data/catalogs"
  mkdir -p "${REF}"
  for f in OMDBv2.0_NT_G_NR95/OMDBv2.0_NT_G_NR95.fna.gz OMDBv2.0_NT_G_NR95/OMDBv2.0_NT_G_NR95.cluster.tsv.gz OMDBv2.0_data.tsv.gz; do
    name="$(basename "${f}")"
    [ -f "${REF}/${name}.ok" ] && continue
    curl -fsSL --retry 10 --retry-delay 30 --retry-connrefused -C - -o "${REF}/${name}.part" "${B}/${f}" \
      && mv "${REF}/${name}.part" "${REF}/${name}" && curl -fsSL -o "${REF}/${name}.md5" "${B}/${f}.md5" \
      && touch "${REF}/${name}.ok" && echo "OMDB: ${name}"
  done
  # The published fasta MD5 is of the uncompressed file
  if [ -f "${REF}/OMDBv2.0_NT_G_NR95.fna.gz.ok" ] && [ ! -f "${REF}/fasta_md5_checked" ]; then
    want="$(cut -d' ' -f1 "${REF}/OMDBv2.0_NT_G_NR95.fna.gz.md5")"
    got="$(gzip -dc "${REF}/OMDBv2.0_NT_G_NR95.fna.gz" | md5sum | cut -d' ' -f1)"
    if [ "${want}" = "${got}" ]; then touch "${REF}/fasta_md5_checked"; echo "OMDB fasta MD5 ok"; else echo "OMDB fasta MD5 MISMATCH" >&2; fi
  fi
fi

# Reads, in the order of SETS; PARALLEL files at a time
xargs -a "${LIST}" -P "${PARALLEL}" -L 1 bash -c 'fetch "$0" "$1" "$2" "$3"'

n_ok=$(wc -l < "${DONE_LOG}"); n_fail=$(wc -l < "${FAIL_LOG}"); n_all=$(wc -l < "${LIST}")
echo "Files verified so far: ${n_ok} (this list: ${n_all}); failed this run: ${n_fail}"
df -h "${REPO_ROOT}" | tail -1
echo "Job finished: $(date)"
echo "If any failed or the job hit 48 h, re-submit the same command; finished files are skipped."
