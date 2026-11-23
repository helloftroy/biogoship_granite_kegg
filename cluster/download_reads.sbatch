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
# sizes, MD5s). With OMDB=1 also the OMDB v2.0 NR95 gene catalog (no longer needed).
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
#   sbatch --account=191001-364393 --export=ALL,SETS=tara_pacific,PARALLEL=2 cluster/download_reads.sbatch
#   # if https keeps failing: --export=ALL,PROTO=ftp
#   # only the first 20 M read pairs per run (Bio-GO-SHIP libraries: median 14 M pairs,
#   # 10-90% 5.7-41 M), about 0.63 TB instead of 3.5 TB; already complete files are cut
#   # locally and the originals moved to fastq_full/:
#   sbatch --account=191001-364393 --export=ALL,SETS=emo_bon,SUBSAMPLE_PAIRS=20000000 cluster/download_reads.sbatch
# Progress: grep -c . logs/download_reads_done.txt; failures in logs/download_reads_failed.txt

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
SETS="${SETS:-tara_pacific jgi_flinkstrom emo_bon}"
SETS="${SETS//,/ }"
PARALLEL="${PARALLEL:-3}"
PROTO="${PROTO:-https}"  # or ftp, if ENA keeps cutting https transfers
MAX_STALLS="${MAX_STALLS:-8}"
SUBSAMPLE_PAIRS="${SUBSAMPLE_PAIRS:-0}"  # >0: only the first N read pairs per run (see fetch_sub)
OMDB="${OMDB:-0}"  # 1 = also fetch OMDB v2.0 NR95 (not needed now: the Bio-GO-SHIP reference is in cluster/download_reference.sbatch)
DONE_LOG="${REPO_ROOT}/logs/download_reads_done.txt"
FAIL_LOG="${REPO_ROOT}/logs/download_reads_failed.txt"
touch "${DONE_LOG}"
: > "${FAIL_LOG}"

echo "Job started: $(date)"
df -h "${REPO_ROOT}" | tail -1

# One line per file: url, md5, bytes, destination
LIST="${TMPDIR}/download_reads_${SLURM_JOB_ID:-local}.tsv"
python - "${LIST}" "${PROTO}" ${SETS} <<'EOF'
import sys
import pandas as pd
out, proto, sets = sys.argv[1], sys.argv[2], sys.argv[3:]
rows = []
for s in sets:
    runs = pd.read_csv(f"data/external/{s}/reads/runs.csv", dtype=str)
    for _, r in runs.iterrows():
        for url, md5, size in zip(r["fastq_ftp"].split(";"), r["fastq_md5"].split(";"), r["fastq_bytes"].split(";")):
            rows.append((f"{proto}://" + url, md5, size, f"data/external/{s}/reads/fastq/" + url.rsplit("/", 1)[1]))
    total = sum(int(x) for v in runs["fastq_bytes"] for x in v.split(";"))
    print(f"{s}: {len(runs)} runs, {total / 1e12:.2f} TB", flush=True)
pd.DataFrame(rows).to_csv(out, sep="\t", header=False, index=False)
EOF

# ENA closes long HTTPS transfers part-way ("curl: (18) transfer closed") and answers
# 403 to bursts of requests. So: resume (-C -) for as long as the file keeps growing,
# back off between attempts (longer after a 403), and give up on a file only after
# MAX_STALLS attempts in a row without progress. A full-size file is MD5-checked.
fetch() {
  local url="$1" md5="$2" size="$3" dest="$4"
  [ -f "${dest}.ok" ] && return 0
  mkdir -p "$(dirname "${dest}")"
  local part="${dest}.part" stalls=0 bad_md5=0 before after rc
  while [ "${stalls}" -lt "${MAX_STALLS}" ]; do
    before=$(stat -c %s "${part}" 2>/dev/null || echo 0)
    if [ "${before}" -lt "${size}" ]; then
      rc=0
      curl -fsS -L --connect-timeout 60 --speed-limit 50000 --speed-time 600 -C - -o "${part}" "${url}" 2> "${part}.err" || rc=$?
      after=$(stat -c %s "${part}" 2>/dev/null || echo 0)
    else
      rc=0; after="${before}"
    fi
    if [ "${after}" -ge "${size}" ]; then
      if [ "${after}" -eq "${size}" ] && [ "$(md5sum "${part}" | cut -d' ' -f1)" = "${md5}" ]; then
        mv "${part}" "${dest}" && touch "${dest}.ok" && rm -f "${part}.err"
        echo "${dest}" >> "${DONE_LOG}"
        return 0
      fi
      echo "MD5 or size mismatch, starting over: ${dest} (${after} of ${size} bytes)" >&2
      rm -f "${part}"
      bad_md5=$((bad_md5 + 1))
      [ "${bad_md5}" -ge 2 ] && break
      continue
    fi
    if [ "${after}" -gt "${before}" ]; then
      stalls=0
      sleep 10
    else
      stalls=$((stalls + 1))
      if grep -q "error: 403" "${part}.err" 2>/dev/null; then
        echo "403 (throttled?), waiting $((300 * stalls)) s: ${url}" >&2
        sleep $((300 * stalls))
      else
        echo "no progress (curl exit ${rc}: $(tr '\n' ' ' < "${part}.err")), waiting $((60 * stalls)) s: ${url}" >&2
        sleep $((60 * stalls))
      fi
    fi
  done
  echo -e "${url}\t${dest}" >> "${FAIL_LOG}"
  return 0
}
# SUBSAMPLE_PAIRS=N: keep only the first N reads of each file (N pairs per run). The
# stream is cut after 4N lines and re-compressed, so ENA's MD5 cannot apply. Checks
# instead: exactly 4N lines, or, for a file with fewer reads, the whole stream ended
# cleanly (curl exit 0 and gzip exit 0, so the gzip CRC and length passed). A file
# that is already complete on disk is cut locally and the full file moved to
# fastq_full/; a .part from an earlier full download is used if it is long enough.
# Done files get <file>.ok plus <file>.sub (holds N), so the mapping job sees them.
cut_stream() {  # stdin: fastq.gz bytes; $1 lines wanted, $2 output, $3 count file
  local zip="gzip -1"
  command -v pigz > /dev/null && zip="pigz -1 -p 2"
  gzip -dc 2> /dev/null | awk -v n="$1" -v cf="$3" '{print} NR>=n{exit} END{print NR > cf}' | ${zip} > "$2"
  local st=("${PIPESTATUS[@]}")
  echo "${st[0]}" > "$3.gz"  # gzip exit: 0 only if the stream ended with a valid CRC
}
fetch_sub() {
  local url="$1" md5="$2" size="$3" dest="$4"
  local lines=$((SUBSAMPLE_PAIRS * 4)) tmp="${dest}.subtmp" cnt="${dest}.subcount"
  [ -f "${dest}.ok" ] && [ -f "${dest}.sub" ] && [ "$(cat "${dest}.sub")" = "${SUBSAMPLE_PAIRS}" ] && return 0
  mkdir -p "$(dirname "${dest}")"
  local src="" full_ok=0
  if [ -f "${dest}.ok" ] && [ ! -f "${dest}.sub" ]; then
    local fulldir="$(dirname "$(dirname "${dest}")")/fastq_full"
    mkdir -p "${fulldir}"
    mv "${dest}" "${fulldir}/" && mv "${dest}.ok" "${fulldir}/"
    src="${fulldir}/$(basename "${dest}")"; full_ok=1
  elif [ -f "${dest}.part" ]; then
    src="${dest}.part"
  fi
  if [ -n "${src}" ]; then
    rm -f "${cnt}"
    cut_stream "${lines}" "${tmp}" "${cnt}" < "${src}"
    local got="$(cat "${cnt}" 2> /dev/null || echo 0)"
    # a complete file with fewer reads than N passes whole (its MD5 was checked)
    if [ "${got}" -eq "${lines}" ] || { [ "${full_ok}" = 1 ] && [ $((got % 4)) -eq 0 ] && [ "${got}" -gt 0 ] && [ "$(cat "${cnt}.gz")" = 0 ]; }; then
      mv "${tmp}" "${dest}" && echo "${SUBSAMPLE_PAIRS}" > "${dest}.sub" && touch "${dest}.ok"
      rm -f "${cnt}" "${cnt}.gz" "${dest}.part" "${dest}.part.err"
      echo "${dest}" >> "${DONE_LOG}"
      return 0
    fi
    rm -f "${tmp}" "${dest}.part" "${dest}.part.err"
  fi
  local tries=0 rc
  while [ "${tries}" -lt "${MAX_STALLS}" ]; do
    rm -f "${cnt}" "${tmp}"
    curl -fsS -L --connect-timeout 60 --speed-limit 50000 --speed-time 600 "${url}" 2> "${tmp}.err" \
      | cut_stream "${lines}" "${tmp}" "${cnt}"
    rc=("${PIPESTATUS[@]}")
    local got="$(cat "${cnt}" 2> /dev/null || echo 0)"
    # cut at 4N lines: curl and gzip are stopped by the closed pipe, which is expected
    if [ "${got}" -eq "${lines}" ] || { [ "${got}" -gt 0 ] && [ $((got % 4)) -eq 0 ] && [ "${rc[0]}" -eq 0 ] && [ "$(cat "${cnt}.gz")" = 0 ]; }; then
      mv "${tmp}" "${dest}" && echo "${SUBSAMPLE_PAIRS}" > "${dest}.sub" && touch "${dest}.ok"
      rm -f "${cnt}" "${cnt}.gz" "${tmp}.err"
      echo "${dest}" >> "${DONE_LOG}"
      return 0
    fi
    tries=$((tries + 1))
    if grep -q "error: 403" "${tmp}.err" 2> /dev/null; then
      echo "403 (throttled?), waiting $((300 * tries)) s: ${url}" >&2; sleep $((300 * tries))
    else
      echo "incomplete (${got} of ${lines} lines; curl exit ${rc[0]}: $(tr '\n' ' ' < "${tmp}.err")), waiting $((60 * tries)) s: ${url}" >&2
      sleep $((60 * tries))
    fi
  done
  rm -f "${tmp}" "${cnt}" "${cnt}.gz" "${tmp}.err"
  echo -e "${url}\t${dest}" >> "${FAIL_LOG}"
  return 0
}
export -f fetch fetch_sub cut_stream
export DONE_LOG FAIL_LOG MAX_STALLS SUBSAMPLE_PAIRS

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
FETCH=fetch
[ "${SUBSAMPLE_PAIRS}" -gt 0 ] && FETCH=fetch_sub && echo "Subsample mode: first ${SUBSAMPLE_PAIRS} read pairs per run"
xargs -a "${LIST}" -P "${PARALLEL}" -L 1 bash -c "${FETCH}"' "$0" "$1" "$2" "$3"'

n_ok=$(wc -l < "${DONE_LOG}"); n_fail=$(wc -l < "${FAIL_LOG}"); n_all=$(wc -l < "${LIST}")
echo "Files verified so far: ${n_ok} (this list: ${n_all}); failed this run: ${n_fail}"
df -h "${REPO_ROOT}" | tail -1
echo "Job finished: $(date)"
echo "If any failed or the job hit 48 h, re-submit the same command; finished files are skipped."
