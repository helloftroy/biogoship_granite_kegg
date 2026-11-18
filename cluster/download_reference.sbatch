#!/bin/bash
#SBATCH --job-name=bgs_download_reference
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# The Bio-GO-SHIP mapping reference (RESULTS.md R30, R32) from the team's Google Drive
# folder to data/reference/biogoship_2026/, with rclone:
#   biogoship_2026_database_v2.mmi   minimap2 index (the one config.yaml names;
#                                    the *_32G.mmi variant is not needed)
#   biogoship_2026_taxonomy.tsv      sequence -> method, domain ... species
#   biogoship_eggnogmapper.tsv       eggNOG-mapper annotations per sequence
# rclone checks each file's MD5 against Drive after the transfer. Finished files are
# skipped when re-submitted.
#
# One-time setup (rclone needs your Google login; the token stays in
# ~/.config/rclone/rclone.conf on Morrill):
#   1. On Morrill:  bash cluster/download_reference.sbatch --install   (puts rclone in
#      /scratch/morrill/users/$USER/bin)
#   2. On Morrill:  /scratch/morrill/users/$USER/bin/rclone config
#        n (new remote) -> name: bgs_drive -> storage: drive -> client_id/secret: blank
#        -> scope: 2 (drive.readonly) -> root_folder_id: 1c4wPBOLKLz7ZIIccmnKc8Rfz08r7vJXB
#        -> service_account_file: blank -> advanced: n -> auto config: n
#      It prints a command `rclone authorize "drive" "..."`. Run that on your Mac
#      (brew install rclone, or the zip from rclone.org), sign in in the browser,
#      and paste the token it prints back into the Morrill prompt. Shared drive: n.
#   3. Check:  /scratch/morrill/users/$USER/bin/rclone lsl bgs_drive:
#
# Then:
#   cd /scratch/morrill/users/$USER/biogoship_granite_kegg
#   sbatch --account=191001-364393 cluster/download_reference.sbatch
# Other file names: --export=ALL,FILES="a.mmi b.tsv c.tsv"

set -euo pipefail
BIN="/scratch/morrill/users/${USER}/bin"
RCLONE="${BIN}/rclone"
if [ "${1:-}" = "--install" ] || [ ! -x "${RCLONE}" ]; then
  mkdir -p "${BIN}"
  tmp="$(mktemp -d)"
  curl -fsSL -o "${tmp}/rclone.zip" https://downloads.rclone.org/rclone-current-linux-amd64.zip
  unzip -q -o "${tmp}/rclone.zip" -d "${tmp}"
  cp "${tmp}"/rclone-*-linux-amd64/rclone "${RCLONE}" && chmod +x "${RCLONE}"
  rm -rf "${tmp}"
  "${RCLONE}" version | head -1
  [ "${1:-}" = "--install" ] && exit 0
fi

source "${SLURM_SUBMIT_DIR:-$(pwd)}/cluster/env.sh"
REMOTE="${REMOTE:-bgs_drive:}"
FILES="${FILES:-biogoship_2026_database_v2.mmi biogoship_2026_taxonomy.tsv biogoship_eggnogmapper.tsv}"
DEST="${REPO_ROOT}/data/reference/biogoship_2026"
mkdir -p "${DEST}"

echo "Job started: $(date)"
df -h "${REPO_ROOT}" | tail -1
echo "Files in ${REMOTE} (size, date, name):"
"${RCLONE}" lsl "${REMOTE}" --max-depth 2 | sort -k4 || true

includes=()
for f in ${FILES}; do includes+=(--include "${f}" --include "${f}.gz"); done
"${RCLONE}" copy "${REMOTE}" "${DEST}" "${includes[@]}" --max-depth 2 \
  --transfers 3 --checkers 4 --drive-chunk-size 128M --retries 10 --low-level-retries 20 \
  --stats 10m --stats-one-line -v

echo "Downloaded:"
ls -la "${DEST}"
for f in "${DEST}"/*.tsv "${DEST}"/*.tsv.gz; do
  [ -f "${f}" ] || continue
  echo "== ${f}: first lines"
  if [[ "${f}" == *.gz ]]; then gzip -dc "${f}" | head -3 | cut -c1-400 || true; else head -3 "${f}" | cut -c1-400 || true; fi
done
df -h "${REPO_ROOT}" | tail -1
echo "Job finished: $(date)"
echo "Paste back: this .out file (file list and the first lines of each table)."
