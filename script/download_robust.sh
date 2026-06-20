#!/usr/bin/env bash
###############################################################################
# Robust, resumable, collapse-proof downloader for large sequencing data.
#
# WHY: big single transfers fail on unstable links (SSH drop, timeout, stall).
# This downloads PER FILE with: resume (-c), multi-connection, auto-retry,
# stall-detection, and md5 verification. Re-running it just resumes — it never
# re-downloads a file that is already complete & correct. Safe to run 100x.
#
# USAGE (on your cloud server):
#   1) Run inside tmux so it survives disconnects:
#        tmux new -s dl
#   2a) From an ENA/SRA study accession (auto-fetches the file list + md5s):
#        bash download_robust.sh ena PRJNA574438 ./fastq
#   2b) From your own list of URLs (one per line), e.g. GEO supplementary:
#        bash download_robust.sh urls my_urls.txt ./data
#   3) Detach with Ctrl-b then d.  Reattach later: tmux attach -t dl
#
# TIP for servers in regions with slow NCBI access: ENA (EBI, Europe) mirrors
# all SRA data and gives direct resumable FASTQ over HTTPS — prefer it.
#
# Needs: aria2c (preferred) or wget; curl. Install on Ubuntu:
#   apt-get update && apt-get install -y aria2 curl
###############################################################################
set -uo pipefail

MODE="${1:-}"; SRC="${2:-}"; OUT="${3:-./download}"
mkdir -p "$OUT"
MAX_TRIES=20            # per-file attempts before giving up (re-run to continue)

have() { command -v "$1" >/dev/null 2>&1; }
md5_of() { if have md5sum; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }

# Download one URL to $OUT, resuming + retrying until md5 matches (md5 optional).
fetch() {
  local url="$1" want_md5="${2:-}" fname
  fname="$OUT/$(basename "${url%%\?*}")"
  for a in $(seq 1 "$MAX_TRIES"); do
    if [[ -f "$fname" && -n "$want_md5" ]]; then
      [[ "$(md5_of "$fname")" == "$want_md5" ]] && { echo "[OK] $fname"; return 0; }
    elif [[ -f "$fname" && -z "$want_md5" && -s "$fname" ]]; then
      # no md5 to check; assume complete if aria2 left no .aria2 control file
      [[ ! -f "${fname}.aria2" ]] && { echo "[OK?] $fname (no md5 to verify)"; return 0; }
    fi
    echo "[try $a/$MAX_TRIES] $url"
    if have aria2c; then
      aria2c -c -x8 -s8 -k1M --retry-wait=10 --max-tries=3 \
             --lowest-speed-limit=10K --summary-interval=0 --console-log-level=warn \
             -d "$OUT" -o "$(basename "$fname")" "$url" && continue
    else
      wget -c --tries=3 --timeout=60 --waitretry=10 --read-timeout=120 \
           -O "$fname" "$url" && continue
    fi
    echo "  (attempt failed; will retry)"; sleep 8
  done
  echo "[FAIL] $fname after $MAX_TRIES tries — re-run the script to keep trying." ; return 1
}

case "$MODE" in
  ena)
    ACC="$SRC"
    echo "Fetching ENA file report for $ACC ..."
    REPORT="$OUT/_ena_report.tsv"
    curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv" -o "$REPORT" \
      || { echo "Could not fetch ENA report. Check accession/network."; exit 1; }
    n=$(($(wc -l < "$REPORT") - 1)); echo "Runs found: $n"
    fails=0
    # columns: run_accession  fastq_ftp  fastq_md5  fastq_bytes  (URLs/md5 are ';'-separated)
    tail -n +2 "$REPORT" | while IFS=$'\t' read -r run ftp md5 bytes; do
      [[ -z "$ftp" ]] && { echo "[skip] $run has no fastq_ftp"; continue; }
      IFS=';' read -ra U <<< "$ftp"; IFS=';' read -ra M <<< "$md5"
      for i in "${!U[@]}"; do
        fetch "https://${U[$i]}" "${M[$i]:-}" || fails=$((fails+1))
      done
    done
    echo "Done. Any [FAIL] lines above mean: just re-run this exact command to resume."
    ;;
  urls)
    [[ -f "$SRC" ]] || { echo "URL list file not found: $SRC"; exit 1; }
    while IFS= read -r url; do
      [[ -z "$url" || "$url" =~ ^# ]] && continue
      fetch "$url" ""
    done < "$SRC"
    ;;
  *)
    echo "Usage:"
    echo "  bash download_robust.sh ena  <ACCESSION> [outdir]     # e.g. ena PRJNA574438 ./fastq"
    echo "  bash download_robust.sh urls <urls.txt>   [outdir]     # one URL per line"
    exit 1
    ;;
esac
