#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./bsub.sh <conf> <final_root> <macro> <jobs> [total_events]

Splits a large Geant4 run into multiple LSF array jobs and merges the ROOT files.
- <conf>         : Path to the .conf file (e.g. param/conf/e90_woHD.conf)
- <final_root>   : Target merged ROOT file path (e.g. rootfile/SigmaNCusp.root)
- <macro>        : Macro with /run/beamOn (e.g. test.mac)
- <jobs>         : Maximum number of array jobs to submit
- [total_events] : Optional override if /run/beamOn is missing in the macro

Environment variables:
- QUEUE=<queue>      : LSF queue name
- LOG_DIR=<dir>      : Log directory (default: logs/<basename>)
- AUTO_MERGE=1       : Wait for completion and run hadd automatically
EOF
  exit 1
}

[[ $# -lt 4 ]] && usage

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONF=$1
FINAL_ROOT=$2
MACRO=$3
REQ_JOBS=$4
TOTAL_EVENTS=${5:-}

[[ -f "$CONF" ]] || { echo "Config not found: $CONF" >&2; exit 1; }
[[ -f "$MACRO" ]] || { echo "Macro not found: $MACRO" >&2; exit 1; }
[[ "$REQ_JOBS" =~ ^[0-9]+$ ]] && ((REQ_JOBS > 0)) || { echo "Jobs must be a positive integer." >&2; exit 1; }

if [[ -z "$TOTAL_EVENTS" ]]; then
  TOTAL_EVENTS=$(awk '/^[[:space:]]*\/run\/beamOn[[:space:]]+/ {print $2; exit}' "$MACRO" || true)
fi
[[ -n "$TOTAL_EVENTS" && "$TOTAL_EVENTS" =~ ^[0-9]+$ && "$TOTAL_EVENTS" -gt 0 ]] || {
  echo "Total events not found. Provide it as arg5 or ensure the macro has '/run/beamOn <N>'." >&2
  exit 1
}

G4S2S_EXE="${script_dir}/bin/G4S2S"
[[ -x "$G4S2S_EXE" ]] || { echo "Executable not found: $G4S2S_EXE" >&2; exit 1; }

CHUNK=$(( (TOTAL_EVENTS + REQ_JOBS - 1) / REQ_JOBS ))
JOBS=$(( (TOTAL_EVENTS + CHUNK - 1) / CHUNK ))
if (( JOBS < REQ_JOBS )); then
  echo "Requested $REQ_JOBS jobs, but only $JOBS contain events (total: $TOTAL_EVENTS)." >&2
fi

OUT_DIR=$(dirname "$FINAL_ROOT")
mkdir -p "$OUT_DIR"
BASE=$(basename "$FINAL_ROOT")
BASE=${BASE%.root}

MACRO_BASE=$(basename "$MACRO")
MACRO_STEM=${MACRO_BASE%.*}
MACRO_DIR="$(pwd)/job_macros_${BASE}"
mkdir -p "$MACRO_DIR"

EVENTS_LEFT=$TOTAL_EVENTS
for i in $(seq 1 "$JOBS"); do
  EVENTS=$CHUNK
  if (( EVENTS_LEFT < CHUNK )); then
    EVENTS=$EVENTS_LEFT
  fi
  EVENTS_LEFT=$(( EVENTS_LEFT - EVENTS ))
  MACRO_PATH="${MACRO_DIR}/${MACRO_STEM}_part${i}.mac"
  awk -v n="$EVENTS" '
    BEGIN {replaced=0}
    /^[[:space:]]*\/run\/beamOn[[:space:]]+/ {print "/run/beamOn " n; replaced=1; next}
    {print}
    END {if(!replaced) print "/run/beamOn " n}
  ' "$MACRO" > "$MACRO_PATH"
done

LOG_DIR="${LOG_DIR:-logs/${BASE}}"
mkdir -p "$LOG_DIR"

JOB_NAME="${BASE}_split"
queue_name="${QUEUE:-s}"
QUEUE_OPT=(-q "$queue_name")

bsub "${QUEUE_OPT[@]}" -J "${JOB_NAME}[1-${JOBS}]" \
  -oo "${LOG_DIR}/${JOB_NAME}.%I.out" \
  -eo "${LOG_DIR}/${JOB_NAME}.%I.err" \
  bash -lc "set -euo pipefail; idx=\${LSB_JOBINDEX}; out=\"${OUT_DIR}/${BASE}_part\${idx}.root\"; macro=\"${MACRO_DIR}/${MACRO_STEM}_part\${idx}.mac\"; cd \"${script_dir}\"; \"${G4S2S_EXE}\" \"${CONF}\" \"\${out}\" \"\${macro}\""

echo "Submitted ${JOBS} jobs as array ${JOB_NAME}[1-${JOBS}] (total events: ${TOTAL_EVENTS}, chunk: ${CHUNK})."
echo "Partial ROOT files: ${OUT_DIR}/${BASE}_part*.root"
echo "To merge after completion:"
echo "  bwait -w \"ended(${JOB_NAME})\""
echo "  hadd -f \"${FINAL_ROOT}\" ${OUT_DIR}/${BASE}_part*.root"

if [[ "${AUTO_MERGE:-0}" == "1" ]]; then
  command -v bwait >/dev/null 2>&1 || { echo "bwait not found; skipping AUTO_MERGE." >&2; exit 0; }
  command -v hadd  >/dev/null 2>&1 || { echo "hadd not found; skipping AUTO_MERGE." >&2; exit 0; }
  echo "Waiting for ${JOB_NAME} to finish..."
  bwait -w "ended(${JOB_NAME})"
  shopt -s nullglob
  PARTS=( ${OUT_DIR}/${BASE}_part*.root )
  shopt -u nullglob
  (( ${#PARTS[@]} > 0 )) || { echo "No partial ROOT files found; merge skipped." >&2; exit 1; }
  hadd -f "$FINAL_ROOT" "${PARTS[@]}"
  echo "Merged into $FINAL_ROOT"
  rm -f "${PARTS[@]}"
  echo "Removed partial files: ${PARTS[*]}"
  rm -rf "$MACRO_DIR"
  echo "Cleaned macro directory: ${MACRO_DIR}"
fi

# Schedule cleanup of macro directory after jobs finish (covers AUTO_MERGE=0 case).
if [[ "${AUTO_MERGE:-0}" != "1" ]]; then
  bsub "${QUEUE_OPT[@]}" -w "ended(${JOB_NAME})" -J "${JOB_NAME}_cleanup" \
    -oo "${LOG_DIR}/${JOB_NAME}_cleanup.out" \
    -eo "${LOG_DIR}/${JOB_NAME}_cleanup.err" \
    bash -lc "rm -rf \"${MACRO_DIR}\""
  echo "Cleanup job scheduled to remove ${MACRO_DIR} after array completion."
fi
