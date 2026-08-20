#!/usr/bin/env bash
set -euo pipefail
CASE_DIR=${1:?'Usage: run_cmc_sweep.sh CASE_DIR [VAMPIRE_EXECUTABLE] [MAX_JOBS]'}
VAMPIRE=${2:-${VAMPIRE_BIN:-vampire-serial}}
MAX_JOBS=${3:-8}
if ! command -v "$VAMPIRE" >/dev/null 2>&1 && [[ ! -x "$VAMPIRE" ]]; then
  echo "Cannot find VAMPIRE executable: $VAMPIRE" >&2
  exit 2
fi
cd "$CASE_DIR"
active=0
for d in angle_*deg; do
  (cd "$d" && "$VAMPIRE" > run.log 2>&1) &
  active=$((active+1))
  if (( active >= MAX_JOBS )); then wait -n; active=$((active-1)); fi
done
wait
echo "All CMC angle jobs completed."
