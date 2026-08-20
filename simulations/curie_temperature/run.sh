#!/usr/bin/env bash
set -euo pipefail
VAMPIRE="${1:-${VAMPIRE_BIN:-vampire-parallel}}"
if ! command -v "$VAMPIRE" >/dev/null 2>&1 && [[ ! -x "$VAMPIRE" ]]; then
  echo "Cannot find VAMPIRE executable: $VAMPIRE" >&2
  echo "Pass it as argument, set VAMPIRE_BIN, or add vampire-parallel to PATH." >&2
  exit 2
fi
"$VAMPIRE" > run.log 2>&1
