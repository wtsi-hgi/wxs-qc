#!/usr/bin/env bash
set -euo pipefail

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

HAS_FAILURE=0

run_check() {
  local name="$1"
  shift
  local out_file="$TMP_DIR/${name}.log"

  echo "[validate] Running ${name}..."
  if "$@" 2>&1 | tee "$out_file"; then
    echo "[validate] ${name}: PASS"
  else
    echo "[validate] ${name}: FAIL"
    HAS_FAILURE=1
  fi
  echo
}

run_check "check" make check
run_check "typecheck" make typecheck

if [[ "$HAS_FAILURE" -ne 0 ]]; then
  echo "[validate] One or more checks failed or were unavailable."
  exit 1
fi

echo "[validate] Baseline checks passed."
