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

if command -v pre-commit >/dev/null 2>&1; then
  run_check "pre-commit-all-files" pre-commit run --all-files
  run_check "mypy" bash -c "scripts/stage_mypy_numbered_scripts.sh && mypy --config-file=pyproject.toml"
else
  echo "[validate] pre-commit is not available; skipping hook checks."
  HAS_FAILURE=1
fi

if [[ "$HAS_FAILURE" -ne 0 ]]; then
  echo "[validate] One or more checks failed or were unavailable."
  exit 1
fi

echo "[validate] Baseline checks passed."
