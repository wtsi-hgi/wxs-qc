#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(git rev-parse --show-toplevel)"
STAGE_DIR="$ROOT_DIR/.mypy_stage"

NUMBERED_DIRS=(
  "0-resource_preparation"
  "1-import_data"
  "2-sample_qc"
  "3-variant_qc"
  "4-genotype_qc"
)

sanitize_package_name() {
  local name="$1"
  name="${name//-/_}"
  if [[ "$name" =~ ^([0-9]+)_(.*)$ ]]; then
    printf "stage_%s_%s" "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
  else
    printf "%s" "$name"
  fi
}

sanitize_module_name() {
  local name="$1"

  if [[ "$name" == "__init__.py" ]]; then
    printf "%s" "$name"
    return
  fi

  name="${name%.py}"
  name="${name//-/_}"
  if [[ "$name" =~ ^[0-9] ]]; then
    name="step_${name}"
  fi
  printf "%s.py" "$name"
}

link_file() {
  local source="$1"
  local target="$2"

  if [[ -e "$target" && ! -L "$target" ]]; then
    echo "Refusing to replace non-symlink mypy stage file: $target" >&2
    exit 1
  fi

  ln -sfn "$source" "$target"
}

mkdir -p "$STAGE_DIR"
find "$STAGE_DIR" -type l -delete

for source_dir in "${NUMBERED_DIRS[@]}"; do
  package_dir="$STAGE_DIR/$(sanitize_package_name "$source_dir")"
  mkdir -p "$package_dir"

  while IFS= read -r -d "" source_file; do
    source_name="$(basename "$source_file")"
    target_name="$(sanitize_module_name "$source_name")"
    link_file "$source_file" "$package_dir/$target_name"
  done < <(find "$ROOT_DIR/$source_dir" -maxdepth 1 -type f -name "*.py" -print0)
done
