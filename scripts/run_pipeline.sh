#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-/data}"
WORKFLOW_STEP="${1:-all}"

SCRIPTS=(
  "1_UMAP_filter_ordered_PG.R"
  "2_Annotation_macro_ordered_PG.R"
  "3_RNA_positive_cells_annotation_ordered_PG.R"
  "4_DEGs_Scorring_ordered_PG.R"
  "5_UMAP_Scorring_Plots_ordered_PG.R"
  "6_Comparison_profiles_ordered_PG.R"
  "7_RNAplus_comparison_ordered_PG.R"
)

if [[ ! -d "$PROJECT_ROOT" ]]; then
  echo "PROJECT_ROOT does not exist: $PROJECT_ROOT" >&2
  exit 1
fi

run_script() {
  local script_name="$1"
  local src="/workspace/${script_name}"
  local tmp="/tmp/${script_name}"

  if [[ ! -f "$src" ]]; then
    echo "Missing script: $src" >&2
    exit 1
  fi

  sed "s|{Your_Folder}|${PROJECT_ROOT}|g" "$src" > "$tmp"
  echo "Running ${script_name} with PROJECT_ROOT=${PROJECT_ROOT}"
  Rscript "$tmp"
}

case "$WORKFLOW_STEP" in
  all)
    for s in "${SCRIPTS[@]}"; do
      run_script "$s"
    done
    ;;
  1|2|3|4|5|6|7)
    run_script "${SCRIPTS[$((WORKFLOW_STEP-1))]}"
    ;;
  *.R)
    run_script "$WORKFLOW_STEP"
    ;;
  *)
    echo "Invalid step: $WORKFLOW_STEP" >&2
    echo "Use one of: all, 1..7, or an exact script filename." >&2
    exit 1
    ;;
esac
