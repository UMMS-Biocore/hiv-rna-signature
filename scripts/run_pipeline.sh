#!/usr/bin/env bash
# run_pipeline.sh — orchestrator for the 7-step R workflow.
#
# Contract:
#   - PROJECT_ROOT must be set (absolute path to the project directory on disk,
#     or /mnt/project inside the Docker image). Scripts refuse to run without it.
#   - WORKFLOW_STEP picks which step(s) to run: all | 1..7 | <script filename>.
#
# Usage:
#   PROJECT_ROOT=/abs/path bash scripts/run_pipeline.sh         # runs all 7
#   PROJECT_ROOT=/abs/path bash scripts/run_pipeline.sh 4       # runs step 4 only
#   PROJECT_ROOT=/abs/path bash scripts/run_pipeline.sh 1_UMAP_filter_ordered_PG.R

set -euo pipefail

# Default only inside the Docker image. On host, PROJECT_ROOT must be set explicitly.
PROJECT_ROOT="${PROJECT_ROOT:-/mnt/project}"
WORKFLOW_STEP="${1:-all}"

# Repo root = parent of this script's dir. No hardcoded absolute paths.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

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
  echo "ERROR: PROJECT_ROOT does not exist: $PROJECT_ROOT" >&2
  echo "Set PROJECT_ROOT to the absolute path of your project directory." >&2
  exit 1
fi

export PROJECT_ROOT

run_script() {
  local script_name="$1"
  local src="${REPO_ROOT}/${script_name}"
  if [[ ! -f "$src" ]]; then
    echo "ERROR: missing script: $src" >&2
    exit 1
  fi
  echo "[pipeline] Running ${script_name}  (PROJECT_ROOT=${PROJECT_ROOT})"
  Rscript "$src"
}

case "$WORKFLOW_STEP" in
  all)
    for s in "${SCRIPTS[@]}"; do run_script "$s"; done
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
