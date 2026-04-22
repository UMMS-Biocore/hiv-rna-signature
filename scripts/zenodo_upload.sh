#!/usr/bin/env bash
# zenodo_upload.sh — create a Zenodo deposition and upload this project's
# dataset payload via the REST API. Does NOT publish — leaves the draft so
# you can review it in the web UI and hit "Publish" manually.
#
# Requirements:
#   - curl, python3, jq
#   - Personal access token from https://zenodo.org/account/settings/applications/tokens/new/
#     with scopes: deposit:write, deposit:actions
#
# Usage:
#   ZENODO_TOKEN=xxxxx bash scripts/zenodo_upload.sh [--sandbox] [--new-version]
#
# Flags:
#   --sandbox       target sandbox.zenodo.org for a dry run (strongly recommended first)
#   --new-version   create a new version of an existing deposition (reads .zenodo_state.json)
#
# What it uploads (from repo root):
#   - Top-level files as-is:
#     bc_RNA_HIV.txt, Clayton_DEG_List.xlsx, Article_DEG_list.xlsx,
#     MANIFEST.sha256, README.md, LICENSE, .zenodo.json, CITATION.cff
#   - One uncompressed .tar per data directory (preserves folder structure
#     inside the archive; no gzip since .rds/.mtx.gz/.tsv.gz payloads are
#     already compressed):
#     GSE239909_RAW.tar, Outputs_Code1_20251130.tar, ...
#
# State:
#   On first run it writes .zenodo_state.json with { deposition_id, bucket_url, html_url, sandbox }.
#   That file is gitignored.

set -euo pipefail

BASE="https://zenodo.org"
SANDBOX=0
NEW_VERSION=0
for arg in "$@"; do
  case "$arg" in
    --sandbox)     SANDBOX=1; BASE="https://sandbox.zenodo.org" ;;
    --new-version) NEW_VERSION=1 ;;
    -h|--help)     sed -n '2,30p' "$0"; exit 0 ;;
    *)             echo "Unknown flag: $arg" >&2; exit 1 ;;
  esac
done

# Auto-load token from .env if present (do NOT echo contents).
ROOT_GUESS="$(cd "$(dirname "$0")/.." && pwd)"
if [[ -f "$ROOT_GUESS/.env" ]]; then
  # shellcheck disable=SC1090,SC1091
  set -a; . "$ROOT_GUESS/.env"; set +a
fi

: "${ZENODO_TOKEN:?Set ZENODO_TOKEN (in .env or environment) to a Zenodo personal access token}"
command -v curl >/dev/null    || { echo "curl missing"    >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 missing" >&2; exit 1; }
command -v jq >/dev/null      || { echo "jq missing — install with 'brew install jq'" >&2; exit 1; }

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

STATE_FILE=".zenodo_state.json"
API="$BASE/api"
AUTH=(-H "Authorization: Bearer $ZENODO_TOKEN")

# ---- 1. Create deposition (or new version) -----------------------------------
if [[ $NEW_VERSION -eq 1 ]]; then
  [[ -f "$STATE_FILE" ]] || { echo "--new-version requires existing $STATE_FILE" >&2; exit 1; }
  PARENT_ID=$(jq -r '.deposition_id' "$STATE_FILE")
  echo "[zenodo] Creating new version of deposition $PARENT_ID"
  NEW=$(curl -fsS -X POST "${AUTH[@]}" "$API/deposit/depositions/$PARENT_ID/actions/newversion")
  DRAFT_URL=$(echo "$NEW" | jq -r '.links.latest_draft')
  DEPOSITION=$(curl -fsS "${AUTH[@]}" "$DRAFT_URL")
else
  echo "[zenodo] Creating new deposition on $BASE"
  DEPOSITION=$(curl -fsS -X POST "${AUTH[@]}" \
    -H "Content-Type: application/json" \
    -d '{"metadata":{}}' \
    "$API/deposit/depositions")
fi

DEP_ID=$(echo "$DEPOSITION" | jq -r '.id')
BUCKET=$(echo "$DEPOSITION" | jq -r '.links.bucket')
HTML=$(echo   "$DEPOSITION" | jq -r '.links.html')

if [[ -z "$DEP_ID" || "$DEP_ID" == "null" ]]; then
  echo "[zenodo] Failed to create deposition:" >&2
  echo "$DEPOSITION" >&2
  exit 1
fi

echo "[zenodo] deposition_id = $DEP_ID"
echo "[zenodo] draft URL     = $HTML"

python3 - "$DEP_ID" "$BUCKET" "$HTML" "$SANDBOX" "$STATE_FILE" <<'PY'
import json, sys
dep_id, bucket, html, sandbox, state = sys.argv[1:]
json.dump({
    "deposition_id": int(dep_id),
    "bucket_url": bucket,
    "html_url": html,
    "sandbox": bool(int(sandbox)),
}, open(state, "w"), indent=2)
PY

# ---- 2. Push metadata from .zenodo.json --------------------------------------
if [[ -f ".zenodo.json" ]]; then
  echo "[zenodo] Pushing metadata from .zenodo.json"
  META=$(python3 -c 'import json,sys; print(json.dumps({"metadata": json.load(open(".zenodo.json"))}))')
  curl -fsS -X PUT "${AUTH[@]}" \
    -H "Content-Type: application/json" \
    -d "$META" \
    "$API/deposit/depositions/$DEP_ID" > /dev/null
fi

# ---- 3. Stage tar archives (uncompressed, one per folder) --------------------
# Each line in $TMP_LIST is: <local_path>\t<zenodo_key>
TMP_LIST=$(mktemp)
TAR_DIR=$(mktemp -d)
trap 'rm -f "$TMP_LIST"; rm -rf "$TAR_DIR"' EXIT

# 3a. Top-level files, shipped as-is.
for f in bc_RNA_HIV.txt Clayton_DEG_List.xlsx Article_DEG_list.xlsx \
         MANIFEST.sha256 README.md LICENSE .zenodo.json CITATION.cff; do
  [[ -f "$f" ]] && printf '%s\t%s\n' "$f" "$f" >> "$TMP_LIST"
done

# 3b. One uncompressed .tar per directory, staged in $TAR_DIR.
tar_dirs=()
[[ -d "GSE239909_RAW" ]] && tar_dirs+=("GSE239909_RAW")
for d in Outputs_Code1_*/ Outputs_Code2_*/ Outputs_Code3_*/ Outputs_Code4_*/ \
         Outputs_Code5_*/ Outputs_Code6_*/ Outputs_Code7_*/; do
  d="${d%/}"
  [[ -d "$d" ]] && tar_dirs+=("$d")
done

for d in "${tar_dirs[@]}"; do
  tarball="$TAR_DIR/${d}.tar"
  echo "[zenodo] Packing $d/ → ${d}.tar"
  tar -cf "$tarball" "$d"
  actual=$(stat -f%z "$tarball" 2>/dev/null || stat -c%s "$tarball")
  printf '[zenodo]   %s.tar  = %d bytes (%.2f GB)\n' "$d" "$actual" "$(python3 -c "print($actual/1e9)")"
  printf '%s\t%s\n' "$tarball" "${d}.tar" >> "$TMP_LIST"
done

COUNT=$(wc -l < "$TMP_LIST" | tr -d ' ')
if [[ "$COUNT" -eq 0 ]]; then
  echo "[zenodo] No files found to upload — aborting." >&2
  exit 1
fi
echo "[zenodo] Will upload $COUNT files"

# ---- 4. Upload via bucket API (PUT, streaming) -------------------------------
upload_one() {
  local path="$1"
  local key="$2"
  local encoded
  encoded=$(python3 -c 'import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1], safe=""))' "$key")
  local size
  size=$(stat -f%z "$path" 2>/dev/null || stat -c%s "$path")
  printf '[zenodo]   %s  (%s bytes)\n' "$key" "$size"
  curl -fS --retry 5 --retry-delay 10 \
    -X PUT "${AUTH[@]}" \
    --upload-file "$path" \
    "${BUCKET}/${encoded}" > /dev/null
}

I=0
while IFS=$'\t' read -r path key; do
  I=$((I+1))
  printf '[zenodo] [%d/%d] ' "$I" "$COUNT"
  upload_one "$path" "$key"
done < "$TMP_LIST"

echo
echo "[zenodo] ============================================================"
echo "[zenodo] Upload complete. Draft is NOT yet published."
echo "[zenodo] Review at: $HTML"
echo "[zenodo] When ready, click 'Publish' in the web UI, or run:"
echo "[zenodo]   curl -X POST -H \"Authorization: Bearer \$ZENODO_TOKEN\" $API/deposit/depositions/$DEP_ID/actions/publish"
echo "[zenodo] ============================================================"
