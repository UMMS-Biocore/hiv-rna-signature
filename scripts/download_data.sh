#!/usr/bin/env bash
# download_data.sh — reconstruct the data tree for a fresh clone.
#
# Fetches:
#   1. Raw 10x matrices from GEO accession GSE239909 (GSE239909_RAW/)
#   2. Derived outputs + DEG inputs + HIV RNA+ table from the Zenodo record
#
# Usage:
#   bash scripts/download_data.sh             # all
#   bash scripts/download_data.sh geo         # GEO only
#   bash scripts/download_data.sh zenodo      # Zenodo only
#
# Environment variables:
#   ZENODO_RECORD_ID   Zenodo record ID (e.g. 1234567). Required for "zenodo".
#   PROJECT_ROOT       Destination root. Defaults to current directory.

set -euo pipefail

MODE="${1:-all}"
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
GEO_ACC="GSE239909"
GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/${GEO_ACC}/suppl/${GEO_ACC}_RAW.tar"

need() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: required tool '$1' not found on PATH" >&2; exit 1; }
}

fetch_geo() {
  need curl
  need tar
  local dest="${PROJECT_ROOT}/${GEO_ACC}_RAW"
  if [[ -d "$dest" && -n "$(ls -A "$dest" 2>/dev/null || true)" ]]; then
    echo "[geo] ${dest} already populated — skipping. Remove to re-download."
    return 0
  fi
  mkdir -p "$dest"
  echo "[geo] Downloading ${GEO_URL}"
  local tmp_tar="${dest}/${GEO_ACC}_RAW.tar"
  curl -fL --retry 5 --retry-delay 5 -o "$tmp_tar" "$GEO_URL"
  echo "[geo] Extracting to ${dest}"
  tar -xf "$tmp_tar" -C "$dest"
  rm -f "$tmp_tar"
  echo "[geo] Done: $(ls "$dest" | wc -l) files"
}

fetch_zenodo() {
  need curl
  need python3
  need tar
  if [[ -z "${ZENODO_RECORD_ID:-}" ]]; then
    echo "ERROR: set ZENODO_RECORD_ID to the published Zenodo record id (e.g. 1234567)" >&2
    echo "       Example: ZENODO_RECORD_ID=1234567 bash scripts/download_data.sh zenodo" >&2
    exit 1
  fi
  local api="https://zenodo.org/api/records/${ZENODO_RECORD_ID}"
  echo "[zenodo] Querying ${api}"
  local listing
  listing=$(curl -fsSL "$api")
  # Fetch every file into $PROJECT_ROOT, verifying MD5s from the Zenodo API.
  python3 - <<'PY' "$listing" "$PROJECT_ROOT"
import json, sys, os, hashlib, subprocess, pathlib
payload = json.loads(sys.argv[1])
root = pathlib.Path(sys.argv[2])
files = payload.get("files", [])
if not files:
    print("[zenodo] No files in record", file=sys.stderr); sys.exit(1)
for f in files:
    name = f["key"]
    url = f["links"]["self"]
    expected = f.get("checksum", "")  # format: "md5:<hex>"
    out = root / name
    out.parent.mkdir(parents=True, exist_ok=True)
    print(f"[zenodo] -> {name}")
    subprocess.check_call(["curl", "-fL", "--retry", "5", "--retry-delay", "5", "-o", str(out), url])
    if expected.startswith("md5:"):
        want = expected.split(":", 1)[1]
        h = hashlib.md5()
        with open(out, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        got = h.hexdigest()
        if got != want:
            print(f"[zenodo] CHECKSUM MISMATCH for {name}: got {got}, want {want}", file=sys.stderr)
            sys.exit(2)
PY
  # Data directories ship as uncompressed .tar archives. Extract them in-place
  # at $PROJECT_ROOT, then remove the archive files.
  shopt -s nullglob
  for t in "${PROJECT_ROOT}"/*.tar; do
    echo "[zenodo] Extracting $(basename "$t")"
    tar -xf "$t" -C "$PROJECT_ROOT"
    rm -f "$t"
  done
  shopt -u nullglob
  echo "[zenodo] Done."
}

case "$MODE" in
  all)    fetch_geo; fetch_zenodo ;;
  geo)    fetch_geo ;;
  zenodo) fetch_zenodo ;;
  *)      echo "Usage: $0 [all|geo|zenodo]" >&2; exit 1 ;;
esac
