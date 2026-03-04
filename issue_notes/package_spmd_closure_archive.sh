#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="/Users/jracine/Development/np-npRmpi"
BASE_TMP="/tmp/spmd_canonical_20260304_0001"
OUT_DIR="${1:-$BASE_TMP/phase37_closure_archive_$(date +%Y%m%d_%H%M%S)}"
EXPORT_DIR="$OUT_DIR/export"

mkdir -p "$OUT_DIR"

"$REPO_DIR/issue_notes/export_spmd_closure_bundle.sh" "$EXPORT_DIR" > "$OUT_DIR/export_stdout.log" 2>&1
EXPORT_TOKEN="$(cat "$EXPORT_DIR/closure_bundle_export.token")"
if [ "$EXPORT_TOKEN" != "CLOSURE_BUNDLE_EXPORT_OK" ]; then
  echo "export failed: token=$EXPORT_TOKEN" >&2
  exit 1
fi

ARCHIVE="$OUT_DIR/spmd_closure_bundle.tgz"
tar -czf "$ARCHIVE" -C "$EXPORT_DIR" closure_bundle
ARCHIVE_SHA="$(shasum -a 256 "$ARCHIVE" | awk '{print $1}')"

{
  echo "generated_at=$(date '+%Y-%m-%dT%H:%M:%S%z')"
  echo "repo_head=$(git -C "$REPO_DIR" rev-parse --short HEAD)"
  echo "export_dir=$EXPORT_DIR"
  echo "export_token=$EXPORT_TOKEN"
  echo "archive=$ARCHIVE"
  echo "archive_sha256=$ARCHIVE_SHA"
} > "$OUT_DIR/closure_archive_summary.env"

echo "$ARCHIVE_SHA  $(basename "$ARCHIVE")" > "$OUT_DIR/closure_archive_sha256.txt"
echo "CLOSURE_ARCHIVE_PACKAGE_OK" > "$OUT_DIR/closure_archive_package.token"
echo "OUT_DIR=$OUT_DIR"
