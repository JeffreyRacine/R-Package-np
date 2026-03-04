#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="/Users/jracine/Development/np-npRmpi"
BASE_TMP="/tmp/spmd_canonical_20260304_0001"
OUT_DIR="${1:-$BASE_TMP/phase36_full_audit_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$OUT_DIR"

VERIFY_DIR="$OUT_DIR/verify"
BUNDLE_DIR="$OUT_DIR/bundle"

"$REPO_DIR/issue_notes/verify_spmd_closure_manifest.sh" "$VERIFY_DIR" > "$OUT_DIR/verify_stdout.log" 2>&1
"$REPO_DIR/issue_notes/export_spmd_closure_bundle.sh" "$BUNDLE_DIR" > "$OUT_DIR/bundle_stdout.log" 2>&1

VERIFY_TOKEN="$(cat "$VERIFY_DIR/closure_manifest_verify.token")"
BUNDLE_TOKEN="$(cat "$BUNDLE_DIR/closure_bundle_export.token")"

{
  echo "verify_dir=$VERIFY_DIR"
  echo "bundle_dir=$BUNDLE_DIR"
  echo "verify_token=$VERIFY_TOKEN"
  echo "bundle_token=$BUNDLE_TOKEN"
} > "$OUT_DIR/full_audit_summary.env"

if [ "$VERIFY_TOKEN" = "CLOSURE_MANIFEST_VERIFY_OK" ] && [ "$BUNDLE_TOKEN" = "CLOSURE_BUNDLE_EXPORT_OK" ]; then
  echo "SPMD_CLOSURE_FULL_AUDIT_OK" > "$OUT_DIR/full_audit.token"
fi

echo "OUT_DIR=$OUT_DIR"
