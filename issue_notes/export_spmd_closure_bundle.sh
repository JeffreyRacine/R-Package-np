#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="/Users/jracine/Development/np-npRmpi"
BASE_TMP="/tmp/spmd_canonical_20260304_0001"
OUT_DIR="${1:-$BASE_TMP/phase35_closure_bundle_$(date +%Y%m%d_%H%M%S)}"
BUNDLE_DIR="$OUT_DIR/closure_bundle"

mkdir -p "$BUNDLE_DIR"

copy_if_exists() {
  local src="$1"
  local dst="$2"
  if [ -e "$src" ]; then
    cp -R "$src" "$dst"
  fi
}

# Repo-resident closeout records
copy_if_exists "$REPO_DIR/REMEDIATION_PLAN_SPMD_CANONICAL_2026-03-03.md" "$BUNDLE_DIR/"
copy_if_exists "$REPO_DIR/SPMD_REMEDIATION_FINAL_HANDOFF_2026-03-04.md" "$BUNDLE_DIR/"
copy_if_exists "$REPO_DIR/SPMD_REMEDIATION_CLOSURE_CERTIFICATE_2026-03-04.md" "$BUNDLE_DIR/"
copy_if_exists "$REPO_DIR/SPMD_REMEDIATION_CLOSURE_MANIFEST_2026-03-04.json" "$BUNDLE_DIR/"
copy_if_exists "$REPO_DIR/issue_notes/verify_spmd_closure_manifest.sh" "$BUNDLE_DIR/"
copy_if_exists "$REPO_DIR/issue_notes/run_spmd_final_gate_pack.sh" "$BUNDLE_DIR/"

# Tmp evidence index (pointer-rich, not immutable storage)
copy_if_exists "$BASE_TMP/FINAL_EVIDENCE_INDEX_20260304.md" "$BUNDLE_DIR/"

# Bundle metadata
{
  echo "generated_at=$(date '+%Y-%m-%dT%H:%M:%S%z')"
  echo "repo_head=$(git -C "$REPO_DIR" rev-parse --short HEAD)"
  echo "tag_closure=$(git -C "$REPO_DIR" rev-parse --short spmd-remediation-complete-2026-03-04^{} 2>/dev/null || true)"
  echo "tag_post_certificate=$(git -C "$REPO_DIR" rev-parse --short spmd-remediation-post-cert-closeout-2026-03-04^{} 2>/dev/null || true)"
  echo "phase31=/tmp/spmd_canonical_20260304_0001/phase31_manual_gatepack_20260304_005852"
  echo "phase32=/tmp/spmd_canonical_20260304_0001/phase32_final_gatepack_20260304_010940"
  echo "phase33=/tmp/spmd_canonical_20260304_0001/phase33_archival_gatepack_20260304_011716"
  echo "phase34=/tmp/spmd_canonical_20260304_0001/phase34_manifest_verify_20260304_013133"
} > "$BUNDLE_DIR/BUNDLE_METADATA.env"

# Checksums for all files in bundle
(
  cd "$BUNDLE_DIR"
  find . -type f -maxdepth 2 | LC_ALL=C sort | while read -r f; do
    shasum -a 256 "$f"
  done
) > "$OUT_DIR/closure_bundle_sha256.txt"

echo "CLOSURE_BUNDLE_EXPORT_OK" > "$OUT_DIR/closure_bundle_export.token"
echo "OUT_DIR=$OUT_DIR"
