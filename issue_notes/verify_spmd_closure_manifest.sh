#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="/Users/jracine/Development/np-npRmpi"
MANIFEST="$REPO_DIR/SPMD_REMEDIATION_CLOSURE_MANIFEST_2026-03-04.json"
OUT_DIR="${1:-/tmp/spmd_canonical_20260304_0001/phase34_manifest_verify_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "$OUT_DIR"

REPORT="$OUT_DIR/closure_manifest_verify.log"
KV="$OUT_DIR/manifest_kv.tsv"
: >"$REPORT"

python - "$MANIFEST" >"$KV" <<'PY'
import json
import sys

manifest_path = sys.argv[1]
m = json.load(open(manifest_path, "r", encoding="utf-8"))

def put(k, v):
    print(f"{k}\t{v}")

put("manifest_head", m["head"])
put("tag_closure_name", m["tags"]["closure"]["name"])
put("tag_closure_target", m["tags"]["closure"]["target"])
put("tag_post_name", m["tags"]["post_certificate"]["name"])
put("tag_post_target", m["tags"]["post_certificate"]["target"])
final_tag = m["tags"].get("final_closure", {})
put("tag_final_name", final_tag.get("name", ""))
put("tag_final_target", final_tag.get("target", ""))

phases = []
for phase in sorted(m["artifacts"].keys()):
    a = m["artifacts"][phase]
    if not phase.startswith("phase"):
        continue
    if not isinstance(a, dict):
        continue
    if not {"path", "token", "gate_matrix_sha256", "ascran_log_sha256"} <= set(a):
        continue
    phases.append(phase)
    put(f"{phase}_path", a["path"])
    put(f"{phase}_token", a["token"])
    put(f"{phase}_gate_sha", a["gate_matrix_sha256"])
    put(f"{phase}_ascran_rel", a.get("ascran_log_relpath", "ascran_clean/rcmd_check_as_cran.log"))
    put(f"{phase}_ascran_sha", a["ascran_log_sha256"])

put("phase_list", ",".join(phases))
PY

getv() {
  local key="$1"
  awk -F'\t' -v k="$key" '$1==k {print $2; exit}' "$KV"
}

record_ok() {
  printf "OK: %s\n" "$1" | tee -a "$REPORT"
}

record_fail() {
  printf "FAIL: %s\n" "$1" | tee -a "$REPORT"
  return 1
}

check_eq() {
  local lhs="$1"
  local rhs="$2"
  local desc="$3"
  if [ "$lhs" = "$rhs" ]; then
    record_ok "$desc"
  else
    record_fail "$desc (expected '$rhs', got '$lhs')"
  fi
}

check_file() {
  local p="$1"
  local desc="$2"
  if [ -f "$p" ]; then
    record_ok "$desc"
  else
    record_fail "$desc (missing: $p)"
  fi
}

check_dir() {
  local p="$1"
  local desc="$2"
  if [ -d "$p" ]; then
    record_ok "$desc"
  else
    record_fail "$desc (missing: $p)"
  fi
}

manifest_head="$(getv manifest_head)"
current_head="$(git -C "$REPO_DIR" rev-parse --short HEAD)"
if git -C "$REPO_DIR" merge-base --is-ancestor "$manifest_head" HEAD; then
  record_ok "manifest head $manifest_head is ancestor of current head $current_head"
else
  record_fail "manifest head $manifest_head is not ancestor of current head $current_head"
fi

tag_closure_name="$(getv tag_closure_name)"
tag_closure_target="$(getv tag_closure_target)"
tag_post_name="$(getv tag_post_name)"
tag_post_target="$(getv tag_post_target)"
tag_final_name="$(getv tag_final_name)"
tag_final_target="$(getv tag_final_target)"

actual_closure_target="$(git -C "$REPO_DIR" rev-parse --short "${tag_closure_name}^{}")"
actual_post_target="$(git -C "$REPO_DIR" rev-parse --short "${tag_post_name}^{}")"

check_eq "$actual_closure_target" "$tag_closure_target" "closure tag target matches manifest"
check_eq "$actual_post_target" "$tag_post_target" "post-certificate tag target matches manifest"
if [ -n "$tag_final_name" ] && [ -n "$tag_final_target" ]; then
  actual_final_target="$(git -C "$REPO_DIR" rev-parse --short "${tag_final_name}^{}")"
  check_eq "$actual_final_target" "$tag_final_target" "final-closure tag target matches manifest"
fi

phase_list="$(getv phase_list)"
if [ -z "$phase_list" ]; then
  record_fail "manifest phase list is empty"
fi

IFS=',' read -r -a phase_arr <<< "$phase_list"
for phase in "${phase_arr[@]}"; do
  phase_path="$(getv "${phase}_path")"
  expected_token="$(getv "${phase}_token")"
  expected_gate_sha="$(getv "${phase}_gate_sha")"
  expected_ascran_rel="$(getv "${phase}_ascran_rel")"
  expected_ascran_sha="$(getv "${phase}_ascran_sha")"

  check_dir "$phase_path" "$phase artifact directory exists"
  check_file "$phase_path/final_gate_pack.token" "$phase token file exists"
  check_file "$phase_path/gate_matrix.log" "$phase gate matrix exists"
  check_file "$phase_path/$expected_ascran_rel" "$phase as-cran log exists"

  actual_token="$(cat "$phase_path/final_gate_pack.token")"
  check_eq "$actual_token" "$expected_token" "$phase token matches manifest"

  actual_gate_sha="$(shasum -a 256 "$phase_path/gate_matrix.log" | awk '{print $1}')"
  actual_ascran_sha="$(shasum -a 256 "$phase_path/$expected_ascran_rel" | awk '{print $1}')"

  check_eq "$actual_gate_sha" "$expected_gate_sha" "$phase gate matrix hash matches manifest"
  check_eq "$actual_ascran_sha" "$expected_ascran_sha" "$phase as-cran hash matches manifest"

  rg -n '^TOK_.*=0$' "$phase_path/gate_matrix.log" >/dev/null
  rg -n '^FAILNZ_.*=[[:space:]]*0$' "$phase_path/gate_matrix.log" >/dev/null
  record_ok "$phase gate matrix zero-token checks passed"
done

if ps -Ao command | rg 'slavedaemon\.R|Rslaves\.sh' | rg -v rg >/dev/null; then
  record_fail "orphan worker scan clean"
else
  record_ok "orphan worker scan clean"
fi

echo "CLOSURE_MANIFEST_VERIFY_OK" >"$OUT_DIR/closure_manifest_verify.token"
record_ok "verification token written to $OUT_DIR/closure_manifest_verify.token"

echo "OUT_DIR=$OUT_DIR"
