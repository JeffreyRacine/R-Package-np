#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/Users/jracine/Development/np-npRmpi"
OUT_DIR="${1:-/tmp/nprmpi_modernization_micro_scan_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "$OUT_DIR"
echo "OUT_DIR=$OUT_DIR"

run_scan() {
  local name="$1"
  local pattern="$2"
  shift 2
  rg -n "$pattern" "$@" > "$OUT_DIR/${name}.log" || true
  wc -l < "$OUT_DIR/${name}.log" > "$OUT_DIR/${name}.count"
}

# Wickham-aligned static scans (R layer + selected C/runtime guardrails)
run_scan "scan_eval" 'eval\(' "$ROOT_DIR/R"
run_scan "scan_parse" 'parse\(' "$ROOT_DIR/R"
run_scan "scan_sapply" '\bsapply\(' "$ROOT_DIR/R"
run_scan "scan_range_one_length" '1:length\(' "$ROOT_DIR/R"
run_scan "scan_runtime_library_require" '\blibrary\(|\brequire\(' "$ROOT_DIR/R"
run_scan "scan_dotC" '\.C\(' "$ROOT_DIR/R"
run_scan "scan_nslaves_zero_runtime" 'nslaves\s*==\s*0|nslaves\s*<=\s*0|npRmpi\.init\(.*nslaves\s*=\s*0' "$ROOT_DIR/R"
run_scan "scan_nslaves_zero_tests" 'nslaves\s*==\s*0|nslaves\s*<=\s*0|npRmpi\.init\(.*nslaves\s*=\s*0' "$ROOT_DIR/tests"
run_scan "scan_demo_masking" '(^|[^#])\btry\s*\(|(^|[^#])\btryCatch\s*\(' "$ROOT_DIR/demo"
run_scan "scan_manual_distributed_call" 'manual_distributed_call\(' "$ROOT_DIR/R" "$ROOT_DIR/tests"
run_scan "scan_stable_helper_symbol" 'np_reg_cv_ls_stable_ll_glp|stable_ll_glp' "$ROOT_DIR/src" "$ROOT_DIR/R"

# Current git status snapshots for audit context
git -C "$ROOT_DIR" status --short > "$OUT_DIR/git_status_np_nprmpi.log" || true

{
  echo "MODERNIZATION_MICRO_SCAN_SUMMARY"
  echo "repo=$ROOT_DIR"
  echo "out_dir=$OUT_DIR"
  for c in "$OUT_DIR"/*.count; do
    b="$(basename "$c" .count)"
    v="$(cat "$c")"
    printf '%s=%s\n' "$b" "$v"
  done | sort
} > "$OUT_DIR/summary.log"

cat "$OUT_DIR/summary.log"
echo "MODERNIZATION_MICRO_SCAN_OK"

