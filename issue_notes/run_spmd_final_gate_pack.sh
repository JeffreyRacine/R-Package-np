#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/Users/jracine/Development"
NPRMPI_DIR="$ROOT_DIR/np-npRmpi"
NPMASTER_DIR="$ROOT_DIR/np-master"
OUT_DIR="${1:-/tmp/spmd_final_gatepack_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "$OUT_DIR"
echo "OUT_DIR=$OUT_DIR"

# Run a command in the background, logging all output, while printing heartbeat
# lines so long-running steps remain observable in constrained terminals.
run_logged() {
  local step_name="$1"
  local logfile="$2"
  shift 2

  : > "$logfile"
  ("$@" > "$logfile" 2>&1) &
  local pid=$!
  while kill -0 "$pid" 2>/dev/null; do
    printf '[%s] running %s\n' "$step_name" "$(date +%H:%M:%S)"
    sleep 5
  done
  wait "$pid"
  local ec=$?
  printf '[%s] exit=%s\n' "$step_name" "$ec"
  return "$ec"
}

# Preflight snapshots
{
  git -C "$ROOT_DIR" status --short || true
} > "$OUT_DIR/git_status_root.log"
{
  git -C "$NPMASTER_DIR" status --short || true
} > "$OUT_DIR/git_status_np_master.log"
{
  git -C "$NPRMPI_DIR" status --short || true
} > "$OUT_DIR/git_status_np_nprmpi.log"

# Static scans
rg -n --glob 'demo/*.R' '(^|[^#])\btry\s*\(|(^|[^#])\btryCatch\s*\(' "$NPRMPI_DIR" > "$OUT_DIR/demo_mask_scan.log" || true
rg -n 'manual_distributed_call\(' "$NPRMPI_DIR/R" "$NPRMPI_DIR/tests" > "$OUT_DIR/scan_manual_distributed_call.log" || true
rg -n 'nslaves\s*==\s*0|nslaves\s*<=\s*0|npRmpi\.init\(.*nslaves\s*=\s*0' "$NPRMPI_DIR/R" > "$OUT_DIR/scan_nslaves_zero_runtime.log" || true
rg -n 'nslaves\s*==\s*0|nslaves\s*<=\s*0|npRmpi\.init\(.*nslaves\s*=\s*0' "$NPRMPI_DIR/tests" > "$OUT_DIR/scan_nslaves_zero_tests.log" || true
rg -n 'np_reg_cv_ls_stable_ll_glp|stable_ll_glp|O\(n\^2\)' "$NPRMPI_DIR/src" "$NPRMPI_DIR/R" > "$OUT_DIR/scan_on2_helper.log" || true

# Orphan pre-scan
ps aux | rg 'slavedaemon\.R|Rslaves\.sh' | rg -v 'rg ' > "$OUT_DIR/orphan_pre.log" || true

LIB_CAND="$OUT_DIR/lib_candidate"
mkdir -p "$LIB_CAND"
run_logged "install_candidate" "$OUT_DIR/install_candidate.log" \
  R CMD INSTALL "$NPRMPI_DIR" -l "$LIB_CAND"

set +e
pushd "$NPRMPI_DIR" >/dev/null
run_logged "test_spmd_step_contract" "$OUT_DIR/test_spmd_step_contract.log" \
  env R_LIBS="$LIB_CAND" NOT_CRAN=true NP_RMPI_ENABLE_ATTACH_TEST=1 NP_RMPI_ENABLE_PROFILE_TEST=1 NP_RMPI_ENABLE_DENSITY_INID_TEST=1 \
  Rscript --no-save -e 'library(testthat); library(npRmpi); test_file("tests/testthat/test-spmd-step-contract.R")'
EC_STEP=$?

run_logged "test_npsigtest_spmd_contract" "$OUT_DIR/test_npsigtest_spmd_contract.log" \
  env R_LIBS="$LIB_CAND" NOT_CRAN=true NP_RMPI_ENABLE_ATTACH_TEST=1 NP_RMPI_ENABLE_PROFILE_TEST=1 NP_RMPI_ENABLE_DENSITY_INID_TEST=1 \
  Rscript --no-save -e 'library(testthat); library(npRmpi); test_file("tests/testthat/test-npsigtest-spmd-contract.R")'
EC_NPSIG=$?

run_logged "session_smoke" "$OUT_DIR/session_smoke.log" \
  env R_LIBS="$LIB_CAND" \
  Rscript --no-save -e 'suppressPackageStartupMessages(library(npRmpi)); npRmpi.init(nslaves=1, quiet=TRUE); on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE); set.seed(1); x<-runif(80); y<-rnorm(80); bw<-npregbw(y~x, regtype="lc", bwmethod="cv.ls", nmulti=1); fit<-npreg(bws=bw); stopifnot(inherits(fit,"npregression")); cat("SESSION_N1_OK\n")'
EC_SESSION=$?

run_logged "validate_attach" "$OUT_DIR/validate_attach.log" \
  env R_LIBS="$LIB_CAND" FI_TCP_IFACE=en0 \
  mpiexec -n 2 Rscript --no-save "$NPRMPI_DIR/issue_notes/validate_route_attach.R"
EC_ATTACH=$?

run_logged "validate_profile" "$OUT_DIR/validate_profile.log" \
  env R_LIBS="$LIB_CAND" R_PROFILE_USER="$NPRMPI_DIR/inst/Rprofile" FI_TCP_IFACE=en0 \
  mpiexec -n 2 Rscript --no-save "$NPRMPI_DIR/issue_notes/validate_route_profile.R"
EC_PROFILE=$?

run_logged "validate_manual" "$OUT_DIR/validate_manual.log" \
  env R_LIBS="$LIB_CAND" \
  Rscript --no-save "$NPRMPI_DIR/issue_notes/validate_route_manual_broadcast.R"
EC_MANUAL=$?

run_logged "test_session_routing_allflags" "$OUT_DIR/test_session_routing_allflags.log" \
  env R_LIBS="$LIB_CAND" NOT_CRAN=true NP_RMPI_ENABLE_ATTACH_TEST=1 NP_RMPI_ENABLE_PROFILE_TEST=1 NP_RMPI_ENABLE_DENSITY_INID_TEST=1 \
  Rscript --no-save -e 'library(testthat); library(npRmpi); test_file("tests/testthat/test-session-routing-subprocess-contract.R")'
EC_CONTRACT=$?
popd >/dev/null
set -e

LIB_NP="$OUT_DIR/lib_np_master"
mkdir -p "$LIB_NP"
run_logged "install_np_master" "$OUT_DIR/install_np_master.log" \
  R CMD INSTALL "$NPMASTER_DIR" -l "$LIB_NP"
run_logged "np_master_anchor_smoke" "$OUT_DIR/np_master_anchor_smoke.log" \
  env R_LIBS="$LIB_NP" \
  Rscript --no-save /tmp/spmd_canonical_20260304_0001/phase23_finalize_20260303_2330/np_master_anchor_smoke.R

# Deterministic local as-cran harness (note-only target)
run_logged "ascran_clean" "$OUT_DIR/ascran_harness_stdout.log" \
  "$NPRMPI_DIR/issue_notes/run_as_cran_local_clean.sh" "$OUT_DIR/ascran_clean"

# Orphan post-scan
ps aux | rg 'slavedaemon\.R|Rslaves\.sh' | rg -v 'rg ' > "$OUT_DIR/orphan_post.log" || true

# Summaries
rg -n '\[ FAIL [1-9]' "$OUT_DIR/test_spmd_step_contract.log" > "$OUT_DIR/test_spmd_step_contract_failnz.log" || true
rg -n '\[ FAIL [1-9]' "$OUT_DIR/test_npsigtest_spmd_contract.log" > "$OUT_DIR/test_npsigtest_spmd_contract_failnz.log" || true
wc -l "$OUT_DIR/test_spmd_step_contract_failnz.log" > "$OUT_DIR/test_spmd_step_contract_failnz.count"
wc -l "$OUT_DIR/test_npsigtest_spmd_contract_failnz.log" > "$OUT_DIR/test_npsigtest_spmd_contract_failnz.count"

# Gate matrix
rg -n 'SESSION_N1_OK' "$OUT_DIR/session_smoke.log" >/dev/null && TOK_SESSION=0 || TOK_SESSION=1
rg -n 'ATTACH_ROUTE_OK' "$OUT_DIR/validate_attach.log" >/dev/null && TOK_ATTACH=0 || TOK_ATTACH=1
rg -n 'PROFILE_ROUTE_OK' "$OUT_DIR/validate_profile.log" >/dev/null && TOK_PROFILE=0 || TOK_PROFILE=1
rg -n 'MANUAL_BCAST_ROUTE_OK' "$OUT_DIR/validate_manual.log" >/dev/null && TOK_MANUAL=0 || TOK_MANUAL=1
rg -n '\[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 38 \]' "$OUT_DIR/test_session_routing_allflags.log" >/dev/null && TOK_CONTRACT=0 || TOK_CONTRACT=1
rg -n 'NP_MASTER_ANCHOR_OK' "$OUT_DIR/np_master_anchor_smoke.log" >/dev/null && TOK_ANCHOR=0 || TOK_ANCHOR=1
rg -n 'Status: 1 NOTE' "$OUT_DIR/ascran_clean/rcmd_check_as_cran.log" >/dev/null && TOK_ASCRAN=0 || TOK_ASCRAN=1

FAILNZ_STEP=$(awk '{print $1}' "$OUT_DIR/test_spmd_step_contract_failnz.count")
FAILNZ_NPSIG=$(awk '{print $1}' "$OUT_DIR/test_npsigtest_spmd_contract_failnz.count")
ORPH_PRE=$(wc -l < "$OUT_DIR/orphan_pre.log")
ORPH_POST=$(wc -l < "$OUT_DIR/orphan_post.log")
if [ "$ORPH_PRE" -eq 0 ] && [ "$ORPH_POST" -eq 0 ]; then TOK_ORPHAN=0; else TOK_ORPHAN=1; fi

printf 'EC_STEP=%s\nEC_NPSIG=%s\nEC_SESSION=%s\nEC_ATTACH=%s\nEC_PROFILE=%s\nEC_MANUAL=%s\nEC_CONTRACT=%s\nTOK_SESSION=%s\nTOK_ATTACH=%s\nTOK_PROFILE=%s\nTOK_MANUAL=%s\nTOK_CONTRACT=%s\nTOK_ANCHOR=%s\nTOK_ASCRAN=%s\nFAILNZ_STEP=%s\nFAILNZ_NPSIG=%s\nTOK_ORPHAN=%s\n' \
  "$EC_STEP" "$EC_NPSIG" "$EC_SESSION" "$EC_ATTACH" "$EC_PROFILE" "$EC_MANUAL" "$EC_CONTRACT" "$TOK_SESSION" "$TOK_ATTACH" "$TOK_PROFILE" "$TOK_MANUAL" "$TOK_CONTRACT" "$TOK_ANCHOR" "$TOK_ASCRAN" "$FAILNZ_STEP" "$FAILNZ_NPSIG" "$TOK_ORPHAN" \
  > "$OUT_DIR/gate_matrix.log"

if [ "$TOK_SESSION" -eq 0 ] && [ "$TOK_ATTACH" -eq 0 ] && [ "$TOK_PROFILE" -eq 0 ] && [ "$TOK_MANUAL" -eq 0 ] && [ "$TOK_CONTRACT" -eq 0 ] && [ "$TOK_ANCHOR" -eq 0 ] && [ "$TOK_ASCRAN" -eq 0 ] && [ "$FAILNZ_STEP" -eq 0 ] && [ "$FAILNZ_NPSIG" -eq 0 ] && [ "$TOK_ORPHAN" -eq 0 ]; then
  echo FINAL_GATE_PACK_OK > "$OUT_DIR/final_gate_pack.token"
fi

echo "DONE"
