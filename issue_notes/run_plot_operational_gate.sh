#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${1:-/tmp/nprmpi_plot_operational_gate_${STAMP}}"
TEST_LOG="${OUT_DIR}/plot_contract_suite.log"
ROUTE_LOG="${OUT_DIR}/route_runner.log"
SUMMARY="${OUT_DIR}/summary.txt"
PRE_SMOKE_LOG="${OUT_DIR}/session_attach_sanity.log"
PRIMARY_IFACE="${NP_RMPI_IFACE_PRIMARY:-en0}"
FALLBACK_IFACE="${NP_RMPI_IFACE_FALLBACK:-lo0}"
GATE_TIMEOUT_TEST_SEC="${NP_RMPI_PLOT_GATE_TEST_TIMEOUT_SEC:-300}"
GATE_TIMEOUT_ROUTE_SEC="${NP_RMPI_PLOT_GATE_ROUTE_TIMEOUT_SEC:-420}"

mkdir -p "${OUT_DIR}"

run_with_timeout() {
  local timeout_sec="$1"
  shift
  if [ "${timeout_sec}" -le 0 ]; then
    "$@"
    return
  fi
  if command -v timeout >/dev/null 2>&1; then
    timeout "${timeout_sec}" "$@"
  elif command -v gtimeout >/dev/null 2>&1; then
    gtimeout "${timeout_sec}" "$@"
  else
    perl -e '$t=shift; alarm $t; exec @ARGV' "${timeout_sec}" "$@"
  fi
}

kill_stray_mpi_slaves() {
  local pids
  pids="$(pgrep -f 'slavedaemon\.R|Rslaves\.sh' || true)"
  if [ -n "${pids}" ]; then
    echo "[warn] cleaning stray slave daemons: ${pids}" | tee -a "${PRE_SMOKE_LOG}"
    printf '%s\n' "${pids}" | xargs kill >/dev/null 2>&1 || true
  fi
}

cleanup() {
  kill_stray_mpi_slaves || true
}
trap cleanup EXIT

echo "[info] output dir: ${OUT_DIR}" | tee "${PRE_SMOKE_LOG}"
echo "[info] repo: ${ROOT_DIR}" | tee -a "${PRE_SMOKE_LOG}"
echo "[info] interface primary/fallback: ${PRIMARY_IFACE}/${FALLBACK_IFACE}" | tee -a "${PRE_SMOKE_LOG}"

# Tiny session/attach sanity smoke (plot errors none) for quick route health.
cat > "${OUT_DIR}/session_attach_sanity.R" <<'RS'
suppressPackageStartupMessages(library(npRmpi))

npRmpi.init(nslaves = 1, quiet = TRUE)
on.exit(try(npRmpi.quit(), silent = TRUE), add = TRUE)
set.seed(1001)
n <- 30L
x <- runif(n)
y <- sin(2 * pi * x) + rnorm(n, sd = 0.05)
bw <- npregbw(y ~ x, data = data.frame(y = y, x = x), bws = 0.2, bandwidth.compute = FALSE)
png(tempfile(fileext = ".png")); on.exit(dev.off(), add = TRUE)
out <- suppressWarnings(plot(
  bw,
  xdat = data.frame(x = x),
  ydat = y,
  plot.behavior = "data",
  perspective = FALSE,
  plot.errors.method = "none"
))
stopifnot(is.list(out), length(out) > 0)
cat("SESSION_PLOT_SANITY_OK\n")
RS

run_with_timeout 120 Rscript "${OUT_DIR}/session_attach_sanity.R" >>"${PRE_SMOKE_LOG}" 2>&1

cat > "${OUT_DIR}/attach_plot_sanity.R" <<'RS'
suppressPackageStartupMessages(library(npRmpi))
npRmpi.init(mode = "attach", quiet = TRUE, autodispatch = TRUE)
if (mpi.comm.rank(1L) == 0L) {
  set.seed(1002)
  n <- 30L
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.05)
  bw <- npregbw(y ~ x, data = data.frame(y = y, x = x), bws = 0.2, bandwidth.compute = FALSE)
  png(tempfile(fileext = ".png")); on.exit(dev.off(), add = TRUE)
  out <- suppressWarnings(plot(
    bw,
    xdat = data.frame(x = x),
    ydat = y,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "none"
  ))
  stopifnot(is.list(out), length(out) > 0)
  cat("ATTACH_PLOT_SANITY_OK\n")
  npRmpi.quit(mode = "attach")
}
RS

set +e
mpiexec -n 2 env FI_PROVIDER=tcp FI_TCP_IFACE="${PRIMARY_IFACE}" \
  Rscript --vanilla "${OUT_DIR}/attach_plot_sanity.R" >>"${PRE_SMOKE_LOG}" 2>&1
rc=$?
if [ "${rc}" -ne 0 ]; then
  mpiexec -n 2 env FI_PROVIDER=tcp FI_TCP_IFACE="${FALLBACK_IFACE}" \
    Rscript --vanilla "${OUT_DIR}/attach_plot_sanity.R" >>"${PRE_SMOKE_LOG}" 2>&1
  rc=$?
fi
set -e
if [ "${rc}" -ne 0 ]; then
  echo "[error] attach sanity smoke failed" | tee -a "${PRE_SMOKE_LOG}"
  exit "${rc}"
fi

kill_stray_mpi_slaves

echo "[info] running plot contract suite" | tee "${TEST_LOG}"
cd "${ROOT_DIR}"
run_with_timeout "${GATE_TIMEOUT_TEST_SEC}" \
  Rscript -e "devtools::test(filter='plot-helper-option-contract|plot-bootstrap-arg-contract|plot-asymptotic-failfast-contract|plot-coef-option-contract|plot-conditional-gradients-bootstrap-contract|plot-mpi-only-bootstrap-contract|plot-autodispatch|plot-guardrails-contract|session-routing-subprocess-contract', reporter='summary')" \
  >>"${TEST_LOG}" 2>&1

echo "[info] running attach/profile/profile-plot route matrix" | tee "${ROUTE_LOG}"
run_with_timeout "${GATE_TIMEOUT_ROUTE_SEC}" \
  bash "${ROOT_DIR}/issue_notes/run_attach_profile_validations.sh" \
  >>"${ROUTE_LOG}" 2>&1

ROUTE_DIR="$(rg -n '\[ok\] out dir:' "${ROUTE_LOG}" | tail -n 1 | sed -E 's/.*out dir: //')"
ps aux | rg -n "slavedaemon\\.R|Rslaves\\.sh" > "${OUT_DIR}/orphan_scan.log" || true
SESSION_SANITY="$(rg -c 'SESSION_PLOT_SANITY_OK' "${PRE_SMOKE_LOG}" || echo 0)"
ATTACH_SANITY="$(rg -c 'ATTACH_PLOT_SANITY_OK' "${PRE_SMOKE_LOG}" || echo 0)"
PLOT_SUITE_DONE="$(rg -c '══ DONE' "${TEST_LOG}" || echo 0)"
PLOT_SUITE_FAILED="$(rg -c '══ Failed' "${TEST_LOG}" || echo 0)"
ROUTE_RUNNER_PASS="$(rg -c '\[ok\] route smoke checks passed' "${ROUTE_LOG}" || echo 0)"

{
  echo "out_dir=${OUT_DIR}"
  echo "session_sanity=${SESSION_SANITY}"
  echo "attach_sanity=${ATTACH_SANITY}"
  echo "plot_suite_done=${PLOT_SUITE_DONE}"
  echo "plot_suite_failed=${PLOT_SUITE_FAILED}"
  echo "route_runner_pass=${ROUTE_RUNNER_PASS}"
  echo "route_runner_out_dir=${ROUTE_DIR:-NA}"
  if [[ -s "${OUT_DIR}/orphan_scan.log" ]]; then echo "orphan_matches=1"; else echo "orphan_matches=0"; fi
} > "${SUMMARY}"

echo "[ok] plot operational gate complete"
cat "${SUMMARY}"
