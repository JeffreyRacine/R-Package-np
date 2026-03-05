#!/usr/bin/env bash
set -euo pipefail

STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${1:-/tmp/nprmpi_plot_helper_semantics_abort_${STAMP}}"
REPEATS="${NP_RMPI_SEMANTICS_REPEATS:-3}"
TIMEOUT_SEC="${NP_RMPI_SEMANTICS_TIMEOUT_SEC:-180}"

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
    printf '%s\n' "${pids}" | xargs kill >/dev/null 2>&1 || true
  fi
}

cat > "${OUT_DIR}/run_case.R" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
case <- args[[1]]
reps <- as.integer(args[[2]])
if (is.na(reps) || reps < 1L) reps <- 1L

suppressPackageStartupMessages(library(npRmpi))
npRmpi.init(nslaves = 1, quiet = TRUE)
on.exit(try(npRmpi.quit(), silent = TRUE), add = TRUE)
options(npRmpi.autodispatch = TRUE, np.messages = FALSE)

quiet_capture <- function(expr) {
  out <- NULL
  sinkfile <- tempfile()
  on.exit(unlink(sinkfile), add = TRUE)
  invisible(capture.output(out <- eval.parent(substitute(expr)), file = sinkfile))
  out
}

run_bootstrap_plot <- function(bw, xdat, ydat, zdat = NULL, boot_num = 9L) {
  quiet_capture(
    suppressWarnings(plot(
      bw,
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = boot_num
    ))
  )
}

for (i in seq_len(reps)) {
  if (case == "r") {
    set.seed(9401 + i)
    n <- 50L
    x <- runif(n)
    y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
    xdat <- data.frame(x = x)
    bw <- npregbw(
      y ~ x,
      data = data.frame(y = y, x = x),
      bws = 0.2,
      bandwidth.compute = FALSE,
      regtype = "lp",
      basis = "glp",
      degree = 2L,
      bernstein.basis = TRUE,
      bwtype = "fixed",
      ckertype = "gaussian"
    )
    out <- run_bootstrap_plot(bw, xdat = xdat, ydat = y, boot_num = 9L)
    cat(sprintf("CASE_OK case=%s rep=%d nout=%d\n", case, i, length(out)))
  } else if (case == "si") {
    set.seed(9501 + i)
    n <- 60L
    xdat <- data.frame(x1 = runif(n), x2 = runif(n))
    y <- xdat$x1^2 + 0.5 * xdat$x2 + rnorm(n, sd = 0.1)
    bw <- npindexbw(
      xdat = xdat,
      ydat = y,
      bws = c(0.2, 0.2, 0.25),
      bandwidth.compute = FALSE,
      regtype = "lp",
      basis = "glp",
      degree = 2L,
      bernstein.basis = TRUE,
      bwtype = "fixed",
      ckertype = "gaussian"
    )
    out <- run_bootstrap_plot(bw, xdat = xdat, ydat = y, boot_num = 9L)
    cat(sprintf("CASE_OK case=%s rep=%d nout=%d\n", case, i, length(out)))
  } else if (case == "sc") {
    set.seed(9601 + i)
    n <- 55L
    xdat <- data.frame(x = runif(n))
    zdat <- data.frame(z = runif(n))
    y <- sin(2 * pi * xdat$x) * zdat$z + rnorm(n, sd = 0.1)
    bw <- npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      bws = 0.2,
      bandwidth.compute = FALSE,
      regtype = "lp",
      basis = "glp",
      degree = 2L,
      bernstein.basis = TRUE,
      bwtype = "fixed",
      ckertype = "gaussian"
    )
    out <- run_bootstrap_plot(bw, xdat = xdat, ydat = y, zdat = zdat, boot_num = 9L)
    cat(sprintf("CASE_OK case=%s rep=%d nout=%d\n", case, i, length(out)))
  } else {
    stop(sprintf("unsupported case '%s'", case))
  }
}
RS

{
  echo "out_dir=${OUT_DIR}"
  echo "repeats=${REPEATS}"
  echo "timeout_sec=${TIMEOUT_SEC}"
} > "${OUT_DIR}/manifest.txt"

echo "case,exit_code,status" > "${OUT_DIR}/results.csv"

for case in r si sc; do
  set +e
  run_with_timeout "${TIMEOUT_SEC}" \
    Rscript "${OUT_DIR}/run_case.R" "${case}" "${REPEATS}" \
    > "${OUT_DIR}/${case}.log" 2>&1
  ec=$?
  set -e
  if [ "${ec}" -eq 0 ]; then
    status="OK"
  elif rg -q "unsupported for smooth coefficient bootstrap in npRmpi canonical SPMD mode" "${OUT_DIR}/${case}.log"; then
    status="FAIL_FAST"
  elif [ "${ec}" -eq 124 ] || [ "${ec}" -eq 137 ]; then
    status="TIMEOUT"
  elif [ "${ec}" -eq 134 ] || [ "${ec}" -eq 139 ] || [ "${ec}" -eq 6 ] || [ "${ec}" -eq 10 ] || [ "${ec}" -eq 11 ]; then
    status="ABORT"
  else
    status="FAIL"
  fi
  echo "${case},${ec},${status}" >> "${OUT_DIR}/results.csv"
  kill_stray_mpi_slaves
done

{
  echo "[ok] semantics abort repro sweep complete"
  cat "${OUT_DIR}/manifest.txt"
  cat "${OUT_DIR}/results.csv"
} | tee "${OUT_DIR}/summary.txt"
