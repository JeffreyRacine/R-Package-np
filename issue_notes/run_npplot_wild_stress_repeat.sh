#!/usr/bin/env bash
set -euo pipefail

REPO="${1:-/Users/jracine/Development/np-npRmpi}"
OUT_DIR="${2:-/tmp/npplot_wild_repeat_$(date +%Y%m%d_%H%M%S)}"
REPEATS="${REPEATS:-5}"
TIMEOUT_SEC="${TIMEOUT_SEC:-180}"
N="${N:-1000}"
B="${B:-199}"
NSLAVES="${NSLAVES:-1}"
BOOT_METHOD="${BOOT_METHOD:-wild}"
WILD_GUARD="${WILD_MASTER_LOCAL_GUARD:-TRUE}"
R_BIN="${R_BIN:-Rscript}"
R_LIBS_VALUE="${R_LIBS_VALUE:-}"
TIMEOUT_BIN="${TIMEOUT_BIN:-/opt/local/bin/timeout}"
PRE_CLEAN_ORPHANS="${PRE_CLEAN_ORPHANS:-1}"

mkdir -p "${OUT_DIR}"

if [ ! -d "${REPO}/issue_notes" ]; then
  echo "ERROR: REPO does not look like npRmpi repo root: ${REPO}" >&2
  exit 2
fi

if ! command -v "${TIMEOUT_BIN}" >/dev/null 2>&1; then
  if command -v timeout >/dev/null 2>&1; then
    TIMEOUT_BIN="timeout"
  else
    echo "ERROR: timeout binary not found (${TIMEOUT_BIN})" >&2
    exit 2
  fi
fi

if [ "${PRE_CLEAN_ORPHANS}" = "1" ]; then
  pkill -f 'slavedaemon\.R|Rslaves\.sh' || true
fi

SCRIPT="${OUT_DIR}/npplot_wild_stress_once.R"
cat > "${SCRIPT}" <<RSCRIPT
options(np.messages=FALSE)
stage <- function(s) { cat(sprintf("STAGE %s\\n", s)); flush.console() }
options(npRmpi.plot.wild.master_local.guard=${WILD_GUARD})
stage("library")
suppressPackageStartupMessages(library(npRmpi))
stage("init")
npRmpi.init(nslaves=${NSLAVES}, quiet=TRUE)
on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)
set.seed(123)
n <- ${N}
x1 <- runif(n); x2 <- runif(n); x3 <- runif(n); x4 <- rbinom(n, 2, .3)
y <- 1 + x1 + x2 + x3 + x4 + rnorm(n)
X <- data.frame(x1, x2, x3, ordered(x4))
stage("bw")
bw <- npregbw(xdat=X, ydat=y, regtype="ll", bwmethod="cv.aic", nmulti=1)
stage("plot")
out <- plot(
  bw,
  perspective=FALSE,
  plot.behavior="data",
  plot.errors.method="bootstrap",
  plot.errors.boot.method="${BOOT_METHOD}",
  plot.errors.center="bias-corrected",
  plot.errors.type="simultaneous",
  plot.errors.boot.num=${B}
)
stage("done")
stopifnot(is.list(out))
cat("RUN_OK\\n")
RSCRIPT

SUMMARY="${OUT_DIR}/summary.tsv"
echo -e "run\trc\trun_ok\tlast_stage\tlog" > "${SUMMARY}"

for i in $(seq 1 "${REPEATS}"); do
  LOG="${OUT_DIR}/run_${i}.log"
  rc=0
  set +e
  if [ -n "${R_LIBS_VALUE}" ]; then
    R_LIBS="${R_LIBS_VALUE}" "${TIMEOUT_BIN}" "${TIMEOUT_SEC}" "${R_BIN}" --no-save "${SCRIPT}" > "${LOG}" 2>&1
    rc=$?
  else
    "${TIMEOUT_BIN}" "${TIMEOUT_SEC}" "${R_BIN}" --no-save "${SCRIPT}" > "${LOG}" 2>&1
    rc=$?
  fi
  set -e

  run_ok=0
  if rg -q '^RUN_OK$' "${LOG}"; then
    run_ok=1
  fi

  last_stage=""
  if rg -q '^STAGE ' "${LOG}"; then
    last_stage="$(rg '^STAGE ' "${LOG}" | tail -n 1 | sed 's/^STAGE //')"
  fi

  echo -e "${i}\t${rc}\t${run_ok}\t${last_stage}\t${LOG}" >> "${SUMMARY}"
  pkill -f 'slavedaemon\.R|Rslaves\.sh' || true
done

pass_count="$(awk -F'\t' 'NR>1 && $3==1 {c++} END{print c+0}' "${SUMMARY}")"
fail_count="$((REPEATS - pass_count))"

echo "OUT_DIR=${OUT_DIR}"
echo "SUMMARY=${SUMMARY}"
echo "PASS=${pass_count} FAIL=${fail_count} REPEATS=${REPEATS}"
cat "${SUMMARY}"
