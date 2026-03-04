#!/usr/bin/env bash
set -euo pipefail

REPO="${1:-/Users/jracine/Development/np-npRmpi}"
OUT_DIR="${2:-/tmp/build_stage_hardening_phaseA_$(date +%Y%m%d_%H%M%S)}"
JOBS="${JOBS:-npreg,npplot,npcondensitybw}"
TIMEOUT_SEC="${TIMEOUT_SEC:-240}"
NSLAVES="${NSLAVES:-1}"
R_BIN="${R_BIN:-Rscript}"
R_LIBS_VALUE="${R_LIBS_VALUE:-}"
KILL_ORPHANS_ON_FAIL="${KILL_ORPHANS_ON_FAIL:-0}"
PRE_CLEAN_ORPHANS="${PRE_CLEAN_ORPHANS:-1}"
WILD_MASTER_LOCAL_GUARD="${WILD_MASTER_LOCAL_GUARD:-}"
BOOTSTRAP_PHASE_TRACE="${BOOTSTRAP_PHASE_TRACE:-0}"
TRANSPORT_TRACE="${TRANSPORT_TRACE:-0}"

mkdir -p "${OUT_DIR}/jobs"
SUMMARY="${OUT_DIR}/summary.tsv"
TOKENS="${OUT_DIR}/summary_tokens.log"

if [ ! -d "${REPO}/issue_notes" ]; then
  echo "ERROR: REPO does not look like npRmpi repo root: ${REPO}" >&2
  exit 2
fi

echo -e "job\tstatus\texit_code\tlog" > "${SUMMARY}"
echo "OUT_DIR=${OUT_DIR}" > "${TOKENS}"
echo "REPO=${REPO}"
echo "OUT_DIR=${OUT_DIR}"
echo "JOBS=${JOBS}"
echo "TIMEOUT_SEC=${TIMEOUT_SEC}"
if [ -n "${WILD_MASTER_LOCAL_GUARD}" ]; then
  echo "WILD_MASTER_LOCAL_GUARD=${WILD_MASTER_LOCAL_GUARD}"
fi
echo "BOOTSTRAP_PHASE_TRACE=${BOOTSTRAP_PHASE_TRACE}"
echo "TRANSPORT_TRACE=${TRANSPORT_TRACE}"

write_job_script() {
  job="$1"
  phase_file="$2"
  script="${OUT_DIR}/jobs/${job}.R"

  case "${job}" in
    npreg)
      cat > "${script}" <<EOF
options(np.messages=FALSE)
options(npRmpi.bootstrap.phase.file="${phase_file}")
stage <- function(label) { cat(sprintf("STAGE %s\\n", label)); flush.console() }
stage("library")
suppressPackageStartupMessages(library(npRmpi))
stage("init")
npRmpi.init(nslaves=${NSLAVES}, quiet=TRUE)
on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)
set.seed(42)
n <- 120
x <- runif(n)
y <- x + rnorm(n, sd=0.3*sd(x))
bw <- npregbw(y ~ x, regtype="ll", bwmethod="cv.ls", nmulti=1)
stage("fit")
fit <- npreg(bws=bw)
stopifnot(inherits(fit, "npregression"))
cat("JOB_OK npreg\\n")
EOF
      ;;
    npplot)
      cat > "${script}" <<EOF
options(np.messages=FALSE)
options(npRmpi.bootstrap.phase.file="${phase_file}")
EOF
      if [ -n "${WILD_MASTER_LOCAL_GUARD}" ]; then
        cat >> "${script}" <<EOF
options(npRmpi.plot.wild.master_local.guard = ${WILD_MASTER_LOCAL_GUARD})
EOF
      fi
      cat >> "${script}" <<EOF
stage <- function(label) { cat(sprintf("STAGE %s\\n", label)); flush.console() }
stage("library")
suppressPackageStartupMessages(library(npRmpi))
stage("init")
npRmpi.init(nslaves=${NSLAVES}, quiet=TRUE)
on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)
set.seed(123)
n <- 100
x1 <- runif(n); x2 <- runif(n); x3 <- runif(n); x4 <- rbinom(n, 2, .3)
y <- 1 + x1 + x2 + x3 + x4 + rnorm(n)
X <- data.frame(x1, x2, x3, ordered(x4))
bw <- npregbw(xdat=X, ydat=y, regtype="ll", bwmethod="cv.aic", nmulti=1)
stage("plot")
out <- plot(
  bw,
  perspective=FALSE,
  plot.behavior="data",
  plot.errors.method="bootstrap",
  plot.errors.center="bias-corrected",
  plot.errors.type="simultaneous",
  plot.errors.boot.num=29
)
stopifnot(is.list(out))
cat("JOB_OK npplot\\n")
EOF
      ;;
    npplot_wild_stress)
      cat > "${script}" <<EOF
options(np.messages=FALSE)
options(npRmpi.bootstrap.phase.file="${phase_file}")
EOF
      if [ -n "${WILD_MASTER_LOCAL_GUARD}" ]; then
        cat >> "${script}" <<EOF
options(npRmpi.plot.wild.master_local.guard = ${WILD_MASTER_LOCAL_GUARD})
EOF
      fi
      cat >> "${script}" <<EOF
stage <- function(label) { cat(sprintf("STAGE %s\\n", label)); flush.console() }
stage("library")
suppressPackageStartupMessages(library(npRmpi))
stage("init")
npRmpi.init(nslaves=${NSLAVES}, quiet=TRUE)
on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)
set.seed(123)
n <- 1000
x1 <- runif(n); x2 <- runif(n); x3 <- runif(n); x4 <- rbinom(n, 2, .3)
y <- 1 + x1 + x2 + x3 + x4 + rnorm(n)
X <- data.frame(x1, x2, x3, ordered(x4))
bw <- npregbw(xdat=X, ydat=y, regtype="ll", bwmethod="cv.aic", nmulti=1)
stage("plot")
out <- plot(
  bw,
  perspective=FALSE,
  plot.behavior="data",
  plot.errors.method="bootstrap",
  plot.errors.boot.method="wild",
  plot.errors.center="bias-corrected",
  plot.errors.type="simultaneous",
  plot.errors.boot.num=799
)
stopifnot(is.list(out))
cat("JOB_OK npplot_wild_stress\\n")
EOF
      ;;
    npcondensitybw)
      cat > "${script}" <<EOF
options(np.messages=FALSE)
options(npRmpi.bootstrap.phase.file="${phase_file}")
stage <- function(label) { cat(sprintf("STAGE %s\\n", label)); flush.console() }
stage("library")
suppressPackageStartupMessages(library(npRmpi))
stage("init")
npRmpi.init(nslaves=${NSLAVES}, quiet=TRUE)
on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)
set.seed(99)
n <- 100
x <- runif(n)
y <- x + rnorm(n, sd=0.25)
bw <- npcdensbw(xdat=data.frame(x=x), ydat=y, bwmethod="cv.ls", nmulti=1)
stage("fit")
fit <- npcdens(bws=bw)
stopifnot(inherits(fit, "condensity"))
cat("JOB_OK npcondensitybw\\n")
EOF
      ;;
    *)
      echo "ERROR: unknown job '${job}'" >&2
      return 1
      ;;
  esac
}

run_with_timeout() {
  timeout="$1"
  shift
  perl -MPOSIX -e '
    $t = shift @ARGV;
    POSIX::setsid();
    $SIG{ALRM} = sub {
      kill "TERM", -$$;
      exit 124;
    };
    alarm($t);
    exec @ARGV;
    exit 127;
  ' "${timeout}" "$@"
}

orphans_present() {
  if ps aux | rg -q "slavedaemon\\.R|Rslaves\\.sh"; then
    return 0
  fi
  return 1
}

kill_orphans() {
  pkill -f "slavedaemon\\.R|Rslaves\\.sh" || true
}

failures=0
IFS=',' read -r -a JOB_ARR <<< "${JOBS}"
for job in "${JOB_ARR[@]}"; do
  job="$(echo "${job}" | sed 's/^ *//;s/ *$//')"
  [ -z "${job}" ] && continue
  echo "RUNNING_JOB=${job}"
  if [ "${PRE_CLEAN_ORPHANS}" = "1" ]; then
    kill_orphans
  fi
  phase_file=""
  if [ "${BOOTSTRAP_PHASE_TRACE}" = "1" ]; then
    phase_file="${OUT_DIR}/${job}.phase.tsv"
    : > "${phase_file}"
  fi
  transport_file=""
  if [ "${TRANSPORT_TRACE}" = "1" ]; then
    transport_file="${OUT_DIR}/${job}.transport.tsv"
    : > "${transport_file}"
  fi
  write_job_script "${job}" "${phase_file}"
  log="${OUT_DIR}/${job}.log"
  rc=0

  set +e
  if [ -n "${R_LIBS_VALUE}" ] && [ -n "${transport_file}" ]; then
    R_LIBS="${R_LIBS_VALUE}" NP_RMPI_TRANSPORT_TRACE_FILE="${transport_file}" run_with_timeout "${TIMEOUT_SEC}" "${R_BIN}" --no-save "${OUT_DIR}/jobs/${job}.R" > "${log}" 2>&1
    rc=$?
  elif [ -n "${R_LIBS_VALUE}" ]; then
    R_LIBS="${R_LIBS_VALUE}" run_with_timeout "${TIMEOUT_SEC}" "${R_BIN}" --no-save "${OUT_DIR}/jobs/${job}.R" > "${log}" 2>&1
    rc=$?
  elif [ -n "${transport_file}" ]; then
    NP_RMPI_TRANSPORT_TRACE_FILE="${transport_file}" run_with_timeout "${TIMEOUT_SEC}" "${R_BIN}" --no-save "${OUT_DIR}/jobs/${job}.R" > "${log}" 2>&1
    rc=$?
  else
    run_with_timeout "${TIMEOUT_SEC}" "${R_BIN}" --no-save "${OUT_DIR}/jobs/${job}.R" > "${log}" 2>&1
    rc=$?
  fi
  set -e

  if [ "${rc}" -eq 0 ] && rg -q -F "JOB_OK ${job}" "${log}"; then
    status="PASS"
  else
    status="FAIL"
    failures=$((failures + 1))
    if [ "${KILL_ORPHANS_ON_FAIL}" = "1" ]; then
      kill_orphans
    fi
  fi

  if orphans_present; then
    orphan_flag="1"
  else
    orphan_flag="0"
  fi

  echo -e "${job}\t${status}\t${rc}\t${log}" >> "${SUMMARY}"
  echo "JOB_DONE=${job} STATUS=${status} EXIT=${rc} ORPHANS=${orphan_flag}"
  {
    echo "JOB_${job}_STATUS=${status}"
    echo "JOB_${job}_EXIT=${rc}"
    echo "JOB_${job}_ORPHANS=${orphan_flag}"
    if [ -n "${phase_file}" ]; then
      echo "JOB_${job}_PHASE_FILE=${phase_file}"
      if [ -s "${phase_file}" ]; then
        echo "JOB_${job}_PHASE_LINES=$(wc -l < "${phase_file}" | tr -d ' ')"
      else
        echo "JOB_${job}_PHASE_LINES=0"
      fi
    fi
    if [ -n "${transport_file}" ]; then
      echo "JOB_${job}_TRANSPORT_FILE=${transport_file}"
      if [ -s "${transport_file}" ]; then
        echo "JOB_${job}_TRANSPORT_LINES=$(wc -l < "${transport_file}" | tr -d ' ')"
      else
        echo "JOB_${job}_TRANSPORT_LINES=0"
      fi
    fi
  } >> "${TOKENS}"
done

echo "FAILURES=${failures}" >> "${TOKENS}"
cat "${TOKENS}"
echo "SUMMARY=${SUMMARY}"

if [ "${failures}" -ne 0 ]; then
  exit 1
fi

echo "MPI_SUBPROCESS_HARNESS_OK"
