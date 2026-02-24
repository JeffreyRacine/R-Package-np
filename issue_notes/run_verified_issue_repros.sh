#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_LOG="${1:-/tmp/np_issue_notes_repros_${STAMP}.log}"
INSTALL_LOG="${RUN_LOG%.log}.install.log"
TMP_LIB="$(mktemp -d /tmp/np_issue_notes_lib.XXXXXX)"
TIMEOUT_SEC="${NP_ISSUE_REPRO_TIMEOUT_SEC:-480}"

run_with_timeout() {
  local timeout_sec="$1"
  shift
  if command -v timeout >/dev/null 2>&1; then
    timeout "${timeout_sec}" "$@"
  elif command -v gtimeout >/dev/null 2>&1; then
    gtimeout "${timeout_sec}" "$@"
  else
    perl -e '$t=shift; alarm $t; exec @ARGV' "${timeout_sec}" "$@"
  fi
}

cleanup() {
  rm -rf "${TMP_LIB}"
}
trap cleanup EXIT

echo "[info] installing np from ${ROOT_DIR}" | tee "${INSTALL_LOG}"
R CMD INSTALL --preclean -l "${TMP_LIB}" "${ROOT_DIR}" >>"${INSTALL_LOG}" 2>&1

echo "[info] running issue-note verified repros" | tee "${RUN_LOG}"
echo "[info] issue-note repro timeout=${TIMEOUT_SEC}s" | tee -a "${RUN_LOG}"
R_LIBS="${TMP_LIB}${R_LIBS:+:${R_LIBS}}" \
  run_with_timeout "${TIMEOUT_SEC}" \
  Rscript "${ROOT_DIR}/issue_notes/verified_issue_repros.R" >>"${RUN_LOG}" 2>&1

if rg -n "FAIL:" "${RUN_LOG}" >/dev/null 2>&1; then
  echo "[error] one or more issue-note repros failed; inspect ${RUN_LOG}" >&2
  rg -n "FAIL:" "${RUN_LOG}" >&2 || true
  exit 1
fi

echo "[ok] all verified issue-note repros passed"
echo "[ok] run log: ${RUN_LOG}"
echo "[ok] install log: ${INSTALL_LOG}"
