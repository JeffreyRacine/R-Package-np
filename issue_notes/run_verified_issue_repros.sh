#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_LOG="${1:-/tmp/nprmpi_issue_notes_repros_${STAMP}.log}"
INSTALL_LOG="${RUN_LOG%.log}.install.log"
TMP_LIB="$(mktemp -d /tmp/nprmpi_issue_notes_lib.XXXXXX)"

cleanup() {
  rm -rf "${TMP_LIB}"
}
trap cleanup EXIT

echo "[info] installing npRmpi from ${ROOT_DIR}" | tee "${INSTALL_LOG}"
R CMD INSTALL --preclean -l "${TMP_LIB}" "${ROOT_DIR}" >>"${INSTALL_LOG}" 2>&1

echo "[info] running issue-note verified repros" | tee "${RUN_LOG}"
R_LIBS="${TMP_LIB}${R_LIBS:+:${R_LIBS}}" \
  Rscript "${ROOT_DIR}/issue_notes/verified_issue_repros.R" >>"${RUN_LOG}" 2>&1

if rg -n "FAIL:" "${RUN_LOG}" >/dev/null 2>&1; then
  echo "[error] one or more issue-note repros failed; inspect ${RUN_LOG}" >&2
  rg -n "FAIL:" "${RUN_LOG}" >&2 || true
  exit 1
fi

echo "[ok] all verified issue-note repros passed"
echo "[ok] run log: ${RUN_LOG}"
echo "[ok] install log: ${INSTALL_LOG}"
