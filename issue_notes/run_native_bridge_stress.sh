#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_LOG="${1:-/tmp/np_native_bridge_stress_${STAMP}.log}"
INSTALL_LOG="${RUN_LOG%.log}.install.log"
TMP_LIB="$(mktemp -d /tmp/np_native_bridge_lib.XXXXXX)"

cleanup() {
  rm -rf "${TMP_LIB}"
}
trap cleanup EXIT

echo "[info] installing np from ${ROOT_DIR}" | tee "${INSTALL_LOG}"
R CMD INSTALL --preclean -l "${TMP_LIB}" "${ROOT_DIR}" >>"${INSTALL_LOG}" 2>&1

echo "[info] running native bridge stress harness" | tee "${RUN_LOG}"
R_LIBS="${TMP_LIB}${R_LIBS:+:${R_LIBS}}" \
  Rscript "${ROOT_DIR}/issue_notes/native_bridge_stress.R" >>"${RUN_LOG}" 2>&1

if rg -n "NATIVE_BRIDGE_STRESS_OK" "${RUN_LOG}" >/dev/null 2>&1; then
  echo "[ok] native bridge stress passed"
  echo "[ok] run log: ${RUN_LOG}"
  echo "[ok] install log: ${INSTALL_LOG}"
  exit 0
fi

echo "[error] native bridge stress did not report success token" >&2
tail -n 80 "${RUN_LOG}" >&2 || true
exit 1
