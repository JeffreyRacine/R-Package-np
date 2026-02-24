#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_LOG="${1:-/tmp/nprmpi_native_bridge_stress_${STAMP}.log}"
INSTALL_LOG="${RUN_LOG%.log}.install.log"
TMP_LIB="$(mktemp -d /tmp/nprmpi_native_bridge_lib.XXXXXX)"
TIMEOUT_SEC="${NP_RMPI_NATIVE_STRESS_TIMEOUT_SEC:-240}"
PRIMARY_IFACE="${NP_RMPI_IFACE_PRIMARY:-en0}"
FALLBACK_IFACE="${NP_RMPI_IFACE_FALLBACK:-lo0}"

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

kill_stray_mpi_slaves() {
  local pids
  pids="$(pgrep -f 'slavedaemon\.R|Rslaves\.sh' || true)"
  if [ -n "${pids}" ]; then
    echo "[warn] cleaning stray slave daemons: ${pids}" | tee -a "${RUN_LOG}"
    printf '%s\n' "${pids}" | xargs kill >/dev/null 2>&1 || true
  fi
}

# shellcheck disable=SC2329 # Invoked via trap.
cleanup() {
  rm -rf "${TMP_LIB}"
}
trap cleanup EXIT

run_native_stress_once() {
  local iface="$1"
  echo "[info] native bridge stress iface=${iface} timeout=${TIMEOUT_SEC}s" | tee -a "${RUN_LOG}"
  R_LIBS="${TMP_LIB}${R_LIBS:+:${R_LIBS}}" \
    FI_TCP_IFACE="${iface}" FI_PROVIDER=tcp FI_SOCKETS_IFACE="${iface}" \
    run_with_timeout "${TIMEOUT_SEC}" \
    Rscript "${ROOT_DIR}/issue_notes/native_bridge_stress.R" >>"${RUN_LOG}" 2>&1
}

echo "[info] installing npRmpi from ${ROOT_DIR}" | tee "${INSTALL_LOG}"
R CMD INSTALL --preclean -l "${TMP_LIB}" "${ROOT_DIR}" >>"${INSTALL_LOG}" 2>&1

echo "[info] running native bridge stress harness" | tee "${RUN_LOG}"

set +e
run_native_stress_once "${PRIMARY_IFACE}"
rc=$?
set -e

if [ "${rc}" -ne 0 ] && [ "${FALLBACK_IFACE}" != "${PRIMARY_IFACE}" ]; then
  echo "[warn] native bridge stress failed on iface=${PRIMARY_IFACE} rc=${rc}; retrying iface=${FALLBACK_IFACE}" | tee -a "${RUN_LOG}"
  kill_stray_mpi_slaves
  set +e
  run_native_stress_once "${FALLBACK_IFACE}"
  rc=$?
  set -e
fi

if [ "${rc}" -ne 0 ]; then
  kill_stray_mpi_slaves
  echo "[error] native bridge stress command failed rc=${rc}" >&2
  tail -n 80 "${RUN_LOG}" >&2 || true
  exit "${rc}"
fi

if rg -n "NATIVE_BRIDGE_STRESS_NPRMPI_OK" "${RUN_LOG}" >/dev/null 2>&1; then
  echo "[ok] native bridge stress passed"
  echo "[ok] run log: ${RUN_LOG}"
  echo "[ok] install log: ${INSTALL_LOG}"
  exit 0
fi

echo "[error] native bridge stress did not report success token" >&2
tail -n 80 "${RUN_LOG}" >&2 || true
exit 1
