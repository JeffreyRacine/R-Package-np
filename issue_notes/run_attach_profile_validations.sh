#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${1:-/Users/jracine/Development/tmp_ll_lc_lp_alignment/nprmpi_route_validation_${STAMP}}"
INSTALL_LOG="${OUT_DIR}/install.log"
MANUAL_LOG="${OUT_DIR}/manual_broadcast.log"
ATTACH_LOG="${OUT_DIR}/attach.log"
PROFILE_LOG="${OUT_DIR}/profile.log"
TMP_LIB="$(mktemp -d /tmp/nprmpi_route_validation_lib.XXXXXX)"
PRIMARY_IFACE="${NP_RMPI_IFACE_PRIMARY:-en0}"
FALLBACK_IFACE="${NP_RMPI_IFACE_FALLBACK:-lo0}"
ROUTE_TIMEOUT_SEC="${NP_RMPI_ROUTE_TIMEOUT_SEC:-120}"

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
    echo "[warn] cleaning stray slave daemons: ${pids}"
    printf '%s\n' "${pids}" | xargs kill >/dev/null 2>&1 || true
  fi
}

cleanup() {
  kill_stray_mpi_slaves || true
  rm -rf "${TMP_LIB}"
}
trap cleanup EXIT

run_manual_route() {
  echo "[info] manual-broadcast route smoke" | tee -a "${MANUAL_LOG}"
  R_LIBS="${TMP_LIB}${R_LIBS:+:${R_LIBS}}" \
    run_with_timeout "${ROUTE_TIMEOUT_SEC}" \
    Rscript "${ROOT_DIR}/issue_notes/validate_route_manual_broadcast.R" >>"${MANUAL_LOG}" 2>&1
  rg -n "MANUAL_BCAST_ROUTE_OK" "${MANUAL_LOG}" >/dev/null
}

run_mpiexec_route() {
  local label="$1"
  local route_script="$2"
  local token="$3"
  local route_log="$4"
  local profile_env="${5:-0}"
  local mpiexec
  local rc=1
  local iface
  local profile_path
  local rscript_bin
  local r_bin
  local -a cmd
  local batch_out

  mpiexec="$(command -v mpiexec || true)"
  if [ -z "${mpiexec}" ]; then
    echo "[error] mpiexec not found in PATH" | tee -a "${route_log}"
    return 127
  fi

  if [ "${profile_env}" -eq 1 ]; then
    profile_path="${TMP_LIB}/npRmpi/Rprofile"
    if [ ! -f "${profile_path}" ]; then
      echo "[error] npRmpi profile template not found at ${profile_path}" | tee -a "${route_log}"
      return 2
    fi
  fi

  rscript_bin="$(R RHOME)/bin/Rscript"
  r_bin="$(R RHOME)/bin/R"

  for iface in "${PRIMARY_IFACE}" "${FALLBACK_IFACE}"; do
    echo "[info] ${label} route attempt on iface=${iface}" | tee -a "${route_log}"
    if [ "${profile_env}" -eq 1 ]; then
      batch_out="${route_log}.Rout"
      rm -f "${batch_out}"
      cmd=("${mpiexec}"
           -env R_PROFILE_USER "${profile_path}"
           -n 2
           "${r_bin}" CMD BATCH --no-save "${route_script}" "${batch_out}")
    else
      cmd=("${mpiexec}" -n 2 "${rscript_bin}" --no-save "${route_script}")
    fi

    set +e
    (
      export R_LIBS="${TMP_LIB}${R_LIBS:+:${R_LIBS}}"
      export FI_TCP_IFACE="${iface}"
      export FI_PROVIDER="tcp"
      export FI_SOCKETS_IFACE="${iface}"
      if [ "${profile_env}" -eq 1 ]; then
        export NP_RMPI_PROFILE_N="${NP_RMPI_PROFILE_N:-120}"
      fi
      run_with_timeout "${ROUTE_TIMEOUT_SEC}" \
        "${cmd[@]}"
    ) >>"${route_log}" 2>&1
    rc=$?
    if [ "${profile_env}" -eq 1 ] && [ -f "${batch_out}" ]; then
      cat "${batch_out}" >>"${route_log}"
    fi
    set -e

    if [ "${rc}" -eq 0 ] && rg -n "${token}" "${route_log}" >/dev/null; then
      return 0
    fi
    kill_stray_mpi_slaves
    if [ "${FALLBACK_IFACE}" = "${PRIMARY_IFACE}" ]; then
      break
    fi
  done

  return "${rc}"
}

echo "[info] installing npRmpi from ${ROOT_DIR}" | tee "${INSTALL_LOG}"
R CMD INSTALL --preclean -l "${TMP_LIB}" "${ROOT_DIR}" >>"${INSTALL_LOG}" 2>&1

run_manual_route
run_mpiexec_route "attach" \
  "${ROOT_DIR}/issue_notes/validate_route_attach.R" \
  "ATTACH_ROUTE_OK" \
  "${ATTACH_LOG}" \
  0
run_mpiexec_route "profile" \
  "${ROOT_DIR}/issue_notes/validate_route_profile.R" \
  "PROFILE_ROUTE_OK" \
  "${PROFILE_LOG}" \
  1

echo "[ok] route smoke checks passed"
echo "[ok] out dir: ${OUT_DIR}"
