#!/usr/bin/env bash
set -euo pipefail

DEV_ROOT="${DEV_ROOT:-/Users/jracine/Development}"
PKG_DIR="${PKG_DIR:-${DEV_ROOT}/np-npRmpi}"
MAN_DIR="${MAN_DIR:-${PKG_DIR}/man}"
CHECK_ARGS="${CHECK_ARGS:---as-cran}"
FORCE_MPI_EXAMPLES="${FORCE_MPI_EXAMPLES:-0}"
TARBALL="${TARBALL:-}"
DRY_RUN="${DRY_RUN:-0}"

cleanup() {
  if [ -x "${MAN_DIR}/dontrun" ]; then
    (
      cd "${MAN_DIR}"
      ./dontrun
    ) || true
  fi
}

trap cleanup EXIT INT TERM

if [ ! -d "${PKG_DIR}" ]; then
  echo "ERROR: package directory not found: ${PKG_DIR}" >&2
  exit 2
fi

if [ ! -x "${MAN_DIR}/run" ]; then
  echo "ERROR: missing executable ${MAN_DIR}/run" >&2
  exit 2
fi

if [ "${DRY_RUN}" = "1" ]; then
  echo "DRY_RUN=1"
  echo "Would run: (cd ${MAN_DIR} && ./run)"
  echo "Would run: (cd ${DEV_ROOT} && R CMD build ${PKG_DIR})"
  if [ "${FORCE_MPI_EXAMPLES}" = "1" ]; then
    echo "Would run: NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1 R CMD check ${CHECK_ARGS} <tarball>"
  else
    echo "Would run: R CMD check ${CHECK_ARGS} <tarball>"
  fi
  echo "Would run via trap: (cd ${MAN_DIR} && ./dontrun)"
  exit 0
fi

(
  cd "${MAN_DIR}"
  ./run
)

(
  cd "${DEV_ROOT}"
  R CMD build "${PKG_DIR}"
)

if [ -z "${TARBALL}" ]; then
  TARBALL="$(find "${DEV_ROOT}" -maxdepth 1 -type f -name 'npRmpi_*.tar.gz' | sort | tail -n 1)"
fi

if [ -z "${TARBALL}" ]; then
  echo "ERROR: could not locate tarball (${TARBALL})" >&2
  exit 2
fi

if [ -f "${TARBALL}" ]; then
  CHECK_TARGET="${TARBALL}"
elif [ -f "${DEV_ROOT}/${TARBALL}" ]; then
  CHECK_TARGET="${DEV_ROOT}/${TARBALL}"
else
  echo "ERROR: could not locate tarball (${TARBALL})" >&2
  exit 2
fi

echo "CHECK_TARGET=${CHECK_TARGET}"
echo "FORCE_MPI_EXAMPLES=${FORCE_MPI_EXAMPLES}"

read -r -a CHECK_ARGS_ARR <<< "${CHECK_ARGS}"

if [ "${FORCE_MPI_EXAMPLES}" = "1" ]; then
  NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1 R CMD check "${CHECK_ARGS_ARR[@]}" "${CHECK_TARGET}"
else
  R CMD check "${CHECK_ARGS_ARR[@]}" "${CHECK_TARGET}"
fi

echo "CHECK_WITH_MAN_TRAP_OK"
