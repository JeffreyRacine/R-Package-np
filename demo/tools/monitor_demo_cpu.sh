#!/bin/bash
set -uo pipefail

if [ "$#" -lt 4 ]; then
  echo "usage: $0 <run_pid> <run_log> <run_root> <out_csv> [interval_sec]" >&2
  exit 2
fi

RUN_PID="$1"
RUN_LOG="$2"
RUN_ROOT="$3"
OUT_CSV="$4"
INTERVAL="${5:-10}"

mkdir -p "$(dirname "${OUT_CSV}")"

if [ ! -f "${OUT_CSV}" ]; then
  printf 'timestamp,run_pid,run_alive,active_stage,pid,ppid,elapsed,pcpu,role,args\n' > "${OUT_CSV}"
fi

csv_escape() {
  local value="${1//$'\n'/ }"
  value="${value//\"/\"\"}"
  printf '"%s"' "${value}"
}

classify_role() {
  local args="$1"
  case "${args}" in
    *slavedaemon.R*) printf 'worker' ;;
    *mpiexec*) printf 'mpiexec' ;;
    /Library/Frameworks/R.framework/*/exec/R*|/Library/Frameworks/R.framework/*/Rscript*|Rscript\ *) printf 'master_or_serial_R' ;;
    *runall*|*makefile*) printf 'controller' ;;
    *) printf 'other' ;;
  esac
}

while kill -0 "${RUN_PID}" >/dev/null 2>&1; do
  ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  active_stage="$(tail -1 "${RUN_LOG}" 2>/dev/null || true)"
  run_alive=1

  ps -axo pid=,ppid=,etime=,pcpu=,args= |
    grep -F "${RUN_ROOT}" |
    grep -v 'grep -F' |
    while read -r pid ppid elapsed pcpu args; do
      role="$(classify_role "${args}")"
      {
        csv_escape "${ts}"; printf ','
        csv_escape "${RUN_PID}"; printf ','
        csv_escape "${run_alive}"; printf ','
        csv_escape "${active_stage}"; printf ','
        csv_escape "${pid}"; printf ','
        csv_escape "${ppid}"; printf ','
        csv_escape "${elapsed}"; printf ','
        csv_escape "${pcpu}"; printf ','
        csv_escape "${role}"; printf ','
        csv_escape "${args}"; printf '\n'
      } >> "${OUT_CSV}"
    done || true

  ps -axo pid=,ppid=,etime=,pcpu=,args= |
    awk '{
      args = $0
      sub(/^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]+[^[:space:]]+[[:space:]]+[^[:space:]]+[[:space:]]+/, "", args)
      if (args ~ /^mpiexec -pmi_args/ ||
          args ~ /^\/bin\/sh .*\/npRmpi\/Rslaves[.]sh / ||
          args ~ /^\/Library\/Frameworks\/R[.]framework\/.*\/exec\/R .*\/npRmpi\/slavedaemon[.]R([[:space:]]|$)/) print
    }' |
    while read -r pid ppid elapsed pcpu args; do
      role="$(classify_role "${args}")"
      {
        csv_escape "${ts}"; printf ','
        csv_escape "${RUN_PID}"; printf ','
        csv_escape "${run_alive}"; printf ','
        csv_escape "${active_stage}"; printf ','
        csv_escape "${pid}"; printf ','
        csv_escape "${ppid}"; printf ','
        csv_escape "${elapsed}"; printf ','
        csv_escape "${pcpu}"; printf ','
        csv_escape "${role}"; printf ','
        csv_escape "${args}"; printf '\n'
      } >> "${OUT_CSV}"
    done || true

  sleep "${INTERVAL}"
done

ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
{
  csv_escape "${ts}"; printf ','
  csv_escape "${RUN_PID}"; printf ','
  csv_escape "0"; printf ','
  csv_escape "run_finished"; printf ',,,,,,\n'
} >> "${OUT_CSV}"
