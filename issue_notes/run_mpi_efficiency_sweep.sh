#!/usr/bin/env bash
set -euo pipefail

REPO="${1:-/Users/jracine/Development/np-npRmpi}"
OUT_DIR="${2:-/tmp/mpi_efficiency_sweep_$(date +%Y%m%d_%H%M%S)}"

if [ ! -d "${REPO}/src" ] || [ ! -d "${REPO}/R" ]; then
  echo "ERROR: expected npRmpi repo root, got: ${REPO}" >&2
  exit 2
fi

mkdir -p "${OUT_DIR}"
SUMMARY="${OUT_DIR}/summary.md"
TOKENS="${OUT_DIR}/tokens.log"

cd "${REPO}"

{
  echo "OUT_DIR=${OUT_DIR}"
  echo "REPO=${REPO}"
} > "${TOKENS}"

# Raw scans
rg -n "npRmpi\.plot\.wild\.master_local\.guard|master-local|master\.local" R/np.plot.helpers.R > "${OUT_DIR}/scan_master_local_guard.txt" || true
rg -n "suppress\.parallel\s*=\s*TRUE" R/*.R > "${OUT_DIR}/scan_suppress_parallel_true_R.txt" || true
rg -n "if\(\(bwm == RBWM_CVLS\) \|\| ks_tree_use \|\| \(BANDWIDTH_reg == BW_ADAP_NN\)\)" src/jksum.c > "${OUT_DIR}/scan_cvls_forced_gate.txt" || true
rg -n "kernel_bandwidth_mean\(" src/jksum.c src/kernele.c src/np.c > "${OUT_DIR}/scan_kernel_bandwidth_calls.txt" || true
rg -n "do not suppress_parallel|suppress_parallel, // suppress_parallel if requested" src/jksum.c src/kernele.c src/np.c > "${OUT_DIR}/scan_suppress_parallel_c.txt" || true
rg -n "MPI_Allgather\(|MPI_Allreduce\(" src/jksum.c src/kernele.c src/statmods.c > "${OUT_DIR}/scan_collectives_core.txt" || true

# Metrics
count_lines() {
  local f="$1"
  if [ -s "$f" ]; then wc -l < "$f" | tr -d ' '; else echo 0; fi
}

MASTER_LOCAL_LINES=$(count_lines "${OUT_DIR}/scan_master_local_guard.txt")
SUPPRESS_TRUE_R_LINES=$(count_lines "${OUT_DIR}/scan_suppress_parallel_true_R.txt")
CVLS_FORCED_GATE_LINES=$(count_lines "${OUT_DIR}/scan_cvls_forced_gate.txt")
KERNEL_BW_CALL_LINES=$(count_lines "${OUT_DIR}/scan_kernel_bandwidth_calls.txt")
SUPPRESS_C_LINES=$(count_lines "${OUT_DIR}/scan_suppress_parallel_c.txt")
COLLECTIVE_CORE_LINES=$(count_lines "${OUT_DIR}/scan_collectives_core.txt")

{
  echo "MASTER_LOCAL_LINES=${MASTER_LOCAL_LINES}"
  echo "SUPPRESS_TRUE_R_LINES=${SUPPRESS_TRUE_R_LINES}"
  echo "CVLS_FORCED_GATE_LINES=${CVLS_FORCED_GATE_LINES}"
  echo "KERNEL_BW_CALL_LINES=${KERNEL_BW_CALL_LINES}"
  echo "SUPPRESS_C_LINES=${SUPPRESS_C_LINES}"
  echo "COLLECTIVE_CORE_LINES=${COLLECTIVE_CORE_LINES}"
} >> "${TOKENS}"

cat > "${SUMMARY}" <<EOF2
# MPI Efficiency Sweep (Static)

Repo: \
\`${REPO}\`

Artifacts:
- \`${OUT_DIR}/scan_master_local_guard.txt\`
- \`${OUT_DIR}/scan_suppress_parallel_true_R.txt\`
- \`${OUT_DIR}/scan_cvls_forced_gate.txt\`
- \`${OUT_DIR}/scan_kernel_bandwidth_calls.txt\`
- \`${OUT_DIR}/scan_suppress_parallel_c.txt\`
- \`${OUT_DIR}/scan_collectives_core.txt\`
- \`${TOKENS}\`

## Headline Counts
- master-local guard hits in plot helpers: **${MASTER_LOCAL_LINES}**
- R-layer \`suppress.parallel = TRUE\` hits: **${SUPPRESS_TRUE_R_LINES}**
- forced CVLS expensive-gate pattern hits in \`src/jksum.c\`: **${CVLS_FORCED_GATE_LINES}**
- kernel_bandwidth_mean callsites scanned: **${KERNEL_BW_CALL_LINES}**
- C-layer suppress-parallel annotation hits: **${SUPPRESS_C_LINES}**
- core collective callsites scanned (Allgather/Allreduce): **${COLLECTIVE_CORE_LINES}**

## Initial Classification
- \`CVLS_FORCED_GATE_LINES = 0\` indicates no reintroduced LL/LP CVLS forced gate pattern of the recent regression class.
- master-local and R suppress-parallel hits require context review (plot/bootstrap vs estimator/CV hot paths).
EOF2

cat "${TOKENS}"
echo "SUMMARY=${SUMMARY}"
echo "MPI_EFFICIENCY_SWEEP_OK"
