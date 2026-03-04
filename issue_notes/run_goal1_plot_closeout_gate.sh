#!/usr/bin/env bash
set -euo pipefail

REPO="/Users/jracine/Development/np-master"
BASE_REF="${1:-checkpoint-20260304_160029-np-master}"
TS="$(date +%Y%m%d_%H%M%S)"
OUT="/tmp/goal1_closeout_${TS}"
mkdir -p "$OUT"/{tests,matrix,perf}

{
  echo "repo=$REPO"
  echo "base_ref=$BASE_REF"
  echo "timestamp=$TS"
  echo "branch=$(git -C "$REPO" rev-parse --abbrev-ref HEAD)"
  echo "sha=$(git -C "$REPO" rev-parse HEAD)"
} > "$OUT/manifest.txt"

TEST_FILTERS=(
  "plot-bootstrap-inid-fastpath"
  "plot-conditional-gradients-bootstrap"
  "plot-contract"
  "plot-guardrails-contract"
  "plot"
)

TEST_SUMMARY="$OUT/tests/test_summary.csv"
echo "filter,fail,warn,skip,pass,log" > "$TEST_SUMMARY"

for f in "${TEST_FILTERS[@]}"; do
  log="$OUT/tests/test_${f}.log"
  Rscript -e "devtools::test(filter='${f}')" > "$log" 2>&1
  line="$(rg -n "\\[ FAIL [0-9]+ \\| WARN [0-9]+ \\| SKIP [0-9]+ \\| PASS [0-9]+ \\]" "$log" | tail -1 | sed -E 's/^[^\\[]*//')"
  if [[ -z "$line" ]]; then
    echo "missing test summary line for filter=${f}" >&2
    exit 1
  fi
  fail="$(echo "$line" | sed -E 's/.*FAIL ([0-9]+).*/\1/')"
  warn="$(echo "$line" | sed -E 's/.*WARN ([0-9]+).*/\1/')"
  skip="$(echo "$line" | sed -E 's/.*SKIP ([0-9]+).*/\1/')"
  pass="$(echo "$line" | sed -E 's/.*PASS ([0-9]+).*/\1/')"
  echo "${f},${fail},${warn},${skip},${pass},${log}" >> "$TEST_SUMMARY"
  if [[ "$fail" != "0" ]]; then
    echo "test failure for filter=${f}" >&2
    exit 1
  fi
done

rm -f "$REPO/tests/testthat/Rplots.pdf"

CLASSIFY_LOG="$OUT/classify_exit134.log"
bash "$REPO/issue_notes/run_classify_exit134_npmaster.sh" > "$CLASSIFY_LOG" 2>&1
CLASSIFY_CSV="$(awk -F': ' '/Results CSV:/ {print $2}' "$CLASSIFY_LOG" | tail -1)"
CLASSIFY_LINE="$(awk '/CLASSIFICATION:/{print $0}' "$CLASSIFY_LOG" | tail -1)"
if [[ -z "$CLASSIFY_CSV" || -z "$CLASSIFY_LINE" ]]; then
  echo "failed to parse exit134 classification output" >&2
  exit 1
fi
if [[ "$CLASSIFY_LINE" != *"DOES NOT REPRODUCE in np-master"* ]]; then
  echo "unexpected classification outcome: $CLASSIFY_LINE" >&2
  exit 1
fi

MATRIX_LOG="$OUT/matrix/matrix.log"
Rscript "$REPO/issue_notes/goal1_plot_closeout_matrix.R" \
  --repo "$REPO" \
  --out-dir "$OUT/matrix" \
  --boot-num 79 \
  --blocklen 5 \
  > "$MATRIX_LOG" 2>&1

WT="/tmp/np_master_goal1_closeout_baseline_${TS}"
cleanup() {
  set +e
  git -C "$REPO" worktree remove --force "$WT" >/dev/null 2>&1 || true
}
trap cleanup EXIT

git -C "$REPO" worktree add --detach "$WT" "$BASE_REF" > "$OUT/perf/worktree_add.log" 2>&1

SEEDS="$(seq 1 25 | paste -sd, -)"
BASE_CSV="$OUT/perf/baseline_runs.csv"
CAND_CSV="$OUT/perf/candidate_runs.csv"

Rscript "$REPO/issue_notes/goal1_plot_closeout_perf_runs.R" \
  --repo "$WT" \
  --label baseline \
  --out-csv "$BASE_CSV" \
  --boot-num 79 \
  --blocklen 5 \
  --seeds "$SEEDS" \
  > "$OUT/perf/baseline.log" 2>&1

Rscript "$REPO/issue_notes/goal1_plot_closeout_perf_runs.R" \
  --repo "$REPO" \
  --label candidate \
  --out-csv "$CAND_CSV" \
  --boot-num 79 \
  --blocklen 5 \
  --seeds "$SEEDS" \
  > "$OUT/perf/candidate.log" 2>&1

Rscript "$REPO/issue_notes/goal1_plot_closeout_perf_summary.R" \
  --baseline-csv "$BASE_CSV" \
  --candidate-csv "$CAND_CSV" \
  --out-dir "$OUT/perf" \
  > "$OUT/perf/summary.log" 2>&1

cleanup
trap - EXIT

{
  echo "test_summary_csv=$TEST_SUMMARY"
  echo "classify_log=$CLASSIFY_LOG"
  echo "classify_csv=$CLASSIFY_CSV"
  echo "classify_result=$CLASSIFY_LINE"
  echo "matrix_summary=$OUT/matrix/matrix_summary.txt"
  echo "matrix_csv=$OUT/matrix/matrix_results.csv"
  echo "perf_verdict=$OUT/perf/perf_verdict.txt"
  echo "perf_summary_csv=$OUT/perf/perf_summary.csv"
  echo "perf_pairs_csv=$OUT/perf/perf_pairs.csv"
  echo "result=PASS"
} >> "$OUT/manifest.txt"

echo "GOAL1_CLOSEOUT_DIR=$OUT"
echo "GOAL1_CLOSEOUT_MANIFEST=$OUT/manifest.txt"
