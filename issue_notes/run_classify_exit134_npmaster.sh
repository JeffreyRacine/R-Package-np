#!/usr/bin/env bash
set -u
set -o pipefail

REPO="/Users/jracine/Development/np-master"
SCRIPT="${REPO}/issue_notes/classify_exit134_npmaster_once.R"
TS="$(date +%Y%m%d_%H%M%S)"
OUT="/tmp/exit134_classify_${TS}"
mkdir -p "$OUT"

CSV="${OUT}/runs.csv"
echo "n,seed,exit_code,max_rss_kb,elapsed_sec,nonfinite,warn_count,status" > "$CSV"

NS=(100 500 1000)
SEEDS=(42 42 42 42 42 1 2 3 4 5)
BOOT_NUM=399

for n in "${NS[@]}"; do
  for seed in "${SEEDS[@]}"; do
    run_dir="${OUT}/n${n}_s${seed}_$RANDOM"
    mkdir -p "$run_dir"

    /usr/bin/time -l \
      Rscript "$SCRIPT" --repo "$REPO" --n "$n" --seed "$seed" --boot-num "$BOOT_NUM" --out-dir "$run_dir" \
      > "${run_dir}/stdout.log" 2> "${run_dir}/stderr.log"
    ec=$?

    rss=$(awk '/maximum resident set size/{print $1}' "${run_dir}/stderr.log" | tail -1)
    [ -z "${rss}" ] && rss="NA"

    elapsed=$(awk -F= '/^elapsed_sec=/{print $2}' "${run_dir}/stdout.log" | tail -1)
    [ -z "${elapsed}" ] && elapsed="NA"

    nonfinite=$(awk -F= '/^nonfinite=/{print $2}' "${run_dir}/stdout.log" | tail -1)
    [ -z "${nonfinite}" ] && nonfinite="NA"

    warn_count=$(awk -F= '/^warn_count=/{print $2}' "${run_dir}/stdout.log" | tail -1)
    [ -z "${warn_count}" ] && warn_count="NA"

    status=$(awk -F= '/^status=/{print $2}' "${run_dir}/stdout.log" | tail -1)
    [ -z "${status}" ] && status="ABORT"

    echo "${n},${seed},${ec},${rss},${elapsed},${nonfinite},${warn_count},${status}" >> "$CSV"
  done
done

echo "Results CSV: $CSV"

if awk -F, 'NR>1 && ($3 != 0 || $6 == "TRUE" || $8 != "ok"){bad=1} END{exit(bad?0:1)}' "$CSV"; then
  echo "CLASSIFICATION: REPRODUCES in np-master (Goal 1 scope)"
else
  echo "CLASSIFICATION: DOES NOT REPRODUCE in np-master (Goal 2 scope)"
fi
