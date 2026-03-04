#!/usr/bin/env bash
set -u
set -o pipefail

REPO="/Users/jracine/Development/np-master"
SCRIPT="${REPO}/issue_notes/classify_exit134_npmaster_once.R"
TS="$(date +%Y%m%d_%H%M%S)"
OUT="/tmp/exit134_classify_${TS}"
mkdir -p "$OUT"

CSV="${OUT}/runs.csv"
echo "n,seed,exit_code,max_rss_kb,elapsed_sec,nan_inf,warn_count,status" > "$CSV"

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

    elapsed=$(grep -Eo 'elapsed_sec=[0-9]+(\.[0-9]+)?' "${run_dir}/stdout.log" | tail -1 | cut -d= -f2)
    [ -z "${elapsed}" ] && elapsed="NA"

    nan_inf=$(grep -Eo 'nan_inf=(TRUE|FALSE)' "${run_dir}/stdout.log" | tail -1 | cut -d= -f2)
    [ -z "${nan_inf}" ] && nan_inf="NA"

    warn_count=$(grep -Eo 'warn_count=[0-9]+' "${run_dir}/stdout.log" | tail -1 | cut -d= -f2)
    [ -z "${warn_count}" ] && warn_count="NA"

    status=$(grep -Eo 'status=[A-Za-z_]+' "${run_dir}/stdout.log" | tail -1 | cut -d= -f2)
    [ -z "${status}" ] && status="ABORT"

    echo "${n},${seed},${ec},${rss},${elapsed},${nan_inf},${warn_count},${status}" >> "$CSV"
  done
done

echo "Results CSV: $CSV"

if awk -F, 'NR>1 && ($3 != 0 || $6 == "TRUE" || $8 != "ok"){bad=1} END{exit(bad?0:1)}' "$CSV"; then
  echo "CLASSIFICATION: REPRODUCES in np-master (Goal 1 scope)"
else
  echo "CLASSIFICATION: DOES NOT REPRODUCE in np-master (Goal 2 scope)"
fi
