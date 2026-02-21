# Reproduce: npcdens Current-vs-CRAN (np)

This document reproduces the 16-combo comparison for `npcdens` in `np` with `n=100`, `times=50`.

## 1) Build isolated libraries

```bash
mkdir -p /tmp/Rlib_np_current /tmp/Rlib_np_cran20
R CMD INSTALL -l /tmp/Rlib_np_current /Users/jracine/Development/np-master
R CMD INSTALL -l /tmp/Rlib_np_cran20 /Users/jracine/Development/CRAN/np_0.60-20.tar.gz
```

## 2) Run 16 combos on CRAN and current

```bash
R_LIBS=/tmp/Rlib_np_cran20 \
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npcdens/run_npcdens_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=cran20_t50_n100

R_LIBS=/tmp/Rlib_np_current \
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npcdens/run_npcdens_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=current_t50_n100
```

## 3) Build comparison tables

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npcdens/compare_npcdens_versions.R \
  --raw_a=/tmp/RAW_CRAN.csv --label_a=np_0.60-20 \
  --raw_b=/tmp/RAW_CURRENT.csv --label_b=np_current \
  --out_timing=/tmp/npcdens_combo16_t50_timing_compare.csv \
  --out_objective=/tmp/npcdens_combo16_t50_objective_compare.csv \
  --out_combo_timing=/tmp/npcdens_combo16_t50_combo_timing_compare.csv \
  --out_combo_objective=/tmp/npcdens_combo16_t50_combo_objective_compare.csv
```

## 4) Generate markdown report

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npcdens/make_npcdens_report.R \
  --timing_csv=/tmp/npcdens_combo16_t50_timing_compare.csv \
  --objective_csv=/tmp/npcdens_combo16_t50_objective_compare.csv \
  --combo_timing_csv=/tmp/npcdens_combo16_t50_combo_timing_compare.csv \
  --combo_objective_csv=/tmp/npcdens_combo16_t50_combo_objective_compare.csv \
  --out_md=/tmp/npcdens_combo16_t50_report.md \
  --title="npcdens np current vs CRAN (n=100, times=50, 16 combos)"
```
