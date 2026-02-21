# Reproduce: npreg Current-vs-CRAN (npRmpi)

This document reproduces the default full-grid comparison for `npRmpi` with `n=100`, `times=50`, `nslaves=1` (128 combos).

## 1) Build isolated libraries

```bash
mkdir -p /tmp/Rlib_nprmpi_current /tmp/Rlib_nprmpi_cran20

R CMD INSTALL -l /tmp/Rlib_nprmpi_current /Users/jracine/Development/np-master
R CMD INSTALL --no-test-load -l /tmp/Rlib_nprmpi_current /Users/jracine/Development/np-npRmpi

R CMD INSTALL -l /tmp/Rlib_nprmpi_cran20 /Users/jracine/Development/CRAN/np_0.60-20.tar.gz
R CMD INSTALL --no-test-load -l /tmp/Rlib_nprmpi_cran20 /Users/jracine/Development/CRAN/npRmpi_0.60-20.tar.gz
```

## 2) Run default combo grid on CRAN and current

```bash
R_LIBS=/tmp/Rlib_nprmpi_cran20 FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npreg/run_npreg_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --nslaves=1 --tag=cran20_t50_n100

R_LIBS=/tmp/Rlib_nprmpi_current FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npreg/run_npreg_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --nslaves=1 --tag=current_t50_n100
```

Each run prints:

- `raw_csv=...`
- `summary_csv=...`
- `manifest_csv=...`

Use those `raw_csv` paths below.

## 3) Build comparison tables

```bash
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npreg/compare_npreg_versions.R \
  --raw_a=/tmp/RAW_CRAN.csv --label_a=npRmpi_0.60-20 \
  --raw_b=/tmp/RAW_CURRENT.csv --label_b=npRmpi_current \
  --out_timing=/tmp/nprmpi_combo32_t50_timing_compare.csv \
  --out_objective=/tmp/nprmpi_combo32_t50_objective_compare.csv \
  --out_combo_timing=/tmp/nprmpi_combo32_t50_combo_timing_compare.csv \
  --out_combo_objective=/tmp/nprmpi_combo32_t50_combo_objective_compare.csv
```

## 4) Generate markdown report

```bash
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npreg/make_npreg_report.R \
  --timing_csv=/tmp/nprmpi_combo32_t50_timing_compare.csv \
  --objective_csv=/tmp/nprmpi_combo32_t50_objective_compare.csv \
  --combo_timing_csv=/tmp/nprmpi_combo32_t50_combo_timing_compare.csv \
  --combo_objective_csv=/tmp/nprmpi_combo32_t50_combo_objective_compare.csv \
  --out_md=/tmp/nprmpi_combo32_t50_report.md \
  --title="npRmpi current vs CRAN (n=100, times=50, 32 combos)"
```

## 5) One-command summary check

```bash
Rscript -e 'print(read.csv("/tmp/nprmpi_combo32_t50_timing_compare.csv",stringsAsFactors=FALSE)); print(read.csv("/tmp/nprmpi_combo32_t50_objective_compare.csv",stringsAsFactors=FALSE))'
```

## 6) Generate combined np + npRmpi report

After producing the `np` comparison CSV files (see `/Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/REPRODUCE.md`), run:

```bash
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npreg/make_npreg_combined_report.R \
  --np_timing_csv=/tmp/np_combo32_t50_timing_compare.csv \
  --np_objective_csv=/tmp/np_combo32_t50_objective_compare.csv \
  --np_combo_timing_csv=/tmp/np_combo32_t50_combo_timing_compare.csv \
  --np_combo_objective_csv=/tmp/np_combo32_t50_combo_objective_compare.csv \
  --nprmpi_timing_csv=/tmp/nprmpi_combo32_t50_timing_compare.csv \
  --nprmpi_objective_csv=/tmp/nprmpi_combo32_t50_objective_compare.csv \
  --nprmpi_combo_timing_csv=/tmp/nprmpi_combo32_t50_combo_timing_compare.csv \
  --nprmpi_combo_objective_csv=/tmp/nprmpi_combo32_t50_combo_objective_compare.csv \
  --out_md=/tmp/npreg_combo32_t50_combined_report.md \
  --title="np + npRmpi current vs CRAN (n=100, times=50, 32 combos)"
```
