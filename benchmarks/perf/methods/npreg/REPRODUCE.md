# Reproduce: npreg Current-vs-CRAN (np)

This document reproduces the default full-grid comparison for `np` with `n=100`, `times=50` (128 combos).

Optional nearest-neighbor sweep:

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/run_npreg_combos.R \
  --n=40 --times=1 --base_seed=42 --nmulti=1 \
  --bwtypes=generalized_nn,adaptive_nn --max_combos=8 --tag=nn_probe
```

## 1) Build isolated libraries

```bash
mkdir -p /tmp/Rlib_np_current /tmp/Rlib_np_cran20
R CMD INSTALL -l /tmp/Rlib_np_current /Users/jracine/Development/np-master
R CMD INSTALL -l /tmp/Rlib_np_cran20 /Users/jracine/Development/CRAN/np_0.60-20.tar.gz
```

## 2) Run default combo grid on CRAN and current

```bash
R_LIBS=/tmp/Rlib_np_cran20 \
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/run_npreg_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=cran20_t50_n100

R_LIBS=/tmp/Rlib_np_current \
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/run_npreg_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=current_t50_n100
```

Each run prints:

- `raw_csv=...`
- `summary_csv=...`
- `manifest_csv=...`

Use those `raw_csv` paths below.

## 3) Build comparison tables

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/compare_npreg_versions.R \
  --raw_a=/tmp/RAW_CRAN.csv --label_a=np_0.60-20 \
  --raw_b=/tmp/RAW_CURRENT.csv --label_b=np_current \
  --out_timing=/tmp/np_combo32_t50_timing_compare.csv \
  --out_objective=/tmp/np_combo32_t50_objective_compare.csv \
  --out_combo_timing=/tmp/np_combo32_t50_combo_timing_compare.csv \
  --out_combo_objective=/tmp/np_combo32_t50_combo_objective_compare.csv
```

## 4) Generate markdown report

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/make_npreg_report.R \
  --timing_csv=/tmp/np_combo32_t50_timing_compare.csv \
  --objective_csv=/tmp/np_combo32_t50_objective_compare.csv \
  --combo_timing_csv=/tmp/np_combo32_t50_combo_timing_compare.csv \
  --combo_objective_csv=/tmp/np_combo32_t50_combo_objective_compare.csv \
  --out_md=/tmp/np_combo32_t50_report.md \
  --title="np current vs CRAN (n=100, times=50, 32 combos)"
```

## 5) One-command summary check

```bash
Rscript -e 'print(read.csv("/tmp/np_combo32_t50_timing_compare.csv",stringsAsFactors=FALSE)); print(read.csv("/tmp/np_combo32_t50_objective_compare.csv",stringsAsFactors=FALSE))'
```

## 6) Generate combined np + npRmpi report

After producing the `npRmpi` comparison CSV files (see `/Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npreg/REPRODUCE.md`), run:

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/make_npreg_combined_report.R \
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

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
