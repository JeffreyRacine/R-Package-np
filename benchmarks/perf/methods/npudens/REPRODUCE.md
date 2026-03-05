# Reproduce: npudens Current-vs-CRAN (np)

This document reproduces the 16-combo comparison for `npudens` in `np` with `n=100`, `times=50`.

## 1) Build isolated libraries

```bash
mkdir -p /tmp/Rlib_np_current /tmp/Rlib_np_cran20
R CMD INSTALL -l /tmp/Rlib_np_current /Users/jracine/Development/np-master
R CMD INSTALL -l /tmp/Rlib_np_cran20 /Users/jracine/Development/CRAN/np_0.60-20.tar.gz
```

## 2) Run 16 combos on CRAN and current

```bash
R_LIBS=/tmp/Rlib_np_cran20 \
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npudens/run_npudens_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=cran20_t50_n100

R_LIBS=/tmp/Rlib_np_current \
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npudens/run_npudens_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=current_t50_n100
```

## 3) Build comparison tables

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npudens/compare_npudens_versions.R \
  --raw_a=/tmp/RAW_CRAN.csv --label_a=np_0.60-20 \
  --raw_b=/tmp/RAW_CURRENT.csv --label_b=np_current \
  --out_timing=/tmp/npudens_combo16_t50_timing_compare.csv \
  --out_objective=/tmp/npudens_combo16_t50_objective_compare.csv \
  --out_combo_timing=/tmp/npudens_combo16_t50_combo_timing_compare.csv \
  --out_combo_objective=/tmp/npudens_combo16_t50_combo_objective_compare.csv
```

## 4) Generate markdown report

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npudens/make_npudens_report.R \
  --timing_csv=/tmp/npudens_combo16_t50_timing_compare.csv \
  --objective_csv=/tmp/npudens_combo16_t50_objective_compare.csv \
  --combo_timing_csv=/tmp/npudens_combo16_t50_combo_timing_compare.csv \
  --combo_objective_csv=/tmp/npudens_combo16_t50_combo_objective_compare.csv \
  --out_md=/tmp/npudens_combo16_t50_report.md \
  --title="npudens np current vs CRAN (n=100, times=50, 16 combos)"
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
