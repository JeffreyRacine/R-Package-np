# npreg Benchmark Harness (np)

This folder provides a parameterized benchmark harness for `npregbw()` + `npreg()` in serial `np` using `microbenchmark`.

## Files

- `bench_npreg_param.R`: run one benchmark configuration with `microbenchmark` repetitions.
- `run_npreg_combos.R`: run all 32 option combinations (`2x2x2x2x2`) with a shared setup.
- `compare_npreg_versions.R`: compare two raw outputs and emit overall and combo-specific comparison tables.
- `make_npreg_report.R`: generate a markdown report from comparison CSV outputs.
- `make_npreg_combined_report.R`: generate a single markdown report combining `np` and `npRmpi` comparison CSV outputs.
- `REPRODUCE.md`: end-to-end current-vs-CRAN replication steps.

## One-Configuration Run

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/bench_npreg_param.R \
  --n=100 --times=50 --base_seed=42 \
  --regtype=lc --bwmethod=cv.ls --nmulti=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npreg_one_raw.csv --out_summary=/tmp/npreg_one_summary.csv
```

## Full Combination Run (32 combos)

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/run_npreg_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=myrun
```

Outputs are written to `/tmp` with run IDs in filenames:

- `npreg_combo_raw_...csv`
- `npreg_combo_summary_...csv`
- `npreg_combo_manifest_...csv`

## Important Defaults

- `n=100`
- `times=50`
- `base_seed=42`
- `nmulti=1`
- `regtype=ll`
- `bwmethod=cv.ls`
- `ckertype=gaussian`
- `np_tree=FALSE`
- `seed_policy=varying`

Change sample size and repetitions via `--n=` and `--times=`.

## Timing/Seed Semantics

- `times` is passed to `microbenchmark(times=...)`.
- One case seed is selected (`--base_seed` unless `--seeds` is given).
- For `seed_policy=fixed`, each microbenchmark repetition uses the same seed.
- For `seed_policy=varying`, repetition `i` uses `seed + (i-1)` (deterministic and reproducible across versions).

## Compatibility Note

`num_fval` is extracted from `bw$num.fval` when available. If missing (older versions), it is reported as `NA`.

## Version Comparison Workflow (np)

1. Install each target version into its own library.

```bash
mkdir -p /tmp/Rlib_np_current /tmp/Rlib_np_cran20
R CMD INSTALL -l /tmp/Rlib_np_current /Users/jracine/Development/np-master
R CMD INSTALL -l /tmp/Rlib_np_cran20 /Users/jracine/Development/CRAN/np_0.60-20.tar.gz
```

2. Run canonical benchmark with each library.

```bash
R_LIBS=/tmp/Rlib_np_cran20 Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/bench_npreg_param.R \
  --n=100 --times=50 --out_raw=/tmp/np_cran20_raw.csv --out_summary=/tmp/np_cran20_summary.csv

R_LIBS=/tmp/Rlib_np_current Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/bench_npreg_param.R \
  --n=100 --times=50 --out_raw=/tmp/np_current_raw.csv --out_summary=/tmp/np_current_summary.csv
```

3. Compare outputs.

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npreg/compare_npreg_versions.R \
  --raw_a=/tmp/np_cran20_raw.csv --label_a=np_0.60-20 \
  --raw_b=/tmp/np_current_raw.csv --label_b=np_current \
  --out_timing=/tmp/np_timing_compare.csv \
  --out_objective=/tmp/np_objective_compare.csv \
  --out_combo_timing=/tmp/np_combo_timing_compare.csv \
  --out_combo_objective=/tmp/np_combo_objective_compare.csv
```

Comparison outputs:

- Timing by function (`npregbw`, `npreg`, `npreg_total`) with mean/median and percent change.
- Objective diagnostics (`fval`, `ifval`, `num_fval`, bandwidth match rate, `ok` match rate).
- Combo timing table by (`regtype`, `bwmethod`, `ckertype`, `np_tree`, `seed_policy`, `nmulti`) with mean/median percent change.
- Combo objective diagnostics by the same combination keys.
