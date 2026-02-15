# npreg Benchmark Harness (npRmpi)

This folder provides a parameterized benchmark harness for `npregbw()` + `npreg()` in `npRmpi` using `microbenchmark`.

## Files

- `bench_npreg_param_nprmpi.R`: run one benchmark configuration with `microbenchmark` repetitions.
- `run_npreg_combos.R`: run all 32 option combinations (`2x2x2x2x2`) with a shared setup.

## One-Configuration Run

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npreg/bench_npreg_param_nprmpi.R \
  --n=100 --times=1 --base_seed=42 \
  --regtype=lc --bwmethod=cv.ls --nmulti=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --rslaves=1 \
  --out_raw=/tmp/npreg_one_mpi_raw.csv --out_summary=/tmp/npreg_one_mpi_summary.csv
```

## Full Combination Run (32 combos)

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npreg/run_npreg_combos.R \
  --n=100 --times=1 --base_seed=42 --nmulti=1 --rslaves=1 --tag=myrun
```

Outputs are written to `/tmp` with run IDs in filenames:

- `npreg_combo_raw_...csv`
- `npreg_combo_summary_...csv`
- `npreg_combo_manifest_...csv`
- Per-combo logs: `npreg_combo_<runid>_<combo>.log`

## Important Defaults

- `n=100`
- `times=1`
- `base_seed=42`
- `nmulti=1`
- `nslaves=1` (`--rslaves` alias supported)
- `regtype=lc`
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

## Version Comparison Workflow (npRmpi)

1. Install each target stack into its own library.

```bash
mkdir -p /tmp/Rlib_nprmpi_current /tmp/Rlib_nprmpi_cran20
R CMD INSTALL -l /tmp/Rlib_nprmpi_current /Users/jracine/Development/np-master
R CMD INSTALL --no-test-load -l /tmp/Rlib_nprmpi_current /Users/jracine/Development/np-npRmpi

R CMD INSTALL -l /tmp/Rlib_nprmpi_cran20 /Users/jracine/Development/CRAN/np_0.60-20.tar.gz
R CMD INSTALL --no-test-load -l /tmp/Rlib_nprmpi_cran20 /Users/jracine/Development/CRAN/npRmpi_0.60-20.tar.gz
```

2. Run canonical benchmark with each stack.

```bash
R_LIBS=/tmp/Rlib_nprmpi_cran20 FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npreg/bench_npreg_param_nprmpi.R \
  --n=100 --times=5 --rslaves=1 \
  --out_raw=/tmp/nprmpi_cran20_raw.csv --out_summary=/tmp/nprmpi_cran20_summary.csv

R_LIBS=/tmp/Rlib_nprmpi_current FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npreg/bench_npreg_param_nprmpi.R \
  --n=100 --times=5 --rslaves=1 \
  --out_raw=/tmp/nprmpi_current_raw.csv --out_summary=/tmp/nprmpi_current_summary.csv
```

3. Compare outputs.

```bash
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npreg/compare_npreg_versions.R \
  --raw_a=/tmp/nprmpi_cran20_raw.csv --label_a=npRmpi_0.60-20 \
  --raw_b=/tmp/nprmpi_current_raw.csv --label_b=npRmpi_current \
  --out_timing=/tmp/nprmpi_timing_compare.csv \
  --out_objective=/tmp/nprmpi_objective_compare.csv
```

Comparison outputs:

- Timing by function (`npregbw`, `npreg`, `npreg_total`) with mean/median and percent change.
- Objective diagnostics (`fval`, `ifval`, `num_fval`, bandwidth match rate, `ok` match rate).
