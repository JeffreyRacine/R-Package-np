# npudens Benchmark Harness (npRmpi)

This folder provides a parameterized benchmark harness for `npudensbw()` + `npudens()` in `npRmpi` using `microbenchmark`.

## Files

- `bench_npudens_param_nprmpi.R`: run one benchmark configuration with `microbenchmark` repetitions.
- `run_npudens_combos.R`: run all 16 option combinations (`2x2x2x2`) with a shared setup.
- `compare_npudens_versions.R`: compare two raw outputs and emit overall and combo-specific comparison tables.
- `make_npudens_report.R`: generate a markdown report from comparison CSV outputs.
- `make_npudens_combined_report.R`: generate a single markdown report combining `np` and `npRmpi` comparison CSV outputs.
- `REPRODUCE.md`: end-to-end current-vs-CRAN replication steps.

## One-Configuration Run

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npudens/bench_npudens_param_nprmpi.R \
  --n=100 --times=50 --base_seed=42 \
  --bwmethod=cv.ml --nmulti=1 --ckertype=gaussian \
  --np_tree=FALSE --seed_policy=fixed --rslaves=1 \
  --out_raw=/tmp/npudens_one_mpi_raw.csv --out_summary=/tmp/npudens_one_mpi_summary.csv
```

## Full Combination Run (16 combos)

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npudens/run_npudens_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --rslaves=1 --tag=myrun
```

## Important Defaults

- `n=100`
- `times=50`
- `base_seed=42`
- `nmulti=1`
- `nslaves=1` (`--rslaves` alias supported)
- `bwmethod=cv.ml`
- `ckertype=gaussian`
- `np_tree=FALSE`
- `seed_policy=varying`

Notes:

- `bwmethod` choices used in combos are `cv.ml` and `cv.ls`.
- `ckertype` choices used in combos are `gaussian` and `epanechnikov`.

## DGP and Formula

The benchmark uses:

```r
set.seed(seed)
x <- runif(n)
z <- factor(rbinom(n, 3, .5))
```

and fits `npudensbw(~ x + z, ...)` then `npudens(bws=bw, data=...)`.
