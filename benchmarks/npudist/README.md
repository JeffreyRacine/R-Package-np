# npudist Benchmark Harness (npRmpi)

This folder provides a parameterized benchmark harness for `npudistbw()` + `npudist()` in `npRmpi` using `microbenchmark`.

## Files

- `bench_npudist_param_nprmpi.R`: run one benchmark configuration with `microbenchmark` repetitions.
- `run_npudist_combos.R`: run all 48 option combinations (`2x3x2x2x2`) with a shared setup.
- `compare_npudist_versions.R`: compare two raw outputs and emit overall and combo-specific comparison tables.
- `make_npudist_report.R`: generate a markdown report from comparison CSV outputs.
- `make_npudist_combined_report.R`: generate a single markdown report combining `np` and `npRmpi` comparison CSV outputs.
- `REPRODUCE.md`: end-to-end current-vs-CRAN replication steps.

## One-Configuration Run

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npudist/bench_npudist_param_nprmpi.R \
  --n=100 --times=50 --base_seed=42 \
  --bwmethod=cv.cdf --bwtype=fixed --nmulti=1 --ckertype=gaussian \
  --np_tree=FALSE --seed_policy=fixed --rslaves=1 \
  --out_raw=/tmp/npudist_one_mpi_raw.csv --out_summary=/tmp/npudist_one_mpi_summary.csv
```

## Full Combination Run (48 combos)

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npudist/run_npudist_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --rslaves=1 --tag=myrun
```

## Important Defaults

- `n=100`
- `times=50`
- `base_seed=42`
- `nmulti=1`
- `nslaves=1` (`--rslaves` alias supported)
- `bwmethod=cv.cdf`
- `ckertype=gaussian`
- `np_tree=FALSE`
- `seed_policy=varying`

Notes:

- `bwmethod` choices used in combos are `cv.cdf` and `normal-reference`.
- `ckertype` choices used in combos are `gaussian` and `epanechnikov`.

## DGP and Formula

The benchmark uses:

```r
set.seed(seed)
x <- runif(n)
z <- ordered(rbinom(n, 3, .5))
```

and fits `npudistbw(~ x + z, ...)` then `npudist(bws=bw, data=...)`.
