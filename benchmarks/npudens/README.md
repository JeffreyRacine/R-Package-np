# npudens Benchmark Harness (np)

This folder provides a parameterized benchmark harness for `npudensbw()` + `npudens()` in serial `np` using `microbenchmark`.

## Files

- `bench_npudens_param.R`: run one benchmark configuration with `microbenchmark` repetitions.
- `run_npudens_combos.R`: run all 16 option combinations (`2x2x2x2`) with a shared setup.
- `compare_npudens_versions.R`: compare two raw outputs and emit overall and combo-specific comparison tables.
- `make_npudens_report.R`: generate a markdown report from comparison CSV outputs.
- `make_npudens_combined_report.R`: generate a single markdown report combining `np` and `npRmpi` comparison CSV outputs.
- `REPRODUCE.md`: end-to-end current-vs-CRAN replication steps.

## One-Configuration Run

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/npudens/bench_npudens_param.R \
  --n=100 --times=50 --base_seed=42 \
  --bwmethod=cv.ml --nmulti=1 --ckertype=gaussian \
  --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npudens_one_raw.csv --out_summary=/tmp/npudens_one_summary.csv
```

## Full Combination Run (16 combos)

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/npudens/run_npudens_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=myrun
```

## Important Defaults

- `n=100`
- `times=50`
- `base_seed=42`
- `nmulti=1`
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
