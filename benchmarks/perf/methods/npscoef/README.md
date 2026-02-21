# npscoef Benchmark Harness (np)

Parameterized microbenchmark harness for `npscoefbw()` + `npscoef()` (serial `np`).

## Grid

`run_npscoef_combos.R` runs 8 combinations (`2x2x2`):

- `ckertype`: `gaussian`, `epanechnikov`
- `np_tree`: `TRUE`, `FALSE`
- `seed_policy`: `fixed`, `varying`

There is no `bwmethod` option for this harness.

## DGP

```r
x <- runif(n)
z <- runif(n, min = -2, max = 2)
y <- x * exp(z) * (1.0 + rnorm(n, sd = 0.2))
```

Bandwidth call:

```r
npscoefbw(y ~ x | z, ...)
```

## One Case

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npscoef/bench_npscoef_param.R \
  --n=250 --times=50 --base_seed=42 --nmulti=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npscoef_one_raw.csv --out_summary=/tmp/npscoef_one_summary.csv
```

## Combo Run

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npscoef/run_npscoef_combos.R \
  --n=250 --times=50 --base_seed=42 --nmulti=1 --tag=myrun
```
