# npindex Benchmark Harness (np)

Parameterized microbenchmark harness for `npindexbw()` + `npindex()` (serial `np`).

## Grid

`run_npindex_combos.R` runs 16 combinations (`2x2x2x2`):

- `method`: `ichimura`, `kleinspady`
- `ckertype`: `gaussian`, `epanechnikov`
- `np_tree`: `TRUE`, `FALSE`
- `seed_policy`: `fixed`, `varying`

## DGP

`method="ichimura"`:

```r
x1 <- runif(n, min = -1, max = 1)
x2 <- runif(n, min = -1, max = 1)
y <- x1 - x2 + rnorm(n)
```

`method="kleinspady"` (binary response stays integer, not factor):

```r
x1 <- runif(n, min = -1, max = 1)
x2 <- runif(n, min = -1, max = 1)
y <- ifelse(x1 + x2 + rnorm(n) > 0, 1L, 0L)
```

## One Case

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npindex/bench_npindex_param.R \
  --n=100 --times=50 --base_seed=42 --method=ichimura --nmulti=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npindex_one_raw.csv --out_summary=/tmp/npindex_one_summary.csv
```

## Combo Run

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npindex/run_npindex_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=myrun
```
