# npscoef Benchmark Harness (npRmpi)

Parameterized microbenchmark harness for `npscoefbw()` + `npscoef()` (`npRmpi`).

## Grid

`run_npscoef_combos.R` runs 8 combinations (`2x2x2`):

- `ckertype`: `gaussian`, `epanechnikov`
- `np_tree`: `TRUE`, `FALSE`
- `seed_policy`: `fixed`, `varying`

No `bwmethod` is used.

## DGP

```r
x <- runif(n)
z <- runif(n, min = -2, max = 2)
y <- x * exp(z) * (1.0 + rnorm(n, sd = 0.2))
```

## One Case

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npscoef/bench_npscoef_param_nprmpi.R \
  --n=250 --times=50 --base_seed=42 --nmulti=1 --nslaves=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npscoef_one_mpi_raw.csv --out_summary=/tmp/npscoef_one_mpi_summary.csv
```

## Combo Run

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npscoef/run_npscoef_combos.R \
  --n=250 --times=50 --base_seed=42 --nmulti=1 --nslaves=1 --tag=myrun
```
