# npplreg Benchmark Harness (npRmpi)

Parameterized microbenchmark harness for `npplregbw()` + `npplreg()` (`npRmpi`).

## Grid

`run_npplreg_combos.R` runs 8 combinations (`2x2x2`):

- `ckertype`: `gaussian`, `epanechnikov`
- `np_tree`: `TRUE`, `FALSE`
- `seed_policy`: `fixed`, `varying`

No `bwmethod` is used.

## DGP

```r
x1 <- rnorm(n)
x2 <- factor(rbinom(n, 1, .5))
z1 <- factor(rbinom(n, 1, .5))
z2 <- rnorm(n)
y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)
```

## One Case

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npplreg/bench_npplreg_param_nprmpi.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --nslaves=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npplreg_one_mpi_raw.csv --out_summary=/tmp/npplreg_one_mpi_summary.csv
```

## Combo Run

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/npplreg/run_npplreg_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --nslaves=1 --tag=myrun
```
