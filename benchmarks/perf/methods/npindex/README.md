# npindex Benchmark Harness (npRmpi)

Parameterized microbenchmark harness for `npindexbw()` + `npindex()` (`npRmpi`).

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

`method="kleinspady"` (binary response stays integer):

```r
x1 <- runif(n, min = -1, max = 1)
x2 <- runif(n, min = -1, max = 1)
y <- ifelse(x1 + x2 + rnorm(n) > 0, 1L, 0L)
```

## One Case

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npindex/bench_npindex_param_nprmpi.R \
  --n=100 --times=50 --base_seed=42 --method=ichimura --nmulti=1 --nslaves=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npindex_one_mpi_raw.csv --out_summary=/tmp/npindex_one_mpi_summary.csv
```

## Combo Run

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npindex/run_npindex_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --nslaves=1 --tag=myrun
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
