# npplreg Benchmark Harness (np)

Parameterized microbenchmark harness for `npplregbw()` + `npplreg()` (serial `np`).

## Grid

`run_npplreg_combos.R` runs 8 combinations (`2x2x2`):

- `ckertype`: `gaussian`, `epanechnikov`
- `np_tree`: `TRUE`, `FALSE`
- `seed_policy`: `fixed`, `varying`

There is no `bwmethod` option for this harness.

## DGP

```r
x1 <- rnorm(n)
x2 <- factor(rbinom(n, 1, .5))
z1 <- factor(rbinom(n, 1, .5))
z2 <- rnorm(n)
y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)
```

Bandwidth call:

```r
npplregbw(y ~ x1 + x2 | z1 + z2, ...)
```

## One Case

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npplreg/bench_npplreg_param.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 \
  --ckertype=gaussian --np_tree=FALSE --seed_policy=fixed \
  --out_raw=/tmp/npplreg_one_raw.csv --out_summary=/tmp/npplreg_one_summary.csv
```

## Combo Run

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/methods/npplreg/run_npplreg_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --tag=myrun
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
