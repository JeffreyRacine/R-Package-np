# npcdens Benchmark Harness (npRmpi)

This folder provides a parameterized benchmark harness for `npcdensbw()` + `npcdens()` in `npRmpi` using `microbenchmark`.

## Files

- `bench_npcdens_param_nprmpi.R`: run one benchmark configuration with `microbenchmark` repetitions.
- `run_npcdens_combos.R`: run all 16 option combinations (`2x2x2x2`) with a shared setup.
- `compare_npcdens_versions.R`: compare two raw outputs and emit overall and combo-specific comparison tables.
- `make_npcdens_report.R`: generate a markdown report from comparison CSV outputs.
- `make_npcdens_combined_report.R`: generate a single markdown report combining `np` and `npRmpi` comparison CSV outputs.
- `REPRODUCE.md`: end-to-end current-vs-CRAN replication steps.

## One-Configuration Run

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npcdens/bench_npcdens_param_nprmpi.R \
  --n=100 --times=50 --base_seed=42 \
  --bwmethod=cv.ml --nmulti=1 \
  --cxkertype=gaussian --cykertype=gaussian \
  --np_tree=FALSE --seed_policy=fixed \
  --rslaves=1 \
  --out_raw=/tmp/npcdens_one_mpi_raw.csv --out_summary=/tmp/npcdens_one_mpi_summary.csv
```

## Full Combination Run (16 combos)

```bash
R_LIBS=/tmp/Rlib_npRmpi_post FI_TCP_IFACE=en0 \
Rscript /Users/jracine/Development/np-npRmpi/benchmarks/perf/methods/npcdens/run_npcdens_combos.R \
  --n=100 --times=50 --base_seed=42 --nmulti=1 --rslaves=1 --tag=myrun
```

## Important Defaults

- `n=100`
- `times=50`
- `base_seed=42`
- `nmulti=1`
- `nslaves=1` (`--rslaves` alias supported)
- `bwmethod=cv.ml`
- `cxkertype=gaussian`
- `cykertype=gaussian`
- `np_tree=FALSE`
- `seed_policy=varying`

Notes:

- `bwmethod` choices used in combos are `cv.ml` and `cv.ls`.
- paired-kernel rule is enforced: `cxkertype == cykertype`, with pair choices `gaussian` or `epanechnikov`.

## DGP

The benchmark uses:

```r
set.seed(seed)
x <- runif(n)
z <- factor(rbinom(n, 3, .5))
y <- x + as.numeric(as.character(z)) + rnorm(n, sd = .25)
```

and fits `npcdensbw(y ~ x + z, ...)` then `npcdens(bws=bw, data=...)`.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
