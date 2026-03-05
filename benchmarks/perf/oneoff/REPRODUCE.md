# Reproduce: one-off benchmarks (np)

Example synthetic run (`npunitest`):

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/oneoff/bench_oneoff_param.R \
  --fun=npunitest --n=1000 --times=5 --seed_policy=varying --base_seed=42 \
  --out_raw=/tmp/np_oneoff_npunitest_raw.csv \
  --out_summary=/tmp/np_oneoff_npunitest_summary.csv
```

Example real-data run (`npqreg`, ignores `--n`):

```bash
Rscript /Users/jracine/Development/np-master/benchmarks/perf/oneoff/bench_oneoff_param.R \
  --fun=npqreg --n=1000 --times=5 --seed_policy=fixed --base_seed=42 \
  --out_raw=/tmp/np_oneoff_npqreg_raw.csv \
  --out_summary=/tmp/np_oneoff_npqreg_summary.csv
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
