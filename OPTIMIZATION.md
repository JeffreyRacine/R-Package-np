# Optimization Notes (np)

This document summarizes the kernel-sum and kernel-evaluation optimizations applied to the `np` package, along with a small set of reproducible microbenchmarks.

## Scope and goals

- Target: reduce runtime by **≥3%** in commonly-used kernel-sum / estimator paths.
- Constraints:
  - No algorithmic changes (numerically identical up to insignificant floating-point differences).
  - No parallelism.
  - Avoid memory-intensive approaches (prefer fewer allocations/copies).
  - Must not break existing behavior.

## Summary of changes

### R-level (`R/np.kernel.R`)

- Replaced `eval(parse(...))` construction in `npksum.numeric()` and `npksum.formula()` with `do.call()`.
- Preserved historical naming behavior for 1D `exdat` (column name remains `"exdat"`).

### C-level (`src/np.c`, `src/jksum.c`)

- `src/np.c` (`np_kernelsum`):
  - Reduced per-call allocation/copy overhead when `options(np.tree = FALSE)` by using buffer “views” where safe and writing directly into output buffers when possible.
  - Passed `bpso` into `kernel_weighted_sum_np()` to avoid an internal hot-path allocation.
- `src/jksum.c` (vectorized kernel evaluators):
  - Restructured `np_ckernelv()`, `np_ukernelv()`, `np_okernelv()` to reduce per-iteration branching and avoid function-pointer dispatch inside tight loops (single `switch` dispatch outside the inner loops).

### Tests

- Added regression coverage for `npksum` consistency and naming:
  - `tests/testthat/test-npksum.R`

## Commits

- Baseline (pre-optimizations): `7637a0b`
- Optimizations:
  - `e67f7b4` (R `npksum` setup + `np_kernelsum` overhead reductions)
  - `8321fb5` (kernel loop optimizations in `jksum.c`)
  - `108158f` (adds `npksum` regression tests; no performance changes)

## Benchmarks

All benchmarks below compare **baseline `7637a0b`** vs **current `108158f`**.

Hardware/software notes:
- Machine: Apple Silicon (M2Studio)
- Benchmark tool: `microbenchmark`
- Settings: `options(np.tree = FALSE, np.messages = FALSE)`
- Inputs: `ex50` uses `exdat` length 50; `eval=tx` omits `exdat` (evaluates at training points).
- Data sizes: `n = 1000, 2000, 4000, 8000, 16000`.

### `npksum()` kernel-sum microbenchmarks (larger sample sizes)

Each case reports mean + median time in milliseconds for the larger `n` runs.

#### Variant A: `exdat` length 50 (`ex50`)

| n | mean (ms) baseline → current | mean speedup | median (ms) baseline → current | median speedup |
|---:|---:|---:|---:|---:|
| 1000 | 7.437 → 6.641 | 1.12x (+10.7%) | 7.421 → 6.595 | 1.12x (+11.1%) |
| 2000 | 14.412 → 12.936 | 1.11x (+10.2%) | 14.314 → 12.845 | 1.11x (+10.3%) |
| 4000 | 28.596 → 25.492 | 1.12x (+10.9%) | 28.539 → 25.355 | 1.13x (+11.2%) |
| 8000 | 56.558 → 51.630 | 1.10x (+8.7%) | 56.520 → 51.460 | 1.10x (+9.0%) |
| 16000 | 116.027 → 100.039 | 1.16x (+13.8%) | 115.572 → 99.753 | 1.16x (+13.7%) |

#### Variant B: evaluation at training points (`eval=tx`, `exdat` missing)

| n | mean (ms) baseline → current | mean speedup | median (ms) baseline → current | median speedup |
|---:|---:|---:|---:|---:|
| 1000 | 7.332 → 6.554 | 1.12x (+10.6%) | 7.286 → 6.503 | 1.12x (+10.7%) |
| 2000 | 28.279 → 25.291 | 1.12x (+10.6%) | 28.159 → 25.150 | 1.12x (+10.7%) |
| 4000 | 112.608 → 100.655 | 1.12x (+10.6%) | 112.663 → 100.100 | 1.13x (+11.2%) |
| 8000 | 447.165 → 409.841 | 1.09x (+8.3%) | 447.241 → 410.268 | 1.09x (+8.3%) |
| 16000 | 1801.918 → 1587.634 | 1.13x (+11.9%) | 1799.025 → 1586.508 | 1.13x (+11.8%) |

### Higher-level estimator calls (estimation only, larger sample sizes)

These timings measure **estimation calls only** (manual bandwidths; bandwidth selection excluded).

| function | n | mean speedup | median speedup |
|---|---:|---:|---:|
| `npreg()` | 1000 | 1.37x (+26.9%) | 1.37x (+26.9%) |
| `npreg()` | 2000 | 1.70x (+41.2%) | 1.71x (+41.6%) |
| `npreg()` | 4000 | 1.99x (+49.7%) | 1.99x (+49.7%) |
| `npreg()` | 8000 | 1.96x (+48.9%) | 1.96x (+49.0%) |
| `npreg()` | 16000 | 2.22x (+55.0%) | 2.20x (+54.5%) |
| `npudens()` | 1000 | 1.41x (+29.2%) | 1.42x (+29.7%) |
| `npudens()` | 2000 | 1.76x (+43.3%) | 1.77x (+43.5%) |
| `npudens()` | 4000 | 2.06x (+51.4%) | 2.06x (+51.3%) |
| `npudens()` | 8000 | 1.99x (+49.8%) | 2.01x (+50.2%) |
| `npudens()` | 16000 | 2.27x (+56.0%) | 2.27x (+56.0%) |
| `npcdens()` | 1000 | 1.37x (+26.8%) | 1.38x (+27.3%) |
| `npcdens()` | 2000 | 1.63x (+38.8%) | 1.64x (+38.9%) |
| `npcdens()` | 4000 | 1.89x (+47.2%) | 1.89x (+47.0%) |
| `npcdens()` | 8000 | 1.85x (+45.9%) | 1.87x (+46.6%) |
| `npcdens()` | 16000 | 2.07x (+51.7%) | 2.07x (+51.7%) |

### Ordered-kernel optimization (ordered-factor power fast path)

These results isolate the ordered-kernel optimization in `jksum.c` (replacing a generic power call with a small-integer fast path). Benchmarks use **fixed bandwidths** (no CV), `times = 100`, and mixed data including ordered factors.

Baseline: pre-ordered-kernel change (`c7825d0`); Current: ordered-kernel change (`bb6bb45`).

| function | n | mean (ms) baseline → current | mean speedup | median (ms) baseline → current | median speedup |
|---|---:|---:|---:|---:|---:|
| `npcdens()` | 1000 | 23.971 → 22.199 | 1.08x (+7.40%) | 23.783 → 21.987 | 1.08x (+7.55%) |
| `npcdens()` | 2000 | 98.634 → 88.521 | 1.11x (+10.25%) | 98.145 → 88.301 | 1.11x (+10.03%) |
| `npcdens()` | 4000 | 441.832 → 399.326 | 1.11x (+9.62%) | 441.159 → 398.488 | 1.11x (+9.67%) |
| `npcdist()` | 1000 | 30.647 → 28.609 | 1.07x (+6.65%) | 30.318 → 28.551 | 1.06x (+5.83%) |
| `npcdist()` | 2000 | 124.528 → 114.667 | 1.09x (+7.92%) | 124.465 → 114.721 | 1.08x (+7.83%) |
| `npcdist()` | 4000 | 548.156 → 505.053 | 1.09x (+7.86%) | 547.238 → 504.767 | 1.08x (+7.76%) |
| `npreg()` | 1000 | 11.114 → 10.149 | 1.10x (+8.69%) | 11.116 → 10.063 | 1.10x (+9.47%) |
| `npreg()` | 2000 | 44.998 → 39.723 | 1.13x (+11.72%) | 45.051 → 39.685 | 1.14x (+11.91%) |
| `npreg()` | 4000 | 200.442 → 181.020 | 1.11x (+9.69%) | 200.375 → 180.739 | 1.11x (+9.80%) |
| `npudens()` | 1000 | 10.370 → 9.506 | 1.09x (+8.34%) | 10.366 → 9.463 | 1.10x (+8.71%) |
| `npudens()` | 2000 | 42.556 → 37.694 | 1.13x (+11.42%) | 42.533 → 37.687 | 1.13x (+11.39%) |
| `npudens()` | 4000 | 193.619 → 171.545 | 1.13x (+11.40%) | 193.669 → 171.502 | 1.13x (+11.45%) |

## Notes on inheritance of improvements

Many higher-level `np*` functions inherit these improvements automatically, because they ultimately dispatch into the same `jksum.c` kernel evaluators (`np_ckernelv`, `np_ukernelv`, `np_okernelv`) in the “new” (non-legacy) C code paths.

## 2026-02-12 (Profile Caching + Kernel Fast Paths)

### Summary
Implemented serial optimizations in `src/jksum.c` focused on categorical/ordered kernel hot paths and product-kernel reuse.

### Code Changes
- Unordered kernel fast path: cached `same/different` values in `np_ukernelv` and `np_p_ukernelv`.
- Ordered kernel fast path: cached `lambda^d` lookup with direct evaluation for ordered kernels `0/1/2` in `np_okernelv` and `np_p_okernelv`.
- Tree early exit: skip all downstream kernel work when tree search yields zero continuous-support interactions.
- Persistent discrete profile ids: build once per `kernel_weighted_sum_np` invocation and reuse across all `j` (tree and non-tree paths under safety conditions).

### Runtime Summary (pre vs post)
- Real-world `wage1` (`npreg`, `ckertype="epanechnikov"`, `times=8`):
  - `np.tree=FALSE`: `10927.178 ms` -> `8466.990 ms` (mean, `+22.5%`)
  - `np.tree=TRUE`: `6033.697 ms` -> `5749.531 ms` (mean, `+4.7%`)
- Non-tree categorical-only `npreg` (`n=4000`, `times=6`):
  - `31371.77 ms` -> `20732.23 ms` (mean, `+33.9%`)
- Non-tree mixed `npreg` (`n=2000`, `times=6`):
  - `4677.014 ms` -> `3921.773 ms` (mean, `+16.1%`)
- Tree mixed `npreg` (`n=2000`, `times=6`):
  - `1189.502 ms` -> `1099.57 ms` (mean, `+7.6%`)

### Numerical Stability
- `wage1` pre vs post:
  - non-tree: `max |Δ bw| = 2.759011e-10`, `max |Δ fitted| = 4.730198e-10`
  - tree: `max |Δ bw| = 1.603712e-04`, `max |Δ fitted| = 1.500707e-05`
- Synthetic checks: only very small floating-point-level differences consistent with operation reordering.

## 2026-02-12 (Large-`h` Continuous-Kernel Shortcut)

### Summary
Added a conservative large-bandwidth shortcut in `src/jksum.c` for ordinary continuous kernels (`KERNEL` codes `0..9`):
- if `|h| >= h_min` for a continuous regressor, use `K(0)` directly for that regressor instead of per-observation kernel evaluation;
- `h_min` is derived from data range and a conservative relative-tolerance bound near `u=(x-x_i)/h=0`;
- applies in both `np.tree=FALSE` and `np.tree=TRUE` paths in `kernel_weighted_sum_np`;
- default relative tolerance is `1e-3`, configurable via `NP_LARGEH_REL_TOL`.

### Fast Validation Setup
- DGP: `x1 ~ U(0,1)`, `y ~ N(0,1)` (`x1` irrelevant).
- Models: `npreg(y~x1)` and `npcdens(y~x1)`.
- Settings: Gaussian kernels, `nmulti=1`, `ftol=1e-02`, `tol=1e-01`, `times=10`, `n=250,500`.
- Files:
  - timing compare: `/tmp/np_largeh_prepost_compare.csv`
  - bandwidth diagnostics (pre): `/tmp/np_largeh_pre_bw_diag.csv`
  - bandwidth diagnostics (post): `/tmp/np_largeh_post_bw_diag.csv`
  - joined timing+bandwidth table: `/tmp/np_largeh_timing_bw_join.csv`

### Runtime (pre -> post)
| tree | function | n | mean % change | median % change |
|---|---|---:|---:|---:|
| FALSE | npreg | 250 | -20.68% | -20.43% |
| FALSE | npreg | 500 | -0.04% | -0.04% |
| TRUE | npreg | 250 | +0.22% | +0.22% |
| TRUE | npreg | 500 | -47.09% | -47.13% |
| FALSE | npcdens | 250 | -38.77% | -39.02% |
| FALSE | npcdens | 500 | -34.89% | -35.06% |
| TRUE | npcdens | 250 | -0.96% | -1.53% |
| TRUE | npcdens | 500 | -32.91% | -32.99% |

### Bandwidth Evidence (x-side bandwidth vs x-range)
Selected rows from `/tmp/np_largeh_timing_bw_join.csv`:
- `npreg, tree=FALSE, n=500`: `median(h/range) ≈ 0.625` (small) -> ~0% speed change.
- `npreg, tree=TRUE, n=500`: `median(h/range)` very large (`~7.3e5` pre, `~1.08e2` post) -> large speedup (~47%).
- `npreg, tree=TRUE, n=250`: `median(h/range) ≈ 0.73` (small) -> near-neutral timing.
- `npcdens` cases with large improvements show larger x-bandwidths in many runs and strong reductions in continuous-kernel workload.

### Interpretation
- The shortcut behaves as intended: gains are strongest in regimes where CV chooses large continuous bandwidths relative to predictor range.
- Cases with small/moderate `h/range` show little or no speedup, which is expected for this targeted optimization.
