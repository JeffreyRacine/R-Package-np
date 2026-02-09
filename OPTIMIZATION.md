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

### `npksum()` kernel-sum microbenchmarks

Each case uses `times=200` and reports mean + median time in milliseconds.

#### Variant A: `exdat` length 50 (`ex50`)

| n | mean (ms) baseline → current | mean speedup | median (ms) baseline → current | median speedup |
|---:|---:|---:|---:|---:|
| 100 | 0.422 → 0.390 | 1.08x (+8.2%) | 0.399 → 0.366 | 1.09x (+9.0%) |
| 200 | 0.468 → 0.429 | 1.09x (+9.1%) | 0.447 → 0.407 | 1.10x (+9.8%) |
| 400 | 0.552 → 0.499 | 1.11x (+10.6%) | 0.524 → 0.474 | 1.11x (+10.5%) |

#### Variant B: evaluation at training points (`eval=tx`, `exdat` missing)

| n | mean (ms) baseline → current | mean speedup | median (ms) baseline → current | median speedup |
|---:|---:|---:|---:|---:|
| 100 | 0.423 → 0.374 | 1.13x (+13.1%) | 0.394 → 0.355 | 1.11x (+11.0%) |
| 200 | 0.687 → 0.563 | 1.22x (+22.0%) | 0.626 → 0.547 | 1.14x (+14.4%) |
| 400 | 1.567 → 1.348 | 1.16x (+16.2%) | 1.528 → 1.335 | 1.14x (+14.5%) |

### Higher-level estimator calls (estimation only)

These timings measure **estimation calls only** (manual bandwidths; bandwidth selection excluded), using `times=100`:

| function | n | mean speedup | median speedup |
|---|---:|---:|---:|
| `npreg()` | 100 | 1.05x (+4.7%) | 1.05x (+4.6%) |
| `npreg()` | 200 | 1.07x (+6.9%) | 1.07x (+6.6%) |
| `npreg()` | 400 | 1.13x (+12.5%) | 1.14x (+13.7%) |
| `npudens()` | 100 | 1.04x (+4.2%) | 1.05x (+5.1%) |
| `npudens()` | 200 | 1.11x (+10.6%) | 1.09x (+9.4%) |
| `npudens()` | 400 | 1.17x (+17.4%) | 1.17x (+17.4%) |
| `npcdens()` | 100 | 1.04x (+4.1%) | 1.06x (+6.4%) |
| `npcdens()` | 200 | 1.10x (+10.5%) | 1.11x (+10.9%) |
| `npcdens()` | 400 | 1.09x (+9.0%) | 1.12x (+12.0%) |

## Notes on inheritance of improvements

Many higher-level `np*` functions inherit these improvements automatically, because they ultimately dispatch into the same `jksum.c` kernel evaluators (`np_ckernelv`, `np_ukernelv`, `np_okernelv`) in the “new” (non-legacy) C code paths.

