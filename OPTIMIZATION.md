# Optimization Summary (npRmpi)

## Scope
This summarizes the npRmpi performance work (ported from np) and the verification steps.

## Numerical Identity Verification
- Script: `/tmp/check_nprmpi_numeric.R`
- Mode: baseline write, current check
- Reference: `/tmp/nprmpi_numeric_ref.rds`
- Tolerance: `1e-12`
- Result: no differences detected

## Benchmark Method
- Script: `/tmp/bench_nprmpi_suite.R`
- Tool: `microbenchmark`
- Runs per case: `NPRMPI_BENCH_TIMES=100`
- Sizes: `n = 1000, 2000, 4000`
- Metrics: mean and median (ms)
- Baseline CSV: `/tmp/nprmpi_bench_baseline.csv`
- Current CSV: `/tmp/nprmpi_bench_current.csv`

## Results
Speedups are reported as `baseline / current` (so `>1` is faster) with percent change in parentheses.

| Function | Variant | n | Mean (ms) Baseline | Mean (ms) Current | Mean Speedup | Median (ms) Baseline | Median (ms) Current | Median Speedup |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| npcdens | est | 1000 | 207.362 | 201.775 | 1.028 (+2.7%) | 216.761 | 210.182 | 1.031 (+3.0%) |
| npcdens | est | 2000 | 244.587 | 240.566 | 1.017 (+1.6%) | 267.603 | 243.081 | 1.101 (+9.2%) |
| npcdens | est | 4000 | 480.371 | 355.959 | 1.350 (+25.9%) | 501.012 | 356.104 | 1.407 (+28.9%) |
| npksum | eval=tx | 1000 | 177.088 | 157.343 | 1.125 (+11.1%) | 202.148 | 200.713 | 1.007 (+0.7%) |
| npksum | eval=tx | 2000 | 196.564 | 167.965 | 1.170 (+14.5%) | 214.406 | 210.355 | 1.019 (+1.9%) |
| npksum | eval=tx | 4000 | 241.301 | 209.249 | 1.153 (+13.3%) | 259.370 | 250.150 | 1.037 (+3.6%) |
| npksum | ex50 | 1000 | 184.007 | 156.685 | 1.174 (+14.8%) | 202.627 | 200.653 | 1.010 (+1.0%) |
| npksum | ex50 | 2000 | 174.767 | 161.799 | 1.080 (+7.4%) | 206.176 | 204.010 | 1.011 (+1.1%) |
| npksum | ex50 | 4000 | 195.427 | 188.331 | 1.038 (+3.6%) | 214.164 | 211.325 | 1.013 (+1.3%) |
| npreg | est | 1000 | 200.368 | 197.160 | 1.016 (+1.6%) | 208.537 | 208.758 | 0.999 (-0.1%) |
| npreg | est | 2000 | 214.573 | 211.928 | 1.012 (+1.2%) | 231.588 | 217.606 | 1.064 (+6.0%) |
| npreg | est | 4000 | 311.843 | 247.301 | 1.261 (+20.7%) | 337.926 | 268.918 | 1.257 (+20.4%) |
| npudens | est | 1000 | 173.687 | 179.049 | 0.970 (-3.1%) | 205.960 | 204.871 | 1.005 (+0.5%) |
| npudens | est | 2000 | 216.519 | 192.817 | 1.123 (+10.9%) | 229.672 | 214.996 | 1.068 (+6.4%) |
| npudens | est | 4000 | 322.556 | 209.382 | 1.541 (+35.1%) | 335.685 | 168.338 | 1.994 (+49.9%) |

## Notes
- With larger sample sizes, the improvements are consistent for most cases, and the communication overhead is less dominant.

## 2026-02-12: jksum Port Verification (MPI, nslaves=1)
- Source change under test: `src/jksum.c` port from `np-master` (unordered/ordered kernel fast paths, cached ordered powers, tree empty-support short-circuit, persistent discrete profile cache).
- Benchmark driver: case-by-case MPI runs using documented `npRmpi.start()/mpi.bcast.cmd(..., caller.execute=TRUE)` pattern.
- Settings: synthetic mixed data, `n = 500, 1000, 2000`, `times = 10` for each case, `ftol = 1e-02`, `tol = 1e-01`.
- MPI transport note: `FI_TCP_IFACE=en0` used to avoid intermittent OFI timeout on `utun4` during heavy `npcdens` collectives.

### Pre/Post Results by Sample Size (averaged over `np.tree` on/off)
`rel = post / pre` (so `< 1.0` is faster post-change)

| Function | n | Mean Rel | Mean % Change | Median Rel | Median % Change |
|---|---:|---:|---:|---:|---:|
| npreg | 500 | 0.9153 | -8.47% | 0.8856 | -11.44% |
| npreg | 1000 | 0.8484 | -15.16% | 0.8715 | -12.85% |
| npreg | 2000 | 0.7744 | -22.56% | 0.7883 | -21.17% |
| npcdens | 500 | 0.9584 | -4.16% | 0.9535 | -4.65% |
| npcdens | 1000 | 0.9660 | -3.40% | 0.9663 | -3.37% |
| npcdens | 2000 | 0.9781 | -2.19% | 0.9821 | -1.79% |
| npudens | 500 | 0.9643 | -3.57% | 1.0228 | +2.28% |
| npudens | 1000 | 0.9920 | -0.80% | 0.9692 | -3.08% |
| npudens | 2000 | 0.9883 | -1.17% | 0.9924 | -0.76% |

### Brief Interpretation
- `npreg` shows strong scaling gains with sample size (largest at `n=2000`).
- `npcdens` shows consistent but smaller gains across sample sizes.
- `npudens` is modestly better on mean at all sample sizes; median is near-neutral at `n=500` and slightly improved at larger `n`.
- Across all 18 tree/function/sample cells: 16 mean wins and 15 median wins.

## 2026-02-18: All-Large Regression CV Fast Path Telemetry (Port Sync)

### Summary
- Ported/refined the shared `src/jksum.c` all-large regression CV path from `np-master`.
- Added cumulative bandwidth summary counters to mirror serial package behavior:
  - `num.feval`
  - `num.feval.fast`
  - `num.feval.fallback`
- Summary now prints: `Number of Function Evaluations: N (fast = X, fallback = Y)`.

### Notes
- Fast-path usage is conditional and may be sparse in nonlinear/small-bandwidth regimes.
- Fallback count tracks all-large attempts that reverted to generic evaluation.
- `npRmpi` behavior should be interpreted with MPI lifecycle/broadcast requirements intact.

### Release/Validation
- Current development baseline: `0.70-0`.
- Built and uploaded tarball to win-builder (`R-release` and `R-devel`) from:
  - `/Users/jracine/Development/npRmpi_0.70-0.tar.gz`
