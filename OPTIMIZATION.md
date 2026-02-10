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
- Runs per case: `NPRMPI_BENCH_TIMES=30`
- Sizes: `n = 100, 200, 400`
- Metrics: mean and median (ms)
- Baseline CSV: `/tmp/nprmpi_bench_baseline.csv`
- Current CSV: `/tmp/nprmpi_bench_current.csv`

## Results
Speedups are reported as `baseline / current` (so `>1` is faster) with percent change in parentheses.

| Function | Variant | n | Mean (ms) Baseline | Mean (ms) Current | Mean Speedup | Median (ms) Baseline | Median (ms) Current | Median Speedup |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| npcdens | est | 100 | 173.060 | 182.944 | 0.946 (-5.7%) | 204.746 | 204.730 | 1.000 (0.0%) |
| npcdens | est | 200 | 163.451 | 180.623 | 0.905 (-10.5%) | 202.330 | 205.211 | 0.986 (-1.4%) |
| npcdens | est | 400 | 195.923 | 152.072 | 1.288 (+22.4%) | 207.032 | 176.655 | 1.172 (+14.7%) |
| npksum | eval=tx | 100 | 134.412 | 160.450 | 0.838 (-19.4%) | 104.197 | 197.872 | 0.527 (-89.9%) |
| npksum | eval=tx | 200 | 160.263 | 161.022 | 0.995 (-0.5%) | 197.840 | 200.249 | 0.988 (-1.2%) |
| npksum | eval=tx | 400 | 168.223 | 168.106 | 1.001 (+0.1%) | 203.068 | 199.821 | 1.016 (+1.6%) |
| npksum | ex50 | 100 | 144.214 | 157.802 | 0.914 (-9.4%) | 104.315 | 200.200 | 0.521 (-91.9%) |
| npksum | ex50 | 200 | 143.692 | 160.547 | 0.895 (-11.7%) | 104.260 | 197.908 | 0.527 (-89.8%) |
| npksum | ex50 | 400 | 174.698 | 181.937 | 0.960 (-4.1%) | 202.707 | 202.825 | 0.999 (-0.1%) |
| npreg | est | 100 | 193.400 | 183.349 | 1.055 (+5.2%) | 204.499 | 204.344 | 1.001 (+0.1%) |
| npreg | est | 200 | 186.291 | 185.583 | 1.004 (+0.4%) | 204.500 | 203.380 | 1.006 (+0.5%) |
| npreg | est | 400 | 170.122 | 183.153 | 0.929 (-7.7%) | 201.285 | 204.970 | 0.982 (-1.8%) |
| npudens | est | 100 | 165.546 | 175.907 | 0.941 (-6.3%) | 203.146 | 203.533 | 0.998 (-0.2%) |
| npudens | est | 200 | 155.286 | 159.563 | 0.973 (-2.8%) | 169.722 | 201.145 | 0.844 (-18.5%) |
| npudens | est | 400 | 169.646 | 163.187 | 1.040 (+3.8%) | 204.354 | 201.330 | 1.015 (+1.5%) |

## Notes
- Results are mixed; some cases improve, others regress. We should re-check the benchmark harness (MPI setup and synchronization) if consistency is required before merging conclusions into user-facing claims.
