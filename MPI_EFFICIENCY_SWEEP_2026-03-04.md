# MPI Efficiency Sweep (`npRmpi`) - 2026-03-04

## Objective
Identify any remaining paths that pull expensive work off MPI or off efficient hot-path routing.

## Scope
1. Static scan of `src/` and `R/` for serializing guards and suppress-parallel toggles.
2. Targeted runtime probe of plot/bootstrap wild path with guard toggles.
3. Focus on current `np-npRmpi` head after LL/LP CVLS gate restoration.

## Artifacts
1. Static sweep root: `/tmp/mpi_efficiency_sweep_20260304_1`
2. Subprocess harness (guard on): `/tmp/mpi_efficiency_guard_true_20260304_1`
3. Subprocess harness (guard off): `/tmp/mpi_efficiency_guard_false_20260304_1`
4. Stress probe logs:
- `/tmp/mpi_eff_true_n300.log`
- `/tmp/mpi_eff_false_n300.log`
- `/tmp/mpi_eff_true_n1000_retry.log`
- `/tmp/mpi_eff_false_n1000_retry.log`
5. Dedicated stress harness roots (`npplot_wild_stress`):
- `/tmp/mpi_efficiency_stress_guard_true_20260304_1`
- `/tmp/mpi_efficiency_stress_guard_false_20260304_1`

## Static Findings
1. `CVLS_FORCED_GATE_LINES=0` for pattern
   - `if((bwm == RBWM_CVLS) || ks_tree_use || (BANDWIDTH_reg == BW_ADAP_NN))`
   indicating no reintroduced LL/LP CVLS forced-gate regression in `src/jksum.c`.
2. No additional estimator/CV-family hot-path callsites were found forcing `suppress_parallel=1` in a way that bypasses outer-rank partitioning.
3. Remaining explicit off-MPI controls are concentrated in plot/bootstrap helpers:
   - `R/np.plot.helpers.R`: `npRmpi.plot.wild.master_local.guard` forces master-local wild bootstrap dispatch when enabled.
   - `R/np.plot.helpers.R`: `.np_plot_kernel_weights_direct()` sets `suppress.parallel = TRUE` for direct `C_np_kernelsum` kernel-weight extraction.

## Runtime Probe Findings
1. Subprocess harness `npplot` job passed with guard both on and off for the small diagnostic job (`plot.errors.boot.num=29`).
2. Medium stress (`n=250/300`, `B=799`) passed with guard on/off; guard-off was slower on this host/config.
3. Larger stress (`n=1000`, `B=799`) aborted (`Abort trap: 6`) with guard both on and off; log files remained empty.
4. Dedicated subprocess stress job (`npplot_wild_stress`) failed with `EXIT=134` for guard both on and off, confirming guard state is not the primary root cause of this failure class.

## Interpretation
1. Core estimator/CV MPI routing is not currently showing additional static evidence of forced serial detours beyond the already-fixed CVLS gate issue.
2. The dominant remaining MPI-efficiency debt is in plot/bootstrap infrastructure, not estimator/CV kernels.
3. The large-`n` abort is not explained solely by the wild master-local guard; there is a deeper plot/bootstrap runtime fragility at higher workload.

## Priority Decision
Highest priority now is **build-stage/runtime hardening**, not core SPMD remediation:
1. Keep core SPMD remediation closed for estimator/CV paths.
2. Continue build-stage hardening to isolate and fix high-workload plot/bootstrap aborts.
3. Defer guard removal until transport/runtime stability is proven under stress workloads.

## Immediate Next Tranche (Recommended)
1. Add targeted diagnostics around bootstrap phases (`dispatch`, `collect`, `done`) and worker return transport for large-`n` wild path.
2. Add dedicated stress job to subprocess harness (`npplot_wild_stress`) with reproducible parameters and timeout.
3. Gate any guard relaxation on:
   - no abort/hang in stress jobs,
   - no orphan workers,
   - no numerical drift for bootstrap outputs.
