# Build-Stage MPI Hardening Plan (`npRmpi`) - 2026-03-04

## Objective
Stabilize build/check-stage MPI validation without masking runtime bugs.

## Scope
In scope:
1. Check/build orchestration hardening.
2. Dedicated subprocess MPI integration harness (one fresh process per job).
3. Legacy-path diagnostics for worker-return protocol fragility.

Out of scope for this tranche:
1. Statistical method changes.
2. Broad MPI runtime rewrites across all legacy helpers.

## Root-Cause Framing
1. `R CMD check --as-cran` forced MPI examples (`NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1`) run many MPI-heavy examples in one long process and expose lifecycle/protocol fragility.
2. This is not primarily shell sequencing; sequencing only affects cleanup hygiene.
3. Existing SPMD opcode ACK/timeout hardening is strong; remaining fragility is concentrated in legacy worker-return paths used by examples.

## Risk Model
Low risk:
1. Trap-safe check orchestration script.
2. Subprocess MPI example/integration harness.

Medium/high risk:
1. Refactoring legacy worker return transport in `R/Rparutilities.R`.

## Hardening Strategy
### Phase A (execute now, low risk)
1. Add trap-safe check runner that always restores `man/dontrun` on failure.
2. Add subprocess MPI integration harness:
   - each job in a fresh `Rscript` process,
   - per-job timeout,
   - per-job logs + status,
   - orphan cleanup checks.
3. Keep `R CMD check --as-cran` default mode as release gate; forced mode only for diagnostics.

### Phase B (next tranche, medium risk)
1. Add diagnostics around legacy return-path shape/type negotiation in `mpi.remote.exec` / `.mpi.worker.exec`.
2. Capture structured failure signatures under subprocess harness.

### Phase C (later, higher risk)
1. Incrementally replace legacy mixed typed gather/gatherv return protocol with single envelope return protocol in targeted paths.
2. Gate each step by session/attach/profile/manual route checks + orphan cleanup + numerical parity.

## Acceptance Criteria (Phase A)
1. Trap-safe check runner restores docs state even on failure.
2. MPI subprocess harness runs multiple jobs and emits machine-readable pass/fail summary.
3. Scripts are `shellcheck`-clean.
4. No runtime/API behavior change in package code.

## Artifacts
1. Script: `issue_notes/run_mpi_examples_subprocess.sh`
2. Script: `issue_notes/run_check_with_man_trap.sh`
3. Validation bundle: `/tmp/build_stage_hardening_phaseA_<timestamp>/`

## Phase A Execution Status (2026-03-04)
Completed:
1. Added trap-safe check runner (`run_check_with_man_trap.sh`) with guaranteed `man/dontrun` cleanup via `trap`.
2. Added subprocess MPI harness (`run_mpi_examples_subprocess.sh`) with:
   - per-job fresh process,
   - per-job timeout,
   - machine-readable summary outputs,
   - orphan-process detection.
3. Hardened harness isolation:
   - fixed exit-code capture bug,
   - isolated each job in a new process group/session (`POSIX::setsid`) to prevent MPI teardown from killing the harness shell.
4. Validation completed:
   - `shellcheck` clean for hardening scripts,
   - trap-safe check runner pass (`CHECK_ARGS='--no-examples --no-tests'`),
   - route sanity smokes pass: session (`nslaves=1`) and attach (`mpiexec -n 2`) for `npreg`.

Observed diagnostics:
1. `npplot` subprocess job still aborts (`EXIT=134`) in forced diagnostic mode; this confirms known runtime fragility remains in that path and is now isolated/reproducible per-job.

Validation artifact roots:
1. `/tmp/build_stage_hardening_tty2_20260304_113247` (PASS: `npreg`, `npcondensitybw`)
2. `/tmp/build_stage_hardening_tty2_20260304_113255` (FAIL: `npplot`, `EXIT=134`)
3. `/tmp/build_stage_hardening_tty2_20260304_113214` (PASS: `npreg`)

## Phase B Incremental Hardening (2026-03-04)
Implemented:
1. Removed forced `gc()` calls from R object transport hot paths:
   - `mpi.bcast.Robj()`
   - `mpi.send.Robj()`
2. Hardened bootstrap fan-out worker serialization by rebinding worker environments to a minimal namespace-backed frame plus declared bindings (`.npRmpi_bootstrap_prepare_worker`).
3. Added targeted runtime guard for known unstable configuration:
   - when `what == "wild"` and `ncol.out <= 3` with active workers, execute that slice on master only with an explicit warning (`npRmpi.plot.wild.master_local.guard`, default `TRUE`).

Validation:
1. Prior crash repro (`npplot` mixed continuous+ordered, wild bootstrap) now completes in patched build.
2. Subprocess harness full set passes against patched library:
   - `/tmp/build_stage_hardening_tty3_20260304_115011` (`npreg`, `npplot`, `npcondensitybw` all PASS)
3. Targeted tests pass:
   - `tests/testthat/test-plot-autodispatch.R`
   - `tests/testthat/test-plot-mpi-only-bootstrap-contract.R`

Follow-up hardening:
1. Added non-interactive progress default-off (`np.plot.progress.noninteractive` opt-in) to avoid check-run console churn.
2. Added robust graphics parameter restore helper (`.np_plot_restore_par`) and wired plot engines to use it.
3. Added conservative MPI-session wild chunk cap (`np.plot.wild.chunk.max.mpi`, default `64`) for stability under repeated bootstrap slices.

Residual blocker:
1. Full forced all-examples `R CMD check` can still abort/hang under strict check internals; current stop-point is in the `np.pairs` example stream (intermittent `MPI_Allgather` truncation / `npudens` ACK mismatch signatures), not the earlier wild-bootstrap `np.plot` crash.
2. Fresh-install subprocess harness remains green for the targeted runtime jobs (`npreg`, `npplot`, `npcondensitybw`), confirming runtime hardening gains while the long monolithic check-stream fragility remains isolated.
3. Next tranche should be docs/example-structure focused for check mode (keep examples lightweight/deterministic in `R CMD check`; keep full MPI lifecycle stress in subprocess integration harnesses).

## Phase B.1 Efficiency Sweep Addendum (2026-03-04)
Static + runtime sweep artifacts:
1. Static sweep: `/tmp/mpi_efficiency_sweep_20260304_1`
2. Guard-on subprocess `npplot`: `/tmp/mpi_efficiency_guard_true_20260304_1`
3. Guard-off subprocess `npplot`: `/tmp/mpi_efficiency_guard_false_20260304_1`
4. Stress probes documented in:
   - `/Users/jracine/Development/np-npRmpi/MPI_EFFICIENCY_SWEEP_2026-03-04.md`

Findings:
1. No reintroduced forced CVLS gate pattern in `src/jksum.c` (`CVLS_FORCED_GATE_LINES=0`).
2. Remaining explicit off-MPI controls are concentrated in plot/bootstrap:
   - `npRmpi.plot.wild.master_local.guard`,
   - `.np_plot_kernel_weights_direct()` with `suppress.parallel = TRUE`.
3. Larger stress (`n=1000`, `plot.errors.boot.num=799`) aborted with guard both on and off (`Abort trap: 6`), indicating a deeper large-workload plot/bootstrap fragility beyond the guard itself.
4. Dedicated subprocess stress job `npplot_wild_stress` is now available in `run_mpi_examples_subprocess.sh` and currently fails in both modes:
   - guard on: `/tmp/mpi_efficiency_stress_guard_true_20260304_1` (`EXIT=134`)
   - guard off: `/tmp/mpi_efficiency_stress_guard_false_20260304_1` (`EXIT=134`)

Priority update:
1. Build-stage/runtime hardening remains the highest-priority active track.
2. Core estimator/CV SPMD remediation remains closed unless new core-family evidence appears.

## Explicit TODO (Deferred Until Current Regression Is Fixed)
1. Restore distributed MPI wild-bootstrap fan-out for `np.plot` paths (remove temporary master-local guard once transport stability is proven).
2. Required acceptance before re-enabling:
   - no crash/hang in targeted wild-bootstrap subprocess soak runs,
   - route sanity pass in session and attach modes,
   - no numerical drift relative to current guarded behavior (within declared tolerance),
   - no material runtime regression vs intended parallel fan-out baseline.
3. Use `npplot_wild_stress` as a mandatory gate for any wild-path transport or guard change; do not relax guard policy until this job is stable with `EXIT=0`.
4. Add phase-level diagnostics in worker return transport (`mpi.remote.exec` / `.mpi.worker.exec`) and bootstrap fanout (`.npRmpi_bootstrap_run_fanout`) to pinpoint `EXIT=134` origin under stress.

## Phase B.2 Forensic Breadcrumb Hardening (2026-03-04)
Implemented:
1. Added opt-in persistent bootstrap phase breadcrumbs in `R/np.plot.helpers.R`:
   - new trace path resolver (`option npRmpi.bootstrap.phase.file` or env `NP_RMPI_BOOTSTRAP_PHASE_FILE`),
   - phase append hook wired into `.npRmpi_bootstrap_phase_mark(...)`.
2. Extended subprocess harness `issue_notes/run_mpi_examples_subprocess.sh`:
   - `BOOTSTRAP_PHASE_TRACE=1` emits per-job phase files (`*.phase.tsv`),
   - stage breadcrumbs (`STAGE library/init/plot`) added to generated R jobs for coarse crash localization.
3. Extended static sweep `issue_notes/run_mpi_efficiency_sweep.sh`:
   - added non-plot scan for `suppress.parallel = TRUE` to catch accidental off-MPI drift outside plot helpers.

Validation:
1. Static sweep rerun:
   - `/tmp/mpi_efficiency_sweep_20260304_2`
   - confirms `SUPPRESS_TRUE_NONPLOT_R_LINES=0` and `CVLS_FORCED_GATE_LINES=0`.
2. Traced subprocess pass (`npplot`):
   - `/tmp/build_stage_hardening_phaseB2_npplot_src_20260304_2` (`PASS`, `PHASE_LINES=16`).
3. Traced stress repro (`npplot_wild_stress`) still fails guard on/off:
   - guard on: `/tmp/build_stage_hardening_phaseB2_wildstress_guardtrue_src_20260304_2` (`EXIT=134`, `PHASE_LINES=12`)
   - guard off: `/tmp/build_stage_hardening_phaseB2_wildstress_guardfalse_src_20260304_2` (`EXIT=134`, `PHASE_LINES=12`)
4. Phase tails in both failing stress runs end at repeated `phase=done` records, indicating fanout dispatch/collect completed for multiple chunks before abort.
5. Route validation on patched build:
   - manual validator pass,
   - attach validator pass,
   - profile validator pass,
   - attach plot smoke pass (`ATTACH_PLOT_OK`).

Interpretation:
1. The new breadcrumbs narrow failure location: aborts are not occurring at fanout entry/dispatch/collect mismatch; they occur after multiple completed fanout cycles.
2. Guard state does not change failure class at stress scale (`EXIT=134` both modes).
3. Next highest-value hardening remains transport/lifecycle diagnostics in legacy return paths, now with stronger phase evidence.

## Phase B.3 Transport Trace Deepening (2026-03-04)
Implemented:
1. Added opt-in transport trace hooks in legacy transport paths (`R/Rparutilities.R`):
   - `mpi.remote.exec` master-side start/type/path/done breadcrumbs,
   - `.mpi.worker.exec` worker-side start/type/send/done breadcrumbs,
   - `mpi.applyLB` master start/done breadcrumbs,
   - `.mpi.worker.applyLB` worker start/malformed/done breadcrumbs.
2. Added plot-bootstrap transport breadcrumbs in `R/np.plot.helpers.R` fanout loop:
   - `fanout.start`, `fanout.master_local.start/done`,
   - `fanout.master_assist.start`, per-message `send/recv`, `fanout.master_assist.done`,
   - terminal `fanout.done` and `fanout.error`.
3. Extended subprocess harness to emit optional per-job transport logs:
   - `TRANSPORT_TRACE=1` writes `*.transport.tsv`,
   - existing phase logs (`*.phase.tsv`) retained.
4. Added trace-hook coverage tests:
   - `tests/testthat/test-transport-trace-hook.R`,
   - extended `tests/testthat/test-plot-bootstrap-phase-trace.R` for bootstrap transport trace.

Validation:
1. Targeted tests pass on patched install:
   - `test-plot-bootstrap-phase-trace.R`,
   - `test-transport-trace-hook.R`,
   - `test-plot-mpi-only-bootstrap-contract.R`.
2. Route validators pass on patched install:
   - `MANUAL_BCAST_ROUTE_OK`,
   - `ATTACH_ROUTE_OK`,
   - `PROFILE_ROUTE_OK`,
   - attach plot smoke `ATTACH_PLOT_OK`.
3. Traced subprocess runs:
   - `/tmp/build_stage_hardening_phaseB3_npplot_src_20260304_2` (`PASS`, `PHASE_LINES=16`, `TRANSPORT_LINES=16`)
   - `/tmp/build_stage_hardening_phaseB3_wildstress_guardtrue_src_20260304_2` (`PASS`, `PHASE_LINES=16`, `TRANSPORT_LINES=16`)
   - `/tmp/build_stage_hardening_phaseB3_wildstress_guardfalse_src_20260304_2` (`FAIL EXIT=134`, `PHASE_LINES=12`, `TRANSPORT_LINES=57`)

Interpretation:
1. New transport logs show complete master-assist send/recv cycles (`fanout.master_assist.done` + `fanout.done`) before the abort in failing runs.
2. Current failure signature remains post-cycle/intermittent (not a direct collectives-cadence divergence).
3. Immediate next hardening should target post-cycle memory/lifecycle pressure points in repeated wild slices (allocation churn, object lifetime, and worker loop teardown timing), using new transport+phase traces as gates.

## Phase B.4 Wild Stress Isolation + Rejected Fixes (2026-03-04)
Repro matrix updates (clean `HEAD`, fresh temp install):
1. Harness stress (`npplot_wild_stress`, `n=1000`) with trace:
   - guard on: `FAIL`, `EXIT=134` (3/3)
   - guard off: `FAIL`, `EXIT=134` (3/3)
   - artifacts: `/tmp/phaseD1_eval_head_20260304_143245`
2. Direct repeated runs (`Rscript`, session mode, `n=1000`, wild):
   - unstable, failing at `STAGE plot` (mixed `rc=6/137/0`)
   - artifacts:
     - `/tmp/direct_repeat_stage_matrix_20260304_143654.tsv`
     - logs: `/tmp/direct_repeat_stage_*.log`
3. Direct repeated runs (`Rscript`, session mode, `n=1000`, inid):
   - stable (`5/5` success)
   - artifacts:
     - `/tmp/direct_repeat_inid_stage_matrix_20260304_143759.tsv`
4. Direct repeated runs (`Rscript`, session mode, `n=100`, wild):
   - stable (`5/5` success)
   - artifact:
     - `/tmp/direct_repeat_wildn100_stage_matrix_20260304_144003.tsv`

Interpretation:
1. The active instability is wild-bootstrap specific at larger workload (`n=1000`) and not a general `plot(..., plot.errors.method="bootstrap")` failure.
2. Guard on/off does not resolve the failure class under stress on current `HEAD`.
3. Failure boundary remains inside wild plot path (`STAGE plot`), consistent with prior post-cycle abort signatures.

Attempted fixes (rejected; not merged):
1. Factor-slice local-hat forcing in wild path:
   - changed failure class from abort to deterministic timeout/hang (`EXIT=142`), therefore rejected.
2. Wild worker prebinding of `H`/inputs:
   - no robust stability gain in repeated stress; rejected as inconclusive.
3. Wild worker GEMM-to-GEMV loop rewrite:
   - no stability gain (worse in sampled run); rejected.
4. All above runtime edits were reverted; branch runtime remains at clean `HEAD` behavior.

Next safe tranche:
1. Keep runtime unchanged until a deterministic keep-worthy fix is proven.
2. Added dedicated deterministic stress gate script for wild (`n=1000`, ordered-factor slice, repeated subprocess runs):
   - script: `issue_notes/run_npplot_wild_stress_repeat.sh`
   - sample artifact: `/tmp/npplot_wild_gate_20260304_1` (`PASS=1 FAIL=2`, failed runs stop at `STAGE plot`).
3. Make this gate mandatory for wild-path changes.
4. Investigate low-risk observability in wild slice internals (without behavior change) to localize abort source before transport/protocol refactor.
