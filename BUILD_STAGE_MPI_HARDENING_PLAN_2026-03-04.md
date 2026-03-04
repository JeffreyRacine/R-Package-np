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
