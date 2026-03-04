# npRmpi CVLS Remediation Plan: Non-Helper Parity with Oracle (2026-03-04)

## Goal
Make the LL/LP non-helper CVLS path (the hoisted/tree/gather-scatter core path) numerically match the CVLS oracle behavior while preserving SPMD safety and keeping performance-oriented structure.

## Decision
1. Keep the non-helper path as the long-term production path for LL/LP CVLS.
2. Use the stable helper path as temporary oracle during remediation.
3. Remove temporary helper-routing dependence only after objective-level parity gates pass.

## Scope
In scope:
1. `src/jksum.c` CVLS LL/LP objective assembly and leave-one-out algebra.
2. LL/LP CVLS routing gates and test-only toggles used to compare helper vs non-helper.
3. Tests/artifacts proving objective parity at fixed bandwidths and optimizer parity at selected seeds.

Out of scope:
1. New statistical methods.
2. API changes.
3. CVAIC changes except regression-risk checks.

## Non-negotiable correctness contract
1. For fixed data + fixed bandwidth vector, helper and non-helper CVLS objective (`fval`) must match within tight tolerance.
2. No route asymmetry: session/attach/profile must stay rank-symmetric for migrated path.
3. No regression to deadlock-prone dispatch behavior.

## Current evidence summary
1. Seed-paired optimizer runs show LL/LP objective divergence between helper on/off (~1-4% relative in sampled runs).
2. LC shows parity in sampled runs.
3. Divergence appears structured (not random floating-point jitter), so root cause likely objective assembly mismatch.
4. Fixed-candidate comparator now confirms divergence at identical candidate bandwidths (non-tree default): objective mismatch exists prior to optimizer-path differences.

## Execution phases

## Phase 0: Oracle lock and reproducible baseline
Status: completed

Tasks:
1. Lock A/B control using process-level `NP_CVLS_STABLE_HELPER=0/1` (avoid rank-local option desync).
2. Produce baseline artifact comparing LL/LP/LC on same seeds.
3. Save artifacts under `/tmp/spmd_cvls_oracle_parity_YYYYMMDD_*`.

Acceptance:
1. Artifacts contain paired on/off optimizer metrics and objective deltas.
2. Completed artifact roots:
   - `/tmp/spmd_ab_switch_20260304_05`
   - `/tmp/spmd_cvls_oracle_parity_20260304_01`

## Phase 1: Fixed-bandwidth objective isolation (remove optimizer confound)
Status: completed

Tasks:
1. Build harness to evaluate CVLS at fixed `bws` for helper and non-helper.
2. Compare `fval` deltas per `(regtype, seed, bw_eval)`.
3. Partition by mode flags (tree on/off, gather-scatter branch relevant cases).

Acceptance:
1. If fixed-bandwidth deltas persist, mismatch is in objective evaluation math.
2. If fixed-bandwidth deltas vanish, mismatch is optimizer-path sensitivity only.
3. Outcome: fixed-candidate deltas persist for LL/LP, therefore mismatch is in objective assembly math (not merely optimizer drift).

## Phase 2: Branch-localize mismatch in non-helper LL/LP core
Status: completed

Tasks:
1. Isolate tree branch vs non-tree branch.
2. Isolate drop-one vs gather-scatter subpath behavior.
3. Instrument intermediate moment terms for a tiny deterministic case (single-seed, small `n`) in debug harness.
4. Implemented one-shot in-core comparator instrumentation in `src/jksum.c`:
   - env `NP_CVLS_COMPARE_HELPER=1` to compare core vs helper on the same candidate call,
   - env `NP_CVLS_COMPARE_ONCE=1` for single diagnostic emission,
   - reports total CV delta and max per-`j` term delta (`j` index + both terms).

Acceptance:
1. One minimal branch/subpath identified as first divergence point.
2. Outcome:
   - First divergence point is the legacy non-tree gather/scatter CVLS subpath, not optimizer drift.
   - Branch-local trigger:
     - LP: `src/jksum.c` gather path around `if (...) else { ... gather-scatter ... }` and post-call `RW` reconstruction loop (near lines ~8333 and ~8486 in current file).
     - LL: analogous block (near lines ~9244 and ~9419 in current file).
   - With `NP_CVLS_ROUTE_DROPONE=1` (default), fixed-candidate helper delta is machine precision.
   - With `NP_CVLS_ROUTE_DROPONE=0`, fixed-candidate helper delta reappears (large, structured).
   - Disabling `RW` correction (`NP_CVLS_DISABLE_RW_FIX=1`) makes mismatch far worse, confirming divergence is tied to gather/reconstruction algebra and not a simple sign flip.
   - CVAIC smoke (LL/LP) is identical under `NP_CVLS_ROUTE_DROPONE=0/1`, confirming this issue is CVLS-specific.

Artifacts:
1. `/tmp/spmd_cvls_oracle_parity_20260304_01/recheck_default_afterexp.log`
2. `/tmp/spmd_cvls_oracle_parity_20260304_01/recheck_gather_afterexp.log`
3. `/tmp/spmd_cvls_oracle_parity_20260304_01/recheck_gather_rw_on.log`
4. `/tmp/spmd_cvls_oracle_parity_20260304_01/recheck_gather_rw_off.log`
5. `/tmp/spmd_cvls_oracle_parity_20260304_01/cvaic_route1.log`
6. `/tmp/spmd_cvls_oracle_parity_20260304_01/cvaic_route0.log`

## Phase 3: Patch non-helper objective assembly to oracle equivalence
Status: pending

Tasks:
1. Apply surgical fix to identified LL/LP non-helper subpath.
2. Preserve existing hoists/caches/tree gates unless directly implicated.
3. Keep MPI collective cadence unchanged unless strictly required.

Acceptance:
1. Fixed-bandwidth parity passes for LL/LP across test grid.
2. No new hangs/timeouts in session/attach/profile smoke.

## Phase 4: Validation gates and cleanup
Status: pending

Tasks:
1. Re-run optimizer-level parity on seed panel (LL/LP/LC).
2. Run existing smoke/tests for touched paths (`npregbw`, `npreg`, route matrix subset).
3. Update remediation docs and mark checkpoints.

Acceptance:
1. Objective-level parity satisfied.
2. Route smokes pass.
3. Existing regression tests pass for touched scope.

## Phase 5: Optional helper deactivation decision
Status: pending

Tasks:
1. Decide whether helper remains fallback-oracle only or is fully disabled for LL/LP CVLS.
2. If disabling, keep a test-only override for forensic comparisons.

Acceptance:
1. Decision documented with evidence and reproducible artifacts.

## Required artifacts per checkpoint
1. Repro script path and exact invocation.
2. Seed list and DGP definitions.
3. Table of `fval` deltas (abs/rel) at fixed `bws`.
4. Route smoke outcomes (`session`, `attach`, `profile` where feasible).
5. Artifact root path under `/tmp`.

## Rollback triggers
1. Any new hang/deadlock.
2. Objective drift worsens after patch.
3. Route-mode failure introduced in touched functions.
