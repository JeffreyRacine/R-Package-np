# Remediation Patch Design: Parallelize Stable CVLS Helper Under SPMD (2026-03-04)

## Status update (2026-03-04, superseding helper-parallelization intent)
1. `np_reg_cv_ls_stable_ll_glp` is currently absent in both repos (`np-master`, `np-npRmpi`).
2. No active call sites route through a stable-helper CVLS path in current `src/jksum.c`.
3. Therefore, a helper-parallelization code tranche is not currently actionable.
4. This document now serves as:
   - historical rationale for why helper-parallelization had been considered,
   - active guidance for adjacent low-risk remediation (dead-scaffold cull + modernization checks + route gates).

## Scope
Target only the CVLS stable-helper path used by LL/LC/LP in `npRmpi`:
- `src/jksum.c`: `np_reg_cv_ls_stable_ll_glp(...)`
- caller sites in `np_kernel_estimate_regression_categorical_ls_aic(...)`

Out of scope for this patch:
- `zzz.R` check-mode lifecycle changes (`npRmpi.reuse.slaves` / unload hard teardown)
- non-regression families (`npudens*`, `npcdens*`, etc.)
- algorithmic/statistical changes

## Why this patch now
Session mode is now SPMD-orchestrated, so rank-symmetric collective entry is enforced at the dispatch layer. That removes the original blocker for parallelizing this helper's outer loop.

The helper remains active and currently runs as per-rank serial work (with inner kernel call using `suppress_parallel=1`).

## Current state (verified)
1. Stable helper symbol/call path is absent in current heads of both repos.
2. Dead scaffolding cull completed for unused `.npRmpi_spmd_next_seq` in `np-npRmpi/R/np.autodispatch.R`.
3. Runtime `nslaves=0` remains rejected in `npRmpi.init` (active contract); only legacy positive-integer guard remains in `R/Rcomm.R::mpi.comm.spawn`.
4. Large-memory kernel-weight path:
   - `return.kernel.weights=TRUE` allocates `n.train * n.eval` storage.
   - In dominant plot/bootstrap usage, `n.eval` is typically fixed/small (often ~100), so growth is effectively linear in `n.train` for fixed evaluation grids.
   - This path is not the same class of risk as unrestricted `O(n^2)` over `n.train * n.train`.

## Design goals
1. Preserve estimator semantics and API.
2. Keep strict collective symmetry under MPI.
3. Avoid nested/asymmetric collective patterns.
4. Improve wall-time for large `n` under `iNum_Processors > 1`.
5. Keep serial behavior identical when MPI is absent.

## Core design

### A) Keep inner kernel call serial-safe
Do **not** change helper's inner `kernel_weighted_sum_np_ctx` call contract:
- retain `suppress_parallel=1`

Rationale: the helper owns outer-loop parallelization; inner parallelization would reintroduce nested collective risk and cadence complexity.

### B) Parallelize only outer `j` loop
Use strided ownership in MPI builds:
- rank `r` computes `j = r, r + iNum_Processors, ...`
- each rank accumulates `cv_local`
- reduce once at end:
  - `MPI_Allreduce(&cv_local, &cv, 1, MPI_DOUBLE, MPI_SUM, comm[1])`

Serial build behavior:
- existing `for (j=0; j<num_obs; j++)` path unchanged.

### C) Failure contract
Introduce rank-local failure flag for singular/solver abort paths inside helper loop:
1. `local_fail` set on irrecoverable per-`j` failure.
2. `MPI_Allreduce(..., MPI_MAX, ...)` to compute `any_fail`.
3. if `any_fail`: return `DBL_MAX` (uniform across ranks), skipping value reduce.

This preserves identical collective order across ranks (fail-check reduce always called in MPI branch).

### D) Maintain parity with `np-master`
- Keep math/formula logic aligned with `np-master` helper body.
- `npRmpi`-specific delta should be limited to MPI loop partition + reductions.
- If desired for textual parity, add the same structural refactor in `np-master` guarded by non-MPI path only.

## Candidate implementation shape

Option preferred (minimal duplication):
1. Extract per-`j` contribution logic into a static helper:
   - computes one `j` SSE contribution
   - returns status (ok/fail) + contribution
2. Existing stable helper becomes orchestration wrapper:
   - serial loop in non-MPI
   - strided loop + allreduce in MPI

Alternative (faster edit, less clean):
- Keep single function and add `#ifdef MPI2` loop partition in place.

## Call-site behavior
No call-site API changes required if orchestration remains internal to `np_reg_cv_ls_stable_ll_glp(...)`.

If introducing a second symbol (e.g., `_spmd` variant), route both current call sites to it:
- `src/jksum.c:7922`
- `src/jksum.c:8616`

## Collective safety checklist (non-negotiable)
1. Exactly one helper-level `MPI_Allreduce` for fail flag (MPI builds).
2. Conditional second `MPI_Allreduce` for `cv` only when `any_fail == 0` (branch is rank-symmetric).
3. No new collectives inside per-`j` worker body.
4. No rank-conditional early returns before helper-level reduction points.

## Validation gates

### Numerical parity
1. Fixed-seed parity for LL/LC/LP CVLS:
   - `np-master` serial vs `npRmpi` session `nslaves=1`
   - compare `fval`, `num.feval`, selected bandwidths
2. Varying-seed smoke for stability (no divergence, no hang).

### Route matrix
1. `np` serial (`np-master`)
2. `npRmpi` session (`npRmpi.init(nslaves=1)`)
3. `npRmpi` attach (`mpiexec -n 2`)
4. `npRmpi` profile/manual-broadcast (`mpiexec -n 2` + profile)

### Performance
1. Screening `times=25` paired interleaved seeds on LL/LC/LP CVLS workloads.
2. Escalate `times>=100` if signal ambiguous.
3. Require non-regression in unaffected families.

### Safety
1. Historical session hang repro passes.
2. No orphan `slavedaemon.R`/`Rslaves.sh` after failures/timeouts.
3. `R CMD check --as-cran` unchanged or improved for this path.

## Risks and mitigations
1. Floating-point reduction-order drift:
   - Mitigation: tolerance-based parity for `fval`, document expected ulp-scale differences.
2. Hidden inner collectives from future edits:
   - Mitigation: comment + test asserting helper calls kernel with `suppress_parallel=1`.
3. Regression masked by session-only tests:
   - Mitigation: enforce attach/profile route gates.
4. Memory pressure in kernel-weight path:
   - Mitigation: keep current path (needed for plotting/bootstrap workflows), but add explicit telemetry/doc notes and avoid broad shape changes in this tranche.

## Updated critique (ROI / risk)
### High ROI / low risk (do now)
1. Remove dead/near-dead legacy scaffolding tied to retired semantics:
   - remove unused `.npRmpi_spmd_next_seq` if confirmed unreachable,
   - remove stale `nslaves==0` scaffolding where still present.
2. Keep stable-helper patch narrow (outer-loop MPI orchestration only), with no API or estimator-math changes.
3. Add static modernization checks on touched files (Wickham-aligned hygiene: scalar safety, eval surface, namespace hygiene, lifecycle cleanup checks).

### Medium ROI / low risk (do in same stream with gates)
1. Add explicit comments/tests around `return.kernel.weights=TRUE` memory shape (`n.train * n.eval`) so behavior is intentional and discoverable.
2. Add targeted memory-smoke scripts for representative plot/bootstrap sizes to catch accidental shape explosions.

### Medium ROI / medium risk (separate tranche, gated; currently N/A)
1. Helper outer-loop parallelization is deferred unless a stable-helper path is intentionally reintroduced.

### Lower ROI / higher risk (defer)
1. Broad helper rewrites or redesign of plot/bootstrap allocation behavior in this same patch.
2. Any refactor that mixes algorithmic changes with orchestration changes.

## Clarification-Driven Decision Matrix (Pro / Con)
### Item A: `return.kernel.weights=TRUE` (`n.train * n.eval`) gateway
1. Pro:
   - required by plot/bootstrap workflows,
   - for common plotting defaults (`n.eval` bounded/small), effective growth is linear in `n.train`,
   - not the same failure class as unrestricted `n.train * n.train` estimator/CV storage.
2. Con:
   - still quadratic in the general two-axis shape and can spike memory if `n.eval` is user-expanded,
   - can be misread as a core-path anti-pattern without explicit documentation.
3. Decision:
   - retain as intentional,
   - add explicit documentation + telemetry checks,
   - do not refactor shape in this tranche.

### Item B: Remove dead/near-dead legacy scaffolding (`.npRmpi_spmd_next_seq`, retired `nslaves==0`)
1. Pro:
   - immediate maintainability gain with near-zero algorithm risk,
   - reduces route ambiguity and stale-contract hazards.
2. Con:
   - tiny risk of deleting a latent callsite if reachability is not verified first.
3. Decision:
   - do in a dedicated first tranche with static reachability proof + targeted lifecycle tests.

### Item C: Wickham-guided modernization sweep integrated into this remediation stream
1. Pro:
   - reduces future technical debt while code context is fresh,
   - catches fragile patterns early (`eval` surface, scalar controls, condition-contract drift).
2. Con:
   - scope creep risk if mixed with behavior-changing patches.
3. Decision:
   - run as gated micro-tranches only, separated from algorithm/orchestration deltas.

## Modernization sweep integration (Wickham-guided)
Treat this remediation as a dual track: correctness + modernization hygiene.

References:
1. `/Users/jracine/Development/WICKHAM_FINAL_MODERNIZATION_SWEEP_PLAN_20260302.md`
2. `/Users/jracine/Development/archive/R_AUDIT_WICKHAM_2026-02-22_20260226.md`
3. `/Users/jracine/Development/VALIDATION_GATES_CANONICAL_2026-02-28.md`
4. `/Users/jracine/Development/archive/MODERNIZATION_REMEDIATION_2026-03-01.md`
5. `/Users/jracine/Development/archive/R_LAYER_MODERNIZATION_BACKLOG_20260226.md`
6. `/Users/jracine/Development/archive/SESSION_CONSOLIDATION_MODERNIZATION_2026-02-24_20260226.md`

Embedded modernization checks for each tranche:
1. R-layer hygiene scan on touched files (`eval(`, `do.call(`, `sapply(`, scalar-control guards, `on.exit(..., add=TRUE)` in lifecycle mutations).
2. Condition-contract stability (errors/warnings/messages unchanged unless explicitly documented).
3. Tarball-first validation and MPI route matrix evidence before checkpoint acceptance.
4. No performance-claim acceptance without paired-seed gate evidence.

## Interaction with current check-mode work
This patch is independent of the check-mode lifecycle patch in `zzz.R`.
Recommended order:
1. Land/validate check-mode lifecycle fix first (to stabilize check harness).
2. Land stable-helper MPI patch as a separate tranche.

## Acceptance criteria
1. No deadlocks across session/attach/profile.
2. LL/LC/LP CVLS stable-helper route numerically matches baseline within declared tolerance.
3. Large-`n` wall-time improves or is neutral under MPI.
4. Parity contract with `np-master` is preserved at estimator-output level.
5. Modernization checks for touched files pass with no condition-contract regressions.

## Execution intent (no-breakage sequence)
1. First tranche: dead-code cull only (`.npRmpi_spmd_next_seq`, residual `nslaves==0` scaffolding), no behavior changes.
2. Validate with full gate pack (serial + session + attach + profile/manual, plus existing smoke/tests).
3. Second tranche: modernization micro-sweep on touched files only (no estimator/CV algorithm changes).
4. Re-run parity + route gates; rollback immediately on any regression/hang.
5. Third tranche: documentation closure updates and artifact capture.

## Intended next actions (explicit)
1. Preserve helper-absence invariant unless explicitly requested to reintroduce a helper path.
2. Keep route-smoke + contract tests mandatory after any runtime/lifecycle cull.
3. Continue Wickham micro-cleanups in isolated tranches with no behavior/default drift.
4. Use full gate-pack evidence before marking the remediation stream complete.

## Execution log (2026-03-04)
1. Completed dead-scaffold cull tranche:
   - removed unreachable `.npRmpi_spmd_next_seq` from `R/np.autodispatch.R`,
   - retained active runtime `nslaves<1` rejection guard in `npRmpi.init`.
2. Validation artifact root:
   - `/tmp/spmd_tranche1_deadscaffold_20260304_103024`
   - tokens: `SESSION_OK=1`, `ATTACH_OK=1`, `PROFILE_OK=1`, `MANUAL_OK=1`, `NPSIGTEST_FAILNZ=0`, orphan pre/post empty.
3. Checkpoint commits:
   - `635ddbb` (dead-scaffold cull + remediation-doc updates),
   - `c950f14` (status correction: helper-parallelization tranche marked N/A on current heads).
4. Completed modernization micro-sweep tranche:
   - hardened `R/Rcomm.R::mpi.comm.spawn` `nslaves` validation to require scalar positive integer while preserving existing error text (`\"Choose a positive number of slaves.\"`),
   - added coverage in `tests/testthat/test-rcomm-arg-contract.R` for invalid `nslaves` forms (`0`, `NA`, vector length > 1, non-integer).
5. Tranche-2 validation artifact root:
   - `/tmp/spmd_tranche2_modsweep_20260304_103653`
   - tokens: `SESSION_OK=1`, `ATTACH_OK=1`, `PROFILE_OK=1`, `MANUAL_OK=1`, `RCOMM_FAILNZ=0`, `NPSIGTEST_FAILNZ=0`, `ORPHAN_PRE_EMPTY=1`, `ORPHAN_POST_EMPTY=1`.
6. Completed modernization-tooling tranche:
   - added scan harness `issue_notes/run_modernization_micro_scan.sh` (shellcheck-clean),
   - generates Wickham-aligned static inventory for R-layer hygiene, stale-branch checks, and helper-symbol presence.
7. Tranche-3 scan artifact root:
   - `/tmp/spmd_tranche3_modscan_20260304_104115`
   - summary highlights: `scan_demo_masking=0`, `scan_dotC=0`, `scan_manual_distributed_call=0`, `scan_stable_helper_symbol=0`, `scan_nslaves_zero_runtime=1` (expected guard in `R/Rcomm.R`), `scan_nslaves_zero_tests=2` (expected rejection coverage).
8. Completed scan-harness refinement tranche:
   - tightened regexes in `issue_notes/run_modernization_micro_scan.sh` to reduce false positives (`deparse`, `tikh.eval`, commented `library/require`).
9. Tranche-4 scan artifact root:
   - `/tmp/spmd_tranche4_modscan_refine_20260304_104220`
   - summary highlights: `scan_eval=4`, `scan_parse=3`, `scan_runtime_library_require=0`, `scan_demo_masking=0`, `scan_dotC=0`, `scan_manual_distributed_call=0`, `scan_stable_helper_symbol=0`.
10. Completed low-risk comment cleanup tranche:
   - removed stale commented `parse(...)` fragments in `R/Rcoll.R` that were no longer part of active execution paths.
11. Tranche-5 validation + scan artifacts:
   - route gate root: `/tmp/spmd_tranche5_gate_20260304_104321`
   - gate tokens: `SESSION_OK=1`, `ATTACH_OK=1`, `PROFILE_OK=1`, `MANUAL_OK=1`, `NPSIGTEST_FAILNZ=0`, `ORPHAN_PRE_EMPTY=1`, `ORPHAN_POST_EMPTY=1`.
   - refreshed scan root: `/tmp/spmd_tranche5_modscan_cleanup_20260304_104312` (`scan_parse=0`).
12. Completed IV-family type-stability tranche:
   - replaced classifier `sapply(..., is.numeric)` calls with fixed-type `vapply(..., logical(1))` in:
     - `R/npregiv.R` (`X.col.numeric`, `xdat.numeric`, `z.numeric`, `w.numeric`)
     - `R/npregivderiv.R` (`z.numeric`, `w.numeric`)
   - no API/default/condition-message changes; classification semantics preserved.
13. Tranche-6 validation artifact root:
   - `/tmp/spmd_tranche6_iv_vapply_20260304_104740`
   - tokens: `SESSION_IV_OK=1`, `ATTACH_OK=1`, `PROFILE_OK=1`, `MANUAL_OK=1`, `NPREGIV_FAILNZ=0`, `NPREGIV_STATE_FAILNZ=0`, `NPSIGTEST_FAILNZ=0`, `ORPHAN_PRE_EMPTY=1`, `ORPHAN_POST_EMPTY=1`.
14. Tranche-6 modernization-scan artifact root:
   - `/tmp/spmd_tranche6_modscan_postvapply_20260304_104841`
   - summary highlights: `scan_sapply=67` (inventory reduced), `scan_eval=4`, `scan_parse=0`, `scan_runtime_library_require=0`, `scan_demo_masking=0`, `scan_dotC=0`, `scan_manual_distributed_call=0`, `scan_stable_helper_symbol=0`.
15. Completed `W.lp` type-stability tranche:
   - replaced `xdat.col.numeric` classifier in `R/util.R` from `sapply(..., is.numeric)` to typed `vapply(..., logical(1))`.
   - no API/default/condition-message changes.
16. Tranche-7 validation artifact root:
   - `/tmp/spmd_tranche7_wlp_vapply_20260304_105424`
   - gate tokens (`summary_tokens_gate.log`): `ORPHAN_PRE_EMPTY=1`, `NPREG_GLP_FAILNZ=0`, `NPREG_CORE_FAILNZ=0`, `NPSIGTEST_FAILNZ=0`, `SESSION_WLP_OK=1`, `ATTACH_OK=1`, `PROFILE_OK=1`, `MANUAL_OK=1`, `ORPHAN_POST_EMPTY=1`.
   - note: `test-npreg-arg-contract.R` is not a valid installed-package gate file because it references unexported symbol `npreg.rbandwidth` by bare name; recorded as `NPREG_ARG_INSTALLED_MODE_NA=1`.
17. Tranche-7 modernization scan:
   - `run_modernization_micro_scan.sh` summary recorded in `/tmp/spmd_tranche7_wlp_vapply_20260304_105424/modscan_run.log`.
   - key delta: `scan_sapply=66` (down from `67`).
18. Completed modernization-tooling I/O hardening tranche:
   - `issue_notes/run_modernization_micro_scan.sh` now accepts explicit positional args `ROOT_DIR` and `OUT_DIR` (`$1`, `$2`) while preserving defaults.
   - added fail-fast root guard (`ROOT_DIR/R` must exist).
19. Tranche-8 tooling validation artifact root:
   - `/tmp/spmd_tranche8_scanio_20260304_105725`
   - `shellcheck` clean; scan run confirms `out_dir=/tmp/spmd_tranche8_scanio_20260304_105725/scan` and `TRANCHE8_OK=1`.
20. Completed plot-engine factor-classifier tranche:
   - replaced `all.isFactor` builders from `sapply(..., is.factor)` to typed `vapply(..., logical(1))` in:
     - `R/np.plot.engine.plbandwidth.R`
     - `R/np.plot.engine.scbandwidth.R`
     - `R/np.plot.engine.condbandwidth.R`
     - `R/np.plot.engine.conbandwidth.R`
   - includes removal of legacy `unlist(sapply(zdat, is.factor))` shape coercion in favor of direct typed vector output.
21. Tranche-9 validation artifact root:
   - `/tmp/spmd_tranche9_plot_vapply_20260304_110304`
   - tokens: `ORPHAN_PRE_EMPTY=1`, `PLOT_AUTODISPATCH_FAILNZ=0`, `PLOT_GUARDRAILS_FAILNZ=0`, `PLOT_COEF_FAILNZ=0`, `PLOT_MPI_BOOT_FAILNZ=0`, `SESSION_PLOT_OK=1`, `ATTACH_OK=1`, `PROFILE_OK=1`, `MANUAL_OK=1`, `ORPHAN_POST_EMPTY=1`.
22. Tranche-9 modernization scan:
   - `/tmp/spmd_tranche9_plot_vapply_20260304_110304/modscan/summary.log`
   - key delta: `scan_sapply=62` (down from `66`).
