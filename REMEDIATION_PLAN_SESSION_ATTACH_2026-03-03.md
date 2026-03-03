# npRmpi Remediation Plan (2026-03-03)

## 1) Candid Assessment

You are right: the session/attach speed-equality objective was not met in the latest tranche.  
The work increased complexity in autodispatch/hot paths and produced measurable regressions in most tested families.

What this plan does:

1. Revert to baseline hot paths anchored at `318e89a9400e568b04d95147b2842f581b0e30be`.
2. Retain only hardening that is outside steady-state compute hot paths (or effectively zero-cost when disabled).
3. Run a strict A/B validation protocol before any commit.
4. Resume optimization as one-risk-at-a-time micro-tranches with stop-loss gates.

## 2) Non-Negotiable Pass Criteria

### A) Runtime parity to baseline commit

Target baseline: `318e89a9400e568b04d95147b2842f581b0e30be`.

Pass gates:

1. Per family (`npreg`, `npudens`, `npudist`, `npcdens`, `npcdist`) and per mode (`session_auto`, `attach_auto`):
   - median elapsed delta vs baseline `<= +2.0%`
   - mean elapsed delta vs baseline `<= +3.0%`
2. Aggregate across the 5 families:
   - geometric mean delta for `session_auto` `<= 0%`
   - geometric mean delta for `attach_auto` `<= 0%`
3. `profile_manual` must not regress beyond noise:
   - median delta `<= +1.0%`

### B) No regressions

1. Numerical parity at fixed seed and varying seed.
2. No new errors/warnings/condition drift on touched APIs.
3. No hangs (including one-rank-spin signature).
4. 4-path route gate: serial + session + attach + profile.
5. Full demo smoke (`n=100`, all `*_serial.R`, `*_npRmpi_attach.R`, `*_npRmpi_profile.R`).

### C) Full report artifact

Must include raw logs + summaries under dated `/tmp` path and a pre/post table.

### D) Commit checkpoint

Only after all gates pass:

`revert hot-path experiments; retain fail-fast MPI hardening; verify parity/perf vs 318e89`

## 3) What Was Attempted and Failed (Current Tranche)

Observed net outcome from tested snapshot:

1. `npreg`: session `+9.55%`, attach `+15.88%`
2. `npudens`: session `+38.21%`, attach `+10.75%`
3. `npudist`: session `+11.15%`, attach `+2.12%`
4. `npcdens`: session `+33.09%`, attach `+23.08%`
5. `npcdist`: session `+34.98%`, attach `-14.45%` (single attach improvement)

Attempt groups introduced:

1. Autodispatch argument caching / prepublish / remote-reference reuse machinery in `R/np.autodispatch.R`.
2. Additional pre-dispatch validation and gating inside estimator entry paths (`R/np.regression.R`, `R/np.density*.R`, `R/np.distribution*.R`, `R/np.condensity.R`, `R/np.condistribution.R`).
3. Session/attach/shutdown worker hardening in `R/session.R`, `inst/slavedaemon.R`, `R/Rparutilities.R`.
4. Expanded docs/tests for new behavior.

Critical interpretation:

1. Hardening work is useful.
2. Performance work was over-broad (too many simultaneous risk axes).
3. We lacked tranche isolation, so root-cause attribution became slow and noisy.

## 4) Revert-and-Retain Scope

### Revert (hot paths)

Revert to baseline for:

1. `R/np.autodispatch.R`
2. `R/np.regression.R`
3. `R/np.density.R`
4. `R/np.density.bw.R`
5. `R/np.distribution.R`
6. `R/np.distribution.bw.R`
7. `R/np.condensity.R`
8. `R/np.condistribution.R`

### Retain (hardening)

Retain (or re-apply if needed) only:

1. Profile startup hard-fail guard in `inst/Rprofile` (dual `R_PROFILE_USER`/`R_PROFILE` conflict).
2. Worker receive-timeout fail-fast controls and diagnostics:
   - `R/session.R`
   - `inst/slavedaemon.R`
3. Safer shutdown/teardown warning path:
   - `R/Rparutilities.R`
4. Documentation + tests strictly tied to above hardening:
   - `man/npRmpi.session.Rd`
   - `man/np.options.Rd`
   - focused test files only.

## 5) Validation Execution Plan (Fast + Comprehensive)

Artifacts root:

`/tmp/nprmpi_remediation_20260303_<timestamp>/`

### Phase 0: Pre-flight safety

1. Save patch snapshot of current dirty state.
2. Confirm no masking wrappers in `demo/*.R`.
3. Ensure no stale MPI workers before runs.

### Phase 1: Build A/B baselines

1. `baseline`: clean worktree at `318e89...`
2. `candidate`: reverted hot paths + retained hardening
3. Install each to isolated `R_LIBS`.

### Phase 2: Route and functional gates

1. Route validators: `manual-broadcast`, `attach`, `session n=1`, `session n=0`.
2. Full demo sweep in `demo` with `n=100`, serial/attach/profile triplets.
3. `testthat` full run (with subprocess route tests enabled).

### Phase 3: Numerical parity gates

1. Fixed seed and varying seed.
2. Compare objective values, bandwidths, predictions.
3. Fail on unexplained drift.

### Phase 4: Performance gates

1. Run identical matrix harness used in current timing snapshot.
2. Collect per family/mode elapsed for baseline and candidate.
3. Compute mean/median deltas; enforce Section 2 gates.
4. Repetition policy (required for accept/reject):
   - `times=5` is smoke only (no keep/drop decisions).
   - decision screening uses interleaved paired-seed `times=25`.
   - predeclare MEI (absolute seconds + relative percent) for each gated method/mode.
   - compute paired deltas and 95% CI for mean delta.
   - allow early accept at `times=25` only for obvious wins (mean+median favorable, CI excludes zero, effect exceeds MEI across all gated methods/modes).
   - otherwise escalate to interleaved paired-seed `times>=100`; if pilot-implied required sample is above `100`, run higher `times` or mark inconclusive for that MEI.
   - performance-claim acceptance requires paired mean and paired median agreement with claimed direction at decision tier.
5. Cross-method agreement policy:
   - all gated methods/modes must align in direction for a performance-claim commit;
   - anomalies must be explicitly rationalized and bounded relative to MEI, otherwise reject/revert.
6. MPI performance mode default:
   - use `nslaves=1` as primary performance gate;
   - run `nslaves=0` only for explicit master-only overhead diagnostics.
7. Reasoning:
   - large effects are detectable at low `times`, while small effects require larger `times` since uncertainty contracts at about `1/sqrt(n)`.

### Phase 5: Report + checkpoint

1. Produce markdown summary with command lines, tables, and verdict.
2. Commit only if all gates pass.

## 6) Micro-Optimization Plan (After Stabilization)

Rule: one micro-change per tranche, each independently benchmarked and revertable.

Candidate sequence:

1. Dispatch-call overhead profiling only (no behavior changes).
2. Optimize only top single hotspot in autodispatch call preparation.
3. Re-run full gates; keep only if neutral/better.
4. Repeat.

Hard stop rules:

1. Two consecutive failed tranches => pause and redesign.
2. Any hang regression => immediate rollback of that tranche.
3. No commit without route + parity + perf gates all green.

## 7) Success Likelihood (Critical)

### For revert-and-retain stabilization

High (about 85-95%) if scope is kept strict and hot-path files are fully reset.

### For session/attach reaching profile-level runtime broadly

Moderate-to-low (about 35-55%) without route-specific API changes, because profile/manual-broadcast is structurally lower overhead.

### For meaningful session/attach improvements without profile regressions

Moderate (about 60-75%) via small, measured micro-optimizations and strict stop-loss gates.

## 8) Additional Controls to Avoid Another Time Sink

1. Pre-register each tranche hypothesis before coding:
   - expected hotspot
   - expected gain
   - explicit rollback condition
2. Maximum changed files per tranche: 1-2 runtime files.
3. Mandatory before/after benchmark run on same host, same seed policy, same harness.
4. No “optimistic” claims until numbers clear gates.
5. Daily checkpoint summary: kept/reverted/gains/losses with raw artifact links.

## 9) Operational Guarantee Language

Absolute mathematical guarantee is impossible.  
Operational guarantee is enforced by hard gates: if any gate fails, the tranche is reverted and not committed.

## 10) External Diagnosis Audit (Critical)

An external diagnosis was provided. Assessment below is evidence-graded.

### A) Autodispatch argument-cache overhead is a primary culprit

Verdict: **high confidence (valid)**.

Evidence in current diff:

1. `R/np.autodispatch.R` has large net growth (`+316/-38`).
2. Per-argument path now includes cache-key and cache-env lookups plus full-object work:
   - `format(caller_env)`
   - `object.size(val)`
   - `identical(val, ent$value)`
   - `exists/get/assign` churn in cache env
3. This executes in the dispatch path repeatedly, so overhead compounds in CV/iterative calls.

### B) Pre-dispatch validation duplication contributes materially

Verdict: **high confidence (valid)**.

Evidence:

1. Added `toFrame`, `%~%`, and `npKernelBoundsCheckEval`/`npKernelBoundsResolve` before autodispatch returns in multiple estimator entry points.
2. Equivalent checks and data handling already exist downstream in the standard estimator/bw execution path.
3. Net effect is duplicated validation work in the same call lifecycle.

### C) `bcast_robj_by_name` rewrite likely increases overhead

Verdict: **moderate-to-high confidence (likely valid)**.

Evidence:

1. Baseline used `mpi.bcast.Robj2slave` path.
2. Current path wraps values in `substitute(assign(...))` and uses command broadcast execution.
3. For large objects, this changes serialization/execution behavior and is plausibly slower.

Note:

Direct microbenchmark attribution per payload size is still required for final proof, but this is a credible regression source and should be reverted in stabilization.

### D) Formula dispatch hook additions are a universal culprit

Verdict: **conditional (partially valid)**.

Evidence:

1. New formula-method autodispatch guards were added in multiple files.
2. Impact depends on workload mix: only calls entering formula methods pay this added layer.
3. This is likely a secondary contributor, not the sole primary driver.

## 11) Plan Updates from Audit

The following are now explicit plan constraints.

1. Revert all argument-cache machinery from hot path in stabilization phase.
2. Revert pre-dispatch duplicate validation from estimator entry points in stabilization phase.
3. Revert `bcast_robj_by_name` back to baseline broadcast path in stabilization phase.
4. Keep hardening-only worker/profile guards (`inst/Rprofile`, `R/session.R`, `inst/slavedaemon.R`, `R/Rparutilities.R`) as retained scope.

## 12) Micro-Optimization Guardrails (Revised)

If cache optimization is retried later, it must satisfy all of the following.

1. Opt-in only, default off.
2. Single-file tranche (`R/np.autodispatch.R`) with no estimator-entry modifications.
3. No full-object equality checks (`identical` on large frames) in hot path.
4. No full-tree size scans (`object.size`) in hot path.
5. No new package dependency solely for identity probing.
6. No reliance on non-API internals such as `.Internal(inspect())`.
7. Profile route must remain neutral (no measurable regression against baseline).

## 13) Rejected/Unsafe Suggestions

The following suggestions are explicitly rejected for this codebase unless separately justified.

1. `data.table::address()` as a required mechanism:
   - introduces a new dependency and still does not guarantee semantic immutability.
2. `.Internal(inspect())` for identity:
   - unsupported internal API, brittle across R versions.
3. Broad multi-file optimization passes:
   - violates one-risk-axis tranche rule and obscures attribution.
