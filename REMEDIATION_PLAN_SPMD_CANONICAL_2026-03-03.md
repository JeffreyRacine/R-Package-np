# npRmpi Canonical Remediation Plan: Session SPMD Unification (2026-03-03)

## Canonical Status
1. This is the active remediation plan for session/attach/profile unification.
2. `nslaves=0` is removed permanently from `npRmpi`.
3. `npRmpi` runtime contract is MPI-only (`nslaves>=1`).

## Execution Status (Updated 2026-03-03/04)
Completed checkpoint tranches:
1. Phase 1 (`nslaves=0` removal): complete (`2fc9317`).
2. Phase 2 (SPMD control-plane scaffolding): complete (`905ffd9`, `6b6be81`, `1fece7e`, `3bf4914`).
3. Phase 3 (`npregbw` LL/LP CV migration + `O(n^2)` helper excision): complete (`24dbd4b` and subsequent route/contract gates).
4. Phase 4 (core-family migration progress): complete for `npudens*`, `npudist*`, `npcdens*`, `npcdist*`, `npscoef*`, `npindex*`, `npplreg*`.
5. Phase 5 (hot-path legacy asymmetry decommission): complete for core estimator/CV hot paths via typed locked opcode handlers and guard-validated route parity.
6. Phase 6 (non-core migration): tranche A complete for `npqreg*`, `npconmode*`, `npksum*`, `npregiv*`, `npregivderiv*`, `npcmstest`, `npqcmstest`, `npdeneqtest`, `npdeptest`, `npsdeptest`, `npsymtest`, and `npunitest` via typed locked opcodes.

Latest checkpoint commits:
1. `8f5e959` (density/distribution opcode classification)
2. `47bbc49` (locked typed handlers for migrated core opcodes)
3. `388d643` (contract test realignment)
4. `9d163f6` (lock typed opcode for `npindexbw`)
5. `f4c83f0` (lock typed opcode for `npindex`)
6. `c34ff68` (lock typed opcodes for remaining core estimators)
7. `138401a` (fix forwarded-dot argument resolution in autodispatch)
8. `519955b` (lock non-core autodispatch families to typed SPMD opcodes)

Latest artifact roots:
1. `/tmp/spmd_canonical_20260304_0001/phase10_density_opcode_timeout_20260303_201934`
2. `/tmp/spmd_canonical_20260304_0001/phase11_locked_opcode_handlers_20260303_202742`
3. `/tmp/spmd_canonical_20260304_0001/phase12_individual_tests_20260303_204052`
4. `/tmp/spmd_canonical_20260304_0001/phase19_preflight_20260303_210733`
5. `/tmp/spmd_canonical_20260304_0001/phase20_noncore_opcode_locks_20260303_211050`

## Objective
Keep user-facing workflow unchanged (`npreg(...)`, `npregbw(...)`, etc.) while making internal execution rank-symmetric SPMD for MPI-sensitive paths in all modes:
1. `session/spawn`
2. `attach`
3. `profile/manual-broadcast`

Core family coverage for this plan:
1. `npreg*`
2. `npudens*`
3. `npcdens*`
4. `npudist*`
5. `npcdist*`
6. `npscoef*`
7. `npindex*`
8. `npplreg*`

Important scope clarification:
1. The core-family list is a migration priority list, not an allowlist.
2. All exported `npRmpi` functions (including non-core families such as `npsigtest`, `npsymtest`, `npcopula`) are expected to run.
3. Non-core families are currently supported under the same runtime contract and are migrated to typed locked opcodes in Phase 6.

## Non-Negotiable Invariants
1. No silent serial fallback when MPI route is selected.
2. No runtime `np::` bridge from `npRmpi`.
3. No masked failures in demos/tests.
4. No default drift unless explicitly approved.
5. One risk axis per tranche.
6. No `O(n^2)` helper-path reintroduction in estimation/CV internals.
7. Helper policy: prefer direct distributed kernels in MPI-sensitive CV/estimation hot paths.

## Baseline and Artifacts
1. Baseline commit for active branch: `45f0e2cbe81dcb4fb40e4d8b01ce8c54e23eeaf6`.
2. All tranche artifacts go under:
   - `/tmp/spmd_canonical_<timestamp>/phase0`
   - `/tmp/spmd_canonical_<timestamp>/phase1`
   - `/tmp/spmd_canonical_<timestamp>/phase2`
   - `/tmp/spmd_canonical_<timestamp>/phase3`
   - `/tmp/spmd_canonical_<timestamp>/phase4`
   - `/tmp/spmd_canonical_<timestamp>/phase5`
   - `/tmp/spmd_canonical_<timestamp>/phase6`

## Required Pre-Flight (Before Every Tranche)
1. `git status --short` for:
   - `/Users/jracine/Development`
   - `/Users/jracine/Development/np-master`
   - `/Users/jracine/Development/np-npRmpi`
2. Demo masking scan:
   - `rg -n --glob 'demo/*.R' '(^|[^#])\\btry\\s*\\(|(^|[^#])\\btryCatch\\s*\\('`
3. Orphan worker scan:
   - `ps aux | rg 'slavedaemon\\.R|Rslaves\\.sh'`
4. Route sanity (`npRmpi`, tiny smokes):
   - session/spawn (`nslaves=1`)
   - attach (`mpiexec -n 2`)

## Route/Smoke Commands (Existing Tools Only)
1. Manual-broadcast validator:
   - `Rscript --no-save /Users/jracine/Development/np-npRmpi/issue_notes/validate_route_manual_broadcast.R`
2. Attach validator:
   - `FI_TCP_IFACE=en0 mpiexec -n 2 Rscript --no-save /Users/jracine/Development/np-npRmpi/issue_notes/validate_route_attach.R`
3. Profile validator:
   - `R_PROFILE_USER=/Users/jracine/Development/np-npRmpi/inst/Rprofile FI_TCP_IFACE=en0 mpiexec -n 2 Rscript --no-save /Users/jracine/Development/np-npRmpi/issue_notes/validate_route_profile.R`
4. Session tiny smoke (`nslaves=1`) token check:
   - `Rscript --no-save -e 'suppressPackageStartupMessages(library(npRmpi)); npRmpi.init(nslaves=1, quiet=TRUE); on.exit(try(npRmpi.quit(), silent=TRUE), add=TRUE); set.seed(1); x<-runif(80); y<-rnorm(80); bw<-npregbw(y~x, regtype=\"lc\", bwmethod=\"cv.ls\", nmulti=1); fit<-npreg(bws=bw); stopifnot(inherits(fit,\"npregression\")); cat(\"SESSION_N1_OK\\n\")'`

## Target Internal Contract (SPMD Step Engine)
For migrated MPI-sensitive opcodes:
1. Master creates envelope: `{seq_id, opcode, args_ref, timeout_class}`.
2. Master broadcasts envelope once per step.
3. All ranks enter same opcode handler.
4. All ranks execute matching collective cadence for that opcode.
5. Completion uses explicit ACK/ERR semantics.
6. Sequence mismatch or timeout is hard-fail with route/opcode diagnostics.

## Phase Plan

## Phase 0: Baseline Freeze
Scope:
1. Snapshot current clean baseline and execution environment.
2. Pin reference logs for route, parity, and performance.

Validation:
1. Route validators all pass.
2. Targeted MPI-sensitive tests pass:
   - `tests/testthat/test-session-routing-subprocess-contract.R`
   - `tests/testthat/test-session-arg-contract.R`
   - `tests/testthat/test-jksum-gating-smoke.R`
   - `tests/testthat/test-ll-lp-degree1-parity.R`

Acceptance:
1. `/tmp` phase0 bundle includes commands, raw logs, and summary table.

## Phase 1: Remove `nslaves=0` Permanently
Scope:
1. Enforce `nslaves>=1` in runtime lifecycle.
2. Remove remaining active `nslaves==0` branches in runtime paths.
3. Update docs/tests/messages to reflect MPI-only contract.

Validation:
1. Negative contract test:
   - `npRmpi.init(nslaves=0)` hard-fails with explicit remediation (`use np for serial`).
2. Positive route check:
   - session smoke with `nslaves=1` passes.
3. Attach/profile/manual validators pass.
4. No orphan worker processes.

Acceptance:
1. No active runtime path supports `nslaves=0`.
2. User-facing message contract is intentional and tested.

## Phase 2: Control-Plane Scaffolding (No Method Logic Change)
Scope:
1. Add envelope structure and `seq_id` management.
2. Add typed opcode registry + dispatcher shell.
3. Add ACK/ERR completion scaffolding.
4. Add timeout/divergence diagnostics.

Validation:
1. One tiny test opcode works in session/attach/profile.
2. Deliberate sequence mismatch fails fast with explicit diagnostics.
3. Existing validators and targeted tests still pass.

Acceptance:
1. Infrastructure exists without altering estimator math.
2. No route regressions.

## Phase 3: Migrate `npregbw` LL/LP CVLS/CVAIC First
Scope:
1. Move highest-risk CV route to typed opcode path.
2. Ensure strict collective-order symmetry in session path.
3. Remove legacy fallback for migrated calls only.
4. Excise known `O(n^2)` helper-path usage from migrated LL/LP CVLS routes.

Validation:
1. Historical session hang repro no longer hangs.
2. Attach/profile behavior preserved.
3. Numerical parity:
   - objective value parity (`fval`) where applicable
   - `num.feval` parity/expected behavior
4. Decision-tier performance gate:
   - interleaved paired `times=25` screening
   - escalate to `times>=100` when not obvious
5. Static helper-excision gate on migrated paths:
   - `rg -n 'np_reg_cv_ls_stable_ll_glp|stable_ll_glp|O\\(n\\^2\\)' /Users/jracine/Development/np-npRmpi/src /Users/jracine/Development/np-npRmpi/R`
   - only acceptable matches are comments/archival notes, not active CV route logic.

Acceptance:
1. Session route becomes stable for LL/LP CVLS/CVAIC path.
2. No numerical regression beyond declared tolerance.
3. No active `O(n^2)` helper dependency remains in the migrated LL/LP route.

## Phase 4: Migrate Remaining Core Families
Scope:
1. Migrate sensitive paths in:
   - `npudens*`
   - `npudist*`
   - `npcdens*`
   - `npcdist*`
   - `npscoef*`
   - `npindex*`
   - `npplreg*`
2. Keep migration opcode-by-opcode.

Validation:
1. Route matrix pass for each migrated family.
2. Family-level parity/performance evidence retained per tranche.

Acceptance:
1. No migrated hot path relies on legacy asymmetric dispatch.

## Phase 5: Decommission Legacy Asymmetry in Hot Paths
Scope:
1. Remove free-form worker eval for estimator/CV hot paths.
2. Keep only controlled opcode execution for migrated families.
3. Finalize docs and route guidance.

Validation:
1. Full route matrix pass:
   - `np` serial anchor (`np-master`)
   - `npRmpi` session/spawn (`nslaves=1` minimum + `nslaves>1` where relevant)
   - `npRmpi` attach (`mpiexec -n 2`)
   - `npRmpi` profile/manual-broadcast (`mpiexec -n 2`)
2. Demo triplet fail-fast parity checks pass.

Acceptance:
1. Session is operationally SPMD-equivalent to attach/profile for migrated hot paths.

## Phase 6: Extend to Remaining Non-Core Function Families
Scope:
1. Migrate additional MPI-sensitive non-core/orchestrator/helper routes not in the core family set.
2. Preserve orchestrator boundary policy:
   - orchestrators remain local unless whole-wrapper dispatch is explicitly route-validated.
3. Keep one-risk-axis micro-tranches and retain route-safe behavior in attach/profile/session.

Validation:
1. Route matrix pass for each touched non-core family.
2. Numerical parity and condition-contract stability for each touched API.
3. Performance gate applied where touched path is hot/iterative.

Acceptance:
1. Touched non-core routes are SPMD-safe under session/attach/profile without masked fallback behavior.

## Gate Pack Required At Every Phase
1. Route gate:
   - session, attach, profile/manual validators.
2. Numerical gate:
   - objective and bandwidth parity (or bounded/justified drift).
3. Performance gate:
   - paired interleaved seeds,
   - declared MEI (absolute + relative),
   - mean + median + tail + spread checks,
   - multiplicity-adjusted beneficiary superiority + unaffected equivalence.
4. Demo/test gate:
   - demo masking scan clean,
   - targeted tests pass,
   - subprocess timeout guards used where appropriate.
5. Cleanup gate:
   - no orphan `slavedaemon.R`/`Rslaves.sh`.
6. Complexity gate:
   - no active introduction of `O(n^2)` helper storage/paths in estimation/CV internals.

## Commit and Rollback Policy
1. Commit only after full phase gates pass.
2. Commit message must include:
   - changed scope,
   - validated routes/methods,
   - artifact path.
3. Immediate rollback of tranche if any hard gate fails:
   - hang/deadlock,
   - numerical drift beyond tolerance,
   - performance regression beyond MEI,
   - condition-message regression not explicitly accepted.

## Immediate Next Action
1. Continue Phase 6 non-core migration tranche-by-tranche for remaining autodispatch/manual-dispatch paths (notably `npsigtest` and any uncategorized call heads found by static scan).
2. Add dedicated source-loaded subprocess coverage for new non-core locked opcodes to complement skip-on-CRAN unit guards.
3. Keep route gate (`session`, `attach`, `profile/manual`) and orphan-cleanup checks mandatory on every tranche.
