# Semiparametric Scale-Factor Search-Bound Repair Plan

Date: 2026-04-27

Scope:

- Source of truth: `/Users/jracine/Development/np-master`
- Parity target after `np` is validated: `/Users/jracine/Development/np-npRmpi`
- Confirmed affected families: `npindexbw` and `npscoefbw`
- Included policy target: uniform continuous scale-factor search defaults and semantics across core `np*` selectors unless Phase 0 finds a concrete estimator-specific reason to retain a documented exception
- Out of scope for this campaign: broad helper refactors, manuscript files, simulation campaign files, unrelated scale-factor policy changes

## Goal

Repair the remaining semiparametric mismatch between the initialization lower bound and the optimizer search floor.

The scientific contract is:

- Users should be able to expect all core `np*` bandwidth selectors to share a common continuous scale-factor search policy unless a selector-specific exception is mathematically or computationally necessary.
- `scale.factor.init` is the deterministic first continuous scale-factor start.
- `scale.factor.init.lower` and `scale.factor.init.upper` define the random multistart initialization range.
- `scale.factor.search.lower` defines the admissible continuous fixed-bandwidth search floor.
- A user may set `scale.factor.search.lower` below `scale.factor.init.lower` to permit search below the random-start range.
- For generated search starts, `scale.factor.init` and `scale.factor.init.lower` may be raised to the floor by the existing `max(..., scale.factor.search.lower)` policy.
- For explicit user bandwidths used for evaluation rather than search starts, the search floor must not rewrite the user-supplied bandwidth.
- The default continuous scale-factor search policy should be globally consistent across core selectors unless a selector has an explicit, documented mathematical or computational reason to differ.

Default target:

- `scale.factor.init = 0.5`
- `scale.factor.init.lower = 0.1`
- `scale.factor.init.upper = 2.0`
- `scale.factor.search.lower = NULL`, resolved through the existing fallback policy to `0.1`

## Execution Discipline

Use a one-risk-axis-per-tranche workflow:

1. Repair the confirmed semiparametric search-floor bug in `np-master`.
2. Validate and commit that bug repair.
3. Harmonize semiparametric defaults in `np-master`.
4. Validate and commit that default harmonization.
5. Port the bound repair to `npRmpi`.
6. Validate and commit the `npRmpi` bound repair.
7. Port the default harmonization to `npRmpi`.
8. Validate and commit the `npRmpi` default harmonization.

All validation must use private installed libraries rooted under `/Users/jracine/Development/tmp/scale_factor_semiparam_bound_repair_20260427`, not the user's ordinary site library. Source-loaded probes may be used only for diagnosis and must not be used for positive claims.

If any tranche fails a hard gate, stop that tranche, preserve the logs, and repair or revert only that tranche before continuing. Do not proceed to `npRmpi` with an unvalidated `np-master` change.

## Confirmed Defect

Claude's audit correctly identified that two semiparametric paths still use the effective initialization lower bound as a hard optimizer lower bound.

Affected source-truth sites:

- `R/np.singleindex.bw.R`: fixed-bandwidth NOMAD lower bound uses `h.start.controls$scale.factor.init.lower`.
- `R/np.smoothcoef.bw.R`: `optim()`-based fixed-bandwidth admissibility lower bound uses `start.controls$scale.factor.init.lower`.

Affected `npRmpi` mirror sites:

- `R/np.singleindex.bw.R`: the same fixed-bandwidth NOMAD lower-bound pattern.
- `R/np.smoothcoef.bw.R`: the same `optim()`-based fixed-bandwidth admissibility lower-bound pattern.

This is wrong when, for example, `scale.factor.search.lower = 0.1` and `scale.factor.init.lower = 0.5`: the optimizer may be unintentionally prevented from exploring the interval `[0.1, 0.5)`.

## Initial Focused Repair Plan

1. Establish read-only diagnosis in `np-master`.
   - Confirm every affected site with line references.
   - Confirm whether each helper object currently carries `scale.factor.search.lower`.
   - Confirm whether any additional semiparametric bound/admissibility sites still use `scale.factor.init.lower` as a search floor.

2. Add narrow source-level instrumentation or unit probes only if needed.
   - Prefer tiny R-level sentinels that inspect bound construction behavior or force candidate admissibility checks.
   - Avoid long stochastic optimization as the first proof.
   - Save raw logs under `/Users/jracine/Development/tmp/scale_factor_semiparam_bound_repair_20260427`.

3. Patch the confirmed `np-master` bound bug surgically.
   - Thread `scale.factor.search.lower` through the relevant semiparametric start-control return objects.
   - Use `scale.factor.search.lower` for hard optimizer lower bounds and admissibility floors.
   - Preserve `scale.factor.init.lower` for deterministic and random start generation.
   - Do not change `npContinuousSearchStartControls()` soft-clip policy.

4. Validate the `np-master` bound repair from an installed private build.
   - Run focused sentinels showing that `search.lower < init.lower` is allowed for search.
   - Run tightening sentinels showing that `search.lower > init.lower` remains enforced.
   - Run representative `npindexbw` and `npscoefbw` smoke examples.
   - Run adjacency checks for object creation, summary, and any lightweight predict/plot route that is already standard for these families.
   - Run `R CMD check` at the tier appropriate for this tranche before any positive claim.

5. Commit the validated `np-master` bound-repair checkpoint.
   - The commit should mention semiparametric scale-factor search floors.
   - Do not port to `npRmpi` before `np-master` passes the installed-build gates.

6. Harmonize semiparametric defaults in `np-master` as a separate tranche if Phase 0 finds no concrete reason for an exception.
   - Target the global defaults listed above.
   - Update `.Rd` usage and argument text for the touched selectors.
   - Validate default-call behavior and user-specified override behavior from an installed build.
   - Commit separately from the bound repair.

7. Port the bound repair surgically to `np-npRmpi`.
   - Apply the same semantic changes by hand, respecting `npRmpi` route-specific code.
   - Do not copy R files between repos.
   - Preserve MPI-specific wrappers, validation harnesses, and lifecycle behavior.

8. Validate the `npRmpi` bound repair from an installed private build.
   - Run the same focused semiparametric sentinels in serial library mode where applicable.
   - Run `npRmpi` attach/profile gates for touched route families if they are part of the current validation harness.
   - Confirm no lingering MPI processes after the gates.
   - Run `R CMD check` at the tier appropriate for this tranche.

9. Commit the validated `npRmpi` bound-repair checkpoint.
   - The commit should explicitly state parity with the validated `np-master` bound repair.

10. Port, validate, and commit the `npRmpi` default-harmonization tranche separately.
    - Keep this separate from the `npRmpi` bound repair so any default-related regression is attributable.

## Criticism Of The Initial Plan

The initial plan is directionally correct, but it is not yet tight enough for a low-regression repair.

First, the naive audit fix says to replace `start.controls$scale.factor.init.lower` with `start.controls$scale.factor.search.lower`. That replacement is not safe until the helper return objects are audited and extended. At present, the canonical start-control helper returns the effective start controls, not necessarily the raw search floor. A blind substitution could introduce `NULL` lower bounds or inconsistent vector scaling.

Second, the plan needs stronger proof that `scale.factor.init.lower` remains in the initialization role. A repair that makes the search floor correct but accidentally changes multistart generation would create a subtle behavioral drift. The acceptance tests must explicitly check the separation between the random-start bracket and the admissible search bracket.

Third, the plan needs to distinguish fixed-bandwidth search from evaluation more explicitly. The package contract permits user-supplied bandwidths below the floor when they are not search starts. Any sentinel that only checks lower-bound construction could miss a regression where explicit evaluation bandwidths are newly rejected.

Fourth, the plan should avoid broad documentation/default churn in the same tranche. The audit raises real follow-up questions about defaults and `.Rd` wording, but mixing those changes with the bug fix would create a wider review surface and make regressions harder to localize.

Fifth, `npRmpi` validation should not be described generically. Because recent work exposed lifecycle and MPI harness sensitivity, the `npRmpi` phase needs explicit clean-field preflight and post-run process checks.

Sixth, acceptance should include negative controls. The best evidence is not merely that the formerly blocked interval opens, but also that existing tightened floors still reject candidates below the requested floor and that default behavior remains unchanged for representative calls.

Seventh, the first draft treated semiparametric default differences as a follow-up design question. That is too weak given the package goal: one coherent public scale-factor search policy across the core `np*` API. The better framing is to require uniform defaults and semantics unless the inventory uncovers a specific estimator contract that would be harmed. Even then, any exception should be explicit in code comments, documentation, and validation artifacts.

## Reformulated Repair Plan

### Phase 0: Preflight And Bound Inventory

1. Confirm clean working trees for both repos.
2. Create the persistent artifact root:
   - `/Users/jracine/Development/tmp/scale_factor_semiparam_bound_repair_20260427`
3. Record the starting commit SHA, `sessionInfo()`, R version, compiler context, and private-library path for every installed-build gate.
4. In `np-master`, inventory all semiparametric occurrences of:
   - `scale.factor.init.lower`
   - `scale.factor.search.lower`
   - `fixed.lower`
   - `fixed.h.lower`
   - `h.lower`
   - candidate admissibility checks
5. Classify each occurrence as one of:
   - deterministic first-start construction,
   - random multistart bracket,
   - optimizer hard lower bound,
   - candidate admissibility floor,
   - user-supplied evaluation bandwidth handling,
   - documentation or object metadata.
6. Do the same read-only classification in `np-npRmpi` before any porting.
7. Record the classification in the validation artifact root.

Acceptance for Phase 0:

- Every affected semiparametric use of `scale.factor.init.lower` has an assigned role.
- The inventory explicitly states whether any hard-bound/admissibility use of `scale.factor.init.lower` remains legitimate.
- No code is edited before the classification is complete.

### Phase 1: Focused Reproduction Sentinels For `np-master`

Create tiny installed-build sentinels that exercise the contract without relying on long optimization.

Sentinels should prefer deterministic `nmulti = 1` and fixed seeds. If an internal helper must be inspected, use the installed namespace via `pkg:::` from the private library and record that the sentinel is checking internal contract plumbing.

Required sentinel behavior:

1. `npindexbw`, permissive floor:
   - set `scale.factor.search.lower < scale.factor.init.lower`;
   - verify the fixed-bandwidth optimizer lower-bound construction uses the search floor scale, not the initialization lower scale.

2. `npindexbw`, tightened floor:
   - set `scale.factor.search.lower > scale.factor.init.lower`;
   - verify the effective lower bound is the tightened floor.

3. `npscoefbw`, permissive floor:
   - set `scale.factor.search.lower < scale.factor.init.lower`;
   - verify fixed-bandwidth candidate admissibility uses the search floor, not the initialization lower scale.

4. `npscoefbw`, tightened floor:
   - set `scale.factor.search.lower > scale.factor.init.lower`;
   - verify candidates below the tightened floor remain inadmissible.

5. Initialization separation:
   - verify generated deterministic/random starts still respect the effective initialization lower bound and are not generated below the floor.

6. Evaluation separation:
   - verify an explicit user bandwidth below `scale.factor.search.lower` is not rejected merely because it is an evaluation bandwidth rather than a search start.

Acceptance for Phase 1:

- At least one sentinel reproduces the current defect before patching or records a direct source-level trace proving the mismatch.
- The sentinels are narrow enough to run quickly and deterministically.
- Each sentinel records whether its evidence is source trace, installed internal-helper trace, or public API behavior.

### Phase 2: `np-master` Surgical Bound Patch

Patch only the semiparametric bound plumbing needed to repair the defect.

Implementation rules:

1. Extend `.npindexbw_h_start_controls()` to retain `scale.factor.search.lower` in its returned object after validation.
2. Extend `.npscoefbw_start_controls()` to retain `scale.factor.search.lower` in its returned object after validation.
3. Replace hard optimizer lower-bound uses of the initialization lower bound with the retained search floor:
   - fixed-bandwidth NOMAD lower bound in `npindexbw`;
   - fixed-bandwidth Powell/NOMAD handoff admissibility lower bound if inventory shows it is affected;
   - `optim()`-based fixed-bandwidth admissibility lower bound in `npscoefbw`.
4. Preserve `scale.factor.init.lower` everywhere its role is start generation.
5. Preserve existing defaults in this bound-repair tranche only.
6. Preserve existing soft-clip semantics in `npContinuousSearchStartControls()`.
7. Do not touch broad documentation wording in this tranche unless a usage signature is directly wrong for the touched arguments.

Acceptance for Phase 2:

- Diff is limited to the two affected R files plus narrow tests/sentinels if they are kept in-tree.
- No default value changes.
- No unrelated helper consolidation.
- Static post-patch search confirms no semiparametric hard optimizer lower bound or candidate admissibility floor still uses `scale.factor.init.lower`.
- Static post-patch search confirms start-generation sites still use `scale.factor.init.lower`.

### Phase 3: `np-master` Installed-Build Validation

Run validation from a private installed build, not source-loaded code.

Required gates:

1. Focused sentinel suite from Phase 1.
2. Representative `npindexbw` smoke:
   - default call;
   - permissive `scale.factor.search.lower < scale.factor.init.lower`;
   - tightened `scale.factor.search.lower > scale.factor.init.lower`.
3. Representative `npscoefbw` smoke:
   - default call;
   - permissive `scale.factor.search.lower < scale.factor.init.lower`;
   - tightened `scale.factor.search.lower > scale.factor.init.lower`.
4. Lightweight adjacency:
   - object summary for both families;
   - prediction or fitted-object route only where the standard example is stable and fast.
5. `R CMD check` at tranche tier.
6. Static post-patch search assertions from Phase 2.

Acceptance for Phase 3:

- The formerly incorrect lower-bound behavior is corrected.
- Tightened floors remain enforced.
- Default representative calls remain unchanged within expected numerical tolerance.
- Explicit old-value overrides remain accepted.
- Logs and summaries are retained under `/Users/jracine/Development/tmp/scale_factor_semiparam_bound_repair_20260427`.

### Phase 4: Commit `np-master` Bound Repair

Commit only after Phase 3 passes.

Commit message target:

`Repair semiparametric scale-factor search floors`

The commit body should mention:

- `npindexbw`;
- `npscoefbw`;
- search floor versus initialization lower-bound separation;
- installed-build validation artifact root.

### Phase 5: `np-master` Core-Selector Default Harmonization Tranche

Run this only after the bound-repair commit is validated.

Goal:

- Bring `npindexbw` and `npscoefbw` into the same default continuous scale-factor policy as the other core selectors, unless Phase 0 produced a concrete estimator-specific exception.
- Confirm that the rest of the core `np*` selectors already share the same user-facing semantics for the four continuous scale-factor search controls.

Implementation target:

- `scale.factor.init = 0.5`
- `scale.factor.init.lower = 0.1`
- `scale.factor.init.upper = 2.0`
- `scale.factor.search.lower = NULL`

Implementation rules:

1. Change only semiparametric public defaults and directly coupled helper defaults.
2. Preserve explicit user overrides exactly.
3. Preserve object metadata and summary reporting of actual selected bandwidths/scale factors.
4. Update `.Rd` usage blocks and argument text for `np.singleindex.bw.Rd` and `np.smoothcoef.bw.Rd`.
5. Document that these four arguments now share the same default policy as the other core continuous fixed-bandwidth selectors.
6. Do not change categorical initialization controls, direction-set controls, or unrelated `*.init` arguments.

Validation gates:

1. Installed-build default-call smoke for `npindexbw` and `npscoefbw`.
2. Installed-build override smoke confirming the old values can still be requested explicitly:
   - `npindexbw(..., scale.factor.init = 1.059224, scale.factor.init.lower = 0.5, scale.factor.init.upper = 1.5)`;
   - `npscoefbw(..., scale.factor.init = 1.0, scale.factor.init.lower = 0.5, scale.factor.init.upper = 1.5)`.
3. Installed-build sentinel confirming `nmulti = 1` uses the global deterministic first-start policy.
4. Documentation check for touched `.Rd` files.
5. Read-only API audit confirming no other core selector exposes divergent defaults or semantics for these four controls without an explicit documented exception.
6. Before/after ledger distinguishing intended default changes from accidental drift:
   - default calls may change because the default policy intentionally changes;
   - explicit old-value override calls should match the pre-harmonization installed build within the tranche tolerance.
7. `R CMD check` at tranche tier.

Commit message target:

`Harmonize semiparametric scale-factor defaults`

### Phase 6: Surgical `npRmpi` Bound-Repair Port

Only after the `np-master` bound-repair tranche is validated and committed:

1. Re-run clean-field preflight:
   - no live `mpiexec`;
   - no live `hydra_pmi`;
   - no live `slavedaemon.R`;
   - no stale `Rslaves.sh` children tied to validation.
2. Apply the same bound-repair semantic changes by hand in `np-npRmpi`.
3. Preserve MPI-specific code and route structure.
4. Do not copy R files between repos.
5. Re-run the occurrence classification after patching to ensure no semiparametric hard-bound site still uses `scale.factor.init.lower`.

Acceptance for Phase 6:

- The `npRmpi` diff is semantically parallel to the accepted `np-master` diff.
- No internal runtime bridge to `np::` is introduced.
- No MPI lifecycle or validation harness behavior is changed.

### Phase 7: `npRmpi` Bound-Repair Installed-Build Validation

Required gates:

1. Focused semiparametric sentinels from the installed `npRmpi` build.
2. `npindexbw` and `npscoefbw` representative smoke examples in ordinary library mode.
3. Attach/profile gates for touched families if supported by the current validation harness.
4. Clean-field postflight:
   - no lingering MPI worker/slave processes;
   - no stale validation sessions.
5. `R CMD check` at tranche tier.

Acceptance for Phase 7:

- `npRmpi` matches the validated `np-master` behavior for the focused contract.
- MPI-specific gates do not regress.
- Logs and summaries are retained under the same dated artifact root or a clearly named `nprmpi` subdirectory.

### Phase 8: Commit `npRmpi` Bound Repair

Commit only after Phase 7 passes.

Commit message target:

`Align semiparametric scale-factor search floors`

The commit body should mention:

- parity with the validated `np-master` checkpoint;
- focused installed-build validation;
- attach/profile status if run.

### Phase 9: Surgical `npRmpi` Default-Harmonization Port

Only after the `np-master` default-harmonization tranche is validated and committed and the `npRmpi` bound-repair tranche is validated and committed:

1. Re-run clean-field preflight.
2. Apply the default-harmonization changes by hand in `np-npRmpi`.
3. Update only directly touched `npRmpi` `.Rd` usage and argument text.
4. Preserve explicit old-value overrides.
5. Re-run read-only API audit for core-selector default consistency.

Acceptance for Phase 9:

- `npRmpi` default signatures and docs match the validated `np-master` policy, modulo package-name-specific text.
- No MPI route behavior is changed except through the public default values.

### Phase 10: `npRmpi` Default-Harmonization Installed-Build Validation

Required gates:

1. Installed-build default-call smoke for `npindexbw` and `npscoefbw`.
2. Installed-build explicit old-value override smoke.
3. Installed-build `nmulti = 1` deterministic-start sentinel.
4. Documentation check for touched `.Rd` files.
5. Attach/profile gates for touched families if supported by the current validation harness.
6. Clean-field postflight.
7. `R CMD check` at tranche tier.

Acceptance for Phase 10:

- `npRmpi` matches the validated `np-master` default policy.
- Explicit old-value overrides remain accepted.
- MPI-specific gates do not regress.

### Phase 11: Commit `npRmpi` Default Harmonization

Commit only after Phase 10 passes.

Commit message target:

`Harmonize npRmpi semiparametric scale-factor defaults`

## Explicit Non-Goals

Do not do any of the following in this repair campaign:

- change the soft-clip policy for initialization controls;
- refactor duplicated helpers unless the patch cannot be made safely without doing so;
- broaden `.Rd` rewrites beyond the touched scale-factor search-control arguments;
- modify manuscript files;
- modify simulation campaign files;
- change C-side floor machinery;
- change `npcdens`, `npcdist`, `npreg`, `npudens`, or `npudist` unless Phase 0 finds the same specific defect there.

## Follow-Up Items After This Campaign

These should be considered only after the focused repair is validated and committed:

1. Run a broader documentation-only pass clarifying `scale.factor.init`, `scale.factor.init.lower`, `scale.factor.init.upper`, and `scale.factor.search.lower` across all selectors after the semiparametric docs are repaired.
2. Consider helper consolidation where duplicate scale-factor validation helpers remain.
3. Consider defensive hardening around global C scaling state only if a concrete risk or failing sentinel is identified.
