# npcdens Scale-Factor Floor Repair Plan

Date: 2026-04-26

Scope:
- Primary repo: `/Users/jracine/Development/np-master`
- Parity repo, only after `np-master` is diagnosed and gated: `/Users/jracine/Development/np-npRmpi`
- Baseline reference tarball: `/Users/jracine/Development/Apr_9_pre_claude_forensic_audit_remediation/np_0.70-1.tar.gz`
- Baseline-matching source commit for relevant files: `fe983631`

## Goal

Correctly implement a user-facing continuous scale-factor lower bound for bandwidth search without breaking established `npcdens(..., bwmethod = "cv.ls")` behavior.

The repaired implementation must:

1. Preserve Apr 9 behavior when `scale.factor.lower.bound` is not active or is set to the legacy equivalent.
2. Correctly enforce the requested lower bound during continuous fixed-bandwidth Powell/NOMAD search.
3. Never confuse the default starting scale factor `cfac.init = 0.5` with the search lower bound.
4. Never accept an invalid objective penalty such as `-1e7` as an optimum.
5. Preserve valid local-constant and local-linear `npcdens` behavior on the patent example and on the simulation sentinels.
6. Preserve `np`/`npRmpi` behavioral parity, with `np-master` remaining the source of truth.

## Known Evidence

The Apr 9 tarball does not expose `scale.factor.lower.bound` in `npcdensbw.default()` and does not contain the later bounded-CVLS quadrature/floor-control machinery.

Current `np-master` prints patent-example scale factors equal to:

- `0.5 * sd(x) * n^(-1/6)` for `hx`
- `0.5 * mad(y) * n^(-1/6)` for `hy`

This means the displayed `Scale Factor: 0.5` is consistent with the starting value `cfac.init = 0.5`; by itself it does not prove that the lower bound was reset to `0.5`.

The more serious symptom is the local-linear patent call returning objective `-1e7` very quickly. In current C code, this is consistent with the invalid-penalty sentinel:

- invalid baseline objective,
- invalid doubled fallback objective,
- `bwm_penalty_value = penalty.multiplier * 1e6`,
- default `penalty.multiplier = 10`,
- returned objective `fval = -1e7`.

Most suspicious recent surfaces:

1. `148635c0` (`Remediate bandwidth floor and nmulti contracts`)
   - changed `nmulti` validation from nonnegative to positive,
   - changed conditional bandwidth search to always set `iMultistart = IMULTI_TRUE`,
   - changed the C history length from `as.integer(max(1, nmulti))` to `as.integer(nmulti)`,
   - added stricter C-side `nmulti` checks,
   - expanded fixed-bandwidth feasibility/floor enforcement.
2. `840394a0` / `be45c280` / `f5685926`
   - introduced and propagated scale-factor lower-bound controls.
3. `7429b605` and following bounded-CVLS commits
   - added floor/feasibility behavior and bounded CVLS quadrature paths that can change invalid-objective behavior.

## Non-Negotiable Engineering Rules

1. Diagnose before patching.
2. Do not bypass Powell, NOMAD, or feasibility checks merely to make the symptom disappear.
3. Do not change public defaults or estimator semantics without explicit acceptance.
4. Do not let an invalid sentinel become a selected optimum.
5. Keep repair axes isolated:
   - floor contract,
   - invalid-objective handling,
   - `nmulti` semantics,
   - bounded-CVLS quadrature route,
   - patent/manuscript application behavior.
6. Use installed-build proofs for positive claims.
7. Port to `npRmpi` only after `np-master` is locally repaired, installed, and gated.

## Initial Diagnostic Plan

### 1. Establish A Clean Repro Matrix

Create a persistent diagnostic workspace under:

`/Users/jracine/Development/tmp/npcdens_scale_floor_repair_20260426`

Run the same scripts against:

- Apr 9 installed tarball,
- current `np-master`,
- selected source commits or detached worktrees during bisection.

Primary repros:

1. Patent local-constant:
   - `npcdens(y ~ x, bwmethod = "cv.ls", data = patent_df)`
2. Patent local-linear:
   - `npcdens(y ~ x, bwmethod = "cv.ls", regtype = "ll", proper = TRUE, data = patent_df)`
3. Small synthetic continuous `x`, continuous `y`, fixed seed, unbounded support.
4. One bounded-response CVLS sentinel from the failed production campaign:
   - chi-square BPA empirical `n = 400` failure case.
5. One Gaussian failed production case:
   - Gaussian BPA empirical or known `n = 800` failure case.

For each run, record:

- package path and version,
- source commit or tarball SHA,
- call,
- elapsed time,
- `fval`, `ifval`,
- `num.feval`, `num.feval.fast`,
- `eval.history`,
- `invalid.history`,
- `fval.history`,
- `xbw`, `ybw`,
- printed scale factors,
- requested `scale.factor.lower.bound` if available,
- computed physical floor values for x/y,
- whether the final candidate is feasible under the declared floor.

### 2. Bisect The Regression Window

Use `fe983631..HEAD` and test only the focused repro matrix.

Decision target:

- identify the first commit where patent LL changes from a real optimization to the rapid `-1e7` sentinel behavior,
- identify whether that same commit also affects production BPA failure cases,
- separate a runtime break from expected bounded-CVLS changes.

Do not infer root cause from commit names. Use run evidence.

### 3. Audit The Floor Contract In R And C

Trace one continuous fixed-bandwidth `npcdensbw(..., bwmethod = "cv.ls")` call from R to C:

1. R default/formal resolution:
   - `npcdensbw.default()`
   - `npcdensbw.conbandwidth()`
2. Effective floor calculation:
   - scale-factor mode (`bwscaling = TRUE`): floor is scale-factor space,
   - physical bandwidth mode (`bwscaling = FALSE`): floor is `scale.factor.lower.bound * EssDee(.) * nconfac`.
3. C argument packing:
   - `myoptd[CBW_SFLOORD]`
   - history length,
   - `iMultistart`,
   - `iNum_Multistart`.
4. C feasibility checks:
   - valid scale-factor check,
   - fixed-continuous floor check,
   - transformed/untransformed candidate handling.
5. Objective handling:
   - raw objective finite/infinite,
   - invalid penalty,
   - best-candidate selection,
   - final objective returned to R.

Expected contract:

- `cfac.init = 0.5` remains a start, not a floor.
- `lbc.init` is only a lower initialization range/control, not a substitute for the user floor.
- `scale.factor.lower.bound = 0.1` is the default lower bound when the new option is active.
- `scale.factor.lower.bound = 0.01` remains valid for legacy-style comparison.
- a candidate below the requested floor is infeasible.
- if all candidates are invalid, the function must fail clearly or keep searching within valid bounds; it must not silently return the invalid sentinel as a fitted optimum.

### 4. Audit `nmulti` Semantics Separately

The previous Apr 9 path allowed `nmulti = 0` at the R level while passing `max(1, nmulti)` as the C history length and setting `iMultistart = IMULTI_FALSE`.

Current code validates `nmulti` as positive and sets `iMultistart = IMULTI_TRUE`.

Decide and document the intended contract:

- if `nmulti = 0` means "no random multistart but one deterministic Powell search", preserve that semantics explicitly;
- if public `nmulti` must be positive, prove no internal hot-start path relies on `0` to disable multistart;
- keep C history allocation separate from the meaning of multistart.

Test cases:

- default `nmulti`,
- explicit `nmulti = 1`,
- internal NOMAD/Powell hot-start path that currently uses "disable_multistart",
- any path that previously used `nmulti = 0`.

### 5. Design The Patch Only After Diagnosis

Candidate repair directions, to be chosen only after the above evidence:

1. Restore a distinct "single deterministic Powell search" state:
   - allow internal no-random-multistart semantics without passing an invalid history length,
   - avoid forcing `iMultistart = IMULTI_TRUE` when the intended mode is no random restarts.
2. Keep the floor check but make final-candidate acceptance explicit:
   - accept only finite objective and feasible candidate,
   - never mark a penalty sentinel as an optimum,
   - emit a clear diagnostic error if no feasible finite candidate exists.
3. If transformed bounds are active, ensure every feasibility check is applied in the same scale as the candidate being checked.
4. Keep invalid penalties as optimizer guidance only, not as accepted final objectives.
5. Preserve old unbounded convolution CVLS behavior for HRL/FYT-like calls unless the call explicitly requests bounded support/quadrature behavior.

### 6. Validation Gates

Run gates in order.

Smoke gates:

- patent LC and LL complete with finite non-sentinel objective values,
- forcing `scale.factor.lower.bound = 0.1` and `0.01` changes only admissible floor behavior, not start-value identity,
- no `-1e7` selected objective in successful fits.

Contract gates:

- existing `testthat` focused floor tests,
- existing `npcdensbw cv.ls` bounded/unbounded tests,
- new regression tests for:
  - patent-style local-linear CVLS not accepting invalid sentinel,
  - `cfac.init = 0.5` not being mistaken for the lower bound,
  - `scale.factor.lower.bound` respected in final fixed-bandwidth candidates,
  - internal hot-start/no-random-multistart semantics if retained.

Production sentinels:

- chi-square BPA empirical `n = 400` representative failed seed,
- Gaussian BPA failed seed,
- small all-DGP smoke covering BPA/HRL/FYT, known/empirical, at least one `n`.

Installed-build gates:

- build tarball,
- install into private library,
- rerun focused repro matrix from installed package,
- run focused `testthat` files from installed package,
- only then claim "validated".

Release gate:

- `R CMD build`,
- `R CMD check --as-cran`,
- retain logs under the diagnostic workspace.

Parity gate for `npRmpi`:

- port only the diagnosed patch,
- do not copy whole files,
- run matching focused tests under serial and MPI-aware setup as appropriate,
- verify no runtime `np::` bridge dependency is introduced.

## Critique Of The Initial Plan

1. The initial plan is broad enough to catch the likely problem, but it risks conflating three changes: floor control, invalid-penalty handling, and bounded-CVLS quadrature changes.
2. Bisection across every post-Apr-9 commit could be slow if each patent LL run takes minutes; a smaller synthetic sentinel is needed for fast iteration.
3. The plan says "preserve Apr 9 behavior" but the current code intentionally added bounded CVLS quadrature behavior. The real standard should be "preserve Apr 9 behavior on unaffected unbounded calls while preserving accepted bounded-CVLS fixes."
4. The plan initially treats `scale.factor.lower.bound` as if it is only a `npcdens` concern, but the global floor-control commits touched multiple bandwidth selectors. A patch in shared helpers can create cross-family regressions.
5. The plan does not yet define exact failure behavior when no finite feasible candidate exists. Silently returning a sentinel is wrong, but a hard error can also disrupt simulation workflows unless the harness can capture and backfill deterministically.
6. The plan does not explicitly distinguish optimizer exploration penalties from final candidate acceptance. This distinction is central.
7. The plan does not yet require comparison of `regtype = "lc"` versus `regtype = "ll"` on the same data and same lower-bound settings. The patent example suggests this split matters.
8. The plan says to port to `npRmpi` after `np-master`, but it should also require checking whether the same bug already exists in `npRmpi` and whether users could be running it before the serial patch is ready.
9. The plan could create collateral damage by changing `nmulti` semantics globally if the fix is too broad.
10. The plan lacks a rollback criterion for candidate patches.

## Reformulated Plan

### Phase 0: Freeze Evidence And Define Patch Boundaries

Create:

`/Users/jracine/Development/tmp/npcdens_scale_floor_repair_20260426`

Save:

- Apr 9 tarball SHA256,
- current `np-master` commit,
- current installed package path/version,
- exact patent LC/LL outputs from Apr 9 and current builds,
- exact synthetic sentinel outputs from Apr 9 and current builds.

Patch boundary:

- start in `np-master`,
- do not touch `npRmpi` until `np-master` has a diagnosed candidate,
- do not touch manuscript files or simulation campaign files during package repair,
- do not modify unrelated bandwidth selectors unless evidence shows the shared helper is wrong.
- do not rely on the systemwide R library for diagnosis or claims; every
  compared build must be installed into a private library under the diagnostic
  workspace, with `.libPaths()` set explicitly in the test script.
- do not source-load R files for positive evidence; source-loaded probes are
  allowed only for local debugging and must be labelled as such.
- do not patch `master` directly for exploratory instrumentation; use a
  detached worktree or scratch copy for temporary tracing.

### Phase 1: Build A Fast Sentinel Before Bisection

Construct a small reproducible R script that fails or exposes the same `-1e7` sentinel behavior quickly.

Minimum required cases:

1. Patent subset or reduced patent call if it reproduces.
2. Full patent LL call if no reduced case reproduces.
3. Synthetic unbounded continuous `x,y` LL CVLS call.
4. One bounded BPA chi-square case if needed to cover quadrature.

The fast sentinel must record objective, bandwidths, evaluation counts, and invalid counts.

Only after the sentinel is reliable should bisection begin.

Additional sentinel design rules:

- test `npcdensbw()` directly before testing full `npcdens()` so bandwidth
  selection is separated from estimation and `proper = TRUE` normalization;
- include both `regtype = "lc"` and `regtype = "ll"` on the same data so the
  local-linear-only failure mode is visible;
- include a no-bound/unbounded case and a bounded-response/quadrature case as
  separate sentinel families;
- make the script self-report `find.package("np")`, package version, `.libPaths()`,
  R version, platform, and the SHA/tarball label being tested;
- define failure codes explicitly:
  - `PASS_FINITE`: finite non-sentinel objective and feasible final bandwidths,
  - `FAIL_SENTINEL`: selected objective equals or is effectively the invalid
    sentinel,
  - `FAIL_ERROR_EXPECTED`: clear error from no feasible finite candidate,
  - `FAIL_ERROR_UNEXPECTED`: unrelated error,
  - `FAIL_TIMEOUT`: runtime exceeds a predeclared timeout.

### Phase 2: Bisect With Two Questions

Question A:

Which commit first causes unbounded patent-style LL CVLS to return the rapid invalid sentinel or unchanged initial bandwidths?

Question B:

Which commit first causes the bounded production BPA failures?

Possible outcomes:

- same commit: one root cause likely;
- different commits: repair separately and keep patches isolated.

Implementation rule:

- each tested commit must be built and installed into its own private library;
- if a commit does not build for reasons unrelated to the sentinel, mark it as
  `SKIP_BUILD` and continue using adjacent commits rather than editing that
  commit;
- retain raw logs for every tested commit so the bisection can be audited.

### Phase 3: Diagnose The First Bad Commit At Source Level

For the first bad commit, trace exactly:

- input bandwidth/start vector,
- transformed/untransformed vector,
- floor coefficient,
- physical floor,
- raw baseline objective,
- penalty objective,
- Powell output,
- final accepted vector,
- final returned objective.

Add temporary diagnostic instrumentation only in a detached worktree or scratch copy, not in the live repo commit path.

Root-cause classification must be one of:

1. wrong floor scale,
2. wrong final-candidate feasibility check,
3. invalid penalty accepted as final optimum,
4. `nmulti`/hot-start semantic regression,
5. bounded-CVLS quadrature objective regression,
6. other, with evidence.

Mandatory source audit checks:

- verify every R-side `myopti`/`myoptd` element is aligned with the matching C
  enum/index in `src/headers.h`;
- verify added C arguments and history lengths cannot shift or truncate data
  passed from R;
- verify the raw objective and penalized wrapper objective are not being
  compared interchangeably during final candidate selection;
- verify transformed-bound candidates are converted exactly once before
  objective evaluation and feasibility checks;
- verify `regtype = "ll"` canonical LP handling is not accidentally routed
  through a bounded-quadrature or degree-search path intended only for
  `regtype = "lp"`/NOMAD.

### Phase 4: Patch Only The Proven Root Cause

Preferred patch shape if evidence confirms the current suspicion:

1. Keep `scale.factor.lower.bound` public API and default.
2. Keep floor enforcement in candidate feasibility checks.
3. Restore a clear separation between:
   - deterministic single Powell search,
   - random multistart count,
   - C history allocation length.
4. Treat invalid penalty values as exploration penalties only.
5. Final candidate acceptance must require:
   - finite raw objective,
   - feasible fixed-bandwidth candidate,
   - objective is not the configured invalid sentinel.
6. If no such candidate exists:
   - return a clear error with enough context for reproducibility,
   - do not return a fitted bandwidth object with `fval = -1e7`.

Avoid:

- changing `cfac.init`,
- increasing the default floor to hide failures,
- disabling Powell/NOMAD,
- weakening floor checks,
- broad rewrites of bounded-CVLS quadrature,
- cross-repo file copying.

Patch discipline:

- one patch per root-cause class;
- if the first bad commit is an `nmulti` semantic change, patch only `nmulti`
  semantics and rerun both sentinel families before touching floor code;
- if the first bad commit is an invalid-sentinel acceptance issue, patch final
  acceptance logic without changing floor values;
- if the first bad commit is a quadrature objective regression, patch only the
  bounded-CVLS route and prove unbounded patent-style CVLS is unaffected.

### Phase 5: Add Narrow Regression Tests

Tests must be small and deterministic.

Required tests:

1. Scale-factor floor API:
   - rejects invalid values,
   - accepts `0`, `0.01`, `0.1`, `0.5`,
   - final continuous fixed bandwidths respect the requested floor.
2. Start-versus-floor distinction:
   - with `cfac.init = 0.5` and floor `0.1`, printed/initial `0.5` is not interpreted as lower-bound enforcement.
3. Invalid sentinel rejection:
   - no successful `npcdensbw` object may return `fval = -1e7` from invalid-penalty selection.
4. `nmulti` contract:
   - deterministic single-start behavior works,
   - random multistart behavior works,
   - NOMAD/Powell hot-start path works.
5. Patent-style local-linear sentinel:
   - local-linear CVLS returns a finite non-sentinel objective or a clear intentional error, never a bogus fitted optimum.

### Phase 6: Installed-Build Validation

Use private libraries, not source-loaded code, for claims.

Validate:

- Apr 9 baseline behavior remains documented,
- current repaired build fixes the bad behavior,
- bounded production failure sentinels pass or fail with an improved, diagnosed message,
- all focused `testthat` files pass,
- `R CMD build` passes,
- `R CMD check --as-cran` passes before any release-ready claim.

Validation must include a before/after table with:

- build label,
- package path,
- sentinel family,
- call,
- elapsed seconds,
- `fval`,
- `ifval`,
- `num.feval`,
- `num.feval.fast`,
- `invalid` count,
- `hx`,
- `hy`,
- computed lower floors,
- status code.

### Phase 7: Port And Gate `npRmpi`

Only after `np-master` passes:

1. Port the minimal patch surgically.
2. Preserve `npRmpi` runtime independence.
3. Run focused parity tests:
   - `nslaves = 1`,
   - at least one MPI setup if the touched path reaches MPI code.
4. Confirm the same patent/sentinel behavior where applicable.

### Phase 8: Simulation Backfill Or Relaunch Decision

Only after repaired installed builds pass:

1. rebuild tarballs with SHA256 provenance,
2. install campaign-local libraries,
3. run installed-build proof on all hosts,
4. backfill failed rows or relaunch only if the campaign plan says so,
5. monitor for worker counts and block-file growth.

## Acceptance Criteria

A candidate patch is acceptable only if:

1. Root cause is identified with evidence.
2. Patch touches the smallest necessary surface.
3. Patent LC and LL behavior is no worse than Apr 9 for the intended unbounded call path.
4. `scale.factor.lower.bound` is genuinely enforced and user-selectable.
5. Invalid penalty sentinels cannot be silently selected as successful optima.
6. Existing bounded-CVLS repairs are not regressed.
7. Focused tests pass from an installed build.
8. `R CMD check --as-cran` passes before tarball/release use.
9. `npRmpi` parity is restored with a surgical port.
10. The repair plan can explain, with evidence, why the patent local-linear
    `-1e7` behavior occurred and why the patch prevents it from recurring.

## Rollback Criteria

Reject or revert a candidate patch if:

1. It fixes the patent LL symptom by disabling Powell/NOMAD or weakening validation.
2. It changes default estimator semantics beyond the scale-floor contract.
3. It lets `-1e7` or another invalid sentinel remain a successful fit.
4. It regresses Apr 9-equivalent unbounded CVLS behavior.
5. It regresses accepted bounded-CVLS sentinels.
6. It requires broad unrelated rewrites to pass focused tests.
7. It cannot pass installed-build gates.
