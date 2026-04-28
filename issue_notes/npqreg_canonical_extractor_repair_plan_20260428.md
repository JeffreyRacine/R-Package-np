# npqreg Documentation, Interface, And Legacy-Path Removal Plan

Date: 2026-04-28

Status: Planning only; no implementation performed here.

## Goal

Remove all legacy `npqreg` quantile-extraction paths according to the one-canonical-route rule, then update the `npqreg` `.Rd` documentation and R-facing function interface so they faithfully expose that single implementation.

The immediate defect is that `npqreg` still contains legacy Powell extraction paths and the documentation/public controls still describe Powell search as though it is the governing quantile-extraction method. The repair must make the canonical one-dimensional quantile solver the only production route. After that, the `.Rd` file must remove Powell references that no longer belong and must document every user-exposed control that remains for the canonical one-dimensional solver.

## Current Working Hypothesis

Subject to confirmation during diagnosis:

- `npqreg()` dispatches to `npqreg.condbandwidth()`.
- `npqreg.condbandwidth()` calls `C_np_quantile_conditional`.
- The default C route constructs a finite support grid, chooses a candidate window, and refines the quantile by minimizing the squared residual between the target quantile level and the estimated conditional CDF.
- Powell-related controls remain in the R formal arguments and `.Rd` text.
- A C-level legacy/Powell route remains selectable through `NP_QREG_EXTRACT_MODE`.
- Powell may also remain as an internal fallback if the candidate route reports failure or an invalid result.

The plan must hold three requirements simultaneously:

1. remove all legacy/Powell extraction branches, environment-selected modes, and production fallback paths from `npqreg`;
2. align documentation and user-facing arguments with the single canonical route;
3. document every control exposed to users for the canonical one-dimensional solver.

## Non-Goals

- Do not redesign `npqreg` numerics during this documentation/interface alignment and legacy-removal campaign.
- Do not replace the current default extractor with a new root solver unless Jeffrey explicitly requests that separate numerical change.
- Do not touch unrelated Powell/NOMAD controls in `npreg*`, `npcdens*`, `npcdist*`, `npscoef*`, `npindex*`, or other estimators where Powell is genuinely part of the live search path.
- Do not leave any legacy `npqreg` quantile-extraction branch reachable after the repair.
- Do not leave Powell references in `npqreg` `.Rd` unless a live `npqreg` production route still uses Powell, which is not the intended end state.
- Do not touch manuscript or simulation campaign files.
- Do not port to `npRmpi` until `np` passes the installed-build gates.
- Do not copy R, C, or `.Rd` files between repos.

## Critique Of The Prior Plan

The prior plan was useful, but not yet sharp enough for this release-stage repair.

Main weaknesses:

1. It did not center the release-facing defect strongly enough: `.Rd` and function arguments must tell users what route they are actually invoking.
2. It correctly identified hidden legacy paths as a problem, but it needed a stricter diagnose-first step to prove whether any sentinel currently relies on them.
3. It did not sharply distinguish user-facing formals from internal implementation scaffolding.
4. It could cause collateral damage by removing arguments that may still be passed through existing examples or tests, without first proving they are unused by the canonical route.
5. It risked broadening from a route-alignment repair into an unrequested numerical-method redesign.

Highest-standard correction:

- Start from what the user sees and can control.
- Audit exactly which arguments affect current `npqreg` output.
- Update `.Rd` wording to be mathematically accurate.
- Remove, reject, or de-document stale Powell-only controls only after tracing whether they are used.
- Remove all legacy `npqreg` extractor branches and environment switches after proving the current canonical route covers the installed sentinels.
- Preserve current default numerical behavior unless a diagnosed mismatch proves a code change is required.

## Reformulated Scientific And Engineering Contract

1. The `npqreg` help page must describe the estimator route the user actually gets by default.
2. The `npqreg.condbandwidth()` formal arguments must not invite users to tune an optimizer that is not governing the default route.
3. There must be one production `npqreg` quantile-extraction path.
4. No environment variable, hidden option, or fallback branch may select a legacy extractor.
5. Any argument retained in the signature must:
   - affect the current route,
   - be needed for documented dispatch or examples,
   - and be documented in the `.Rd` file when it is user-facing.
6. Any argument removed from documentation must also be checked against examples, tests, demos, and reverse internal uses.
7. Ordinary fitted values and summaries from `npqreg` must remain unchanged on deterministic sentinels.
8. Internal cleanup is mandatory for legacy extractor paths, but it must be sequenced after diagnosis so the canonical route is repaired rather than silently replacing a fallback with a new failure.
9. `np` is repaired and validated first; `npRmpi` receives a surgical parity update only after `np` is green.

## Minimal Blast-Radius Rules

This campaign must be surgical. Removing legacy `npqreg` paths must not become a broad optimizer, bandwidth, or quantile-regression refactor.

Hard limits:

1. Touch only the files required for `npqreg` documentation, `npqreg` R argument handling, and the `npqreg` C extraction helper.
2. Do not alter shared Powell, NOMAD, bandwidth-search, density, distribution, regression, single-index, or smooth-coefficient code except for read-only audit.
3. Do not change the general-purpose `powell()` implementation or shared optimizer headers unless diagnosis proves an `npqreg`-only declaration is obsolete.
4. Avoid changing `.Call` signatures or option-vector layouts unless an obsolete `npqreg`-only control cannot otherwise be removed safely.
5. If a lower-risk path exists, prefer leaving internal argument slots reserved but unused over reshaping shared option vectors, provided the public API and documentation do not expose stale controls and no legacy route remains reachable.
6. Keep all helper renames local to `npqreg` and only where they reduce ambiguity; do not perform cosmetic renames in adjacent routines.
7. Preserve fitted-object structure, summary output, plot/predict compatibility, and dispatch behavior unless a diagnosed mismatch proves one is part of the stale interface.
8. Every code edit must have a direct line back to one of these goals:
   - remove a reachable legacy `npqreg` path;
   - remove/reject stale user-facing Powell controls for `npqreg`;
   - document the canonical one-dimensional solver accurately;
   - preserve validation after the first three changes.

Stop condition:

- If removal of the legacy path requires touching shared optimizer architecture or adjacent estimator families, stop and split that into a separately reviewed plan rather than broadening this campaign.

## Phase 1: Baseline And Scope Lock

Create a dated scratch root:

- `/Users/jracine/Development/tmp/npqreg_doc_interface_alignment_20260428`

Record:

- `np-master` and `np-npRmpi` commits and dirty status;
- R version and compiler details;
- installed package provenance;
- exact current `npqreg` formals and `.Rd` usage output.

Primary scope:

- `R/np.qregression.R`
- `man/np.qregression.Rd`
- `src/kernele.c`
- `src/kernelcv.c` if needed for names, comments, or status propagation
- `src/headers.h` if signatures change
- directly required examples/tests for `npqreg`
- the corresponding `npRmpi` files after `np` is green

C files are edited only to remove legacy `npqreg` extractor paths, environment switching, stale comments/names, and obsolete signatures. Unrelated Powell users must not be touched.

Before implementation, write a short scope ledger in the scratch root listing each intended file edit and the exact reason it is necessary. Any file not on that ledger is out of scope unless the ledger is updated with a concrete `npqreg` reason before editing.

## Phase 2: Diagnose The Live Interface And Legacy Reachability

Before patching, produce an argument-use table for `npqreg.condbandwidth()`:

- `tau`
- `ftol`
- `tol`
- `small`
- `itmax`
- `zero`
- `lbc.dir`
- `dfc.dir`
- `cfac.dir`
- `initc.dir`
- `lbd.dir`
- `hbd.dir`
- `dfac.dir`
- `initd.dir`
- any other passed-through `...` values

For each argument, classify it as:

- live canonical extractor control;
- live numerical tolerance/control but poorly documented;
- stale Powell-only control;
- dispatch/compatibility argument;
- unused or misleading argument.

Required evidence:

- R call path from `npqreg()` to `.Call("C_np_quantile_conditional", ...)`;
- C-side use of each option in the default extractor route;
- whether `ftol`, `tol`, `small`, and `itmax` affect the golden-section path, Powell path, or both;
- whether direction-set controls affect ordinary default output;
- whether examples, demos, tests, or package internals pass Powell-only controls to `npqreg`;
- whether `NP_QREG_EXTRACT_MODE` or any other switch can select a legacy path;
- whether any installed sentinel reaches the legacy/Powell fallback under ordinary default usage;
- raw trace for any fallback-triggering sentinel: support bounds, grid best, bracket, refined candidate, objective value, feasibility status, and failure status.
- a symbol/reference map proving the legacy helper being removed is used only by `npqreg`;
- a call-neighborhood map showing adjacent functions in the same C file are not affected by the patch.

No patch is acceptable until the diagnosis explains whether the fallback is dead code or covering a real canonical-route edge case.

## Phase 3: Decide The Minimum Function-Surface Patch

Use the Phase 2 table to choose the narrowest truthful R-level repair.

Preferred outcome:

1. Keep live canonical controls in the function signature.
2. Remove stale Powell-only controls from `npqreg.condbandwidth()` formals if no examples/tests/internal calls require them.
3. Remove stale Powell-only controls from the accepted public interface. If unavoidable transition protection is needed for one release-stage tranche, reject user-supplied non-default values with a precise error and then remove that transition scaffold before final acceptance.
4. Preserve `...` dispatch behavior for `npqreg.formula()` and `npqreg.default()`.
5. Preserve all ordinary numerical output for calls that do not supply stale controls.

Hard constraints:

- Do not silently ignore a user-supplied argument that appears to control optimization. Either it is live and documented, or it is rejected clearly.
- Do not retain a user-selectable or environment-selectable legacy extractor.
- Do not leave production fallback to Powell behind the canonical route.
- Do not leave any final accepted `npqreg` API control undocumented in the `.Rd` file.

## Phase 4: Update `np-master` Documentation

Revise `man/np.qregression.Rd` so the user can understand the actual route.

Required changes:

1. Replace the “Powell Search Controls” subsection with a heading such as:
   - `Quantile Extraction Controls`
   - or `Quantile Extraction Tolerances`
2. Describe the canonical one-dimensional solver plainly:
   - the conditional CDF is estimated at candidate response values;
   - the target quantile is found by the package's one-dimensional residual-solving/refinement route;
   - the exact implementation is described according to the final C route after legacy removal.
3. Document `tau` first among extraction controls.
4. Document every live tolerance/control argument exposed to users for the canonical one-dimensional solver.
5. Remove Powell direction-set descriptions and all other Powell references from `npqreg` if Powell is no longer part of the production route.
6. Preserve the grouped `.Rd` structure and user mental model developed in the documentation campaign.
7. Check examples so they do not mention or imply stale Powell control use.

Avoid overpromising:

- Do not call the implementation a pure root solver unless the final C route confirms that is exactly what it is.
- If the final route remains residual minimization/refinement rather than bracketed root finding, document that accurately while still explaining the one-dimensional quantile-solving problem.

## Phase 5: Update `np-master` Function Interface

Patch `R/np.qregression.R` to match the accepted Phase 3 decision.

Likely changes:

1. Remove stale Powell-only direction controls from `npqreg.condbandwidth()` formals:
   - `lbc.dir`
   - `dfc.dir`
   - `cfac.dir`
   - `initc.dir`
   - `lbd.dir`
   - `hbd.dir`
   - `dfac.dir`
   - `initd.dir`
2. Remove validation and `myoptd` plumbing for controls that are no longer accepted.
3. Keep live tolerance controls only if the C default extractor actually consumes them.
4. Ensure all remaining canonical one-dimensional solver controls are validated, passed intentionally, and documented.
5. Do not change `npqreg()` dispatch, formula handling, or fitted-object structure unless diagnosis proves they are part of the mismatch.

## Phase 6: Remove `np-master` Legacy Extractor Paths

Patch only the `npqreg` C route.

Required C-side changes:

1. Remove `NP_QREG_EXTRACT_MODE` behavior from production source.
2. Remove the legacy/Powell `npqreg` extractor helper if it is not used by any non-`npqreg` route.
3. Remove production fallback from the canonical extractor to Powell.
4. Replace fallback with explicit failure/status propagation if the canonical extractor cannot produce a finite in-support result.
5. If Phase 2 finds a sentinel currently depends on fallback, repair the canonical extractor mechanism directly before removing the fallback.
6. Update helper names/comments where needed so the source describes the single canonical route.
7. Confirm the general-purpose `powell()` implementation remains untouched because other estimators genuinely use it.

Patch discipline:

- Keep diffs small and reviewable.
- Prefer deleting the `npqreg`-specific legacy helper and its local call sites over rewriting surrounding extraction machinery.
- Do not reindent or restyle neighboring C blocks.
- Do not change any shared constants, macros, or optimizer tolerances unless they are proven `npqreg`-local.
- If option-vector slots must remain for ABI/layout stability, mark them internally unused and remove them from the public R interface and `.Rd`.

## Phase 7: Validate `np-master`

Static gates:

- `git diff --check`
- `tools::checkRd()` for `man/np.qregression.Rd`
- `R CMD build` parsing of `.Rd`
- grep gate showing no stale “Powell Search Controls” wording in `man/np.qregression.Rd`
- grep gate showing removed/de-documented arguments are not used in examples/tests/demos
- grep gate showing no `NP_QREG_EXTRACT_MODE` in production source
- grep gate showing no reachable `npqreg` legacy/Powell extractor branch remains
- grep gate showing shared Powell users outside `npqreg` remain intact
- `nm` or source-symbol audit, if useful, showing no unintended exported-symbol churn

Behavioral installed-build sentinels:

1. simple continuous example with deterministic seed and several `tau` values;
2. existing example from `man/np.qregression.Rd`;
3. a formula interface call and an object interface call;
4. summary/print behavior for the resulting `npqregression` object;
5. if stale arguments are rejected, one negative test showing the message is clear;
6. an edge-case sentinel that would have triggered fallback if one exists from Phase 2; it must either pass through the canonical extractor or fail with the new precise error.
7. one adjacent-function smoke sentinel that exercises a nearby non-`npqreg` route but is not expected to change, to catch accidental collateral damage.

Numerical acceptance:

- For ordinary calls, fitted values and summaries match the pre-patch current default within tight tolerance.
- No ordinary sentinel requires the removed legacy path.
- No new warnings/errors occur for documented examples.
- Runtime does not materially regress on sentinels.

Package gates:

- `R CMD build --no-build-vignettes --no-manual`
- `R CMD check --no-manual --ignore-vignettes`

## Phase 8: Port Surgically To `npRmpi`

After `np-master` is green:

1. Repeat the read-only argument-use and legacy-reachability table in `np-npRmpi`.
2. Apply the same R-interface, `.Rd`, and C-route decisions surgically.
3. Preserve any MPI-specific wrapper behavior.
4. Do not introduce runtime `np::` calls.
5. Do not copy files from `np-master`; patch matching sections manually.

## Phase 9: Validate `npRmpi`

Static gates:

- `git diff --check`
- `tools::checkRd()` for the touched `npRmpi` `.Rd`
- grep gates matching `np-master`

Installed-build sentinels:

- `npRmpi.init(nslaves=1)` route for deterministic `npqreg`;
- `npRmpi.init(nslaves=2)` route if feasible;
- same documented example class as `np`;
- clean `npRmpi.quit()`.

Package gates:

- `R CMD build --no-build-vignettes --no-manual`
- `R CMD check --no-manual --ignore-vignettes`

Parity:

- Compare `np` and `npRmpi` results on deterministic sentinels within documented tolerance.
- Confirm no unexpected MPI slowdown or worker idle pattern appears from the route cleanup.

## Phase 10: Release Notes

Update `CHANGELOG` only after both repos pass.

Suggested wording:

- Clarified `npqreg` documentation and argument handling to describe the current quantile-extraction route accurately.
- Removed stale Powell-only `npqreg` controls and legacy quantile-extraction paths so `npqreg` has one documented production route.

Do not claim a new numerical algorithm unless one is actually implemented and validated.

## Acceptance Criteria

This campaign is complete only when:

1. `npqreg` documentation accurately describes the canonical one-dimensional quantile solver.
2. The accepted function signature contains no Powell-only controls.
3. All legacy `npqreg` extractor paths, environment switches, and hidden fallback branches are removed.
4. Every user-exposed control for the canonical one-dimensional solver is documented in `.Rd`.
5. Existing ordinary `npqreg` results are unchanged on deterministic installed-build sentinels.
6. Examples, tests, demos, and help pages are internally consistent.
7. `np` passes source, Rd, build, check, and installed sentinels.
8. `npRmpi` receives the same user-facing and internal-route alignment and passes matching gates.
9. All artifacts are retained under `/Users/jracine/Development/tmp`.

## Suggested Checkpoint Commits

1. `Plan npqreg documentation and legacy-path alignment`
2. `Align npqreg documentation with extractor route`
3. `Align npqreg interface with documented controls`
4. `Remove npqreg legacy extractor paths`
5. `Port npqreg canonical extractor alignment to npRmpi`
6. `Note npqreg canonical extractor alignment in changelogs`
