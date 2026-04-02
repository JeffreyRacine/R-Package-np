# `np-master` Audit Remediation Plan

Date: 2026-04-02
Scope: `/Users/jracine/Development/np-master`
Status: proposed

## Objective

Address the audit findings in `np-master` with the highest probability of eliminating the real risks while minimizing regression, behavioral drift, and collateral damage to estimator correctness, package interfaces, and release hygiene.

This plan is explicitly conservative:

1. one narrow risk axis per tranche,
2. source-loaded proof before installed/tarball proof,
3. no cross-repo mutation without explicit approval,
4. no default drift,
5. no silent semantic remapping,
6. no broad refactor justified only by code style,
7. no live-repo code mutation until a detached validation workspace has re-earned keep status.

## Audit Findings To Remediate

1. Licensing conflict risk:
   - `DESCRIPTION` declares `License: GPL`.
   - `src/nr.c` contains a Numerical Recipes notice that appears incompatible with GPL redistribution.
2. Exported `npseed()` accepts invalid inputs:
   - `NA`, `Inf`, and coercion edge cases are not rejected before reaching the C RNG path.
3. Unchecked native allocations:
   - user-reachable paths in `src/jksum.c` dereference `malloc` results without always checking for `NULL`.
4. Overflow-prone allocation sizing:
   - several allocation sites multiply signed counts before `malloc`, risking wraparound and undersized buffers.
5. Error-path cleanup leakage:
   - helper functions allocate multiple heap blocks and may call `error()` after partial allocation, leaking memory for the life of the R session.
6. Coverage gaps:
   - many `.Rd` examples remain under `\\dontrun{}`.
   - some tests depend on optional packages such as `crs` and therefore do not always run in the default local/CI lane.

## First-Pass Remediation Plan

### Phase 1: Fix the legal blocker

1. Inventory every Numerical Recipes-derived function in `src/nr.c`.
2. Determine which are still used by package-shipped runtime paths.
3. Replace the used routines with GPL-compatible implementations.
4. Update `inst/COPYRIGHTS`, `DESCRIPTION`, and any provenance notes.

### Phase 2: Harden `npseed()`

1. Add strict R-side validation:
   - scalar,
   - finite,
   - non-missing,
   - integer-representable.
2. Mirror the validation in `C_np_set_seed`.
3. Add tests for valid and invalid seeds.
4. Update `man/npseed.Rd`.

### Phase 3: Sweep native allocation safety

1. Add checked allocation helpers.
2. Replace raw `malloc`/`calloc`/`realloc` patterns across `src/jksum.c`, `src/tree.c`, `src/linalg.c`, and `src/mat_vec.c`.
3. Add overflow-safe size multiplication helpers.
4. Replace transient heap allocations with `R_alloc` where feasible.

### Phase 4: Improve cleanup discipline

1. Refactor partial-allocation code paths to use a single cleanup block.
2. Eliminate `error()` calls that bypass owned heap cleanup.
3. Re-test touched native entry points.

### Phase 5: Improve coverage

1. Convert safe `.Rd` examples from `\\dontrun{}` to runnable examples.
2. Add or strengthen tests for native error branches.
3. Add a CI lane or local gate with optional dependencies installed.

## Critique Of The First-Pass Plan

The first-pass plan is directionally correct but not yet good enough by the highest engineering standards.

### 1. It is too broad in the native-code phases

"Sweep native allocation safety" bundles four different risk axes:

1. introducing helper APIs,
2. changing allocation semantics,
3. changing overflow behavior,
4. changing cleanup strategy.

That is too much blast radius for a package whose hot paths are in native code and whose correctness matters more than aesthetic cleanup. A broad sweep would maximize merge size, reduce attribution when a regression appears, and make rollback harder.

### 2. It treats the licensing issue like an ordinary code fix

The `nr.c` licensing issue is not just a patch task. It is a release-governance and provenance blocker. A code-first replacement before a usage and compatibility inventory risks unnecessary churn, and a partial replacement could still leave incompatible code in-tree. The plan needs an explicit stop/go decision point before source mutation.

### 3. It does not separate user-facing risk from infrastructure cleanup

`npseed()` is an exported surface with a likely low-regression, high-value fix. Native allocation infrastructure work is broader and riskier. They should not share a tranche or even the same acceptance logic.

### 4. It does not define user-reachable hotspot order

`src/jksum.c` is large. A credible plan cannot say "replace raw allocation patterns across jksum.c" without prioritizing by entry point and reachability. The plan should begin with the exported entry points most likely to exercise unsafe allocations, not with a file-wide mechanical sweep.

### 5. It under-specifies validation

For native safety work, "re-test touched native entry points" is too vague. We need:

1. source-loaded seam proofs for each touched entry point,
2. installed package proofs,
3. tarball-first confirmation at the end,
4. targeted malformed-input/error-path tests,
5. compiler-assisted evidence where feasible.

Without that, the plan leaves too much room for optimistic interpretation.

### 6. It overvalues `R_alloc` as a universal remedy

`R_alloc` is useful for request-scoped temporary buffers, but it is not automatically the right replacement everywhere. Some buffers are cached, shared across helper calls, or intentionally released by package-owned teardown logic. Blind conversion would risk lifetime bugs and performance regressions.

### 7. It does not distinguish "fix now" from "raise coverage later"

Coverage and documentation improvements are worthwhile but should not delay the real safety and release blockers. They belong after the highest-risk correctness and provenance items are controlled.

### 8. It lacks explicit no-regression boundaries

For a package like `np`, every native hardening tranche must preserve:

1. estimator values,
2. bandwidth decisions,
3. user-facing errors/warnings/messages unless intentionally changed,
4. performance in hot paths unless the user explicitly accepts a tradeoff.

The initial plan implies this but does not operationalize it.

## Reformulated Plan

This is the recommended plan.

## Phase -1: Live Repo Protection And Promotion Gate

Goal: ensure the live `np-master` tree is not used as the experimental proving ground for native hardening.

Rules:

1. all implementation happens first in a detached worktree or scratch clone rooted under `/Users/jracine/Development/tmp`,
2. the live repo may receive code only after the detached workspace passes its full tranche gate,
3. promotion to the live repo must use the exact kept patch from the detached workspace, not a fresh hand-reimplementation,
4. if a tranche is inconclusive, regressing, or only partially validated, it does not get promoted.

Minimum promotion gate before any live-repo code update:

1. source-loaded seam proof in the detached workspace,
2. installed-package proof for the named sentinels in the detached workspace,
3. numerical parity evidence on named sentinels,
4. no unexpected hot-path slowdown on the touched route,
5. no new warnings/notes/errors in tranche-level installed validation,
6. retained raw artifacts under `/Users/jracine/Development/tmp`.

Final package-level gate before any "ready" claim:

1. `R CMD INSTALL` + relevant installed tests,
2. final tarball-first `R CMD check --as-cran`,
3. comparison against accepted warning/note baseline.

Why this is safer:

1. the live repo stays rebuild-safe,
2. implementation and promotion are decoupled,
3. only evidence-backed patches reach the working tree the user may rely on.

## Phase 0: Evidence Freeze And Scope Lock

Goal: create a stable baseline before any mutation.

1. Record the current `np-master` head, test baseline, and audit findings in a dated artifact directory under `/Users/jracine/Development/tmp`.
2. Build a function-level inventory for:
   - `src/nr.c` exported/internal consumers,
   - exported R functions that reach the touched native paths,
   - the specific `jksum.c` entry points that are user-reachable.
3. Define the named sentinel set for this campaign:
   - one simple sentinel per touched family,
   - one representative heavier sentinel per touched family.
4. Do not edit code in this phase.

Exit criteria:

1. exact list of touched runtime routes,
2. named validation set,
3. artifact location recorded.

## Phase 1: Licensing/Provenance Decision Tranche

Goal: resolve the `nr.c` blocker with the least risky path.

1. Perform a precise provenance inventory:
   - which routines in `src/nr.c` are still compiled and used,
   - whether any compatible relicensing evidence already exists outside the file header,
   - whether package-shipped distribution under GPL remains defensible.
2. If licensing is clearly incompatible:
   - stop release-readiness claims immediately for the package,
   - choose one of two explicit paths:
     - preferred: replace only the actually-used Numerical Recipes routines with GPL-compatible implementations behind behavior-locking tests,
     - fallback: if replacement scope is larger than one safe tranche, document the blocker and do not mutate unrelated code under the guise of progress.
3. Do not mix this tranche with allocation or RNG changes.

Why this is safer:

1. it isolates a governance blocker from engineering cleanup,
2. it avoids a speculative full-file rewrite,
3. it preserves rollback clarity.

Exit criteria:

1. clear legal disposition documented,
2. if replacement is required, a bounded replacement map exists before implementation begins.

## Phase 2: Exported `npseed()` Hardening

Goal: remove a user-reachable correctness/safety defect with minimal blast radius.

Implementation:

1. Add strict R-side checks in `npseed()`:
   - length exactly `1`,
   - not `NULL`,
   - not `NA`,
   - finite,
   - integer-representable.
2. Add matching C-side checks in `C_np_set_seed` so direct `.Call` misuse also fails fast.
3. Keep accepted semantics unchanged for valid seeds.
4. Update `man/npseed.Rd` to describe rejected inputs explicitly.

Validation:

1. source-loaded tests for:
   - valid scalar integer,
   - numeric scalar convertible to integer exactly if that remains supported,
   - `NA`,
   - `NaN`,
   - `Inf`,
   - length `0`,
   - length `>1`.
2. installed-package confirmation.
3. no other package behavior changes in this tranche.

Why this is safer:

1. exported surface,
2. easy to test exhaustively,
3. tiny rollback surface,
4. little chance of performance collateral damage.

## Phase 3: Allocation-Overflow Guard Infrastructure

Goal: add the minimum shared infrastructure needed for later native hardening without changing many call sites yet.

Implementation:

1. Introduce a small internal helper layer for:
   - checked multiplication to `size_t`,
   - checked `malloc`/`calloc` wrappers for request-scoped allocations.
2. Keep the helper API narrow and local to native code.
3. Do not migrate the entire tree in this tranche.

Validation:

1. compile cleanly,
2. unit-style tests where reachable,
3. no estimator-path behavior changes.

Why this is safer:

1. the helper contract is stabilized before broad adoption,
2. later tranches become smaller and more mechanical,
3. performance impact can be measured once.

## Phase 4: User-Reachable Native Hotspot Hardening By Entry Point

Goal: fix real safety risk where it matters first.

Order:

1. `C_np_kernelsum`
2. `C_np_regression`
3. `C_np_density`
4. conditional density/distribution entry points only after the simpler routes are green

For each entry point, run a separate tranche:

1. replace unchecked allocations only in the call graph needed for that entry point,
2. apply overflow-safe sizing only to that slice,
3. preserve allocation lifetime semantics unless there is a proven leak/error-path defect,
4. do not opportunistically rewrite neighboring code.

Validation for each tranche:

1. source-loaded targeted seam proof,
2. explicit malformed-input or forced-error-path test where feasible,
3. installed-package sentinel proof,
4. numerical parity check on named sentinels,
5. no accepted performance regression in the touched route.

Why this is safer:

1. small diffs,
2. route-level attribution,
3. fast rollback,
4. avoids turning `jksum.c` into a monolithic refactor.

## Phase 5: Error-Path Cleanup Conversion

Goal: remove leak-prone `malloc` + `error()` patterns where ownership is clear.

Rules:

1. convert only truly request-scoped scratch buffers to `R_alloc`,
2. use explicit cleanup labels for owned heap objects with nontrivial lifetime,
3. do not convert static caches or stateful buffers blindly,
4. do not combine cleanup conversion with algorithmic changes.

Initial targets:

1. `src/mat_vec.c`
2. `src/linalg.c`
3. `src/tree.c`

Only after those are green should similar patterns in `jksum.c` be revisited.

Why this is safer:

1. these files are smaller and easier to reason about than `jksum.c`,
2. they provide a pattern library for later migration,
3. they reduce longjmp leak risk without destabilizing the heaviest kernels first.

## Phase 6: Coverage And Documentation Follow-Through

Goal: raise long-term confidence after the real blockers are controlled.

Implementation:

1. audit the `46` `\\dontrun{}` pages and classify them:
   - truly too heavy for examples,
   - convertible to a smoke example,
   - better moved to vignette/test coverage.
2. add a non-default lane that installs optional dependencies such as `crs`.
3. reduce skip-driven blind spots for the touched audit routes first.
4. add regression tests for any newly hardened error paths.

Why this is safer:

1. it avoids delaying safety work,
2. it improves future confidence using evidence from completed fixes,
3. it keeps default-lane runtime growth intentional.

## Non-Negotiable Acceptance Rules

For every implementation tranche above:

1. validate `np-master` first and only `np-master` unless explicit approval is given to port elsewhere,
2. keep one risk axis per tranche,
3. no positive "fixed" claim until the named installed sentinel set passes,
4. no release-ready claim until final tarball-first `R CMD check --as-cran`,
5. if source-loaded and installed results disagree, stop and resolve that mismatch first,
6. if any hardening change regresses hot-path performance beyond noise, replace it with a leaner equivalent or revert it.

## Recommended Execution Order

1. Phase 0
2. Phase 1
3. Phase 2
4. Phase 3
5. Phase 4 entry-point tranches one by one
6. Phase 5
7. Phase 6

## Additional Refinements For Maximum Safety

These refinements further optimize for success probability and minimum collateral damage.

### A. Detached-Worktree Execution Rule

All implementation tranches after Phase 0 should begin in a detached worktree or equivalent scratch copy rooted under `/Users/jracine/Development/tmp` until that tranche re-earns keep status.

Why this is safer:

1. it preserves the live repo as rebuild-safe,
2. it makes rollback trivial,
3. it prevents partially-hardened native work from contaminating unrelated local validation.

### B. Reproducer-Before-Fix Rule

No native safety tranche should start with code edits unless at least one of the following exists first:

1. a concrete malformed-input or forced-failure test for the defect class,
2. a static proof tied to a user-reachable call path and a specific allocation site,
3. a compiler/runtime diagnostic captured in artifacts from the touched route.

Why this is safer:

1. it prevents speculative cleanup,
2. it keeps each change tied to evidence,
3. it improves post-fix confidence because the pre-fix failure mode is known.

### C. Tranche Diff-Size Budget

For Phases 2 through 5, each tranche should target:

1. one primary runtime file,
2. at most one small shared helper addition when truly necessary,
3. at most one documentation/test file set needed to validate the tranche.

If a proposed change exceeds that shape, split it before implementation.

Why this is safer:

1. smaller reviews,
2. better attribution,
3. lower merge conflict risk,
4. easier revert boundaries.

### D. No Mixed Semantic And Safety Changes

A safety tranche must not also:

1. alter estimator algebra,
2. change kernel defaults,
3. change bandwidth ranking behavior,
4. revise progress-stream design,
5. restructure object layouts for convenience.

If a safety change appears to require semantic movement, stop and split the problem.

Why this is safer:

1. safety evidence stays interpretable,
2. parity failures have fewer candidate causes,
3. user-visible drift is less likely.

### E. Native Hardening Order Inside A Route

When hardening a route, apply this order:

1. add overflow checks,
2. add NULL guards,
3. add cleanup-path discipline,
4. only then consider replacing heap scratch with `R_alloc`.

Do not start with lifetime changes.

Why this is safer:

1. overflow and NULL-guard fixes are usually behavior-preserving,
2. cleanup and lifetime changes are more invasive,
3. it front-loads the highest-value, lowest-regression edits.

### F. Dual Validation For Native Tranches

Each native tranche should have both:

1. behavioral proof:
   - source-loaded seam test,
   - installed sentinel proof,
   - numerical parity,
2. structural proof:
   - static scan confirming the targeted unsafe pattern is removed from the touched slice,
   - where feasible, a diagnostic run under stricter compiler/runtime instrumentation in scratch only.

Instrumentation examples may include:

1. extra warning flags,
2. allocator poisoning tools,
3. address/undefined-behavior sanitizers in a scratch-only build,
4. forced low-memory or injected-allocation-failure harnesses where practical.

These are diagnostic aids, not substitutes for package-level validation.

Why this is safer:

1. behavior can stay green while a structural defect remains,
2. structural cleanup can look complete while behavior regresses,
3. using both reduces false confidence.

### G. Tranche Stop Conditions

Stop and revert the tranche immediately if any of the following occurs:

1. named sentinel numerical drift without an already-declared tolerance basis,
2. hot-path slowdown beyond noise on the touched route,
3. new warning/note/error in installed or tarball validation,
4. source-loaded and installed disagreement,
5. evidence that the edited helper is shared by more routes than the tranche was scoped to protect.

Why this is safer:

1. it prevents sunk-cost continuation,
2. it preserves trust in each checkpoint,
3. it forces re-scoping before broader damage accumulates.

### H. `jksum.c` Sub-Slicing Rule

For `jksum.c`, "one route per tranche" is still not enough if the touched route spans multiple internal kernels. Use this sub-slicing order inside the selected route:

1. wrapper/entry allocation layer,
2. one helper cluster,
3. one cached-state or tree-related cluster only if proven necessary.

Do not modify multiple unrelated helper clusters in one `jksum.c` tranche just because they share allocation style.

Why this is safer:

1. `jksum.c` is too large for style-based sweeping,
2. helper-level containment improves blame assignment,
3. performance regressions become easier to localize.

### I. Commit And Claim Discipline

For this campaign:

1. no checkpoint commit should combine more than one tranche,
2. no tranche should be described as "hardened" or "fixed" until installed sentinels pass,
3. legal/provenance work should never be described as resolved by implication from engineering progress,
4. coverage work should not be used to justify unsafe code churn that has not yet re-earned parity.

Why this is safer:

1. commit messages remain trustworthy,
2. user-facing status remains precise,
3. partial progress is less likely to be mistaken for closure.

### J. Endgame Rule

Do not begin Phase 6 broadly until:

1. Phase 1 has a clear disposition,
2. Phase 2 is green,
3. at least one Phase 4 native tranche has successfully completed end-to-end.

This ensures the documentation/coverage expansion is informed by proven repair patterns instead of guesswork.

### K. Checkpoint Cadence Rule

Create explicit checkpoints at these moments:

1. after Phase 0 evidence freeze,
2. after Phase 1 legal disposition,
3. after Phase 2 `npseed()` hardening,
4. after each kept Phase 4 tranche,
5. before any promotion from detached workspace to live repo.

Each checkpoint should include:

1. git commit in the detached workspace for kept changes,
2. a dated artifact note summarizing scope, validation mode, and outcome,
3. the exact base commit used,
4. the exact sentinel set used.

Why this is safer:

1. it prevents long unsnapshotted native work,
2. it makes rollback and resumption reliable,
3. it reduces ambiguity about what evidence belongs to which tranche.

### L. Live-Repo Update Rule

When a detached tranche is ready to promote:

1. refresh the live repo status first,
2. verify the live repo has not drifted in touched files,
3. apply only the kept patch,
4. rerun at least the installed sentinel set on the live repo after promotion,
5. only then consider the tranche truly landed.

If the live repo has drifted in touched files, stop and re-evaluate rather than force-applying.

Why this is safer:

1. it avoids stale patch promotion,
2. it catches merge-time regressions,
3. it respects the "no breakage of the package" goal at the actual landing point.

## What This Plan Explicitly Avoids

1. a file-wide `jksum.c` rewrite,
2. a blind `malloc` to `R_alloc` conversion sweep,
3. mixing legal, RNG, allocation, and documentation work in one tranche,
4. cross-repo porting before `np-master` re-earns proof,
5. performance tradeoffs disguised as safety cleanup,
6. replacing large amounts of native code before first locking behavior with sentinels.
