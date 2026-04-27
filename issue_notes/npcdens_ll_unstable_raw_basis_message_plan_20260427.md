# Plan: Targeted Local-Linear Raw-Basis Failure Message

Date: 2026-04-27

## Goal

Make the existing `npcdensbw()` and, if reachable, `npcdistbw()` fixed-bandwidth
search failures more useful when canonical local-linear conditional routes fail
because the raw local-polynomial basis is numerically unstable.

The user-facing contract should be:

- preserve current estimator semantics;
- add no data pre-scan and no numerical overhead in successful searches;
- add no per-objective, per-candidate, per-multistart, or per-iteration work;
- do not silently switch `regtype = "ll"` to another route;
- when the existing search failure occurs for canonical `ll`, surface a short
  actionable error telling the user to retry the equivalent degree-1
  local-polynomial route with Bernstein basis, or to center/scale continuous
  regressors.

## Evidence Anchor

Patents diagnostic under `/Users/jracine/Development/bkcde_ER/R`:

- `npcdens(y ~ x, data = patent_df, regtype = "ll", bwmethod = "cv.ls")`
  fails with:
  `C_np_density_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective`
- `regtype = "lp", degree = 1, bernstein.basis = TRUE` succeeds.
- `regtype = "lp", degree = 1, bernstein.basis = FALSE` fails with the same
  invalid raw objective.
- Centering or scaling the continuous regressor also succeeds.

Mechanism identified in `np-master`:

- `R/np.condensity.bw.R` canonicalizes `regtype = "ll"` to the LP engine with
  degree 1 and `bernstein.basis.engine = FALSE`.
- The unbounded `npcdens` LP CVLS route in `src/jksum.c` uses raw full-row
  local-polynomial x weights for the quadratic/I1 term.
- For year-valued predictors around 1983-1991, the raw `[1, x]` weighted
  crossproduct can trip the existing `rcond < 1e-10` rejection, producing an
  invalid raw objective.
- The final fixed-bandwidth raw-objective gate in `src/np.c` correctly rejects
  that candidate rather than returning a sentinel.

`npcdist` is less exposed in the ordinary block path because its LP CVLS route
uses dropped/LOO x weights, but it has the same final raw-objective gate and an
all-large fast path that can encounter raw-basis conditioning. Therefore the
same message should be installed only if the same exact failure reaches the
R boundary for canonical `ll`.

## Non-Goals

1. Do not change objective formulas, bandwidth search, estimator output, or
   default routes.
2. Do not make `ll` automatically use Bernstein basis.
3. Do not add a data-condition scan, condition-number preflight, or retry path.
4. Do not alter explicit `lp(..., bernstein.basis = FALSE)` behavior.
5. Do not edit manuscript, campaign, or simulation files.

## Critical Review Of The Initial Plan

The initial plan was intentionally narrow, but it still needed tightening before
implementation.

1. It was too willing to diagnose the exact numerical mechanism from a generic
   final raw-objective error. The patents evidence supports raw-basis
   instability, but the C error text itself only says that the final objective
   is invalid. A robust user message should state the known facts and give the
   proven remedy without claiming more than the failing junction proves.

2. It did not explicitly protect successful-route behavior with pre/post
   installed-build parity. Since the intended change is only error handling,
   successful `lc`, `ll`, and `lp` calls must have identical bandwidths,
   objectives, and selected metadata before and after the patch.

3. It made `.Rd` updates sound optional. They should be required, because the
   error message points users to `lp`, degree 1, and Bernstein basis. The help
   pages should explain that relationship compactly.

4. It proposed exact C-message matching, which is good for collateral-damage
   control, but did not say how to avoid brittle accidental broadening. The
   helper must match a small explicit map of function-specific C messages, not
   a loose regular expression.

5. It did not specify preservation of error behavior outside the targeted
   condition. Tests must prove that unrelated errors, explicit raw-basis `lp`
   failures, and non-`cv.ls` failures are not rewritten.

6. It did not sufficiently separate CRAN-time tests from slower installed
   sentinels. The patents example is valuable as proof but too slow to become a
   routine unit-test dependency unless a very fast minimized reproducer is
   found.

7. It did not make the no-slowdown contract mechanically explicit enough. The
   implementation must not touch C objective functions, optimizer wrappers,
   candidate feasibility checks, or multistart loops. The only acceptable
   runtime hook is a single R condition handler around the top-level `.Call`,
   which runs once per bandwidth-selection call and does no work unless an
   error is thrown.

## Proposed User Message

Use an error, not a warning, because computation has already failed.

For `npcdensbw()`:

```text
npcdensbw() local-linear cv.ls failed while using the canonical raw degree-1 basis. Try regtype = "lp", degree = 1, bernstein.basis = TRUE, or center/scale continuous regressors.
```

For `npcdistbw()`:

```text
npcdistbw() local-linear cv.ls failed while using the canonical raw degree-1 basis. Try regtype = "lp", degree = 1, bernstein.basis = TRUE, or center/scale continuous regressors.
```

Keep the message short. Avoid a long explanation of internal LP canonicalization
in the primary error text. The phrase "while using" is deliberate: it states the
route that failed without overclaiming that every invalid raw objective has the
same numerical cause as the patents example.

## Reformulated Implementation Plan In `np-master`

### 1. Diagnose Before Patching

Create a retained diagnostic bundle under:

`/Users/jracine/Development/tmp/ll_raw_basis_message_20260427`

Record:

- current `np-master` commit and dirty status;
- installed library path used for proof;
- patents `npcdens` failure;
- patents `lp(degree = 1, bernstein.basis = TRUE)` success;
- one small synthetic shifted-x reproducer if it is deterministic and fast
  enough for tests.
- pre-patch successful-route baseline for several small routes whose behavior
  must remain bit-for-bit or tolerance-identical after the message patch.

No source changes until the reproducer and exact error junction are documented.

### 2. Add A Shared R-Side Error Augmentation Helper

Add a small internal helper, likely in `R/util.R`, that receives:

- function name: `npcdensbw` or `npcdistbw`;
- requested canonical spec;
- bandwidth method;
- caught condition.

It should return the enhanced message only when all are true:

- the caught error message is exactly one of a small function-specific map:
  - `C_np_density_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective`
  - `C_np_distribution_conditional_bw: optimizer returned a fixed-bandwidth candidate with invalid raw objective`
- the public request was `regtype = "ll"`;
- `bwmethod = "cv.ls"`;
- the canonical engine is LP degree 1 with `bernstein.basis.engine = FALSE`;
- at least one continuous regressor is present, so the local-linear raw-basis
  advice is meaningful.

Otherwise it must rethrow the original condition unchanged.

This keeps the check at the existing failure boundary and avoids preflight
numerics. The only routine-path cost is a single R condition handler around the
existing `.Call`; no objective, bandwidth, or data matrix is recomputed.

Implementation details:

- implement the helper entirely on the R side unless diagnosis discovers that
  the public `regtype` cannot be recovered there;
- do not edit C objective code, C optimizer loops, C feasibility checks, or
  candidate scoring;
- do not allocate or inspect data matrices in the success path;
- do not use a broad regex such as `"invalid raw objective"`;
- do not rewrite explicit `regtype = "lp"` failures;
- do not rewrite `cv.ml` failures;
- do not rewrite infeasible-candidate or feasible-candidate failures with
  different C messages;
- do not bury the original C message in a long paragraph. If retained, keep it
  as a short parenthetical only after the actionable message.

Hard performance contract:

- no new code may run inside `bwmfunc`, `bwmfunc_wrapper`, Powell/NOMAD
  candidate evaluation, or any C objective helper;
- successful bandwidth searches may pay only the constant R `tryCatch` /
  condition-handler setup cost around the one existing `.Call`;
- if profiling or timing suggests measurable successful-route slowdown beyond
  noise, the tranche fails and must be redesigned.

### 3. Wire `npcdensbw()`

Wrap the existing `.Call("C_np_density_conditional_bw", ...)` in
`npcdensbw.condbandwidth()` with the helper.

Acceptance criteria:

- ordinary successful `npcdensbw()` calls are unchanged;
- non-`ll` failures retain the previous error text;
- explicit `lp(..., degree = 1, bernstein.basis = FALSE)` retains the previous
  error text unless a separate explicit-LP message is later requested;
- patents canonical `ll` failure gets the new short actionable error.
- the resulting condition is still an error, not a warning, retry, fallback, or
  silently altered search.

### 4. Wire `npcdistbw()` At The Same Boundary

Apply the same helper around `.Call("C_np_distribution_conditional_bw", ...)`
in `npcdistbw.condbandwidth()`.

Because `npcdist` may not hit this failure in the patents route, do not force a
slow or artificial estimator-level failure test. Instead:

- test the shared helper directly for the `npcdistbw` C message;
- run installed-build smoke/parity checks proving ordinary `npcdist` `ll` and
  `lp` routes are unchanged;
- add an estimator-level `npcdist` failure test only if a genuinely fast,
  deterministic reproducer exists.

### 5. Tests

Add focused tests in `np-master/tests/testthat`:

- helper preserves unrelated errors;
- helper augments only canonical `ll` + `cv.ls` + matching invalid raw objective;
- helper does not augment explicit `lp` raw-basis failure;
- helper does not augment the wrong function's C message;
- helper does not augment when no continuous regressor is present;
- patents or synthetic `npcdens` failure emits the new message if runtime is
  acceptable for the test suite;
- `lp(degree = 1, bernstein.basis = TRUE)` succeeds on the same data in an
  installed-build sentinel, even if not used as a CRAN-time unit test.

Prefer a small synthetic test over the full patents data if it reproduces the
same failure quickly and deterministically.

Do not put the full patents runtime in ordinary package tests unless it has been
minimized to a fast deterministic case.

### 6. Documentation

Update `.Rd` documentation:

- `npcdensbw.Rd`: mention that canonical `ll` uses the raw degree-1 local
  polynomial basis, and that users with unstable raw continuous regressors can
  use `regtype = "lp", degree = 1, bernstein.basis = TRUE` or center/scale.
- `npcdistbw.Rd`: same note if the option surface and failure mode are present.

Keep wording compact and avoid presenting Bernstein as a silent replacement for
`ll`; it is an explicit user-selected equivalent-degree LP remedy.

Documentation should also make clear that this is a remedy for search failures
or numerical instability, not a new default or an automatic fallback.

### 7. Installed-Build Validation

Build and install `np` into a private library under the diagnostic bundle.

Required gates:

- `R CMD build`;
- private `R CMD INSTALL`;
- focused testthat tests for the new helper/message;
- pre/post successful-route parity for small `npcdens` and `npcdist` examples:
  compare objective values, bandwidths, `regtype.engine`, degree, basis, and
  Bernstein metadata;
- patents installed-build sentinel:
  - `npcdens(..., regtype = "ll", bwmethod = "cv.ls")` fails with the new
    actionable message;
  - `npcdens(..., regtype = "lp", degree = 1, bernstein.basis = TRUE,
    bwmethod = "cv.ls")` succeeds;
- ordinary `lc` patents route still succeeds;
- small `npcdist` smoke route still succeeds for `ll` and `lp`;
- existing scale-factor/search-floor focused tests remain green.
- a lightweight timing smoke on fast successful `npcdensbw` and `npcdistbw`
  routes before/after the patch, with the expected result being no measurable
  change beyond noise. This is a guardrail, not a performance-claim benchmark.

Retain raw logs and summaries under the diagnostic bundle.

Positive status may be claimed only after installed-build proof. Source-loaded
helper tests are not sufficient by themselves.

### 8. Commit Checkpoint

Commit `np-master` only after installed-build validation passes.

Suggested commit message:

`Improve local-linear raw-basis failure guidance`

## `npRmpi` Port Plan

After `np-master` passes and is committed:

1. Inspect the corresponding `npRmpi` R/C boundary.
2. Port only the R-side helper/message wiring and tests needed for parity.
3. Preserve MPI-specific behavior and avoid any runtime `np::` bridge.
4. Install `npRmpi` into a private library.
5. Validate:
   - focused message tests;
   - serial/session smoke route;
   - `nslaves = 1` and, where applicable, `nslaves > 1` smoke coverage for
     touched conditional routes;
   - ordinary successful routes retain prior behavior.
6. Commit `npRmpi` only after its installed-build gates pass.

## Regression Risks And Mitigations

1. Risk: swallowing unrelated C errors.
   Mitigation: match exact C error text and otherwise rethrow unchanged.

2. Risk: implying automatic semantic equivalence between `ll` and all `lp`
   settings.
   Mitigation: recommend only `regtype = "lp", degree = 1,
   bernstein.basis = TRUE` as an explicit retry, and keep documentation clear.

3. Risk: adding runtime overhead.
   Mitigation: no preflight and no numerical diagnostics; only an error handler
   around an already existing `.Call`.

4. Risk: overfitting the message to `npcdens`.
   Mitigation: share the helper but apply function-specific labels; for
   `npcdist`, augment only the identical invalid-raw-objective failure.

5. Risk: tests become slow.
   Mitigation: use helper-level tests and a fast synthetic reproducer if
   available; keep patents as an installed sentinel rather than a default CRAN
   unit test if runtime is material.

## Done Definition

The tranche is complete only when:

- `np-master` has the targeted message, tests, documentation, private install
  proof, and commit;
- `npRmpi` has the surgical parity port, tests, private install proof, and
  commit;
- no estimator output, default route, or successful search behavior changes;
- all retained logs live under `/Users/jracine/Development/tmp`;
- no manuscript or simulation campaign files are touched.
