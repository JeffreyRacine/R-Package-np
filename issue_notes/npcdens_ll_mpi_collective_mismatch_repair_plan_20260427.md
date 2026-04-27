# npRmpi `npcdens` Local-Linear MPI Collective-Mismatch Repair Plan

Date: 2026-04-27

Scope: `np-npRmpi`, with read-only comparison against `np-master` as the serial
source of truth. This plan does not authorize edits to `np-master` unless a
separate serial bug is found.

## Problem Statement

The patents-data call

```r
npcdens(y ~ x, data = patent_df, regtype = "ll", bwmethod = "cv.ls")
```

fails quickly and correctly in serial `np`, but the corresponding `npRmpi`
route can run for more than 30 minutes without producing the expected failure.
The current diagnostic bundle is:

```text
/Users/jracine/Development/tmp/ll_raw_basis_message_20260427
```

Key observations:

- `np` serial default patents `ll`: fails in `3.669s`.
- `np` serial `nmulti = 1`: fails in `1.868s`.
- `npRmpi`, `autodispatch = FALSE`, `nslaves = 1`: still running at `30s`.
- `npRmpi`, `autodispatch = FALSE`, `nmulti = 1`: still running at `30s`.
- `npRmpi` forced through `.npRmpi_with_local_regression(...)`: fails in
  `4.037s` with the same local-linear raw-basis guidance as `np`.
- Stack sample with autodispatch on and off shows the coordinator in
  `C_np_density_conditional_bw -> np_density_conditional_bw ->
  np_cv_func_con_density_categorical_ls_npksum -> MPI_Allreduce`, while the
  worker is in the Rmpi slave command loop around `mpi_allgather` /
  `MPI_Allgather`.

NOMAD control-path verification on the same installed builds:

- `np` serial `npcdens(y ~ x, data = patent_df, bwmethod = "cv.ls",
  nomad = TRUE)`: success in `66.058s`.
- `npRmpi`, `nslaves = 1`, `autodispatch = TRUE`: success in `34.012s`.
- `npRmpi`, `nslaves = 2`, `autodispatch = TRUE`: success in `34.645s`.
- `npRmpi`, `nslaves = 1`, `autodispatch = FALSE`: success in `34.839s`.
- All NOMAD control runs returned the same objective and bandwidth summary:
  `fval = 0.02706553`, `xbw = 15060.177765257`,
  `ybw = 0.734772233787208`, `regtype.engine = "lp"`, `degree = 1`,
  `basis.engine = "glp"`, `bernstein.basis.engine = TRUE`.

Working diagnosis:

The failure is an MPI collective mismatch in the conditional-density bandwidth
path. It is not explained by missing high-level option forwarding, and it is not
multistart noise.

## Scientific And Engineering Contract

1. `npRmpi` must remain functionally equivalent to `np` for the selected
   estimator semantics.
2. A route that is expected to fail because the raw local-linear basis is
   numerically unstable must fail promptly with the same user-facing guidance,
   not hang in MPI collectives.
3. MPI workers must participate in the same computational collective contract as
   the coordinator, or the coordinator must run the route in explicitly local
   single-process mode. Mixed coordinator-C-collective plus worker-Rmpi-command
   collectives are invalid.
4. No repair may silently change user-selected method semantics, estimator
   semantics, bandwidth defaults, quadrature controls, or search controls.
5. No repair may be accepted on the patents repro alone. The campaign must audit
   adjacent density/distribution bandwidth selectors and any shared MPI
   collective helper path that could share the defect.

## Initial Repair Plan

### Phase 0: Preserve The Failure Evidence

Retain and index the existing proof files:

- `results/np_patents_ll_default_timing_trace.log`
- `results/np_patents_ll_nmulti1_timing_trace.log`
- `results/nprmpi_patents_ll_default_trace_monitor.log`
- `results/nprmpi_patents_ll_autodFALSE_30s.log`
- `results/nprmpi_patents_ll_autodFALSE_nmulti1_30s.log`
- `results/nprmpi_patents_ll_forced_local_30s.log`
- `results/nprmpi_patents_ll_stack_sample.log`
- `results/nprmpi_patents_ll_autodTRUE_stack_sample.log`

Create a short trace summary in the same tmp bundle before patching. The summary
must state exact elapsed times, process IDs sampled, stack locations, and the
conclusion that the observed issue is a collective mismatch.

### Phase 1: Route Inventory And Risk Map

Audit all bandwidth-selection routes that can enter compiled MPI collectives
while invoked from a master/slave Rmpi execution context:

- `npcdensbw()` and `npcdens()`
- `npcdistbw()` and `npcdist()`
- `npudensbw()` and `npudens()`
- `npudistbw()` and `npudist()`
- Regression-family routes using the analogous SPMD/autodispatch machinery:
  `npregbw()`, `npscoefbw()`, `npplregbw()`, and `npindexbw()`

For each route, classify:

- whether the top-level R route uses autodispatch, manual broadcast, local
  shadow, or direct MPI C collectives;
- whether workers execute the same compiled routine as the coordinator;
- whether failures in the compiled routine are collectively propagated;
- whether the route contains an explicit local-mode escape comparable to
  `.npRmpi_with_local_regression(...)`;
- whether the route already has installed sentinels for `nslaves = 1` and
  `nslaves > 1`.

Do not patch until the inventory identifies the minimal shared defect surface.

### Phase 2: Determine The Intended MPI Contract For Conditional Density

For `npcdensbw()` fixed-bandwidth `cv.ls` with local-linear/local-polynomial
continuous regressors, decide from code evidence which of the following is the
correct contract:

1. True compiled MPI path:
   The coordinator and all workers enter the same `C_np_density_conditional_bw`
   calculation and participate in the same `MPI_Allreduce` / `MPI_Allgather`
   sequence.
2. Coordinator-local path:
   The route is intentionally unsafe or not beneficial under Rmpi worker
   dispatch, so it must be evaluated under `.npRmpi_with_local_regression(...)`
   while preserving the public `npRmpi` interface and error behavior.
3. Dedicated shadow path:
   A route-specific shadow evaluator must be used so the worker payload is
   self-contained and collective-safe.

The decision must be evidence-based. It must cite the exact R and C call sites
that prove the route can or cannot safely execute as a true compiled MPI path.

### Phase 3: Candidate Repair

Implement only the smallest repair that enforces the chosen contract.

Likely candidate, pending Phase 2 evidence:

- For conditional-density fixed-bandwidth `cv.ls` local-linear/raw-basis search
  routes that currently create a coordinator-C-collective / worker-Rmpi-command
  mismatch, route the bandwidth computation through explicit local mode before
  entering `C_np_density_conditional_bw`.
- If `npcdistbw()` has the same collective mismatch shape, apply the same
  contract there.
- If `npudensbw()` or `npudistbw()` expose the same stack mismatch, include
  them in the same tranche only if the repair is the same shared helper
  contract; otherwise split into separate tranches.
- Do not change the serial `np` implementation.
- Do not alter optimizer semantics, floors, bounds, quadrature defaults,
  NOMAD/Powell policy, or user argument validation.

### Phase 4: Failure-Stage Instrumentation

Before and after patching, capture:

- selected public call and resolved method/options;
- `regtype`, engine regtype, degree, basis, Bernstein flag;
- `bwmethod`, `bwtype`, `nmulti`, `scale.factor.*`, quadrature controls;
- coordinator and worker rank/size;
- whether the call is inside autodispatch broadcast;
- whether local mode is active;
- exact compiled call entered, if any;
- exact error or success status;
- elapsed time and process CPU pattern.

Instrumentation must be in scratch harnesses or tests, not permanent verbose
runtime output, unless the permanent output is a deliberate user-facing
condition.

### Phase 5: Validation Gates

Minimum installed-build gates for `npRmpi`:

1. Patents failure sentinel:
   - `npcdens(..., regtype = "ll", bwmethod = "cv.ls")`
   - default multistart and `nmulti = 1`
   - `autodispatch = TRUE` and `autodispatch = FALSE`
   - `nslaves = 1`
   - must fail within `2x` serial `np` baseline plus a small absolute MPI
     overhead allowance, and with the expected guidance.
2. Positive patents alternative:
   - `regtype = "lp", degree = 1, bernstein.basis = TRUE`
   - must remain functional.
3. Adjacent conditional-density smoke:
   - `regtype = "lc"` still succeeds/fails as before with no objective drift.
   - bounded quadrature routes remain functional.
4. Conditional distribution parity:
   - run equivalent `npcdistbw()` / `npcdist()` `ll` and success-route sentinels
     if inventory shows shared risk.
5. Unconditional density/distribution guard:
   - run focused `npudensbw()` and `npudistbw()` `cv.ls` MPI sentinels to verify
     no shared-helper regression.
6. MPI mode coverage:
   - at least `nslaves = 1` and `nslaves = 2`.
   - autodispatch on and off where both modes are supported.
7. Package gates:
   - focused testthat sentinel set.
   - installed private build.
   - `R CMD check --no-manual --ignore-vignettes`.
   - `R CMD check --as-cran` only after the focused campaign is green.

All logs must be retained under a dated path in
`/Users/jracine/Development/tmp`.

## Critique Of The Initial Plan

The initial plan is directionally sound but too permissive in four ways.

1. It risks broadening too early.
   Auditing all density, distribution, and regression routes is appropriate,
   but patching across all touched routes in one tranche would create an
   unnecessary blast radius. The observed deterministic defect is currently
   `npcdensbw()` conditional-density `cv.ls`; other routes should be included
   only if they exhibit the same collective mismatch or share the exact helper
   repaired.

2. The candidate repair is still too vague.
   "Route through local mode" might be correct, but it must not be chosen merely
   because it makes the patents repro finish. If a true MPI implementation is
   intended and only one worker-dispatch gate is missing, a local-mode bypass
   would throw away intended parallelism. The plan must require a direct
   contract decision from source evidence before selecting the repair shape.

3. The validation timing rule needs a precise threshold.
   "`2x` serial baseline plus small overhead" is useful as a diagnostic rule,
   but acceptance should use a concrete wall-clock cap recorded in the harness,
   e.g. `max(15 seconds, 2 * serial_elapsed + 5 seconds)` for the failure
   sentinel. This avoids ambiguous pass/fail interpretation.

4. The plan needs a collective-safety invariant test.
   A one-time stack sample proves the bug, but future regressions need an
   automated guard. The final repair should include a test or harness that fails
   if a known-fast error route enters the MPI C collective and stalls instead of
   returning a condition. Where feasible, use a short timeout sentinel and
   explicit process cleanup.

## Reformulated Execution Plan

### Tranche 1: Locked Diagnosis, No Production Patch

1. Create a trace-summary artifact under:

   ```text
   /Users/jracine/Development/tmp/nprmpi_cdens_collective_mismatch_20260427
   ```

2. Re-run only short, bounded probes:
   - serial `np` baseline;
   - `npRmpi` forced-local;
   - `npRmpi` direct MPI with a `30s` timeout;
   - one stack sample confirming `MPI_Allreduce` vs worker `MPI_Allgather`.

3. Record a route map for:
   - `npcdensbw()`
   - `npcdistbw()`
   - `npudensbw()`
   - `npudistbw()`
   - the NOMAD shortcut route for `npcdensbw(..., nomad = TRUE)`

4. Stop if the new trace contradicts the existing diagnosis.

Exit criterion:

- The route map states exactly whether the defect is isolated to
  `npcdensbw()` or shared with another density/distribution selector.
- The route map separately states whether the working NOMAD shortcut path is
  independent of the failing raw local-linear path, and whether it already uses
  a collective-safe worker contract.

Decision discipline:

- Route-map entries are evidence, not repair authorization. A route may be
  patched only when the map identifies either the exact same collective mismatch
  or a shared helper that must be changed to repair the primary defect.
- If a sibling route is merely suspicious but not reproduced, add a sentinel and
  leave production behavior untouched.

### Tranche 2: Minimal Contract Repair

Patch only the smallest R/C surface needed to enforce one collective contract.

Preferred ordering:

1. If source evidence shows workers are supposed to enter the same compiled
   routine, repair the worker execution path so both ranks enter the same
   collective sequence.
2. If source evidence shows `npcdensbw()` cannot safely enter true MPI
   collectives from the current Rmpi worker-command architecture for this route,
   enforce local mode at the narrow R boundary for the affected fixed `cv.ls`
   raw local-linear/local-polynomial failure path only. This local-mode use must
   be explicit, documented in the code comment or test name, and justified by
   route-map evidence.
3. If a dedicated shadow route already exists and can be reused without changing
   estimator semantics, prefer that over broad local-mode routing when it
   preserves intended MPI execution and has a smaller behavioral blast radius.
4. If `npcdistbw()` has the same proven defect, apply the same narrow repair
   there in the same tranche only if the helper and validation surface are
   materially identical.
5. Do not touch `npudensbw()` / `npudistbw()` / regression routes unless the
   route map proves they share the same faulty helper or exhibit the same
   collective mismatch.

Hard prohibitions:

- no silent serial fallback for broad `npRmpi` execution;
- no option remapping;
- no optimizer fallback;
- no NOMAD-only degradation;
- no changes to `np-master`;
- no changes to manuscript, simulation, or campaign files.
- no permanent diagnostic chatter;
- no changes to the currently working `nomad = TRUE` semantics unless the
  route-map evidence proves the NOMAD route shares the repaired helper.

NOMAD preservation rule:

- Before patching, establish the serial `np` and `npRmpi` installed-build timing
  for `npcdens(..., nomad = TRUE)` on the patents data.
- After patching, repeat the same timing. The repair is not acceptable if the
  NOMAD route loses success status, result shape, or its expected speedup
  direction relative to serial `np`, unless a separate root cause is identified
  and repaired in the same campaign.

### Tranche 3: Focused Regression Tests

Add focused tests that exercise the repaired contract using small data and the
patents sentinel when feasible.

Required test classes:

- fast unit-level guard for the chosen dispatch/local-mode condition;
- installed-script sentinel for patents `npcdens ll cv.ls`;
- success-route sentinel for Bernstein local polynomial;
- NOMAD shortcut sentinel for patents `npcdens(..., nomad = TRUE)`;
- adjacent `lc` conditional-density sentinel;
- route-map-driven sentinels for any additionally affected selector.

Timeout behavior:

- For the patents failure sentinel, accepted elapsed time must be no worse than
  `max(15, 2 * serial_np_elapsed + 5)` seconds on the local machine.
- Any timeout must clean worker processes and fail the gate.

### Tranche 4: Installed Validation

Build and install private `npRmpi`, then run:

- focused testthat file(s);
- patents failure/success sentinels with `nslaves = 1`;
- patents failure/success sentinels with `nslaves = 2`;
- patents NOMAD shortcut sentinel with serial `np`, `npRmpi nslaves = 1`, and
  `npRmpi nslaves = 2`;
- autodispatch on/off variants where supported;
- conditional distribution and unconditional density/distribution sentinels only
  as required by the route map;
- `R CMD check --no-manual --ignore-vignettes`.

Only after the above is green, run final `R CMD check --as-cran`.

### Tranche 5: Checkpoint Commit And Closeout

Commit only after installed gates pass.

Closeout must report:

- exact files changed;
- exact repair contract chosen and why;
- serial `np` baseline time;
- `npRmpi` post-repair times;
- pre/post `nomad = TRUE` times and result summaries;
- any routes audited and explicitly excluded from patching;
- any remaining risks or route-specific limitations.

## Acceptance Definition

The campaign is complete only when:

1. `npRmpi` no longer hangs or exceeds the timing cap on the patents
   local-linear raw-basis failure route.
2. The user receives the intended guidance condition.
3. Positive conditional-density routes still work.
4. The patents `nomad = TRUE` route remains functional and preserves its
   expected speedup direction relative to serial `np`.
5. Any affected sibling route identified by the route map is repaired and
   validated.
6. No unrelated estimator family is touched without evidence.
7. The final private installed build passes the named gates.
