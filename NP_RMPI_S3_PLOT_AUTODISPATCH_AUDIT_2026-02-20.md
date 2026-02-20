# npRmpi Audit: S3 Plot + Autodispatch (2026-02-20)

Scope: recent changes in `np-npRmpi` focused on (i) S3 plotting path, (ii) autodispatch replacing explicit `mpi.bcast.*` wrappers, and (iii) production-readiness gaps.

## Executive verdict

Current state is a strong functional improvement for users:

- direct autodispatch works for major calls (no mandatory `mpi.bcast.cmd(...)` wrappers in normal use),
- direct bootstrap plotting now executes in smoke tests,
- `npsigtest()` direct autodispatch also executes.

But the plotting implementation is **not yet modern/modular end-state**. The package still relies on the legacy monolithic `npplot` backend and extensive `eval(parse(...))` chains. So this is **functionally improved but structurally transitional**, not yet the final production architecture.

## What was validated now

Smoke run (`npRmpi.start(nslaves=1)`, `options(npRmpi.autodispatch=TRUE, np.messages=FALSE)`) with direct calls and no `mpi.bcast.*` wrappers:

- `plot(npreg(...), plot.errors.method="bootstrap", plot.errors.type="all", plot.errors.boot.num=9)` -> pass
- `plot(npcdens(...), ... bootstrap ...)` -> pass
- `plot(npcdist(...), ... bootstrap ...)` -> pass
- `plot(npudens(...), ... bootstrap ...)` -> pass
- `plot(npudist(...), ... bootstrap ...)` -> pass
- `npsigtest(npreg_obj, boot.num=19)` -> pass

Observed in `/tmp/audit_plot_autodispatch_smoke2.log`.

## Architecture reality (important)

### S3 layer exists, but it is thin wrappers

- S3 `plot.*` methods are defined in `R/np.plot.methods.R`.
- Most methods immediately delegate to `npplot(...)` by extracting `bws` (or `bw`) and forwarding dots.

Key references:

- `/Users/jracine/Development/np-npRmpi/R/np.plot.methods.R:1`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.methods.R:53`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.methods.R:71`

This means modular S3 dispatch is present at the API boundary, but core plotting logic still lives in monolithic `npplot` methods.

### Monolithic backend remains dominant

- `npplot` generic and class methods are still central and huge (`~7.2k` LOC in `R/np.plot.R`).

Key references:

- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:1452`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:1471`

## Findings (ranked)

### P1: Structural debt remains high in plotting core

The plot engine still contains very large monolithic logic with many repeated `eval(parse(...))` execution paths, which are fragile to refactor and hard to validate exhaustively.

Representative locations:

- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:394`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:2013`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:2660`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:3491`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:4893`

Impact:

- high maintenance risk,
- non-local side effects harder to reason about,
- higher chance of regression in edge plotting paths.

### P1: `npplot` is still required; cannot be safely removed yet

S3 methods are wrappers around `npplot`. Removing `npplot` today would break most plotting behavior.

References:

- `/Users/jracine/Development/np-npRmpi/R/np.plot.methods.R:6`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.methods.R:53`
- `/Users/jracine/Development/np-npRmpi/NAMESPACE:151`
- `/Users/jracine/Development/np-npRmpi/NAMESPACE:390`

### P2: Global environment staging in autodispatch is deliberate but still check-NOTE prone

Autodispatch materializes temp args in `.GlobalEnv` and cleans them afterward; this triggered known `R CMD check` NOTE about assignment to global environment.

References:

- `/Users/jracine/Development/np-npRmpi/R/np.autodispatch.R:216`
- `/Users/jracine/Development/np-npRmpi/R/np.autodispatch.R:230`

Impact:

- likely acceptable operationally,
- but undesirable for strict package-hygiene and future review burden.

### P2: Formula fail-fast currently requires explicit data-like args

Autodispatch intentionally errors if formula interface is used without `data=`/`xdat`/`ydat` style arguments to avoid unresolved symbols on slaves.

Reference:

- `/Users/jracine/Development/np-npRmpi/R/np.autodispatch.R:159`

Impact:

- safer runtime behavior,
- but still a usability constraint versus ideal “it just works from parent env symbols.”

### P3: Legacy guard text/path still implies old bootstrap limitations

`np.plot.R` still contains guard scaffolding originally intended to block direct autodispatch bootstrap in some contexts; current behavior now allows direct bootstrap in key paths via payload.

References:

- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:141`
- `/Users/jracine/Development/np-npRmpi/R/np.plot.R:1590`

Impact:

- mostly documentation/consistency risk now,
- low immediate functional risk if tests stay green.

## What can be safely removed now

Safe immediate removals (low risk):

1. Any remaining old-style `mpi.bcast.cmd(np*...)` example/test code in `np*` man/demo/test assets that has direct autodispatch equivalent and has already been validated.
2. Redundant warning/guard text in docs that says direct bootstrap plotting is unsupported (if still present after doc refresh).
3. Obsolete transitional comments referring to mandatory wrapper workflows.

## What should NOT be removed yet

Do not remove yet:

1. `npplot` generic/method family.
2. `compute.bootstrap.errors.*` S3 methods.
3. `npplot` S3 registrations in `NAMESPACE`.

Rationale: `plot.*` methods currently delegate into this backend; removing these now would break broad behavior.

## Required patch queue for production readiness

### Phase A (stabilization, no behavior change)

1. Replace high-risk `eval(parse(...))` blocks in `R/np.plot.R` with structured `do.call()`/function dispatch for top-level duplicated patterns.
2. Centralize bootstrap payload invocation so each family (`npreg`, `npcdens`, `npcdist`, etc.) does not carry duplicated expression-construction logic.
3. Add targeted unit tests for direct autodispatch bootstrap plotting per class (`npregression`, `condensity`, `condistribution`, `npdensity`, `npdistribution`) with `plot.errors.type='all'` and `plot.errors.boot.num=9`.

### Phase B (modularization)

1. Move class-specific plotting logic out of monolithic `np.plot.R` into dedicated S3 helpers.
2. Keep `npplot` as compatibility shim initially, then progressively route to modular internals.
3. Only after parity and tests, deprecate/remove external reliance on `npplot` in docs/examples.

### Phase C (hygiene)

1. Remove `.GlobalEnv` assignment pattern from autodispatch temp staging if possible (or isolate and annotate if not avoidable due Rmpi constraints).
2. Update `R CMD check` expectation and document rationale if NOTE remains intentional.

## Production-readiness assessment

- Functional readiness for core direct-autodispatch workflows: **Good**.
- Structural/modularity readiness of plotting subsystem: **Not yet complete**.
- Risk profile: **Moderate** (mainly maintainability and regression risk from legacy monolith + dynamic eval).

## Recommended next checkpoint

Implement Phase A in small patches with parity tests after each patch, then re-run:

- targeted plot bootstrap smoke matrix,
- `tests/testthat` suite,
- extracted `np*.Rd` dontrun script.

