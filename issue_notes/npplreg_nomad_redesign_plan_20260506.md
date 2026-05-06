# npplreg NOMAD Redesign Plan

Date: 2026-05-06

## Goal

Repair and modernize `npplreg(..., nomad = TRUE)` in both `npRmpi` and
`np-master` by replacing the inherited joint high-dimensional NOMAD abstraction
with the natural partially linear decomposition:

1. select bandwidths and local-polynomial degree for `y ~ z`;
2. select bandwidths and local-polynomial degree for each `x_j ~ z`;
3. perform the final partially linear solve from the resulting residualized
   components.

The key statistical contract is:

> `npplreg` is a composition of nuisance regressions plus a final partially
> linear solve. NOMAD belongs to each nuisance regression, not to one artificial
> product-space optimizer over all child regressions at once.

This is not merely a performance refactor. It is the natural local-polynomial
generalization of the partially linear estimator:

- `regtype = "lc"` means every child regression uses degree `0`;
- `regtype = "ll"` means every child regression uses degree `1`;
- fixed `regtype = "lp"` accepts a common degree vector as a compatibility
  shorthand, but the full API must also accept a child-indexed degree matrix or
  list in the same order as the bandwidth object:
  `yzbw`, then each `x_j` child;
- `nomad = TRUE` lets each child regression select its own admissible degree and
  bandwidths.

The final partially linear fit is not a NOMAD problem.

## Current Diagnosis

The old `npplreg` NOMAD route optimizes one combined bandwidth/degree point for
the full collection of child regressions. Conditional on the child regression
specification, the objective is a separable sum of child regression criteria.
That makes the product-space optimizer the wrong abstraction:

- it searches a larger space than the estimator requires;
- it imposes a common selected degree where different nuisance regressions
  should be allowed to choose differently;
- it encouraged route-specific MPI scheduling machinery that split children
  across ranks or ran service loops, rather than delegating work to the proven
  `npregbw()` MPI-aware regression kernels;
- it made profile, attach, and session behavior harder to reason about and
  easier to accidentally desynchronize.

Earlier scratch experiments that tried to make the old abstraction faster are
now diagnostic-only. They showed where time and CPU were being lost, but the
right repair is to remove the joint abstraction rather than keep optimizing it.

## Replacement Design

Create a private child-specific orchestration path for `npplregbw(...,
nomad = TRUE)`:

1. Build the standard child response list:
   - child `1`: `y ~ z`;
   - children `2:(p + 1)`: each `x_j ~ z`.
2. For each child, call the already validated `npregbw()` automatic
   local-polynomial degree-search route with the active MPI pool:
   - `regtype = "lp"`;
   - `nomad = TRUE`;
   - `degree.select`, `search.engine`, `degree.min`, `degree.max`,
     `degree.start`, `degree.restarts`, `degree.max.cycles`,
     `nomad.nmulti`, `bernstein.basis`, kernels, and bounds propagated
     faithfully;
   - no explicit `degree`, because `npregbw(..., nomad = TRUE)` owns degree
     selection.
3. Sum the child objective values for the `plbandwidth` public `fval`.
4. Sum `num.feval` and `num.feval.fast` across children.
5. Preserve per-child `rbandwidth` objects as the source of truth for:
   - bandwidths;
   - selected degree;
   - child objective;
   - child timing/evaluation metadata where available.
6. Perform the downstream `npplreg(bws = ...)` fit exactly as a sequential
   composition of these child objects.

The old joint `npplreg` NOMAD functions should first be bypassed, not deleted.
Once the new path passes installed sentinels and release gates, remove the dead
scaffold in a separate cleanup tranche.

## Object And Reporting Contract

`plbandwidth` must report child-level degree information beside the child
bandwidths.

Required summary shape:

```text
y on z:
Degree: ...
Bandwidth: ...

x1 on z:
Degree: ...
Bandwidth: ...

x2 on z:
Degree: ...
Bandwidth: ...
```

For `lc` and `ll`, this should still be explicit:

- `lc`: all child degree blocks print `0`;
- `ll`: all child degree blocks print `1`.

For `nomad = TRUE`, child degrees may differ. The child `rbandwidth` objects are
the source of truth. A legacy/global `bws$degree` field may remain as a
compatibility convenience only; summary, fit, prediction, plotting, and future
validation should not rely on it when child degrees differ.

## Progress And User-I/O Contract

The child-specific design must be visible to the user while it runs. Long
silent stretches are not acceptable.

For `npplregbw(..., nomad = TRUE)`, progress should identify:

1. the total number of nuisance regressions, `p + 1`;
2. the active child regression:
   - `E[y | z]`;
   - `E[x1 | z]`;
   - ...;
   - `E[xp | z]`;
3. the completion index, for example `2/4`;
4. whether the child is in the NOMAD search or Powell refinement;
5. the final partially linear fit stage after all child bandwidth searches
   complete.

The messages should follow the same philosophy as the plot progress messages:
they should reassure the user that real stages are advancing, not merely print
abstract optimizer activity. The child label must be present in both NOMAD and
Powell progress where the underlying `npregbw()` route emits progress.

## Implementation Tranches

### Tranche 0: Scratch-Only Candidate

Start in `npRmpi`, because the active contradiction arose from MPI route
behavior. Use the detached worktree:

`/Users/jracine/Development/tmp/npplreg_nomad_separable_20260506/np-npRmpi-separable`

No live-repo package source edits until this tranche passes installed proof.

Tasks:

1. Add a private child-specific NOMAD helper for `npplregbw.default()`.
2. Route `npplregbw(..., nomad = TRUE)` through that helper.
3. Add `plbandwidth` child-degree metadata derived from the child `rbandwidth`
   objects.
4. Update `print.plbandwidth()` and `summary.plbandwidth()` to display child
   degrees.
5. Add child-aware progress labels for `E[y | z]`, `E[x_j | z]`, NOMAD,
   Powell, and final fit stages.
6. Leave the old joint NOMAD scaffold in place but unused.

### Tranche 1: Correctness Proof

Use fixed seeds and small enough data for quick iteration, then representative
sizes for meaningful confidence.

Required checks:

1. `lc`, `ll`, and fixed `lp` degree routes preserve existing numerical output
   and now report explicit child degrees.
2. `nomad = TRUE` returns valid child-specific degrees and bandwidths.
3. `npplreg(bws = bw)` consumes the returned object without needing a global
   common degree.
4. Fitted values, coefficients, residual lengths, and prediction/evaluation
   paths are valid.
5. Child objective values sum to public `bw$fval`.
6. Progress output has no silent child-search or Powell-refinement gap longer
   than the configured progress delay/throttle window.

Numerical expectations:

- fixed-degree routes should remain numerically identical;
- NOMAD child-specific degree search may legitimately change selected degrees,
  bandwidths, and final objective relative to the old forced-common-degree
  route;
- any changed NOMAD result must be explained by the new estimator contract, not
  by accidental option drift.

### Tranche 2: MPI Route Proof

Validate the same installed scratch build across:

- serial;
- profile;
- attach;
- session.

Milestones:

1. `nslaves = 1` should be practically no worse than `np-master`/serial for a
   sufficiently large workload and should use the proven `npregbw` route.
2. For sufficiently large `n`, `nslaves > 1` should improve elapsed time.
3. Active child-regression phases should reach high CPU utilization when `n` is
   large enough that MPI communication overhead is not the explanation.
4. Profile, attach, and session should show the same qualitative scaling shape;
   large route gaps become defects to diagnose, not acceptable mode semantics.

### Tranche 3: Documentation And Live Migration

Only after Tranches 0-2 pass:

1. Update `.Rd` language to describe child-specific NOMAD degree selection for
   `npplreg`.
2. Port surgically into the live `np-npRmpi` repo.
3. Run installed sentinels for `npplreg` and adjacent semiparametric routes.
4. Run tarball-first package checks.
5. Commit only after green installed proof.

### Tranche 4: Cross-Package Parity

After the `npRmpi` design is green, apply the same estimator contract to
`np-master`. `np-master` is not a valid performance or semantic baseline for
the new NOMAD route until this tranche is complete.

Tasks:

1. Refactor `np-master::npplregbw(..., nomad = TRUE)` to use the same
   child-specific `npregbw()` orchestration.
2. Add the same child-degree metadata/reporting contract to `plbandwidth`.
3. Add the same fixed-`lp` child-degree API in `np-master`:
   common vector shorthand, plus canonical child-indexed matrix/list input.
4. Update `np-master` `.Rd` documentation and examples.
5. Validate `np-master` fixed-degree parity, child-specific fixed-degree API,
   and child-specific NOMAD behavior.
6. Validate `npRmpi` and `np-master` against comparable fixed-seed examples
   after both packages implement the same statistical contract.

## Critical Self-Review

The plan above is attractive, but the risks are real:

1. It changes estimator semantics for `nomad = TRUE` by allowing child-specific
   selected degrees. That is intended, but it must be documented as a deliberate
   statistical correction, not smuggled in as a speed fix.
2. Downstream code may still assume `plbandwidth$degree` is global. The repair
   must audit fit, predict, plot, demo parsing, and summary paths and make child
   objects authoritative where needed.
3. Calling public `npregbw()` inside `npplregbw()` is the right abstraction for
   this estimator, but it must not create nested MPI dispatch or lifecycle
   trouble in attach/session/profile.
4. Child-specific degree search may increase total work for tiny examples,
   because it now runs several full child searches. Acceptance should be based
   on correctness and large-enough scaling, not tiny-`n` noise.
5. We should not delete the old joint scaffold until the new path has passed
   installed sentinels. Removing too much at once would combine semantic,
   performance, and cleanup risk.
6. The first scratch implementation should not attempt a clever custom
   scheduler. If simple child-specific `npregbw()` orchestration is not fast
   enough, diagnose why the child calls differ from standalone `npregbw()`
   before adding any new MPI machinery.
7. The plan must treat progress output as a correctness surface. A numerically
   valid repair that leaves users staring at a static `Refining...` line for a
   long child Powell phase is incomplete.
8. Cross-package parity is required. `npRmpi` should prove the route first
   because the contradiction arose from MPI behavior, but `np-master` must be
   refactored before final closeout so both packages expose the same estimator
   semantics.

## Reformulated Plan

Proceed in one narrow production-shaped tranche:

1. Implement child-specific NOMAD orchestration in scratch only.
2. Store/report child degrees from child `rbandwidth` objects.
3. Preserve fixed-degree common-vector semantics exactly while adding the
   proper child-indexed fixed-`lp` API:
   - `degree = c(...)` means the same degree vector for every child;
   - `degree = matrix(...)` must have `(p + 1)` rows by `ncon(z)` columns in
     `yzbw`, `x1`, ..., `xp` order, with row names accepted and reordered;
   - `degree = list(...)` must have one entry per child, with names accepted and
     reordered;
   - missing or ambiguous child names must fail fast.
4. Validate `nomad = TRUE` as a deliberate new estimator contract with
   child-specific selected degrees.
5. Instrument child-aware progress before timing claims:
   - `E[y | z]`, `E[x_j | z]`, and `i/(p + 1)` must appear;
   - NOMAD and Powell phases must identify the active child;
   - the final partially linear fit must be announced.
6. Audit consumers of `plbandwidth$degree` before live migration:
   - fit/eval/predict must use child `rbandwidth` objects;
   - summary/print must display child degrees;
   - plot/demo/test helpers must not silently collapse child-specific degrees
     into one misleading scalar.
7. Validate profile/attach/session after the installed scratch build proves the
   local object contract.
8. Only then migrate live and document the new contract.
9. Port the settled contract to `np-master`, validate, document, and commit so
   the two packages remain functionally aligned.

## Comfort Level After Refinement

I am comfortable with this plan because it minimizes regression risk in three
ways:

1. it uses `npregbw()` as the proven child kernel instead of inventing a new
   optimizer or MPI scheduler;
2. it keeps fixed-degree `lc`/`ll`/manual `lp` paths on their existing semantics;
3. it makes the intended `nomad=TRUE` semantic change visible in object
   metadata, summary output, progress output, documentation, and validation.

The remaining discomfort is appropriate: this changes `npplreg`'s automatic
degree-selection contract. That is the correct statistical move, but it means
the acceptance proof must focus on estimator validity, child-object integrity,
and route behavior rather than forcing equality to the old joint optimizer.
