# Bootstrap Interface Rewrite Status (npRmpi)

## Scope

This note summarizes what has been learned while migrating `npRmpi` from explicit user-side `mpi.bcast.cmd(...)` wrapping toward a simpler autodispatch interface, with emphasis on bootstrap-heavy paths (`npsigtest`, `plot(..., plot.errors.method="bootstrap")`).

## Current Stable State

1. Core autodispatch workflow works for major `np*` calls.
2. `npsigtest` direct autodispatch now works in smoke tests (`boot.num=19`).
3. Non-bootstrap plotting works in direct autodispatch.
4. `npregression` bootstrap plotting now works in direct autodispatch via payload bootstrap engine.
5. Non-`rbandwidth` bootstrap plotting remains stable via explicit manual broadcast:
   - `mpi.bcast.cmd(plot(..., plot.errors.method="bootstrap"), caller.execute=TRUE)`

## Key Findings (Confirmed)

1. MPI execution model is the dominant constraint:
   - many `np*` fit/predict operations are effectively collective in MPI context.
   - ranks must execute compatible call sequences; divergence causes hangs.

2. Rank-splitting bootstrap replications is unsafe for current `npRmpi` internals:
   - assigning replicate subsets to different ranks and having each rank run independent `npreg` loops can hang.

3. Synchronized-rank bootstrap loops are viable in minimal settings:
   - all ranks run the same replicate loop (`b = 1..B`) with identical ordering.

4. Mixed auto/manual object lifecycle is fragile:
   - objects created under autodispatch can carry state/shape that does not always safely survive transition into manual broadcast bootstrap workflows.

5. Disabling graphics devices on slaves is not the root cause:
   - hangs were reproducible in compute-only bootstrap prototypes with no plotting.

## Problems Still Open

1. Direct autodispatch bootstrap plotting remains unresolved outside the `rbandwidth` path.
2. There is still at least one unresolved bridge issue between:
   - autodispatch-created bandwidth/object state, and
   - manual-style distributed bootstrap compute for non-`rbandwidth` paths.

## Failed / Insufficient Approaches So Far

1. Generic direct-autodispatch enablement for bootstrap plot paths:
   - produced hangs; reverted to guarded behavior.

2. Rank-partitioned bootstrap execution:
   - deadlocks due to collective mismatch in internal MPI calls.

3. Object untagging alone as a cure:
   - necessary for some guard interactions, but not sufficient to eliminate all hangs.

4. Replacing temp-object broadcast strategy (`assign`-literal vs symbol broadcast):
   - did not by itself resolve bootstrap-plot deadlocks.

5. Forcing direct plot routing through experimental wrapper paths:
   - unstable in current form; reverted from unstable variants.

## Changes That Did Work

1. `npsigtest` autodispatch stabilization:
   - route direct call through a manual-style distributed helper.
   - for `npregression` inputs, explicitly recover `xdat`/`ydat` from model frame on master.
   - normalize `bws` in `npsigtest.rbandwidth`.
   - avoid residuals re-entry path that dropped `txdat`.

2. Plot policy refinement:
   - non-bootstrap direct autodispatch supported.
   - bootstrap plot explicitly guarded with clear remediation message for non-`rbandwidth` paths.
3. `rbandwidth` bootstrap payload engine:
   - synchronized payload bootstrap loop now drives direct autodispatch bootstrap for `plot.npregression`,
   - smoke-tested with `plot.errors.boot.num=19` under `npRmpi.autodispatch=TRUE`.

## Best Current Strategy

1. Continue staged migration by method family (`rbandwidth` done first).
2. Keep stable guard for all remaining direct bootstrap plotting until each method gets a payload path.
3. Build a dedicated synchronized bootstrap compute helper (single execution model):
   - one MPI layer
   - same replicate loop on all ranks
   - deterministic replicate seeding
   - master-only collection of outputs and rendering.
4. Integrate this helper into one method at a time, then re-open direct plot autodispatch for that method only.
5. Validate each increment with small deterministic smoke tests (`B=2`, small `n`) before scaling.

## Immediate Practical Guidance

1. For now, use:
   - direct autodispatch for non-bootstrap plotting and estimators,
   - direct autodispatch bootstrap plotting for `npregression`,
   - explicit `mpi.bcast.cmd(...)` for bootstrap plotting in other method families.
2. For significance testing:
   - direct `npsigtest` autodispatch path is now available (smoke-tested).
