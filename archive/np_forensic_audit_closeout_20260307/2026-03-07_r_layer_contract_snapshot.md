# R-Layer Contract Snapshot for Future Consolidation

Date: 2026-03-07
Repo: `np-master`
Status: freeze note only, no code change

## Purpose

Record the concrete duplication seams that must be frozen before any R-layer
consolidation tranche. This note is intentionally narrower than the forensic
audit language: it identifies where the current contracts live, not how to
rewrite them yet.

## Primary Freeze Surface

The highest-value shared shape is the bandwidth-selector family:

1. `npregbw` in [`R/np.regression.bw.R`](/Users/jracine/Development/np-master/R/np.regression.bw.R)
2. `npudensbw` in [`R/np.density.bw.R`](/Users/jracine/Development/np-master/R/np.density.bw.R)
3. `npcdensbw` in [`R/np.condensity.bw.R`](/Users/jracine/Development/np-master/R/np.condensity.bw.R)
4. `npudistbw` in [`R/np.distribution.bw.R`](/Users/jracine/Development/np-master/R/np.distribution.bw.R)
5. `npcdistbw` in [`R/np.condistribution.bw.R`](/Users/jracine/Development/np-master/R/np.condistribution.bw.R)

These files all repeat the same broad pattern:

1. generic dispatch via `.np_bw_dispatch_target(...)`
2. `.formula` entry that builds a `model.frame`
3. `.NULL` entry that normalizes raw inputs and recursively re-enters `.default`
4. class-specific bandwidth entry (`*.bandwidth` / `*.rbandwidth`)
5. `toFrame()` then `toMatrix()` coercion passes
6. call/formula/rows-omitted metadata repair after recursion

## Contracts That Must Not Drift

Before any helper extraction, preserve:

1. formula parsing behavior, including the time-series `predvars` rewrite paths
2. `rows.omit`, `nobs.omit`, `terms`, `formula`, and `call` fields
3. current defaulting of `nmulti`, `remin`, `bandwidth.compute`, and
   `invalid.penalty`
4. type checks and type-specific error messages for factor/ordered/numeric
   inputs
5. current `np.tree` option usage at call-materialization time
6. the exact `.Call(...)` payload shaping for the compiled entry points

## Safe First Extraction Candidates

Only two helper shapes currently look safe enough for a future first refactor:

1. a pure `model.frame`/terms helper for `.formula` selectors
2. a pure metadata finalizer for `call`/`formula`/`rows.omit` repair after the
   recursive handoff

Anything that mixes coercion, validation, and `.Call(...)` argument assembly in
one helper is still too risky for the first consolidation tranche.

## `npRmpi` Relevance

The mirrored selector files in `np-npRmpi` retain the same core duplication but
overlay MPI-specific route guards and autodispatch checks. Example:

1. [`np-npRmpi/R/np.density.bw.R`](/Users/jracine/Development/np-npRmpi/R/np.density.bw.R)

That means:

1. `np-master` remains the source-of-truth for shared contract freezing,
2. helper extraction must be proven serial-safe before any port,
3. the `npRmpi` port must preserve route-guard placement and must not collapse
   MPI-specific wrapper behavior into serial helpers.

## Recommended Next Step

If this tranche is pursued later, start with a proof-only helper extraction on
one selector family (`npregbw` preferred), backed by snapshot tests that freeze:

1. formula-to-default object fields,
2. raw-input-to-default object fields,
3. representative error messages,
4. `.Call(...)` payload dimensions/types,
5. `np-master` serial parity before any `npRmpi` port.
