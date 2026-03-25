# Deferred `npRmpi` Bounded Adaptive Selector Gap

## Scope

This note tracks the remaining bounded nonfixed public gap after the certified
bounded rollout tranches completed on 2026-03-25.

Out of scope here:

- serial `np` bounded rollout
- already-certified `npRmpi` bounded `generalized_nn`
- already-certified `npRmpi` bounded `adaptive_nn` fit-only routes
- semiparametric/index bounded routes already widened

In scope:

- `npudistbw()` with bounded `bwtype = "adaptive_nn"`
- `npcdistbw()` with bounded `bwtype = "adaptive_nn"`

## Current Public Policy

These selector routes remain intentionally blocked in `npRmpi` with explicit
user-facing errors:

- `bounded adaptive_nn remains unsupported for npudistbw() in npRmpi`
- `bounded adaptive_nn remains unsupported for npcdistbw() in npRmpi`

This is deliberate containment. The estimator-side fit-only routes using
precomputed bounded-adaptive bandwidth objects are already certified and should
remain enabled.

## Why It Is Deferred

Fresh proof tranches showed:

- estimator completion is not the main issue
- the unresolved risk sits in selector/autodispatch/MPI routing for the
  distribution-family bounded adaptive selectors

So this is not a broad bounded-kernel correctness question anymore. It is a
narrow MPI selector-orchestration problem.

## Required Proof-First Reopen Sequence

1. `npudistbw()` first
- build a tiny public reproducer for bounded `adaptive_nn`
- compare:
  - session/spawn
  - attach
  - profile/manual-broadcast
- isolate first failing seam:
  - selector objective route
  - autodispatch wrapper
  - collect/broadcast lifecycle
  - worker-local helper payload

2. `npcdistbw()` second
- only after deciding whether it shares the `npudistbw()` seam

3. Keep fit-only baseline frozen
- do not break the certified fit-only bounded-adaptive routes:
  - `npudist()` with precomputed `dbandwidth`
  - `npcdist()` with precomputed `condbandwidth`

## Acceptance Gate For Any Future Keep

No future widening keep is acceptable unless all of the following remain green:

- session/spawn
- attach
- profile/manual-broadcast
- explicit blocked-route tests are replaced with positive-route tests
- tarball-first `R CMD check --as-cran`

If the selector route requires deeper MPI refactoring, stop and defer rather
than forcing pre-release widening.
