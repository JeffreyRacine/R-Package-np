# npRmpi Autodispatch Refactor Plan (Staged, Safe)

## Objective
Replace patch-accumulated autodispatch orchestration with a simpler, maintainable execution core while preserving current working behavior.

## Constraints
1. No regressions for currently passing core workflows (`npreg*`, `npcdens*`, `npcdist*`).
2. No new hangs: fail-fast beats implicit deadlock.
3. Keep manual `mpi.bcast.cmd(..., caller.execute=TRUE)` workflows functional.
4. Do not break existing examples during transition.

## What We Keep
1. Existing preflight guards (`npRmpi.start` required, formula+no-data fail-fast).
2. Existing option sync behavior across ranks.
3. Existing manual-mode safety for `npsigtest`.

## Target Architecture
1. One internal executor for autodispatch calls.
2. Explicit execution mode tagging on returned objects (`auto` vs `manual`).
3. Centralized mode guards in method entry points (predict/plot family).
4. Minimal per-function wrapper logic: wrappers call executor and return.

## Execution Phases

### Phase 1: Stabilize + Instrument (no behavior change)
1. Add internal execution diagnostics (optional) for mode and call class.
2. Add lifecycle mode tags to returned objects for core constructors.
3. Add mode-mismatch fail-fast messages in selected high-risk methods (`plot`, key `predict.*`).

Acceptance:
- Current passing smokes remain green.
- Known failures remain only `np.plot` and `npsigtest` pathways.

### Phase 2: Introduce New Executor Behind Feature Flag
1. Add `options(npRmpi.autodispatch.engine = "legacy"|"payload")`.
2. Implement payload-engine executor path without touching default.
3. Keep legacy path fully intact and default.
4. Route only core constructor wrappers through switchable engine.

Acceptance (payload engine only):
- `npregbw/npreg`, `npcdensbw/npcdens`, `npcdistbw/npcdist` pass fixed smokes.
- direct call parity vs legacy on matched seeds for bw/fit outputs.

### Phase 3: Method-chain Migration
1. Migrate key `predict.*` chain to payload engine under flag.
2. Keep `plot` guarded/manual-safe until predict chain proves stable.
3. Keep `npsigtest` manual-only for now.

Acceptance:
- Core predict-chain smokes pass in payload mode.
- No new deadlocks in mixed call sequences.

### Phase 4: Resolve Plot + npsigtest
1. `plot`: choose one strict policy and enforce consistently.
   - Preferred near-term: manual-only for MPI compute path with explicit fail-fast guidance.
2. `npsigtest`: keep explicit/manual path and harden, or add dedicated bootstrap-safe executor branch.

Acceptance:
- `np*.Rd` dontrun sweep failures are reduced to 0, or documented intentional guards with revised examples.

### Phase 5: Promote Payload Engine
1. If all acceptance gates pass, switch default engine to `payload`.
2. Keep legacy engine one release cycle behind fallback option.
3. Remove legacy path only after full cycle of stable checks.

## Test Matrix (Required at each Phase Gate)
1. Core constructors: `npreg*`, `npcdens*`, `npcdist*`.
2. Method chain: `predict` then `plot` (where policy allows).
3. Start/stop cycles and repeated session runs.
4. Load-order checks (`np`, `Rmpi`, `npRmpi`).
5. `np*.Rd` dontrun harness (subprocess-isolated, per-block timeout).

## Rollback Rule
If any phase introduces hangs or broad regressions:
1. Keep default on legacy engine.
2. Isolate failure with per-block subprocess harness.
3. Revert phase patch; do not stack speculative fixes.

## Immediate Next Implementation Step
Phase 1 only:
1. add mode tagging and mode-mismatch guardrails,
2. add diagnostics toggles,
3. rerun core smokes and dontrun harness subset.
