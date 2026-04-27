# npRmpi NOMAD scaling diagnosis and repair plan

Date: 2026-04-27

## Goal

Repair any real `npRmpi` NOMAD scaling defect surfaced by `npcdens(..., nomad=TRUE)` with bounded kernels, without changing estimator semantics, serial `np` behavior, or unrelated `npRmpi` routes. Confirm whether the same defect affects `npcdist*`; include it only if the evidence shows the same contract is violated there.

## Current evidence

- Serial `np` bounded `npcdens` NOMAD on the exact uniform `n=1000` probe took about 23.3s.
- Installed `npRmpi`, `nslaves=1`, on the same probe took about 13.5s with the same objective, bandwidths, degree, and fast-gate counts. This confirms the large-bandwidth fast gate is active in `npRmpi` on this path.
- A sinusoidal stress probe showed `nslaves=1` and `nslaves=2` both near 90s for `n=2000`, same objective and evaluation counts. CPU sampling showed both coordinator and worker active, so this is not the earlier "worker idle because options were not passed" failure pattern.
- Source inspection shows `npreg` uses an SPMD NOMAD driver: every rank runs the same NOMAD control loop, while each candidate objective evaluation is a collective C evaluation. That can be correct when the C evaluator partitions work and returns synchronized objective values.
- `npcdens` and `npcdist` also run the NOMAD driver on every rank. Therefore, the central question is not "why is the driver duplicated?", but whether each candidate evaluator partitions and synchronizes work correctly across all active ranks, and whether communication/serial work dominates as ranks increase.

## Scientific and engineering contract

1. The selected optimizer semantics remain unchanged: `nomad`/`nomad+powell`, degree search, Bernstein basis, fixed-bandwidth constraints, floor checks, and Powell handoff remain as specified by the user.
2. A given candidate point must produce rank-consistent objective, feasibility, fast-count, and evaluation-count state before the NOMAD driver uses it.
3. When active MPI ranks are present, computational work inside a candidate objective should be partitioned or knowingly collective; ranks must not perform independent full-data work whose results are ignored.
4. Scaling claims must be stated honestly: if the route is limited by fast-gate short-circuiting, rank synchronization, communication, or unavoidable serial phases, the benchmark report must say so.
5. `npcdens` and `npcdist` are checked separately. Shared repair is allowed only where the failing contract is demonstrably shared.
6. No `np-master` mutation is part of this campaign. `np-master` is a read-only reference for semantics and timing.

## Initial implementation plan

1. Preserve raw evidence under `/Users/jracine/Development/tmp/nprmpi_nomad_scaling_20260427`.
2. Build a structural map for `npreg`, `npcdens`, and `npcdist` NOMAD SPMD paths:
   - where the driver runs,
   - whether candidate evaluation is collective,
   - how rank-local objective pieces are reduced or broadcast,
   - how fast/evaluation counts are synchronized.
3. Add a local diagnostic harness that records per-rank participation for `npcdens`/`npcdist` without changing package semantics.
4. If the flaw is rank-result synchronization, repair the C evaluator so rank 0's accepted objective/counts are broadcast consistently, following the established `npreg` shadow pattern.
5. If the flaw is full-data duplication inside a candidate evaluator, repair the evaluator or R dispatch so the candidate work partitions across ranks.
6. If the flaw is not a correctness defect but poor scalability from fast-gate/communication dominance, do not patch semantics; document the measured limitation and preserve the benchmark artifacts.
7. Validate with installed builds:
   - exact uniform `npcdens` fast-gate probe,
   - sinusoidal `npcdens` stress probe with `nslaves=1` and `nslaves=2`,
   - an explicit `nmulti` probe if multistart availability limits observable scaling,
   - corresponding `npcdist` probe if source or timing shows it is affected,
   - a small `npreg` NOMAD sentinel to ensure the known-good SPMD route is not regressed.

## Critique of the initial plan

- It is too quick to treat `npreg` as a direct template. `npreg` broadcasts the rank-0 objective after `bwmfunc_wrapper()`, but that does not prove `npcdens` needs the identical post-wrapper broadcast; some density/distribution CV routines already reduce and broadcast internally.
- A restart-partitioning repair would be attractive but risky. It changes optimizer scheduling, progress accounting, and random-start assignment. It should not be the first repair unless evidence proves candidate-level parallelism is structurally impossible or broken.
- A timing-only diagnosis is insufficient. Equal `nslaves=1` and `nslaves=2` wall time can arise from a real design flaw, a fast-gate-dominated objective, load imbalance, communication overhead, or a small effective parallel fraction.
- `npcdist` should not be swept into the patch merely because it has similar R scaffolding. Its fixed-degree evaluator has a different R/C path and must earn inclusion with source-level or timing evidence.
- The validation set must protect correctness before speed. Matching objective, bandwidths, degree, `num.feval`, and `num.feval.fast` matters as much as elapsed time.

## Reformulated plan

1. Diagnose first, then patch:
   - Compare `npreg`, `npcdens`, and `npcdist` candidate-evaluation synchronization down to the C-level reduction/broadcast contract.
   - Use one lightweight installed-build timing probe per route to separate "active but not scaling" from "idle or serial fallback".
2. Classify the failure:
   - `Correctness/synchronization defect`: ranks can see different candidate objective/count state or rank 0 ignores non-root computed state.
   - `Work-partition defect`: every rank computes the full objective and only one rank's result is used.
   - `Parallel-fraction limitation`: the candidate evaluator is collective and correct, but communication or fast-gate dominance prevents improvement for the tested size.
3. Repair only the classified defect:
   - For synchronization, make `npcdens`/`npcdist` return rank-consistent candidate objective/counts using the narrowest C-level broadcast or allreduce that matches the existing route.
   - For work partitioning, fix the evaluator-level partition/reduction path, not the optimizer's scientific contract.
   - Avoid restart partitioning unless candidate-level parallelism cannot be repaired without larger risk.
4. Validate in gates:
   - Source-level invariants: no `np::` bridge, no serial fallback, no changed defaults.
   - Installed correctness: same fitted result metadata as pre-repair on deterministic seeds.
   - Installed speed: `nslaves=2` should not be slower than `nslaves=1` beyond noise on a route where enough non-fast work remains; if no improvement is expected because the fast gate dominates, record that as a limitation, not a repair failure.
   - `npcdist` included only if affected.
5. Commit only after installed-build gates pass and retained logs clearly support the status claim.

## Comfort and final refinement

I am comfortable proceeding with this refined plan because it narrows the risk axis to candidate-evaluation MPI semantics and explicitly rejects broad optimizer scheduling changes unless the evidence forces them. The final guardrail is to make the smallest measurable repair that restores the existing SPMD contract, then prove it with deterministic installed-build probes before committing.

No further user decision is needed unless diagnosis proves the only meaningful scaling improvement would require a larger redesign, such as objective-subproblem scheduling or restart-level partitioning with changed progress/random-start accounting.

## Execution note

The root cause found during implementation was a candidate-evaluation work-partition defect in the `npcdens` local-polynomial conditional-density CV streams used by NOMAD:

- `cv.ml` LP streaming evaluated every training row on every rank.
- bounded `cv.ls` quadrature LP streaming likewise evaluated every observation block on every rank before any MPI partitioning could apply.
- `npcdist` LP CVLS already had rank-partitioned block streaming with allreduce aggregation and was not the same defect.

The accepted repair should therefore be limited to `npcdens` C-level candidate evaluators: partition the outer row/block work by `my_rank`, aggregate scalar objective terms with `MPI_Allreduce`, and keep serial/local mode unchanged.
