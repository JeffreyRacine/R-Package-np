# Method Suites

This folder contains estimator-family performance suites, one directory per
method.

## Method Index

- `npreg/`: `run_npreg_combos.R`, `bench_npreg_param_nprmpi.R`
- `npcdens/`: `run_npcdens_combos.R`, `bench_npcdens_param_nprmpi.R`
- `npudens/`: `run_npudens_combos.R`, `bench_npudens_param_nprmpi.R`
- `npcdist/`: `run_npcdist_combos.R`, `bench_npcdist_param_nprmpi.R`
- `npudist/`: `run_npudist_combos.R`, `bench_npudist_param_nprmpi.R`
- `npplreg/`: `run_npplreg_combos.R`, `bench_npplreg_param_nprmpi.R`
- `npscoef/`: `run_npscoef_combos.R`, `bench_npscoef_param_nprmpi.R`
- `npindex/`: `run_npindex_combos.R`, `bench_npindex_param_nprmpi.R`

Each method directory also contains report/comparison helpers.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
