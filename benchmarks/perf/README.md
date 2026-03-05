# Performance Benchmarks

Use this folder for timing-focused scripts and summary utilities.

- `methods/`: estimator-family benchmark suites (`npreg`, `npcdens`, etc.).
- `oneoff/`: one-off benchmark harnesses for methods not covered by combo grids.

## Output Conventions

- Write raw and summary artifacts to `/tmp`.
- Include key timing/objective columns (`elapsed`, objective value where relevant,
  selected bandwidths, and `num.feval` when available).
- For post-vs-pre reports, use positive percent values to mean "post is faster".

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
