# AGENTS.md

Use `/Users/jracine/Development/AGENTS.md` as the canonical workflow and benchmarking policy.
Hard default: enforce `Sweep Safety Protocol (Default Across Scope Repos)` from the canonical AGENTS file.

Repo-specific note:
- `np-master` is the source of truth for shared C-core performance changes before porting to `np-npRmpi`.
- If asked for an `npRmpi` counterpart, first verify serial behavior in `np-master` and then port with MPI-specific broadcast/start-stop rules from `/Users/jracine/Development/AGENTS.md` and `/Users/jracine/Development/np-npRmpi/benchmarks/perf/oneoff/run_oneoff_suite.R`.
- For regression bandwidth summary work, keep `num.feval` and `num.feval.fast` semantics aligned with `np-npRmpi`.
- Keep kernel discoverability docs aligned via `man/np.kernels.Rd` and cross-links in method `.Rd` pages.
- `npRmpi` parity ports from `np-master` must preserve `npRmpi` runtime independence: no runtime `np::` bridge calls and no silent serial fallback in MPI-selected plot/bootstrap/helper paths.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
