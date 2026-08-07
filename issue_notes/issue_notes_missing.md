# Issue Notes Index (`np-npRmpi`)

Last refresh: 2026-08-07

## Active / Keep In Place

- `issue_notes/issue_3_local_linear_high_order.md` (open methodological issue)

## Resolved / Historical

Resolved plans and status notes were removed from the live source tree during
the 2026-08-07 pre-release modernization. Their exact contents remain in Git
history and in the dated workspace closeout archives. Canonical current policy
and release gates live at the workspace root; this directory is not a second
status tracker.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
