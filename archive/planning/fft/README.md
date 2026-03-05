# Archived FFT Planning Drafts

Canonical active plan:

1. `/Users/jracine/Development/np-master/FFT_EXECUTION_READINESS_V5_2026-03-02.md`

Archived superseded drafts:

1. `FFT_ACCELERATION_PLAN_2026-03-02.md`
2. `FFT_ACCELERATION_DECISION_LOG_2026-03-02.md`
3. `FFT_ALIGNMENT_PLAN_V4_2026-03-02.md`

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
