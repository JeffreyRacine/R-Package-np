# Issue #4: npcdensbw cv.ls → cv.ml segfault when y is factor

**Status:** Resolved (in commit ae9e41c)

## Repro
```r
library(MASS); data(birthwt)
birthwt$low <- factor(birthwt$low)

library(np)
nmulti <- 1
bw1 <- npcdensbw(low ~ lwt, bwmethod="cv.ls", data=birthwt, nmulti=nmulti)
bw2 <- npcdensbw(low ~ lwt, bwmethod="cv.ml", data=birthwt, nmulti=nmulti)
```

## Root cause
C state for conditional density used stale pointers to categorical structures between calls.

## Fix
Reset conditional‑density externals to NULL at entry and after free in `src/np.c` (np_density_conditional_bw / np_distribution_conditional_bw). Prevents stale pointer reuse.

## Suggested GitHub response (draft)
- Confirm segfault and fix in C.
- Provide reproducible check; confirm `cv.ls` then `cv.ml` works.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
