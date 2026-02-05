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
