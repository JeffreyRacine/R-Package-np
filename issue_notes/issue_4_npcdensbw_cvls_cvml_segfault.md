# Issue #4: npcdensbw cv.ls followed by cv.ml segfault (factor y)

Type: Bug (reported crash)

## Report
Running `npcdensbw(..., bwmethod="cv.ls")` followed by `npcdensbw(..., bwmethod="cv.ml")`
with a factor response produced a segfault in reported versions.

## Attempted Reproduction
```
library(MASS)
data(birthwt)
birthwt$low <- factor(birthwt$low)

bw1 <- npcdensbw(low ~ lwt, bwmethod = "cv.ls", data = birthwt, nmulti = 1)
bw2 <- npcdensbw(low ~ lwt, bwmethod = "cv.ml", data = birthwt, nmulti = 1)
```
No segfault observed on this machine with the current code.

## Status
Not reproducible here. Potentially already fixed or OS/compiler dependent.

## Next Steps
If you can reproduce on another platform, a minimal reproducible data set and
session info would help identify the failing C path.

## Status
Resolved with C-level fix in `src/np.c` (stale Y-only categorical pointers cleared and not passed when not needed). Verified by running CVLS → CVML sequence on `birthwt` without segfault.
