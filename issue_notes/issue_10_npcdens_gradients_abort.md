# Issue #10: npcdens with gradients causes R abort

Type: Bug (reported crash)

## Report
`npcdens(..., gradients=TRUE)` may abort when y is binary and x contains one
continuous and one unordered categorical variable. The report indicates a crash
after bandwidth selection via `npcdensbw(..., bwmethod="cv.ls")`.

## Attempted Reproduction
```
set.seed(1)
DF <- data.frame(
  y = factor(rbinom(50,1,0.5)),
  x1 = runif(50),
  x2 = factor(rbinom(50,1,0.5))
)
BW <- npcdensbw(y ~ x1 + x2, data = DF, bwmethod = "cv.ls", nmulti = 1,
               uykertype = "liracine", cxkertype = "epanechnikov", uxkertype = "liracine")
npcdens(BW)                # ok
npcdens(BW, gradients=TRUE) # ok in this run
```
No crash observed on this machine with the above parameters.

## Status
Not reproducible with the small sample above. Likely data-dependent or related
to particular kernels, sample sizes, or compiled code state.

## Next Steps
If you have a data set that reproduces the crash, we can attach it and
narrow to the specific C routine called by `npcdens(..., gradients=TRUE)`.

## Status
Not reproducible in local tests; no fix applied.
