# Issue #6: npindexbw / npindex bandwidth type option

Type: Bug (argument ignored)

## Report
`npindexbw(..., method="kleinspady", bwtype="adaptive_nn")` returned `type: fixed`.

## Cause
In `R/np.singleindex.bw.R`, `npindexbw.sibandwidth()` rebuilt the bandwidth
object with `sibandwidth(...)` but did not pass through `bwtype`, so it defaulted
back to `fixed`.

## Fix
Pass `bwtype = bws$type` when rebuilding:
```
 bws <- sibandwidth(..., bwtype = bws$type, ...)
```

## Validation Script
```
library(np)
set.seed(123)
n <- 200
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)

bw_fixed <- npindexbw(y ~ x1 + x2, method = "kleinspady", gradients = TRUE, bwtype = "fixed")
bw_adapt <- npindexbw(y ~ x1 + x2, method = "kleinspady", gradients = TRUE, bwtype = "adaptive_nn")
bw_gen <- npindexbw(y ~ x1 + x2, method = "kleinspady", gradients = TRUE, bwtype = "generalized_nn")

stopifnot(bw_fixed$type == "fixed")
stopifnot(bw_adapt$type == "adaptive_nn")
stopifnot(bw_gen$type == "generalized_nn")

cat("OK: bwtype propagation works\n")
```

## Status
Resolved. `bwtype` now propagates into the rebuilt `sibandwidth` object.
