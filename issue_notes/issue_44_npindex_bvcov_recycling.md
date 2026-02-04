# Issue #44: Bvcov in npindex.sibandwidth (matrix recycling)

Type: Bug (potential numerical correctness)

## Report
Matrix/vector recycling in `np.singleindex.R` may multiply the gradient vector
by columns rather than rows, leading to incorrect `Bvcov` when `k > 1`.

## Location
`R/np.singleindex.R` around the `Bvcov` calculation in `npindexbw.sibandwidth`.

## Proposed Fix
Replace column-wise recycling with explicit row-wise scaling:
```
# old:
# dg.db.xmex <- index.tgrad[,1] * xmex
# Sigma <- (uhat * dg.db.xmex) %*% t(uhat * dg.db.xmex)

# new:
# multiply each column (of t(xmex)) by index.tgrad
# then scale rows by uhat before crossprod
 dg.db.xmex <- t(index.tgrad[,1] * t(xmex))
 Sigma <- t(uhat * t(dg.db.xmex)) %*% t(uhat * t(dg.db.xmex))
```

## Verification
Ran a small test:
```
set.seed(123)
n <- 200
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)
npindexbw(y ~ x1 + x2, method = "kleinspady", gradients = TRUE)
```
No error observed.

## Risk / Caveats
- Changes the numerical value of `Bvcov` for `k > 1`.
- If `dg.db.xmex %*% t(dg.db.xmex)` is not positive-definite, `chol()` can fail
  (this is inherent to the prior computation as well, but the change may surface
  it in different cases).
- Consider a guarded fallback (e.g., use `chol()` inside `tryCatch` or add a small
  ridge) if failures are observed.

## Status
Implemented in `R/np.singleindex.R` (row-wise scaling). Needs confirmation that resulting `Bvcov` is numerically stable across larger samples.

## Update
Reverted: The row-wise scaling change caused non-PD errors in `chol()` for some fits
(e.g., `npindex(..., gradients=TRUE)`). Leaving as-is pending co-author review.
Regression cause identified: missing `xmex` vector guard in 0.60-21; restored baseline handling from 0.60-20.
