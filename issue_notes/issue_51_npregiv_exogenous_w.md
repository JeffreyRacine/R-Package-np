# Issue #51: npregiv exogenous covariates error

Type: Bug

## Report
When providing exogenous covariates `w` to `npregiv`, the call fails with:
`Error in deriv < 0 || deriv > degree : 'length = 2' in coercion to 'logical(1)'`

## Reproduction (current repo + installed to /tmp/np_lib)
```
set.seed(123)
y <- rnorm(100)
x <- rnorm(100)
z <- rnorm(100)
w <- rnorm(100)
np::npregiv(y, z, w, x)
```

## Findings
The error occurs inside `glpreg()` because `deriv` can be a vector; the current
check uses scalar logical operators.

## Proposed Fix
Use vector-safe checks in `R/npregiv.R`:
```
if (p > 0 && (any(deriv < 0) || any(deriv > degree)))
  stop("deriv must lie between 0 and degree")
```

## Status
Verified: after this change, the reproduction completes successfully.

## Risk
Low. Only affects argument validation and does not change computation for valid inputs.

## Status
Resolved (vector-safe `deriv` validation).
