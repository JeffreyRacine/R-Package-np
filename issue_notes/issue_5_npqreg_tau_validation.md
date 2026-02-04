# Issue #5: npqreg() does not validate tau

**Type:** Bug

## Report
`npqreg()` accepts invalid `tau` values (e.g., <0 or >1), which can lead to nonsensical results or errors downstream.

## Fix
Add explicit validation in `npqreg.condbandwidth`:
```
if (!is.numeric(tau) || length(tau) != 1 || is.na(tau) || tau <= 0 || tau >= 1)
  stop("'tau' must be a single numeric value in (0,1)")
```

## Risk assessment
Low. This is defensive input validation, only affects invalid inputs.

## Status
Implemented locally in `R/np.qregression.R`. Needs package install/test.

## Status
Resolved (explicit validation for `tau` in (0,1)).
