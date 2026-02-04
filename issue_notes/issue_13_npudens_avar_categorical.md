# Issue #13: npudens variance for categorical-only data

**Type:** Bug

## Report
For purely categorical density estimation with `bws=0`, `npudens()` returns variance `p/n`, but the correct variance is `p(1-p)/n`.

Example from issue:
```
model <- npudens(~factor(X), bws=0, ukertype="aitchisonaitken")
p.hat <- unique(fitted(model))
Avar.p.hat <- unique(se(model))^2
```
Expected: `p.hat*(1-p.hat)/n`. Observed: `p.hat/n`.

## Root cause
The C routine `np_density` returns `derr` that corresponds to `p/n` for categorical-only, zero-bandwidth cases. This is correct for continuous kernels but not for pure sample proportions.

## Fix
In `R/np.density.R` (after the `.C("np_density", ...)` call), adjust `derr` for the special case:
- `ncon == 0`, `nord == 0`, `nuno > 0`
- all unordered bandwidths are zero

Use:
```
myout$derr <- sqrt(p * (1 - p) / n)
```
with `p` clipped to [0,1].

## Risk assessment
Low risk: condition is very narrow (categorical-only + bandwidths zero). Continuous or smoothed cases unchanged.

## Status
Implemented locally in `R/np.density.R`. Needs package install to test end-to-end.

## Status
Resolved (categorical-only variance fixed to p*(1-p)/n when bws=0).
