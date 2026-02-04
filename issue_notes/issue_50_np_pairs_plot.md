# Issue #50: Cross-Validated Pairs Plots

**Type:** Feature request

## Request
Add helper functions to compute bandwidths for all pairs of variables and plot density/conditional mean pairs using cross-validated bandwidths.

## Proposed implementation
Added two functions locally:
- `np.pairs(y_vars, y_dat, ...)` builds pairwise density/regression objects.
- `np.pairs.plot(pair_list)` plots density on diagonal and regression off-diagonal.

Implementation stored in `R/np.pairs.R` (local).

## Reproduction / demo
Using USArrests example from issue:
```
y_vars <- c('Murder','UrbanPop')
names(y_vars) <- c('Murder Arrests per 100K', 'Pop.PercentUrban')
pair_list <- np.pairs(y_vars=y_vars, y_dat=USArrests,
                      ckertype='epanechnikov', bwscaling=TRUE)
np.pairs.plot(pair_list)
```

## Risk assessment
Low risk (new feature). Needs:
- NAMESPACE export entries
- man page (`man/np.pairs.Rd`)
- Any naming conflict checks

## Status
Implemented locally in `R/np.pairs.R` for evaluation. Not wired into NAMESPACE or docs yet.

## Status
Feature prototype added in `R/np.pairs.R`, not yet wired into NAMESPACE/docs.

## Guidance
Would require a new `.Rd` page with the user-provided example (likely under `\dontrun{}`), guardrails for inputs, and tests. Not implemented yet.

## Status
Implemented with new \code{np.pairs} and \code{np.pairs.plot} functions, exported in NAMESPACE, and documented in \code{man/np.pairs.Rd} with a \dontrun{} example.
