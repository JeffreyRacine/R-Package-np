# Issue #7: uocquantile with factor subsets

**Status:** Resolved (ae9e41c)

## Root cause
`table(x)` used full factor levels even after subsetting, causing mismatch with `unique(x)`.

## Fix
Drop unused levels before computing mode in plotting path.
