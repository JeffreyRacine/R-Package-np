# Issue #13: npudens variance for categorical-only data

**Status:** Resolved (ae9e41c)

## Root cause
Returned `p/n` instead of `p(1-p)/n` for categorical-only densities.

## Fix
Adjusted variance calculation in `R/np.density.R`.
