# Issue #51: npregiv with exogenous covariates w

**Status:** Resolved (ae9e41c)

## Root cause
Vector length error in derivative check when w is provided.

## Fix
Make derivative check vector‑safe; allow exogenous covariates.
