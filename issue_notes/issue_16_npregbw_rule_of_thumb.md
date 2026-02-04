# Issue #16: Rule of thumb for npregbw

Type: Feature request

## Report
Request to add a plug-in (rule-of-thumb) bandwidth method for regression
bandwidth selection (npregbw), similar to density plug-in methods.

## Notes
- Current `npregbw()` supports `bwmethod = "cv.ls"` and `"cv.aic"`.
- There is no standard multivariate regression rule-of-thumb in base R.

## Status
Feature request, not a bug. Requires design decision and documentation.

## Status
Feature request; no code changes proposed yet.

## Guidance
Users can supply manual bandwidths via `bws=`. A general rule-of-thumb for multivariate mixed-data regression is not available/standard, so no generic plug-in method is proposed.

## Status
Resolved via guidance (manual bandwidths supported; no general multivariate rule-of-thumb).
