# Issue #44: Bvcov in npindex.sibandwidth

**Status:** Open (no code change retained)

## Notes
A proposed change to multiplication order was tested but caused instability in `npindex(..., gradients=TRUE)` (chol failure). Reverted to baseline pending co‑author review and deeper theoretical check.
