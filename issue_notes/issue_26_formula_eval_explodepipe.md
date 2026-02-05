# Issue #26: Stored/dynamic formula parsing in npplregbw / npcdist

**Status:** Resolved (ae9e41c)

## Root cause
`explodePipe()` used `as.character(formula)` on unevaluated calls.

## Fix
Evaluate formula before splitting (handles stored/dynamic formulas and `~ .`).
