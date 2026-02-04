# Issue #18: Multi-dimensional instrument support for Tikhonov in npregiv

Type: Feature request

## Report
Request to support multi-dimensional instruments in the Tikhonov (npregiv)
implementation.

## Status
Feature request, not a bug. Requires design and likely C-level changes.

## Status
Feature request; no code changes proposed yet.

## Guidance
If multi-dimensional instruments are supported for Landweber–Fridman, extending Tikhonov may be feasible, but requires careful review of the Tikhonov implementation and its interface. No code changes yet.

## Update
Confirmed: Tikhonov path failed with multi-dimensional instruments because `Kmat.lp` assumed scalar `p` and used `if (p < 0)`.
Fixed by normalizing `p` to a scalar inside `Kmat.lp` (uses first element), so calls that pass `p=rep(p, k)` no longer error.

## Status
Resolved. Verified Tikhonov with multi-dimensional `w` now runs (see user repro).
