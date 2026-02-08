# Issue #18: Multi-dimensional instrument support for Tikhonov in npregiv

**Status:** Resolved (ae9e41c)

## Root cause
Vector `p` (order of polynomial) used in scalar checks, triggering length errors when `w` has multiple columns.

## Fix
Normalize `p` to scalar at top of Tikhonov path.
