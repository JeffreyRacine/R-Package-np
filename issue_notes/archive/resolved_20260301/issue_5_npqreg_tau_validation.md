# Issue #5: npqreg should validate tau

**Status:** Resolved (ae9e41c)

## Fix
Added guard to ensure `tau` in (0,1). Returns informative error instead of silent misuse.

## Suggested GitHub response
Acknowledge and note validation added to prevent invalid probability values.
