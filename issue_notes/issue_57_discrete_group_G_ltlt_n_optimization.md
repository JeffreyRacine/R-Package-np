# issue: exploit discrete group structure when G << n

status: open (parking lot)

context:
- In CV loops, when discrete predictors induce only `G` unique combinations with `G << n`, we can reduce repeated work by solving at the group level and expanding back to observations.
- Current gating/profile-cache paths already reduce some repeated kernel work, but there is no dedicated end-to-end `G`-group CV objective pathway yet.

next pass scope:
- Add a strict, correctness-first fast path for regression CV where grouped discrete structure can be used without changing numerics.
- Preserve fallback to existing paths when assumptions fail.
