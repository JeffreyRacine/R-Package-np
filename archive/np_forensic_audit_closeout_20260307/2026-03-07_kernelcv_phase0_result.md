# Phase 0 Result: `kernelcv.c` Reusable-Buffer Candidate Rejected

Date: 2026-03-07
Repo: `np-master`
Status: rejected, no code retained
Candidate: reuse the `mean` workspace in `cv_func_regression_categorical_ls_nn()`

## Decision

Do not keep this patch.

## Reason

The candidate satisfied the targeted NN canaries, but it failed the stricter
collateral-damage standard:

1. beneficiary canaries (`generalized_nn`, `adaptive_nn`) showed exact
   bandwidth/objective parity on the focused regression test,
2. the broader paired screen showed numerical drift in the unaffected family,
   specifically fixed-bandwidth regression selectors that should not have been
   touched by this change,
3. the unaffected-family drift reproduced in fresh sessions, so it was not an
   artifact of the benchmark loop ordering or stale buffer state.

Under the current rules, an optimization candidate that perturbs untouched
selectors must be rejected even if the beneficiary path looks numerically clean.

## Key Evidence

Artifacts:

1. `/tmp/np_kernelcv_phase0_20260307_100045`
2. `/tmp/np_kernelcv_phase0_20260307_100310`

Focused reproducer (`bwtype="fixed"`, same seed/data, fresh session):

1. baseline:
   `0.027840098137775372,0.066279913099310236`
   `fval=0.0095975864557935001`
2. candidate:
   `0.027840098137771493,0.093684831770326035`
   `fval=0.0095975864557934949`

This reproducer remained different after preceding NN calls plus explicit
`C_np_release_static_buffers`, which confirms the candidate is not isolated
enough for acceptance.

## Implication

If this concern is revisited, it should be treated as a deeper investigation
into hidden coupling/undefined behavior in the regression selector stack, not as
a simple allocation micro-optimization.
