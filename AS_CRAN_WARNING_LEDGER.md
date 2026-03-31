# `R CMD check --as-cran` Warning/Note Ledger (`np-npRmpi`)

Tracker State: ACTIVE (canonical: `/Users/jracine/Development/ACTIVE_ISSUES_CANONICAL_2026-02-28.md`)

Last refresh: 2026-03-31 (full tarball-first wrapper-hardening closeout refresh)  
Tarball: `npRmpi_0.70-1.tar.gz`  
Check log: `/Users/jracine/Development/tmp/npRmpi_wrapper_string_contracts_20260331/artifacts/stage1_candidate/check_as_cran/00check.log`  
Validation summary: `/Users/jracine/Development/tmp/npRmpi_wrapper_string_contracts_20260331/artifacts/stage1_candidate/validation_summary_2026-03-31.md`

## Current Status
- `Status: 1 WARNING, 1 NOTE` (local full `--as-cran` closeout with tests enabled)

## WARNINGs / NOTEs
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.1`, `existing: 0.60.20`) and stale `Date`.
- Disposition: accepted for the current modernization release line; not a package-code regression.
2. `top-level files` checkbashisms condition
- Local default path emits WARNING:
  - `A complete check needs the 'checkbashisms' script.`
- Disposition: local tooling condition; no package-code action required.

## Effective Local Commands
1. Build tarball:
- `R CMD build /Users/jracine/Development/np-npRmpi`
2. Run full `--as-cran`:
- `env NP_RMPI_ENABLE_ATTACH_TEST=1 NP_RMPI_ENABLE_PROFILE_TEST=1 R CMD check --as-cran /Users/jracine/Development/npRmpi_0.70-1.tar.gz`

## Latest Local Closeout (2026-03-31)
1. Observed check status:
- `Status: 1 WARNING, 1 NOTE`
2. Validation bundle:
- `/Users/jracine/Development/tmp/npRmpi_wrapper_string_contracts_20260331/artifacts/stage1_candidate/validation_summary_2026-03-31.md`
- `/Users/jracine/Development/tmp/npRmpi_wrapper_string_contracts_20260331/artifacts/stage2_live_promotion/validation_summary_2026-03-31.md`
3. Wrapper/test-closeout outcome:
- full `R CMD check --as-cran` passed with tests enabled on the wrapper-hardening closeout tree;
- live repo was promoted to the exact tested tree and rebuilt successfully.

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
