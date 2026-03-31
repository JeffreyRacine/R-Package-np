# `R CMD check --as-cran` Warning/Note Ledger (`np-master`)

Tracker State: ACTIVE (canonical: `/Users/jracine/Development/ACTIVE_ISSUES_CANONICAL_2026-02-28.md`)

Last refresh: 2026-03-31 (full tarball-first semihat closeout refresh)  
Tarball: `np_0.70-1.tar.gz`  
Check log: `/Users/jracine/Development/tmp/np_semihat_fix_20260331_live_promotion/check_as_cran/check.log`  
Run bundle: `/Users/jracine/Development/tmp/np_semihat_fix_20260331_live_promotion`

## Current Status
- `Status: 1 WARNING, 2 NOTEs` (local full `--as-cran` closeout with tests enabled)

## WARNINGs / NOTEs
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.1`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.
2. `future file timestamps` current-time NOTE
- Detail: `unable to verify current time`
- Disposition: environment/network timing condition; not a package-code regression.
3. `top-level files` checkbashisms condition
- Local default path emits WARNING:
  - `A complete check needs the 'checkbashisms' script.`
- Disposition: local tooling condition; no package-code action required.

## Effective Local Commands
1. Build tarball:
- `R CMD build /Users/jracine/Development/np-master`
2. Run full `--as-cran`:
- `R CMD check --as-cran /Users/jracine/Development/np_0.70-1.tar.gz`

## Latest Local Closeout (2026-03-31)
1. Observed check status:
- `Status: 1 WARNING, 2 NOTEs`
2. Validation bundle:
- `/Users/jracine/Development/tmp/np_semihat_fix_20260331_live_promotion/check_as_cran/check.log`
- `/Users/jracine/Development/tmp/np_semihat_fix_20260331_live_promotion/validation_summary_2026-03-31.md`
3. Touched-suite outcome:
- installed `testthat` harness finishes cleanly in the full check path.

## Win-Builder Submission Evidence (2026-03-04, late)
1. Submitted via canonical devtools route:
- `devtools::check_win_release(pkg='/Users/jracine/Development/np-master', quiet=FALSE)`
- `devtools::check_win_oldrelease(pkg='/Users/jracine/Development/np-master', quiet=FALSE)`
2. Submission artifacts:
- `/tmp/winbuilder_submit_20260304_233529/np_win_release_dir_noemail.log`
- `/tmp/winbuilder_submit_20260304_233529/np_win_oldrelease_dir_noemail.log`
3. Submission status:
- submission command exit codes are `0` for both targets;
- upload acceptance evidence observed via FTP queue visibility/disappearance in timed poll:
  - `/tmp/winbuilder_poll_20260304_234306/summary.txt`
  - tarballs present through poll 8 then absent from poll 9 onward in `R-oldrelease`.
- final win-builder check/result disposition remains pending maintainer email receipt (random result directory links are emailed by win-builder).
4. Result ingestion command (when email links arrive):
- `/Users/jracine/Development/winbuilder_ingest_results.sh`

## Win-Builder Result Disposition (2026-03-05)
1. Ingested result links:
- oldrelease (R 4.4.3): `https://win-builder.r-project.org/wUqA3lnnJu1V`
- release (R 4.5.2): `https://win-builder.r-project.org/3gsWddmTIwas`
2. Ingestion bundle:
- `/tmp/winbuilder_results_20260305_user_emails_055344`
3. Observed check status for both targets:
- `Status: 1 ERROR, 1 NOTE`
4. Blocking `ERROR` root cause (from `00install.out`):
- unresolved LAPACK symbols at link time from `linalg.o`:
  - `dgetrf_`, `dgesv_`, `dgetri_`
- install stops with `ERROR: compilation failed for package 'np'`.
5. NOTE interpretation:
- `CRAN incoming feasibility` (version jump) is expected;
- URL reachability checks can emit transient/host-policy outcomes and are non-blocking for this compile failure.
6. Remediation applied (pending win-builder re-run verification):
- `/Users/jracine/Development/np-master/src/Makevars.win` now links:
  - `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`
7. Compiler warning interpretation:
- warnings such as `-Wunused-but-set-variable` and `-Wmaybe-uninitialized` are diagnostics and did not cause this `ERROR` termination.

## Win-Builder Refresh Disposition (2026-03-23)
1. oldrelease (R 4.4.3):
- URL: `https://win-builder.r-project.org/wUqA3lnnJu1V`
2. release refresh:
- URL: `https://win-builder.r-project.org/yJCf5YLo8jFI`
- disposition: `1 NOTE` only; tests pass; compiled code clean.
3. devel refresh:
- URL: `https://win-builder.r-project.org/SP81r7OJ1ZYV`
- disposition: `1 NOTE` only; tests pass; compiled code clean.
4. Evidence source:
- `/Users/jracine/Development/PRE_RELEASE_ISSUES_2026-03-23.md`
- referenced ingest root:
  - `/tmp/winbuilder_results_20260323_np_nprmpi_refresh_1`

## Cleared In This Refresh
1. `top-level files` warning from non-package artifacts (`archive`, `cran-comments.md`) cleared via `.Rbuildignore`.
2. `hidden files/directories` note from local artifacts (`.DS_Store`, `.Rlib`) cleared via `.Rbuildignore`.
3. `inst/doc` PDF size warning cleared via compact-vignette build path.

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
