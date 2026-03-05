# `R CMD check --as-cran` Warning/Note Ledger (`np-master`)

Tracker State: ACTIVE (canonical: `/Users/jracine/Development/ACTIVE_ISSUES_CANONICAL_2026-02-28.md`)

Last refresh: 2026-03-05 (win-builder result ingestion + Windows link remediation)  
Tarball: `np_0.70-1.tar.gz`  
Check log: `/Users/jracine/Development/np.Rcheck/00check.log`  
Run bundle: `/tmp/release_hygiene_followup_20260304_232456`

## Current Status
- `Status: 1 NOTE`

## Notes
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.1`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.

## Effective Local Commands
1. Build compacted tarball:
- `R CMD build --compact-vignettes=both /Users/jracine/Development/np-master`
2. Run `--as-cran`:
- `R CMD check --as-cran /Users/jracine/Development/np_0.70-1.tar.gz`

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

## Cleared In This Refresh
1. `top-level files` warning from non-package artifacts (`archive`, `cran-comments.md`) cleared via `.Rbuildignore`.
2. `hidden files/directories` note from local artifacts (`.DS_Store`, `.Rlib`) cleared via `.Rbuildignore`.
3. `inst/doc` PDF size warning cleared via compact-vignette build path.

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression.
