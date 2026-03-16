# `R CMD check --as-cran` Warning/Note Ledger (`np-npRmpi`)

Tracker State: ACTIVE (canonical: `/Users/jracine/Development/ACTIVE_ISSUES_CANONICAL_2026-02-28.md`)

Last refresh: 2026-03-16 (no-vignette tarball-first closeout refresh)  
Tarball: `npRmpi_0.70-1.tar.gz`  
Check log: `/tmp/nprmpi_fix_check.log`  
Run bundles:
- default local: `/tmp/nprmpi_fix_check.log`
- helper/static proof: `/tmp/nprmpi_fix5_static.log`
- exact tiny-support proof: `/tmp/nprmpi_fix5_tiny_support.log`

## Current Status
- Default local (`PATH` without `checkbashisms` wrapper): `Status: 1 WARNING, 1 NOTE`
- Wrapper-assisted (`PATH="/tmp:$PATH"` with `/tmp/checkbashisms` shim): `Status: 2 NOTEs`

## WARNINGs / NOTEs
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.1`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.
2. `top-level files` checkbashisms condition:
- Default local path emits WARNING:
  - `A complete check needs the 'checkbashisms' script.`
- Wrapper-assisted path removes that WARNING but may emit a NOTE from shellcheck parsing `configure.ac` m4 syntax (SC1064/SC1065/SC1072/SC1073).
- Disposition: both are local tooling conditions; no package-code action required.

## Effective Local Commands
1. Build no-vignette tarball:
- `R CMD build --no-build-vignettes --no-manual /Users/jracine/Development/np-npRmpi`
2. Run `--as-cran` (default local baseline):
- `R CMD check --as-cran --ignore-vignettes /Users/jracine/Development/npRmpi_0.70-1.tar.gz`
3. Optional wrapper-assisted run:
- `PATH="/tmp:$PATH" R CMD check --as-cran --ignore-vignettes /Users/jracine/Development/npRmpi_0.70-1.tar.gz`
- wrapper path: `/tmp/checkbashisms` (shellcheck-backed compatibility shim with `-p/-d` acceptance).

## Latest Local Closeout (2026-03-16)
1. Observed check status:
- default local no-vignette baseline: `Status: 1 WARNING, 1 NOTE`
2. Validation bundle:
- `/tmp/nprmpi_fix_check.log`
- `/tmp/nprmpi_fix5_static.log`
- `/tmp/nprmpi_fix5_tiny_support.log`
3. Touched-suite outcome:
- touched static contracts passed and the exact tiny-support helper proof returned `TINY_SUPPORT_OK 0e+00`.

## Win-Builder Submission Evidence (2026-03-04, late)
1. Submitted via canonical devtools route:
- `devtools::check_win_release(pkg='/Users/jracine/Development/np-npRmpi', quiet=FALSE)`
- `devtools::check_win_oldrelease(pkg='/Users/jracine/Development/np-npRmpi', quiet=FALSE)`
2. Submission artifacts:
- `/tmp/winbuilder_submit_20260304_233529/nprmpi_win_release_dir_noemail.log`
- `/tmp/winbuilder_submit_20260304_233529/nprmpi_win_oldrelease_dir_noemail.log`
3. Submission status:
- submission command exit codes are `0` for both targets;
- upload acceptance evidence observed via FTP queue visibility/disappearance in timed poll:
  - `/tmp/winbuilder_poll_20260304_234306/summary.txt`
  - tarballs present through poll 8 then absent from poll 9 onward in `R-oldrelease`.
- final win-builder check/result disposition remains pending maintainer email receipt (random result directory links are emailed by win-builder).
4. Result ingestion command (when email links arrive):
- `/Users/jracine/Development/winbuilder_ingest_results.sh`

## Win-Builder Result Disposition (2026-03-05)
1. Ingested result link:
- oldrelease (R 4.4.3): `https://win-builder.r-project.org/whIxNg5rQr8Q`
2. Ingestion bundle:
- `/tmp/winbuilder_results_20260305_user_emails_055344`
3. Observed check status:
- `Status: 1 ERROR, 1 NOTE`
4. Blocking `ERROR` root cause (from `00install.out`):
- unresolved LAPACK symbols at link time from `linalg.o`:
  - `dgetrf_`, `dgesv_`, `dgetri_`
- install stops with `ERROR: compilation failed for package 'npRmpi'`.
5. NOTE interpretation:
- `CRAN incoming feasibility` (version jump) is expected;
- URL checks to `mpich.org`/`sshrc-crsh.gc.ca` can return host-policy/transient responses and are not this compile blocker.
6. Remediation applied (pending win-builder re-run verification):
- `/Users/jracine/Development/np-npRmpi/src/Makevars.win` now appends:
  - `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` to both WIN64 and WIN32 link lines.
7. Compiler warning interpretation:
- warnings such as `-Wmaybe-uninitialized` are diagnostics and did not cause this `ERROR` termination.

## Decoupling Notes
1. `Rmpi` is intentionally removed from `Suggests` in `DESCRIPTION`.
2. Rd cross-reference dependency on external `Rmpi` docs was removed (`man/hosts.Rd` points to local `npRmpi.init` usage context).
3. Session/runtime version reporting resolves through embedded/internal `npRmpi` interfaces (no external `Rmpi` namespace lookup).

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression unless explicitly accepted.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
