# `R CMD check --as-cran` Warning/Note Ledger (`np-npRmpi`)

Tracker State: ACTIVE (canonical: `/Users/jracine/Development/ACTIVE_ISSUES_CANONICAL_2026-02-28.md`)

Last refresh: 2026-03-04 (release-hygiene follow-up)  
Tarball: `npRmpi_0.70-1.tar.gz`  
Check log: `/Users/jracine/Development/npRmpi.Rcheck/00check.log`  
Run bundles:
- default local: `/tmp/release_hygiene_20260304_231716`
- wrapper-assisted: `/tmp/release_hygiene_followup_20260304_232456`

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
1. Build tarball:
- `R CMD build /Users/jracine/Development/np-npRmpi`
2. Run `--as-cran` (default local baseline):
- `R CMD check --as-cran /Users/jracine/Development/npRmpi_0.70-1.tar.gz`
3. Optional wrapper-assisted run:
- `PATH="/tmp:$PATH" R CMD check --as-cran /Users/jracine/Development/npRmpi_0.70-1.tar.gz`
- wrapper path: `/tmp/checkbashisms` (shellcheck-backed compatibility shim with `-p/-d` acceptance).

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

## Decoupling Notes
1. `Rmpi` is intentionally removed from `Suggests` in `DESCRIPTION`.
2. Rd cross-reference dependency on external `Rmpi` docs was removed (`man/hosts.Rd` points to local `npRmpi.init` usage context).
3. Session/runtime version reporting resolves through embedded/internal `npRmpi` interfaces (no external `Rmpi` namespace lookup).

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression unless explicitly accepted.
