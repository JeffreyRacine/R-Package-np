# `R CMD check --as-cran` Warning/Note Ledger (`np-master`)

Last refresh: 2026-02-23 (post compact-vignette + shellcheck-wrapper checkpoint)  
Tarball: `np_0.70-0.tar.gz`  
Check log: `/tmp/np_master_check_ascran_compact_shellcheck_20260223.log`

## Current Status
- `Status: 1 NOTE`

## Notes
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.0`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.

## Effective Local Commands
1. Build compacted tarball:
- `R CMD build --compact-vignettes=both /Users/jracine/Development/np-master`
2. Run `--as-cran` with local wrapper:
- `PATH="/tmp:$PATH" R CMD check --as-cran np_0.70-0.tar.gz`
- wrapper path: `/tmp/checkbashisms` (shellcheck-backed compatibility shim with `-p/-d` support)

## Cleared In This Refresh
1. Prior `npsigtest.rbandwidth` code/documentation mismatch warning is cleared after `man/np.sigtest.Rd` default alignment.
2. Prior `top-level files` `checkbashisms` warning is cleared in local runs by supplying `/tmp/checkbashisms` on `PATH`.
3. Prior `inst/doc` PDF size warning is cleared by building with `--compact-vignettes=both`.

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression.
