# `R CMD check --as-cran` Warning/Note Ledger (`np-npRmpi`)

Last refresh: 2026-02-23 (post `np.*test` + `np.sigtest.Rd` alignment checkpoint)  
Tarball: `npRmpi_0.70-0.tar.gz`  
Check log: `/tmp/nprmpi_check_ascran_refresh_20260223.log`

## Current Status
- `Status: 1 WARNING, 1 NOTE`

## Notes
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.0`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.

## Warning
1. `top-level files` WARNING (`checkbashisms` script unavailable)
- Detail: environment/tooling warning from local check host.
- Disposition: accepted for local runs; re-check on release host/CI with `checkbashisms` installed.

## Cleared In This Refresh
1. Prior `future file timestamps` NOTE is no longer present in current local run.
2. Prior `npsigtest.rbandwidth` code/documentation mismatch warning is cleared after `man/np.sigtest.Rd` default alignment.

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression.
