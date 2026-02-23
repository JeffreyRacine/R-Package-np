# `R CMD check --as-cran` Warning/Note Ledger (`np-npRmpi`)

Last refresh: 2026-02-23 (post assignment-control refactors)  
Tarball: `npRmpi_0.70-0.tar.gz`  
Check log: `/tmp/nprmpi_check_ascran_refresh3_20260223.log`

## Current Status
- `Status: 1 WARNING, 2 NOTEs`

## Notes
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.0`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.

2. `future file timestamps` NOTE
- Detail: local environment time-verification limitation (`unable to verify current time`).
- Disposition: accepted for local runs; re-check on release host/CI.

## Warning
1. `top-level files` WARNING (`checkbashisms` script unavailable)
- Detail: environment/tooling warning from local check host.
- Disposition: accepted for local runs; re-check on release host/CI with `checkbashisms` installed.

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression.
