# `R CMD check --as-cran` Warning/Note Ledger (`np-npRmpi`)

Last refresh: 2026-02-23 (post shellcheck-wrapper checkpoint)  
Tarball: `npRmpi_0.70-0.tar.gz`  
Check log: `/tmp/nprmpi_check_ascran_compact_shellcheck_20260223.log`

## Current Status
- `Status: 1 NOTE`

## Notes
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.0`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.

## Effective Local Commands
1. Build compacted tarball:
- `R CMD build --compact-vignettes=both /Users/jracine/Development/np-npRmpi`
2. Run `--as-cran` with local wrapper and MPI env:
- `PATH="/tmp:$PATH" FI_TCP_IFACE=en0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=en0 R CMD check --as-cran npRmpi_0.70-0.tar.gz`
- fallback: `FI_TCP_IFACE=lo0` if needed
- wrapper path: `/tmp/checkbashisms` (shellcheck-backed compatibility shim with `-p/-d` support)

## Cleared In This Refresh
1. Prior `future file timestamps` NOTE is no longer present in current local run.
2. Prior `npsigtest.rbandwidth` code/documentation mismatch warning is cleared after `man/np.sigtest.Rd` default alignment.
3. Prior `top-level files` `checkbashisms` warning is cleared in local runs by supplying `/tmp/checkbashisms` on `PATH`.

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression.
