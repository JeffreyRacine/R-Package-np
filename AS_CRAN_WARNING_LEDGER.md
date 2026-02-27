# `R CMD check --as-cran` Warning/Note Ledger (`np-npRmpi`)

Last refresh: 2026-02-27 (strict `Rmpi` decoupling refresh)  
Tarball: `npRmpi_0.70-1.tar.gz`  
Check log: `/tmp/nprmpi_check_notimeout_20260227_071702.log`

## Current Status
- `Status: 1 WARNING, 2 NOTEs`

## WARNINGs
1. `top-level files` WARNING
- Detail: local environment does not provide `checkbashisms` by default.
- Message: `A complete check needs the 'checkbashisms' script.`
- Disposition: accepted for local runs in this environment; for full local parity, use the compatibility wrapper path described below.

## NOTEs
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.1`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.
2. `future file timestamps` NOTE
- Detail: `unable to verify current time`.
- Disposition: accepted as environment/transient; no package-code action.

## Effective Local Commands
1. Build tarball:
- `R CMD build /Users/jracine/Development/np-npRmpi`
2. Run `--as-cran` with MPI env:
- `FI_TCP_IFACE=en0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=en0 R CMD check --as-cran npRmpi_0.70-1.tar.gz`
- fallback: `FI_TCP_IFACE=lo0` if needed
3. Optional shell lint parity (`checkbashisms` compatibility):
- `PATH="/tmp:$PATH" ... R CMD check --as-cran ...`
- wrapper path: `/tmp/checkbashisms` delegating to `shellcheck` with `-p/-d` support.

## Decoupling Notes
1. `Rmpi` is intentionally removed from `Suggests` in `DESCRIPTION`.
2. Rd cross-reference dependency on external `Rmpi` docs was removed (`man/hosts.Rd` now points to local `npRmpi.init` usage context).
3. Session/runtime version reporting now resolves through embedded/internal `npRmpi` interfaces (no external `Rmpi` namespace lookup).

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression unless explicitly accepted.
