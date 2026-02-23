# `R CMD check --as-cran` Warning/Note Ledger (`np-master`)

Last refresh: 2026-02-23  
Tarball: `np_0.70-0.tar.gz`  
Check log: `/tmp/np_master_check_ascran_20260223.log`

## Current Status
- `Status: 2 WARNINGs, 1 NOTE`

## Notes
1. `CRAN incoming feasibility` NOTE
- Detail: version jump (`submitted: 0.70.0`, `existing: 0.60.20`).
- Disposition: accepted for modernization release line.

## Warnings
1. `top-level files` WARNING (`checkbashisms` script unavailable)
- Detail: environment/tooling warning from local check host.
- Disposition: accepted for local runs; re-check on release host/CI with `checkbashisms` installed.

2. `inst/doc` PDF size WARNING
- Detail: `gs+qpdf` suggests compaction opportunities for generated vignette PDFs.
- Disposition: accepted for now; optional packaging optimization before final CRAN submission (`tools::compactPDF` or `--compact-vignettes`).

## Gate Policy
- Any *new* warning/note not listed here is treated as a regression.
