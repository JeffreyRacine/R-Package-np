## R CMD check results

Local (macOS, R 4.5.2, --as-cran):

0 errors | 2 warnings | 2 notes

### Notes

1. `checking CRAN incoming feasibility ... NOTE`
- Version jump from CRAN `0.60-20` to submission `0.70-1` is intentional for a major modernization release.

2. `checking for future file timestamps ... NOTE`
- Environment clock verification note (`unable to verify current time`) observed intermittently on local host.

### Warnings

1. `checking top-level files ... WARNING`
- Non-standard top-level directory `archive/` retained intentionally for project provenance/history.

2. `checking sizes of PDF files under inst/doc ... WARNING`
- Vignette PDFs can be further compacted by `tools::compactPDF(...)`; content is unchanged and installed files are valid.

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
