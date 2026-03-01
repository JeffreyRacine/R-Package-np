# Diff vs CRAN np_0.60-19

This report provides a brief narrative of what changed in each file compared to the CRAN posted package.

- CRAN posted: /Users/jracine/Development/np_0.60-19_CRAN_posted
- Current: /Users/jracine/Development/np-master

## Summary

Added files: 24
Removed files: 0
Changed files: 23

## Added

- `BUILD.md`: Build instructions added.
- `CRAN_DIFF.md`: Added in the current codebase (new file).
- `R/np.pairs.R`: Added cross‑validated pairs plotting utilities (np.pairs / np.pairs.plot) for quick bivariate diagnostics. (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)
- `WORKTREES.md`: Added in the current codebase (new file).
- `np_0.60-21.tar.gz`: Added in the current codebase (new file).
- `tests/testthat.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-extra.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-misc.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Expand unit test coverage to include distribution, copula, quantile regression, tests of independence/specification, and boundary/shape-constrained methods)
- `tests/testthat/test-npcdens.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-npcdist.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-npcopula.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Expand unit test coverage to include distribution, copula, quantile regression, tests of independence/specification, and boundary/shape-constrained methods)
- `tests/testthat/test-npdeptest.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Expand unit test coverage to include distribution, copula, quantile regression, tests of independence/specification, and boundary/shape-constrained methods)
- `tests/testthat/test-npindex.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-npplreg.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-npqreg.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Expand unit test coverage to include distribution, copula, quantile regression, tests of independence/specification, and boundary/shape-constrained methods)
- `tests/testthat/test-npreg.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-npregiv.R`: Added in the current codebase (new file). (last change: 2026-02-02 | Update unit tests and fix documentation mismatches for npregiv/npregivderiv)
- `tests/testthat/test-npscoef.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-npsigtest.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-nptests.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Expand unit test coverage to include distribution, copula, quantile regression, tests of independence/specification, and boundary/shape-constrained methods)
- `tests/testthat/test-npudens.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `tests/testthat/test-npudist.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Expand unit test coverage to include distribution, copula, quantile regression, tests of independence/specification, and boundary/shape-constrained methods)
- `tests/testthat/test-npuniden.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Expand unit test coverage to include distribution, copula, quantile regression, tests of independence/specification, and boundary/shape-constrained methods)
- `tests/testthat/test-utils.R`: Added in the current codebase (new file). (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)

## Changed

- `.Rbuildignore`: Updated relative to CRAN (see file for details). (last change: 2014-01-06 | Updated .Rbuildignore for files in man, changes to np.quantile.Rd)
- `DESCRIPTION`: Version/metadata updates. (last change: 2026-02-04 | Bumped up version # to 0.60-21 in anticipation of moving to open issues)
- `NAMESPACE`: Export/method registration updates for new/changed functions. (last change: 2026-02-02 | Implement S3 methods (print, summary, plot) for npregiv and npregivderiv)
- `R/np.condensity.bw.R`: Minor conditional density bandwidth selection fix (consistency with updated inputs). (last change: 2026-02-05 | Fix np package registrations after sync)
- `R/np.density.R`: Density estimation tweaks and guardrails (minor logic adjustments for stability/consistency). (last change: 2026-02-05 | Fix np package registrations after sync)
- `R/np.plot.R`: Plotting updates to align with revised outputs and defaults. (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)
- `R/np.qregression.R`: Quantile regression path adjusted for consistency with updated bandwidth/fit handling. (last change: 2026-02-05 | Fix np package registrations after sync)
- `R/np.sigtest.R`: Adjusted significance test workflow and guardrails (fixes error when training/eval data frames are not aligned; supports updated regression calls). (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)
- `R/np.singleindex.R`: Single‑index model updates (minor logic/guard changes). (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)
- `R/np.singleindex.bw.R`: Bandwidth selection updates for single‑index models (consistency checks and defaults). (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)
- `R/npregiv.R`: Instrumental variables regression updates (fixes/guards in evaluation and fit paths). (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)
- `R/npregivderiv.R`: npregiv derivative computation and summaries aligned with updated argument ordering. (last change: 2026-02-03 | Add bandwidth details to npregiv summaries)
- `R/util.R`: Utility updates and numerical safeguards (NZD helpers and related checks) aligned with newer codepaths. (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)
- `R/zzz.R`: Startup/initialization tweaks (package message and load behavior alignment). (last change: 2026-02-04 | Bumped up version # to 0.60-21 in anticipation of moving to open issues)
- `README.md`: Updated relative to CRAN (see file for details). (last change: 2026-01-27 | # Commit Message)
- `demo/Engel95.R`: Demo updates for Engel95 dataset usage. (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `demo/npregiv.R`: Demo updates for npregiv usage. (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `man/data-Engel95.Rd`: Engel95 dataset documentation updates. (last change: 2026-02-01 | Add unit tests using testthat and clean up phihat to phi in npregiv*)
- `man/np.regressioniv.Rd`: npregiv documentation updates. (last change: 2026-02-02 | Update unit tests and fix documentation mismatches for npregiv/npregivderiv)
- `man/np.regressionivderiv.Rd`: npregivderiv documentation updates. (last change: 2026-02-02 | Reorder npregivderiv arguments alphabetically after fixed ones)
- `npver`: Version tracking file updated. (last change: 2026-02-04 | Bumped up version # to 0.60-21 in anticipation of moving to open issues)
- `src/Makevars.win`: Windows MPI build configuration added/updated. (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)
- `src/np.c`: C core updates (logic tweaks and stability fixes in computation kernels). (last change: 2026-02-05 | Sync non-MPI code fixes from npRmpi)

## Issue Mapping (Best-Effort)

The following issues appear to be addressed by changes in this tree:

- `#4` (npcdensbw cv.ls -> cv.ml segfault): addressed in C/R core updates (see `src/np.c` and related bandwidth logic).
- `#6` (npindexbw bwtype option): addressed in `R/np.singleindex.bw.R` and `R/np.singleindex.R`.
- `#7` (uocquantile subset of factor): addressed in `R/util.R`.
- `#13` (npudens categorical variance): addressed in `R/np.density.R`.
- `#18` (npregiv Tikhonov multi‑instrument): addressed in `R/npregiv.R`.
- `#50` (cross‑validated pairs plot): addressed by `R/np.pairs.R`.

If any of these are not fully resolved, treat this as a pointer to where the relevant changes landed.
