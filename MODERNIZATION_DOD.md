# np-master Modernization Definition of Done

## Goal
Ship a release-candidate-quality `np` that is modern, stable, performance-accountable, and aligned with:

- https://r-pkgs.org/
- https://adv-r.hadley.nz/
- https://r4ds.hadley.nz/

## Scope Priority
1. Core estimator families first: `npreg*`, `npudens*`, `npcdens*`, `npudist*`, `npcdist*`, `npscoef*`, `npindex*`, `npplreg*`.
2. `npregiv*` is explicitly lower priority unless blocking core release readiness.

## Gate Snapshot (2026-02-23)
- [x] R-layer `.C(` callsite retirement complete (`0`).
- [x] Active `<<-` retirement complete in R layer (`0`).
- [x] String-based `do.call("<string>", ...)` retirement complete in active R-layer paths.
- [x] High-risk `eval(parse(...))` absent; residual `eval(...)` centralized in shared helper with contracts.
- [x] Assignment-inside-`if` control flow retired in active R-layer codepaths (including `npregiv`).
- [x] Core modernization checkpoints validated with targeted contract tests + tarball checks.
- [x] Core bandwidth constructor scalar branches (`ifelse` -> scalar `if`) aligned in `dbandwidth/rbandwidth/conbandwidth/condbandwidth/smoothbandwidth/sibandwidth` plus `npcopula` and `gsl_bspline` helper guards (`51e9232`).
- [x] Core bw selector indexing is zero-length-safe (`seq_len`) in distribution/conditional/index bw paths (`0716cd6`).
- [x] Residual `1:n` index patterns retired in core estimator bw/index/smoothcoef/plreg paths (`37a157c`).
- [x] Conditional bw `goodrows` row-indexing now uses `seq_len(nrow(...))` in density/distribution selectors (`6d602bb`), with pre/post perf + parity artifacts recorded.
- [x] Conditional bw column/index reconstruction now uses `seq_len(...)`/safe ranges (`df809f8`) including hardened `gbw` split indexing for `xncon==0` edge paths.
- [x] Verified issue-note repro harness includes `npreg` factor-dispatch guard (`4802530`).
- [x] Native bridge stress harness added and passing for touched `.Call` surfaces (`issue_notes/native_bridge_stress.R`).
- [x] `--as-cran` reports no code/documentation mismatches (`/tmp/np_master_check_ascran_postloadhook_20260223.log`).
- [x] Local `--as-cran` warning closure achieved; only accepted CRAN incoming version-jump NOTE remains (`/tmp/np_master_check_ascran_compact_shellcheck_20260223.log`).
- [ ] Win-builder validation still required before release candidate.

## Mandatory Release Gates

### 1) Interface and Semantics
- [x] Public APIs for core estimators are stable (signature + return-structure contracts).
- [x] Formula/default method pairs have parity tests (including `subset`/`na.action` behavior).
- [x] S3 docs and method signatures match exactly (no codoc mismatches).

### 2) Evaluation and Call Construction
- [x] No `eval(parse(...))` in R layer.
- [x] High-risk `eval(...)` call construction paths are replaced with structured calls (`do.call`, explicit call objects) where appropriate.
- [x] Any remaining `eval(...)` has a documented reason and a regression test.

### 3) Native Interface Safety
- [x] `.C` callsites in R layer are `0`.
- [x] `.Call` interface paths have stress tests for touched entry points.
- [x] PROTECT/UNPROTECT accounting validated in modified C entry points.
- [x] No new stack-imbalance warnings in targeted stress runs.

### 4) Performance Governance
- [ ] Every performance patch includes pre/post comparison with identical script/args.
- [ ] Reports include mean and median percent deltas.
- [ ] Both fixed-seed and varying-seed runs are reported.
- [ ] At least one numerical parity check is included (objective and/or bandwidth).
- [ ] Artifacts are saved in `/tmp` with clear names.

### 5) Documentation and Examples
- [ ] Examples are runnable/minimal; heavy workflows in `\dontrun{}`.
- [ ] CRAN gate assumes `\dontrun{}` policy; full `\dontrun{}` execution is optional and used only for intentional behavior validation.
- [ ] Canonical references are linked from core docs (`np.kernels`, `np.options`, `np.plot`).
- [ ] Any user-facing mode changes have migration notes in NEWS/CHANGELOG.

### 6) Check and Packaging Hygiene
- [x] Tarball-first validation is cleanly run:
  - `R CMD build np-master`
  - `R CMD check --as-cran np_<ver>.tar.gz`
- [x] Issue-note regression sweep is run periodically and after core modernization touches:
  - `./issue_notes/run_verified_issue_repros.sh`
- [x] New warnings/notes are treated as regressions unless explicitly accepted and documented.
- [ ] Windows validation completed via win-builder prior to release candidate.

## Required Benchmark/Validation Record per Checkpoint
Include in commit body or companion note:

1. Workload definition (`DGP`, formula, `n`, `times`, seed policy, serial/parallel mode).
2. Mean and median performance deltas.
3. Numerical parity outcome and tolerance.
4. Known regime-specific caveats.
5. Output artifact paths (`/tmp/...`).
6. If shared `np`/`npRmpi` codepaths are touched, record separate smoke outcomes for each package; do not treat one as a proxy for the other.

## Current Residual Risks (Known)
- Some documentation/check warnings in non-core areas are still present and must be triaged as either:
  - accepted non-target technical debt with rationale, or
  - fixed before release candidate.
- Any change touching formula evaluation semantics remains medium/high risk and requires focused contract tests.

## Conditional BW Column-Index Safety Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced residual `1:n` column reconstruction and `setdiff(1:(...))` forms with `seq_len(...)` in conditional bw selectors.
2. Hardened `gbw` split initialization in `npcdistbw.condbandwidth` to avoid `1:0`/descending range hazards when `xncon==0`.
3. Scope:
   - `R/np.condensity.bw.R`
   - `R/np.condistribution.bw.R`
4. Commit:
   - `np-master`: `df809f8`
5. Validation:
   - parse gates for touched files: `PARSE_OK`
   - targeted tests:
     - `/tmp/np_master_condbw_seqcols_tests_20260223.log` (`TEST_RC:0`)
   - edge smoke (`xncon==0` path) on fresh install:
     - `/tmp/np_condbw_xncon0_smoke_installed_20260223.out` (`NP_CONDBW_XNCON0_SMOKE_OK`)
   - issue-note verified repro sweep:
     - `/tmp/np_issue_notes_repros_seqcols_20260223.log`
   - tarball check:
     - `/tmp/np_master_check_seqcols_20260223.log` (`Status: OK`)

## BW Guarded-Slice `seq_len` Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced guarded `1:gbw` fallback slices with explicit `seq_len(gbw)` indices in bw constructors.
2. Scope:
   - `R/np.condensity.bw.R`
   - `R/np.distribution.bw.R`
3. Commit:
   - `np-master`: `495c8a3`
4. Validation:
   - parse gates: `PARSE_OK`
   - targeted tests:
     - `/tmp/np_master_gbwidx_tests_20260223.log` (`TEST_RC:0`)
   - fresh-install surface smoke:
     - `/tmp/np_gbwidx_surface_smoke3_20260223.out` (`NP_GBWIDX_SURFACE_SMOKE_OK`)
   - issue-note verified repro sweep:
     - `/tmp/np_issue_notes_repros_gbwidx_20260223.log`
   - tarball check:
     - `/tmp/np_master_check_gbwidx_20260223.log` (`Status: OK`)

## DBandwidth `rorder` seq_len Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced fragile `rorder` reconstruction using `(1:ncol)[...]` with zero-length-safe `seq_len(ncol)` indexing.
2. Scope:
   - `R/np.distribution.bw.R`
3. Commit:
   - `np-master`: `4854b5f`
4. Validation:
   - parse gate: `PARSE_OK`
   - targeted tests:
     - `/tmp/np_master_dist_rorder_tests_20260223.log` (`TEST_RC:0`)
   - issue-note verified repro sweep:
     - `/tmp/np_issue_notes_repros_dist_rorder_20260223.log`
   - tarball check:
     - `/tmp/np_master_check_dist_rorder_20260223.log` (`Status: OK`)

## Density/Regression `rorder` seq_len Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced residual `(1:ncol)[...]` `rorder` reconstruction with `seq_len(ncol)` in core density/regression paths.
2. Scope:
   - `R/np.density.bw.R`
   - `R/np.regression.bw.R`
   - `R/np.regression.R`
3. Commit:
   - `np-master`: `a6dcd09`
4. Validation:
   - parse gates: `PARSE_OK`
   - targeted tests:
     - `/tmp/np_master_rorder_regdens_tests_20260223.log` (`TEST_RC:0`)
   - fresh-install smoke:
     - `/tmp/np_regdens_rorder_smoke_20260223.out` (`NP_REGDENS_RORDER_SMOKE_OK`)
   - issue-note verified repro sweep:
   - `/tmp/np_issue_notes_repros_rorder3_20260223.log`
   - tarball check:
     - `/tmp/np_master_check_rorder3_20260223.log` (`Status: OK`)

## Core Estimator `seq_len` Loop/Index Hardening Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced residual `1:n` loop/index forms in core estimator families with zero-length-safe `seq_len(...)`.
2. Scope:
   - `R/np.distribution.bw.R`
   - `R/np.condistribution.bw.R`
   - `R/np.singleindex.bw.R`
   - `R/np.smoothcoef.bw.R`
   - `R/np.smoothcoef.R`
   - `R/np.plregression.R`
3. Commit:
   - `np-master`: `a9b0f40`
4. Validation:
   - parse gates:
     - `/tmp/np_master_seq_len_core2_parse_20260223.log` (`RC:0`)
   - targeted contracts:
     - `/tmp/np_master_seq_len_core2_tests_postfix_20260223.log` (`PASS 29, FAIL 0`)
   - issue-note verified repro sweep:
     - `/tmp/np_issue_notes_repros_seqlen_core2_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_seqlen_core2_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_seqlen_core2_20260223.log` (`CHECK_RC:0`)

## Smoothcoef NDIM `seq_len` Hardening Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced residual `sapply(1:bws$ndim, ...)` with zero-length-safe `sapply(seq_len(bws$ndim), ...)` in core smooth coefficient bandwidth scaling setup.
2. Scope:
   - `R/np.smoothcoef.bw.R`
3. Commit:
   - `np-master`: `e9fd725`
4. Validation:
   - parse gate:
     - `/tmp/np_master_scoef_ndim_parse_20260223.log` (`RC:0`)
   - targeted contracts:
     - `/tmp/np_master_scoef_ndim_tests_20260223.log` (`PASS 19, FAIL 0`)
   - issue-note verified repro sweep:
     - `/tmp/np_issue_notes_repros_scoef_ndim_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_scoef_ndim_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_scoef_ndim_20260223.log` (`CHECK_RC:0`)

## Sigtest Index/Bootstrap Hygiene Checkpoint (2026-02-23)
Completed in `np-master`:
1. Fixed duplicate-index validation bug in `npsigtest` and modernized related bootstrap indexing patterns.
2. Scope:
   - `R/np.sigtest.R`
   - `tests/testthat/test-npsigtest.R`
3. Changes:
   - fixed repeated-index guard typo:
     - `if(length(unique(index)) < length(unique))` -> `if(length(unique(index)) < length(index))`
   - `index` default now uses `seq_len(ncol(xdat))`
   - `goodrows` initialization uses `seq_len(nrow(xdat))`
   - bootstrap loops use `seq_len(boot.num)`
   - pairwise bootstrap resampling uses `sample.int(num.obs, replace = TRUE)`
   - indicator means use logical means directly (`mean(In.vec > In)`).
4. Commit:
   - `np-master`: `6198c04`
5. Validation:
   - parse:
     - `/tmp/np_master_sigtest_parse_20260223.log` (`RC:0`)
   - focused tests:
     - `/tmp/np_master_sigtest_tests_20260223.log` (`PASS 5, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_sigtest_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_sigtest_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_sigtest_20260223.log` (`CHECK_RC:0`)

## Unitest Bootstrap Sampling Hygiene Checkpoint (2026-02-23)
Completed in `np-master`:
1. Modernized null-draw bootstrap sampling in `npunitest` to use `sample.int(...)` instead of index-vector sampling via `sample(1:length(...), ...)`.
2. Scope:
   - `R/np.unitest.R`
3. Commit:
   - `np-master`: `d2ecf8c`
4. Validation:
   - parse:
     - `/tmp/np_master_unitest_parse_20260223.log` (`RC:0`)
   - focused tests:
     - `/tmp/np_master_unitest_tests_20260223.log` (`PASS 9, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_unitest_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_unitest_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_unitest_20260223.log` (`CHECK_RC:0`)
5. Note:
   - `np-master` check warning count in this run reflects pre-existing non-target tree state (top-level/codoc drift), not this `np.unitest.R` change.

## Unitest Loop/P-Value Scalar Hygiene Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced residual scalar bootstrap loop/index and indicator-mean pattern in `npunitest`.
2. Scope:
   - `R/np.unitest.R`
3. Changes:
   - `for(b in 1:boot.num)` -> `for (b in seq_len(boot.num))`
   - `mean(ifelse(resampled.stat > test.stat, 1, 0))` -> `mean(resampled.stat > test.stat)`
4. Commit:
   - `np-master`: `dfb5ff6`
5. Validation:
   - parse:
     - `/tmp/np_master_unitest2_parse_20260223.log` (`RC:0`)
   - focused tests:
     - `/tmp/np_master_unitest2_tests_20260223.log` (`PASS 9, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_unitest2_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_unitest2_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_unitest2_20260223.log` (`CHECK_RC:0`)

## Np*test Bootstrap Loop/Index Hygiene Checkpoint (2026-02-23)
Completed in `np-master`:
1. Modernized bootstrap loop/index idioms in `np.*test` paths and fixed `npsdeptest` bootstrap lag routing semantics for `lag.num > 1`.
2. Scope:
   - `R/np.cmstest.R`
   - `R/np.qcmstest.R`
   - `R/np.deneqtest.R`
   - `R/np.deptest.R`
   - `R/np.sdeptest.R`
3. Changes:
   - `1:boot.num`/`1:lag.num` loops -> `seq_len(...)`
   - indicator p-value means use logical means directly
   - `npsdeptest` bootstrap now computes lagged pairs inside the per-lag loop from a single resampled series (`resampled.ts`), removing stale-`k` behavior.
4. Commit:
   - `np-master`: `ee2b478`
5. Validation:
   - parse gate: `PARSE_OK`
   - focused tests:
     - `/tmp/np_master_nptests_fix_20260223.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_nptestsfix_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_nptestsfix_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_nptestsfix_20260223.log` (`Status: 1 WARNING` pre-fix codoc)

## Sigtest Rd Signature Alignment Checkpoint (2026-02-23)
Completed in `np-master`:
1. Aligned `npsigtest.rbandwidth` Rd default for `index` with code (`seq_len(ncol(xdat))`).
2. Scope:
   - `man/np.sigtest.Rd`
3. Commit:
   - `np-master`: `da2e4a8`
4. Validation:
   - tarball-first:
     - `/tmp/np_master_check_sigdoc_20260223.log` (`Status: OK`)
   - `--as-cran` refresh:
     - `/tmp/np_master_check_ascran_refresh_20260223.log` (`Status: 2 WARNINGs, 1 NOTE`; no codoc warning).

## Core Estimator Performance Governance Checkpoint (2026-02-23)
Completed in `np-master` (serial mode):
1. Pre/post commit comparison on core estimators with identical scripts/args:
   - pre: `5298a49`
   - post: `4158fbb`
2. Workload and policy:
   - script: `/tmp/bench_core_np_lib_20260223.R`
   - estimators: `npreg`, `npudens`, `npcdens`
   - `n=250`, `times=5`
   - seed policy: fixed (`42`) and varying (`42 + i - 1`)
3. Invocation:
   - `Rscript /tmp/bench_core_np_lib_20260223.R <lib> <csv> <metrics.rds> 5 250`
   - compare: `Rscript /tmp/bench_compare_prepost_20260223.R <pre_csv> <post_csv> <pre_rds> <post_rds> <summary_csv> <parity_txt>`
4. Mean/median percent change summary:
   - fixed:
     - `npreg`: mean `-1.31%`, median `0.00%`
     - `npudens`: mean `-1.19%`, median `-1.23%`
     - `npcdens`: mean `+1.23%`, median `+1.05%`
   - varying:
     - `npreg`: mean `-3.54%`, median `-3.92%`
     - `npudens`: mean `-0.76%`, median `-2.82%`
     - `npcdens`: mean `-0.53%`, median `-1.02%`
5. Parity:
   - fixed-seed metric parity max abs diff: `0`
   - note: `npcdens_bw_mean` is non-finite (`NA`) in both runs; fitted-value parity remained exact.
6. Artifacts:
   - pre/post raw:
     - `/tmp/bench_np_pre_full_20260223.csv`
     - `/tmp/bench_np_post_full_20260223.csv`
   - summary/parity:
     - `/tmp/bench_np_full_summary_20260223.csv`
     - `/tmp/bench_np_full_parity_20260223.txt`

## Win-Builder Submission Checkpoint (2026-02-23)
Submitted from this environment:
1. `np-master` release check submitted via:
   - `Rscript -e \"suppressPackageStartupMessages(library(devtools)); check_win_release('/Users/jracine/Development/np-master')\"`
   - log: `/tmp/winbuilder_submit_np_20260223.log`
2. `np-npRmpi` companion submission also launched (for synchronized modernization release readiness):
   - `Rscript -e \"suppressPackageStartupMessages(library(devtools)); check_win_release('/Users/jracine/Development/np-npRmpi')\"`
   - log: `/tmp/winbuilder_submit_nprmpi_20260223.log`
3. Status:
   - submissions accepted by win-builder; final pass/fail results pending email report (~15-30 minutes post submission).

## Plot Engine Scalar-Branch Hygiene Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced high-confidence scalar `ifelse(...)` branches in plot engines with scalar `if` expressions.
2. Scope:
   - `R/np.plot.engine.condbandwidth.R`
   - `R/np.plot.engine.conbandwidth.R`
   - `R/np.plot.engine.scbandwidth.R`
   - `R/np.plot.engine.plbandwidth.R`
3. Changes:
   - `dsf = ifelse(gradients, bws$xndim, 1)` -> `dsf = if (gradients) bws$xndim else 1`
   - `i = ifelse(plot.index <= bws$xndim, ...)` -> scalar `if` index branch.
4. Commit:
   - `np-master`: `2e57883`
5. Validation:
   - parse gate: `PARSE_OK`
   - focused tests:
     - `/tmp/np_master_plot_scalar_ifelse_tests_20260223.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_plotscalar_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_plotscalar_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_plotscalar_20260223.log` (`Status: OK`)

## Plot Engine Scalar-`ifelse` Expansion Checkpoint (2026-02-23)
Completed in `np-master`:
1. Retired additional scalar `ifelse(...)` usage in plot engines where control expressions are scalar (`xi.factor`, `gradients`, `quantreg`, `cdf`).
2. Scope:
   - `R/np.plot.engine.condbandwidth.R`
   - `R/np.plot.engine.conbandwidth.R`
   - `R/np.plot.engine.scbandwidth.R`
   - `R/np.plot.engine.plbandwidth.R`
3. Changes:
   - scalar error-bar style/shape selectors now use scalar `if`
   - scalar gradient/quantile error-name/bias-name selectors now use scalar `if` chains
   - scalar label composition branches (`ylab`, `zlab`) now use scalar `if`
4. Commit:
   - `np-master`: `0871f2b`
5. Validation:
   - parse gate: `PARSE_OK`
   - focused tests:
     - `/tmp/np_master_plot_ifelse2_tests_20260223.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_plotifelse2_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_plotifelse2_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_plotifelse2_20260223.log` (`Status: OK`)

## Plot Bandwidth/Dbandwidth Scalar-Branch Checkpoint (2026-02-23)
Completed in `np-master`:
1. Replaced residual scalar `ifelse(...)` branches in density/regression bandwidth plot engines.
2. Scope:
   - `R/np.plot.engine.bandwidth.R`
   - `R/np.plot.engine.dbandwidth.R`
3. Changes:
   - scalar `xi.factor` error-style/error-bar selectors now use scalar `if`
   - scalar `xi.factor` linetype selector in `dbandwidth` now uses scalar `if`
4. Commit:
   - `np-master`: `d9d8242`
5. Validation:
   - parse gate: `PARSE_OK`
   - focused tests:
     - `/tmp/np_master_plot_ifelse3_tests_20260223.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_plotifelse3_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_plotifelse3_20260223.log` (`BUILD_RC:0`)
     - `/tmp/np_master_check_plotifelse3_20260223.log` (`Status: OK`)

## Plot `rbandwidth` Scalar-Branch + Vector-Safe Error-Bar Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced remaining scalar `ifelse(...)` branches in `np.plot.engine.rbandwidth.R` with scalar `if` for:
   - gradient label prefix,
   - factor plot error style/bar selection,
   - factor/continuous `lty` selection.
2. Modernized scalar helper branches in `np.plot.helpers.R`:
   - `gen.label(...)`,
   - `gen.tflabel(...)`.
3. Guarded vector semantics in `draw.error.bars(...)`:
   - initial scalarization of `ifelse(htest, ...)` was invalid because `htest` is vector-valued,
   - finalized to vector-safe equivalent `hdelta = pmin(yg, hbardist)/2`.
4. Commit:
   - `np-master`: `dfe83c8`
5. Validation:
   - targeted tests:
     - `/tmp/np_master_plot_ifelse4_tests2_20260224.log` (`PASS 26, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_plotifelse4_20260224.log` (`RC:0`)
   - tarball-first:
     - `/tmp/np_master_build_plotifelse4_20260224.log` (`BUILD_RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_plotifelse4_20260224.log` (`Status: OK`)

## Indicator/Scalar `ifelse` Cleanup Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced remaining scalar/indicator `ifelse(...)` in active runtime paths:
   - `R/np.plot.engine.sibandwidth.R`: scalar gradient-name prefix branch.
   - `R/np.cmstest.R`, `R/np.qcmstest.R`, `R/np.symtest.R`: p-value indicator means now use direct logical means.
2. Changes:
   - `paste(ifelse(gradients, ...))` -> scalar `if` branch.
   - `mean(ifelse(condition, 1, 0))` -> `mean(condition)`.
3. Commit:
   - `np-master`: `f8f1ff4`
4. Validation:
   - targeted tests:
     - `/tmp/np_master_nptests_plotifelse5_20260224.log` (`PASS 26, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_ifelse5_20260224.log` (`RC:0`)
   - tarball-first:
     - `/tmp/np_master_build_ifelse5_20260224.log` (`BUILD_RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ifelse5_20260224.log` (`Status: OK`)

## Print/Reporting Scalar-Branch Cleanup Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced scalar `ifelse(...)` in print/reporting paths with explicit scalar `if` branches.
2. Scope:
   - `R/cmstest.R`
   - `R/deneqtest.R`
   - `R/deptest.R`
   - `R/sdeptest.R`
   - `R/symtest.R`
   - `R/unitest.R`
3. Commit:
   - `np-master`: `9ea1c59`
4. Validation:
   - parse gates: `PARSE_OK`
   - targeted tests:
     - `/tmp/np_master_print_ifelse6_tests_20260224.log` (`PASS 13, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_printifelse6_20260224.log` (`RC:0`)
   - tarball-first:
     - `/tmp/np_master_build_printifelse6_20260224.log` (`BUILD_RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_printifelse6_20260224.log` (`Status: OK`)

## `b.star` Elementwise Rounding/Bounds Fix Checkpoint (2026-02-24)
Completed in `np-master`:
1. Fixed a multivariate rounding/bounds bug in `b.star(..., round=TRUE)` where `BstarCB` used scalar `max(1, round(BstarCB))` inside vectorized `ifelse`, which could collapse column-wise outputs.
2. Modernized bounds logic with elementwise vector ops:
   - `M <- min(2 * mhat, mmax)`
   - `pmin(...)` / `pmax(...)` for SB/CB bounds.
3. Added regression coverage for multivariate `round=TRUE` elementwise behavior.
4. Scope:
   - `R/b.star.R`
   - `tests/testthat/test-utils.R`
5. Commit:
   - `np-master`: `c2916ae`
6. Validation:
   - focused tests:
     - `/tmp/np_master_bstar_fix_tests_20260224.log` (`PASS 6, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_bstarfix_20260224.log` (`RC:0`)
   - tarball-first:
     - `/tmp/np_master_build_bstarfix_20260224.log` (`BUILD_RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_bstarfix_20260224.log` (`Status: OK`)

## Issue-Note Hat/Bootstrap Regression Guard Checkpoint (2026-02-24)
Completed in `np-master`:
1. Extended periodic verified-issue repro harness with explicit serial guards for modern hat/bootstrap paths:
   - `npreghat` fit-operator parity versus `npreg` fitted values.
   - `plot(..., plot.errors.boot.method="wild-hat")` smoke on `rbandwidth`.
2. Scope:
   - `issue_notes/verified_issue_repros.R`
3. Commit:
   - `np-master`: `80191c2`
4. Validation:
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_032316.log` (all verified repros passed)

## Quantile/Copula Loop-Index Hygiene Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced residual `1:length(...)` / `1:n...` loop patterns in auxiliary quantile/copula runtime paths with `seq_along(...)` / `seq_len(...)`.
2. Scope:
   - `R/np.quantile.R`
   - `R/np.copula.R`
3. Commit:
   - `np-master`: `b99c0f2`
4. Validation:
   - targeted tests:
     - `/tmp/np_master_copula_quantile_tests_20260224.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_032635.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_loophyg_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_loophyg_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## Utility Loop-Index Hygiene Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced residual `1:length(...)` / `1:ncol(...)` index patterns in shared utility/runtime helpers with `seq_along(...)` / `seq_len(...)`.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-master`: `fd4e13e`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='utils|npcopula|npudist|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_033540.log` (all verified repros passed)

## Univariate Density Loop-Index Hygiene Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced residual `1:length(...)` / `1:n...` loop/index patterns in univariate boundary/shape-constrained density helpers with `seq_along(...)` / `seq_len(...)`.
2. Scope:
   - `R/npuniden.boundary.R`
   - `R/npuniden.sc.R`
3. Commit:
   - `np-master`: `f2ec364`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='npuniden|utils', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_033821.log` (all verified repros passed)

## Plot-Helper Indexing Hygiene Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced residual `1:length(...)` / `1:ncol(...)` indexing in plot helper paths with `seq_along(...)` / `seq_len(...)`.
2. Replaced manual cumulative-loop `sapply(1:length(tq), ...)` with `cumsum(tq)` in quantile helper.
3. Scope:
   - `R/np.plot.helpers.R`
4. Commit:
   - `np-master`: `3ddccc7`
5. Validation:
   - targeted tests:
     - `testthat::test_local(filter='semihat|npreghat|plot|npindex|npplreg|npscoef', reporter='summary')` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_034150.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_posthyg3_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_posthyg3_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## Utility/Test-Stat Loop Finalization Checkpoint (2026-02-24)
Completed in `np-master`:
1. Finalized residual scalar range loops in utility/test-stat helpers:
   - `cast(...)` dataframe-column loop,
   - raw polynomial derivative loop in `mypoly(...)`,
   - `b.star(...)` column loop,
   - `sdeptest` bootstrap/report loops.
2. Scope:
   - `R/util.R`
   - `R/b.star.R`
   - `R/sdeptest.R`
3. Commit:
   - `np-master`: `55aafa4`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='utils|nptests|sdeptest|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_035016.log` (all verified repros passed)

## Plot-Helper Residual Index Cleanup Checkpoint (2026-02-24)
Completed in `np-master`:
1. Removed remaining `1:length(...)` index helpers in low-level plot helper vectors (`jj` and factor plotting separators) using `seq_along(...)`.
2. Scope:
   - `R/np.plot.helpers.R`
3. Commit:
   - `np-master`: `df51a5c`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='plot|semihat|npreghat', reporter='summary')` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_035414.log` (all verified repros passed)

## Consolidated Tarball Gate (2026-02-24)
1. `R CMD build /Users/jracine/Development/np-master`
   - `/tmp/np_master_build_postmodern_final_20260224.log` (`RC:0`, `creating vignettes ... OK`)
2. `R CMD check --as-cran np_0.70-0.tar.gz`
   - `/tmp/np_master_check_ascran_postmodern_final_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## Plot-Engine Loop-Range Hygiene Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced residual `1:...` loop/index ranges in plot-engine runtime paths with `seq_len(...)`-safe forms.
2. Scope:
   - `R/np.plot.engine.bandwidth.R`
   - `R/np.plot.engine.dbandwidth.R`
   - `R/np.plot.engine.conbandwidth.R`
   - `R/np.plot.engine.condbandwidth.R`
3. Commit:
   - `np-master`: `a0bf55b`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='plot|semihat|npreghat', reporter='summary')` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_041140.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_plotengine_seq_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_plotengine_seq_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## `npregiv` Loop-Range Hygiene Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced residual `1:ncol(...)` / `1:length(...)` index ranges in `npregiv` support paths with `seq_len(...)` / `seq_along(...)` to avoid `1:0` edge behavior.
2. Scope:
   - `R/npregiv.R`
3. Commit:
   - `np-master`: `a1a56ab`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='npregiv', reporter='summary')` (`RC:0`; expected pre-existing monotone-stopping warnings)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_042226.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_npregiv_seq_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_npregiv_seq_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## `npregiv` Loop-Range Hygiene Follow-On Checkpoint (2026-02-24)
Completed in `np-master`:
1. Completed residual `1:n...` conversions in `npregiv` iterative/LOO loops and `sapply` index ranges.
2. Scope:
   - `R/npregiv.R`
3. Commit:
   - `np-master`: `de162bf`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='npregiv', reporter='summary')` (`RC:0`; expected pre-existing monotone-stopping warnings)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_043014.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_npregiv_seq2_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_npregiv_seq2_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 1 NOTE`; warning set unchanged from existing top-level/vignette-size debt)

## Bandwidth Metadata / `npregivderiv` Loop-Range Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced residual metadata/index range builders with `seq_len(...)` in bandwidth metadata constructors and `npregivderiv` numeric-column detection.
2. Scope:
   - `R/bandwidth.R`
   - `R/rbandwidth.R`
   - `R/smoothbandwidth.R`
   - `R/npregivderiv.R`
3. Commit:
   - `np-master`: `e5fb784`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='npregiv|bandwidth|bw-dispatch|formula-bw-contract', reporter='summary')` (`RC:0`; expected pre-existing monotone-stopping warnings)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_043811.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_bwmeta_seq_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_bwmeta_seq_20260224_withcrs.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)
5. Environment note:
   - local `--as-cran` checks required `crs` in `R_LIBS`; command used:
     - `R_LIBS=/tmp/crs_check_lib R CMD check --as-cran --no-manual np_0.70-0.tar.gz`

## Bandwidth Metadata Parity Follow-On Checkpoint (2026-02-24)
Completed in `np-master`:
1. Applied matching `seq_len(...)` metadata-index conversion in `R/dbandwidth.R` for cross-file parity.
2. Scope:
   - `R/dbandwidth.R`
3. Commit:
   - `np-master`: `1ff6a07`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='bw-dispatch|formula-bw-contract|bandwidth', reporter='summary')` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_044753.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_bwmeta_seq2_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_bwmeta_seq2_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## `b.star` Index-Range Guard Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced the `1:(mmax-Kn+1)` insignificant-run index builder with bounded `seq_len(max(mmax - Kn + 1L, 0L))` to avoid `1:0`/descending edge behavior while preserving normal positive-range semantics.
2. Scope:
   - `R/b.star.R`
3. Commit:
   - `np-master`: `9e95e51`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='utils|nptests|sdeptest|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_045406.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_bstar_seq_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_bstar_seq_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## `npuniden.boundary` Scalar-Branch Hygiene Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced scalar `ifelse(...)` branches in `npuniden.boundary` CV/optimizer routing with scalar `if` expressions:
   - finite fallback in LS-CV objective,
   - optimizer start-index choice (`cv.ml` max vs `cv.ls` min),
   - optimizer upper-bound choice for `beta2` kernels.
2. Scope:
   - `R/npuniden.boundary.R`
3. Commit:
   - `np-master`: `f1d3f7e`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/npuniden.boundary.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='npuniden|nptests|sdeptest')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
     - log: `/tmp/np_master_npuniden_scalar_ifelse_tests_20260224.log`
   - direct runtime smoke:
     - `/tmp/np_master_npuniden_scalar_ifelse_smoke_20260224.out` (`NPUNIDEN_BOUNDARY_SMOKE_OK`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_050442.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_npuniden_scalar_ifelse_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_npuniden_scalar_ifelse_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## Wild-Bootstrap Draw Allocation Checkpoint (2026-02-24)
Completed in `np-master`:
1. Reduced allocations in wild-bootstrap draw helpers by replacing matrix `ifelse(...)` draw generation with preallocated matrices and logical indexing.
2. Scope:
   - `R/np.plot.helpers.R`
3. Commit:
   - `np-master`: `da2eee9`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/np.plot.helpers.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='plot|semihat|npreghat|wild', reporter='summary')` (`RC:0`)
     - log: `/tmp/np_master_wilddraw_alloc_tests_20260224.log`
   - direct helper contract smoke:
     - `/tmp/np_master_wilddraw_contract_20260224b.out` (`NP_WILDDRAW_CONTRACT_OK`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_051458.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_wilddraw_alloc_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_wilddraw_alloc_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## `dlev` Sequence-Safety Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced `as.numeric(1:nlevels(x))` with `as.numeric(seq_len(nlevels(x)))` in `dlev()` to eliminate `1:0` edge behavior while preserving factor-level mapping.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-master`: `219d062`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='utils|nptests|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
     - log: `/tmp/np_master_dlev_seq_tests_20260224.log`
   - direct smoke:
     - `/tmp/np_dlev_seq_smoke_20260224.out` (`DLEV_SEQ_OK`)
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_052223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_seq_edge_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_seq_edge_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## Bootstrap Wild-Draw Refactor Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced repeated `ifelse(rbinom(...))` bootstrap multipliers with a lightweight preallocated vector helper (`runif` + logical assignment) in test/bootstrap-heavy paths.
2. Scope:
   - `R/np.cmstest.R`
   - `R/np.qcmstest.R`
   - `R/np.sigtest.R`
3. Commit:
   - `np-master`: `8361490`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/np.cmstest.R')); invisible(parse(file='R/np.qcmstest.R')); invisible(parse(file='R/np.sigtest.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='cmstest|qcmstest|sigtest', reporter='summary')` (`RC:0`)
     - log: `/tmp/np_master_bootdraw_refactor_tests_20260224.log`
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_053149.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_bootdraw_refactor_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_bootdraw_refactor_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 1 NOTE`; warning set unchanged from existing top-level/vignette-size debt)

## `npqcmstest` Quantile-Residual Helper Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced repeated `ifelse(model.resid <= 0, 1 - tau, -tau)` logic with a local helper preserving NA semantics, reused across bandwidth selection and test-statistic internals.
2. Scope:
   - `R/np.qcmstest.R`
3. Commit:
   - `np-master`: `80c6fd1`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/np.qcmstest.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='qcmstest|cmstest|sigtest', reporter='summary')` (`RC:0`)
     - log: `/tmp/np_master_qresidual_helper_tests_20260224.log`
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_053830.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_qresidual_helper_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_qresidual_helper_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## `npuniden.boundary` + `npregiv` Edge-Safety Follow-On (2026-02-24)
Completed in `np-master`:
1. Hardened remaining low-risk edge paths:
   - replaced vector `ifelse(...)` finite clamp in CV-ML objective with explicit mask assignment (`f.safe`) to reduce branch allocation and keep invalid values at `.Machine$double.xmin`,
   - replaced debug-only `seq(1:num.bw)` with `seq_len(num.bw)` in `npregiv`.
2. Scope:
   - `R/npuniden.boundary.R`
   - `R/npregiv.R`
3. Commit:
   - `np-master`: `759d749`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/npuniden.boundary.R')); invisible(parse(file='R/npregiv.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='npuniden|npregiv', reporter='summary')` (`RC:0`; expected pre-existing `npregiv` and `npuniden.sc` warnings)
     - log: `/tmp/np_master_uniden_npregiv_edge_tests_20260224.log`
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_054656.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_uniden_npregiv_edge_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_uniden_npregiv_edge_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## `genGofStr` Scalar-Branch Hygiene Checkpoint (2026-02-24)
Completed in `np-master`:
1. Replaced nested scalar `ifelse(...)` in `genGofStr()` with explicit scalar `if` branches for `MSE`/`R2` report fragments.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-master`: `4e48d75`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='utils|nptests|sdeptest|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
     - log: `/tmp/np_master_gofstr_scalar_tests_20260224.log`
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_055424.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_gofstr_scalar_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_gofstr_scalar_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)

## `util.R` Scalar-Branch Follow-On (`NZD`/`genBwScaleStrs`) (2026-02-24)
Completed in `np-master`:
1. Removed additional scalar `ifelse(...)` usage in hot/summary utility paths:
   - `NZD()` small-value assignment now uses direct sign-mask replacement,
   - `genBwScaleStrs()` scalar label/summary string branches now use explicit scalar `if`.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-master`: `b795e3b`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='utils|bandwidth|bw-dispatch|formula-bw-contract', reporter='summary')` (`RC:0`)
     - log: `/tmp/np_master_util_scalar2_tests_20260224.log`
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_060257.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_util_scalar2_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_util_scalar2_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 1 NOTE`; warning set unchanged from existing top-level/vignette-size debt)

## `util.R` Vector-Branch Follow-On (`tgauss`/`pad`/`rpad`) (2026-02-24)
Completed in `np-master`:
1. Replaced additional `ifelse(...)` usages with explicit vector-safe masked assignments in utility hot paths:
   - `nptgauss()` local `tgauss` now uses precomputed vector + mask assignment,
   - `pad()` and `rpad()` now use vector-safe index assignment while preserving names.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-master`: `f08e8d4`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests:
     - `testthat::test_local(filter='utils|nptests|sdeptest|npuniden|bandwidth', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
     - log: `/tmp/np_master_util_scalar3_tests_20260224.log`
   - issue-note repro sweep:
     - `/tmp/np_issue_notes_repros_20260224_061008.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/np_master_build_util_scalar3_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/np_master_check_ascran_util_scalar3_20260224.log` (`RC:0`, `Status: 2 WARNINGs, 2 NOTEs`; warning set unchanged from existing top-level/vignette-size debt)
