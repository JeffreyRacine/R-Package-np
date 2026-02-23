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
- [ ] Full `--as-cran` warning/note closure still required (accepted-warning ledger now tracked in `AS_CRAN_WARNING_LEDGER.md`).
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
