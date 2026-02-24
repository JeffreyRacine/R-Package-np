# np-npRmpi Modernization Definition of Done

## Goal
Ship a release-candidate-quality `npRmpi` that is modern, robust in MPI lifecycle behavior, and aligned with:

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
- [x] Session-mode regression (`npRmpi.init(nslaves=1)` startup hang path) fixed in spawn-report route (`7b5ba36`).
- [x] Session API scalar controls validated early (`npRmpi.init/quit`) to fail before MPI internals on invalid inputs (`3cf92e7`).
- [x] Subprocess timeout routing contracts added for session smoke and skip-init fast-fail (`e6581f3`).
- [x] Subprocess timeout routing contract added for manual-broadcast mode smoke (`0db36a8`).
- [x] Opt-in `mpiexec` attach-mode smoke contract added (`c719107`, enabled by `NP_RMPI_ENABLE_ATTACH_TEST=1`).
- [x] Core modernization checkpoints validated with targeted contract tests + tarball checks in working MPI env.
- [x] Core bandwidth constructor scalar branches (`ifelse` -> scalar `if`) aligned in `dbandwidth/rbandwidth/conbandwidth/condbandwidth/smoothbandwidth/sibandwidth` plus `npcopula` and `gsl_bspline` helper guards (`91ec5dc`).
- [x] Core bw selector indexing is zero-length-safe (`seq_len`) in distribution/conditional/index bw paths (`57222c8`).
- [x] Residual `1:n` index patterns retired in core estimator bw/index/smoothcoef/plreg paths (`531216c`).
- [x] Conditional bw `goodrows` row-indexing now uses `seq_len(nrow(...))` in density/distribution selectors (`29aed7c`), with pre/post perf + parity artifacts recorded.
- [x] Conditional bw column/index reconstruction now uses `seq_len(...)`/safe ranges (`45dc5cc`) including hardened `gbw` split indexing for `xncon==0` edge paths.
- [x] Verified issue-note repro harness includes session `npreg` factor-routing guard (`b5b597a`).
- [x] Load hook now supports source-tree/dev loading without installed-package lookup dependence (`61b4783`).
- [x] Native bridge stress harness added and passing for touched `.Call` surfaces in session mode (`issue_notes/native_bridge_stress.R`).
- [x] `--as-cran` reports no code/documentation mismatches (`/tmp/nprmpi_check_ascran_postloadhook_20260223.log`).
- [x] Local `--as-cran` warning closure achieved; only accepted CRAN incoming version-jump NOTE remains (`/tmp/nprmpi_check_ascran_compact_shellcheck_20260223.log`).
- [ ] Win-builder validation still required before release candidate.

## Mandatory Release Gates

### 1) Interface and Semantics
- [x] Public APIs for core estimators are stable (signature + return-structure contracts).
- [x] Formula/default method pairs have parity tests (including `subset`/`na.action` behavior).
- [x] S3 docs and method signatures match exactly (no codoc mismatches).

### 2) MPI Execution Model Integrity
- [x] Attach mode (`mpiexec` + `.Rprofile`) works for core workflows.
- [x] Session mode (`npRmpi.init/quit`) works for core workflows.
- [x] Manual broadcast high-performance path remains functional.
- [ ] Any new MPI helper path has explicit cleanup guarantees (`npRmpi.stop(force=TRUE)` or `mpi.quit()` route).
- [x] Smoke protocol explicitly requires `npRmpi.init(...)` before estimator calls and guaranteed `npRmpi.quit(...)` cleanup (never run `npRmpi` smoke as plain `np`).

### 3) Evaluation and Call Construction
- [x] No `eval(parse(...))` in R layer.
- [x] High-risk `eval(...)` paths replaced where appropriate with structured calls.
- [x] Any remaining `eval(...)` is documented and covered by tests.

### 4) Native Interface Safety
- [x] `.C` callsites in R layer are `0`.
- [x] `.Call` interface paths have stress tests for touched entry points.
- [x] PROTECT/UNPROTECT accounting validated for modified C entry points.
- [x] No new stack-imbalance warnings in targeted serial + MPI runs.

### 5) Performance Governance
- [ ] Every performance patch includes pre/post comparison with identical script/args.
- [ ] Reports include mean and median percent deltas.
- [ ] Both fixed-seed and varying-seed runs are reported.
- [ ] At least one numerical parity check is included.
- [ ] Artifacts are saved in `/tmp` with clear names.
- [ ] Performance notes distinguish execution model effects (attach vs session vs manual broadcast).

### 6) Documentation and Examples
- [ ] `npRmpi` docs clearly show three usage modes:
  - user-friendly session mode,
  - user-friendly `mpiexec` attach mode,
  - high-performance manual broadcast mode.
- [ ] Demos include working minimal scripts for serial, attach, and profile/manual broadcast paths.
- [ ] Examples are runnable/minimal; heavy workflows in `\dontrun{}`.
- [ ] CRAN gate assumes `\dontrun{}` policy; full `\dontrun{}` execution is optional and used only for intentional behavior validation.

### 7) Check and Packaging Hygiene
- [x] Tarball-first validation is cleanly run:
  - `R CMD build np-npRmpi`
  - `R CMD check --as-cran npRmpi_<ver>.tar.gz`
- [x] Issue-note regression sweep is run periodically and after core modernization touches:
  - `./issue_notes/run_verified_issue_repros.sh`
- [x] New warnings/notes are treated as regressions unless explicitly accepted and documented.
- [ ] Windows validation completed via win-builder prior to release candidate.

## Required Benchmark/Validation Record per Checkpoint
Include in commit body or companion note:

1. Workload definition (`DGP`, formula, `n`, `times`, seed policy, mode).
2. Mean and median performance deltas.
3. Numerical parity outcome and tolerance.
4. MPI lifecycle notes (spawn/attach cleanup behavior).
5. Output artifact paths (`/tmp/...`).
6. Timeout/cleanup notes for aborted runs (including residual `slavedaemon.R` / `Rslaves.sh` handling).

## Current Residual Risks (Known)
- Attach/session mode behavior may differ by environment and MPI interface settings (`FI_TCP_IFACE`), so smoke gates must run in both modes before release.
- Remaining non-target doc warnings and duplicated alias/cross-reference warnings should be triaged and either fixed or explicitly accepted.

## Conditional BW Column-Index Safety Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Replaced residual `1:n` column reconstruction and `setdiff(1:(...))` forms with `seq_len(...)` in conditional bw selectors.
2. Hardened `gbw` split initialization in `npcdistbw.condbandwidth` to avoid `1:0`/descending range hazards when `xncon==0`.
3. Scope:
   - `R/np.condensity.bw.R`
   - `R/np.condistribution.bw.R`
4. Commit:
   - `np-npRmpi`: `45dc5cc`
5. Validation:
   - parse gates for touched files: `PARSE_OK`
   - targeted tests:
     - `/tmp/nprmpi_condbw_seqcols_targeted_tests_20260223.log` (`TEST_RC:0`)
   - edge smoke (`xncon==0` path) on fresh install:
     - `/tmp/nprmpi_condbw_xncon0_smoke_installed_20260223.out` (`NPRMPI_CONDBW_XNCON0_SMOKE_OK`)
   - issue-note verified repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_seqcols_en0_20260223.log`
   - session-mode regression smoke (`npRmpi.init` + `npcdens`):
     - `/tmp/nprmpi_user_npcdens_smoke_20260223.out` (`NPRMPI_USER_NPCDENS_SMOKE_OK`)
   - tarball check (MPI env pinned):
     - `/tmp/nprmpi_check_seqcols_20260223.log` (`Status: OK`)

## BW Guarded-Slice `seq_len` Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Replaced guarded `1:gbw` fallback slices with explicit `seq_len(gbw)` indices in bw constructors.
2. Scope:
   - `R/np.condensity.bw.R`
   - `R/np.distribution.bw.R`
3. Commit:
   - `np-npRmpi`: `6ed24a6`
4. Validation:
   - parse gates: `PARSE_OK`
   - targeted tests:
     - `/tmp/nprmpi_gbwidx_tests_20260223.log` (`includes pre-existing npudist test failure outside touched constructors`)
   - focused fresh-install surface smoke:
     - `/tmp/nprmpi_gbwidx_surface_smoke3_20260223.out` (`NPRMPI_GBWIDX_SURFACE_SMOKE_OK`)
   - issue-note verified repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_gbwidx_20260223.log`
   - tarball check (MPI env pinned):
     - `/tmp/nprmpi_check_gbwidx_20260223.log` (`Status: OK`)

## Autodispatch Tmp-Call Replacement Fix (2026-02-23)
Completed in `np-npRmpi`:
1. Fixed large-argument autodispatch return rewriting so bandwidth `call` objects no longer retain unresolved temporary symbols (e.g. `.__npRmpi_autod_dat_1`).
2. Root cause:
   - return-path replacement only used `prepared$tmpvals`; large arguments are staged in `prepared$prepublish`.
3. Change:
   - in `.npRmpi_distributed_call_impl`, use combined replacement map `c(prepared$tmpvals, prepared$prepublish)` for list/call return rewriting.
4. Scope:
   - `R/np.autodispatch.R`
5. Commit:
   - `np-npRmpi`: `d6f5b6e`
6. Validation:
   - targeted unit file:
     - `/tmp/nprmpi_test_npudist_20260223.log` (`TEST_RC:0`)
   - direct repro (fresh install):
     - `/tmp/nprmpi_npudist_call_repro_20260223.out` (`NPRMPI_NPUDIST_CALL_REPRO_OK`)
   - broader targeted filter:
     - `/tmp/nprmpi_gbwidx_tests_postfix_20260223.log` (`TEST_RC:0`)
   - issue-note verified repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_autodtmpfix_20260223.log`
   - session-mode npcdens smoke:
     - `/tmp/nprmpi_user_npcdens_smoke_fix_20260223.out` (`NPRMPI_USER_NPCDENS_SMOKE_OK`)
   - tarball check:
     - `/tmp/nprmpi_check_autodtmpfix_20260223.log` (`Status: OK`)

## Routing Regression Contract Lock (2026-02-23)
Completed in `np-npRmpi`:
1. Added explicit tests to prevent recurrence of large-argument autodispatch call-symbol regressions in `npudist(bws=...)`.
2. Scope:
   - `tests/testthat/test-autodispatch-call-helpers.R`
3. Commit:
   - `np-npRmpi`: `1ebe1c0`
4. New contracts:
   - structural contract that `.npRmpi_distributed_call_impl` rewrites return objects using `c(prepared$tmpvals, prepared$prepublish)`.
   - runtime contract forcing prepublish path (`npRmpi.autodispatch.arg.broadcast.threshold = 1L`) and asserting:
     - `.np_eval_bws_call_arg(bw, "dat")` resolves to a full data.frame,
     - `npudist(bws = bw)` completes and yields fitted output length equal to training rows.
5. Validation:
   - `/tmp/nprmpi_test_autod_call_helpers_20260223d.log` (`TEST_RC:0`)
   - `/tmp/nprmpi_test_autod_npudist_filter_20260223d.log` (`TEST_RC:0`)
   - `/tmp/nprmpi_issue_notes_repros_autodtests_20260223.log`
   - `/tmp/nprmpi_check_autodtests_20260223.log` (`Status: OK`)

## DBandwidth `rorder` seq_len Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Replaced fragile `rorder` reconstruction using `(1:ncol)[...]` with zero-length-safe `seq_len(ncol)` indexing.
2. Scope:
   - `R/np.distribution.bw.R`
3. Commit:
   - `np-npRmpi`: `8969569`
4. Validation:
   - parse gate: `PARSE_OK`
   - targeted tests:
     - `/tmp/nprmpi_dist_rorder_tests_targeted_20260223.log` (`PASS 34, FAIL 0`)
   - note on broader formula-bw filter (pre-existing unrelated failures):
     - `/tmp/nprmpi_dist_rorder_tests_20260223.log`
   - issue-note verified repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_dist_rorder_20260223.log`
   - tarball check:
     - `/tmp/nprmpi_check_dist_rorder_20260223.log` (`Status: OK`)

## Density/Regression `rorder` seq_len Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Replaced residual `(1:ncol)[...]` `rorder` reconstruction with `seq_len(ncol)` in core density/regression paths.
2. Scope:
   - `R/np.density.bw.R`
   - `R/np.regression.bw.R`
   - `R/np.regression.R`
3. Commit:
   - `np-npRmpi`: `84a75dd`
4. Validation:
   - parse gates: `PARSE_OK`
   - targeted tests:
     - `/tmp/nprmpi_rorder_regdens_tests_20260223.log` (includes one pre-existing `npudens` formula-path failure unrelated to touched lines)
     - `/tmp/nprmpi_dist_rorder_tests_targeted_20260223.log` (`PASS 34, FAIL 0`)
   - fresh-install smoke:
     - `/tmp/nprmpi_regdens_rorder_smoke_20260223.out` (`NPRMPI_REGDENS_RORDER_SMOKE_OK`)
   - issue-note verified repro sweep:
   - `/tmp/nprmpi_issue_notes_repros_rorder3_20260223.log`
   - tarball check:
     - `/tmp/nprmpi_check_rorder3_20260223.log` (`Status: OK`)

## Autodispatch Formula-Default `gdat` Routing Fix (2026-02-23)
Completed in `np-npRmpi`:
1. Fixed autodispatch target-argument materialization to include default-interface `gdat` (in addition to formula `gdata`), preventing unresolved symbols in worker-evaluated default calls.
2. Kept formal-first argument forcing in materialization path to preserve forwarded-call stability for nested dispatch (`..1`/promise cases).
3. Scope:
   - `R/np.autodispatch.R`
   - `tests/testthat/test-autodispatch-call-helpers.R`
4. Commit:
   - `np-npRmpi`: `356c09c`
5. Validation:
   - focused contracts:
     - `/tmp/nprmpi_formula_autod_helpers_post_gdatfix_20260223.log` (`PASS 40, FAIL 0`)
   - issue-note verified repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_nonstop_20260223.log` (all verified repros passed)
   - user session-mode smoke (`npRmpi.init` + `npcdens`):
     - `/tmp/nprmpi_user_session_smoke_20260223.out` (`SESSION_SMOKE_OK`)
   - tarball build/check:
     - `/tmp/nprmpi_build_gdatfix_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_gdatfix_en0_20260223.log` (`CHECK_RC:0`, `--as-cran`, `FI_*_IFACE=en0`)

## Core Estimator `seq_len` Loop/Index Hardening Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Replaced residual `1:n` loop/index forms in core estimator families with zero-length-safe `seq_len(...)`.
2. Scope:
   - `R/np.distribution.bw.R`
   - `R/np.condistribution.bw.R`
   - `R/np.singleindex.bw.R`
   - `R/np.smoothcoef.bw.R`
   - `R/np.smoothcoef.R`
   - `R/np.plregression.R`
3. Commit:
   - `np-npRmpi`: `46f706f`
4. Validation:
   - parse gates:
     - `/tmp/nprmpi_seq_len_core2_parse_20260223.log` (`RC:0`)
   - targeted contracts:
     - `/tmp/nprmpi_seq_len_core2_tests_postfix_20260223.log` (`PASS 64, FAIL 0`)
   - issue-note verified repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_seqlen_core2_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_seqlen_core2_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_seqlen_core2_20260223.log` (`CHECK_RC:0`, `--as-cran`, `FI_*_IFACE=en0`)

## Autodispatch Missing-Arg Return-Rewrite Hardening (2026-02-23)
Completed in `np-npRmpi`:
1. Fixed return-call tmp replacement recursion to safely handle calls containing missing arguments (`data = , ...`) during autodispatch post-processing.
2. Root cause:
   - recursive replacement evaluated missing-argument placeholders and failed with `argument "xi" is missing, with no default` (surfaced by `npplregbw` path).
3. Changes:
   - `R/np.autodispatch.R`
     - added `.npRmpi_is_missing_call_arg(...)` helper,
     - switched call/pairlist traversal to `as.list(...)` reconstruction with missing-arg skip.
   - `tests/testthat/test-autodispatch-call-helpers.R`
     - added contract test for replacement on calls with missing arguments.
4. Commit:
   - `np-npRmpi`: `4415a9b`
5. Validation:
   - focused contracts:
     - `/tmp/nprmpi_autod_missingarg_fix_focus2_20260223.log` (`TEST_RC:0`)
   - broad targeted contracts (includes `npplreg`):
     - `/tmp/nprmpi_seq_len_core2_tests_postfix_20260223.log` (`PASS 64, FAIL 0`)
   - issue-note verified repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_seqlen_core2_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_check_seqlen_core2_20260223.log` (`CHECK_RC:0`)

## Smoothcoef NDIM `seq_len` Hardening Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Replaced residual `sapply(1:bws$ndim, ...)` with zero-length-safe `sapply(seq_len(bws$ndim), ...)` in core smooth coefficient bandwidth scaling setup.
2. Scope:
   - `R/np.smoothcoef.bw.R`
3. Commit:
   - `np-npRmpi`: `cdd8acd`
4. Validation:
   - parse gate:
     - `/tmp/nprmpi_scoef_ndim_parse_20260223.log` (`RC:0`)
   - targeted contracts:
     - `/tmp/nprmpi_scoef_ndim_tests_20260223.log` (`PASS 54, FAIL 0`)
   - issue-note verified repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_scoef_ndim_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_scoef_ndim_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_scoef_ndim_20260223.log` (`CHECK_RC:0`, `FI_*_IFACE=en0`)

## Sigtest Index/Bootstrap Hygiene Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
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
   - `np-npRmpi`: `08bf9fa`
5. Validation:
   - parse:
     - `/tmp/nprmpi_sigtest_parse_20260223.log` (`RC:0`)
   - focused tests:
     - `/tmp/nprmpi_sigtest_tests_20260223.log` (`PASS 7, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_sigtest_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_sigtest_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_sigtest_20260223.log` (`CHECK_RC:0`, `FI_*_IFACE=en0`)

## Unitest Bootstrap Sampling Hygiene Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Modernized null-draw bootstrap sampling in `npunitest` to use `sample.int(...)` instead of index-vector sampling via `sample(1:length(...), ...)`.
2. Scope:
   - `R/np.unitest.R`
3. Commit:
   - `np-npRmpi`: `a918d1c`
4. Validation:
   - parse:
     - `/tmp/nprmpi_unitest_parse_20260223.log` (`RC:0`)
   - focused tests:
     - `/tmp/nprmpi_unitest_tests_20260223.log` (`PASS 9, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_unitest_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_unitest_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_unitest_20260223.log` (`CHECK_RC:0`, `FI_*_IFACE=en0`)
5. Note:
   - check warning count in this run reflects pre-existing non-target tree state (top-level/codoc drift), not this `np.unitest.R` change.

## Unitest Loop/P-Value Scalar Hygiene Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Replaced residual scalar bootstrap loop/index and indicator-mean pattern in `npunitest`.
2. Scope:
   - `R/np.unitest.R`
3. Changes:
   - `for(b in 1:boot.num)` -> `for (b in seq_len(boot.num))`
   - `mean(ifelse(resampled.stat > test.stat, 1, 0))` -> `mean(resampled.stat > test.stat)`
4. Commit:
   - `np-npRmpi`: `6946033`
5. Validation:
   - parse:
     - `/tmp/nprmpi_unitest2_parse_20260223.log` (`RC:0`)
   - focused tests:
     - `/tmp/nprmpi_unitest2_tests_20260223.log` (`PASS 9, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_unitest2_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_unitest2_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_unitest2_20260223.log` (`CHECK_RC:0`, `FI_*_IFACE=en0`)

## Np*test Bootstrap Loop/Index Hygiene Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
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
   - `np-npRmpi`: `91cae32`
5. Validation:
   - parse gate: `PARSE_OK`
   - focused tests:
     - `/tmp/nprmpi_nptests_fix_20260223.log` (`PASS 13, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_nptestsfix_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_nptestsfix_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_nptestsfix_20260223.log` (`Status: 1 WARNING` pre-fix codoc)

## Sigtest Rd Signature Alignment Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Aligned `npsigtest.rbandwidth` Rd default for `index` with code (`seq_len(ncol(xdat))`).
2. Scope:
   - `man/np.sigtest.Rd`
3. Commit:
   - `np-npRmpi`: `6a0f08d`
4. Validation:
   - tarball-first:
     - `/tmp/nprmpi_check_sigdoc_20260223.log` (`Status: OK`)
   - `--as-cran` refresh:
     - `/tmp/nprmpi_check_ascran_refresh_20260223.log` (`Status: 1 WARNING, 1 NOTE`; no codoc warning).

## Core Estimator Performance Governance Checkpoint (2026-02-23)
Completed in `np-npRmpi` (session mode, `nslaves=1`):
1. Pre/post commit comparison on core estimators with identical scripts/args:
   - pre: `3e4b0da`
   - post: `6feca6d`
2. Workload and policy:
   - script: `/tmp/bench_core_nprmpi_lib_20260223.R`
   - estimators: `npreg`, `npudens`, `npcdens`
   - run mode: session (`npRmpi.init(nslaves=1, quiet=TRUE)`)
   - `n=200`, `times=4`
   - seed policy: fixed (`42`) and varying (`42 + i - 1`)
3. Invocation:
   - `Rscript /tmp/bench_core_nprmpi_lib_20260223.R <lib> <csv> <metrics.rds> 4 200`
   - compare: `Rscript /tmp/bench_compare_prepost_20260223.R <pre_csv> <post_csv> <pre_rds> <post_rds> <summary_csv> <parity_txt>`
4. Mean/median percent change summary:
   - fixed:
     - `npreg`: mean `-2.07%`, median `-3.87%`
     - `npudens`: mean `-2.44%`, median `-6.20%`
     - `npcdens`: mean `-4.49%`, median `-5.24%`
   - varying:
     - `npreg`: mean `+3.50%`, median `+2.78%`
     - `npudens`: mean `+4.78%`, median `+0.81%`
     - `npcdens`: mean `+4.36%`, median `+7.14%`
5. Parity:
   - fixed-seed metric parity max abs diff: `0`
   - note: `npcdens_bw_mean` is non-finite (`NA`) in both runs; fitted-value parity remained exact.
6. Artifacts:
   - pre/post raw:
     - `/tmp/bench_nprmpi_pre_full_20260223.csv`
     - `/tmp/bench_nprmpi_post_full_20260223.csv`
   - summary/parity:
     - `/tmp/bench_nprmpi_full_summary_20260223.csv`
     - `/tmp/bench_nprmpi_full_parity_20260223.txt`

## Win-Builder Submission Checkpoint (2026-02-23)
Submitted from this environment:
1. `np-master` release check submitted via:
   - `Rscript -e \"suppressPackageStartupMessages(library(devtools)); check_win_release('/Users/jracine/Development/np-master')\"`
   - log: `/tmp/winbuilder_submit_np_20260223.log`
2. `np-npRmpi` release check submitted via:
   - `Rscript -e \"suppressPackageStartupMessages(library(devtools)); check_win_release('/Users/jracine/Development/np-npRmpi')\"`
   - log: `/tmp/winbuilder_submit_nprmpi_20260223.log`
3. Status:
   - submissions accepted by win-builder; final pass/fail results pending email report (~15-30 minutes post submission).

## Plot Engine Scalar-Branch Hygiene Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
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
   - `np-npRmpi`: `75f3e84`
5. Validation:
   - parse gate: `PARSE_OK`
   - focused tests:
     - `/tmp/nprmpi_plot_scalar_ifelse_tests_20260223.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_plotscalar_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_plotscalar_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_plotscalar_20260223.log` (`Status: OK`)

## Plot Engine Scalar-`ifelse` Expansion Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
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
   - `np-npRmpi`: `f4b177b`
5. Validation:
   - parse gate: `PARSE_OK`
   - focused tests:
     - `/tmp/nprmpi_plot_ifelse2_tests_20260223.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_plotifelse2_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_plotifelse2_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_plotifelse2_20260223.log` (`Status: OK`)

## Plot Bandwidth/Dbandwidth Scalar-Branch Checkpoint (2026-02-23)
Completed in `np-npRmpi`:
1. Replaced residual scalar `ifelse(...)` branches in density/regression bandwidth plot engines.
2. Scope:
   - `R/np.plot.engine.bandwidth.R`
   - `R/np.plot.engine.dbandwidth.R`
3. Changes:
   - scalar `xi.factor` error-style/error-bar selectors now use scalar `if`
   - scalar `xi.factor` linetype selector in `dbandwidth` now uses scalar `if`
4. Commit:
   - `np-npRmpi`: `cc25213`
5. Validation:
   - parse gate: `PARSE_OK`
   - focused tests:
     - `/tmp/nprmpi_plot_ifelse3_tests_20260223.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_plotifelse3_20260223.log` (all verified repros passed)
   - tarball-first:
     - `/tmp/nprmpi_build_plotifelse3_20260223.log` (`BUILD_RC:0`)
     - `/tmp/nprmpi_check_plotifelse3_20260223.log` (`Status: OK`)

## Plot `rbandwidth` Scalar-Branch + Vector-Safe Error-Bar Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
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
   - `np-npRmpi`: `4847719`
5. Validation:
   - targeted tests:
     - `/tmp/nprmpi_plot_ifelse4_tests2_20260224.log` (`PASS 17, FAIL 0, SKIP 1`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_plotifelse4_20260224.log` (`RC:0`)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_plotifelse4_20260224.log` (`BUILD_RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_plotifelse4_20260224.log` (`Status: OK`)

## Indicator/Scalar `ifelse` Cleanup Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced remaining scalar/indicator `ifelse(...)` in active runtime paths:
   - `R/np.plot.engine.sibandwidth.R`: scalar gradient-name prefix branch.
   - `R/np.cmstest.R`, `R/np.qcmstest.R`, `R/np.symtest.R`: p-value indicator means now use direct logical means.
2. Changes:
   - `paste(ifelse(gradients, ...))` -> scalar `if` branch.
   - `mean(ifelse(condition, 1, 0))` -> `mean(condition)`.
3. Commit:
   - `np-npRmpi`: `528c639`
4. Validation:
   - targeted tests:
     - `/tmp/nprmpi_nptests_plotifelse5_20260224.log` (`PASS 18, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_ifelse5_20260224.log` (`RC:0`)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_ifelse5_20260224.log` (`BUILD_RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ifelse5_20260224.log` (`Status: OK`)

## Print/Reporting Scalar-Branch Cleanup Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced scalar `ifelse(...)` in print/reporting paths with explicit scalar `if` branches.
2. Scope:
   - `R/cmstest.R`
   - `R/deneqtest.R`
   - `R/deptest.R`
   - `R/sdeptest.R`
   - `R/symtest.R`
   - `R/unitest.R`
3. Commit:
   - `np-npRmpi`: `513057c`
4. Validation:
   - parse gates: `PARSE_OK`
   - targeted tests:
     - `/tmp/nprmpi_print_ifelse6_tests_20260224.log` (`PASS 13, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_printifelse6_20260224.log` (`RC:0`)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_printifelse6_20260224.log` (`BUILD_RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_printifelse6_20260224.log` (`Status: OK`)

## `b.star` Elementwise Rounding/Bounds Fix Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Fixed a multivariate rounding/bounds bug in `b.star(..., round=TRUE)` where `BstarCB` used scalar `max(1, round(BstarCB))` inside vectorized `ifelse`, which could collapse column-wise outputs.
2. Modernized bounds logic with elementwise vector ops:
   - `M <- min(2 * mhat, mmax)`
   - `pmin(...)` / `pmax(...)` for SB/CB bounds.
3. Added regression coverage for multivariate `round=TRUE` elementwise behavior.
4. Scope:
   - `R/b.star.R`
   - `tests/testthat/test-utils.R`
5. Commit:
   - `np-npRmpi`: `4e78cae`
6. Validation:
   - focused tests:
     - `/tmp/nprmpi_bstar_fix_tests_20260224.log` (`PASS 6, FAIL 0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_bstarfix_20260224.log` (`RC:0`)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_bstarfix_20260224.log` (`BUILD_RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_bstarfix_20260224.log` (`Status: OK`)

## Session Hat/Bootstrap Routing Guard Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Added bounded subprocess contracts to lock session-mode routing for modern hat/bootstrap surfaces:
   - `npreghat` parity smoke (`npRmpi.init` -> `npregbw/npreg/npreghat` -> `npRmpi.quit`).
   - `plot(..., plot.errors.boot.method="wild-hat")` session smoke.
2. Extended periodic verified-issue repro harness with matching session guards so resolved regressions are rechecked by default.
3. Scope:
   - `tests/testthat/test-session-routing-subprocess-contract.R`
   - `issue_notes/verified_issue_repros.R`
4. Commits:
   - `np-npRmpi`: `a550e9b`
   - `np-npRmpi`: `5635b05`
5. Validation:
   - subprocess contracts (`NOT_CRAN=true`):
     - `/tmp` console run (8 pass, attach-mode test skipped unless `NP_RMPI_ENABLE_ATTACH_TEST=1`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_032316.log` (all verified repros passed)

## Quantile/Copula Loop-Index Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced residual `1:length(...)` / `1:n...` loop patterns in auxiliary quantile/copula runtime paths with `seq_along(...)` / `seq_len(...)`.
2. Scope:
   - `R/np.quantile.R`
   - `R/np.copula.R`
3. Commit:
   - `np-npRmpi`: `4a2376a`
4. Validation:
   - targeted tests:
     - `/tmp/nprmpi_copula_quantile_tests_20260224.log` (`RC:0`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_032635.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_loophyg_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_loophyg_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)
