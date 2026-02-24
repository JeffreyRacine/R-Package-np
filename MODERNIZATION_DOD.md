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

## Current Status (2026-02-24)
1. `npRmpi` modernization is in late routing/runtime hardening mode:
   - core estimator paths are stable with subprocess session contracts and periodic issue-note sweeps,
   - recent patches continue to prioritize regression prevention in user-facing `npRmpi.init(...)` workflows.
2. Recent completed slices include:
   - default-quiet `npRmpi.init(nslaves=1)` + `npcdens` subprocess contract lock,
   - GLP/`dim_basis` contract hardening parity with `np-master`,
   - deterministic NA/range guard hardening and helper micro-modernization.
3. Remaining highest-priority work:
   - bounded-kernel/convolution native-path completion and MPI-path parity validation,
   - performance-governance closure with fixed/varying seed comparisons for performance-sensitive changes,
   - win-builder closure before release candidate.

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

## Utility Loop-Index Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced residual `1:length(...)` / `1:ncol(...)` index patterns in shared utility/runtime helpers with `seq_along(...)` / `seq_len(...)`.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-npRmpi`: `3560c2b`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='utils|npcopula|npudist|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_033540.log` (all verified repros passed)

## Univariate Density Loop-Index Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced residual `1:length(...)` / `1:n...` loop/index patterns in univariate boundary/shape-constrained density helpers with `seq_along(...)` / `seq_len(...)`.
2. Scope:
   - `R/npuniden.boundary.R`
   - `R/npuniden.sc.R`
3. Commit:
   - `np-npRmpi`: `019793f`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='npuniden|utils', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_033821.log` (all verified repros passed)

## Plot-Helper Indexing Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced residual `1:length(...)` / `1:ncol(...)` indexing in plot helper paths with `seq_along(...)` / `seq_len(...)`.
2. Replaced manual cumulative-loop `sapply(1:length(tq), ...)` with `cumsum(tq)` in quantile helper.
3. Scope:
   - `R/np.plot.helpers.R`
4. Commit:
   - `np-npRmpi`: `2b6b0c6`
5. Validation:
   - targeted tests:
     - `testthat::test_local(filter='semihat|npreghat|plot|npindex|npplreg|npscoef', reporter='summary')` (`RC:0`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_034151.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_posthyg3_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_posthyg3_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## Utility/Test-Stat Loop Finalization Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
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
   - `np-npRmpi`: `a0cc555`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='utils|nptests|sdeptest|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_035016.log` (all verified repros passed)

## Plot-Helper Residual Index Cleanup Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Removed remaining `1:length(...)` index helpers in low-level plot helper vectors (`jj` and factor plotting separators) using `seq_along(...)`.
2. Scope:
   - `R/np.plot.helpers.R`
3. Commit:
   - `np-npRmpi`: `5bd980b`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='plot-autodispatch|semihat', reporter='summary')` (`RC:0`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_035414.log` (all verified repros passed)

## MPI Helper Loop-Bounds Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced residual `1:...` loop/index ranges in MPI utility helpers with `seq_len(...)`-safe forms to eliminate `1:0` style edge behavior.
2. Scope:
   - `R/Rcoll.R`
   - `R/Rparutilities.R`
3. Commit:
   - `np-npRmpi`: `855b097`
4. Validation:
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='mpi-helpers|mpi-mixed-mode-guards|session-routing-subprocess-contract|plot-autodispatch|semihat|rngstream-contract', reporter='summary')` (`RC:0`; attach smoke intentionally skipped unless `NP_RMPI_ENABLE_ATTACH_TEST=1`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_040132.log` (all verified repros passed)

## Consolidated Tarball Gate (2026-02-24)
1. `R CMD build /Users/jracine/Development/np-npRmpi`
   - `/tmp/nprmpi_build_postmodern_final_20260224.log` (`RC:0`, `creating vignettes ... OK`)
2. `FI_TCP_IFACE=en0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=en0 R CMD check --as-cran npRmpi_0.70-0.tar.gz`
   - `/tmp/nprmpi_check_ascran_postmodern_final_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## Plot-Engine Loop-Range Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced residual `1:...` loop/index ranges in plot-engine runtime paths with `seq_len(...)`-safe forms.
2. Scope:
   - `R/np.plot.engine.bandwidth.R`
   - `R/np.plot.engine.dbandwidth.R`
   - `R/np.plot.engine.conbandwidth.R`
   - `R/np.plot.engine.condbandwidth.R`
3. Commit:
   - `np-npRmpi`: `6fdeaa9`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='plot-autodispatch|semihat', reporter='summary')` (`RC:0`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_041140.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_plotengine_seq_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_plotengine_seq_20260224_en0.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)
     - `/tmp/nprmpi_check_ascran_plotengine_seq_20260224.log` (`RC:1`; `lo0` NIC setting failed during MPI init in install stage, retried successfully with `en0`)

## `npregiv` Loop-Range Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced residual `1:ncol(...)` / `1:length(...)` index ranges in `npregiv` support paths with `seq_len(...)` / `seq_along(...)` to avoid `1:0` edge behavior.
2. Scope:
   - `R/npregiv.R`
3. Commit:
   - `np-npRmpi`: `4a1fe8f`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='npregiv', reporter='summary')` (`RC:0`; expected pre-existing monotone-stopping warnings)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_042226.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_npregiv_seq_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_npregiv_seq_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `npregiv` Loop-Range Hygiene Follow-On Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Completed residual `1:n...` conversions in `npregiv` iterative/LOO loops and `sapply` index ranges.
2. Scope:
   - `R/npregiv.R`
3. Commit:
   - `np-npRmpi`: `a7d61d8`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='npregiv', reporter='summary')` (`RC:0`; expected pre-existing monotone-stopping warnings)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_043014.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_npregiv_seq2_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_npregiv_seq2_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## Bandwidth Metadata / `npregivderiv` Loop-Range Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced residual metadata/index range builders with `seq_len(...)` in bandwidth metadata constructors and `npregivderiv` numeric-column detection.
2. Scope:
   - `R/dbandwidth.R`
   - `R/rbandwidth.R`
   - `R/smoothbandwidth.R`
   - `R/npregivderiv.R`
3. Commit:
   - `np-npRmpi`: `04199a4`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='npregiv|bandwidth|bw-dispatch|formula-bw-contract', reporter='summary')` (`RC:0`; expected pre-existing monotone-stopping warnings)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_043811.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_bwmeta_seq_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_bwmeta_seq_20260224_withcrs.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)
5. Environment note:
   - local `--as-cran` checks required `crs` in `R_LIBS`; commands used:
     - `R_LIBS=/tmp/crs_check_lib FI_TCP_IFACE=en0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=en0 R CMD check --as-cran --no-manual npRmpi_0.70-0.tar.gz`

## Bandwidth Metadata Parity Follow-On Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Applied matching `seq_len(...)` metadata-index conversion in `R/bandwidth.R` for cross-file parity.
2. Scope:
   - `R/bandwidth.R`
3. Commit:
   - `np-npRmpi`: `457098b`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='bw-dispatch|formula-bw-contract|bandwidth', reporter='summary')` (`RC:0`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_044753.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_bwmeta_seq2_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_bwmeta_seq2_20260224.log` (`RC:0`, `Status: 1 WARNING, 1 NOTE`; warning set unchanged from existing top-level-file debt)

## `b.star` Index-Range Guard Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced the `1:(mmax-Kn+1)` insignificant-run index builder with bounded `seq_len(max(mmax - Kn + 1L, 0L))` to avoid `1:0`/descending edge behavior while preserving normal positive-range semantics.
2. Scope:
   - `R/b.star.R`
3. Commit:
   - `np-npRmpi`: `21627e1`
4. Validation:
   - targeted tests:
     - `testthat::test_local(filter='utils|nptests|sdeptest|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_045406.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_bstar_seq_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_bstar_seq_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `npuniden.boundary` Scalar-Branch Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced scalar `ifelse(...)` branches in `npuniden.boundary` CV/optimizer routing with scalar `if` expressions:
   - finite fallback in LS-CV objective,
   - optimizer start-index choice (`cv.ml` max vs `cv.ls` min),
   - optimizer upper-bound choice for `beta2` kernels.
2. Scope:
   - `R/npuniden.boundary.R`
3. Commit:
   - `np-npRmpi`: `03d0e3f`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/npuniden.boundary.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='npuniden|nptests|sdeptest')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
     - log: `/tmp/nprmpi_npuniden_scalar_ifelse_tests_20260224.log`
   - direct runtime smoke:
     - `/tmp/nprmpi_npuniden_scalar_ifelse_smoke_20260224.out` (`NPRMPI_NPUNIDEN_BOUNDARY_SMOKE_OK`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_050442.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_npuniden_scalar_ifelse_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_npuniden_scalar_ifelse_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## Wild-Bootstrap Draw Allocation Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Reduced allocations in wild-bootstrap draw helpers by replacing matrix `ifelse(...)` draw generation with preallocated matrices and logical indexing.
2. Scope:
   - `R/np.plot.helpers.R`
3. Commit:
   - `np-npRmpi`: `6fe2c3e`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/np.plot.helpers.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='plot-autodispatch|semihat', reporter='summary')` (`RC:0`)
     - log: `/tmp/nprmpi_wilddraw_alloc_tests_20260224b.log`
   - direct helper contract smoke:
     - `/tmp/nprmpi_wilddraw_contract_20260224.out` (`NPRMPI_WILDDRAW_CONTRACT_OK`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_051458.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_wilddraw_alloc_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_wilddraw_alloc_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## Sequence-Range Edge-Safety Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Eliminated `1:n` edge behavior in key utility/MPI helpers:
   - `dlev()` now uses `as.numeric(seq_len(nlevels(x)))`,
   - `.splitIndices()` now uses `seq_len/seq.int`, fixing the `nx=0` phantom index bug.
2. Scope:
   - `R/util.R`
   - `R/Rparutilities.R`
3. Commit:
   - `np-npRmpi`: `3fc6792`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); invisible(parse(file='R/Rparutilities.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='mpi-helpers|nptests|utils', reporter='summary')` (`RC:0`)
     - log: `/tmp/nprmpi_splitindices_dlev_tests_20260224.log`
   - direct edge-case smoke:
     - `/tmp/nprmpi_splitindices_zero_smoke_20260224.out` (`NPRMPI_SPLITINDICES_ZERO_OK`)
     - `/tmp/np_dlev_seq_smoke_20260224.out` (`DLEV_SEQ_OK`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_052223.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_seq_edge_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_seq_edge_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## Bootstrap Wild-Draw Refactor Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced repeated `ifelse(rbinom(...))` bootstrap multipliers with a lightweight preallocated vector helper (`runif` + logical assignment) in test/bootstrap-heavy paths.
2. Scope:
   - `R/np.cmstest.R`
   - `R/np.qcmstest.R`
   - `R/np.sigtest.R`
3. Commit:
   - `np-npRmpi`: `b07921f`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/np.cmstest.R')); invisible(parse(file='R/np.qcmstest.R')); invisible(parse(file='R/np.sigtest.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='cmstest|qcmstest|sigtest', reporter='summary')` (`RC:0`)
     - log: `/tmp/nprmpi_bootdraw_refactor_tests_20260224.log`
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_053149.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_bootdraw_refactor_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_bootdraw_refactor_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `npqcmstest` Quantile-Residual Helper Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced repeated `ifelse(model.resid <= 0, 1 - tau, -tau)` logic with a local helper preserving NA semantics, reused across bandwidth selection and test-statistic internals.
2. Scope:
   - `R/np.qcmstest.R`
3. Commit:
   - `np-npRmpi`: `013efc9`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/np.qcmstest.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='qcmstest|cmstest|sigtest', reporter='summary')` (`RC:0`)
     - log: `/tmp/nprmpi_qresidual_helper_tests_20260224.log`
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_053830.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_qresidual_helper_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_qresidual_helper_20260224.log` (`RC:0`, `Status: 1 WARNING, 1 NOTE`; warning set unchanged from existing top-level-file debt)

## `npuniden.boundary` + `npregiv` Edge-Safety Follow-On (2026-02-24)
Completed in `np-npRmpi`:
1. Hardened remaining low-risk edge paths:
   - replaced vector `ifelse(...)` finite clamp in CV-ML objective with explicit mask assignment (`f.safe`) to reduce branch allocation and keep invalid values at `.Machine$double.xmin`,
   - replaced debug-only `seq(1:num.bw)` with `seq_len(num.bw)` in `npregiv`.
2. Scope:
   - `R/npuniden.boundary.R`
   - `R/npregiv.R`
3. Commit:
   - `np-npRmpi`: `968e73c`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/npuniden.boundary.R')); invisible(parse(file='R/npregiv.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='npuniden|npregiv', reporter='summary')` (`RC:0`; expected pre-existing `npregiv` and `npuniden.sc` warnings)
     - log: `/tmp/nprmpi_uniden_npregiv_edge_tests_20260224.log`
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_054656.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_uniden_npregiv_edge_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_uniden_npregiv_edge_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `genGofStr` Scalar-Branch Hygiene Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced nested scalar `ifelse(...)` in `genGofStr()` with explicit scalar `if` branches for `MSE`/`R2` report fragments.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-npRmpi`: `333d64a`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='utils|nptests|sdeptest|npuniden', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
     - log: `/tmp/nprmpi_gofstr_scalar_tests_20260224.log`
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_055424.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_gofstr_scalar_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_gofstr_scalar_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `util.R` Scalar-Branch Follow-On (`NZD`/`genBwScaleStrs`) (2026-02-24)
Completed in `np-npRmpi`:
1. Removed additional scalar `ifelse(...)` usage in hot/summary utility paths:
   - `NZD()` small-value assignment now uses direct sign-mask replacement,
   - `genBwScaleStrs()` scalar label/summary string branches now use explicit scalar `if`.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-npRmpi`: `162b09e`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='utils|bandwidth|bw-dispatch|formula-bw-contract', reporter='summary')` (`RC:0`)
     - log: `/tmp/nprmpi_util_scalar2_tests_20260224.log`
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_060257.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_util_scalar2_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_util_scalar2_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `util.R` Vector-Branch Follow-On (`tgauss`/`pad`/`rpad`) (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced additional `ifelse(...)` usages with explicit vector-safe masked assignments in utility hot paths:
   - `nptgauss()` local `tgauss` now uses precomputed vector + mask assignment,
   - `pad()` and `rpad()` now use vector-safe index assignment while preserving names.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-npRmpi`: `687bb21`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); cat('PARSE_OK\n')"` (`RC:0`)
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='utils|nptests|sdeptest|npuniden|bandwidth', reporter='summary')` (`RC:0`; expected pre-existing `npuniden.sc` warning)
     - log: `/tmp/nprmpi_util_scalar3_tests_20260224.log`
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_061008.log` (all verified repros passed)
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_util_scalar3_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_util_scalar3_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `util.R` Scalar-Branch Follow-On (`genBwKerStrs`/`genBwKerStrsXY`) (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced nested `ifelse(...)` branch construction in kernel-summary string helpers with explicit scalar/vector-safe branching.
2. Scope:
   - `R/util.R`
3. Commit:
   - `np-npRmpi`: `0546ebe`
4. Validation:
   - parse gate:
     - `Rscript -e "invisible(parse(file='R/util.R')); cat('PARSE_OK\n')"` (`RC:0`)
     - log: `/tmp/nprmpi_util_genBwKerStrs_parse_20260224.log`
   - targeted tests (`NOT_CRAN=true`):
     - `testthat::test_local(filter='plot-contract|formula-bw-contract|bw-dispatch-contract|call-bounds-guard-contract|npksum|session-arg-contract')` (`RC:0`)
     - log: `/tmp/nprmpi_util_genBwKerStrs_tests2_20260224.log`
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_util_genBwKerStrs.log` (all verified repros passed)
     - run artifact: `/tmp/nprmpi_issue_notes_repros_20260224_062150.log`
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_util_genBwKerStrs_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_util_genBwKerStrs_withcrs_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## Core NA/Range Guard Hardening (`npquantile`/`npcopula`/`npsigtest`/`npreg`) (2026-02-24)
Completed in `np-npRmpi`:
1. Hardened argument-validation paths to fail deterministically on missing values and fixed one ordering bug in `npquantile`:
   - `tau` validation now rejects missing values explicitly and range-checks with `na.rm = TRUE`,
   - `npquantile` now validates `bws$xnames` only after constructing/validating `bws`,
   - `u`/`prob`/`index` range checks now guard missing values explicitly,
   - Bernstein-support guard in `npreg` now uses `na.rm = TRUE` for deterministic behavior.
2. Scope:
   - `R/np.quantile.R`
   - `R/np.copula.R`
   - `R/np.plot.helpers.R`
   - `R/np.sigtest.R`
   - `R/np.regression.R`
3. Commit:
   - `np-npRmpi`: `07c93cb`
4. Validation:
   - parse gate:
     - log: `/tmp/nprmpi_na_guard_parse_20260224.log` (`PARSE_OK`)
   - targeted tests (`NOT_CRAN=true`):
     - logs:
       - `/tmp/nprmpi_na_guard_tests_20260224.log` (`partial in this environment due long-running MPI `npreg-glp-higher-order` context`)
       - `/tmp/nprmpi_na_guard_tests2_20260224.log` (`FAIL 0`, `PASS 42` on MPI-safe filter)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_na_guard.log` (all verified repros passed)
     - run artifact: `/tmp/nprmpi_issue_notes_repros_20260224_063809.log`
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_na_guard_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_na_guard_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `npuniden*` NA-Safe Range Validation Hardening (2026-02-24)
Completed in `np-npRmpi`:
1. Added deterministic missing-value guards and NA-safe range checks for univariate density helper entry points.
2. Scope:
   - `R/npuniden.boundary.R`
   - `R/npuniden.reflect.R`
   - `R/npuniden.sc.R`
3. Commit:
   - `np-npRmpi`: `565ee2e`
4. Validation:
   - parse gate:
     - log: `/tmp/nprmpi_uniden_na_guard_parse_20260224.log` (`PARSE_OK`)
   - targeted tests (`NOT_CRAN=true`):
     - log: `/tmp/nprmpi_uniden_na_guard_tests_20260224.log` (`FAIL 0`; expected pre-existing `npuniden.sc` warning only)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_uniden_na_guard.log` (all verified repros passed)
     - run artifact: `/tmp/nprmpi_issue_notes_repros_20260224_064531.log`
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_uniden_na_guard_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_uniden_na_guard_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## GLP Validator Typing/NA Contract Hardening (2026-02-24)
Completed in `np-npRmpi`:
1. Hardened GLP validation helpers to reject non-numeric and missing derivative/degree vectors deterministically before numeric comparisons.
2. Scope:
   - `R/util.R`
   - `tests/testthat/test-glp-validator-contract.R`
3. Commit:
   - `np-npRmpi`: `990dbc5`
4. Validation:
   - parse gate:
     - `/tmp/nprmpi_glp_validator_parse_20260224.log` (`PARSE_OK`)
   - targeted tests (`NOT_CRAN=true`):
     - `/tmp/nprmpi_glp_validator_tests_20260224.log` (`FAIL 0`, `PASS 32`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_glp_validator.log` (all verified repros passed)
     - run artifact: `/tmp/nprmpi_issue_notes_repros_20260224_065332.log`
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_glp_validator_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_glp_validator_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## `dim_basis` Integer-Like Contract Hardening (2026-02-24)
Completed in `np-npRmpi`:
1. Hardened `dim_basis()` input contracts to validate numeric/integer-like vectors before coercion for `degree`, `segments`, `include`, and `categories`.
2. This closes a latent coercion hazard where character inputs could be silently coerced and fail later with less precise diagnostics.
3. Scope:
   - `R/util.R`
   - `tests/testthat/test-glp-validator-contract.R`
4. Commit:
   - `np-npRmpi`: `a4302df`
5. Validation:
   - parse gate:
     - `/tmp/nprmpi_dimbasis_parse_20260224.log` (`PARSE_OK`)
   - targeted tests (`NOT_CRAN=true`):
     - `/tmp/nprmpi_dimbasis_tests_20260224.log` (`FAIL 0`, `PASS 36`)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_20260224_dimbasis.log` (all verified repros passed)
     - run artifact: `/tmp/nprmpi_issue_notes_repros_20260224_070017.log`
   - tarball-first (MPI env pinned):
     - `/tmp/nprmpi_build_dimbasis_20260224.log` (`RC:0`, `creating vignettes ... OK`)
     - `/tmp/nprmpi_check_ascran_dimbasis_20260224.log` (`RC:0`, `Status: 1 WARNING, 2 NOTEs`; warning set unchanged from existing top-level-file debt)

## Session Routing Contract Extension (Default-Quiet `npcdens`) (2026-02-24)
Completed in `np-npRmpi`:
1. Added a subprocess regression contract for user-style session initialization with default `npRmpi.init(nslaves=1)` (no explicit `quiet=`), followed by `npcdens` + `plot`.
2. This locks the exact class of session hang regression seen when examples are run interactively without an explicit `quiet` argument.
3. Scope:
   - `tests/testthat/test-session-routing-subprocess-contract.R`
4. Commit:
   - `np-npRmpi`: `4ec0b8c`
5. Validation:
   - focused contracts (`NOT_CRAN=true`, MPI env pinned):
     - `/tmp/nprmpi_session_routing_contract_tests_20260224b.log` (`RC:0`)
   - issue-note repro sweep (post-change):
     - `/tmp/nprmpi_issue_notes_repros_20260224_defaultquiet_contract.log` (all verified repros passed)
     - run artifact: `/tmp/nprmpi_issue_notes_repros_20260224_071246.log`

## `uocquantile` NA-Guard Micro-Modernization (`anyNA`) (2026-02-24)
Completed in `np-npRmpi`:
1. Replaced `any(is.na(x) | is.nan(x))` with `anyNA(x)` in `uocquantile` for clearer semantics and lower overhead in a hot helper path.
2. Scope:
   - `R/np.plot.helpers.R`
3. Commit:
   - `np-npRmpi`: `8c823eb`
4. Validation:
   - targeted tests (`NOT_CRAN=true`, MPI env pinned):
     - `/tmp/nprmpi_plothelpers_anyna_tests_20260224.log` (`RC:0`)
   - issue-note repro sweep:
     - `/tmp/nprmpi_issue_notes_repros_20260224_plothelpers_anyna.log` (all verified repros passed)
     - run artifact: `/tmp/nprmpi_issue_notes_repros_20260224_071545.log`

## Session Plot Bootstrap `inid` Stall Fix (2026-02-24)
Completed in `np-npRmpi`:
1. Fixed user-facing `plot(npreg(...), plot.errors.method="bootstrap", plot.errors.boot.method="inid")` session stalls by removing MPI/distributed refit pressure from bootstrap inner loops.
2. Changes:
   - Added `.npRmpi_with_local_bootstrap(...)` wrapper to disable autodispatch during bootstrap/tsboot loops.
   - Added `.npRmpi_bootstrap_estimator(...)` resolver to use serial `np` estimators (`np::npreg`, `np::npreghat`) when available, with safe fallback to `npRmpi` implementations.
   - Applied local-bootstrap wrapper across all `compute.bootstrap.errors.*` bootstrap branches.
   - Added explicit session subprocess contract for `inid` plot path and mirrored issue-note smoke guard.
3. Scope:
   - `R/np.plot.helpers.R`
   - `tests/testthat/test-session-routing-subprocess-contract.R`
   - `issue_notes/verified_issue_repros.R`
4. Commit:
   - `np-npRmpi`: `3da77d0`
5. Validation:
   - exact user-style repro (`n=1000`, `boot.num=9999`) now completes:
     - `/tmp/repro_plot_boot_npRmpi_fix3_20260224.out`
     - elapsed:
       - plain `plot(g)`: `0.794s`
       - wild-hat bootstrap: `3.286s`
       - inid bootstrap: `9.307s`
   - issue-note sweep with inid plot guard:
     - `/tmp/nprmpi_issue_notes_repros_20260224_inidguard3.log` (`REPRO_RC:0`)
     - run artifact: `/tmp/nprmpi_issue_notes_repros_20260224_073703.log`

## `npRmpi` `inid` Plot Path: Fixed-`lc` Fast Path + Session Guards (2026-02-24)
Completed in `np-npRmpi`:
1. Ported the fixed-`lc` `inid` counts-weighted hat-operator fast path to `compute.bootstrap.errors.rbandwidth` under the same contract gate as `np-master`.
2. Added robust fallback semantics:
   - if `npreghat` is unavailable/fails/returns incompatible shape, code falls back to legacy bootstrap-refit route.
3. Strengthened routing regression guards to lock this exact route:
   - session subprocess contract now builds `bw <- npregbw(..., regtype='lc', bws=0.25, bandwidth.compute=FALSE)` before `plot(..., plot.errors.boot.method='inid')`,
   - issue-note verified repro `#inidplot` now mirrors the same fixed-`lc` setup.
4. Scope:
   - `R/np.plot.helpers.R`
   - `tests/testthat/test-session-routing-subprocess-contract.R`
   - `issue_notes/verified_issue_repros.R`
5. Validation:
   - install/smoke in temp library:
     - `/tmp/npfast_install_nprmpi.log` (`RC:0`)
     - `/tmp/npfast_nprmpi_validate.out` (`NPRMPI_VALIDATE_OK`)
   - subprocess routing contracts (`NOT_CRAN=true`, MPI env pinned):
     - `/tmp/nprmpi_test_session_routing_20260224.log` (`FAIL 0`, one expected attach test skip)
   - issue-note regression sweep:
     - `/tmp/nprmpi_issue_notes_repros_fastinid_20260224_080637.log` (`all verified repros passed`)
6. Environment note:
   - in this sandbox, `FI_TCP_IFACE=en0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=en0` was required for stable MPI init.

## Plot Bootstrap Modernization: `wild` Clean Break + Chunked `inid` Fast Path (2026-02-24)
Completed in `np-npRmpi`:
1. Standardized plot-bootstrap method naming on `plot.errors.boot.method="wild"` as a clean break (deprecated `wild-hat` removed).
2. Ported chunked fixed-`lc` `inid` fast-path execution in `compute.bootstrap.errors.rbandwidth` to bound memory for large bootstrap jobs:
   - new option: `np.plot.inid.chunk.size` (positive integer).
3. Updated session subprocess contract to exercise `wild` alias in user-mode routing tests.
4. Scope:
   - `R/np.plot.helpers.R`
   - `tests/testthat/test-session-routing-subprocess-contract.R`
   - `tests/testthat/test-semihat.R`
   - `issue_notes/verified_issue_repros.R`
5. Validation:
   - targeted session-routing subprocess contracts (`NOT_CRAN=true`, MPI env pinned):
     - `/tmp/npwildchunk_test_nprmpi_20260224.log` (`FAIL 0`, one expected attach-mode skip)
   - issue-note repro sweep (MPI env pinned):
     - `/tmp/nprmpi_issue_notes_repros_wildchunk_20260224_083729.log` (all verified repros passed)
   - bounded performance benchmark (`n=10000`, `boot.num=999`, session `npRmpi.init(nslaves=1)`, fixed+varying seeds):
     - raw: `/tmp/bench_nprmpi_inid_n10000_b999_bounded_raw_20260224.csv`
     - summary: `/tmp/bench_nprmpi_inid_n10000_b999_bounded_summary_20260224.txt`
     - fixed-seed mean/median delta: `-94.83% / -94.83%`
     - varying-seed mean/median delta: `-94.54% / -94.54%`
6. Runtime note:
   - sandbox MPI execution used `FI_TCP_IFACE=en0 FI_PROVIDER=tcp FI_SOCKETS_IFACE=en0`.

## Plot Bootstrap MPI Fan-Out Probe (2026-02-24)
Completed in `np-npRmpi`:
1. Added chunk-fanout scaffolding for bootstrap `wild`/`inid` paths in `R/np.plot.helpers.R`.
2. Found reproducible hangs in session-mode (`npRmpi.init(nslaves=1)`) when using apply/remote-exec style fanout primitives in this runtime.
3. Implemented a safety gate so fanout is disabled by default:
   - `options(np.plot.bootstrap.mpi.experimental = TRUE)` is now required to activate the experimental path.
4. Default behavior remains the known-safe local/chunked execution path.
5. Validation:
   - session script (`npreg` plot timing with `wild` + `inid`) completes without hang:
     - `/tmp/repro_plot_npreg_session_20260224_gated.out`
   - session `npcdens` plot example completes without hang:
     - `/tmp/repro_user_npcdens_session_20260224.out`

## Plot Bootstrap `npksum` Inid Fast Path (Staged, Default-Off) (2026-02-24)
Completed in `np-npRmpi`:
1. Ported `np-master` helper tranche for inid `npksum` acceleration into `R/np.plot.helpers.R`:
   - `.np_boot_matrix_from_ksum`
   - `.np_inid_boot_from_ksum_unconditional`
   - `.np_inid_boot_from_ksum_conditional`
   - supporting conditional kbandwidth helper constructors.
2. Wired routes in:
   - `compute.bootstrap.errors.bandwidth`
   - `compute.bootstrap.errors.dbandwidth`
   - `compute.bootstrap.errors.conbandwidth`
   - `compute.bootstrap.errors.condbandwidth`
3. Runtime finding:
   - in this session-mode runtime, direct `npksum(...)` usage under density inid plot paths can stall/hang.
4. Safety gate:
   - fast path is now explicit opt-in and disabled by default:
     - `options(np.plot.inid.ksum.fastpath.nprmpi = TRUE)`.
   - opt-in path uses serial `np` namespace kernels (`np::npksum` + `np` kbandwidth helpers)
     via `.npRmpi_bootstrap_estimator(...)` to avoid session-mode `npRmpi::npksum` stalls.
5. Contract strategy:
   - default subprocess routing contracts stay in the stable path,
   - density inid session smoke added as opt-in contract:
     - `NP_RMPI_ENABLE_DENSITY_INID_TEST=1`.
6. Commit:
   - `np-npRmpi`: `896ca4b` (staging + guards), `91132b8` (serial-`np` routing for opt-in path).
7. Validation:
   - `/tmp/nprmpi_session_routing_contract_20260224_101550.log` (`PASS 20, FAIL 0, SKIP 2`)
   - `/tmp/nprmpi_issue_notes_repros_20260224_102108.log` (`all verified issue-note repros passed`)

## GLP Gradient Mask Alignment Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Updated `gradients.npregression(...)` GLP masking logic to preserve any requested derivative order up to polynomial degree (instead of first-order-only masking).
2. Removed stale warning text claiming higher-order GLP derivatives were unavailable.
3. Scope:
   - `R/regression.R`
4. Commit:
   - `np-npRmpi`: `8be3fc1`
5. Validation:
   - direct session smoke: `/tmp/nprmpi_glp_gradmask_repro_20260224.out` (`NPRMPI_GLP_GRADMASK_OK`)
   - issue-note sweep: `/tmp/nprmpi_issue_notes_repros_20260224_afterpatch.log` (`all verified issue-note repros passed`)

## Plot Bootstrap Wild-Selector Wiring Checkpoint (2026-02-24)
Completed in `np-npRmpi`:
1. Wired `plot.errors.boot.wild` through MPI plot engines to fully propagate draw-selector choice into bootstrap helper calls:
   - `plot.rbandwidth`, `plot.scbandwidth`, `plot.plbandwidth`, `plot.sibandwidth`.
2. Preserved `npRmpi`-specific session/autodispatch guards and rank behavior while adding selector propagation.
3. Scope:
   - `R/np.plot.engine.rbandwidth.R`
   - `R/np.plot.engine.scbandwidth.R`
   - `R/np.plot.engine.plbandwidth.R`
   - `R/np.plot.engine.sibandwidth.R`
4. Commit:
   - `np-npRmpi`: `bada05c`
5. Validation:
   - user-style session repro: `/tmp/repro_nprmpi_plot_hang_check_20260224_en0_afterpatch.out` (`NPRMPI_PLOT_REPRO_OK`)
   - small session repro: `/tmp/nprmpi_plot_repro_small_20260224.out` (`NPRMPI_PLOT_SMALL_OK`)
   - issue-note sweep: `/tmp/nprmpi_issue_notes_repros_20260224_afterpatch.log` (`all verified issue-note repros passed`)
   - tarball build: `/tmp/nprmpi_build_20260224_checkpoint.log` (`creating vignettes ... OK`)

## `npRmpi` Ksum Fast-Path Default-On Promotion Probe (2026-02-24)
Status: deferred (no code change).
1. Ran bounded session-mode probes to evaluate promoting `np.plot.inid.ksum.fastpath.nprmpi` from opt-in to default-on.
2. In this sandbox/runtime, method-level probes did not complete within practical time bounds for all density-family cases, so promotion was not enabled in this checkpoint.
3. Safety decision: keep current default-off gate, retain fallback behavior, and require explicit opt-in until bounded reproducible performance/parity evidence is complete.
