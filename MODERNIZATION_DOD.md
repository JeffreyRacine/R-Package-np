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
- [ ] Full `--as-cran` warning/note closure still required (accepted-warning ledger now tracked in `AS_CRAN_WARNING_LEDGER.md`).
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
