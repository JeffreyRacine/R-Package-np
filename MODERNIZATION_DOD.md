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
