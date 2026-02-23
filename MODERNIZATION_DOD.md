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
- [x] Verified issue-note repro harness includes session `npreg` factor-routing guard (`b5b597a`).
- [x] Load hook now supports source-tree/dev loading without installed-package lookup dependence (`61b4783`).
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
- [ ] High-risk `eval(...)` paths replaced where appropriate with structured calls.
- [ ] Any remaining `eval(...)` is documented and covered by tests.

### 4) Native Interface Safety
- [x] `.C` callsites in R layer are `0`.
- [ ] `.Call` interface paths have stress tests for touched entry points.
- [ ] PROTECT/UNPROTECT accounting validated for modified C entry points.
- [ ] No new stack-imbalance warnings in targeted serial + MPI runs.

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
