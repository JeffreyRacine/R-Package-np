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
- [x] Core modernization checkpoints validated with targeted contract tests + tarball checks.
- [ ] Remaining codoc/doc-signature drift in `man/*.Rd` must be reconciled.
- [ ] Full `--as-cran` warning/note closure and accepted-warning ledger still required.
- [ ] Win-builder validation still required before release candidate.

## Mandatory Release Gates

### 1) Interface and Semantics
- [ ] Public APIs for core estimators are stable (signature + return-structure contracts).
- [ ] Formula/default method pairs have parity tests (including `subset`/`na.action` behavior).
- [ ] S3 docs and method signatures match exactly (no codoc mismatches).

### 2) Evaluation and Call Construction
- [ ] No `eval(parse(...))` in R layer.
- [ ] High-risk `eval(...)` call construction paths are replaced with structured calls (`do.call`, explicit call objects) where appropriate.
- [ ] Any remaining `eval(...)` has a documented reason and a regression test.

### 3) Native Interface Safety
- [ ] `.C` callsites in R layer are `0`.
- [ ] `.Call` interface paths have stress tests for touched entry points.
- [ ] PROTECT/UNPROTECT accounting validated in modified C entry points.
- [ ] No new stack-imbalance warnings in targeted stress runs.

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
- [ ] Tarball-first validation is cleanly run:
  - `R CMD build np-master`
  - `R CMD check --as-cran np_<ver>.tar.gz`
- [ ] New warnings/notes are treated as regressions unless explicitly accepted and documented.
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
