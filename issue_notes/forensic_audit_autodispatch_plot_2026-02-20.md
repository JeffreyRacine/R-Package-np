# Forensic Audit: Autodispatch + S3 Plot Stack (2026-02-20)

## Scope
- Repo: `/Users/jracine/Development/np-npRmpi`
- Focus: `npRmpi.autodispatch`, S3 plotting migration, bootstrap plot chain, legacy MPI wrapping assumptions.

## Confirmed Architecture State
1. S3 methods are registered for bandwidth and fitted object classes (`plot.*` in `NAMESPACE`).
2. `plot.*` methods now route through `R/np.plot.methods.R` and then into specialized `npplot.*` backends.
3. The monolithic legacy entry (`npplot`) still exists as compatibility orchestrator (not removed yet).
4. Bootstrap-capable plot paths are centralized via helpers in `R/np.plot.helpers.R`, including payload-based bootstrap execution.

## Findings

### A. GlobalEnv temp-object dependency in autodispatch (investigated; not yet removable)
- Previous behavior relied on assigning transient objects into `.GlobalEnv` before broadcasting.
- This created side effects and was fragile under nested/manual dispatch boundaries.
- Attempted refactor to private package env was rolled back after smoke failures:
  - `npregbw`/formula and call-environment paths require symbol visibility patterns tied to existing evaluation chains.
  - Regression observed: unresolved temporary symbols in downstream call reconstruction/evaluation.
  - Conclusion: `.GlobalEnv` transient staging is currently still required for compatibility; removal deferred.

### B. `eval(parse())` usage remains extensive (partially improved only)
- Replaced one broadcast-by-name parse path with direct call construction.
- Large residual `eval(parse())` blocks remain in legacy `npplot.*` backend files and several non-plot legacy files.
- These should be reduced incrementally, not in one wholesale patch, due behavior risk.

### C. Legacy plotting backend still mixed-mode
- S3 facade is in place.
- Backend still uses old dynamic-string invocation style in multiple files:
  - `R/np.plot.npplot.bandwidth.R`
  - `R/np.plot.npplot.rbandwidth.R`
  - `R/np.plot.npplot.conbandwidth.R`
  - `R/np.plot.npplot.condbandwidth.R`
  - `R/np.plot.npplot.dbandwidth.R`
  - `R/np.plot.npplot.plbandwidth.R`
  - `R/np.plot.npplot.scbandwidth.R`
  - `R/np.plot.npplot.sibandwidth.R`
- Result: modular entry points are improved, but internals are still legacy-heavy.

### D. Safe removals identified
1. Dead helper removed in this patch:
   - `.npRmpi_bcast_robj_by_name()`
2. Additional safe-removal candidates (pending verification in next pass):
   - stale/duplicated plotting argument munging fragments once specific `eval(parse())` regions are rewritten.

## Production-Risk Notes
1. Plot/bootstrap pathways are now operational with payload strategy but still sensitive to legacy internals in `npplot.*` backends.
2. Any broad refactor replacing string-built expressions with `do.call` should be done per backend family with smoke tests after each file.
3. Keep `npplot` generic until each `npplot.*` backend is proven equivalent under both autodispatch and manual MPI contexts.

## Next Patch Queue (ordered)
1. Rebuild/install `npRmpi` after current autodispatch private-env patch.
2. Smoke tests:
   - `npregbw` + `npreg` + `gradients`
   - `npsigtest`
   - `plot(..., plot.errors.method="bootstrap")` for `npreg`, `npcdens`, `npcdist`
3. Run `man/run` then `R CMD build` + `R CMD check --as-cran` for `npRmpi`.
4. Start staged `eval(parse())` retirement in plotting internals:
   - first target `R/np.plot.npplot.bandwidth.R` and `R/np.plot.npplot.rbandwidth.R`.
5. After parity confirmation, mirror compatible plot modularization cleanups into `np-master`.
