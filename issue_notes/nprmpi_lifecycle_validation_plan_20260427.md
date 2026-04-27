# npRmpi Lifecycle Validation Plan (2026-04-27)

## Goal

Validate the lifecycle-hardening changes in `npRmpi` with the highest chance of
catching regressions and the lowest chance of collateral package changes.

This plan covers commit `e4f5dc2e` and its immediate successors if validation
requires a narrow follow-up. It does not authorize estimator, bandwidth-search,
or compiled numerical changes unless a reproducible production-code root cause
is identified.

## Scope Lock

Validate only these change classes:

1. Production runtime:
   - `.onUnload()` finalizes MPI when needed but no longer explicitly calls
     `library.dynam.unload("npRmpi", ...)`.
2. Route validation harness:
   - `validate_route_attach.R` uses the canonical attach-mode init return value
     to gate master work.
   - `validate_route_attach.R` uses explicit `data=` for the `npregbw()` /
     `npreg()` smoke, matching shipped attach-demo style.
   - `run_attach_profile_validations.sh` installs into a retained private
     library under its output directory, runs attach/profile/profile-plot before
     manual-broadcast, and disables spawned-slave reuse for the one-off manual
     route.

Out of scope unless a later root-cause proof requires it:

- estimator semantics,
- bandwidth search logic,
- C-level objective or collective code,
- public API defaults,
- `np-master`.

## Risk Hypotheses

1. The `.onUnload()` change could improve unload/finalize stability while
   accidentally changing load/unload behavior on this host or CRAN-like systems.
2. The route validator could pass because the harness was softened rather than
   because the underlying route remains healthy.
3. MPI state contamination from back-to-back launcher modes could be mistaken
   for an estimator defect.
4. Direct attached demos could remain healthy while a bundled validator fails,
   or vice versa; both surfaces must be distinguished.
5. A stress-style runner that repeatedly mixes spawn, attach, and profile
   launches can expose host MPI fragility without proving a package regression.
   Such evidence is useful, but it must not override cleaner direct-route
   evidence unless the same failure reproduces in a canonical shipped route.

## Clean-Field Preflight

Before every route batch:

1. Record commit SHA:
   - `git -C /Users/jracine/Development/np-npRmpi rev-parse --short HEAD`
2. Record worktree state:
   - `git -C /Users/jracine/Development/np-npRmpi status --short`
3. Confirm no stale MPI processes:
   - scan for `mpiexec`, `hydra_pmi`, `slavedaemon.R`, `Rslaves.sh`,
     `validate_route`, and `npRmpi`.
4. Save all logs under:
   - `/Users/jracine/Development/tmp/nprmpi_lifecycle_validation_20260427`
5. Use fresh output subdirectories for each route batch.
6. Record the validation mode for each artifact:
   - installed/live library,
   - retained private library,
   - tarball-first check,
   - copied demo harness.
7. For profile/manual-broadcast launches, set `R_LIBS` at process start; do not
   rely on late `.libPaths()` changes after `inst/Rprofile` has loaded
   `npRmpi`.

If stale route processes are found, stop and identify their owner before
killing anything.

Do not run overlapping MPI route batches. If Jeffrey is running demos or local
MPI probes, wait for those to finish or coordinate a clean validation window.

## Gate 1: Production Runtime

Purpose: validate the actual runtime change in `R/zzz.R`.

Run:

1. Install the current tree:
   - `R CMD INSTALL --preclean --no-multiarch --with-keep.source .`
2. Repeated load/unload subprocess loop:
   - repeatedly run `library(npRmpi); unloadNamespace("npRmpi")`
   - capture warnings, errors, and exit status.
3. Tarball-first check:
   - from `/Users/jracine/Development`,
   - `R CMD build np-npRmpi`,
   - `R CMD check --no-manual --ignore-vignettes npRmpi_0.70-1.tar.gz`.

Pass criteria:

- load/unload loop exits cleanly,
- `R CMD check` reports `Status: OK`,
- no new unload or MPI-finalize warnings appear.
- `checking whether the package can be unloaded cleanly` remains `OK`.

Rollback trigger:

- if this gate fails in a way attributable to `.onUnload()`, revert only the
  `.onUnload()` production change and rerun the gate before touching anything
  else.

## Gate 2: Canonical Route Validators

Purpose: validate lifecycle route surfaces without relying on demos alone.

Run the direct route validators first:

1. direct attach validator slices:
   - run `issue_notes/validate_route_attach.R` in a fresh `mpiexec` world for
     each `NP_RMPI_ATTACH_VALIDATE_ROUTE` slice:
     `npreg`, `npscoef`, `npplreg`, `npindex`, `npcopula`, `npconmode`;
   - do not accumulate all attach-route families in one attach world, because
     shipped attach demos are one route per process and repeated mixed-family
     attach worlds are a stress probe rather than a user-facing contract.
2. direct profile validator:
   - launch with exactly one profile source, using `inst/Rprofile`.
3. direct manual-broadcast validator:
   - run with `NP_RMPI_NO_REUSE_SLAVES=1` because it is a one-off subprocess
     route gate, not an interactive reused-slave session.

Then run the full validator:

4. full validator:
   - `issue_notes/run_attach_profile_validations.sh <fresh-output-dir>`
   - run once as a hard acceptance gate.
   - a second immediate repeat is a stress probe, not a hard acceptance gate,
     unless the same failure also appears in a direct route or shipped demo.

Pass criteria:

- `ATTACH_ROUTE_OK`,
- `ATTACH_NPREG_ROUTE_OK`,
- `ATTACH_NPSCOEF_ROUTE_OK`,
- `ATTACH_NPPLREG_ROUTE_OK`,
- `ATTACH_NPINDEX_ROUTE_OK`,
- `ATTACH_NPCOPULA_ROUTE_OK`,
- `ATTACH_NPCONMODE_ROUTE_OK`,
- `PROFILE_ROUTE_OK`,
- `PROFILE_NPCOPULA_ROUTE_OK`,
- `PROFILE_NPCONMODE_ROUTE_OK`,
- `PROFILE_PLOT_ROUTE_OK`,
- `MANUAL_BCAST_ROUTE_OK`,
- no orphan MPI processes afterward.

Failure triage:

- if direct attach/profile pass but full validator fails, treat the runner as
  suspect first;
- if direct route validators fail, compare against shipped demos before
  touching package code;
- if a multi-family attach validator fails but sliced attach validators and
  shipped one-route demos pass, classify the multi-family run as a stress
  harness, not as production-route evidence;
- if a direct route fails only after a previous mixed-mode failure, reset the
  MPI field and reproduce from a clean-field preflight before classifying it;
- do not add sleeps, broad retries, or estimator fallbacks as fixes.
- do not alter production attach/session/profile code unless a canonical
  shipped route fails reproducibly from a clean field.

## Gate 3: Demo Triplet Sentinels

Purpose: ensure shipped user-facing demo launch patterns still work.

Use a copied demo harness under the validation artifact root so tracked demo
outputs are not overwritten.

Run at `NP_DEMO_N=100`, `NP=2`:

1. `npreglcls`
   - sentinel for the previously fragile attach route surface.
2. `npcdensml`
   - sentinel for recent `npcdens` ML / NOMAD MPI work.
3. `npcdensls`
   - sentinel for recent `npcdens` LS / bounded-quadrature MPI work.

For each sentinel, run:

- serial,
- attach,
- profile/manual-broadcast.

Run from a copied demo tree under the artifact root. Preserve the canonical
layout (`makefile`, scripts, `serial`, `n_2_attach`, `n_2_profile`) and set an
explicit absolute `RPROFILE` for copied profile runs.

Pass criteria:

- each demo command exits zero,
- each `.Rout` has no `Error:`, `Execution halted`, `Abort(`, or demo
  `FAILED` marker,
- no orphan MPI processes afterward.

## Gate 4: Negative Controls And Claim Discipline

Retain prior failure logs showing:

- temporary-library runner failures,
- formula/no-`data=` attach validator failures,
- back-to-back mixed-mode failures.

These are evidence for why the harness was changed. They are not acceptance
evidence for production estimator changes.

Claim rules:

- call the validator changes "harness hardening";
- call the `.onUnload()` change "production runtime lifecycle hardening";
- do not claim estimator, bandwidth-search, or C collective repair from this
  tranche;
- report any remaining MPICH/OpenMPI instability as route-launcher evidence
  unless a direct shipped route reproduces deterministically.
- distinguish accepted gates from stress probes in the closeout. A stress probe
  failure can motivate a separate plan, but it cannot by itself invalidate a
  clean direct-route and demo-triplet acceptance set.

## Gate 5: CRAN-Specific Follow-Up

Purpose: keep the CRAN warning context separate from local MPICH route stress.

If local gates pass, preserve a short CRAN-focused note with:

1. The exact `.onUnload()` delta.
2. Hao's stated rationale: avoid explicit dynamic-library unload racing MPI.
3. Local evidence:
   - unload-clean check is `OK`,
   - namespace unload check is `OK`,
   - tarball-first `R CMD check --no-manual --ignore-vignettes` is `Status: OK`.
4. Any remaining limitation:
   - local MPICH version differs from CRAN's upgraded stack;
   - local evidence is supportive, not a substitute for CRAN's Debian check.

## Final Acceptance

This tranche is acceptable only if:

1. Gate 1 passes.
2. Gate 2 passes.
3. Gate 3 passes.
4. Gate 5 note is prepared if CRAN communication is needed.
5. No stale MPI processes remain.
6. The git diff remains limited to lifecycle/harness files.
7. All artifact roots and logs are named in the final closeout.

If any hard gate fails, stop and classify the failure before adding more code.

## Non-Goals

Do not:

- tune MPI provider settings as a package fix;
- change estimator behavior to placate a harness;
- add timing sleeps as production synchronization;
- broaden route validators until the current lifecycle surface is stable;
- port anything to `np-master`, which is out of scope.
