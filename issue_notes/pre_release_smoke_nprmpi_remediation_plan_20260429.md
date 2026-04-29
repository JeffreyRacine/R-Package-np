# npRmpi Pre-Release Smoke Remediation Plan

Date: 2026-04-29

## Post-Repair Status

Resolved in the live tree on 2026-04-29.

Two issues were repaired:

1. The native estimator-global alias cleanup from `np` was surgically ported to
   `npRmpi`.
2. The full test harness now exits nonzero when `testthat` reports failures,
   and stale full-lane expectations were reconciled with current installed
   check context, progress output, and scale-factor search policy.

Validation artifacts:

- `/Users/jracine/Development/tmp/pre_release_np_abort_20260429/logs/nprmpi_alias_cleanup_sentinel_live.log`
- `/Users/jracine/Development/tmp/pre_release_np_abort_20260429/logs/nprmpi_reconciled_targeted_tests.log`
- `/Users/jracine/Development/tmp/pre_release_np_abort_20260429/logs/nprmpi_full_lane_reconciled_check.log`
- `/Users/jracine/Development/tmp/pre_release_np_abort_20260429/logs/nprmpi_default_as_cran_reconciled_check.log`

The reconciled full lane is now green under `NP_CHECK_FULL=1`.

## Scope

This plan covers issues found while exercising the submitted
`/Users/jracine/Development/npRmpi_0.70-1.tar.gz` tarball through pre-release
smoke surfaces.

The primary release and MPI smoke lanes are green:

- local `R CMD check --as-cran`: `1 WARNING, 1 NOTE`
- win-builder R-release: `1 NOTE`
- win-builder R-devel: `1 NOTE`
- `NP_CHECK_FULL=1 R CMD check --no-manual --no-vignettes --ignore-vignettes`:
  top-level check status `OK`
- one-off MPI benchmark smoke, `n = 50`, `times = 1`, `nslaves = 1`: all 11
  functions passed
- copied demo matrix with `NP_DEMO_N = 100`: completed successfully across
  serial, attach, profile, and `NP = 2, 3, 4`
- generated 175 demo `.Rout` files; scan found no `Error`, `Execution halted`,
  `MPI_Abort`, timeout, or fatal tokens
- no stale `slavedaemon.R`, `Rslaves.sh`, demo, or `mpiexec` processes remained
  after the run

Primary artifacts:

- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/logs/nprmpi_full_check_no_vignettes.log`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/npRmpi.Rcheck`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/results/nprmpi_oneoff_manifest.csv`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/logs/nprmpi_demo_runall_n100.log`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/nprmpi_demo`

## Issue NR-1: Full Test Output Contains Failures Despite Top-Level OK

Observed command:

```sh
NP_CHECK_FULL=1 R CMD check --no-manual --no-vignettes --ignore-vignettes \
  /Users/jracine/Development/npRmpi_0.70-1.tar.gz
```

Observed result:

- top-level `R CMD check` status is `OK`;
- however, `npRmpi.Rcheck/tests/testthat.Rout` contains a `Failed tests`
  section with `FAIL 13`;
- failures include:
  - source-file reading tests unavailable in installed check context,
  - a source-test direct call to unqualified internal helper
    `updateBwNameMetadata`,
  - NOMAD progress expectation drift for several MPI families,
  - `npindexbw` and `npscoefbw` start-control expectations that still refer to
    legacy defaults rather than the current scale-factor search policy.

This is primarily a test-lane integrity issue: the check harness must not allow
visible test failures to coexist with a top-level `OK` status.

## Highest-Standards Repair Strategy

1. First repair the test harness semantics, not package runtime behavior:
   - make `NP_CHECK_FULL=1` fail hard when any test fails;
   - confirm `R CMD check` returns nonzero when a deliberate sentinel failure is
     inserted in a detached scratch tree;
   - preserve normal default-lane behavior.
2. Split full-lane tests by execution context:
   - tests that read source files must either run only from a source-tree
     context with explicit paths, or skip clearly in installed check context;
   - tests that inspect internals must use `getFromNamespace()` or helper setup
     explicitly rather than relying on unqualified names from source loading.
3. Reassess the reported behavioral expectations after the harness is honest:
   - NOMAD progress-output tests should reflect the current intended MPI
     progress contract; if the output is wrong, repair progress I/O at the
     route; if the test is obsolete, update the test expectation with a short
     rationale.
   - `npindexbw` and `npscoefbw` start-control tests must be reconciled with
     the accepted global `scale.factor.*` policy; update tests only after
     verifying runtime behavior matches the current scientific contract.
4. Do not mark failing full-lane tests as skipped merely to get `OK`.
   Every skip must document why the route cannot be evaluated in that context
   and where the same contract is covered instead.
5. Validation after candidate repair:
   - `R CMD check --no-manual --no-vignettes --ignore-vignettes` default lane,
   - `NP_CHECK_FULL=1 R CMD check --no-manual --no-vignettes --ignore-vignettes`,
   - one-off MPI benchmark smoke with `nslaves=1`,
   - demo matrix with `NP_DEMO_N=100`,
   - stale-worker process scan after completion.

## Non-Issues From This Pass

- The win-builder install-log `WARNING: this package has a configure script`
  is informational and not a source compiler warning.
- The MPI demo matrix itself is healthy under the documented smoke settings.

## Release Risk Classification

Medium.

Runtime release surfaces are green, but the extended test lane has governance
problems. Before final release, the full lane should either be made honest and
green or explicitly moved out of the release-hardening gate with a documented
replacement coverage plan.
