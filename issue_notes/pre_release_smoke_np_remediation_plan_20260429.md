# np Pre-Release Smoke Remediation Plan

Date: 2026-04-29

## Post-Repair Status

Resolved in the live tree on 2026-04-29.

Root cause was stale native estimator-global alias pointers left non-NULL after
bandwidth-selector cleanup. Later shadow/reset proof paths treated those stale
aliases as live owned state and could free invalid memory. The repair is a
narrow native cleanup helper that clears estimator-global aliases after the
affected bandwidth/estimator routes release their borrowed state.

Validation artifacts:

- `/Users/jracine/Development/tmp/pre_release_np_abort_20260429/logs/repro_shadow_all_large_pair_live.log`
- `/Users/jracine/Development/tmp/pre_release_np_abort_20260429/logs/test_file_shadow_live.log`
- `/Users/jracine/Development/tmp/pre_release_np_abort_20260429/logs/np_default_as_cran_live_tarball_check.log`

The full `NP_CHECK_FULL=1` lane no longer aborts, but it exposes unrelated
full-only stale test expectations. Those are not part of this native repair.

## Scope

This plan covers issues found while exercising the submitted
`/Users/jracine/Development/np_0.70-1.tar.gz` tarball through pre-release smoke
surfaces.

The default release lanes are green:

- local `R CMD check --as-cran`: `1 WARNING, 1 NOTE`
- win-builder R-release: `1 NOTE`
- win-builder R-devel: `1 NOTE`
- default source smoke filter `check-core-smoke`: passed
- demos `Engel95`, `npregiv`, and `tree`: passed
- one-off benchmark smoke, `n = 50`, `times = 1`: all 11 functions passed
- tiny serial synthetic microbenchmark smoke: passed

Primary artifacts:

- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/logs/np_full_check_no_vignettes.log`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/np.Rcheck`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/test_logs/np_source_full_testthat.log`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/results/np_oneoff_manifest.csv`
- `/Users/jracine/Development/tmp/pre_release_smoke_20260429/results/np_serial_microbench_synth`

## Issue NP-1: Full Extended Test Lane Aborts In Native Shadow Proof Surface

Observed command:

```sh
NP_CHECK_FULL=1 R CMD check --no-manual --no-vignettes --ignore-vignettes \
  /Users/jracine/Development/np_0.70-1.tar.gz
```

Observed result:

- `Status: 1 ERROR`
- `tests/testthat.R` aborts with `Abort trap: 6`
- output stops while running the native conditional-density shadow proof
  surface, after `npc-cv-shadow-proof-contract` progress and repeated
  degree-search progress lines.

This is not a CRAN/default-lane failure, but it is a real pre-release hardening
signal because the extended lane should be runnable and interpretable.

## Highest-Standards Repair Strategy

1. Reproduce in a detached worktree or scratch source tree, never in the live
   release tree first.
2. Build a minimal failing reproducer from
   `tests/testthat/test-npc-cv-shadow-proof-contract.R`:
   - run file sections in order,
   - identify the first `test_that()` block that aborts,
   - reduce that block to the smallest `.Call()`/R wrapper sequence that still
     aborts.
3. Capture native failure evidence before editing:
   - exact `.Call` symbol,
   - bandwidth object metadata,
   - `regtype`, `degree`, `basis`, `bwtype`,
   - tree flag,
   - criterion,
   - shadow reset state before/after,
   - backtrace from `lldb` or equivalent.
4. Inspect the touched native helper and its state lifecycle:
   - PROTECT/UNPROTECT balance,
   - allocation sizes and zero-column matrix handling,
   - shadow-state reset guarantees,
   - ownership of returned SEXP objects,
   - assumptions about degree vector length and kernel-code layout.
5. Decide whether the abort exposes:
   - a real native helper defect, which must be repaired in C with a focused
     contract test; or
   - a test-only misuse of an internal native diagnostic helper, which must be
     repaired in the test harness without weakening production contracts.
6. Do not suppress, skip, or remove the shadow proof to obtain a green lane
   unless root-cause evidence proves the test is invalid and a smaller valid
   replacement preserves the same contract.
7. After any candidate repair, validate in this order:
   - minimal reproducer,
   - full `test-npc-cv-shadow-proof-contract.R`,
   - `NP_CHECK_FULL=1 R CMD check --no-manual --no-vignettes --ignore-vignettes`,
   - default `R CMD check --as-cran`,
   - a focused installed-build sentinel using `npcdensbw()` / `npcdistbw()`
     LP shadow-adjacent routes.

## Non-Issues From This Pass

- The `PXTKX` compiler warning is retained as accepted false-positive context.
- The local `checkbashisms` warning is host tooling, not package behavior.
- The microbenchmark comparator script requires pre/post directories; invoking
  it without those arguments was a harness misuse, not a package issue.

## Release Risk Classification

Medium.

The default CRAN and win-builder lanes are green, but an extended native proof
lane aborts. This should be repaired or explicitly dispositioned before final
release if time permits, because native aborts are never acceptable as lingering
unknowns.
