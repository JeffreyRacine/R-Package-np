# Phase 0 Design: `kernelcv.c` Follow-On Investigation

Date: 2026-03-07
Repo: `np-master`
Status: design only, no code change approved by this note

## Purpose

Define the smallest safe investigation plan for the next tranche after the
`statmods.c` safety fix. This note does not authorize any `kernelcv.c`
modification by itself. It exists to force predeclaration before touching a
timing-sensitive path.

## Scope Candidate

The only follow-on candidate currently worth investigating is:

1. repeated allocation in `cv_func_regression_categorical_ls_nn()`

Candidates explicitly not authorized by this note:

1. removing or weakening `check_valid_scale_factor_cv(...)`,
2. changing the public `timing` contract,
3. broad CV extern-state refactoring,
4. multi-candidate edits in one tranche.

## Change Classification

Blast-radius class:
- `P1` hot-path micro-optimization

Reason:
- this function is in the bandwidth-selection objective path;
- even a small internal change can alter timing and, if done poorly, numerical
  behavior.

## Hypothesis

Reusing the `mean` workspace in `cv_func_regression_categorical_ls_nn()` may
reduce repeated allocation/free overhead without changing:

1. selected bandwidth,
2. objective value,
3. invalid-parameter behavior,
4. public object fields,
5. `npRmpi` route behavior.

## Canary Set

Primary canary:
1. `npregbw` on continuous regression data with `bwtype="generalized_nn"`,
   `bwmethod="cv.ls"`, `nmulti=1`

Secondary canary:
1. `npregbw` on continuous regression data with `bwtype="adaptive_nn"`,
   `bwmethod="cv.ls"`, `nmulti=1`

Reason:
- these are the narrowest public consumers most likely to exercise the exact
  objective path under investigation.

## Beneficiary Family

Predeclared beneficiary family for this candidate:

1. `npregbw` continuous-data regression selectors using the NN CVLS path
2. serial `np`
3. `npRmpi` session/attach route smokes after serial approval

## Unaffected Family

Predeclared unaffected family for this candidate:

1. `npregbw` fixed-bandwidth CVLS regression selectors
2. regression selectors not using the NN CVLS objective
3. density/conditional-density/distribution selectors

This family is intentionally narrow enough to validate, but broad enough to
catch accidental spillover.

## Numerical Parity Contract

For any candidate patch, all of the following must hold:

1. exact equality where feasible for canary outputs,
2. otherwise explicit `all.equal`/max-abs-diff evidence,
3. bandwidth parity must be exact or explained with documented floating-point
   justification,
4. invalid-history and eval-history shape must not drift,
5. warnings/errors/messages must remain intentional.

## Timing Decision Design

This candidate is not a performance-claim tranche. It is a no-regression
micro-optimization investigation.

Required design:

1. interleaved paired-seed screening at `times=25`,
2. fixed-seed and varying-seed runs,
3. mean and median absolute and percent deltas,
4. tail and spread checks,
5. escalation to `times>=100` if the screening signal is not clearly favorable,
6. no keep decision on mixed or ambiguous evidence.

## MEI

Because this is a hot-path candidate, any accepted regression above noise is a
fail. For screening, use:

1. beneficiary superiority interest:
- at least `max(0.01 sec, 1.5%)` improvement in mean and median on the medium
  canary

2. unaffected-family equivalence margin:
- within `max(0.01 sec, 1.5%)` in absolute/relative terms

3. tail tolerance:
- no `p90` or `p95` regression beyond the unaffected-family margin

4. spread tolerance:
- no meaningful `sd`/`IQR`/`MAD` inflation beyond ordinary run-to-run noise

These are starting margins only; if pilot variance shows they are too loose or
too tight, revise them before the confirmatory tier and record the revision.

## Scripts and Artifacts

Preferred harnesses:

1. `/Users/jracine/Development/np-master/benchmarks/np_serial_microbench.R`
2. `/Users/jracine/Development/np-master/benchmarks/oneoff/run_oneoff_suite.R`
3. `/Users/jracine/Development/np-master/benchmarks/validation/hotspot_wp0_wp1_pre.R`
4. `/Users/jracine/Development/np-master/benchmarks/validation/hotspot_wp2_compare_prepost.R`

If none directly fit the candidate path, a tranche-specific `/tmp` harness is
allowed, provided it records:

1. exact invocation,
2. package/library path,
3. seed policy,
4. workload size,
5. raw timings,
6. summary statistics,
7. parity outputs.

## `npRmpi` Follow-On Requirement

No `npRmpi` port happens unless the `np-master` candidate survives:

1. serial parity,
2. serial performance gate,
3. installed-package targeted regression coverage.

After a port, required evidence is:

1. session/spawn tiny smoke,
2. attach tiny smoke,
3. profile/manual route if the touched family is represented there and the host
   environment allows it.

No `npRmpi` timing claim is required for this candidate unless the serial
results are strong enough to justify a larger follow-on study.

## Stop-Loss Rules

Drop the candidate immediately if any of the following occur:

1. bandwidth or objective parity drifts,
2. mean and median disagree in sign at screening tier,
3. affected and unaffected families disagree in direction,
4. tail or spread regress materially,
5. the gain is too small to clear the proof burden,
6. route behavior in `npRmpi` becomes less stable.

## Acceptable Outcomes

This note explicitly allows all of the following conclusions:

1. keep and implement the candidate,
2. inconclusive, do not implement,
3. valid concern but not worth patching,
4. defer until a stronger trigger appears.

The existence of a plausible micro-optimization is not enough to force a code
change.
