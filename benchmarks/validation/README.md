# Validation Harnesses

Use this folder for numerical-parity checks, option probes, and
smoke-validation scripts.

Scripts in this folder should prioritize coverage and diagnostics over timing.

## `nmulti` Default Study

- `nmulti_default_study.R`
  - paired simulation harness for investigating whether smaller
    `nmulti` defaults are scientifically defensible,
  - compares runtime, bandwidth drift, evaluation-surface drift, and
    truth-based error across matched seeds and scenarios,
  - writes raw runs, evaluation-point outputs, summaries, and a markdown
    report to a chosen `/tmp` directory.

- `nmulti_default_study_README.md`
  - reproduction notes and interpretation guidance for the `nmulti`
    study.

- `nmulti_core_surface_study.R`
  - focused default-choice harness for `npreg` and `npcdens`,
  - spans `lc`, `ll`, and `lp` (`degree=1,2`) with one continuous and
    one categorical regressor,
  - records runtime, objective diagnostics, and fitted-surface drift
    relative to a larger `nmulti` reference.

- `nmulti_core_surface_study_README.md`
  - reproduction notes and design rationale for the focused core study.

- `nmulti_wage1_start_stress.R`
  - real-data robustness study for `npregbw()` on `wage1`,
  - varies LP starting bandwidth vectors across categorical and
  continuous dimensions,
  - compares `nmulti = 1, 2, 5` using `fval` and fitted-surface drift
  relative to the same-start `nmulti=5` run.

- `nmulti_borderline_density_dist_study.R`
  - targeted synthetic follow-up for borderline `npcdens`/`npcdist`
    routes,
  - restricts attention to `ll` and `lp(degree=2)` under two
    mixed-data DGPs,
  - compares `nmulti = 1, 2, 5` using truth-based surface error and
    drift relative to `nmulti=5`.

- `nmulti_wage1_nomad_stress.R`
  - real-data start-value stress test for automatic LP degree search on
    `wage1`,
  - uses `search.engine='nomad+powell'` with user-supplied bandwidth
    starts,
  - compares `nmulti = 1, 2, 5` on selected degree, `fval`, and fitted
    drift relative to `nmulti=5`.

- `nmulti_default_cap_summary.md`
  - closeout note for the empirical default-choice decision,
  - summarizes the synthetic, real-data, and NOMAD evidence supporting
    the new `nmulti` cap of `min(2, p)`.

## Hotspot Remediation Baseline

- `hotspot_wp0_wp1_pre.R`
  - short-run pre-baseline harness for remediation WP0/WP1,
  - captures elapsed timing, `num.feval`, objective metadata, bandwidth signatures,
    output signatures, and GC peak memory proxies,
  - writes raw + summary + pre/post scaffold CSV artifacts to `/tmp` by default.

- `hotspot_wp2_compare_prepost.R`
  - strict pre/post comparator for WP2 edits,
  - enforces CV path-lock checks on `bw` scenarios (`num.feval`, objective history, final `fval`, bandwidth signature),
  - emits scenario-level parity classes: `PASS_STRICT`, `PATH_DIVERGENCE`, `FAIL`,
  - writes row-level and summary CSV outputs to a chosen `/tmp` directory.

## Triage First Steps

- Confirm argument parity between compared runs.
- Check seed policy (`fixed` vs `varying`) and keep it explicit.
- Verify outputs include enough diagnostics (objective values, selected
  bandwidths, and max-abs diffs when relevant).

## Canonical Implementation Directive (2026-03-05)

This repository follows a strict canonical execution rule:

1. One canonical implementation per method (outside explicit `np.tree` branching).
2. Unsupported configurations must fail fast with explicit `stop(...)` diagnostics.
3. No silent remap/coercion of user-selected options (for example `bwmethod`, `regtype`, kernels, `cv.iterate`, or bounds transforms).
4. No hidden alternate execution paths for the same method semantics.
5. All fit-defining options (for example `degree`, `basis`, `bernstein.basis`, kernels, and bounds) must be propagated and used by the canonical path.
6. `np.tree=FALSE` is the default; when `np.tree=TRUE`, behavior must remain semantics-preserving and option-compatible with the canonical path.
7. Remove or reject legacy/debug compatibility branches that add redundant runtime overhead once canonical behavior exists.
