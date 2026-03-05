# Validation Harnesses

Use this folder for numerical-parity checks, option probes, and
smoke-validation scripts.

Scripts in this folder should prioritize coverage and diagnostics over timing.

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
