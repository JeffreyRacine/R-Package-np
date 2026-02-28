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
