# Validation Harnesses

Use this folder for numerical-parity checks, option probes, and
smoke-validation scripts.

Scripts in this folder should prioritize coverage and diagnostics over timing.

## Triage First Steps

- Confirm argument parity between compared runs.
- Check seed policy (`fixed` vs `varying`) and keep it explicit.
- Verify outputs include enough diagnostics (objective values, selected
  bandwidths, and max-abs diffs when relevant).
- For MPI failures, inspect startup/interface settings before method logic.
