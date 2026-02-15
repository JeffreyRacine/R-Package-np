# Reproduce: npscoef (npRmpi)

1. Build old/current `np` + `npRmpi` libs.
2. Run `run_npscoef_combos.R` for each library with identical `n/times/base_seed/nmulti/nslaves`.
3. Compare raw outputs with `compare_npscoef_versions.R`.
4. Render markdown via `make_npscoef_report.R`.
