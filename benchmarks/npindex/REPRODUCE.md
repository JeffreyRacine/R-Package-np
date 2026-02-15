# Reproduce: npindex (npRmpi)

1. Build old/current `np` + `npRmpi` libs.
2. Run `run_npindex_combos.R` for each library with identical `n/times/base_seed/nmulti/nslaves`.
3. Compare raw outputs with `compare_npindex_versions.R`.
4. Render markdown via `make_npindex_report.R`.
