# `nmulti` Default-Cap Decision (Mirror)

`np-master` is the source of truth for the empirical `nmulti`
default-cap study. The serial study harnesses live in:

- `/Users/jracine/Development/np-master/benchmarks/validation/nmulti_core_surface_study.R`
- `/Users/jracine/Development/np-master/benchmarks/validation/nmulti_borderline_density_dist_study.R`
- `/Users/jracine/Development/np-master/benchmarks/validation/nmulti_wage1_start_stress.R`
- `/Users/jracine/Development/np-master/benchmarks/validation/nmulti_wage1_nomad_stress.R`
- `/Users/jracine/Development/np-master/benchmarks/validation/nmulti_default_cap_summary.md`

Mirrored conclusion for `npRmpi`:

- adopt the same default multistart cap of `min(2, p)` for the core
  bandwidth selectors
- keep `nmulti=1` for the univariate boundary density helper
- mirror the same NOMAD fallback behavior when `nmulti` is not supplied

The point of this note is to keep the rationale visible inside the
`npRmpi` benchmark tree without maintaining duplicate empirical harness
files in two repos.
