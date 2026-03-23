# npRmpi 0.70-1

- The default multistart cap for bandwidth selection now follows
  `min(2, p)` across the mirrored estimator families, replacing the
  older `min(5, p)` cap. This includes automatic LP degree-search calls
  when `search.engine="nomad"` or `"nomad+powell"` and `nmulti` is not
  supplied explicitly.
- The univariate boundary density helper `npuniden.boundary()` now
  defaults to `nmulti=1`.
- The empirical studies supporting this mirror change are documented in
  `np-master/benchmarks/validation/`, with a summary note kept in this
  repository's `benchmarks/validation/` folder.
- LP-capable front ends now accept `nomad=TRUE` as a documented
  convenience preset for the recommended automatic NOMAD
  local-polynomial route, mirroring the serial package defaults and
  help-page guidance.
