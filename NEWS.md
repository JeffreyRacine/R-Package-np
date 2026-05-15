# npRmpi 0.70-2

- `npqreg()` is now a fully fledged MPI-aware quantile-regression front
  end. It supports the formula/data workflow, internally computes
  `npcdistbw()` bandwidths when a bandwidth object is not supplied,
  accepts scalar or vector `tau`, reuses selected bandwidths for
  additional quantiles in `plot()`, and exposes the usual S3 surface:
  `fitted()`, `predict()`, `predict(..., se.fit=TRUE)`, `se()`,
  `gradients()`, `summary()`, `print()`, `quantile()`, and `plot()`.
- `npqreg()` prediction now honors the standard `newdata` workflow while
  preserving native `exdat` precedence for compatibility with existing
  `npRmpi` call surfaces. Formula-based prediction validates that new
  data contain the required right-hand-side variables.
- `npqreg()` plotting has been expanded for vector quantiles,
  level/gradient displays, ordered predictors, user-specified legends,
  and object-fed plotting of additional `tau` values without recomputing
  cross-validation. The fixed-bandwidth gradient path now uses the
  MPI-aware helper route.
- `npconmode()` is now a first-class conditional-mode estimator. It
  supports formula/data and bandwidth-object workflows, forwards
  bandwidth-selection options to `npcdensbw()`, propagates local
  polynomial and NOMAD metadata, and exposes `fitted()`, `predict()`,
  `summary()`, `print()`, `gradients()`, and `plot()` methods.
- `npconmode()` now supports optional class-probability matrices and
  level-specific probability gradients. For non-local-constant fits,
  probabilities are normalized to be non-negative and to sum to one
  across the discrete response support before modal classification.
- `npconmode()` now fails early for non-categorical responses and
  validates formula-based `newdata` against the original right-hand-side
  variables.
- `npconmode()` plotting now supports object-fed class-probability slices
  and two-dimensional probability surfaces, optional `rgl` rendering, and
  probability-level asymptotic intervals where defined. Surface bootstrap
  intervals for class probabilities remain intentionally deferred.
- `npcopula()` is now a first-class copula estimator. It supports
  formula/data and bandwidth-object workflows, automatic two-dimensional
  probability grids, explicit `u` evaluation grids, and ordinary
  extractable object components including `$bws`.
- `npcopula()` now provides `fitted()`, `predict()`, `predict(...,
  se.fit=TRUE)`, `se()`, `summary()`, `print()`, `as.data.frame()`, and
  richer `plot()` methods. Plotting supports base `persp`, `image`, and
  optional `rgl` rendering, with asymptotic and MPI-fanned bootstrap
  intervals for copula surfaces where defined.
- `npcopula()` explicit-grid evaluation now uses the direct estimator
  route, preserving numerical results while avoiding the severe runtime
  growth of the previous expanded-grid path when users request larger
  probability grids.
- The automatic local-polynomial NOMAD controls have been split into
  explicit restart toggles: `powell.remin` for Powell restarts and
  `nomad.remin` for the second NOMAD hot start. This preserves the
  Powell Numerical Recipes restart default while allowing NOMAD hot
  starts to be controlled separately.
- Deprecated legacy `remin` remains accepted by `npregbw()` and `npreg()`
  with a warning and is mapped to the modern `powell.remin`/`nomad.remin`
  controls where appropriate, preserving downstream compatibility while
  documenting the new spelling.
- Hat-operator helpers now support an additional constraint-oriented
  output route for objects needed by shape-constrained quadratic
  programming workflows, avoiding reimplementation of local-polynomial
  hat-matrix construction in user examples.
- Local-polynomial derivative support has been broadened across the
  conditional estimator family. `npreg()`, `npcdens()`, and `npcdist()`
  now honor `gradient.order` more consistently for fitted, evaluated,
  predicted, and plotted objects when the selected polynomial degree is
  high enough, including vector derivative orders over continuous
  predictors and tensor/additive/Bernstein local-polynomial bases. The
  MPI implementation dispatches the corresponding conditional hat-apply
  helper work across the active worker pool where applicable.
- Core and semiparametric S3 prediction paths have been hardened around
  `newdata`, native evaluation-argument precedence, formula RHS
  validation, and `se.fit` handling while preserving `npRmpi` route
  independence.
- Front-end/bandwidth argument hygiene has been tightened so
  estimator-only controls such as `proper` are not forwarded into
  bandwidth selectors that do not accept them.
- MPI lifecycle and plotting routes received additional hardening,
  including soft `npRmpi.quit()` behavior, local object-fed plot
  computation where required, and explicit fanout of applicable
  bootstrap workloads.
- Documentation has been refreshed for the promoted `npqreg()`,
  `npconmode()`, and `npcopula()` workflows, including the
  local-polynomial NOMAD route, probability/gradient outputs, plot
  controls, and examples that use the streamlined interfaces.
- The pre-release validation suite was expanded with focused hostile
  argument tests, S3 contract tests, installed/tarball proof scripts,
  route-aware MPI probes, and serial/MPI parity checks for the newly
  promoted estimator families.

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
