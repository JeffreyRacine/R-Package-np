# npRmpi 0.70-3

- Local-polynomial regression cross-validation now uses a leaner hot
  symmetric weighted-sum loop. Fixed-bandwidth `npregbw(..., regtype = "lp",
  bwmethod = "cv.ls")` objective probes in active MPI sessions match serial
  `np` objective values to numerical precision while substantially reducing
  local-polynomial CV evaluation time.
- Unconditional density least-squares cross-validation now uses a leaner
  fixed-bandwidth Gaussian convolution loop. Fixed-bandwidth
  `npudensbw(..., bwmethod = "cv.ls")` objective probes preserve objective
  values exactly in the focused validation rows while materially reducing the
  convolution portion of the objective calculation. Conditional-density
  least-squares objective probes inherit the same fixed-bandwidth Gaussian
  convolution improvement.
- Conditional density and conditional distribution least-squares
  cross-validation now use a size-aware row-block policy for local-polynomial
  objective evaluation. The accepted route keeps the bounded-quadrature cap
  unchanged, bounds transient memory by sample size, and preserves objective
  values to numerical precision while materially reducing evaluator overhead
  for fixed-bandwidth CVLS probes in serial and MPI sessions.
- Large-sample categorical-only regression now uses the MPI-safe
  profile-compressed route under `options(np.categorical.compress = TRUE)`,
  which is enabled by default. The legacy `options(np.tree = TRUE)` switch
  also enables the same route for compatibility. Repeated predictor profiles
  are compressed before bandwidth search, fitting, prediction/evaluation,
  standard errors, hat-helper use, and plot bootstrap helpers, preserving the
  established dense-route numerical contract while reducing repeated work.
- Categorical-only unconditional density routes now use the same
  profile-compression idea when `options(np.categorical.compress = TRUE)` is
  enabled. The fixed-bandwidth fit/evaluation route preserves dense-route
  fitted/evaluation values while avoiding repeated computation over identical
  categorical profiles, and the bandwidth-search route now uses the same
  compressed support representation for all-categorical data. As with other
  flat categorical search surfaces, selected smoothing parameters may drift by
  optimizer-path amounts while preserving the objective scale.
  Very fast
  compressed routes may remain overhead-floor limited, so MPI acceleration is
  most useful once the uncompressed work would be genuinely long-running.
- Categorical-only conditional density and conditional distribution bandwidth
  searches now honor `options(np.categorical.compress = TRUE)`. The promoted
  route preserves the objective value to numerical precision while allowing
  harmless optimizer-path drift in selected smoothing parameters, especially
  near upper-bound or large-bandwidth regions where the objective is flat.
- Ordered-only unconditional distribution bandwidth search and fit/evaluation
  routes also use profile compression when
  `options(np.categorical.compress = TRUE)` is enabled. The bandwidth-search
  route preserves the objective value to numerical precision while allowing
  harmless optimizer-path drift in selected smoothing parameters; fitted
  distribution values and standard errors are preserved while avoiding repeated
  computation over identical ordered profiles.
- Fixed-bandwidth local-constant `npscoef()` fits now use categorical-profile
  compression when all `Z` variables are categorical and
  `options(np.categorical.compress = TRUE)` is enabled. The route preserves
  fitted means, coefficient surfaces, asymptotic mean standard errors, and
  coefficient/gradient standard errors for training and evaluation fits while
  avoiding repeated work over duplicate `Z` profiles. The corresponding
  `npscoefhat(output = "apply")` path and count-based plot-bootstrap helper
  use the same profile compression without changing the explicit full-matrix
  `output = "matrix"` contract.
- Internal categorical-profile and large-bandwidth caches are now cleared at
  the relevant top-level density, distribution, conditional-density,
  conditional-distribution, and regression cleanup points. These caches are
  keyed by call-local row pointers, so clearing them per `.Call` prevents stale
  same-process state from leaking across unrelated data sets or MPI dispatch
  modes.
- Fixed `npcdens()` and `npcdist()` formula calls with explicit numeric
  smoothing parameters, such as `npcdist(y ~ x, data = dat, bws = c(.25,.25))`,
  so `npRmpi` preserves the established formula-to-bandwidth-object rewrite
  before MPI autodispatch.
- Hardened the `npudist()` formula route so formula calls are handled before
  MPI autodispatch.
- `npplreg()` now activates the already validated categorical regression
  compression path for its internal all-categorical `Z` regressions when
  `options(np.categorical.compress = TRUE)` is enabled, without requiring users
  to turn on the broader legacy `options(np.tree = TRUE)` switch.
- Formula variables whose names contain dots, such as `y.irr ~ x`, are no
  longer mistaken for the formula wildcard `.` in conditional density and
  conditional distribution bandwidth routes. The conditional-density bandwidth
  formula route also now expands the actual wildcard form `y ~ .` using the
  supplied `data` frame, matching the conditional-distribution route.
- Fixed MPI conditional-density and conditional-distribution NOMAD degree-search
  routes so Powell refinement and promoted wrappers such as `npconmode()` reach
  the intended bandwidth-object construction path rather than the pre-search
  autodispatch preflight used by non-degree-search routes.

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
