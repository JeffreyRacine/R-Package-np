# np 0.70-3

- Added `nplsqreg()`/`nplsqregbw()` as a location-scale quantile-regression
  front end with formula/data and bandwidth-object workflows, scalar/vector
  `tau`, prediction, residual extraction, summaries, and plot routes built on
  the shared quantile plotting engine.
- Supported MADS/NOMAD-backed bandwidth-search routes now use the final native
  `crs` NOMAD C API rather than the retired legacy `snomadr()` fallback.
  This covers the promoted regression, density, distribution, conditional
  density, conditional distribution, smooth-coefficient, single-index,
  partially linear, and location-scale quantile search surfaces where those
  routes support native NOMAD/MADS. The runtime dependency on `crs`
  is now declared in `Imports`, while `LinkingTo` remains for the native
  header.
- Native NOMAD routes now preserve progress best-record reporting, expose
  cache/evaluation diagnostics, honor explicit start and option controls, and
  reject unsupported or indeterminate cache-off settings before solver entry.
  Inadmissible GLP degree candidates are guarded before expensive evaluator
  work.
- `npindexbw(..., method = "ichimura", regtype = c("ll", "lp"))` now reuses
  the established local-polynomial regression objective evaluator for
  fixed-degree and NOMAD degree-search routes. Focused sentinel runs preserved
  selected objective payloads while materially reducing runtime for
  local-linear and local-polynomial Ichimura single-index bandwidth searches.
- `options(np.tree = "auto")` is now the default tree mode. In auto mode,
  continuous kd-tree routes are enabled only for bounded-support continuous
  kernels (`"epanechnikov"` and `"uniform"`); `np.tree = TRUE` remains the
  explicit force-on override and `np.tree = FALSE` remains the force-off
  diagnostic path.
- Powell bandwidth searches now expose package-side repeated-candidate
  objective caching through `options(np.objective.cache = TRUE/FALSE)`. The
  cache remains enabled by default and is scoped to one bandwidth solve, so it
  can reuse exact candidates across Powell restarts without carrying state
  across datasets or later calls. Continuous-only generalized/adaptive
  nearest-neighbor routes also retain their integer nearest-neighbor objective
  cache under the same switch; NOMAD solver caching and extended-NN distance
  reuse remain separate mechanisms.
- Continuous large-bandwidth shortcut evaluations can now be disabled with
  `options(np.largeh = FALSE)`, and discrete near-upper-bandwidth shortcut
  evaluations can now be disabled with `options(np.largelambda = FALSE)`.
  Both remain enabled by default. These switches are intended for diagnostic
  timing and reproducibility studies that need to separate tree effects from
  large-bandwidth and large-lambda fast paths without changing the canonical
  dense/tree objective machinery.
- Local-polynomial regression cross-validation now uses a leaner hot
  symmetric weighted-sum loop. Fixed-bandwidth `npregbw(..., regtype = "lp",
  bwmethod = "cv.ls")` objective probes show substantially faster
  local-polynomial CV evaluation while preserving objective values to
  numerical precision; adjacent density bandwidth probes preserve their
  objective values as well.
- Shared weighted outer-product accumulation in `npksum()` now uses a guarded
  BLAS `dgemm` route when the operation is dense, non-permuted, and
  memory-bounded. Focused fixed-bandwidth probes preserve objective values to
  numerical precision while substantially accelerating high-basis
  local-polynomial regression and smooth-coefficient objective rows; small and
  scalar routes remain on the established loop path.
- Unconditional density least-squares cross-validation now uses a leaner
  fixed-bandwidth Gaussian convolution loop. Fixed-bandwidth
  `npudensbw(..., bwmethod = "cv.ls")` objective probes preserve objective
  values exactly in the focused validation rows while materially reducing the
  convolution portion of the objective calculation. Conditional-density
  least-squares objective probes inherit the same fixed-bandwidth Gaussian
  convolution improvement.
- Non-Gaussian scalar-bandwidth convolution helpers now hoist the response
  bandwidth power outside the inner loop, improving fixed-bandwidth
  least-squares density cross-validation with compact-support kernels while
  preserving objective values exactly in focused probes.
- Continuous-kernel vector helpers now reuse the loop-invariant signed inverse
  bandwidth scale inside their inner loops. Focused density, conditional
  density, and regression objective probes preserved objective values exactly
  while reducing repeated scaling work in shared C hot paths.
- Conditional density and conditional distribution least-squares
  cross-validation now use a size-aware row-block policy for local-polynomial
  objective evaluation. The accepted route keeps the bounded-quadrature cap
  unchanged, bounds transient memory by sample size, and preserves objective
  values to numerical precision while materially reducing evaluator overhead
  for fixed-bandwidth CVLS probes.
- Local-polynomial conditional density maximum-likelihood cross-validation now
  uses the same bounded-memory block machinery for fixed and generalized
  nearest-neighbor bandwidths. Focused `npcdensbw(..., bwmethod = "cv.ml",
  regtype = "lp")` probes preserve objective values and selected bandwidths to
  numerical precision while reducing objective and full-search runtime.
- Large-sample categorical-only regression now has a profile-compressed
  execution route controlled by `options(np.categorical.compress = TRUE)`,
  which is enabled by default. This categorical route is independent of
  `options(np.tree)`. For local constant categorical regression, repeated
  predictor profiles are compressed before fitting, prediction/evaluation,
  standard errors, gradients where meaningful, bandwidth search, hat-helper
  use, and plot bootstrap helpers.
  This preserves the established dense-route numerical contract while greatly
  reducing work for large samples with many repeated
  factor/ordered predictor combinations.
- Categorical-only unconditional density routes now use the same
  profile-compression idea when `options(np.categorical.compress = TRUE)` is
  enabled. The fixed-bandwidth fit/evaluation route preserves dense-route
  fitted/evaluation values while avoiding repeated computation over identical
  categorical profiles, and the bandwidth-search route now uses the same
  compressed support representation for all-categorical data. As with other
  flat categorical search surfaces, selected smoothing parameters may drift by
  optimizer-path amounts while preserving the objective scale.
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
  same-process state from leaking across unrelated data sets.
- Formula variables whose names contain dots, such as `y.irr ~ x`, are no
  longer mistaken for the formula wildcard `.` in conditional density and
  conditional distribution bandwidth routes. The conditional-density bandwidth
  formula route also now expands the actual wildcard form `y ~ .` using the
  supplied `data` frame, matching the conditional-distribution route.

# np 0.70-2

- `npqreg()` is now a fully fledged quantile-regression front end. It
  supports the formula/data workflow, internally computes
  `npcdistbw()` bandwidths when a bandwidth object is not supplied,
  accepts scalar or vector `tau`, reuses selected bandwidths for
  additional quantiles in `plot()`, and exposes the usual S3 surface:
  `fitted()`, `predict()`, `predict(..., se.fit=TRUE)`, `se()`,
  `gradients()`, `summary()`, `print()`, `quantile()`, and `plot()`.
- `npqreg()` prediction now honors the standard `newdata` workflow while
  preserving native `exdat` precedence for compatibility with existing
  `np` call surfaces. Formula-based prediction validates that new data
  contain the required right-hand-side variables.
- `npqreg()` plotting has been expanded for vector quantiles,
  level/gradient displays, ordered predictors, user-specified legends,
  and object-fed plotting of additional `tau` values without recomputing
  cross-validation.
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
  optional `rgl` rendering, with asymptotic and bootstrap intervals for
  copula surfaces where defined.
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
  predictors and tensor/additive/Bernstein local-polynomial bases.
- Core and semiparametric S3 prediction paths have been hardened around
  `newdata`, native evaluation-argument precedence, formula RHS
  validation, and `se.fit` handling.
- Front-end/bandwidth argument hygiene has been tightened so
  estimator-only controls such as `proper` are not forwarded into
  bandwidth selectors that do not accept them.
- Documentation has been refreshed for the promoted `npqreg()`,
  `npconmode()`, and `npcopula()` workflows, including the
  local-polynomial NOMAD route, probability/gradient outputs, plot
  controls, and examples that use the streamlined interfaces.
- The pre-release validation suite was expanded with focused hostile
  argument tests, S3 contract tests, installed/tarball proof scripts,
  and cross-package parity checks for the newly promoted estimator
  families.

# np 0.70-1

- The default multistart cap for bandwidth selection now follows
  `min(2, p)` across the core estimator families, replacing the older
  `min(5, p)` cap. This includes automatic LP degree-search calls when
  `search.engine="nomad"` or `"nomad+powell"` and `nmulti` is not
  supplied explicitly.
- The univariate boundary density helper `npuniden.boundary()` now
  defaults to `nmulti=1`.
- The empirical studies supporting this change are documented under
  `benchmarks/validation/`.
- LP-capable front ends now accept `nomad=TRUE` as a documented
  convenience preset for the recommended automatic NOMAD local-polynomial
  route. Missing settings expand to the same long-form LP/NOMAD defaults
  documented in the bandwidth help pages, and regression formula calls
  such as `npreg(y ~ x, nomad = TRUE)` now carry that shortcut through
  the internally computed bandwidth path.
