# np 0.70-6

* `npksum()` now validates its public logical controls and `kernel.pow` before
  native dispatch. This prevents malformed scalar controls from being masked
  by short-circuit evaluation and prevents empty, missing, non-finite, or
  vector-valued kernel powers from reaching C; an empty kernel power could
  previously terminate the R process. Documented scalar inputs and valid
  integer powers retain their established behavior.

* Repaired `npksum()` score and OCG output for mixed data. `compute.score=TRUE`
  now allocates its native result buffer, and packed score, OCG, and continuous
  permutation blocks (including derivative kernel weights) are returned in
  original data-column order. Categorical base-kernel factors are now retained
  in mixed-data continuous permutation blocks, and score-only ordered kernels
  no longer dereference OCG-only state. The defects could previously cause a
  process-level crash, fail during R reconstruction, or return an incomplete
  mixed-data derivative.

* Corrected Landweber--Fridman state coherence in `npregiv()`: state `N`
  now means exactly `N` completed updates, and `phi.mat[, N]`,
  `norm.stop[N]`, `norm.index`, the returned curve, derivatives, and optional
  weights refer to that same state. This can change the selected iteration and
  returned estimate. Bandwidth replay is now boundary-safe at `norm.index = 1`
  and retains recomputed stopping diagnostics. Documented numeric
  `starting.values` now initialize a complete level-and-derivative state.

* Corrected multivariate `npregiv()` derivatives to return one named pure
  coordinate partial per continuous structural coordinate. Local-constant
  Tikhonov derivative weights now use the derivative operator and bandwidth
  for the requested coordinate instead of combining coordinates or recycling
  bandwidth divisors. Univariate Tikhonov results are unchanged.

* Separated `npregivderiv()` training and evaluation state. A fit-time
  `zeval`/formula `newdata` grid now affects only `phi.prime.eval` and its path;
  training operators, residuals, centering, stopping, selected state, fitted
  values, and training gradients remain unchanged. Arbitrary positive grid
  sizes are supported, while a different `weval` is rejected because the
  inverse problem is defined on the training instruments. The default Issue
  57 local-linear training route is numerically unchanged.

* Hardened the IV public surface: ordinary IV regression stages now own
  `bandwidth.divide`, `ukertype`, and `okertype` without duplicate-argument
  failures and use kernel-appropriate categorical search bounds;
  `npregivderiv(random.seed=)` now controls its internal bandwidth searches;
  scalar controls fail early; summaries report continuous and categorical
  counts and derivative residual bandwidths; and `fitted()`, `residuals()`,
  `gradients()`, and plot methods consistently distinguish training from
  evaluation fields. Post-fit `predict()` remains deliberately unavailable.

* The omitted regression-smoothing choice in `npregivderiv()` is now local
  linear (`regtype = "ll"`, degree one), matching the longstanding `p = 1`
  default of `npregiv()`. Explicit `regtype = "lc"` reproduces the former
  omitted-default computation, while `regtype = "lp", degree = ...` remains
  available. Local-polynomial order applies only to continuous predictors;
  categorical-only internal stages use the equivalent local-constant route.
  This intentional default change can alter derivative paths, stopping states,
  and fitted structural functions for calls that previously omitted
  `regtype`.

* `npregivderiv()` now forwards user-supplied unordered and ordered
  categorical regression kernels to its internal regression stages without
  colliding with the private Equation (14) adjoint. The adjoint continues to
  own its required Li--Racine categorical kernels and ordinary-CDF
  normalization; previously successful calls and estimator defaults are
  unchanged.

* `npregivderiv()` progress now identifies the mathematical object being
  computed, including `E[y|w]`, initialization derivatives, conditional
  residuals, and the derivative adjoint `T*{E[y-phi(z)|w]}`. State-zero work
  is distinguished from work at completed iteration `N`; estimator results,
  stopping behavior, and smoothing choices are unchanged.

* Modernized the public `npregiv()` and `npregivderiv()` interfaces while
  preserving their established numerical defaults. Both now support explicit
  IV formulas (`y ~ z | w`, with optional `| x` where the estimator supports
  it), `data`, `subset`, `na.action`, fit-time `newdata`, fixed
  `regtype = "lc"`/`"ll"`/`"lp"` and scalar `degree` controls, structured
  summaries, and `fitted()`, `gradients()`, and training-row `residuals()`.
  Objects retain compact bandwidth, smoothing, and stage metadata. Legacy
  `p` remains supported by `npregiv()` with unchanged default behavior.
  Unsupported post-fit prediction, derivative `x`, and automatic-degree
  NOMAD routes are documented and fail explicitly instead of being silently
  accepted or approximated.

* Beta associated kernels now interpret \code{ckerbound="range"} using
  outer half-spacing bounds based on the two smallest and two largest distinct
  training values. This keeps every observation strictly inside the fitted
  support, removes empirical-extremum jumps and infinite raw-extremum
  derivatives, and preserves tied-extremum multiplicity. Explicit fixed beta
  bounds remain literal, while every non-beta range route continues to use
  exact sample minima and maxima.

* Beta range resolution now validates distinct-extrema metadata coordinate by
  coordinate. This fixes an unreleased 0.70-6 development defect in which two
  or more continuous beta coordinates failed during support resolution because
  vector and scalar logical operators were mixed.

* Automatic univariate second-order beta range searches now certify fixed-
  bandwidth density CVLS and distribution CDF solutions against the
  double-precision uniform-limit objective. A material improvement triggers
  one deterministic Powell refinement from the resolved support width, and
  the best material candidate is retained. The user's ordinary search and
  \code{nmulti} setting are unchanged, and cold starts and restarts from
  bandwidth objects share the same contract; summaries report the additional
  certification evaluations. Explicit fixed bounds, CVML, higher beta orders,
  nearest-neighbour modes, multivariate searches, and non-Powell solvers are
  unchanged.

* Corrected recursive fitted-value centering in `npregivderiv()` so both
  empirical terms of the Equation (14) adjoint use the same fitted
  conditional-residual vector. The recursive path had instead centered its
  second term on the raw residual even though its first term used the fitted
  residual. This can change derivative trajectories, stopping states, and
  fitted curves. The `npregiv()` and `npregivderiv()` examples now both use
  `n = 500`, for which the documented seed produces stable estimates.

* Corrected Landweber-Fridman state indexing in `npregivderiv()`. Iteration
  `N` now consistently denotes `N` completed derivative updates: column `N`
  of `phi.prime.mat` and `phi.mat`, `norm.stop[N]`, `num.iterations`, and the
  returned derivative/function now identify the same state. The initialization
  remains state zero, and evaluated stopping-rule overshoots remain available
  in the iteration matrices while `num.iterations` identifies the selected
  state. This can change the selected iteration and estimates because the
  stopping criterion is now evaluated as `N` times the residual norm at state
  `N`, as described by Florens, Centorrino, and Racine.

* Corrected the empirical adjoint in `npregivderiv()` to use the ordinary
  kernel CDF required by Equation (14) of Florens, Centorrino, and Racine.
  The integral kernel sum had retained an extra continuous-bandwidth factor,
  which could make the Landweber-Fridman stopping norm increase and bias the
  recovered conditional mean. The private adjoint now owns its required
  normalization, while public `npksum()` defaults and regression argument
  forwarding remain unchanged. The iteration guard also examines only the
  computed prefix of the preallocated stopping vector. This resolves issue
  #57.

* Retired the unused experimental truncated-Gaussian continuous kernel and
  its `nptgauss()` configuration helper. The supported continuous kernels are
  Gaussian, Epanechnikov, and uniform; their public interfaces, native codes,
  and numerical behavior are unchanged.

* Corrected heterogeneous generalized local-polynomial (GLP) construction to
  include every coordinate-capped term through the declared total degree. For
  example, `degree = c(2, 1)` now includes the `x1*x2` term. Raw and
  `bernstein.basis = TRUE` GLP fits now use exact representations of the same
  complete polynomial space; the latter uses a deterministic degree-graded,
  orthonormal shifted-Legendre representation. Raw heterogeneous GLP results
  can therefore differ from earlier versions when terms are restored;
  multivariate `bernstein.basis = TRUE` GLP and automatic degree-search results
  can also differ because the former fixed-degree Bernstein columns did not
  span the corrected GLP space once interactions entered. Univariate and
  maximum-degree-one GLP retain the same polynomial space, although their
  floating-point paths can differ. Dimension guards, higher-order
  derivatives, fitted values, standard errors, and the native evaluator use
  the same term definition. Local constant, local linear, univariate,
  additive, tensor, and non-local-polynomial semantics are unchanged.

* Hardened native NOMAD observer and interrupt handling. Explicit user
  interruption is now reported as an R `interrupt` condition only after
  native cleanup, ordinary observer errors remain fail-open, native callables
  are resolved per solve rather than retained across package reloads, and the
  declared R 3.5 compatibility floor is preserved.

* Restored timely progress updates during native NOMAD bandwidth and degree
  searches, including long compiled objective evaluations. Iteration, current
  degree, and accepted-best details now follow the package-wide progress
  interval without changing optimization results or evaluation accounting.
  This requires `crs` 0.15-46 or later.

* Restored partial mixed-degree local-polynomial gradient evaluation and
  plotting across regression, conditional density/distribution, least-squares
  quantile regression, conditional quantiles, and conditional modes. A
  continuous derivative requested above that coordinate's fitted degree is
  retained as `NA` and shown as an empty panel while other supported
  derivatives and categorical first-difference effects continue to plot.
  Requests with no available component now fail early with the affected
  predictor orders and degrees and actionable refitting guidance.

* Corrected common-scale uncertainty-band rendering for partially available
  gradient plots. Empty derivative panels now retain their list positions and
  finite-range calculations no longer emit spurious `min`/`max` warnings;
  supported estimates and interval bands are unchanged.

# np 0.70-5

* Corrected the heteroskedasticity-robust Ichimura index-coefficient
  covariance returned by `vcov(npindex(..., gradients = TRUE))` when the model
  has more than one free index coefficient. Link-gradient and residual weights
  are now applied observation by observation. Fits with one free coefficient
  remain exactly unchanged, as do beta, bandwidth, objective, fitted-value,
  residual, and gradient results.

* Repaired finite-support bounded continuous-kernel normalization for all
  supported Gaussian, Epanechnikov, uniform, and truncated-Gaussian kernels.
  The centered normalization now avoids tail cancellation, uses the platform
  C99 error function where appropriate, and preserves the uniform-density
  limit as the bandwidth tends to infinity.

* Replaced the bounded distribution operator with the mathematically correct
  observation-centered truncated-kernel CDF. Analytic centered interval
  primitives and hoisted invariant work preserve exact support endpoints,
  numerical accuracy, and the finite-support large-bandwidth uniform limit.

* Corrected `transform.bounds = TRUE` initialization for regression initial
  starts and multistarts and for later conditional-distribution multistarts.
  External bandwidth starts are now inverse-mapped before transformed Powell
  search; the public default remains `FALSE`.

* Clarified optimizer summaries so objective-cache lookups, NOMAD point
  lookups, family-native R-`optim` refinement, and total function-evaluation
  accounting are labelled consistently without implying that their
  denominators are interchangeable.

* Clarified bandwidth-object replay documentation across estimator families.
  Stored bandwidth objects retain search and fit metadata but do not silently
  materialize or expand training data; callers must provide data again where
  the documented replay route requires it.

* Standardized the public `nomad.opts` contract across supported bandwidth
  selectors and hardened native NOMAD callback cleanup, nearest-neighbor cache
  release, and native allocation-dimension guards.

* Reworked `npindexbw()` / `npindex()` internals after 0.70-4. Ichimura and
  Klein-Spady single-index objectives now reuse the established `npreg`
  leave-one-out backend where applicable, preserving the public estimator
  contract while materially improving high-dimensional and local-polynomial
  objective evaluation. The single-index formula/materialization routes were
  tightened so formula-selected variables, explicit bandwidth objects, and
  direct formula calls with explicit smoothing parameters reenter consistently.

* Added user documentation for `npindex()` optimizer choice and beta
  interpretation. The documentation now gives practical guidance on when
  derivative-free Nelder-Mead remains a reasonable low-dimensional default and
  when BFGS is useful for higher-dimensional index searches, and it explains
  the relative interpretation of normalized single-index beta coefficients.

* Hardened single-index fit, evaluation, plotting, and summary behavior.
  Bounded continuous-kernel options are now carried consistently through
  objective, fit, evaluation, variance, and bootstrap routes; Klein-Spady
  confusion-matrix output is guarded against out-of-range fitted values; and
  large-bandwidth shortcuts are restricted to kernel/order combinations whose
  constant-weight approximation was validated.

* Repaired single-index plot bootstrap memory usage. Wild-bootstrap
  self-maps now avoid public-facing `O(n^2)` kernel-weight allocations for the
  ordinary large-sample plot route, while preserving the public fitted and
  plotting contracts.

* Repaired generalized-nearest-neighbor local-polynomial derivative ownership
  across `npreghat()` and related public routes. Mixed-degree local-polynomial
  fits now route available derivative components through the correct owner,
  preserve `apply == H %*% y` contracts, and report unavailable derivative
  components consistently rather than silently applying the wrong operator.

* Repaired local-constant derivative ownership in `npreghat()` so degree-zero
  local-polynomial derivative requests use the analytic local-constant
  derivative contract only when that contract is mathematically available.
  Scalar, matrix, and multi-column apply routes were validated separately.

* Improved conditional mixed-degree local-polynomial gradients for
  `npcdens()` / `npcdist()` and public conditional-gradient accessors. The
  partial-availability contract now matches the regression-family policy:
  available components are returned, unavailable components are represented as
  `NA`, and incoherent metadata fails clearly.

* Repaired formula/data reentry contracts across density, distribution,
  conditional density, conditional distribution, single-index, smooth
  coefficient, partially linear, quantile, conditional-mode, copula, and
  significance-test routes. Explicit estimator `data=` now overrides stored
  bandwidth-object data where that public call shape is supported, formula
  `newdata` is validated against fitted RHS variables, and direct formula
  calls with numeric smoothing parameters no longer misroute formula objects
  as native data.

* Repaired unconditional density and distribution edge contracts. Ordered
  kernel-code selection is now consistent between conditional bandwidth
  selection and fitting where normalization is required, all-NA input is
  rejected before native calls, training and evaluation omission metadata are
  retained separately where needed, categorical zero-bandwidth standard-error
  handling is consistent across categorical configurations, and unsupported
  bandwidth-selection method codes fail clearly.

* Clarified normal-reference bandwidth documentation for density and
  distribution routes. The rule-of-thumb formulas are documented as fast
  exploratory Silverman-style heuristics, not production substitutes for
  cross-validation or likelihood-based selection.

* Hardened `npcdistbw()` normal-reference method handling so the R method code
  and native C method code stay aligned, avoiding accidental fall-through to
  an unintended bandwidth-selection branch.

* Repaired and hardened public estimator contracts found during adversarial
  audits of `npscoef`, `npplreg`, `npqreg`, `nplsqreg`, `npconmode`,
  `npcopula`, `npreg`, `npcdens`, `npcdist`, `npudens`, and `npudist`.
  Repairs include `npscoef` iterated backfit behavior, partially-linear fit
  reentry, quantile inversion/clamping contracts, conditional-mode
  probability/tie handling, copula sample reentry, regression tree/large-h
  predicate alignment, and density/distribution formula/native argument
  consistency.

* Plot-bootstrap memory hardening now covers single-index, conditional-mode,
  and partially-linear plot routes. The default plot evaluation grids remain
  linear in the training sample size, and explicit `neval == ntrain` style
  requests remain user-controlled.

* Categorical-gradient and derivative workspace handling were repaired so
  no-gradient and gradient calls can be mixed safely in one session without
  stale native derivative workspace state affecting later calls.

* Native shadow-object and cache lifecycles were hardened. Conditional-density
  shadow pointers, native objective-cache state, and regression large-h /
  large-lambda caches are cleared at the appropriate top-level lifecycle
  boundary so pointer-keyed helper state cannot leak across independent calls
  or datasets.

* `npksum()` numeric-bandwidth dispatch now constructs data-aware bandwidth
  objects for mixed continuous/categorical data instead of falling back to a
  default-typed object. The serial and MPI packages now share the same public
  behavior for numeric smoothing parameters.

* The local-polynomial regression CVKS low-support objective path now uses a
  named bandwidth-method code and aligns the R and C routing contracts,
  reducing the risk of drift between search metadata and native objective
  handling.

* Expanded focused tests, demos, benchmarks, and release-protocol sentinels
  across the public exported surface. The release protocol now requires an
  explicit public-exported-surface inventory, estimator-family sentinels,
  documentation/demo/benchmark smoke coverage, and installed package proof
  before release-ready claims.

# np 0.70-4

* Hardened proactive C cleanup paths by clearing extended nearest-neighbor
  alias state in the shared estimator cleanup helper, routing selected
  bandwidth-constructor allocation and unsupported-method failures through
  existing cleanup labels, and guarding large `np_kernelsum()` allocation-size
  products before calling the legacy `alloc_vecd(int)` allocator.

* Added `nomad = "auto"` for local-polynomial degree searches. The automatic
  policy uses exhaustive Powell-style degree search for small
  one-dimensional degree lattices where evidence showed it is more reliable
  than heuristic NOMAD restarts, while preserving NOMAD for larger or
  explicitly requested search surfaces.

* Bandwidth and fit summaries now report cumulative search diagnostics more
  clearly. NOMAD cache output distinguishes repeated point lookups avoided
  by NOMAD from true objective computations; Powell summaries expose
  repeated objective lookups avoided by the package-side cache; hybrid
  NOMAD+Powell and exhaustive Powell timing labels are reported consistently.

* Base-graphics plot scaling and legends have been hardened for multi-panel
  displays. Plots now honor active `mfrow`/`cex` behavior more consistently,
  use role-appropriate legend sizes, and draw factor legends with point glyphs
  matching the plotted estimates.

* Fixed-bandwidth bias-corrected bootstrap plot intervals have been reworked
  across the supported plot families. Bias-corrected centers and intervals
  now share a common centering contract for regression, unconditional density
  and distribution, conditional density and distribution, single-index,
  partially linear, and smooth-coefficient plots where supported.

* Pair/block/geometric bootstrap intervals in regression-style routes now use
  smooth-bootstrap bias correction when `center = "bias-corrected"` is
  requested; wild-bootstrap regression intervals retain the standard
  wild-bootstrap correction. Density and distribution routes use
  perturbation-based smooth-bootstrap bias correction. Gaussian, uniform, and
  second-order Epanechnikov perturbation kernels are supported; higher-order
  signed perturbation kernels fail closed.

* Bias-corrected plot support now covers mixed-data unconditional density and
  distribution, mixed-data conditional density and distribution, and
  conditional gradient displays. Non-fixed/adaptive/generalized bandwidth
  bias correction and empirically unsupported `npqreg` bias correction remain
  fail-closed with explicit messages.

* Bias-corrected plot payloads returned by `output = "data"` and
  `output = "both"` now consistently expose fitted values, bias-corrected
  values, gradients, gradient bias corrections, and interval payloads with the
  same centering contract used for rendering.

* Derivative-order validation has been tightened across plot/gradient routes
  so unsupported derivative requests fail early or are represented
  consistently rather than silently plotting a lower-order derivative.

* Single-index formula dispatch with explicit bandwidths has been repaired,
  and single-index bias-corrected plot centers, output payloads, and legends
  have been tightened.

* The native `loadNamespace("crs")` call construction used by CRS-backed
  native search routes is now protected, resolving the RCHK protection finding
  without changing the native NOMAD search contract.

# np 0.70-3

* Added `nplsqreg()`/`nplsqregbw()` as a location-scale quantile-regression
  front end with formula/data and bandwidth-object workflows, scalar/vector
  `tau`, prediction, residual extraction, summaries, and plot routes built on
  the shared quantile plotting engine.

* Supported MADS/NOMAD-backed bandwidth-search routes now use the final native
  `crs` NOMAD C API rather than the retired legacy `snomadr()` fallback.
  This covers the promoted regression, density, distribution, conditional
  density, conditional distribution, smooth-coefficient, single-index,
  partially linear, and location-scale quantile search surfaces where those
  routes support native NOMAD/MADS. The runtime dependency on `crs`
  is now declared in `Imports`, while `LinkingTo` remains for the native
  header.

* Native NOMAD routes now preserve progress best-record reporting, expose
  cache/evaluation diagnostics, honor explicit start and option controls, and
  reject unsupported or indeterminate cache-off settings before solver entry.
  Inadmissible GLP degree candidates are guarded before expensive evaluator
  work.

* `npindexbw(..., method = "ichimura", regtype = c("ll", "lp"))` now reuses
  the established local-polynomial regression objective evaluator for
  fixed-degree and NOMAD degree-search routes. Focused sentinel runs preserved
  selected objective payloads while materially reducing runtime for
  local-linear and local-polynomial Ichimura single-index bandwidth searches.

* `options(np.tree = "auto")` is now the default tree mode. In auto mode,
  continuous kd-tree routes are enabled only for bounded-support continuous
  kernels (`"epanechnikov"` and `"uniform"`); `np.tree = TRUE` remains the
  explicit force-on override and `np.tree = FALSE` remains the force-off
  diagnostic path.

* Powell bandwidth searches now expose package-side repeated-candidate
  objective caching through `options(np.objective.cache = TRUE/FALSE)`. The
  cache remains enabled by default and is scoped to one bandwidth solve, so it
  can reuse exact candidates across Powell restarts without carrying state
  across datasets or later calls. Continuous-only generalized/adaptive
  nearest-neighbor routes also retain their integer nearest-neighbor objective
  cache under the same switch; NOMAD solver caching and extended-NN distance
  reuse remain separate mechanisms.

* Continuous large-bandwidth shortcut evaluations can now be disabled with
  `options(np.largeh = FALSE)`, and discrete near-upper-bandwidth shortcut
  evaluations can now be disabled with `options(np.largelambda = FALSE)`.
  Both remain enabled by default. These switches are intended for diagnostic
  timing and reproducibility studies that need to separate tree effects from
  large-bandwidth and large-lambda fast paths without changing the canonical
  dense/tree objective machinery.

* Local-polynomial regression cross-validation now uses a leaner hot
  symmetric weighted-sum loop. Fixed-bandwidth `npregbw(..., regtype = "lp",
  bwmethod = "cv.ls")` objective probes show substantially faster
  local-polynomial CV evaluation while preserving objective values to
  numerical precision; adjacent density bandwidth probes preserve their
  objective values as well.

* Shared weighted outer-product accumulation in `npksum()` now uses a guarded
  BLAS `dgemm` route when the operation is dense, non-permuted, and
  memory-bounded. Focused fixed-bandwidth probes preserve objective values to
  numerical precision while substantially accelerating high-basis
  local-polynomial regression and smooth-coefficient objective rows; small and
  scalar routes remain on the established loop path.

* Unconditional density least-squares cross-validation now uses a leaner
  fixed-bandwidth Gaussian convolution loop. Fixed-bandwidth
  `npudensbw(..., bwmethod = "cv.ls")` objective probes preserve objective
  values exactly in the focused validation rows while materially reducing the
  convolution portion of the objective calculation. Conditional-density
  least-squares objective probes inherit the same fixed-bandwidth Gaussian
  convolution improvement.

* Non-Gaussian scalar-bandwidth convolution helpers now hoist the response
  bandwidth power outside the inner loop, improving fixed-bandwidth
  least-squares density cross-validation with compact-support kernels while
  preserving objective values exactly in focused probes.

* Continuous-kernel vector helpers now reuse the loop-invariant signed inverse
  bandwidth scale inside their inner loops. Focused density, conditional
  density, and regression objective probes preserved objective values exactly
  while reducing repeated scaling work in shared C hot paths.

* Conditional density and conditional distribution least-squares
  cross-validation now use a size-aware row-block policy for local-polynomial
  objective evaluation. The accepted route keeps the bounded-quadrature cap
  unchanged, bounds transient memory by sample size, and preserves objective
  values to numerical precision while materially reducing evaluator overhead
  for fixed-bandwidth CVLS probes.

* Local-polynomial conditional density maximum-likelihood cross-validation now
  uses the same bounded-memory block machinery for fixed and generalized
  nearest-neighbor bandwidths. Focused `npcdensbw(..., bwmethod = "cv.ml",
  regtype = "lp")` probes preserve objective values and selected bandwidths to
  numerical precision while reducing objective and full-search runtime.

* Large-sample categorical-only regression now has a profile-compressed
  execution route controlled by `options(np.categorical.compress = TRUE)`,
  which is enabled by default. This categorical route is independent of
  `options(np.tree)`. For local constant categorical regression, repeated
  predictor profiles are compressed before fitting, prediction/evaluation,
  standard errors, gradients where meaningful, bandwidth search, hat-helper
  use, and plot bootstrap helpers.
  This preserves the established dense-route numerical contract while greatly
  reducing work for large samples with many repeated
  factor/ordered predictor combinations.

* Categorical-only unconditional density routes now use the same
  profile-compression idea when `options(np.categorical.compress = TRUE)` is
  enabled. The fixed-bandwidth fit/evaluation route preserves dense-route
  fitted/evaluation values while avoiding repeated computation over identical
  categorical profiles, and the bandwidth-search route now uses the same
  compressed support representation for all-categorical data. As with other
  flat categorical search surfaces, selected smoothing parameters may drift by
  optimizer-path amounts while preserving the objective scale.

* Categorical-only conditional density and conditional distribution bandwidth
  searches now honor `options(np.categorical.compress = TRUE)`. The promoted
  route preserves the objective value to numerical precision while allowing
  harmless optimizer-path drift in selected smoothing parameters, especially
  near upper-bound or large-bandwidth regions where the objective is flat.

* Ordered-only unconditional distribution bandwidth search and fit/evaluation
  routes also use profile compression when
  `options(np.categorical.compress = TRUE)` is enabled. The bandwidth-search
  route preserves the objective value to numerical precision while allowing
  harmless optimizer-path drift in selected smoothing parameters; fitted
  distribution values and standard errors are preserved while avoiding repeated
  computation over identical ordered profiles.

* Fixed-bandwidth local-constant `npscoef()` fits now use categorical-profile
  compression when all `Z` variables are categorical and
  `options(np.categorical.compress = TRUE)` is enabled. The route preserves
  fitted means, coefficient surfaces, asymptotic mean standard errors, and
  coefficient/gradient standard errors for training and evaluation fits while
  avoiding repeated work over duplicate `Z` profiles. The corresponding
  `npscoefhat(output = "apply")` path and count-based plot-bootstrap helper
  use the same profile compression without changing the explicit full-matrix
  `output = "matrix"` contract.

* Internal categorical-profile and large-bandwidth caches are now cleared at
  the relevant top-level density, distribution, conditional-density,
  conditional-distribution, and regression cleanup points. These caches are
  keyed by call-local row pointers, so clearing them per `.Call` prevents stale
  same-process state from leaking across unrelated data sets.

* Formula variables whose names contain dots, such as `y.irr ~ x`, are no
  longer mistaken for the formula wildcard `.` in conditional density and
  conditional distribution bandwidth routes. The conditional-density bandwidth
  formula route also now expands the actual wildcard form `y ~ .` using the
  supplied `data` frame, matching the conditional-distribution route.

# np 0.70-2

* `npqreg()` is now a fully fledged quantile-regression front end. It
  supports the formula/data workflow, internally computes
  `npcdistbw()` bandwidths when a bandwidth object is not supplied,
  accepts scalar or vector `tau`, reuses selected bandwidths for
  additional quantiles in `plot()`, and exposes the usual S3 surface:
  `fitted()`, `predict()`, `predict(..., se.fit=TRUE)`, `se()`,
  `gradients()`, `summary()`, `print()`, `quantile()`, and `plot()`.

* `npqreg()` prediction now honors the standard `newdata` workflow while
  preserving native `exdat` precedence for compatibility with existing
  `np` call surfaces. Formula-based prediction validates that new data
  contain the required right-hand-side variables.

* `npqreg()` plotting has been expanded for vector quantiles,
  level/gradient displays, ordered predictors, user-specified legends,
  and object-fed plotting of additional `tau` values without recomputing
  cross-validation.

* `npconmode()` is now a first-class conditional-mode estimator. It
  supports formula/data and bandwidth-object workflows, forwards
  bandwidth-selection options to `npcdensbw()`, propagates local
  polynomial and NOMAD metadata, and exposes `fitted()`, `predict()`,
  `summary()`, `print()`, `gradients()`, and `plot()` methods.

* `npconmode()` now supports optional class-probability matrices and
  level-specific probability gradients. For non-local-constant fits,
  probabilities are normalized to be non-negative and to sum to one
  across the discrete response support before modal classification.

* `npconmode()` now fails early for non-categorical responses and
  validates formula-based `newdata` against the original right-hand-side
  variables.

* `npconmode()` plotting now supports object-fed class-probability slices
  and two-dimensional probability surfaces, optional `rgl` rendering, and
  probability-level asymptotic intervals where defined. Surface bootstrap
  intervals for class probabilities remain intentionally deferred.

* `npcopula()` is now a first-class copula estimator. It supports
  formula/data and bandwidth-object workflows, automatic two-dimensional
  probability grids, explicit `u` evaluation grids, and ordinary
  extractable object components including `$bws`.

* `npcopula()` now provides `fitted()`, `predict()`, `predict(...,
  se.fit=TRUE)`, `se()`, `summary()`, `print()`, `as.data.frame()`, and
  richer `plot()` methods. Plotting supports base `persp`, `image`, and
  optional `rgl` rendering, with asymptotic and bootstrap intervals for
  copula surfaces where defined.

* `npcopula()` explicit-grid evaluation now uses the direct estimator
  route, preserving numerical results while avoiding the severe runtime
  growth of the previous expanded-grid path when users request larger
  probability grids.

* The automatic local-polynomial NOMAD controls have been split into
  explicit restart toggles: `powell.remin` for Powell restarts and
  `nomad.remin` for the second NOMAD hot start. This preserves the
  Powell Numerical Recipes restart default while allowing NOMAD hot
  starts to be controlled separately.

* Deprecated legacy `remin` remains accepted by `npregbw()` and `npreg()`
  with a warning and is mapped to the modern `powell.remin`/`nomad.remin`
  controls where appropriate, preserving downstream compatibility while
  documenting the new spelling.

* Hat-operator helpers now support an additional constraint-oriented
  output route for objects needed by shape-constrained quadratic
  programming workflows, avoiding reimplementation of local-polynomial
  hat-matrix construction in user examples.

* Local-polynomial derivative support has been broadened across the
  conditional estimator family. `npreg()`, `npcdens()`, and `npcdist()`
  now honor `gradient.order` more consistently for fitted, evaluated,
  predicted, and plotted objects when the selected polynomial degree is
  high enough, including vector derivative orders over continuous
  predictors and tensor/additive/Bernstein local-polynomial bases.

* Core and semiparametric S3 prediction paths have been hardened around
  `newdata`, native evaluation-argument precedence, formula RHS
  validation, and `se.fit` handling.

* Front-end/bandwidth argument hygiene has been tightened so
  estimator-only controls such as `proper` are not forwarded into
  bandwidth selectors that do not accept them.

* Documentation has been refreshed for the promoted `npqreg()`,
  `npconmode()`, and `npcopula()` workflows, including the
  local-polynomial NOMAD route, probability/gradient outputs, plot
  controls, and examples that use the streamlined interfaces.

* The pre-release validation suite was expanded with focused hostile
  argument tests, S3 contract tests, installed/tarball proof scripts,
  and cross-package parity checks for the newly promoted estimator
  families.

# np 0.70-1

* The default multistart cap for bandwidth selection now follows
  `min(2, p)` across the core estimator families, replacing the older
  `min(5, p)` cap. This includes automatic LP degree-search calls when
  `search.engine="nomad"` or `"nomad+powell"` and `nmulti` is not
  supplied explicitly.

* The univariate boundary density helper `npuniden.boundary()` now
  defaults to `nmulti=1`.

* The empirical studies supporting this change are documented under
  `benchmarks/validation/`.

* LP-capable front ends now accept `nomad=TRUE` as a documented
  convenience preset for the recommended automatic NOMAD local-polynomial
  route. Missing settings expand to the same long-form LP/NOMAD defaults
  documented in the bandwidth help pages, and regression formula calls
  such as `npreg(y ~ x, nomad = TRUE)` now carry that shortcut through
  the internally computed bandwidth path.
