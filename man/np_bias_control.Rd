\name{np_bias_control}
\alias{np_bias_control}
\title{Auxiliary Controls for Bias-Corrected Regression Plot Intervals}
\description{
  Constructs the auxiliary settings for \code{plot.npregression} and
  regression bandwidth plots with \code{center="bias-corrected"}.
  The fitted model's bandwidths, polynomial degree and basis remain fixed.
}
\usage{
np_bias_control(construction = c("iterated", "legacy"),
                iterations = 16L, degree = NULL, bernstein.basis = TRUE,
                p0.scale = 4, p1.scale = 1.5,
                operator.scale = 4/3, variance.scale = 4/3,
                p0.bandwidth = NULL, p1.bandwidth = NULL,
                operator.bandwidth = NULL, variance.bandwidth = NULL,
                adaptation = c("auto", "none", "local-to-global"))
}
\arguments{
  \item{adaptation}{Coordinatewise local-to-global auxiliary correction.
    \code{"auto"} uses the validated transition for continuous-only fixed
    Gaussian GLP models when P0, P1 and operator bandwidth controls are omitted.
    Local coordinates remain unchanged. Explicit geometry bandwidth or scale
    overrides retain the original construction under \code{"auto"}; iterations,
    auxiliary degree and variance controls remain effective.
    \code{"none"} retains the frozen pre-adaptation construction.
    \code{"local-to-global"} explicitly applies the transition, treating
    supplied P0, P1 and operator bandwidths as reference values that may grow.
    It errors for unsupported model classes. Older serialized version-1/2
    controls retain \code{"none"}. Main model bandwidths and degrees never change.}

  \item{construction}{\code{"iterated"} selects the corrected observation
    operator and its full covariance. \code{"legacy"} explicitly requests
    the previous bootstrap centering construction; it requires bootstrap
    errors and cannot be combined with auxiliary overrides.}
  \item{iterations}{Positive integer number of residual correction iterations.
    The default is 16. Larger values increase computation and need not improve
    coverage monotonically.}
  \item{degree}{Auxiliary polynomial degree, distinct from the fitted model's
    degree. The default \code{NULL} resolves coordinate by coordinate to
    \code{pmax(5L, main.degree)}. A positive integer scalar applies to every
    continuous predictor; an unnamed vector follows continuous-predictor order,
    and a named vector must name every continuous predictor exactly. Explicit
    values are literal overrides. The auxiliary GLP space must contain the
    main polynomial space. Profile-specific auxiliaries require positive
    residual degrees of freedom in every profile; pooled auxiliaries require
    this for each profile main space and for their joint auxiliary space.
    Actual local rank checks also apply.}
  \item{bernstein.basis}{Logical. Must be \code{TRUE} for auxiliary degree
    greater than one. The canonical stable GLP representation (degree-graded shifted Legendre)
    is used for the auxiliary polynomial space, including derivative rows.
    This control never changes the fitted main model basis.}
  \item{p0.scale, p1.scale, operator.scale, variance.scale}{Positive auxiliary
    bandwidth scale factors. Each may be a scalar, an unnamed vector in
    continuous-predictor order, or a vector named for every continuous
    predictor. The default rule is given below.}
  \item{p0.bandwidth, p1.bandwidth, operator.bandwidth, variance.bandwidth}{Optional
    positive bandwidths in predictor units, with the same vector conventions.
    For each role, an explicit bandwidth and an explicit scale are mutually
    exclusive. Raw bandwidths must be transformed when predictor units change.}
}
\details{
  Supply \code{bias.control=np_bias_control(...)} to the regression plot.
  Omitting this argument with \code{center="bias-corrected"} uses the iterated
  defaults. \code{center="estimate"} remains the ordinary plot default.

  For \eqn{d} continuous predictors, let \eqn{s_j} be the package's
  \code{EssDee} spread on the retained training sample and let \eqn{n_{min}}
  be the smallest observed joint categorical-profile size (the full sample
  size without categorical predictors). Each role uses
  \deqn{h_j = c_j s_j n_{min}^{-1/(4+d)}.}
  Thus \eqn{c_j} is dimensionless. These are profile-adjusted auxiliary
  scales, not the main model's bandwidth scale factors. Effective full-sample
  package scale factors are also recorded. Centering and scaling the internal
  coordinates preserves the fitted model; derivative units are restored.

  The P0 local polynomial fits the response on the complete continuous
  polynomial basis separately within each categorical profile. Its projection
  geometry supplies polynomial correction and its residuals supply variance
  estimation. The stable Gaussian operator applies the specified number of
  residual correction iterations. P1 supplies an additional local polynomial
  projection; its complete mean or derivative row is included in the final
  corrected observation row. Variance estimation and calibration use that
  complete row, including correction covariance. No training-square smoothing
  matrix is retained. Training preparation is reused and query rows are blocked.

  For an active local-to-global transition, each continuous coordinate uses
  the dimensionless ratio \eqn{r_j=\max(1,h_j/b_{0j})}, where \eqn{h_j} is
  the fitted bandwidth and \eqn{b_{0j}} the reference P0 bandwidth in the same
  units. P0, P1 and residual-operator bandwidths are multiplied by \eqn{r_j};
  the variance-smoothing bandwidth is unchanged. A fixed coordinatewise
  polynomial penalty reduces auxiliary components beyond the main model's
  polynomial space as that coordinate becomes global. The frozen residual
  iteration normalization makes the same transition. This is an auxiliary
  inference construction, not degree selection or a change to the main fit.
  Adapted P0 row-variance-normalized residuals supply variance estimation and
  residual-score calibration. Actual fitted observation rows are reused,
  including each complete plot panel's large-bandwidth weight decisions.
  Returned settings distinguish requested reference and effective bandwidths,
  effective scale factors, activation reason and transition geometry.
  Categorical and NN models retain their existing correction under
  \code{"auto"}. Continuous dimension and main degree are not capped by this
  activation rule; coverage evidence for the transition concentrates on
  one/two continuous predictors and LC/LL/quadratic models, with numerical
  validation of higher and mixed degrees. This does not establish general
  high-dimensional or high-degree coverage. Transition coverage was assessed
  with Gaussian and Mammen calibration; other laws retain their stated
  finite-sample qualifications.

  The implementation admits an unbounded second-order Gaussian kernel,
  LC, LL and package-admissible GLP main degree vectors, including mixed
  degrees, with any positive number of continuous predictors. Fixed,
  adaptive-neighbour and generalized-neighbour bandwidths are supported.
  Fit higher-degree models with \code{bernstein.basis=TRUE} for numerical
  stability; plotting preserves the stored basis and selected degrees.
  Larger auxiliary spaces require more observations in every profile and
  more computation; infeasible or unidentified systems fail explicitly, with
  no automatic degree reduction. Explicit auxiliary vectors can describe
  smaller anisotropic spaces that still contain the main model.
  The original profile-specific construction admits at most one binary
  unordered and one binary ordered predictor, with every joint profile observed.
  It supports means, first partial derivatives and categorical first
  differences. Categorical differences have
  numerical operator validation; coverage evidence must be interpreted for
  the configurations actually studied. Higher derivatives, other kernels,
  additive and tensor main bases remain outside this construction.

  Other categorical layouts use a partially pooled auxiliary for fixed-bandwidth
  means with Gaussian or Mammen calibration. This includes multiple unordered
  and ordered predictors and factors with more than two levels. Higher-order
  continuous curvature is shared across observed joint profiles; each profile
  retains its own polynomial terms through the fitted main degree. P0 and P1
  use this same joint space. Existing profile-specific routes are unchanged;
  no failure silently selects a different construction. Evaluation profiles
  must occur in the training sample. Unused factor levels do not create fitted
  profiles. Local-to-global adaptation is not activated for this pooled path.
  Its qualified workspace contains at most 128 joint auxiliary columns, and
  every profile main space must leave positive residual support. Ill-conditioned
  pilots fail explicitly, even with stable bases; no profile is dropped or
  degree reduced. The joint space, profile sizes and rank criterion are recorded
  in \code{inference$settings$auxiliary_space}. The existing auxiliary controls
  remain effective, including literal raw-bandwidth overrides.

  Bounded pooled-mean simulations with two continuous, two binary unordered and
  two binary ordered predictors, degree-one main fits, and sample sizes 400/800
  gave encouraging overall Gaussian/Mammen coverage under balanced and sparse
  profile counts. This is narrower evidence than numerical support for other
  counts, levels and degrees. Sparse weighted pilot geometry, including the
  retained wage1 stress case, can still fail numerical identification.
  Nearest-neighbor bandwidths, derivatives, categorical contrasts and other
  bootstrap laws remain unqualified for the pooled path.

  Outside the pooled categorical path, higher pure coordinate derivatives use
  the existing \code{gradient.order} argument. Each requested order greater than
  one must satisfy \code{gradient.order[j] <= degree[j] - 1}; the fitted main
  degree is never increased. Scalar orders are recycled across continuous
  coordinates and vectors retain their coordinate-specific meaning. These are
  pure partial derivatives, not cross partials. Newly admitted higher orders
  use Gaussian or Mammen calibration. The complete weighted P1 pilot is
  differentiated analytically, including its weight and coefficient dependence,
  and all correction covariance terms and derivative units are retained.
  The effective vector is recorded in \code{inference$settings$gradient.order}.
  Numerical support does not imply well-determined high derivatives at small
  sample sizes; intervals can be wide and sparse pilot geometry can fail.

  Use \code{errors="asymptotic"} for Gaussian calibration or
  \code{errors="bootstrap"} with the existing \code{bootstrap} and
  \code{\link{np_boot_control}} choices. Mammen is the default wild law;
  explicit Rademacher remains available. In the iterated construction, Gaussian
  and Mammen use a residual-pilot studentized multiplier process. Each draw
  perturbs the fitted P0 pilot by Gaussian or Mammen innovations scaled by the
  observed smoothed variance. The same normalized P0 residual operator and
  variance smoother then provide that draw's denominator. Numerators use the
  complete corrected observation rows. The corrected center and observed
  standard errors are unchanged; replicate-specific denominators calibrate all
  three interval types. No bandwidth or degree is reselected. Rademacher retains
  its fixed-factor multiplier calibration. The iterated inid, fixed and geom
  constructions resample corrected residual scores rather than refitting each
  row resample. Inid uses the smoothed variance refinement; fixed and geom use
  their score-bootstrap scale. Automatic block selection uses raw P0 residuals
  with \code{b.star(..., round=TRUE)}: BstarCB for fixed and BstarSB for geom.
  The legacy \code{nonfixed="frozen"} resampling switch is not applicable.
  Explicit block lengths are retained as overrides.

  \code{band="pointwise"}, \code{"simultaneous"}, \code{"bonferroni"}
  and \code{"all"} share the corrected center and scale. Simultaneous
  coverage refers to the finite displayed grid within each scalar target
  family. Bands are nested and coincide for a singleton family. Gaussian
  and bootstrap marginal limits use absolute studentized draws and an outward
  Bonferroni envelope. \code{band="pmzsd"} retains
  its normal-quantile times standard-error meaning. Use \code{gaussian.draws}
  or \code{B} to control calibration precision. Increasing the number of
  draws can be necessary for remote empirical Bonferroni tails on large grids.

  The original one/two-continuous, LC/LL/quadratic defaults have encouraging
  finite-sample evidence, especially for Gaussian and Mammen calibration.
  Independent confirmation of the studentized construction used fixed
  one-predictor models at sample sizes 100 and 250, degrees 0, 2 and 4, with
  quadratic and sinusoidal signals. Broader mixed-data and nearest-neighbor
  support has numerical validation; this is not evidence of comparable coverage
  throughout those classes. Regularized generalized-neighbor main
  fits with a stored Bernstein basis use their complete canonical observation
  rows, including the numerical ridge and intercept restoration. Polynomial
  correction includes their regularization bias; covariance uses the complete
  corrected rows. This does not establish coverage for unidentified pilot fits
  or a uniform honesty guarantee. These newly admitted regularized rows
  currently support Gaussian calibration only. Bounded fixed-specification
  simulations with Gaussian errors gave adequate coverage, but a two-predictor
  small-neighbor mean case had substantially inadequate Mammen coverage despite
  negligible residual bias and accurate variance. Some corrected rows were
  dominated by one observation. Such concentration also limits justification
  of the Gaussian approximation for non-Gaussian errors; this validation does
  not establish arbitrary-error or uniform coverage. Bootstrap calibration of
  these regularized rows is explicitly refused, without changing the selected
  law. Bootstrap support for previously admitted zero-ridge rows is unchanged.
  Other regularized main rows remain unqualified.
  Per-target/frame numerical ridge values are returned in
  \code{inference$main.regularization}, and \code{inference$covariance}
  describes the variance and calibration construction.
  Auxiliary fits require adequate local
  support throughout the retained training sample, including observations
  outside the displayed grid. Skewed or heavy-tailed predictors can cause
  explicit auxiliary-rank failures; trimming the plot grid does not fix these.
  Inspect the reported pilot role and its bandwidth or degree controls. A
  transformed-predictor model is a new fitted specification, not a plotting
  adjustment. Stable bases also do not guarantee well-conditioned normal
  equations for very wide predictor ranges.

  The profile-specific correction can remove much of the main estimator's
  categorical borrowing and substantially increase interval widths, especially
  in small profiles. Computation includes streamed quadratic-time training
  operators; absence of a training-square matrix does not imply linear time.
  MPI distributes calibration draws, including their residual-pilot projections
  and variance calculations. Observed-data preparation retains its existing
  computational ownership.
  These results are not a uniform honesty guarantee
  over arbitrary regression functions or smoothing choices. Difficult
  localized bumps, residual bias and dependent-data mean coverage remain
  qualifications. Bandwidth selection and degree search uncertainty are not
  included. Overrides change the auxiliary construction and require judgment;
  they are not selected separately for each data-generating process.

  In \pkg{npRmpi}, invoke the plot from the coordinator after initializing
  the pool. Whole-plot \code{mpi.bcast.cmd(...)} execution is unsupported.
  A profile-created worker pool can be adopted with
  \code{npRmpi.init(mode="spawn")}, which reuses the existing communicator;
  the profile's workers continue to execute the bounded calibration tasks.

  For historical wild-bootstrap reproduction, specify both
  \code{np_bias_control(construction="legacy")} and
  \code{np_boot_control(wild="rademacher")}. Unsupported iterated settings
  never silently select the legacy construction.
}
\value{
  A classed list consumed by the plot method. With \code{output="data"} or
  \code{"plot-data"}, each returned regression component includes
  \code{inference}: original and corrected estimates, correction, standard
  errors, all three bands, evaluation grid and target mapping, requested and
  effective controls (including the named auxiliary degree vector, its rule
  and term count in \code{settings$degree_spec}), main-model settings, training-row metadata, calibration
  law and draw count, block selection and RNG state. For Gaussian/Mammen
  studentized calibration, \code{inference$studentization} records the
  construction and bounded training/draw block sizes; \code{inference$covariance}
  records the denominator and calibration methods. The historical \code{bias}
  and \code{gbias} sign conventions are preserved. With \code{band="all"},
  \code{band.errors} also retains the displayed error widths. The original
  fitted object is unchanged.
}
\seealso{\code{\link{plot.np}}, \code{\link{np_boot_control}}, \code{\link{b.star}}}
\examples{
\dontrun{
npRmpi.init(nslaves=2)
plot(fit, errors="asymptotic", center="bias-corrected",
     band="simultaneous", bias.control=np_bias_control(), output="data")
plot(fit, errors="bootstrap", center="bias-corrected", bootstrap="wild",
     band="all", bias.control=np_bias_control(iterations=32L), output="data")
npRmpi.quit()
}
}
\keyword{hplot}
