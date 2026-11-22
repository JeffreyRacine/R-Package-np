\name{nppltest}
\alias{nppltest}
\alias{print.pltest}
\alias{summary.pltest}
\title{Projection Test of Partial Linearity}
\description{
  Tests whether the conditional mean has the form
  \eqn{E(Y\mid X,Z)=X^T\beta+g(Z)} against an unrestricted conditional
  mean. A rejection can indicate curvature in the linear regressors or an
  interaction between the linear and nuisance regressors.
}
\usage{
nppltest(object, bw.method = c("sun", "cv", "object"),
         boot.method = c("wild-rademacher", "wild"), B = 399L,
         random.seed = 42L, projection.scale = c("original", "standardized"),
         ..., boot.correction = c("none", "cv-pilot", "pilot-hc2", "cv-raw"),
         integration = c("sobol", "exact"), directions = 2048L)
}
\arguments{
  \item{object}{A fitted \code{\link{npplreg}} or \code{\link{npplregbw}}
    object with recoverable training data. The response must be numeric.
    The test uses the aligned complete training observations, including when
    the fitted object contains separate evaluation data.}
  \item{bw.method}{\code{"sun"} uses local-linear undersmoothing for the
    continuous nuisance coordinates and exact matching on nuisance factors.
    \code{"cv"} selects least-squares CV bandwidths separately for the
    conditional means of the response and each linear regressor.
    \code{"object"} uses the bandwidths and regression specification of
    \code{object}. The selected test bandwidths are held fixed across draws.}
  \item{boot.method}{Wild multiplier law. \code{"wild-rademacher"} uses
    independent equiprobable signs. \code{"wild"} uses Mammen's asymmetric
    two-point law with mean zero and unit variance.}
  \item{B}{Integer number of draws, at least 9. Ties count as exceedances.}
  \item{random.seed}{An integer for temporary seeding. The caller's RNG
    state is restored on return or error. \code{NULL} uses the current stream
    and allows it to advance. A single seed covers bandwidth search and
    bootstrap in sequence. Changing a randomized search may change the
    multiplier stream as well as the selected bandwidths. Temporary seeding
    requires a restorable built-in RNG, excluding Box--Muller normals and
    user-supplied generators. See \code{\link{np_rng_scope}} for supported
    temporary seed states.}
  \item{projection.scale}{Use original predictor units or standardize numeric
    coordinates and ordered-factor scores for the projection geometry.
    In the null fit, linear factors use the slope coding of
    \code{\link{npplreg}}. In the projection, unordered levels are encoded so two different levels are distance one.
    Scaling can affect power.}
  \item{\dots}{Named controls for \code{bw.method = "cv"} only:
    \code{bwtype}, \code{bwscaling}, \code{ckertype}, \code{ckerorder},
    \code{ckerbound}, \code{ckerlb}, \code{ckerub}, \code{ukertype},
    \code{okertype}, \code{nmulti}, \code{itmax}, \code{ftol},
    \code{tol}, \code{small}, \code{powell.remin}, and
    \code{scale.factor.search.lower}. These control the test smoother.
    The separate CV bootstrap pilot uses fixed bandwidths and one start.}
  \item{boot.correction}{\code{"none"} uses the fitted null mean and its raw
    residuals. \code{"cv-pilot"} uses a separate CV partially linear fit for
    the center and divides its residuals by their residual-operator row
    norms (generalized HC2). \code{"cv-raw"} uses the same pilot and its raw
    residuals. Select \code{boot.method = "wild"} as well to reproduce the
    raw CV/Mammen recipe. \code{"pilot-hc2"} uses a leave-one-out-selected
    restricted polynomial pilot of degree at most three and its HC2 residuals.}
  \item{integration}{The default, \code{"sobol"}, averages over a fixed,
    shifted Sobol sequence of projection directions. Select \code{"exact"}
    for the analytical angular integration of Sun et al. (2021).}
  \item{directions}{Positive integer number of Sobol directions. The default
    is 2,048. Smaller counts reduce projection time. Larger counts extend
    the same sequence. Omit this argument with \code{integration = "exact"}.}
}
\details{
  Sun, Chen, Liang and Ruppert (2021) integrate a squared marked residual
  process over directions of the joint predictors and empirical thresholds.
  In Sobol mode, the package shifts each coordinate of the Sobol points
  modulo one, applies the standard normal quantile function and normalizes
  each vector to unit length. The fixed shift was generated independently
  of the data and bootstrap (Mersenne--Twister seed 514907). The observed
  and bootstrap statistics use identical directions and empirical thresholds.
  Tied projected values enter the residual sum together. Unordered factors
  use orthonormal contrasts scaled to unit distance between levels. Sobol
  integration supports up to 1,111 projection coordinates, counting one per
  numeric or ordered predictor and one fewer than the number of observed
  levels per unordered factor. Exact integration uses the original angular
  formula. Kernel smoothing
  estimates the null conditional means. A bootstrap response is its chosen
  center plus a residual amplitude times an independent multiplier. Each
  draw is refitted with the test bandwidths. The test forms the training
  operator from the conditional smoother weights and a single factorization of
  the residualized linear covariates. Factor nuisance variables retain the
  categorical kernel conventions of the fitted regression.

  For a single continuous nuisance variable the default Gaussian bandwidth
  is \eqn{1.06\min(s,0.75\,IQR)n^{-1/3}}. For \eqn{q=2,3} the
  test uses coordinate-specific scales and exponent \eqn{-3/(q+8)}.
  With factors, the default fits within observed joint cells. An all-factor
  nuisance uses cell means. Zero scales, unidentified local designs, absent
  residual variation and undefined requested HC2 adjustments produce errors.
  For \eqn{q\ge4} the proposed local-linear Sun rule is unavailable.

  The published theory assumes a scalar continuous nuisance. The
  multivariate rate, factor encoding, matching and pilot calibrations are
  research extensions. The sampling unit is an independent observation.
  In the independent twelve-cell study the original bootstrap rejected
  15.5 percent of true nulls in the four-slope, binary-factor cell at
  \eqn{n=120}. CV-pilot HC2 and raw CV/Mammen were conservative across
  those designs; nominal power also differed. The test returns no automatic
  choice based on dimension or a preliminary test.

  Finite-direction integration is a computational approximation to the
  published statistic. In two mixed continuous/binary designs at
  \eqn{n=240} and \eqn{800}, with 1,000 samples per design and size and
  399 draws, the largest observed differences from exact integration at
  2,048 directions were 0.4 percentage points for size and 0.5 for power
  against curvature or interactions, across fixed and freshly shifted
  sequences. All paired 95 percent intervals were within two percentage
  points. Exact rejection rates at nominal five percent were 8.1 to 12.7
  percent on these designs. Other designs can give different calibration.

  Sobol integration sorts the projected observations and scores the
  bootstrap residuals in blocks. With a fixed number of directions and
  draws, sorting for one block takes order \eqn{n\log n}. The null fit
  still stores a dense training operator and reuses it across bootstrap
  responses. Exact angular integration stores an additional dense matrix
  and takes cubic time at fixed predictor dimension. The CV-pilot HC2 option also computes the
  squared row norms of the separate pilot residual operator.
}
\value{
  \code{total.time} gives wall-clock seconds for this test call, including
  preparation and bootstrap calibration. The summary reports this time.
  With \code{getOption("show.signif.stars")} enabled, the printed bootstrap
  p-value carries the significance codes used by \code{\link{npsigtest}}.
  Each symbol requires a p-value strictly below its threshold. A zero is displayed as \samp{0 (0 of B)}, with the actual number of draws.

  A \code{pltest} list. \code{statistic} is the observed quadratic form.
  \code{P} is the proportion of bootstrap statistics greater than or
  equal to the observed statistic.
  \code{exceedances}, \code{bootstrap} and \code{critical.values}
  retain the count, bootstrap statistics in replication order and empirical 90, 95 and 99 percent
  quantiles. \code{bw.method}, \code{bws}, \code{bandwidths},
  \code{boot.method}, \code{boot.correction} and \code{pilot} identify the
  calibration. \code{n}, \code{continuous.z}, \code{cell.counts},
  \code{projection.scale}, \code{integration}, \code{directions} (\code{NULL} for exact integration),
  \code{boot.num} and \code{random.seed}
  describe the sample and computation. \code{call} and \code{model.call}
  retain the test and input model calls. \code{x.names}, \code{z.names} and
  \code{factor.levels} identify the predictors and nuisance factor levels.
  \code{rank.counts} is NULL for serial execution.
  \code{print} and \code{summary} return the object invisibly.
}
\section{Comparisons in the literature}{
  Section 6.1 of Sun et al. directly compares Fan and Li's kernel
  U-statistic, Zhu and Ng's residual empirical process, and Xia's adaptive
  projection test. Sun et al. report the highest empirical power in their
  two-regressor interaction design. With 20 linear regressors, two peer
  statistics frequently degenerated and Xia's had low power. These are
  comparisons under the paper's designs and calibrations. They do not
  establish a mixed-data ranking.
}
\references{
  Sun, Z., Chen, F., Liang, H. and Ruppert, D. (2021).
  A projection-based consistent test incorporating dimension-reduction in
  partially linear models. \emph{Statistica Sinica}, 31, 1489--1508.
  \doi{10.5705/ss.202019.0238}.

  Fan, Y. and Li, Q. (1996). Consistent model specification tests:
  omitted variables and semiparametric functional forms.
  \emph{Econometrica}, 64, 865--890. \doi{10.2307/2171848}.

  Zhu, L.-X. and Ng, K. W. (2003). Checking the adequacy of a partial
  linear model. \emph{Statistica Sinica}, 13, 763--781.

  Xia, Y. (2009). Model checking in regression via dimension reduction.
  \emph{Biometrika}, 96, 133--148. \doi{10.1093/biomet/asn074}.

  Mammen, E. (1993). Bootstrap and wild bootstrap for high dimensional
  linear models. \emph{Annals of Statistics}, 21, 255--285.
  \doi{10.1214/aos/1176349025}.
}
\seealso{\code{\link{npplreg}}, \code{\link{npplregbw}}, \code{\link{npcmstest}}}
\examples{
\dontrun{
# set.seed() controls the simulated data. The test's random.seed argument
# controls the bootstrap draws and defaults to 42 in both calls below.
set.seed(42)
n <- 120L
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z = runif(n, -1, 1))
noise <- with(d, (0.45 + 0.30*z^2)*rnorm(n))

# Under the null, x1 and x2 have constant slopes and z enters smoothly.
d$y0 <- with(d, 0.8*x1 - 0.6*x2 + sin(pi*z) + noise)
bw0 <- npplregbw(y0 ~ x1 + x2 | z, data = d, regtype = "ll")
fit0 <- npplreg(bw0)
summary(nppltest(fit0))
# P Value: 0.4862155 (194 of 399) at seed 42. Do not reject partial linearity at level 0.05.

# For the alternative, keep the predictors and noise and add curvature in x1.
# The nuisance function of z cannot absorb this term.
d$y1 <- with(d, y0 + 1.5*(x1^2 - 1))
bw1 <- npplregbw(y1 ~ x1 + x2 | z, data = d, regtype = "ll")
summary(nppltest(bw1))
# P Value: 0 (0 of 399) at seed 42. Reject partial linearity at level 0.05.

# To change the integration, use nppltest(fit0, directions = 512L)
# or nppltest(fit0, integration = "exact").

# With another data seed, even a true null can be rejected, and a sample
# drawn under the alternative can fail to reject. With the data held fixed,
# changing random.seed can change the p-value through bootstrap variation.
} % enddontrun
}
\keyword{htest}
