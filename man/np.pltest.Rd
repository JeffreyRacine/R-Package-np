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
         ..., boot.correction = c("none", "cv-pilot", "pilot-hc2", "cv-raw"))
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
  \item{B}{Integer number of draws, at least 9. The smallest attainable
    p-value is \eqn{1/(B+1)}.}
  \item{random.seed}{An integer for temporary seeding. The caller's RNG
    state is restored on return or error. \code{NULL} uses the current stream
    and allows it to advance. A single seed covers bandwidth search and
    bootstrap in sequence. Changing a randomized search may change the
    multiplier stream as well as the selected bandwidths. Temporary seeding
    requires a restorable built-in RNG, excluding Box--Muller normals and
    user-supplied generators. A manual MPI collective also requires matching,
    restorable generator kinds when the seed is NULL. All its ranks then
    adopt the master's current stream and leave with its advanced state. See \code{\link{np_rng_scope}} for supported
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
}
\details{
  Sun, Chen, Liang and Ruppert (2021) integrate a squared marked residual
  process over directions of the joint predictors and empirical thresholds.
  The angle calculation integrates directions analytically. Kernel smoothing
  estimates the null conditional means. A bootstrap response is its chosen
  center plus a residual amplitude times an independent multiplier. Each
  draw is refitted with the test bandwidths. The test forms the training
  operator from shared fitted smoother weights for continuous Sun nuisance variables
  or from public unit-response fits for other configurations, preserving the fitted
  estimator's categorical kernel conventions.

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

  Computation stores dense training matrices and takes cubic time for the
  angular calculation at fixed predictor dimension. The continuous Sun test reuses
  one training smoother and the residualized linear projection across all draws.
  Other configurations construct the operator using public unit-response fits.
  The CV-pilot HC2 option constructs a second such operator. In \pkg{npRmpi},
  rank zero calculates the angular matrix. In a collective session it shares
  this matrix with the other ranks. To use MPI,
  initialize an \pkg{npRmpi} session or attach pool and call the test on
  the master. Under the package startup profile, make the same training data
  and object available on every rank and use, for example,
  \code{mpi.bcast.cmd(ans <- nppltest(fit), caller.execute = TRUE)}.
  Every rank returns the result. To print one summary, print it on rank zero.
  Workers compute disjoint sets of replications. The returned vector is in
  serial replication order. Floating-point reductions in the fitted model
  can produce small numerical differences between package runs.
}
\value{
  \code{total.time} gives wall-clock seconds for this test call, including
  preparation and bootstrap calibration. The summary reports this time.
  With \code{getOption("show.signif.stars")} enabled, the printed bootstrap
  p-value carries the significance codes used by \code{\link{npsigtest}}.
  Each symbol requires a p-value strictly below its threshold.

  A \code{pltest} list. \code{statistic} is the observed quadratic form.
  \code{P} is one plus the number of bootstrap statistics greater than or
  equal to the observed statistic, divided by \eqn{B+1}.
  \code{exceedances}, \code{bootstrap} and \code{critical.values}
  retain the count, bootstrap statistics in replication order and empirical 90, 95 and 99 percent
  quantiles. \code{bw.method}, \code{bws}, \code{bandwidths},
  \code{boot.method}, \code{boot.correction} and \code{pilot} identify the
  calibration. \code{n}, \code{continuous.z}, \code{cell.counts},
  \code{projection.scale}, \code{boot.num} and \code{random.seed}
  describe the sample and computation. \code{call} and \code{model.call}
  retain the test and input model calls. \code{x.names}, \code{z.names} and
  \code{factor.levels} identify the predictors and nuisance factor levels.
  In MPI, \code{rank.counts} is a named integer vector counting draws
  computed by each participating rank. It is \code{NULL} in serial R.
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
# Run in an interactive master session.
force.run <- nzchar(Sys.getenv("NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK"))

in.check <- !force.run && (
            nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
            nzchar(Sys.getenv("_R_CHECK_INTERNALS2_")) ||
            nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_")) ||
            nzchar(Sys.getenv("_R_CHECK_EXAMPLE_TIMING_THRESHOLD_"))
          )

if (!in.check) {
  npRmpi.init(nslaves = 1L)
  tryCatch({
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
    # p = 0.4875 at seed 42. Do not reject partial linearity at level 0.05.

    # For the alternative, keep the predictors and noise and add curvature in x1.
    # The nuisance function of z cannot absorb this term.
    d$y1 <- with(d, y0 + 1.5*(x1^2 - 1))
    bw1 <- npplregbw(y1 ~ x1 + x2 | z, data = d, regtype = "ll")
    summary(nppltest(bw1))
    # p = 0.0025 at seed 42. Reject partial linearity at level 0.05.

    # With another data seed, even a true null can be rejected, and a sample
    # drawn under the alternative can fail to reject. With the data held fixed,
    # changing random.seed can change the p-value through bootstrap variation.
  }, finally = npRmpi.quit())
}
} % enddontrun
}
\keyword{htest}
