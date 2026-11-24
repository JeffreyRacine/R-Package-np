\name{npindextest}
\alias{npindextest}
\alias{indextest}
\alias{print.indextest}
\alias{summary.indextest}
\title{Bootstrap Specification Test for a Single-Index Model}
\description{
Tests whether the conditional mean of a continuous response, or the success
probability of a binary response, depends on the covariates through one
unknown linear index with an unrestricted link (Maistre and Patilea, 2019).
}
\usage{
npindextest(object, bw.method = c("published", "maximum", "separate"),
            B = 399L, random.seed = 42L)
\method{print}{indextest}(x, ...)
\method{summary}{indextest}(object, ...)
}
\arguments{
  \item{object}{An \code{\link{npindex}} fit or \code{\link{npindexbw}}
    object obtained by joint estimation of the direction and link bandwidth
    on the retained training data. The default index fit is local constant.
    Local-linear and prespecified local-polynomial fits use the degree and
    basis stored in the supplied object. A fitted object with replacement training
    data must reproduce the stored search criterion. For \code{summary}, an \code{indextest}
    result.}
  \item{bw.method}{Test bandwidth choice. \code{"published"} uses the
    authors' empirical factor one. \code{"maximum"} jointly calibrates three
    factors. \code{"separate"} sets a separate residual bandwidth and a
    narrower pairwise bandwidth. Choose one procedure before inspecting its
    p-value.}
  \item{B}{Number of bootstrap draws, an integer of at least 9.}
  \item{random.seed}{A nonnegative integer for temporary seeding. The
    caller's RNG state is restored on exit. \code{NULL} advances the current
    stream. See \code{\link{np_rng_scope}}.}
  \item{x}{An \code{indextest} result.}
  \item{...}{Further arguments to the print or summary method.}
}
\details{
Use independent observations, at least two continuous numeric covariates of
full centered rank, and a fixed, unbounded second-order Gaussian
index fit. Klein--Spady requires both numeric response values 0
and 1. The hypothesis concerns existence of a single-index representation.
Off-index curvature and interactions are examples of departures. No link
function is specified under the null.

The statistic combines off-diagonal products of density-weighted,
leave-one-out residuals with Gaussian weights in the estimated index and
its orthogonal complement. All test bandwidths are on the standardized,
unit-direction index. The published convention sets \eqn{h=n^{-2/9}} and
uses the fitted CV link bandwidth for the residuals of a local-constant fit. The joint maximum uses
the largest studentized statistic over factors 0.5, 1 and 2 and recomputes
that maximum on every bootstrap draw. The experimental separate-bandwidth
version uses \eqn{g=n^{-2/9}} for residuals and \eqn{h=n^{-17/36}} for pairs.
These last two choices extend the paper's empirical procedure. The separate
pair was conservative and had lower power in small two-covariate simulations.

For a local-linear or fixed-degree polynomial index fit, each draw keeps
the fitted degree, basis and Bernstein setting. These fits extend the
procedure of Maistre and Patilea. On the observed standardized index,
a separate Gaussian local-constant least-squares CV bandwidth supplies the
pilot, including valid Bernoulli probabilities. With \code{"published"} or
\code{"maximum"} it also supplies the residual bandwidth and is reselected
on each refitted index. With \code{"separate"} the residual bandwidth follows
the fixed rule above. The polynomial degree must be prespecified.

In simulations at \eqn{n=192} with two continuous covariates and fixed
Bernstein degrees 0 through 5, higher degrees often reduced power against
off-index interactions. On one curved Klein--Spady null, degrees 2 and 3
rejected 12 and 13 times in 128 replications at the 5 percent level.
Degrees outside this range were not studied.

Every draw refits the direction and CV bandwidth using the original model's
search settings. Ichimura uses leave-one-out CV fitted means and Mammen wild
weights. The Klein--Spady adaptation draws Bernoulli responses from inclusive
CV fitted probabilities. The one-sided p-value is the proportion of bootstrap statistics
greater than or equal to the observed statistic. Predictor geometry is reusable, while the
estimated direction and its kernel weights depend on each bootstrap response.
For Ichimura fits, a common change of response units cancels from the
studentized score at a given direction and bandwidth. The calculation rescales
the response before forming pair differences. If refitting in different units,
the starting direction and absolute stopping tolerance of \code{\link{npindexbw}}
can affect the selected fit and hence the bootstrap results. For very small
or large units, fit and test in numerically convenient response units.

For local-constant bootstrap Nelder--Mead searches with default convergence tolerances,
at least 160 observations and automatic or disabled trees, Gaussian
cross-validation uses interpolation and FFT convolution on the trial index.
Sparse rows can be summed directly. Extreme response scales, unsuitable
grids or insufficient donor mass use the native criterion. Other optimizer
settings and explicitly enabled trees use the native evaluator throughout.
The search checks the winning objective with the native calculation.
A full native search repeats the original starts and seed if the
discrepancy exceeds \eqn{2\times10^{-10}\max(1,|CV|)}, or if the two best
final start values are within this tolerance. Trial rounding near competing minima
can still change individual draws. In a paired study with two covariates and
240 observations, the two evaluators gave identical 5 percent decisions
on 3,360 datasets across the two estimators, ordinary and weak-index nulls,
off-index curvature and interaction alternatives.

The pilot and statistic process pairs in blocks. For a fixed number of
covariates their workspace grows linearly with the sample size.

The fit and test should use the same package defaults. Some explicit
optimizer controls or starting values survive only as expressions in the
original call. If their values cannot be recovered, fit with literal values
or use \code{do.call(npindexbw, args)} with evaluated controls.
}
\value{
  \code{total.time} gives wall-clock seconds for this test call, including
  preparation and bootstrap calibration. The summary reports this time.
  With \code{getOption("show.signif.stars")} enabled, the printed bootstrap
  p-value carries the significance codes used by \code{\link{npsigtest}}.
  Each symbol requires a p-value strictly below its threshold. A zero is displayed as \samp{0 (0 of B)}, with the actual number of draws.

An object of class \code{indextest}. \code{statistic}, \code{P},
\code{exceedances}, \code{critical.values} and \code{bootstrap} report the
observed statistic, p-value, exceedance count, upper critical values and
ordered bootstrap statistics. \code{components} gives the observed values
for \code{h} in the same order. \code{g.cv}, \code{g.pilot} and \code{g.test} are the
observed model, bootstrap pilot and score residual bandwidths on the
standardized index. \code{beta}, \code{method},
\code{bw.method}, \code{n}, \code{refit.controls}, \code{boot.num} and
\code{random.seed} describe the fit and calibration. \code{call} and
\code{model.call} retain the calls. \code{cv.counts} records the numbers of grid and native bootstrap CV
trial evaluations and complete native search replays. Evaluations used for
final certification are excluded. \code{rank.counts} is NULL in serial
execution. Print and summary return the result invisibly.
}
\references{
Maistre, S. and Patilea, V. (2019). Nonparametric model checks of
single-index assumptions. \emph{Statistica Sinica}, 29, 113--138.
\url{https://www3.stat.sinica.edu.tw/statistica/oldpdf/A29n17.pdf}.

Mammen, E. (1993). Bootstrap and wild bootstrap for high dimensional linear
models. \emph{Annals of Statistics}, 21, 255--285.
}
\seealso{\code{\link{npindex}}, \code{\link{npindexbw}},
\code{\link{nppltest}}, \code{\link{npcmstest}}}
\examples{
\dontrun{
# set.seed() controls the simulated sample. Each test uses the published bandwidth
# choice and the default bootstrap seed, random.seed = 42.
set.seed(42)
n <- 160L
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
z <- with(d, x1 + 0.6*x2)
w <- with(d, x2 - 0.6*x1)

# Ichimura null. The mean depends only on z. The variance may depend on w.
d$y0 <- z + 0.45*sin(z) + (0.6 + 0.2*abs(w))*rnorm(n)
# Ichimura alternative. Keep the same noise and add an off-index interaction.
d$y1 <- d$y0 + 1.8*z*w

# Klein-Spady alternative. The interaction in the log odds makes the
# success probability depend on both z and its orthogonal direction w.
d$choice1 <- rbinom(n, 1, plogis(z + 1.8*z*w))
# Klein-Spady null. The success probability depends only on z.
d$choice0 <- rbinom(n, 1, plogis(z))

bw0 <- npindexbw(y0 ~ x1 + x2, data = d, method = "ichimura")
summary(npindextest(npindex(bw0)))
# P Value: 0.7619048 (304 of 399) at seed 42. Do not reject the single-index null at level 0.05.

bw1 <- npindexbw(y1 ~ x1 + x2, data = d, method = "ichimura")
summary(npindextest(npindex(bw1)))
# P Value: 0 (0 of 399) at seed 42. Reject the single-index null at level 0.05.

binary0 <- npindexbw(choice0 ~ x1 + x2, data = d, method = "kleinspady")
summary(npindextest(binary0))
# P Value: 0.06766917 (27 of 399) at seed 42. Do not reject the single-index null at level 0.05.

binary1 <- npindexbw(choice1 ~ x1 + x2, data = d, method = "kleinspady")
summary(npindextest(binary1))
# P Value: 0 (0 of 399) at seed 42. Reject the single-index null at level 0.05.

# A reported zero means none of the 399 bootstrap draws reached the observed statistic.
# A different data seed gives another sample. At level 0.05 a true null
# can be rejected, and a departure can go undetected. Changing only
# random.seed varies the bootstrap p-value for the same observations.
} % enddontrun
}
\keyword{htest}
