\name{npscoeftest}
\alias{npscoeftest}
\alias{print.scoeftest}
\alias{summary.scoeftest}
\title{Kernel Specification Test for a Smooth-Coefficient Regression}
\description{Tests whether the conditional mean is affine in the slope regressors at each value of the modifiers. Both the intercept and the slopes may vary freely with the modifiers.}
\usage{
npscoeftest(object, method = c("wide", "jlwy"), B = 399L,
            random.seed = 42L, h = NULL,
            jlwy.bandwidth = c("application", "simulation"), trim = NULL, ...)
\method{print}{scoeftest}(x, ...)
\method{summary}{scoeftest}(object, ...)
}
\arguments{
\item{object}{A training-sample \code{\link{npscoef}} fit or \code{\link{npscoefbw}} object retaining separate continuous X and Z and a finite numeric response. Supply a non-iterated local-constant fit with fixed, unbounded, second-order Gaussian bandwidths. The summary method takes a test result.}
\item{method}{\code{"wide"} uses the supplied null bandwidths, normally selected by cross-validation, and a signed, unstudentized off-diagonal score. \code{"jlwy"} uses the studentized statistic of Jansen, Li, Wang and Yang (2008) and refits the null at their bandwidth rule.}
\item{B}{Number of fixed-design Mammen wild bootstrap draws, at least 9. The default P-value resolution is 1/399.}
\item{random.seed}{Temporary integer seed in the caller's RNG kind. \code{NULL} consumes the current stream.}
\item{h}{Optional positive pair bandwidths in physical units, one per coordinate of X followed by Z. A fully named vector is matched by unique column names. The defaults are given below.}
\item{jlwy.bandwidth}{Bandwidth convention for \code{"jlwy"}. The \code{"application"} setting reproduces the rule of the paper's two-modifier application (equations 15--16). The \code{"simulation"} setting uses the rule from its scalar simulation (Section 3.3). Other dimensions extrapolate these rules. An explicit \code{h} replaces the pair bandwidths.}
\item{trim}{Optional logical vector for \code{"jlwy"}. TRUE retains an observation in the pair sums. Use a fixed trimming region in covariate space. All observations still enter each null fit. The default retains every observation.}
\item{x}{A test result.}
\item{...}{Unused arguments to the test produce an error.}
}
\details{
Write \eqn{W=(X,Z)} and let \eqn{e_i} denote the residuals under the smooth-coefficient null. With the product Gaussian \eqn{K_{ij}=\exp(-\sum_k[(W_{ik}-W_{jk})/h_k]^2/2)}, the default reports \eqn{Q=\sum_{i\ne j} K_{ij} e_i e_j}. This score can have either sign. A large positive value counts against the null. Each coordinate uses its own sample standard deviation \eqn{s_k}. The default pair widths are \eqn{h_k=\sqrt{2}\,s_k n^{-1/5}}. The supplied null bandwidths are held fixed. A manually specified input uses those manual bandwidths.

The JLWY statistic divides Q by \eqn{[2\sum_{i\ne j}K_{ij}^2 e_i^2 e_j^2]^{1/2}} and recomputes this denominator for every bootstrap response. If trimming is supplied both sums use the retained pairs. Its null widths are \eqn{s_k n^{-1/6}} on Z. Its pair widths on W are \eqn{s_k n^{-1/6}} for the application and \eqn{s_k n^{-1/5}} for the simulation. Positive kernel normalizing constants cancel under studentization.

Both methods generate responses from the fitted null mean plus the observed null residual times an independent Mammen multiplier. Each bootstrap response is refitted on the original predictors at the same null bandwidths. The P value is the proportion of bootstrap scores at least as large as the observed score. To use the paper's strict empirical-percentile decision, compare the observed statistic to the corresponding returned critical value.

The wider unstudentized test is a modification of JLWY supported by finite-sample simulation. At n = 240 with two slopes and one, two or three modifiers, size at the nominal 5 percent level was 7.0, 5.5 and 8.6 percent respectively (128 null samples each). More slopes or modifiers can markedly reduce power. The paper's asymptotic results concern its studentized statistic under its assumptions. The wild bootstrap assumes an appropriate conditional error structure. Local moment systems must be numerically nonsingular.}
\value{An object of class \code{scoeftest} with the signed observed \code{statistic}, \code{bootstrap} scores, \code{P}, \code{exceedances}, \code{critical.values} (empirical quantiles of type 1), method, null and pair bandwidths, Monte Carlo seed and draw count, retained sample count, effective degrees of freedom, elapsed seconds and \code{rank.counts} (owned score panels). Normally the raw scores are in squared response units. When this conversion would overflow, underflow or change an exceedance, multiply all returned scores by \eqn{\exp(L)} with L given by \code{log.multiplier} to recover original units. The multiplier is common and positive.}
\references{Jansen, D. W., Li, Q., Wang, Z. and Yang, J. (2008). Fiscal policy and asset markets: A semiparametric analysis. \emph{Journal of Econometrics}, 147, 141--150. \doi{10.1016/j.jeconom.2008.09.007}.}
\seealso{\code{\link{npscoef}}, \code{\link{npscoefbw}}, \code{\link{nppltest}}, \code{\link{npindextest}}}
\examples{
\dontrun{
npRmpi.init(nslaves = 2L)
set.seed(42)
n <- 240L
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z = runif(n, -1, 1))
# Null: both slopes vary with z. Most null samples should not reject.
d$y <- with(d, sin(z) + x1 * cos(z) + x2 * z + rnorm(n))
bw <- npscoefbw(y ~ x1 + x2 | z, data = d)
summary(npscoeftest(bw))
# The paper's application-bandwidth version, on the same data.
summary(npscoeftest(bw, method = "jlwy"))
# Alternative: an interaction between the slopes violates the null.
# Strong signals tend to give small P values. Results vary across seeds.
d$y <- with(d, y + 1.5 * x1 * x2)
summary(npscoeftest(npscoef(y ~ x1 + x2 | z, data = d)))
npRmpi.quit()
}
}
