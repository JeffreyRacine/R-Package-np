\name{np_conditional_bias_control}
\alias{np_conditional_bias_control}
\title{Experimental Bias-Corrected Conditional Density and Distribution Bands}
\description{Selects explicitly the experimental local-constant conditional
inference constructions. Existing plots and methods are unchanged when this
control is omitted.}
\usage{
np_conditional_bias_control(construction = c("analytic", "jackknife", "jackknife-mp1"))
}
\arguments{
\item{construction}{The corrected functional. Use \code{"analytic"} for
conditional CDF inference and conditional density bootstrap inference;
use \code{"jackknife"} for conditional density asymptotic inference.
Use \code{"jackknife-mp1"} for the separate asymmetric studentized
conditional density bootstrap described below.}
}
\value{A validated control object for the \code{bias.control} argument of
\code{plot} applied to \code{npcdist} or \code{npcdens} fits or bandwidths.}
\details{
These are experimental opt-ins for local constant, one continuous explanatory
variable and one continuous response, fixed unbounded second-order Gaussian
kernels. Gradients, quantiles, shape projections, nonfixed bandwidths and mixed
variables are not admitted. Original fits and selected bandwidths remain fixed.

Specify \code{center="bias-corrected"} and an explicit \code{band}:
\code{"pointwise"}, \code{"simultaneous"}, \code{"bonferroni"} or \code{"all"}.
With \code{errors="bootstrap"}, explicitly specify
\code{bootstrap="conditional-pilot"} and
\code{np_boot_control(calibration="symmetric-sd")} or
\code{np_boot_control(calibration="asymmetric-sd")}.
This law keeps explanatory observations fixed and draws each response from its
conditional Gaussian donor mixture, with pilot bandwidths multiplied by
\eqn{\sqrt{2}}. It is not ordinary pairs resampling and does not replace
\code{inid}, \code{fixed}, \code{geom}, residual or wild resampling.

The analytic correction contains both explanatory-weight and response-kernel
curvature. Each bootstrap root includes the complete corrected functional,
centered by its conditional-pilot expectation. A common marginal bootstrap SD,
maximized within each explanatory slice, scales the simultaneous signed or
absolute extremes. The asymmetric method inverts the signed extremes; it does
not use two unrelated lower/upper SD estimates.

With \code{errors="asymptotic"}, complete conditional-pilot covariances feed
Gaussian finite-grid calibration. For CDFs the analytic correction and slice
scales are used. For densities the corrected estimate is twice the estimate at
the fitted bandwidths minus the estimate at bandwidths enlarged by
\eqn{\sqrt{2}} in both coordinates. All cross-covariance terms are retained.
The density asymptotic and bootstrap corrections therefore differ deliberately.

The separate density construction \code{"jackknife-mp1"} requires
\code{errors="bootstrap"}, \code{bootstrap="conditional-pilot"} and
\code{calibration="asymmetric-sd"} in \code{np_boot_control()}.
It uses the same jackknife corrected functional as density asymptotic inference.
Each complete corrected bootstrap root is divided by its own MP1 standard
error. MP1 averages the original and resample pilot covariance matrices and
adds one quarter of the sum of per-observation mean-difference outer products.
This is the covariance of independent equal pilot mixtures for each training
observation, retaining every correction cross term. Its variance is at least
half the original pilot variance; this is an algebraic property, not a variance
floor. Nonpositive or nonfinite component variances cause an error.
Signed studentized roots are inverted using the original standard errors.
The selector is currently restricted to asymmetric calibration; it does not
change either existing analytic SD bootstrap construction.

A nominal marginal sharp-score safeguard uses the certified range of each
signed contribution and its heterogeneous variance bound. It augments the
Gaussian or bootstrap intervals. The density asymptotic and jackknife-MP1 constructions use a
marginal hull; the other constructions also propagate required widths within
explanatory slices. Companion hulls protect interval nesting. These are
approximate confidence constructions, not finite-sample honesty guarantees.
No clipping, positivity projection or forced properness is applied to corrected
centers or endpoints.

The family consists of the actual requested finite evaluation points within
one panel or surface. Separate panels are not a joint family; a finite grid is
not the entire continuous domain. Initial admission is at most 256 target points
and 25 distinct response targets. Choose, for example, \code{neval=9} for an
81-point surface; the ordinary default 50-by-50 surface exceeds this new method's
scope. Inputs exceeding the numerical workspace budget of 256 MiB fail before
resampling; the grid and draw count are never silently reduced. The budget is
not a bound on total R-process memory. Bootstrap \code{B} and asymptotic
\code{gaussian.draws} must be integers at least two; defaults are 9999.
Low finite-draw tail resolution is reported. Alpha is propagated throughout
in \eqn{(0,0.5)}; population coverage evidence is currently for nominal 95\% only.

The returned object's \code{inference} component records absolute endpoints for
all three bands, correction, pilot law, calibration, physical bandwidths,
requested grid, scales and covariance provenance. Jackknife-MP1 also records
its construction, studentization, mixture weight and attained scale-ratio range. For these new routes,
\code{bias.corrected + conderr[,1]} and
\code{bias.corrected + conderr[,2]} recover the selected absolute endpoints.
For \code{band="all"}, \code{conderr} represents the simultaneous band;
all endpoints remain in \code{inference$bands}. Corrected centers may be signed.

Saved Gaussian-design studies support full-grid CDF coverage at n=100 and 400,
and density coverage at n=400 and 800. CDF plot subfamilies and marginal intervals
can be conservative above 99\%. Density bootstrap small-n studies failed lower
coverage gates. These failures remain limitations; they do not establish
uniform coverage, robustness to arbitrary DGPs, or adjustment for bandwidth
selection uncertainty. Fresh paired Gaussian studies of asymmetric MP1 at n=400 and 800 retained
481/512 and 488/512 simultaneous coverages, respectively, on the studied
81-point random grids. These results qualify this experimental opt-in, not
a general coverage or optimality guarantee. The options are intended for
further experimentation.
}
\examples{
\dontrun{
set.seed(42)
x <- rnorm(100)
y <- 0.5*x + rnorm(100)
model <- npcdist(y ~ x, regtype="lc")
p <- plot(model, errors="asymptotic", center="bias-corrected",
          bias.control=np_conditional_bias_control("analytic"),
          band="all", neval=9, output="plot-data")
p <- plot(model, errors="bootstrap", bootstrap="conditional-pilot",
          center="bias-corrected", band="all", neval=9,
          bias.control=np_conditional_bias_control("analytic"),
          boot.control=np_boot_control(calibration="asymmetric-sd"),
          output="plot-data")

model.d <- npcdens(y ~ x, regtype="lc")
p.d <- plot(model.d, errors="bootstrap", bootstrap="conditional-pilot",
            center="bias-corrected", band="all", neval=9,
            bias.control=np_conditional_bias_control("jackknife-mp1"),
            boot.control=np_boot_control(calibration="asymmetric-sd"),
            output="plot-data")
}
}
\seealso{\code{\link{plot.np}}, \code{\link{np_boot_control}},
\code{\link{npcdist}}, \code{\link{npcdens}}}
\keyword{nonparametric}
