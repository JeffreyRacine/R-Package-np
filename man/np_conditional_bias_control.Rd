\name{np_conditional_bias_control}
\alias{np_conditional_bias_control}
\title{Experimental Bias-Corrected Conditional Density and Distribution Bands}
\description{
\strong{Bias-corrected bounds are experimental (pre-beta research), incomplete,
and not recommended for general-purpose inference.} See the experimental-status
section below.
Selects explicitly the experimental local-constant conditional
inference constructions. Existing plots and methods are unchanged when this
control is omitted.}
\usage{
np_conditional_bias_control(construction = c("analytic", "jackknife",
  "jackknife-mp1", "pairs-score"), pointwise = NULL)
}
\arguments{
\item{construction}{The corrected functional. Use \code{"analytic"} for
conditional CDF inference and conditional density bootstrap inference;
use \code{"jackknife"} for conditional density asymptotic inference.
Use \code{"jackknife-mp1"} for the separate asymmetric studentized
conditional density bootstrap described below. Use \code{"pairs-score"} for the conditional CDF pairs bootstrap.}
\item{pointwise}{For pairs-score intervals, the default uses sharp-score pointwise protection.
With fixed Gaussian conditioning and a second-order Gaussian, Epanechnikov or
uniform response kernel, or a beta2 response on known finite support, explicitly
select \code{"local-wilson"} to use a positive local-probability Wilson
companion. The jackknife center and basic pairs interval are included.
The common simultaneous band and Bonferroni endpoints expand where required
for pointwise nesting. The Wilson scale assumes locally homogeneous conditional
probabilities. Coverage is experimental; joint inference is for the finite
requested grid. The strongest evidence is on trimmed explanatory grids.}
}
\value{A validated control object for the \code{bias.control} argument of
\code{plot} applied to \code{npcdist} or \code{npcdens} fits or bandwidths.}
\details{
These are experimental opt-ins for local constant, one continuous explanatory
variable and one continuous response, fixed unbounded second-order Gaussian
kernels. The explicit pairs-score construction also admits an unbounded
second-order Epanechnikov or uniform response kernel, and a beta2 response on known finite
support with Gaussian explanatory smoothing. For beta responses the donor range is
enclosed using monotonicity of the integrated kernels at both bandwidth scales.
Gradients, quantiles, shape projections, nonfixed bandwidths and mixed
variables are not admitted. Original fits and selected bandwidths remain fixed.

Specify \code{center="bias-corrected"} and an explicit \code{band}:
\code{"pointwise"}, \code{"simultaneous"}, \code{"bonferroni"} or \code{"all"}.
For the analytic or studentized bootstrap, explicitly specify
\code{bootstrap="conditional-pilot"} and
\code{interval="symmetric-sd"} or
\code{interval="asymmetric-sd"}.
This law keeps explanatory observations fixed and draws each response from its
conditional Gaussian donor mixture, with pilot bandwidths multiplied by
\eqn{\sqrt{2}}. It is not ordinary pairs resampling and does not replace
\code{inid}, \code{fixed}, \code{geom}, residual or wild resampling.

With \code{construction="pairs-score"}, specify \code{errors="bootstrap"},
\code{bootstrap="inid"} and \code{interval="percentile"}.
For a uniform response, both ratios use its integrated rectangular kernel.
The signed donor range is computed at the four support knots. The widened
response bandwidth must also preserve the fitted kernel. Sparse explanatory
samples can produce wide simultaneous intervals even on smooth distributions.
For an Epanechnikov response, both ratios and the signed donor range use the
integrated Epanechnikov kernel at the fitted and widened response bandwidths.
On 32 saved fits with 9,999 draws and 50-by-50 targets this construction covered
30 trimmed and all 32 untrimmed surfaces. Mean simultaneous widths were 0.53
and 1.87 CDF units respectively. Pointwise intervals were conservative on both
grids. Evidence supports trimmed grids at sample sizes 150 and 400; misses in the five-by-five
comparison concentrated in the mixture design. On 64 separate dense grids, trimmed
simultaneous coverage was 62 of 64; one original of each design missed at n=400. Untrimmed joint widths commonly exceed
one CDF unit. Pointwise coverage is near 99 percent, and the dense-grid Bonferroni tails
are resolution-limited. The same seed supplies common bootstrap counts across fits of each size.
The ordinary estimate uses the native kernel constants; the corrected ratios use the
exact integrated polynomial (about 1.1e-10 maximum primitive difference).
The original kernel, bandwidth and finite target family are held fixed.

This is a fixed-X sharp-score envelope around the jackknife center with a
basic inid pairs companion. Either term can determine the displayed limits.
Each multinomial sample recomputes both conditional CDF ratios in the
two-bandwidth jackknife. Pointwise limits combine inverted signed empirical
root quantiles and a marginal sharp-score interval. The simultaneous radius
covers the empirical absolute maximum quantile and all marginal intervals.
Bonferroni limits contain both the tail-adjusted basic limits and the joint
band. The fitted bandwidths remain fixed while predictor and response observations
are resampled together. The marginal score companion conditions on the
observed predictor weights. Its certified range concerns the smoothed
functional and does not bound residual bias for the population CDF.

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
\code{interval="asymmetric-sd"} in the plot call.
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
marginal hull. The pairs-score construction uses one absolute radius across
all targets. The analytic pilot bootstrap and analytic CDF asymptotic bands
propagate required widths within explanatory slices. Companion hulls protect interval nesting. These are
approximate confidence constructions, not finite-sample honesty guarantees.
No clipping, positivity projection or forced properness is applied to corrected
centers or endpoints.

The family consists of the actual requested finite evaluation points within
one panel or surface. Separate panels are not a joint family; a finite grid is
not the entire continuous domain. The analytic and conditional-pilot methods admit up to 256 target
points and 25 distinct response coordinates. The explicit \code{"pairs-score"}
method admits up to 2,500 points and 50 response coordinates when its numerical
workspace fits. For requests up to 9,999 draws, target tiles can share the same pairs draws
when storing the full draw matrix would exceed this budget. The complete grid
still determines each simultaneous maximum. At most twelve target tiles reuse the counts. On a 50-by-50 grid the
native-center workspace admits at most 556 observations on this tiled route.
For 2,500 targets and 9,999 draws, nominal 95\\% Bonferroni limits
have only 0.10 expected draws per tail. For example, \code{neval=9} requests an 81-point surface.
On 32 saved cross-validated samples with 25 targets, nominal 95\% simultaneous
coverage was 31/32 for \code{"pairs-score"} and 0/32 for the ordinary corrected
inid percentile construction selected by omitting \code{bias.control}. On the corresponding
50-by-50 grids, the earlier untiled calculation used 797 draws at n=150 and
283 at n=400. The simultaneous band covered all 32 untrimmed surfaces, with
a width exceeding one CDF unit in every sample (mean 1.81). Trimming the
predictor by 10\% at each end gave 30/32 and mean width 0.52. The empirical
Bonferroni tails on these grids reach the extreme observed draws; their
envelope also includes the simultaneous band.
Inputs exceeding the numerical workspace budget of 256 MiB fail before
resampling; the grid and draw count are never silently reduced. The budget is
not a bound on total R-process memory. Bootstrap \code{B} and asymptotic
\code{gaussian.draws} must be integers at least two; defaults are 9999.
Low finite-draw tail resolution is reported. Alpha is propagated throughout
in \eqn{(0,0.5)}; population coverage evidence is currently for nominal 95\% only.

Data-returning plots use \code{\link{np_plot_data}}. The value tables contain
original and corrected estimates and available SEs; interval tables contain
absolute endpoints for every requested band. Compact inference records identify
the construction, law, calibration and available studentization diagnostics.
Corrected centers may be signed and need not be interval midpoints. Internal
pilot samples, covariance workspaces and RNG states are not returned.

For the earlier analytic and jackknife constructions, saved Gaussian-design studies support full-grid CDF coverage at n=100 and 400,
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
npRmpi.init(nslaves=1)
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
          interval="asymmetric-sd",
          output="plot-data")

model.d <- npcdens(y ~ x, regtype="lc")
p.d <- plot(model.d, errors="bootstrap", bootstrap="conditional-pilot",
            center="bias-corrected", band="all", neval=9,
            bias.control=np_conditional_bias_control("jackknife-mp1"),
            interval="asymmetric-sd",
            output="plot-data")
npRmpi.quit()
}
}

\section{Experimental status of bias-corrected inference}{
\strong{All bias-corrected confidence intervals and bands are experimental
(pre-beta research), incomplete, and not recommended for general-purpose
inference.} This status remains in effect unless the documentation explicitly
states otherwise for a particular method and scope. It applies wherever
bias-corrected bounds are available: asymptotic and bootstrap constructions,
including legacy methods, pointwise intervals, Bonferroni intervals and
simultaneous bands, for levels and supported derivatives.

These methods form an ongoing research project. Implementation support,
numerical checks and favorable simulation results in selected settings do not
establish generally valid coverage. Coverage can fail outside the studied
settings, and available options do not constitute a general solution to
bias-corrected inference. Method-specific assumptions, restrictions and known
failures must be considered when experimenting with these bounds.
}

\seealso{\code{\link{plot.np}}, \code{\link{np_boot_control}},
\code{\link{npcdist}}, \code{\link{npcdens}}}
\keyword{nonparametric}
