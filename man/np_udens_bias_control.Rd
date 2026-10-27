\name{np_udens_bias_control}
\alias{np_udens_bias_control}
\title{Controls for Bias-Corrected Unconditional Density Bands}
\description{Selects joint-scale bias-corrected density intervals, the earlier positive-covariance construction, or legacy bootstrap centering.}
\usage{np_udens_bias_control(construction = c("joint-scale", "positive-covariance", "legacy"))}
\arguments{\item{construction}{The joint-scale envelope, earlier positive-covariance construction, or existing smooth-bootstrap centering route.}}
\details{
With \code{center="bias-corrected"}, omitting \code{bias.control} selects the
joint-scale construction. The default \code{center="estimate"} is unchanged.
Use \code{errors="asymptotic"} for Gaussian calibration or
\code{errors="bootstrap"}, \code{bootstrap="inid"} for literal independent
resampling. Explicit \code{construction="legacy"} retains the previous
bootstrap correction. Explicit \code{construction="positive-covariance"}
creates a version-1 control retaining the earlier F4 bootstrap construction
and its mixed-data scope and limitations. Saved positive-covariance controls
with versions 1, 2 and 3 retain their original constructions; version 2 is the
fixed reference and version 3 is the joint-scale envelope.
New \code{construction="joint-scale"} controls have version 3 and add a joint log-scale
uncertainty envelope with fixed scale-tail probability 0.05. The F3 center
and F4 marginal scales remain unchanged; the complete influence of both
pilot variances supplies the additional scale uncertainty. This is a
practical envelope with bounded empirical qualification, not a finite-sample
coverage guarantee.
Regression iterations and residual/pilot assumptions do not apply.

The joint-scale construction supports fixed, unbounded, second-order Gaussian continuous
kernels for one or two numeric variables. Mixed-data qualification of the
new envelope remains deferred; saved controls retain their previous scope.
Selected bandwidths remain fixed. The corrected center and complete
observation contributions are unchanged from the previous density correction.
The covariance combines the previous positive-pilot marginal scales with the
complete empirical corrected correlation. It intentionally differs from the
previous pilot correlation.

The independent bootstrap draws literal multinomial counts. Its calibrated
statistic is the centered corrected observation process transported to the
observed marginal scales, using a fixed observed normalization. It is not the
previous replicate-studentized bootstrap and is not the untransformed
bootstrap distribution of the fitted estimator. Zero or nonfinite observed
scales fail explicitly; there are no variance floors, deleted draws or redraws.

Gaussian grids of at most 512 targets use the existing covariance calibrator.
Larger grids use bounded observation-factor tiles with common innovations
across the entire grid. The default two-dimensional 50 by 50 surface retains
all targets and corners. No training-size-squared matrix or full grid-by-draw
matrix is retained. The inference-storage forecast remains limited to 256 MiB
summed over participating processes, excluding baseline R/package memory.

Bootstrap Bonferroni requests fail when \code{alpha * B / m < 1}, where
\code{m} is the panel grid size. In particular the default 2,500-target grid
and 9,999 draws do not admit the new bootstrap Bonferroni or \code{band="all"}
route. Pointwise and simultaneous requests remain distinct admitted routes.
No law or draw count is silently changed. Raw empirical tail quantiles do not
establish a finite-sample Bonferroni guarantee; raw and nesting-adjusted
criticals are retained separately. Gaussian Bonferroni uses analytic normal
tails. Singleton bands coincide.

Nearest-neighbor bandwidths, dependent bootstrap laws, higher-order or bounded
continuous kernels, \code{pmzsd}, factor boxplots and derivatives remain outside
the joint-scale construction. Unsupported corrected requests fail clearly. Existing estimate
and explicit legacy routes retain their previous contracts and limitations.

Full-sample-range tests use ordinary per-dataset bandwidth selection and the
public plot calls, with 50 points per axis. At 95 percent nominal coverage,
normal and smooth-mixture univariate cases had simultaneous coverage from
96.9 to 99.2 percent; a correlated bivariate normal case with 800 observations
had 97.7 percent for both calibrations. The continuing cells used 128 datasets.
Smaller bivariate normal and mixture cases stopped at 64 datasets because
bands exceeded the prespecified width criteria. Local pointwise coverage can
remain weak: the bivariate mixture minimum was 81.3 percent despite favorable
simultaneous coverage. These finite-grid results do not establish arbitrary
sample-size, distribution or continuum coverage. At 90 percent nominal, the
univariate normal 100-observation Gaussian case covered 97.7 percent on 128
retained datasets, so substantial conservatism remains possible.
Full bivariate calibration can take materially longer than a univariate plot.
Earlier interior-grid evidence is not used to qualify sample-range plots.
The fitted density is unchanged; the corrected center and confidence endpoints
may be negative and are not clipped or renormalized. For version-2 and version-3 outputs,
\code{bias.corrected + derr[,1]} and \code{bias.corrected + derr[,2]} recover
lower and upper endpoints for curves and surfaces. Legacy/version-1 surfaces
retain their historical positive lower-width convention.
}
\value{A family-specific control object. Corrected plot-data objects include an
\code{inference} component with the center, marginal scales, absolute bands,
complete covariance or its factor-representation description, raw and adjusted
criticals, effective calibration settings and random-stream provenance.}
\seealso{\code{\link{npudens}}, \code{\link[=plot.np]{npplot}}, \code{\link{np_bias_control}}}
\examples{
npRmpi.init(nslaves = 1)
set.seed(123)
x <- data.frame(x = rnorm(80))
bw <- npudensbw(dat = x, bws = 0.45, bandwidth.compute = FALSE)
z <- plot(bw, xdat = x, errors = "bootstrap", bootstrap = "inid",
          center = "bias-corrected", band = "all", B = 399,
          neval = 3, output = "data")
z[[1]]$inference$settings$primary
npRmpi.quit()
}
\keyword{nonparametric}
