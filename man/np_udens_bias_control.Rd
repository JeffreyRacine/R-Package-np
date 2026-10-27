\name{np_udens_bias_control}
\alias{np_udens_bias_control}
\title{Controls for Bias-Corrected Unconditional Density Bands}
\description{Selects the fixed-bandwidth corrected density construction or explicit legacy bootstrap centering. This detached research candidate has not been qualified for public activation.}
\usage{np_udens_bias_control(construction = c("positive-covariance", "legacy"))}
\arguments{\item{construction}{The corrected contribution and covariance construction, or the existing smooth-bootstrap centering route.}}
\details{
With \code{center="bias-corrected"}, omitting \code{bias.control} selects the
new construction. The default \code{center="estimate"} is unchanged.
Use \code{errors="asymptotic"} for Gaussian calibration or
\code{errors="bootstrap"}, \code{bootstrap="inid"} for literal independent
resampling. Explicit \code{construction="legacy"} retains the previous
bootstrap correction. Previously saved version-1 controls retain their F4
bootstrap construction and its limitations. Newly created controls have version 2.
Regression iterations and residual/pilot assumptions do not apply.

The candidate supports fixed, unbounded, second-order Gaussian continuous
kernels with the supported Aitchison--Aitken unordered and Li--Racine ordered
kernels, including nonbinary categories. Selected bandwidths and categorical
smoothing parameters remain fixed. The corrected center and complete
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
this candidate. Unsupported corrected requests fail clearly. Existing estimate
and explicit legacy routes retain their previous contracts and limitations.

The full-sample-range empirical coverage of this recovery candidate remains
under investigation. Earlier interior-grid evidence does not qualify it.
The fitted density is unchanged; the corrected center and confidence endpoints
may be negative and are not clipped or renormalized. For version-2 outputs,
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
set.seed(123)
x <- data.frame(x = rnorm(80))
bw <- npudensbw(dat = x, bws = 0.45, bandwidth.compute = FALSE)
z <- plot(bw, xdat = x, errors = "bootstrap", bootstrap = "inid",
          center = "bias-corrected", band = "all", B = 399,
          neval = 3, output = "data")
z[[1]]$inference$settings$primary
}
\keyword{nonparametric}
