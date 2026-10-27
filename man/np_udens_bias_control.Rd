\name{np_udens_bias_control}
\alias{np_udens_bias_control}
\title{Controls for Bias-Corrected Unconditional Density Bands}
\description{Selects an optional positive-covariance construction for fixed-bandwidth unconditional density plotting. Existing plot defaults are unchanged.}
\usage{np_udens_bias_control(construction = c("positive-covariance", "legacy"))}
\arguments{\item{construction}{The new corrected contribution and covariance construction, or the existing smooth-bootstrap centering route.}}
\details{
Pass the result as \code{bias.control} to \code{plot} or \code{npplot}, with
\code{errors="bootstrap"}, \code{bootstrap="inid"},
\code{center="bias-corrected"} and an explicit \code{band="pointwise"},
\code{"simultaneous"}, \code{"bonferroni"} or \code{"all"}.
Omitting \code{bias.control} retains the existing construction, as does an
explicit \code{construction="legacy"}. Regression controls and iterations
do not apply to this density construction.

The new construction supports fixed, unbounded, second-order Gaussian
continuous kernels with the supported Aitchison--Aitken unordered and
Li--Racine ordered kernels. Numeric, factor and ordered variables can be
combined, including nonbinary categories. It corrects continuous and categorical
smoothing contributions, preserves the selected bandwidths, and uses complete
replicate studentization under literal independent resampling. Its observed
positive covariance uses the same-bandwidth pilot and a positive marginal
congruence correction. Absolute studentized pivots calibrate all three bands;
there is no normal-critical floor, draw deletion or redraw.

The initial construction is limited to 512 evaluation points per panel and a
256 MiB resident-process memory forecast. For a two-dimensional surface,
\code{neval=20} uses 400 points; the default 50 by 50 surface exceeds the cap.
Use an explicit smaller grid or slices. No grid is silently reduced. The
\code{pmzsd} band and bootstrap factor boxplots are not available for this
construction. Existing large-grid and boxplot routes remain available with
the legacy construction.

Gaussian calibration, nearest-neighbor bandwidths, dependent resampling,
higher-order kernels, bounded continuous kernels and derivatives are not
admitted to the new public construction. Unsupported requests fail rather
than changing their law or bandwidth type. Existing plot routes for those
settings are unchanged.

Finite-grid coverage improved substantially in bounded fixed-bandwidth,
data-selected-bandwidth experiments with one or two numeric variables and
mixed data. This is not a uniform honesty claim. Simultaneous and Bonferroni
bands can be conservative and wide, especially for sparse mixed regions.
A zero replicate scale with nonzero centered numerator produces an infinite
pivot; undefined zero-over-zero pivots or uncertified zero scales fail the
request. Infinite endpoints are reported, not discarded or replaced by a
variance floor. The fitted density itself is unchanged; the corrected center
can be signed and is not clipped or renormalized.
}
\value{A family-specific control object. Returned plot-data density objects
include an \code{inference} component with the corrected center, covariance,
absolute bands, bandwidths, effective calibration settings, degeneracy counts
and random-stream provenance.}
\seealso{\code{\link{npudens}}, \code{\link{npplot}}, \code{\link{np_bias_control}}}
\examples{
set.seed(123)
x <- data.frame(x = rnorm(80))
bw <- npudensbw(dat = x, bws = 0.45, bandwidth.compute = FALSE)
z <- plot(bw, xdat = x, errors = "bootstrap", bootstrap = "inid",
          center = "bias-corrected", band = "all", B = 399,
          neval = 3, output = "data", bias.control = np_udens_bias_control())
z[[1]]$inference$settings$primary
}
\keyword{nonparametric}
