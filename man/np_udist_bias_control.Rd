\name{np_udist_bias_control}
\alias{np_udist_bias_control}
\title{Private Controls for Bias-Corrected Unconditional Distribution Bands}
\description{Selects the research CDF construction or legacy bootstrap centering.}
\usage{np_udist_bias_control(construction = c("cdf", "legacy"))}
\arguments{\item{construction}{The joint corrected CDF process or legacy correction.}}
\details{
This private research implementation is undergoing full empirical-range
qualification. It does not establish a finite-sample coverage guarantee.
With \code{center="bias-corrected"}, omitting \code{bias.control} selects
\code{"cdf"}. The default estimate center remains unchanged. Use
\code{errors="asymptotic"} for Gaussian calibration or \code{errors="bootstrap"}
with \code{bootstrap="inid"} for literal independent resampling. Explicit
\code{construction="legacy"} retains previous bootstrap correction.

The new construction currently admits fixed unbounded second-order Gaussian
continuous kernels and normalized Li--Racine ordered smoothing, preserving
integer lattice distances and declared levels. Unordered factors are excluded.
Selected bandwidths and ordered smoothing parameters stay fixed within a fit.
Joint continuous bias correction and the complete corrected observation
covariance are used. Small Gaussian families use bounded covariance calibration;
larger families use target tiles with common observation innovations. No
training-size-squared or unrestricted target covariance matrix is retained.

Independent bootstrap draws are literal multinomial resamples. Each replicate
uses its own complete corrected scale. Certified zero replicate scale with
nonzero displacement gives an infinite pivot; undefined zero over zero fails
the intact family. Draws are never removed or replaced and scales are not floored.
Bootstrap quantiles use type one. Bonferroni empirical tail resolution is limited
by the requested draw count; the raw and nesting-adjusted critical values are
retained. Gaussian pointwise and Bonferroni quantiles are analytic. Bands coincide
for singleton families and nest pointwise, simultaneous, Bonferroni.

Use an explicit \code{band} of \code{"pointwise"}, \code{"simultaneous"},
\code{"bonferroni"} or \code{"all"}. Nearest-neighbor bandwidths, fixed/geom
bootstrap laws, bounded or higher-order kernels, large-bandwidth limits and
factor boxplots remain unqualified for this private corrected route and fail
clearly. Their legacy routes remain available. No explicit law is substituted.
}
\value{A versioned control object for the plot \code{bias.control} argument.
For data-returning corrected plots, each panel's \code{inference$bands} contains
absolute endpoints; \code{bias.corrected + derr[,1]} and
\code{bias.corrected + derr[,2]} give matching signed-offset endpoints.
For \code{band="all"}, \code{derr} refers to the pointwise band.}
\seealso{\code{\link{npudist}}, \code{\link[=plot.np]{npplot}}}
