\name{np_udist_bias_control}
\alias{np_udist_bias_control}
\title{Private Controls for Bias-Corrected Unconditional Distribution Bands}
\description{Selects the research CDF construction or legacy bootstrap centering.}
\usage{np_udist_bias_control(construction = c("joint-process", "cdf", "legacy"))}
\arguments{\item{construction}{The joint corrected CDF process or legacy correction.}}
\details{
This private research implementation is undergoing full empirical-range
qualification. It does not establish a finite-sample coverage guarantee.
With \code{center="bias-corrected"}, omitting \code{bias.control} selects
\code{"joint-process"}. The default estimate center remains unchanged. Use
\code{errors="asymptotic"} for Gaussian calibration or \code{errors="bootstrap"}
with \code{bootstrap="inid"}, \code{"fixed"} or \code{"geom"} for literal
independent, circular fixed-length or stationary geometric-block resampling. Explicit
\code{construction="legacy"} retains previous bootstrap correction.

The default joint-process construction admits fixed and generalized-nearest-neighbor
unbounded second-order Gaussian kernels with at most two numeric and two ordered variables. Ordered
smoothing uses the normalized Li--Racine kernel and the complete joint corrected
contribution, retaining declared levels and integer gaps. It calibrates the maximum
absolute corrected CDF displacement, with complete covariance and common
observation innovations across target tiles. For fixed bandwidths, Gaussian and literal multinomial
inid processes use the same centered observation factor. Critical values are in
CDF probability units, not standardized units. Selected bandwidths or neighbor counts and categorical smoothing parameters stay fixed.
With fixed bandwidths, fixed/geom processes use literal selected-law block counts and the corresponding
finite-sample covariance; their reported scales are not IID scales. The automatic
block length remains \code{b.star(xdat, round=TRUE)[1,1]}, selected once per plot.
The first input series may not represent dependence elsewhere in multivariate
data. Explicit \code{boot.control=np_boot_control(blocklen=...)} overrides are
honored. Neither method divides by a replicate studentizer. Full-length fixed
blocks with zero process variance are unavailable.
No local or replicate studentizer division, variance floor, trimming, endpoint
clipping or resampling-law substitution is used. An entirely zero observed
process is unavailable; zero marginal scales are reported explicitly.

For generalized nearest-neighbor bandwidths, Gaussian calibration retains the
complete influence of the sample-specific query radii. Bootstrap calibration
recomputes every physical radius and corrected center using the literal selected
law, and uses their absolute displacements. The displayed corrected center retains
its original summation; paired evaluations of the same algebraic sums ensure
zero displacement for an unchanged sample. Reported bootstrap scales are empirical
standard deviations under that selected law, including radius randomness.
This construction has no replicate studentizing divisor and does not qualify
the preceding spacing-pilot bootstrap-t method. A zero physical radius is
unavailable, including when produced by ordinary CV and literal resampling;
no draw is discarded, replaced or floored. Pure ordered NN inputs retain the
fixed-geometry identity. Numeric adaptive NN remains parked at its population
bias-remainder gate. These extensions require their own public-call coverage
screen and do not inherit coverage from fixed bandwidths.

Gaussian marginal quantiles are analytic. Joint critical values use type-one
quantiles bounded by the analytic marginal and Bonferroni maximum bounds.
Inid marginal and joint quantiles use type one. Bonferroni intervals envelope
the joint band. Singleton bands coincide. Standalone Gaussian Bonferroni uses
the default joint draws for its envelope; request \code{band="all"} to set
\code{gaussian.draws} and obtain all bands from the same process. Finite-draw
Bonferroni tail resolution and marginal sparse-tail coverage remain limitations.
The maximum-process draws are retained as \code{inference$process.maximum}.

Explicit \code{construction="cdf"} retains the preceding private studentized
construction, including normalized Li--Racine ordered smoothing and bounded
small-grid covariance calibration. Its bootstrap uses replicate scales and
retains the extended degenerate-pivot contract. This explicit construction does
not inherit qualification from the absolute process. Explicit \code{"legacy"}
retains the historical bootstrap center. All constructions remain research
candidates; full empirical-range screening is distinct from confirmation.

Use an explicit \code{band} of \code{"pointwise"}, \code{"simultaneous"},
\code{"bonferroni"} or \code{"all"}. Numeric adaptive nearest-neighbor bandwidths, unordered factors,
bounded or higher-order kernels, large-bandwidth limits and
factor boxplots remain unqualified for this private corrected route and fail
clearly. Their legacy routes remain available. No explicit law is substituted.

The finite family is the returned grid within each panel. Above two input
variables, ordinary plots return slices; their bands do not establish simultaneous
coverage over the union of panels or a full multidimensional rectangle. Ordered
and block extensions have separate qualification gates and do not inherit
coverage evidence from the numeric IID construction.
}
\value{A versioned control object for the plot \code{bias.control} argument.
For data-returning corrected plots, each panel's \code{inference$bands} contains
absolute endpoints; \code{bias.corrected + derr[,1]} and
\code{bias.corrected + derr[,2]} give matching signed-offset endpoints.
For \code{band="all"}, \code{derr} refers to the pointwise band.}
\seealso{\code{\link{npudist}}, \code{\link[=plot.np]{npplot}}}
