\name{np_udist_bias_control}
\alias{np_udist_bias_control}
\title{Private Controls for Bias-Corrected Unconditional Distribution Bands}
\description{Selects the research CDF construction or legacy bootstrap centering.}
\usage{np_udist_bias_control(construction = c("joint-process", "cdf", "legacy"),
                      calibration = c("default", "score-envelope"))}
\arguments{
\item{construction}{The joint corrected CDF process or legacy correction.}
\item{calibration}{The default construction or the explicit, side-specific
score envelope for supported fixed numeric Gaussian asymptotic bands.}
}
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
smoothing retains the canonical categorical representation: nonnumeric ordered
labels use ordinal ranks, while numeric labels retain their metric gaps. It
uses the normalized Li--Racine kernel and the complete joint corrected
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
the preceding spacing-pilot bootstrap-t method. If a bootstrap query radius
collapses to zero, the Gaussian integrated kernel uses its exact one-sided limit:
one below the query, one half at a tie, zero above it; its correction moment is
zero. Positive-radius coordinates retain their original arithmetic. This is an
extension of the bootstrap CDF statistic, with no draw discarded, replaced or
floored. Observed zero radii and invalid Gaussian spacing pilots remain unavailable. Pure ordered NN inputs retain the
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
Every positive draw count is supported when the provider arithmetic is defined.
Poorly resolved empirical quantiles return all requested bands with one aggregated
warning and retained \code{inference$resolution} diagnostics. No draw count is
increased and no requested band is omitted for resolution alone. Numeric GNN
bootstrap scales are empirical standard deviations across draws, so fewer than
two draws remains a distinct undefined-arithmetic error. Analytic Gaussian
marginal quantiles do not require empirical tail resolution.
The maximum-process draws are retained as \code{inference$process.maximum}.

For numeric-only fixed bandwidths with one or two variables, Gaussian and inid
pointwise bands additionally envelope a signed-contribution score interval.
The contribution has a certified range: writing \eqn{t=\sqrt{1+2/q}}, its
lower bound is \eqn{a=\Phi(-t)-(q/2)t\phi(t)} and upper bound is \eqn{1-a}.
A Wilson score interval for the rescaled contribution mean is mapped back to
this range. This is approximate inference for the fixed-bandwidth contribution
expectation, not an exact confidence interval for the unsmoothed CDF.
Data-selected bandwidths and evaluation ranges require separate empirical
qualification. The inid interval is a hybrid marginal envelope: literal draws
and quantiles are unchanged, but the final endpoints can be wider and asymmetric.
The corrected center is unchanged. Simultaneous bands retain their raw process
calibration and expand only by the minimum uniform amount required to contain
the revised pointwise bands. Bonferroni bands retain nesting; singleton bands
coincide. \code{inference$reference.bands} and
\code{inference$reference.critical} retain the preceding intervals and criticals.
Final \code{critical} values report maximum one-sided extents around the
preserved center; use \code{bands} for the asymmetric endpoints. Settings record
the contribution range and nesting expansion. Mixed, NN and block routes do
not receive this marginal envelope. An entirely zero observed process remains
unavailable.

For supported corrected empirical bootstrap bands, explicit
\code{boot.control=np_boot_control(calibration="separate-tails")} uses signed
lower and upper empirical tails in CDF units. The simultaneous roots are
separate quantiles of the family minimum and maximum; endpoints invert those
signed estimate displacements. \code{calibration="tail-envelope"} additionally
takes the interval hull with the unchanged absolute construction. The default
is unchanged, and \code{calibration="absolute"} requests it explicitly.
These options preserve the selected inid, fixed or geom law, the complete
contribution process, block choice, corrected center and displayed variance.
They currently admit fixed and purely ordered NN geometry; numeric NN
asymmetric calibration requires separate qualification. Density's
\code{"percentile-tails"} choice is not a CDF method.

For fixed numeric Gaussian2 inputs with one or two variables, the inid signed
bands retain the existing marginal score protection. Its two sides are nested
outward separately, preserving asymmetry instead of imposing the larger
uniform displacement on both sides. An explicit asymptotic option,
\code{bias.control=np_udist_bias_control(calibration="score-envelope")},
uses this same side-specific composition with the unchanged Gaussian process.
Gaussian draws remain reflection-symmetric; this option's possible asymmetry
comes from the deterministic score envelope, and some bands can remain
symmetric. Default pointwise numeric bands already include score protection.
The new Gaussian option changes how that protection extends to the simultaneous
and Bonferroni bands; it does not simulate a skewed Gaussian law.

Explicit asymmetric options require \code{band="pointwise"},
\code{"simultaneous"}, \code{"bonferroni"} or \code{"all"}; they do not
silently alter normal-scale \code{"pmzsd"} intervals. Singleton bands coincide,
and the declared family remains the full returned grid within each panel.
Results retain \code{absolute.intervals}, \code{raw.process.intervals},
\code{raw.intervals}, signed tail orientation where applicable and side-specific
nesting diagnostics. No [0,1] clipping or monotonicity projection is applied.
These are private candidates requiring their own coverage evidence; containing
a baseline interval is not a proof of nominal coverage.

Explicit \code{construction="cdf"} retains the preceding private studentized
construction, including normalized Li--Racine ordered smoothing and bounded
small-grid covariance calibration. Its bootstrap uses replicate scales and
retains the extended degenerate-pivot contract. This explicit construction does
not inherit qualification from the absolute process. Explicit \code{"legacy"}
retains the historical bootstrap center. All constructions remain research
candidates; full empirical-range screening is distinct from confirmation.

Omitting \code{band}, or selecting \code{band="pmzsd"}, uses the corrected
center plus or minus \eqn{z_{1-\alpha/2}} times \code{inference$se}.
The existing pointwise provider supplies the selected-law scale, retaining
block dependence and exact GNN resampling where applicable. Its computation
and RNG contract are retained. Zero marginal scales give zero normal width.
These normal intervals do not apply the score or joint envelopes and do not
inherit their coverage results. Endpoints are in \code{inference$bands$pmzsd};
\code{inference$normal.band} records the rule. Explicit \code{"pointwise"},
\code{"simultaneous"}, \code{"bonferroni"} and \code{"all"} retain their
calibrated constructions. Numeric adaptive nearest-neighbor bandwidths, unordered factors,
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
