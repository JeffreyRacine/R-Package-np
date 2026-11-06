\name{np_udens_bias_control}
\alias{np_udens_bias_control}
\title{Controls for Bias-Corrected Unconditional Density Bands}
\description{
\strong{Bias-corrected bounds are experimental (pre-beta research), incomplete,
and not recommended for general-purpose inference.} See the experimental-status
section below.
Selects joint-scale bias-corrected density intervals, the earlier positive-covariance construction, or legacy bootstrap centering.}
\usage{np_udens_bias_control(construction = c("joint-scale", "positive-covariance", "legacy"))}
\arguments{\item{construction}{The joint-scale envelope, earlier positive-covariance construction, or existing smooth-bootstrap centering route.}}
\details{
With \code{center="bias-corrected"}, \code{errors="asymptotic"}, and
omitted or \code{NULL} \code{bias_control}, fixed unbounded second-order
Gaussian density plots with one to three continuous numeric variables use the
D8 positive-pilot construction described below, for \eqn{10^{-4}<\alpha<0.5}.
Other previously supported requests retain their existing construction and
admission rules. Supplying \code{np_udens_bias_control()} explicitly selects
the joint-scale construction, including for asymptotic inference. Omitted
controls for unbounded Gaussian bootstrap inference continue to select joint-scale.
Known-bound beta plots use the separate construction described below.
The default \code{center="estimate"} is unchanged.
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

For the joint-scale construction, omitting \code{band}, or selecting
\code{band="pmzsd"}, uses the existing
normal-quantile convention: the corrected center plus or minus
\eqn{z_{1-\alpha/2}} times the returned \code{std_error}. The family provider supplies
the corrected scale; its pointwise calibration work and RNG contract are
retained. These normal intervals do not apply the joint-scale uncertainty
envelope and do not inherit its empirical coverage results. Their endpoints
are in the interval rows with \code{band="pmzsd"}. Explicit \code{"pointwise"}, \code{"simultaneous"},
\code{"bonferroni"} and \code{"all"} retain their calibrated constructions.


\subsection{Known-bound beta density inference}{
With omitted or \code{NULL} \code{bias_control}, fixed second-order beta
kernels support \code{errors="asymptotic"} and independent bootstrap
\code{calibration="symmetric-sd"} or \code{"asymmetric-sd"}. Request
\code{center="bias-corrected"} and explicitly select \code{band="pointwise"}
or \code{"simultaneous"}. This plot-only construction supports one or two
continuous variables, at least ten observations and known finite coordinatewise
\code{ckerbound="fixed"} bounds. Observations and evaluation points must lie
strictly inside those bounds. Estimated \code{ckerbound="range"} bounds and
literal support endpoints are rejected; no jittering or target deletion occurs.
Bandwidths, bounds and targets stay fixed during inference. Core fitted values,
estimator standard errors, bandwidth selection and \code{se=FALSE} are unchanged.
An object fitted with \code{se=FALSE} can receive explicitly requested plot bands.

On normalized support, write \eqn{L_r=4K_r-6K_{2r}+4K_{3r}-K_{4r}} and
\eqn{\Psi=L_h+L_b-L_hL_b}. Products across coordinates are formed within each
scale term. The corrected plot center is
\eqn{C=\widehat R+g_{\widehat\theta}}, where
\eqn{\widehat R=P_n\Psi-E_{g_{\widehat\theta}}\Psi}. The fitted guide has two
edge-shape parameters in one dimension, and four edge-shape parameters plus
an interaction in two dimensions. Half the nominal error probability is
assigned to a nonlinear guide confidence region and half to the residual
process. The residual covariance includes the guide-estimation projection
and cross terms. A positive ratio pilot supplies its transported scale.
This scale is not a complete normal standard error for the composite interval.

Asymptotic pointwise intervals use \eqn{z_{1-\alpha/4}}; simultaneous bands
calibrate the absolute Gaussian maximum under the complete residual covariance
on the finite requested grid. For both SD variants, every literal independent
multinomial resample refits the guide in the residual; the original empirical
residual scale is held fixed. This is not replicate studentization or a
bootstrap of the complete corrected center followed by a second guide region.
\code{"asymmetric-sd"} inverts the signed residual roots. \code{"symmetric-sd"}
uses their absolute roots, then takes the smallest outward-rounded interval
symmetric about the corrected center that contains the composite interval.
The guide region is added once in each case. These beta-specific rules do not
change the existing Gaussian density bootstrap constructions.

All slices in one call share one calibration family. Returned endpoints are
absolute limits; they and the corrected center remain signed, without clipping
or renormalization. Stored plotting reuses those limits. The selected RNG kind
and literal resampling law are retained, including Box--Muller continuation.
Low tail-draw resolution warns and returns finite limits; failed numerical
certificates and workspace limits are errors. The requested \code{B} or
\code{gaussian_draws} is never reduced. The current implementation uses bounded
resident workspaces and invocation-owned temporary files, with limits of
256 MiB aggregate inference workspace and 512 MiB temporary storage. These
limits are distinct from total R-process memory.

The retained Beta-family, uniform, mixture and dependent two-dimensional
experiments show generally conservative finite-grid coverage. They do not
establish continuum coverage, honesty over unknown boundary smoothness,
uncertainty in estimated bounds or bandwidth-selection uncertainty. Very wide
near-edge intervals, especially for singular densities, remain possible.
The experimental status of bias-corrected inference applies to this method too.
Other bounded kernels, nearest-neighbour bandwidths, mixed data, dimensions
above two, derivatives, other bootstrap laws and the \code{pmzsd},
\code{bonferroni} and \code{all} bands are outside this new beta construction.
Explicit saved density controls keep their original construction and scope.
}

\subsection{Default D8 asymptotic density inference}{
Let \eqn{H=\mathrm{diag}(h_1^2,\ldots,h_q^2)} and
\eqn{D_H=\sum_j h_j^2\partial_{jj}}. The plotted corrected curve averages
\eqn{L_8=(I-D_H/2+D_H^2/8-D_H^3/48)G_H} over the observations.
This total-degree radial correction includes multivariate cross terms; it is
not a product of univariate eighth-order kernels. The fitted density,
bandwidths and estimator standard errors are unchanged.

The working covariance is \eqn{\mathrm{Cov}_P(L_8(g-X),L_8(g'-X))/n},
where the positive pilot \eqn{P} is the Gaussian mixture centered at the
observations with the same fitted bandwidths. Both pilot means are subtracted.
Pointwise and \code{pmzsd} intervals use normal critical values and these
complete corrected scales. Simultaneous bands calibrate the standardized
absolute Gaussian maximum on the complete requested grid. All slices from
one plot call form one union; a surface uses its full Cartesian grid.
Bonferroni bands use normal tails and an outward nesting adjustment.
Corrected curves and endpoints remain signed, without clipping or renormalization.

Bounded resident panels and invocation-owned temporary factor files avoid an
unrestricted resident grid-square matrix. A row-major draw workspace reduces
repeated file seeks; peak temporary disk can be twice the completed factor
size. The resident allocation forecast is distinct from total R heap/RSS.
In \pkg{npRmpi}, this conservative inference-storage forecast is summed
across the coordinator and workers. Increasing the worker count can therefore
exceed the 256 MiB forecast even for a grid that a smaller pool admits.
The error does not imply that the observed R resident memory reached that
limit. An explicitly selected smaller pool can retain the same complete grid;
the package does not reduce the grid or switch methods automatically.
A residual covariance bound is included in maximum calibration. Requested
Gaussian draw counts are retained, with the existing precision warning when
small. Seeded calls preserve ordinary caller RNG continuation, including
Box--Muller via an isolated base-R normal stream; no RNG kind is remapped.

Coverage evidence is strongest for one and two continuous variables. In the
retained ordinary cases, nominal 95 percent simultaneous coverage ranged from
93.4 to 97.2 percent for univariate cases (500 replications per cell) and from
92.8 to 98.0 percent for bivariate cases (250 per cell). Narrow-spike and
sparse-tail limitations remain. Three-dimensional support has arithmetic and
engineering checks, not comparable broad coverage qualification. These are
finite-grid, fixed-bandwidth approximations, not uniform honesty guarantees or
uncertainty for bandwidth selection. Bootstrap and explicit saved controls
retain the distinct constructions documented below.
}

The joint-scale construction supports fixed, unbounded, second-order Gaussian
continuous kernels, including multivariate and mixed-data layouts with unordered
Aitchison--Aitken and ordered Li--Racine kernels. Mathematical and implementation
support is distinct from statistical coverage qualification: one or two continuous
variables remain the primary assessment, and mixed/higher-dimensional evidence
is limited. Saved controls retain their construction semantics.
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

Defined version-2/3 bootstrap bands are returned at every positive draw count,
including Bonferroni bands on the full surface. Coarse location or fixed-beta
scale quantiles produce one aggregate warning and retained
\code{inference[[1]]$resolution} diagnostics. Small draw counts do not establish
statistical qualification. Invalid inputs and undefined arithmetic remain errors.
No law or draw count is silently changed. Raw empirical tail quantiles do not
establish a finite-sample Bonferroni guarantee; raw and nesting-adjusted
criticals are retained separately. Gaussian Bonferroni uses analytic normal
tails. Singleton bands coincide.

Bootstrap \code{inid}, \code{fixed} and \code{geom} share the same corrected
provider. Explicit block length one reproduces independent resampling exactly;
longer fixed/geometric blocks remain different laws. Mixed-data block requests
require an explicit \code{np_boot_control(blocklen=...)} while automatic mixed
block-length selection remains unqualified. Dependent-data coverage has not been
established by the independent-data evidence below.

Select bootstrap calibration through \code{\link{np_boot_control}}:
\code{"symmetric-sd"} retains absolute-root bands;
\code{"asymmetric-sd"} retains transported asymmetric tail bands;
\code{"smoothed-percentile"} uses outer-calibrated percentiles of complete
unshifted corrected smoothed-bootstrap curves, with companion-band nesting and
certified full-support completion. The last method is the refined OP1+SC1
construction with multivariate/mixed-data and categorical-kernel extensions;
it replaces the earlier private centered-smoothed experiment as the public
percentile candidate. Its empirical SD is descriptive and does not scale
endpoints. The original fit and bandwidth are unchanged; the 39 fresh outer
samples receive their own CV bandwidth selection. Returned settings distinguish
the construction and retain pilot, support, resolution and stream provenance.
This method requires \pkg{ks}, at least one numeric coordinate, and inid or
unit-length fixed/geometric blocks. Nonunit-block calibration and NN are not
implemented. Its conditional bootstrap curves can be signed. Public surface
and slice family definitions are preserved, and no grid is trimmed implicitly.
Full-support intervals may be wide. A still-unattainable calibration remains
explicitly diagnostic-only. The q=1 Gaussian comparison is an outcome-informed
paired assessment; multivariate/mixed capability proof is broader than current
coverage evidence. The SD-method evidence below does not qualify percentile
coverage in those other domains.

An unused pure categorical target may have zero empirical variance, making SD
studentization undefined. Such requests retain their explicit arithmetic error;
this is different from a small-draw resolution warning. Direct percentile bands
remain available when their corrected curves are finite, including constant
curves. No variance floor or silent method substitution is used.

Nearest-neighbor bandwidths, higher-order or bounded
continuous kernels, factor boxplots and derivatives remain outside
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
may be negative and are not clipped or renormalized. Absolute endpoints are returned consistently for all
constructions, including legacy and version-1 surfaces; their internal width
conventions do not appear in the public plot result.
}
\value{A family-specific control object consumed by plot methods.
Data-returning plots use \code{\link{np_plot_data}}: original and corrected
values, available standard errors, absolute endpoints and compact method and
resolution diagnostics. Select a panel, scalar target and band with
\code{as.data.frame(result, panel=1, target=1, band="simultaneous")}.
Full fitted state, covariance workspaces, bootstrap draws and RNG states are
internal computation records, not public plot fields. The fitted object is unchanged.}

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

\seealso{\code{\link{npudens}}, \code{\link[=plot.np]{npplot}}, \code{\link{np_bias_control}}}
\examples{
npRmpi.init(nslaves = 1)
set.seed(123)
x <- data.frame(x = rnorm(80))
bw <- npudensbw(dat = x, bws = 0.45, bandwidth.compute = FALSE)
z <- plot(bw, xdat = x, errors = "bootstrap", bootstrap = "inid",
          center = "bias-corrected", band = "all", B = 399,
          neval = 3, output = "data")
z$inference[[1]]$primary
npRmpi.quit()
}
\keyword{nonparametric}
