\name{np_udens_bias_control}
\alias{np_udens_bias_control}
\title{Controls for Bias-Corrected Unconditional Density Bands}
\description{Selects joint-scale bias-corrected density intervals, the earlier positive-covariance construction, or legacy bootstrap centering.}
\usage{np_udens_bias_control(construction = c("joint-scale", "positive-covariance", "legacy"))}
\arguments{\item{construction}{The joint-scale envelope, earlier positive-covariance construction, or existing smooth-bootstrap centering route.}}
\details{
With \code{center="bias-corrected"}, \code{errors="asymptotic"}, and
omitted or \code{NULL} \code{bias.control}, fixed unbounded second-order
Gaussian density plots with one to three continuous numeric variables use the
D8 positive-pilot construction described below, for \eqn{10^{-4}<\alpha<0.5}.
Other previously supported requests retain their existing construction and
admission rules. Supplying \code{np_udens_bias_control()} explicitly selects
the joint-scale construction, including for asymptotic inference. Omitted
controls for bootstrap inference continue to select joint-scale.
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

Omitting \code{band}, or selecting \code{band="pmzsd"}, uses the existing
normal-quantile convention: the corrected center plus or minus
\eqn{z_{1-\alpha/2}} times \code{inference$se}. The family provider supplies
the corrected scale; its pointwise calibration work and RNG contract are
retained. These normal intervals do not apply the joint-scale uncertainty
envelope and do not inherit its empirical coverage results. Their endpoints
are in \code{inference$bands$pmzsd}; \code{inference$normal.band} records the
rule. Explicit \code{"pointwise"}, \code{"simultaneous"},
\code{"bonferroni"} and \code{"all"} retain their calibrated constructions.


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
\code{inference$resolution} diagnostics. Small draw counts do not establish
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
may be negative and are not clipped or renormalized. For version-2 and version-3 outputs,
\code{bias.corrected + derr[,1]} and \code{bias.corrected + derr[,2]} recover
lower and upper endpoints for curves and surfaces. Legacy/version-1 surfaces
retain their historical positive lower-width convention.
}
\value{A family-specific control object. Corrected plot-data objects include an
\code{inference} component with the center, marginal scales where applicable, bands,
complete covariance or its factor-representation description, raw and adjusted
criticals, effective calibration settings and random-stream provenance.
Smoothed percentiles retain the selected pilot, outer calibration records,
rank/content/support diagnostics and separate stream provenance; they do not retain the full draw matrix.}
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
