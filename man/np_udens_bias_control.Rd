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
Known-bound beta and Gaussian plots use the separate constructions described below.
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
or \code{"simultaneous"}, \code{"bonferroni"} or \code{"all"}. This plot-only construction supports one or two
continuous variables, at least ten observations and known finite coordinatewise
\code{ckerbound="fixed"} bounds. Observations and evaluation points must lie
strictly inside those bounds. Estimated \code{ckerbound="range"} bounds and
literal support endpoints are rejected; no jittering or target deletion occurs.
Bandwidths, bounds and targets stay fixed during inference. Core fitted values,
estimator standard errors, bandwidth selection and \code{se=FALSE} are unchanged.
An object fitted with \code{se=FALSE} can receive explicitly requested plot bands.

On normalized support, write \eqn{L_r=4K_r-6K_{2r}+4K_{3r}-K_{4r}}.
Here \eqn{r} is the dimensionless Chen kernel scale, not a distance
bandwidth. On \eqn{[0,1]},
\deqn{K_r(x,v)=\mathrm{Beta}(v;1+x/r,1+(1-x)/r).}
If the fitted physical bandwidth in coordinate \eqn{j} is \eqn{h_j^{phys}}
and the known support width is \eqn{w_j}, the operator uses
\eqn{\rho_{h,j}=(h_j^{phys}/w_j)^2}. Its auxiliary scale is
\eqn{b=n^{-7/17}} for one coordinate and \eqn{b=n^{-13/54}}
for two coordinates. In the two-coordinate formula, \eqn{K_{sr}}
is the product of the coordinate kernels at the same multiplier \eqn{s};
the four-term linear combination is taken after this product.

The product \eqn{L_{\rho_h}L_b} denotes operator composition:
\deqn{(L_{\rho_h}L_b)(x,v)=\int L_{\rho_h}(x,u)L_b(u,v)\,du.}
The fitted-bandwidth operator is outer and the auxiliary operator is inner.
It is neither a pointwise product nor the reversed composition. Thus the
residual kernel is \eqn{\Psi=L_{\rho_h}+L_b-L_{\rho_h}L_b}.

The corrected plot center is
\eqn{C=\widehat R+g_{\widehat\theta}}, where
\eqn{\widehat R=P_n\Psi-E_{g_{\widehat\theta}}\Psi}. The fitted guide has two
edge-shape parameters in one dimension, and four edge-shape parameters plus
an interaction in two dimensions. More explicitly,
\deqn{g_\theta(u)\propto\prod_j u_j^{\theta_{2j-1}}(1-u_j)^{\theta_{2j}}}
in one dimension, with the additional factor
\eqn{\exp\{\theta_5(2u_1-1)(2u_2-1)\}} in two dimensions.
Each edge parameter exceeds \eqn{-1}; the interaction is unconstrained.
The guide region uses \eqn{\eta_j=\log(1+\theta_j)} for edge parameters
and \eqn{\eta_5=\theta_5} for the interaction. With the empirical sandwich
covariance \eqn{\widehat V_\eta}, it takes the minimum and maximum guide
density over the Wald ellipsoid
\deqn{(\eta-\widehat\eta)'\widehat V_\eta^{-1}(\eta-\widehat\eta)
       \leq \chi^2_{p,1-\alpha/2},\qquad p=2\ \mathrm{or}\ 5.}
Half the nominal error probability is
assigned to a nonlinear guide confidence region and half to the residual
process. The residual covariance includes the guide-estimation projection
and cross terms. A positive ratio pilot supplies its transported scale.
This scale is not a complete normal standard error for the composite interval.

The positive ratio pilot is proportional to
\deqn{g_{\widehat\theta}(u)\frac{1}{n}\sum_i
       \frac{K_{\rho_p}(u,X_i)}{g_{\widehat\theta}(X_i)},}
normalized by numerical integration. Its Chen scale is
\eqn{\rho_p=n^{-2/5}} in one dimension and \eqn{\rho_p=n^{-1/3}}
in each coordinate in two dimensions. This is distinct from the auxiliary
operator scale \eqn{b}. Density values and interval limits are transformed
back to physical units by dividing by the product of support widths.

Asymptotic pointwise intervals use \eqn{z_{1-\alpha/4}}; simultaneous bands
calibrate the absolute Gaussian maximum under the complete residual covariance
on the finite requested grid. Asymptotic Bonferroni bands use
\eqn{z_{1-\alpha/(4m)}} for the residual at each of the \eqn{m} distinct
targets, retaining the same joint guide region at risk \eqn{\alpha/2}.
The guide risk is not divided again. With \code{band="all"}, one simultaneous
preparation supplies pointwise, simultaneous and Bonferroni limits; duplicate
targets do not inflate multiplicity and slices retain the whole family.
For both SD variants, every literal independent
multinomial resample refits the guide in the residual; the original empirical
residual scale is held fixed. This is not replicate studentization or a
bootstrap of the complete corrected center followed by a second guide region.
\code{"asymmetric-sd"} inverts the signed residual roots. \code{"symmetric-sd"}
uses their absolute roots, then takes the smallest outward-rounded interval
symmetric about the corrected center that contains the composite interval.
For \eqn{B} residual draws, the signed lower rank is
\eqn{k=\lfloor B\alpha/4\rfloor+1} and the upper rank is \eqn{B+1-k}.
Inversion subtracts the upper root at the lower endpoint and the lower root
at the upper endpoint. The absolute rank is
\eqn{k_a=\lceil B(1-\alpha/2)\rceil}. Pointwise ranks use each target's
roots; simultaneous signed ranks use the minimum and maximum across the
supplied finite family, and absolute ranks use its maximum absolute root.
The guide region is added once in each case. These beta-specific rules do not
change the existing Gaussian density bootstrap constructions.

For SD Bonferroni bands, let \eqn{m} be the number of distinct targets in the
whole supplied family, and \eqn{B} the unchanged number of bootstrap draws.
The symmetric absolute-root order statistic has rank
\eqn{\lceil B(1-\alpha/(2m))\rceil}. The asymmetric signed lower rank is
\eqn{k=\lfloor B\alpha/(4m)\rfloor+1}, with upper rank \eqn{B+1-k};
the signed interval is inverted before multiplying by the residual scale.
The guide region retains risk \eqn{\alpha/2}. Each complete marginal
Bonferroni interval is then united with its same-draw simultaneous interval.
This unconditional simultaneous hull cannot reduce coverage relative to the
existing simultaneous band, but can widen the result. It is not an estimate
of unseen bootstrap tails or a finite-sample population coverage guarantee.
The returned construction and resolution diagnostics identify this hull.

With \code{band="all"}, the existing pointwise and simultaneous endpoints are
retained exactly, together with this safeguarded Bonferroni companion. The
same preparation, guide region, resamples and guide refits are reused; no
additional bootstrap samples are generated. Very small Bonferroni tail
probabilities may select the observed extremes. Expected tail counts below
the existing resolution threshold produce one warning, while retaining
finite endpoints and the simultaneous hull. Requested draw counts, the
resampling law and existing standalone methods are unchanged.


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

Both the beta construction and the fixed-bound Gaussian FG-D construction
require supported built-in uniform and normal RNG kinds when they simulate
calibration draws. User-supplied uniform or normal generators are refused
before stream work; they are not silently replaced. To use these draw-based
routes, explicitly select a supported built-in combination, for example
\code{RNGkind("Mersenne-Twister", "Inversion", "Rejection")}, and set the
intended seed. Deterministic asymptotic pointwise inference does not create a
calibration stream. This statement does not certify arbitrary custom state
for other plotting constructions.

For a workspace or owned-file budget refusal, reduce the requested grid size
and/or the explicit \code{B} or \code{gaussian_draws}, as applicable to the
reported allocation. Reducing draws changes simulation resolution, not the
requested nominal level. No draw count, grid or method is reduced automatically.
Temporary-file creation also requires a writable temporary directory with
sufficient free space; the per-invocation storage budget is not a promise
about total R-process memory or other concurrent jobs.

The retained Beta-family, uniform, mixture and dependent two-dimensional
experiments show generally conservative finite-grid coverage. They do not
establish continuum coverage, honesty over unknown boundary smoothness,
uncertainty in estimated bounds or bandwidth-selection uncertainty. Very wide
near-edge intervals, especially for singular densities, remain possible.
The experimental status of bias-corrected inference applies to this method too.
Other bounded kernels, nearest-neighbour bandwidths, mixed data, dimensions
above two, derivatives, other bootstrap laws and \code{pmzsd} are outside
this beta construction.
Explicit saved density controls keep their original construction and scope.
}

\subsection{Known-bound univariate normal intervals}{
With omitted or \code{NULL} \code{bias_control}, known-fixed-bound order-2 beta
and Gaussian density plots support \code{errors="asymptotic"},
\code{center="bias-corrected"} and omitted \code{band} or \code{band="pmzsd"}
for one continuous variable. Fixed bandwidths, at least ten observations,
strictly interior data and targets, and \code{0 < alpha < 0.5} are required.
The original corrected center, observation operator and fitted bandwidth are
unchanged. The normal standard error includes the full guide-estimation and
cross-covariance terms, integrated under the normalized positive ratio pilot.
The guide fitted to that reference law is used only for the variance, not to
replace the original curve. Intervals are the corrected center plus or minus
\eqn{z_{1-\alpha/2}} times this complete scale; no guide-region risk split,
maximum calibration, variance floor, envelope or bootstrap resampling is used.

These are experimental marginal normal intervals. They do not inherit the
simultaneous coverage claims of the separate calibrated constructions. The
fitted pilot can substantially understate uncertainty even with complete
covariance calculations; sparse-tail deficits and residual mixture bias remain.
See the experimental-status section of \code{\link{plot.np}}. Its warning also
applies to the separate one-dimensional uniform and Epanechnikov normal
intervals described there. The calculation conditions on the selected
bandwidth, fixed bounds and targets; it does not add bandwidth-selection or
estimated-bound uncertainty. Negative centers and endpoints are retained.
The complete scale and its distinct normal/reference-law provenance are
available in the plot inference record. Core fitted values, estimator standard
errors and \code{se=FALSE} are unchanged.

Existing pointwise, simultaneous, Bonferroni, all-band, estimate-centered,
unbounded and explicit-control methods retain their constructions.
This univariate construction does not support bootstrap pmzsd or Gaussian
two-dimensional pmzsd. For the separate two-dimensional beta and compact normal constructions,
see also \code{\link{plot.np}}. It does not support nearest-neighbour bandwidths, mixed data, range-derived bounds or literal
support endpoints. In particular, the reference-law scale is not the standard
deviation of iid bootstrap refits.
}

\subsection{Known-bound compact normal intervals in two dimensions}{
With omitted or \code{NULL} \code{bias_control}, fixed second-order uniform and
Epanechnikov kernels support corrected asymptotic normal intervals for exactly
two continuous coordinates. Request omitted \code{band} or \code{band="pmzsd"}.
The fixed bandwidth, finite coordinatewise bounds, strictly interior data and
evaluation targets, at least ten observations and \code{0 < alpha < 0.5} are
required. The product local-cubic observation operator and corrected center are
unchanged. The variance integrates the complete influence function under the
normalized positive joint ratio pilot, including guide-score dependence and
all guide/kernel cross terms. A guide fitted to this reference law affects
only the variance, not the original corrected center.

Intervals use ordinary normal critical values without simulated draws, a
nonlinear guide region, a variance floor or a maximum calibration. They are
marginal intervals conditional on the bandwidth, fixed bounds and supplied
targets. The pilot can substantially understate uncertainty; sparse-target
information and residual bias can impair coverage even when numerical
calculations succeed. See the experimental-status section of
\code{\link{plot.np}}. There is no simultaneous coverage claim. Negative centers and endpoints are retained.

The calculation processes targets in bounded blocks, recomputing target-specific
compact-support quadrature. It retains the whole requested grid and checks the
192 MiB plot-workspace limit before allocation; this is not a bound on total R
process memory. The separate calibrated and SD q2 addition is described below.
Explicit controls, range-derived bounds, literal support endpoints, mixed data,
nearest-neighbour bandwidths, higher kernel orders and bootstrap pmzsd are not
admitted by this two-dimensional normal addition. Core estimates, estimator
standard errors and \code{se=FALSE} are unchanged.
}

\subsection{Two-dimensional compact calibrated and SD bands}{
With omitted density bias controls, fixed second-order uniform and Epanechnikov
kernels support corrected asymptotic \code{band="pointwise"} and
\code{"bonferroni"}, and independent \code{"symmetric-sd"} or
\code{"asymmetric-sd"} calibration with explicit \code{"pointwise"},
\code{"simultaneous"}, \code{"bonferroni"} or \code{"all"} bands.
The SD route requires explicit \code{bootstrap="inid"} and calibration in
\code{np_boot_control()}; an omitted or different bootstrap law is not remapped.
Use two continuous variables, known finite coordinatewise fixed bounds,
fixed bandwidths, at least ten strictly interior observations and interior
evaluation points, and \code{0 < alpha < 1}. Range bounds, literal support
endpoints, mixed variables, nearest-neighbour bandwidths, higher kernel orders,
explicit density controls and additional bootstrap laws are outside this addition.

The product compact correction and fitted center are unchanged. A joint
five-score guide retains cross-coordinate dependence. Residual uncertainty uses
both the complete ratio-pilot and empirical influence variances; endpoints take
their containing envelope and include a separate nonlinear guide region. The
risk allocation remains alpha/2 for the guide and alpha/2 for the residual.
Analytic pointwise and Bonferroni bands do not require Gaussian draws or a full
covariance matrix. Asymptotic simultaneous/all calibration is not admitted here.

SD bands use literal shared multinomial counts and one actual joint-guide refit
per draw. Roots are standardized by the original empirical residual scale;
both original scales enter the endpoint envelope. Symmetric-sd retains the
absolute containing hull; asymmetric-sd inverts signed residual tails.
Bonferroni companions use empirical tails with a simultaneous containing hull.
Small finite-draw tail resolution is reported explicitly, without extrapolating,
deleting draws or increasing the requested B. All-band output shares the same
preparation, draws and finite family as its individual components.

A complete common quadrature includes every supplied target's compact-support
cuts. Storage blocks do not define independent inference families or refit guides.
The 192 MiB additional inference-array and 512 MiB invocation-storage bounds are
checked; neither is a bound on the complete R process. Fitted values, estimator
standard errors, core \code{se=FALSE}, normal pmzsd and other existing methods
are unchanged. Negative returned centers and interval endpoints remain available.

Retained evidence is conservative rather than proof of universal nominal
coverage. In the retained six law/sample-size cells, calibrated simultaneous and
Bonferroni intervals covered all 128 originals per cell; SD evidence used only
16 originals per cell. Sampling uncertainty is correspondingly larger for SD.
Claims concern the supplied finite grid and do not establish continuum-wide
coverage or protect arbitrary sparse-tail distributions.
}

\subsection{Known-bound Gaussian density inference}{
With omitted or \code{NULL} \code{bias_control}, fixed second-order Gaussian
kernels with known \code{ckerbound="fixed"} bounds use the FG-D plot
construction for corrected asymptotic and independent symmetric-sd/asymmetric-sd
intervals. Request explicit pointwise, simultaneous, Bonferroni or
\code{band="all"}. As for the
known-bound beta construction, one or two continuous coordinates and at least
ten observations are required; data and evaluation points must be strictly
interior. Range-derived bounds, support endpoints, nearest-neighbour bandwidths,
higher kernel orders, derivatives, other bootstrap laws and pmzsd
are outside this construction. Explicit density controls, including the
legacy calibration alias, retain their existing meanings.

The query-normalized truncated Gaussian is corrected by a cubic donor-moment
projection on each normalized support coordinate. Coordinate operators are
multiplied in two dimensions. Write this operator as \eqn{A}; the corrected
plot center is \eqn{g_{\widehat\theta}+P_n A-E_{g_{\widehat\theta}}A}.
The guide and its nonlinear region have the same parameterization and
half-guide/half-residual risk allocation described above. Residual uncertainty
includes the guide-estimation projection, cross terms and mean products.

FG-D constructs two complete intervals, using respectively the positive
ratio-pilot residual covariance and the observation-based residual covariance,
then takes their outward union. Asymptotic simultaneous inference uses each
full covariance law on the entire finite target family, with common Gaussian
innovations; pointwise inference uses \eqn{z_{1-\alpha/4}}. Analytic Bonferroni
companions apply \eqn{z_{1-\alpha/(4m)}} separately to both residual scales
and combine their complete intervals. The guide region is unchanged. All-band
requests reuse one simultaneous preparation. SD inference refits
the guide in each actual multinomial sample and holds the original empirical
residual denominator fixed. The same residual roots are transported under each
scale. The asymmetric rule inverts their signed ranks; the symmetric rule forms
a containing symmetric hull for each complete interval before the union.
Taking the maximum of the two scales is not the asymmetric construction.
The guide region is added once per constituent interval. SD Bonferroni uses
the same marginal ranks and unconditional simultaneous hull described for
beta, separately for each complete constituent before their union; SD
\code{band="all"} preserves its existing pointwise and simultaneous components.

Returned limits and the corrected center are signed, without clipping or
renormalization. Bounds, evaluation grid and selected bandwidth are conditioned
on, not re-estimated. The full finite target family, requested draw count and
selected RNG are preserved; numerical or workspace failures are explicit.
Near-edge intervals may be very wide. Retained coverage experiments are
conservative over the investigated Beta-family and mixture cases, but do not
establish continuum or universal coverage. Core estimators and \code{se=FALSE}
are unchanged; an existing fitted object can request these plot intervals.
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
\dontrun{
npRmpi.init(nslaves = 1)
set.seed(123)
x <- data.frame(x = rnorm(80))
bw <- npudensbw(dat = x, bws = 0.45, bandwidth.compute = FALSE)
z <- plot(bw, xdat = x, errors = "bootstrap", bootstrap = "inid",
          center = "bias-corrected", band = "all", B = 399,
          neval = 3, output = "data")
z$inference[[1]]$primary
npRmpi.quit()
} % enddontrun
}
\keyword{nonparametric}
