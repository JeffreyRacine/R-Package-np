\name{np_plot_controls}
\alias{np_boot_control}
\alias{np_grid_control}
\alias{np_render_control}
\title{Validated plot control helpers}
\description{
\strong{Bias-corrected bounds are experimental (pre-beta research), incomplete,
and not recommended for general-purpose inference.} See the experimental-status
section below.

  Constructors for detailed \code{plot} controls used by the redesigned
  \pkg{npRmpi} plotting interface. Common choices such as \code{errors},
  \code{band}, \code{alpha}, \code{bootstrap}, and \code{B} remain direct
  \code{plot} arguments; these helpers collect less common options.
}
\usage{
np_boot_control(nonfixed = c("exact", "frozen"),
                wild = c("mammen", "rademacher"),
                blocklen = NULL,
                calibration = c("default", "absolute", "separate-tails",
                                "tail-envelope", "percentile-tails",
                                "symmetric-sd", "asymmetric-sd",
                                     "smoothed-percentile", "legacy"))

np_grid_control(xtrim = NULL, xq = NULL, slices = NULL)

np_render_control(style = c("band", "bar"),
                  bar = c("|", "I"),
                  bar.num = NULL)
}
\arguments{
  \item{nonfixed}{Bootstrap behavior for non-fixed bandwidth routes.
    \code{"exact"} reconstructs each duplicate-row resample and recomputes
    nearest-neighbour bandwidths; \code{"frozen"} reuses the bandwidths from
    the fitted object.}
  \item{wild}{Wild-bootstrap multiplier distribution. Mammen is the default;
    Rademacher remains available explicitly.}
  \item{blocklen}{Optional block length for block bootstrap methods.
  Corrected density providers reject fixed blocks whose length equals the
  training sample size: each resample then contains every observation exactly
  once and has no resampling variation. Geometric blocks are unaffected.}
  \item{calibration}{For bias-corrected bootstrap empirical bands,
    \code{"symmetric-sd"} selects the existing absolute-root construction and
    \code{"asymmetric-sd"} selects separate empirical tails with the family's
    existing scale treatment. These names apply to regression, unconditional
    density and unconditional distribution estimation. The experimental option
    \code{"smoothed-percentile"} selects calibrated unshifted corrected
    bootstrap curves, with family-specific resampling described below. Existing defaults and older names retain their
    meanings. See Details for supported configurations.}
  \item{xtrim}{Optional two-element trim range for evaluation grids.}
  \item{xq}{Optional evaluation quantiles.}
  \item{slices}{Optional slice specification for routes that support it.}
  \item{style}{Interval display style.}
  \item{bar}{Bar style when \code{style = "bar"}.}
  \item{bar.num}{Optional number of bars.}
}
\details{
The descriptive bootstrap choices require \code{errors="bootstrap"},
\code{center="bias-corrected"} and \code{band="pointwise"},
\code{"simultaneous"}, \code{"bonferroni"}, or \code{"all"}.
\code{"symmetric-sd"} retains \code{"absolute"}.
\code{"asymmetric-sd"} retains \code{"separate-tails"} for regression and
CDFs, and transported \code{"percentile-tails"} for density. The suffix
\code{sd} does not imply normal critical values: these empirical bands use
bootstrap calibration. Existing \code{"tail-envelope"} remains the conservative
hull with the absolute-root interval; it is distinct from the raw asymmetric
choice. Omitted/default controls preserve their earlier behavior.

The \code{"legacy"} selection retains the pre-iterated regression and
pre-joint density/CDF bias-corrected bootstrap constructions. It is shorthand
for the corresponding family \code{bias_control} with
\code{construction="legacy"}; it is not a fourth new calibration algorithm.
It requires corrected bootstrap inference, preserves that legacy route's
supported bands, laws and returned-error conventions, and rejects an explicitly
supplied nonlegacy bias control. Omitted calibration remains unchanged.

For known-fixed-bound beta density plots with omitted \code{bias_control},
\code{"symmetric-sd"} and \code{"asymmetric-sd"} select the separate refitted
residual construction described in \code{\link{np_udens_bias_control}}.
It requires independent resampling and explicit pointwise or simultaneous
bands. The asymmetric rule inverts signed residual roots; the symmetric rule
adds the guide region and takes a containing symmetric hull once. Gaussian
density tail transport and all existing calibration aliases retain their
meanings; omitted calibration does not silently select a beta SD variant.

For fixed-bound order-2 Gaussian density plots with omitted \code{bias_control},
the same two SD selections use FG-D: construct complete signed or symmetric
intervals under the pilot and empirical residual scales, then take their outward
union. The original empirical residual denominator remains fixed across actual
guide refits. Known fixed bounds, strictly interior observations/targets and
explicit pointwise/simultaneous bands are required. Range bounds are rejected;
\code{"legacy"} retains its existing route. See
\code{\link{np_udens_bias_control}}.

For density, \code{"smoothed-percentile"} uses complete unshifted corrected
curves, a data-selected outer smoothing pilot and 39 outer calibration samples.
The complete density and its correction vary jointly. The original fitted
object and selected bandwidth remain unchanged. Fresh outer samples repeat the
original bandwidth rule when its selection provenance is available: manual
physical bandwidths stay fixed, manual scale factors retain their units,
normal-reference bandwidths are recomputed, and CV repeats its retained
objective and search settings. Newly selected density and distribution
bandwidth objects retain these evaluated settings in \code{selection.procedure};
no stored call or calling environment is evaluated to recover them.
Older manual and normal-reference objects retain enough information to recover
those rules. Older CV objects without a complete descriptor retain historical
default-CV outer calibration with an aggregated warning and
\code{legacy-unresolved} provenance diagnostics. Refit such objects to obtain
faithful replay. Coverage for nondefault selection procedures remains
unqualified; matching the procedure is not itself a coverage theorem.
Inclusive family ranks select the
simultaneous content. Companion hulls preserve pointwise/simultaneous/Bonferroni
nesting and singleton equality. For a family with more than one grid point,
pointwise bands are nominal marginal percentile bands; the outer calibration
chooses simultaneous content, with a pointwise floor. The Bonferroni companion
contains both the nominal Bonferroni and calibrated simultaneous bands. For a
singleton all three coincide, including any support-completed endpoint.
Certified conditional bootstrap support supplies
the full-content endpoint when needed. No centering shift, SD transport,
zero clipping or forced point-estimate containment is applied.

This option requires the suggested package \pkg{ks}, fixed unbounded
second-order Gaussian continuous kernels, at least one continuous coordinate,
and admitted unordered Aitchison--Aitken or normalized ordered Li--Racine
coordinates. Ordered smoothing uses the integer lattice, including unobserved
gaps. Near-pooling categorical operators and their provenance are recorded.
It supports \code{inid}, or \code{fixed}/\code{geom} with block length one;
nonunit-block outer calibration is not implemented. The requested draw count
is retained, up to the prototype's 9999-draw limit. Low resolution returns every
requested band with diagnostics. Unattainable outer calibration is explicitly
flagged as diagnostic-only; support completion can yield wide intervals.
A support-completed band is distinct from empirical-tail calibration and can
be effectively uninformative. It is retained with a warning and diagnostics;
it does not mean that calibration unattainability has been silently discarded.
The Gaussian q=1 assessment is stronger than the current multivariate/mixed
coverage evidence. No uniform coverage or dominance claim is made.

Symmetric/asymmetric density SD methods retain their existing admitted
\code{inid}/\code{fixed}/\code{geom} laws and scope. Mixed-data block requests
require explicit block lengths. Asymmetric CDF bands support admitted fixed
bandwidth joint-process \code{inid}/\code{fixed}/\code{geom} inference;
numeric nearest-neighbor asymmetric CDF inference remains unqualified.
The experimental CDF \code{"smoothed-percentile"} option uses a
Gaussian derivative-three smoothing pilot, 39 fresh outer-sample fits using
the same bandwidth-procedure policy described above,
complete bias-corrected CDF curves and calibrated family ranks. Its current
scope is one numeric variable, fixed unbounded Gaussian2 bandwidths,
\code{inid} or literal block length one, and the complete untrimmed sample-range
evaluation grid. All three companion bands are returned when requested, with the same nominal
pointwise and calibrated simultaneous companion convention as density.
Useful population coverage for this CDF prototype has not been established;
it is exposed for experimentation. Fresh outer CV may select slightly different
bandwidths across serial and MPI fits, as with the existing estimator.

For regression, the common \code{"smoothed-percentile"} label selects the
E1 mean / E6 first-derivative calibrated percentile construction. This is a
fixed-design standardized empirical-residual bootstrap: residual innovations
are centered and normalized, then resampled with replacement and multiplied by
the fitted error scale. No kernel jitter is added to the residual distribution.
The smooth fitted pilot generates responses, and the complete correction
operator acts on every response draw. Curves are not recentered on the observed
corrected estimate. The 39 outer responses refresh the pilot, residual shape
and variance using the original selected geometry; they do not reselect
bandwidths. Derivative marginal and whole-family contents are calibrated
separately. For regression means with more than one grid point the pointwise
band is the nominal marginal percentile band, not an outer-calibrated marginal
band; it remains a floor for the calibrated simultaneous band. Derivative
pointwise content is separately outer-calibrated. Bonferroni companions include
the simultaneous envelope; singleton equality is enforced.

Selecting this regression option with \code{bootstrap} omitted uses
\code{bootstrap="residual"}; explicitly choosing a different bootstrap law is
an error. The residual law is exclusive to this option and does not replace
wild/Mammen or the existing pairs/block laws. With \code{bias_control} omitted,
the option uses \code{np_bias_control(adaptation="none")}. Initial scope is
fixed Gaussian2, one or two numeric predictors, local constant, local linear,
or Bernstein local polynomial degrees at most three, first derivatives,
16 correction iterations, \code{alpha=0.05} and \code{B<=9999}. Explicit
incompatible controls are rejected. Resident work is bounded by a forecast;
no training-square matrix or outer-by-draw-by-grid bank is retained.

This regression prototype has known localized-signal coverage failures, in
addition to conservative q=2 cells and heteroskedastic sensitivity/control-loss
failures. Applying the correction operator to responses generated by the smooth
pilot can place the bootstrap curves away from the observed corrected estimate.
The empirical residual-shape normalization can also be numerically unstable
when residuals are at floating-point resolution. This experimental opt-in is
retained for investigation; its public availability is not a general coverage
guarantee or a recommendation over the existing SD calibrations. All percentile methods
preserve support-completion and unattainability diagnostics. Percentile bounds
need not be symmetric about or contain the corrected estimate. Existing
symmetric/asymmetric SD methods and omitted calibration defaults are unchanged.
See \code{\link{np_udens_bias_control}} and \code{\link{np_udist_bias_control}}.

Gaussian simultaneous variability bands use the Gaussian-multiplier framework
of Chernozhukov, Chetverikov and Kato (2013, 2015). For the density/distribution
families, related simultaneous-inference research is due to Cattaneo, Jansson
and Ma (2024) and Cattaneo, Chandak, Jansson and Ma (2024). See
\code{\link{plot.np}} for the full published references, implementation scope
and estimator-specific qualifications.


  Asymmetric regression calibration requires \code{errors="bootstrap"},
  \code{center="bias-corrected"}, the iterated construction and
  \code{band="pointwise"}, \code{"simultaneous"}, \code{"bonferroni"}, or
  \code{"all"}. It supports the admitted \code{wild} (Mammen or explicit
  Rademacher), \code{inid}, \code{fixed}, and \code{geom} routes. Existing
  kernel, regularization and auxiliary-space restrictions still apply.
  Estimate-centered, legacy, normal-scale \code{pmzsd} and Gaussian asymptotic
  inference do not accept an explicit calibration selector. Families other
  than the regression, density and distribution routes described here retain their existing
  selection contracts. The Gaussian approximation is reflection-symmetric.
  Omitted selectors in previously saved controls retain their meaning.

  The raw method uses marginal lower/upper order statistics for pointwise
  bands and separate lower/upper family-extremum quantiles for simultaneous
  bands. Bonferroni is nested outward to contain the simultaneous band.
  Endpoints invert signed estimation errors, retaining within-panel and
  within-functional families. This is not the experimental marginal-rank
  balancing construction. The envelope is the interval hull with the same
  call's empirical absolute-root rectangle; it can be wider and does not
  repair poor raw tail calibration. Neither option uses a normal critical
  value or changes the selected bootstrap law.

  Corrected empirical \code{fixed}/\code{geom} bands use the pilot scale
  times the square root of the exact selected-block to independent score
  variance ratio, including with default absolute calibration. Literal draws
  and block selection are unchanged; pooled score studentization and
  \code{pmzsd} retain their existing owners. Low draw counts return every
  defined requested band with one aggregated resolution warning and retained
  \code{inference[[1]]$resolution} diagnostics. Undefined variances still error.
  For smoothed-percentile bands, \code{order_index} is the minimum empirical
  index actually used over the family after companion-band enlargement;
  \code{resolution_limited} indicates that this index is one.
  \code{nominal_tail_draws} retains its nominal meaning. Internal regression rank
  details also retain point-specific \code{used.indices}; mean-family
  \code{content} counts original curves contained by the final simultaneous
  band, which can differ from the selected calibration score. Analytic support
  enlargement is recorded separately and does not create additional draws.

  Raw asymmetric bands remain an optional method with limited finite-sample
  evidence. Exposed regression replays show improvements from the block-scale
  refinement, but remaining mean undercoverage and localized width inflation.
  The retained Mammen envelope has separate independent confirmation evidence
  for its original numeric settings; this does not qualify raw bands or
  other laws. No uniform coverage guarantee or dominance is asserted.
  The returned \code{\link{np_plot_data}} records selected absolute endpoints
  and compact method/resolution diagnostics. Signed-root workspaces and
  alternative protected/root arrays remain internal provider records. Certified deterministic targets, including a
  categorical contrast against its own reference, retain zero-width intervals
  and remain in the simultaneous family. The object does not retain full draws.

  Regression auxiliary bias controls are supplied separately through
  \code{\link{np_bias_control}}.

  These helpers validate names and values early. Raw lists are deliberately not
  accepted as substitutes, so misspelled detailed controls fail before plot
  computation or MPI fanout begins.

  Compound public plot arguments use snake_case. Control constructors retain
  their existing arguments; only plot-facing names change. Returned plot data
  use the common \code{\link{np_plot_data}} schema. Common rendering
  controls such as \code{data_overlay}, \code{data_rug}, \code{layout}, and
  \code{output} are supplied directly to \code{plot(...)} rather than through
  \code{np_render_control()}.

  Gaussian simultaneous variability bands use the direct plot argument
  \code{gaussian_draws} (default 9999), including supported density and
  distribution plots. This controls Gaussian calibration only; \code{B} and
  \code{np_boot_control()} continue to control bootstrap resampling. See
  \code{\link{plot.np}} for supported targets, covariance contracts and
  the references to Chernozhukov, Chetverikov and Kato (2013, 2015).

  Conditional density and distribution level bootstraps permit different
  supported continuous kernel families and orders on the X and Y sides,
  including beta on either side. Beta-containing routes remain limited to
  continuous X and Y with local-constant fitting and retain the resolved
  support bounds stored in the bandwidth object. The \code{nonfixed} choice
  changes nearest-neighbour bandwidth replay, not those fitted support bounds.
}
\value{
  A small classed list consumed by \code{plot} methods.
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

\examples{
\dontrun{
plot(fit,
     errors = "bootstrap",
     bootstrap = "inid",
     B = 9999,
     band = "pointwise",
     boot_control = np_boot_control(nonfixed = "exact"))
}
}
\keyword{hplot}
