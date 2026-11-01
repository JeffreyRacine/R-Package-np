\name{np_plot_controls}
\alias{np_boot_control}
\alias{np_grid_control}
\alias{np_render_control}
\title{Validated plot control helpers}
\description{
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
                                "symmetric-sd", "asymmetric-sd", "smoothed-percentile", "legacy"))

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
  \item{blocklen}{Optional block length for block bootstrap methods.}
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
for the corresponding family \code{bias.control} with
\code{construction="legacy"}; it is not a fourth new calibration algorithm.
It requires corrected bootstrap inference, preserves that legacy route's
supported bands, laws and returned-error conventions, and rejects an explicitly
supplied nonlegacy bias control. Omitted calibration remains unchanged.

For density, \code{"smoothed-percentile"} uses complete unshifted corrected
curves, a data-selected outer smoothing pilot and 39 outer calibration samples.
The complete density and its correction vary jointly. Each fresh outer sample
gets its own ordinary CV bandwidth selection; the original fitted object and
selected bandwidth remain unchanged. The outer fits currently use the default
CV procedure even when the original fit used a manually supplied bandwidth,
normal-reference selection or a nondefault CV procedure; coverage under these
nondefault selection procedures has not been qualified. Inclusive family ranks select the
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
Gaussian derivative-three smoothing pilot, 39 fresh outer-sample CV fits,
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
wild/Mammen or the existing pairs/block laws. With \code{bias.control} omitted,
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
  \code{inference$resolution} diagnostics. Undefined variances still error.

  Raw asymmetric bands remain an optional method with limited finite-sample
  evidence. Exposed regression replays show improvements from the block-scale
  refinement, but remaining mean undercoverage and localized width inflation.
  The retained Mammen envelope has separate independent confirmation evidence
  for its original numeric settings; this does not qualify raw bands or
  other laws. No uniform coverage guarantee or dominance is asserted.
  Returned \code{inference} records requested/effective calibration, signed
  root orientation, raw/protected root limits where applicable, and the
  block-scale refinement. Certified deterministic targets, including a
  categorical contrast against its own reference, retain zero-width intervals
  and remain in the simultaneous family. The object does not retain full draws.

  Regression auxiliary bias controls are supplied separately through
  \code{\link{np_bias_control}}.

  These helpers validate names and values early. Raw lists are deliberately not
  accepted as substitutes, so misspelled detailed controls fail before plot
  computation or MPI fanout begins.

  Compound public plot arguments use dotted names, consistently with the
  established package API. Constructor function names and returned object
  fields are unchanged. Common rendering
  controls such as \code{data.overlay}, \code{data.rug}, \code{layout}, and
  \code{output} are supplied directly to \code{plot(...)} rather than through
  \code{np_render_control()}.

  Gaussian simultaneous variability bands use the direct plot argument
  \code{gaussian.draws} (default 9999), including supported density and
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
\examples{
\dontrun{
plot(fit,
     errors = "bootstrap",
     bootstrap = "inid",
     B = 9999,
     band = "pointwise",
     boot.control = np_boot_control(nonfixed = "exact"))
}
}
\keyword{hplot}
