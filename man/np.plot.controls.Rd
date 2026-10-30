\name{np_plot_controls}
\alias{np_boot_control}
\alias{np_grid_control}
\alias{np_render_control}
\title{Validated plot control helpers}
\description{
  Constructors for detailed \code{plot} controls used by the redesigned
  \pkg{np} plotting interface. Common choices such as \code{errors},
  \code{band}, \code{alpha}, \code{bootstrap}, and \code{B} remain direct
  \code{plot} arguments; these helpers collect less common options.
}
\usage{
np_boot_control(nonfixed = c("exact", "frozen"),
                wild = c("mammen", "rademacher"),
                blocklen = NULL,
                calibration = c("default", "absolute", "separate-tails",
                                "tail-envelope", "percentile-tails",
                                "symmetric-sd", "asymmetric-sd", "smoothed-percentile"))

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
    density and unconditional distribution estimation. For density only,
    \code{"smoothed-percentile"} selects calibrated unshifted corrected
    smoothed-bootstrap curves. Existing defaults and older names retain their
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

For density, \code{"smoothed-percentile"} uses complete unshifted corrected
curves, a data-selected outer smoothing pilot and 39 outer calibration samples.
The complete density and its correction vary jointly. Each fresh outer sample
gets its own ordinary CV bandwidth selection; the original fitted object and
selected bandwidth remain unchanged. Inclusive family ranks select the
simultaneous content. Companion hulls preserve pointwise/simultaneous/Bonferroni
nesting and singleton equality. Certified conditional bootstrap support supplies
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
The Gaussian q=1 assessment is stronger than the current multivariate/mixed
coverage evidence. No uniform coverage or dominance claim is made.

Symmetric/asymmetric density SD methods retain their existing admitted
\code{inid}/\code{fixed}/\code{geom} laws and scope. Mixed-data block requests
require explicit block lengths. Asymmetric CDF bands support admitted fixed
bandwidth joint-process \code{inid}/\code{fixed}/\code{geom} inference;
numeric nearest-neighbor asymmetric CDF inference remains unqualified.
\code{"smoothed-percentile"} is not yet implemented for CDFs or regression.
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
  computation begins.

  Compound public plot arguments use dots, consistently with the package API.
  Constructor function names and returned object fields are unchanged. Common rendering
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
