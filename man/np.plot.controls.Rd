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
                calibration = c("default", "absolute", "separate-tails", "tail-envelope"))

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
  \item{calibration}{For iterated bias-corrected regression bootstrap empirical
    bands, \code{"separate-tails"} selects raw asymmetric empirical lower and
    upper tails; \code{"tail-envelope"} additionally protects with the
    absolute-root interval. \code{"default"} preserves existing dispatch and
    \code{"absolute"} explicitly requests absolute-root calibration. See Details.}
  \item{xtrim}{Optional two-element trim range for evaluation grids.}
  \item{xq}{Optional evaluation quantiles.}
  \item{slices}{Optional slice specification for routes that support it.}
  \item{style}{Interval display style.}
  \item{bar}{Bar style when \code{style = "bar"}.}
  \item{bar.num}{Optional number of bars.}
}
\details{
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
  Estimate-centered, legacy, normal-scale \code{pmzsd}, other estimator
  families and Gaussian asymptotic inference do not accept an explicit
  calibration selector. The Gaussian approximation is reflection-symmetric.
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
