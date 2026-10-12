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
                wild = c("rademacher", "mammen"),
                blocklen = NULL)

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
  \item{wild}{Wild-bootstrap multiplier distribution.}
  \item{blocklen}{Optional block length for block bootstrap methods.}
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
