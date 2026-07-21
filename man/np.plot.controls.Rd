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
                  bar_num = NULL)
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
  \item{bar_num}{Optional number of bars.}
}
\details{
  These helpers validate names and values early. Raw lists are deliberately not
  accepted as substitutes, so misspelled detailed controls fail before plot
  computation or MPI fanout begins.

  New public plot arguments use \code{snake_case}. Dotted argument names are
  reserved for S3 methods and internal plot-engine plumbing. Common rendering
  controls such as \code{data_overlay}, \code{data_rug}, \code{layout}, and
  \code{output} are supplied directly to \code{plot(...)} rather than through
  \code{np_render_control()}.

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
     B = 1999,
     band = "pointwise",
     boot_control = np_boot_control(nonfixed = "exact"))
}
}
\keyword{hplot}
