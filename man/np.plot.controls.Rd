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
                  bar.num = NULL,
                  rug = FALSE,
                  overlay = FALSE,
                  rotate = FALSE)
}
\arguments{
  \item{nonfixed}{Bootstrap refit behavior for non-fixed bandwidth routes.}
  \item{wild}{Wild-bootstrap multiplier distribution.}
  \item{blocklen}{Optional block length for block bootstrap methods.}
  \item{xtrim}{Optional two-element trim range for evaluation grids.}
  \item{xq}{Optional evaluation quantiles.}
  \item{slices}{Optional slice specification for routes that support it.}
  \item{style}{Interval display style.}
  \item{bar}{Bar style when \code{style = "bar"}.}
  \item{bar.num}{Optional number of bars.}
  \item{rug}{Logical; request rug display where supported.}
  \item{overlay}{Logical; request data overlay where supported.}
  \item{rotate}{Logical; request rotation where supported.}
}
\details{
  These helpers validate names and values early. Raw lists are deliberately not
  accepted as substitutes, so misspelled detailed controls fail before plot
  computation or MPI fanout begins.

  Historical \code{plot.errors.*} arguments remain accepted as compatibility
  inputs in the 0.70-2 replacement tranche and are translated at the plot
  boundary. New code should prefer the direct scalar arguments and these
  helpers.
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
