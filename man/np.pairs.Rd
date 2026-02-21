\name{np.pairs}
\alias{np.pairs}
\alias{np.pairs.plot}
\title{Cross-Validated Pairs Plot (Helper Functions)}
\description{
Compute pairwise nonparametric regressions and densities for a set of
variables, then plot a pairs-style display with fitted smoothers.
}
\usage{
np.pairs(y_vars, y_dat, ...)
np.pairs.plot(pair_list)
}
\arguments{
  \item{y_vars}{
    character vector of column names in \code{y_dat}. If \code{y_vars}
    is named, the names are used as plot labels.
  }
  \item{y_dat}{
    data frame containing the variables listed in \code{y_vars}.
  }
  \item{...}{
    additional arguments passed to \code{\link{npudens}} and
    \code{\link{npreg}}.
  }
  \item{pair_list}{
    list returned by \code{np.pairs}.
  }
}
\value{
  \code{np.pairs} returns a list with components
  \code{y_vars}, \code{pair_names}, and \code{pair_kerns}.
  \code{np.pairs.plot} returns \code{NULL} (invisibly).
}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, and \code{\link{plot}} for plotting options.

On the diagonal, \code{npudens} is used to compute kernel density
estimates. Off-diagonal panels use \code{npreg} with residuals to draw
scatterplots and smoothers.
}
\examples{
\dontrun{
data("USArrests")
y_vars <- c("Murder", "UrbanPop")
names(y_vars) <- c("Murder Arrests per 100K", "Pop. Percent Urban")

pair_list <- np.pairs(
  y_vars = y_vars,
  y_dat = USArrests,
  ckertype = "epanechnikov",
  bwscaling = TRUE
)

np.pairs.plot(pair_list)
}
}
\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}\code{\link{npudens}}, \code{\link{npreg}}}
\keyword{hplot}
