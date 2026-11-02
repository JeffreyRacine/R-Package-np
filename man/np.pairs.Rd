\name{np.pairs}
\alias{np.pairs}
\alias{np.pairs.plot}
\title{Cross-Validated Pairs Plot (Helper Functions)}
\description{
Compute pairwise nonparametric regressions and densities for a set of
variables, then plot a pairs-style display with fitted smoothers.
}
\usage{
np.pairs(y.vars, y.dat, ...)
np.pairs.plot(pair.list)
}
\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the variables, data, and pair specifications to plot.
}

  \item{pair.list}{
    list returned by \code{np.pairs}.
  }


  \item{y.dat}{
    data frame containing the variables listed in \code{y.vars}. With two or
    more distinct selected variables, each must be a numeric (integer or double)
    response vector, because each appears as an off-diagonal regression outcome.
    Factors, ordered factors and logical responses are rejected before fitting;
    their integer codes are not treated as numeric outcomes.
  }


  \item{y.vars}{
    character vector of column names in \code{y.dat}. If \code{y.vars}
    is named, the names are used as plot labels.
  }


\subsection{Additional Arguments}{
Further graphical arguments are passed through to plotting methods.
}

  \item{...}{
    additional arguments passed to \code{\link{npudens}} and
    \code{\link{npreg}}.
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
y.vars <- c("Murder", "UrbanPop")
names(y.vars) <- c("Murder Arrests per 100K", "Pop. Percent Urban")

pair.list <- np.pairs(
  y.vars = y.vars,
  y.dat = USArrests,
  ckertype = "epanechnikov",
  bwscaling = TRUE
)

np.pairs.plot(pair.list)
}
}
\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}\code{\link{npudens}}, \code{\link{npreg}}}
\keyword{hplot}
