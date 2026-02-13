\name{np.options}
\alias{np.options}
\title{Global Package Options for \pkg{npRmpi}}
\description{
  Global options controlling selected computational and display behavior
  for the \pkg{npRmpi} package.
}
\details{
  The following options are recognized by \pkg{npRmpi}.

  \itemize{
    \item \code{np.messages} (logical): controls console/progress output.
    Default is \code{TRUE}.

    \item \code{np.tree} (logical): enables kd-tree acceleration when
    supported by the selected kernel/operator combination. Default is
    \code{FALSE}.

    \item \code{np.largeh.rel.tol} (numeric): relative tolerance used by the
    continuous large-\eqn{h} shortcut. When all standardized distances for a
    continuous predictor are sufficiently close to zero, the corresponding
    kernel factor is approximated by \eqn{K(0)} to reduce repeated kernel
    evaluations. Default is \code{1e-3}. Valid range is \code{(0, 0.1)}.

    \item \code{np.disc.upper.rel.tol} (numeric): relative tolerance used by
    the discrete upper-bound shortcut for bandwidths near their feasible upper
    bounds. When same/different-category kernel values are numerically close,
    the corresponding discrete kernel factor is treated as constant to reduce
    repeated category comparisons. Default is \code{1e-2}. Valid range is
    \code{(0, 0.5)}.

    \item \code{plot.par.mfrow} (logical): used by \code{\link{npplot}} to
    determine whether plotting layout is automatically managed via
    \code{par(mfrow=...)}. If \code{NULL} (default behavior), \pkg{npRmpi} uses
    its internal plotting defaults.
  }

  Option values can be set globally via \code{\link{options}} and restored
  with \code{\link{on.exit}} in scripts/functions for reproducibility.
}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
  \code{\link{npRmpi}}, \code{\link{npplot}}, \code{\link{options}}
}
\examples{
\dontrun{
old <- options(
  np.tree = TRUE,
  np.messages = FALSE,
  np.largeh.rel.tol = 1e-3,
  np.disc.upper.rel.tol = 1e-2
)
on.exit(options(old), add = TRUE)

## ... run bandwidth selection / estimation ...
}
}
\keyword{misc}
