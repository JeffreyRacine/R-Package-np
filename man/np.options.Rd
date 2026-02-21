\name{np.options}
\alias{np.options}
\title{Global Package Options for \pkg{npRmpi}}
\description{
  Global options controlling selected computational and display behavior
  for the \pkg{npRmpi} package.
}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels and
  \code{\link{npplot}} for plotting options.

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
    bounds. The near-upper check is applied relative to each kernel's own
    feasible upper bound (e.g., Aitchison-Aitken depends on category
    cardinality), with a tiny machine-precision floor for numerical robustness.
    When same/different-category kernel values are numerically close, the
    corresponding discrete kernel factor is treated as constant to reduce
    repeated category comparisons. Default is \code{1e-2}. Valid range is
    \code{(0, 0.5)}.

    \item \code{np.groupcv.fast} (logical): enables the grouped
    \eqn{G \ll n}{G << n} cross-validation fast path for regression bandwidth
    selection when the all-large-\eqn{h} gate is active and grouped discrete
    structure can be exploited. This can substantially reduce computation in
    discrete-heavy settings, while potentially yielding tiny objective-level
    floating-point differences relative to the legacy accumulation order.
    Default is \code{TRUE}.

    \item \code{plot.par.mfrow} (logical): used by \code{\link{npplot}} to
    determine whether plotting layout is automatically managed via
    \code{par(mfrow=...)}. If \code{NULL} (default behavior), \pkg{npRmpi} uses
    its internal plotting defaults.

    \item \code{npRmpi.autodispatch} (logical): when \code{TRUE}, eligible
    \code{np*} calls are auto-dispatched across MPI ranks without explicit
    \code{mpi.bcast.cmd(...)} wrapping. Default is \code{FALSE}.
    For formula interfaces under autodispatch, provide explicit
    \code{data=} (or explicit \code{xdat}/\code{ydat}-style arguments) to
    avoid unresolved-symbol failures on slave ranks.
  }

  Option values can be set globally via \code{\link{options}} and restored
  with \code{\link{on.exit}} in scripts/functions for reproducibility.
}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
  \code{\link{np.kernels}}, \code{\link{npplot}}
  \code{\link{npRmpi}}, \code{\link{npplot}}, \code{\link{options}}
}
\examples{
\dontrun{
npRmpi.init(nslaves=1)
options(npRmpi.autodispatch=TRUE, np.messages=FALSE)

old <- options(
  np.tree = TRUE,
  np.messages = FALSE,
  np.largeh.rel.tol = 1e-3,
  np.disc.upper.rel.tol = 1e-2,
  np.groupcv.fast = TRUE
)
on.exit(options(old), add = TRUE)
on.exit(npRmpi.quit(force=TRUE), add = TRUE)

## ... run bandwidth selection / estimation ...
}
}
\keyword{misc}
