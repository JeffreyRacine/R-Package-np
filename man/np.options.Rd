\name{np.options}
\alias{np.options}
\title{Global Package Options for \pkg{npRmpi}}
\description{
  Global options controlling selected computational and display behavior
  for the \pkg{npRmpi} package.
}
\details{
  For MPI startup/performance guidance (including message-passing tradeoffs and the manual-broadcast template), see \code{\link{npRmpi.init}} details and \code{inst/Rprofile}.
  Documentation guide: see \code{\link{np.kernels}} for kernels and
  \code{\link{plot}} for plotting options.

  The following options are recognized by \pkg{npRmpi}.

  \itemize{
    \item \code{np.messages} (logical): controls console/progress output.
    Default is \code{TRUE}.

    \item \code{np.plot.progress} (logical): controls bounded
    plot/bootstrap progress heartbeats on the master rank. Default is
    \code{TRUE}.

    \item \code{np.plot.progress.start.grace.sec} (numeric): delay before
    the first plot/bootstrap progress line is shown. Default is
    \code{0.75}.

    \item \code{np.plot.progress.interval.sec} (numeric): minimum elapsed
    time between plot/bootstrap heartbeat lines once progress reporting has
    started. Default is \code{0.5}.

    \item \code{np.plot.progress.max.intermediate} (integer): maximum
    number of mid-run plot/bootstrap heartbeat lines emitted between the
    initial start notice and final completion line. Default is \code{3}.

    \item \code{np.tree} (logical or \code{"auto"}): controls kd-tree
    acceleration for continuous variables. \code{TRUE} forces eligible tree
    routes on, \code{FALSE} disables them, and \code{"auto"} enables trees
    only when active continuous kernels are bounded support
    (\code{"epanechnikov"} or \code{"uniform"}). Default is \code{"auto"}.

    \item \code{np.categorical.compress} (logical): enables exact
    unique-profile compression for eligible all-categorical routes. This is
    distinct from \code{np.tree}: no continuous kd-tree is requested merely
    because categorical compression is enabled. The default is \code{TRUE}.
    Set \code{options(np.categorical.compress = FALSE, np.tree = FALSE)} to
    force the historical dense categorical route for debugging or timing
    comparisons.

    \item \code{np.powell.cache} (logical): enables exact caching of repeated
    integer nearest-neighbor candidates in Powell bandwidth searches for
    eligible continuous-only generalized/adaptive nearest-neighbor routes.
    Default is \code{TRUE}. Set \code{options(np.powell.cache = FALSE)} to
    force every repeated Powell candidate through the objective evaluator for
    diagnostic timing comparisons. This option is synchronized to MPI workers
    in autodispatch sessions. MADS/NOMAD solver caching is separate and is not
    controlled by this option.

    \item \code{np.largeh} (logical): enables continuous large-\eqn{h}
    shortcut evaluations. Default is \code{TRUE}. Set
    \code{options(np.largeh = FALSE)} to keep tree/dense routing unchanged
    while forcing continuous large-bandwidth objective rows through the
    ordinary kernel-summation machinery for diagnostic timing comparisons.

    \item \code{np.extendednn} (logical): enables extended nearest-neighbor
    bandwidth search and manually supplied nearest-neighbor bandwidths above
    the ordinary \eqn{n-1} neighbor-count limit for eligible
    generalized/adaptive nearest-neighbor routes. When \code{TRUE}, a
    nearest-neighbor search parameter above \eqn{n-1} scales the saturated
    nearest-neighbor bandwidth. Default is \code{TRUE}. This option is
    synchronized to MPI workers in autodispatch sessions.

    \item \code{np.largeh.rel.tol} (numeric): relative tolerance used by the
    continuous large-\eqn{h} shortcut. When all standardized distances for a
    continuous predictor are sufficiently close to zero, the corresponding
    kernel factor is approximated by \eqn{K(0)} to reduce repeated kernel
    evaluations. Default is \code{1e-3}. Valid range is \code{(0, 0.1)}.

    \item \code{np.largelambda} (logical): enables discrete near-upper
    bandwidth shortcut evaluations. Default is \code{TRUE}. Set
    \code{options(np.largelambda = FALSE)} to keep categorical routing
    unchanged while forcing near-upper discrete bandwidth rows through the
    ordinary category-comparison machinery for diagnostic timing comparisons.

    \item \code{np.disc.upper.rel.tol} (numeric): relative tolerance used by
    the discrete upper-bound shortcut for bandwidths near their feasible upper
    bounds. The near-upper check is applied relative to each kernel's own
    feasible upper bound (e.g., Aitchison-Aitken depends on category
    cardinality), with a tiny machine-precision floor for numerical robustness.
    When same/different-category kernel values are numerically close, the
    corresponding discrete kernel factor is treated as constant to reduce
    repeated category comparisons. Default is \code{1e-2}. Valid range is
    \code{(0, 0.5)}.

    \item \code{plot.par.mfrow} (logical): used by \code{\link{plot}} to
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
\seealso{\code{\link{npRmpi.init}}. 
  \code{\link{np.kernels}}, \code{\link{plot}}
  \code{\link{npRmpi}}, \code{\link{plot}}, \code{\link{options}}
}
\examples{
\dontrun{
npRmpi.init(nslaves=1)

old <- options(
  np.tree = "auto",
  np.categorical.compress = TRUE,
  np.powell.cache = TRUE,
  np.largeh = TRUE,
  np.extendednn = TRUE,
  np.largelambda = TRUE,
  np.messages = FALSE,
  np.largeh.rel.tol = 1e-3,
  np.disc.upper.rel.tol = 1e-2
)
on.exit(options(old), add = TRUE)
on.exit(npRmpi.quit(force=TRUE), add = TRUE)

## ... run bandwidth selection / estimation ...
}
}
\keyword{misc}
