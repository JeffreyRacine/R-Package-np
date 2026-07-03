\name{np.options}
\alias{np.options}
\title{Global Package Options for \pkg{np}}
\description{
  Global options controlling selected computational and display behavior
  for the \pkg{np} package.
}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels and
  \code{\link{plot}} for plotting options.

  The following options are recognized by \pkg{np}.

  \subsection{Output and Plot Progress}{
    \itemize{
      \item \code{np.messages} (logical): controls console/progress output.
      Default is \code{TRUE}.

      \item \code{np.plot.progress} (logical): controls bounded
      plot/bootstrap progress heartbeats. Default is \code{TRUE}.

      \item \code{np.plot.progress.start.grace.sec} (numeric): delay before
      the first plot/bootstrap progress line is shown. Default is
      \code{0.75}.

      \item \code{np.plot.progress.interval.sec} (numeric): minimum elapsed
      time between plot/bootstrap heartbeat lines once progress reporting has
      started. Default is \code{0.5}.

      \item \code{np.plot.progress.max.intermediate} (integer): maximum
      number of mid-run plot/bootstrap heartbeat lines emitted between the
      initial start notice and final completion line. Default is \code{3}.
    }
  }

  \subsection{Kernel Routing and Acceleration}{
    \itemize{
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

      \item \code{np.macMseries.accelerate} (logical or \code{"auto"}):
      controls macOS arm64 Accelerate/vDSP/vForce fast paths for eligible
      Gaussian-kernel and weighted-reduction computations. \code{"auto"}
      enables these paths only when they were compiled for macOS arm64 with
      Accelerate support; \code{TRUE} requests them when available; and
      \code{FALSE} disables them exactly. Default is \code{"auto"}. On
      non-macOS, non-arm64, or non-Accelerate builds, this option falls back to
      the ordinary portable code path. Accelerated and portable paths may use
      different math-library kernels or summation orders, so platform or option
      comparisons should use numerical tolerances; reserve exact equality for
      controlled runs using the same acceleration setting.
    }
  }

  \subsection{Search and Objective Caching}{
    \itemize{
      \item \code{np.objective.cache} (logical): enables package-side objective
      caching for repeated Powell bandwidth-search candidates. For all Powell
      routes using the shared C optimizer wrapper, exact repeated raw candidate
      vectors within one bandwidth solve reuse the previous objective value,
      including repeated candidates reached by later Powell restarts. For
      eligible continuous-only generalized/adaptive nearest-neighbor routes,
      this also enables the existing integer nearest-neighbor objective cache.
      Mixed continuous/categorical rows use only the exact full-candidate cache,
      so categorical smoothing parameters remain part of the cache key. Default
      is \code{TRUE}. Set \code{options(np.objective.cache = FALSE)} to force
      repeated Powell candidates through the objective evaluator for diagnostic
      timing comparisons. MADS/NOMAD solver caching is separate and is not
      controlled by this option, and the saturated extended-NN distance cache is
      governed by \code{np.extendednn}. Bandwidth summaries keep
      \code{Number of Function Evaluations} as the established cumulative
      optimizer evaluation count and report NOMAD, Powell, nearest-neighbor,
      degree-search, guarded-evaluation, and other fast-CV diagnostics on
      separate lines when the corresponding counters are available. These
      diagnostics describe their own optimizer/cache layer and are not intended
      to decompose the function-evaluation count arithmetically.

      \item \code{np.nomad.degree.start.policy} (character): controls automatic
      polynomial-degree starting values for NOMAD local-polynomial searches when
      \code{degree.start} is not supplied. The default,
      \code{"low_first_full_random"}, uses \code{degree.min} for the first start
      and draws later starts over the full user-specified
      \code{degree.min}/\code{degree.max} range. Other accepted values are
      \code{"mid_first_full_random"}, \code{"anchor_then_random"},
      \code{"spread_then_random"}, and \code{"random_full_only"}. Explicit
      \code{degree.start} values remain the first start.
    }
  }

  \subsection{Bandwidth Shortcuts}{
    \itemize{
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
      nearest-neighbor bandwidth. Default is \code{TRUE}.

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
      cardinality), with a tiny machine-precision floor for numerical
      robustness. When same/different-category kernel values are numerically
      close, the corresponding discrete kernel factor is treated as constant to
      reduce repeated category comparisons. Default is \code{1e-2}. Valid range
      is \code{(0, 0.5)}.
    }
  }

  \subsection{Plot Layout}{
    \itemize{
      \item \code{plot.par.mfrow} (logical): used by \code{\link{plot}} to
      determine whether plotting layout is automatically managed via
      \code{par(mfrow=...)}. If \code{NULL} (default behavior), \pkg{np} uses
      its internal plotting defaults.
    }
  }

  Option values can be set globally via \code{\link{options}} and restored
  with \code{\link{on.exit}} in scripts/functions for reproducibility.
}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\seealso{
  \code{\link{np.kernels}}, \code{\link{plot}}
  \code{\link{np}}, \code{\link{plot}}, \code{\link{options}}
}
\examples{
\dontrun{
old <- options(
  np.tree = "auto",
  np.categorical.compress = TRUE,
  np.macMseries.accelerate = "auto",
  np.objective.cache = TRUE,
  np.largeh = TRUE,
  np.extendednn = TRUE,
  np.largelambda = TRUE,
  np.messages = FALSE,
  np.largeh.rel.tol = 1e-3,
  np.disc.upper.rel.tol = 1e-2
)
on.exit(options(old), add = TRUE)

## ... run bandwidth selection / estimation ...
}
}
\keyword{misc}
