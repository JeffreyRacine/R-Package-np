%
\name{npcdens}
\alias{npcdens}
\alias{npcdens.conbandwidth}
\alias{npcdens.default}
\alias{npcdens.formula}
\title{Kernel Conditional Density Estimation with Mixed Data Types}

\description{
  \code{npcdens} computes kernel conditional density estimates on
  \eqn{p+q}-variate evaluation data, given a set of training data (both
  explanatory and dependent) and a bandwidth specification (a
  \code{conbandwidth} object or a bandwidth vector, bandwidth type, and
  kernel type) using the method of Hall, Racine, and Li (2004).
  The data may be continuous, discrete (unordered and ordered
  factors), or some combination thereof.
}

\usage{
npcdens(bws, \dots)

\method{npcdens}{formula}(bws, data = NULL, newdata = NULL, \dots)


\method{npcdens}{conbandwidth}(bws,
        txdat = stop("invoked without training data 'txdat'"),
        tydat = stop("invoked without training data 'tydat'"),
        exdat,
        eydat,
        gradients = FALSE,
        proper = FALSE,
        proper.method = c("project"),
        proper.control = list(),
        \dots)

\method{npcdens}{default}(bws, txdat, tydat, nomad = FALSE, \dots)

}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface, and training data.
}

  \item{bws}{
    a bandwidth specification. This can be set as a \code{conbandwidth}
    object returned from a previous invocation of
    \code{\link{npcdensbw}}, or as a \eqn{p+q}-vector of bandwidths,
    with each element \eqn{i} up to \eqn{i=q} corresponding to the
    bandwidth for column \eqn{i} in \code{tydat}, and each element
    \eqn{i} from \eqn{i=q+1} to \eqn{i=p+q} corresponding to the
    bandwidth for column \eqn{i-q} in \code{txdat}. If specified as a
    vector, then additional arguments will need to be supplied as
    necessary to specify the bandwidth type, kernel types, training
    data, and so on.
  }


  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
   \code{\link{npcdensbw}} was called.
  }


  \item{txdat}{
    a \eqn{p}-variate data frame of sample realizations of explanatory
    data (training data). Defaults to the training data used to
    compute the bandwidth object.
  }


  \item{tydat}{
    a \eqn{q}-variate data frame of sample realizations of dependent
    data (training data). Defaults to the training data used to
    compute the bandwidth object.
  }


\subsection{Bandwidth Search Shortcut}{
This argument passes the recommended automatic local-polynomial NOMAD preset to \code{\link{npcdensbw}} when bandwidths are computed inside \code{npcdens}.
}

  \item{nomad}{
    logical shortcut passed through to \code{\link{npcdensbw}} when
    bandwidths are computed inside \code{npcdens}. When \code{TRUE},
    the bandwidth route fills any missing values among \code{regtype},
    \code{search.engine}, \code{degree.select},
    \code{bernstein.basis}, \code{degree.min}, \code{degree.max},
    \code{degree.verify}, and \code{bwtype} with the recommended
    automatic LP NOMAD preset documented in \code{\link{npcdensbw}}.
    Additional NOMAD tuning arguments such as \code{nomad.nmulti} may
    also be supplied through \code{\dots}; \code{nmulti} remains the
    outer restart count while \code{nomad.nmulti} controls inner
    \code{crs::snomadr()} multistarts within each outer restart.
    After fitting, inspect \code{fit$bws$nomad.shortcut} on the
    returned object \code{fit} to see the normalized shortcut metadata.
  }


\subsection{Evaluation Data And Returned Quantities}{
These arguments control where the fitted conditional density is evaluated and which estimates are returned.
}

  \item{exdat}{
    a \eqn{p}-variate data frame of explanatory data on which
    conditional densities will be evaluated. By default,
    evaluation takes place on the data provided by \code{txdat}.
  }


  \item{eydat}{
    a \eqn{q}-variate data frame of dependent data on which conditional
    densities will be evaluated. By default,
    evaluation takes place on the data provided by \code{tydat}.
  }


  \item{gradients}{
    a logical value specifying whether to return estimates of the
    gradients at the evaluation points. Defaults to \code{FALSE}.
  }


  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.
  }



\subsection{Fit Properization Controls}{
These arguments control optional post-estimation properization of the fitted conditional density.
}

  \item{proper}{
    a logical value specifying whether to post-process the estimated
    conditional density so that it integrates to one over the
    evaluation grid. Defaults to \code{FALSE}.
  }


  \item{proper.control}{
    a named list of control parameters for properization. Supported
    entries are \code{tol}, \code{grid.check}, \code{store.raw}, and
    \code{fail.on.unsupported}.
  }



  \item{proper.method}{
    the properization method. Currently only
    \code{"project"} is supported.
  }


\subsection{Additional Arguments}{
Further arguments are passed to \code{\link{npcdensbw}} when bandwidths are computed internally, or used to interpret a numeric \code{bws} vector.
}

  \item{\dots}{
    additional arguments supplied to \code{\link{npcdensbw}} when
    \code{npcdens} computes bandwidths internally, or arguments needed
    to interpret a numeric \code{bws} vector. This is where bandwidth
    selection controls such as \code{bwmethod}, \code{bwtype},
    kernel/support controls such as \code{cxkertype},
    \code{cykertype}, \code{cxkerbound}, and \code{cykerbound},
    search controls such as \code{nmulti},
    \code{scale.factor.search.lower}, and \code{nomad.nmulti}, and
    local-polynomial controls such as \code{regtype}, \code{degree},
    \code{basis}, and \code{bernstein.basis} are supplied. See
    \code{\link{npcdensbw}} for the complete bandwidth-selection
    argument surface.

  }

}
\details{
  Documentation guide: see \code{\link{npcdensbw}} for bandwidth selection and search controls, \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, \code{\link{plot}}, \code{\link{plot.np}} for plotting options, and \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See \code{\link{npRmpi.init}} details for performance tradeoffs (message passing/startup mode) and the \code{inst/Rprofile} manual-broadcast template.

  When \code{bws} is omitted, the formula and default methods call
  \code{\link{npcdensbw}} first and pass bandwidth-selection arguments
  from \code{\dots} to that call. When \code{bws} is already a
  \code{conbandwidth} object, \code{npcdens} estimates with the stored
  bandwidth metadata in that object.

  Argument groups for bandwidth selection are documented on
  \code{\link{npcdensbw}}. The most common workflow is to initialize MPI
  execution if needed, choose data and bandwidth inputs, then bandwidth
  criterion and representation, then kernel/support controls, numerical
  search controls, bounded \code{cv.ls} quadrature controls if relevant,
  and finally local-polynomial/NOMAD controls for polynomial-adaptive
  fits.

  For S3 plotting help, use \code{methods("plot")} and query
  class-specific help topics such as \code{?plot.npregression} and
  \code{?plot.rbandwidth}. You can inspect implementations with
  \code{getS3method("plot","npregression")}.

  \code{npcdens} implements a variety of methods for estimating
  multivariate conditional distributions (\eqn{p+q}-variate) defined
  over a set of possibly continuous and/or discrete (unordered, ordered)
  data. The approach is based on Li and Racine (2004) who employ
  \sQuote{generalized product kernels} that admit a mix of continuous
  and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating
  the density at the point \eqn{x}. Generalized nearest-neighbor
  bandwidths change with the point at which the density is estimated,
  \eqn{x}. Fixed bandwidths are constant over the support of \eqn{x}.

  Training and evaluation input data  may be a
  mix of continuous (default), unordered discrete (to be specified in
  the data frames using \code{\link{factor}}), and ordered discrete (to be
  specified in the data frames using \code{\link{ordered}}). Data can be
  entered in an arbitrary order and data types will be detected
  automatically by the routine (see \code{\link{npRmpi}} for details).

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

  For practitioners who want the recommended automatic LP NOMAD route
  without spelling out all LP tuning arguments,
  \code{npcdens(..., nomad=TRUE)} and \code{npcdensbw(..., nomad=TRUE)}
  expand missing settings to the same documented preset. Explicit
  incompatible settings fail fast rather than being silently rewritten.

}
\value{

  \code{npcdens} returns a \code{condensity} object. The generic
  accessor functions \code{\link{fitted}}, \code{\link{se}}, and
  \code{\link{gradients}}, extract estimated values, asymptotic standard
  errors on estimates, and gradients, respectively, from the returned
  object. Furthermore, the functions \code{\link{predict}},
  \code{\link{summary}} and \code{\link{plot}} support objects of both
  classes. The returned objects have the following components:

  \item{xbw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    explanatory data, \code{txdat} }
  \item{ybw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    dependent data, \code{tydat} }
  \item{xeval}{ the evaluation points of the explanatory data }
  \item{yeval}{ the evaluation points of the dependent data }
  \item{condens}{ estimates of the conditional density at the evaluation
    points }
  \item{conderr}{ standard errors of the conditional density
    estimates }
  \item{congrad}{ if invoked with \code{gradients = TRUE}, estimates of
    the gradients at the evaluation points }
  \item{congerr}{ if invoked with \code{gradients = TRUE}, standard
    errors of the gradients at the evaluation points }
  \item{log_likelihood}{ log likelihood of the conditional density estimate }

}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and J.S. Racine and Q. Li (2004), \dQuote{Cross-validation and the
    estimation of conditional probability densities,} Journal of the
  American Statistical Association, 99, 1015-1026.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,} Cambridge
  University Press.

  Scott, D.W. (1992), \emph{Multivariate Density Estimation. Theory,
  Practice and Visualization,} New York: Wiley.

  Silverman, B.W. (1986), \emph{Density Estimation,} London: Chapman and
  Hall.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
    for discrete distributions,}  Biometrika, 68, 301-309.

}
\author{
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
  \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.
}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}, \code{\link{plot.np}}
  \code{\link{npudens}}
}
\examples{
\dontrun{
## Not run in checks: excluded to keep MPI examples stable and check times short.
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi_getting_started", package = "npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.quit(force=TRUE)` then restart.

force.run <- nzchar(Sys.getenv("NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK"))

in.check <- !force.run && (
            nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
            nzchar(Sys.getenv("_R_CHECK_INTERNALS2_")) ||
            nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_")) ||
            nzchar(Sys.getenv("_R_CHECK_EXAMPLE_TIMING_THRESHOLD_"))
          )

if (!in.check) {
npRmpi.init(nslaves=1)

data("Italy")

bw <- npcdensbw(formula=gdp~ordered(year), data=Italy)

fhat <- npcdens(bws=bw)

summary(fhat)

## Variations on local polynomial conditional density estimation with
## proper = TRUE.

Italy2 <- within(Italy, {
  year <- as.numeric(as.character(year))
})

## Plot only: make the plotted surface proper on the plot evaluation grid.

fhat <- npcdens(gdp ~ year, data = Italy2,
                regtype = "lp", degree = 3, nmulti = 1)

plot(fhat, proper = TRUE)

## Fit an object whose fitted values are themselves proper.

ctrl_fit <- list(
  mode = "slice",
  apply = "fitted",
  slice.grid.size = 101L,
  slice.extend.factor = 0.1
)

fhat_fit <- npcdens(
  gdp ~ year,
  data = Italy2,
  regtype = "lp",
  degree = 3,
  nmulti = 1,
  proper = TRUE,
  proper.control = ctrl_fit
)

fit_proper <- fitted(fhat_fit)
fit_raw <- fhat_fit$condens.raw

## Display the repaired and raw fitted values for cases where the raw
## fitted density is negative.

head(cbind(fit_proper, fit_raw)[which(fit_raw < 0), ])

## Predict on a common explicit y-grid for several years, and render
## those predictions proper.

g.grid <- seq(min(Italy2$gdp), max(Italy2$gdp), length.out = 200)

nd_grid <- expand.grid(
  gdp = g.grid,
  year = c(1955, 1975, 1995)
)

pred_grid <- predict(fhat, newdata = nd_grid, proper = TRUE)

## Predict on paired rows with different gdp grids by year, and still
## make the predictions proper via slice mode.

g1 <- seq(quantile(Italy2$gdp, 0.10),
          quantile(Italy2$gdp, 0.60), length.out = 60)
g2 <- seq(quantile(Italy2$gdp, 0.30),
          quantile(Italy2$gdp, 0.90), length.out = 35)

nd_slice <- rbind(
  data.frame(gdp = g1, year = rep(1960, length(g1))),
  data.frame(gdp = g2, year = rep(1985, length(g2)))
)

pred_slice <- predict(
  fhat,
  newdata = nd_slice,
  proper = TRUE,
  proper.control = list(mode = "slice")
)

## One object that carries properization for fitted values and for later
## predict() calls.

ctrl_both <- list(
  mode = "slice",
  apply = "both",
  slice.grid.size = 101L,
  slice.extend.factor = 0.1
)

fhat_both <- npcdens(
  gdp ~ year,
  data = Italy2,
  regtype = "lp",
  degree = 3,
  nmulti = 1,
  proper = TRUE,
  proper.control = ctrl_both
)

fit_both <- fitted(fhat_both)
pred_both <- predict(
  fhat_both,
  newdata = nd_slice,
  proper.control = ctrl_both
)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

## Note: on some systems (notably macOS+MPICH), repeatedly spawning and
## tearing down slaves in the same R session can lead to hangs/crashes.
## npRmpi may therefore keep slave daemons alive by default and
## `npRmpi.quit()` performs a "soft close". Use `force=TRUE` to
## actually shut down the slaves.
##
## You can disable reuse via `options(npRmpi.reuse.slaves=FALSE)` or by
## setting the environment variable `NP_RMPI_NO_REUSE_SLAVES=1` before
## loading the package.

npRmpi.quit()               ## soft close (may keep slaves alive)
## npRmpi.quit(force=TRUE)  ## hard close
} else {
  message("Skipping MPI spawn in check context (set NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1).")
}
} % enddontrun
}
\keyword{nonparametric}
