\name{npudist}
\alias{npudist}
\alias{npudist.dbandwidth}
\alias{npudist.default}
\alias{npudist.formula}

\title{ Kernel Distribution Estimation with Mixed Data Types }
\description{

  \code{npudist} computes kernel unconditional cumulative distribution
  estimates on evaluation data, given a set of training data and a
  bandwidth specification (a \code{dbandwidth} object or a bandwidth
  vector, bandwidth type, and kernel type) using the method of Li, Li
  and Racine (2017).

}
\usage{
npudist(bws, \dots)

\method{npudist}{formula}(bws, 
        data = NULL, 
        newdata = NULL, 
        \dots)

\method{npudist}{dbandwidth}(bws,
        tdat = stop("invoked without training data 'tdat'"),
        edat,
        \dots)

\method{npudist}{default}(bws, 
        tdat, 
        \dots)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface, and training data.
}

  \item{bws}{
    a \code{dbandwidth} specification. This can be set as a \code{dbandwidth}
    object returned from an invocation of \code{\link{npudistbw}}, or as a
    \eqn{p}-vector of bandwidths, with an element for each variable in the
    training data. If specified as a vector, then additional arguments
    will need to be supplied as necessary to change them from the
    defaults to specify the bandwidth type, kernel types, training data,
    and so on.
  }



  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npudistbw}} was called.
  }


  \item{tdat}{
    a \eqn{p}-variate data frame of sample realizations (training data)
    used to estimate the cumulative distribution. Defaults to the corresponding
    training data recoverable from the bandwidth object in the current
    session.
  }


\subsection{Evaluation Data}{
These arguments control where the fitted cumulative distribution is evaluated.
}

  \item{edat}{
    a \eqn{p}-variate data frame of cumulative distribution evaluation points. By default,
    evaluation takes place on the data provided by \code{tdat}.
  }


  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.
  }


\subsection{Additional Arguments}{
Further arguments are passed to \code{\link{npudistbw}} when bandwidths are computed internally, or used to interpret a numeric \code{bws} vector.
}

  \item{\dots}{
    additional arguments supplied to \code{\link{npudistbw}} when
    \code{npudist} computes bandwidths internally, or arguments needed
    to interpret a numeric \code{bws} vector. This is where bandwidth
    selection controls such as \code{bwmethod}, \code{bwtype},
    \code{bwscaling}, kernel/support controls such as \code{ckertype},
    \code{ckerorder}, and \code{ckerbound}, categorical kernel controls
    such as \code{okertype}, and search controls such as
    \code{nmulti}, \code{ngrid}, and
    \code{scale.factor.search.lower} are supplied. See
    \code{\link{npudistbw}} for the complete bandwidth-selection
    argument surface.
  }
}
\details{
  Documentation guide: see \code{\link{npudistbw}} for bandwidth
  selection and search controls, \code{\link{np.kernels}} for kernels,
  \code{\link{np.options}} for global options, \code{\link{plot}},
  \code{\link{plot.np}} for plotting options, and
  \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See
  \code{\link{npRmpi.init}} details for performance tradeoffs (message
  passing/startup mode) and the \code{inst/Rprofile} manual-broadcast
  template.

  When \code{bws} is omitted, the formula and default methods call
  \code{\link{npudistbw}} first and pass bandwidth-selection arguments
  from \code{\dots} to that call. When \code{bws} is already a
  \code{dbandwidth} object, \code{npudist} estimates with the stored
  bandwidth metadata in that object.

  With finite continuous bounds, \code{npudist} uses the
  observation-centered truncated-kernel CDF documented in
  \code{\link{np.kernels}}. It is exactly zero and one at the declared
  lower and upper endpoints and has the linear uniform-CDF limit as the
  bandwidth diverges.

  The beta associated CDF is available at orders 2, 4, 6, and 8 with manual
  or cross-validated metric bandwidths. Fixed, generalized-NN, and adaptive-NN
  modes are supported. Set \code{ckertype="beta"} and use either
  \code{ckerbound="fixed"} with finite \code{ckerlb}/\code{ckerub}, or
  \code{ckerbound="range"} for the exact training-sample range. The data must
  be continuous and all training and evaluation values must lie within the
  resulting supports. Categorical products and normal-reference selection are
  not supported.

  For observation \eqn{x_i}, target \eqn{x}, support \eqn{[a,b]},
  \eqn{L=b-a}, \eqn{z_i=(x_i-a)/L}, \eqn{t=(x-a)/L}, and
  \eqn{\tau=(L/h)^2}, the contribution is
  \deqn{I_t(1+z_i\tau,1+(1-z_i)\tau).}
  The returned estimate is the sample mean of these contributions; its
  reported standard error is the sample standard deviation of the
  multivariate contribution divided by \eqn{\sqrt{n}}. The implementation
  accumulates this variance in centered form to avoid cancellation near zero
  or one. This beta CDF is observation-centered and is not the target-variable
  integral of Chen's evaluation-centered beta density kernel. See
  \code{\link{np.kernels}} for the design rationale and references.

  Fixed and adaptive order-2 beta CDF estimates are proper CDFs on the
  declared support. A generalized-NN order-2 curve retains exact endpoints
  and contribution-wise range but may be locally nonmonotone because its
  bandwidth changes with the target. Orders 4, 6, and 8 are signed and may
  overshoot \eqn{[0,1]} or be locally nonmonotone in any bandwidth mode.
  Returned values are not silently clipped or monotonized.

  Argument groups for bandwidth selection are documented on
  \code{\link{npudistbw}}. The most common workflow is to choose data
  and bandwidth inputs first, then bandwidth criterion and
  representation, then kernel/support controls, distribution-specific
  integral/grid controls, and numerical search controls if defaults
  need to be changed.

  For S3 plotting help, see \code{\link{plot.np}}. You can list
  available plot methods with \code{methods("plot")}.


  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    Usage 1: first compute the bandwidth object via npudistbw and then
    compute the cumulative distribution:

    ## Start npRmpi for interactive execution. If slaves are already running and
    ## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
    ## reuse the existing pool instead of respawning. To change the number of
    ## slaves, call `npRmpi.quit(force=TRUE)` then restart.
    npRmpi.init(nslaves=1)
    bw <- npudistbw(formula = ~y, data = mydat)
    Fhat <- npudist(bw)
    npRmpi.quit()

    Usage 2: alternatively, compute the bandwidth object indirectly:

    ## Start npRmpi for interactive execution. If slaves are already running and
    ## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
    ## reuse the existing pool instead of respawning. To change the number of
    ## slaves, call `npRmpi.quit(force=TRUE)` then restart.
    npRmpi.init(nslaves=1)
    Fhat <- npudist(formula = ~y, data = mydat)
    npRmpi.quit()

    Usage 3: modify the default kernel and order:

    ## Start npRmpi for interactive execution. If slaves are already running and
    ## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
    ## reuse the existing pool instead of respawning. To change the number of
    ## slaves, call `npRmpi.quit(force=TRUE)` then restart.
    npRmpi.init(nslaves=1)
    Fhat <- npudist(formula = ~y, data = mydat, ckertype="epanechnikov", ckerorder=4)
    npRmpi.quit()

    Usage 4: use the data frame interface rather than the formula
    interface:

    ## Start npRmpi for interactive execution. If slaves are already running and
    ## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
    ## reuse the existing pool instead of respawning. To change the number of
    ## slaves, call `npRmpi.quit(force=TRUE)` then restart.
    npRmpi.init(nslaves=1)
    Fhat <- npudist(tdat = y, ckertype="epanechnikov", ckerorder=4)
    npRmpi.quit()
  }

  \code{npudist} implements a variety of methods for estimating
  multivariate cumulative distributions (\eqn{p}-variate) defined over a
  set of possibly continuous and/or discrete (ordered) data. The
  approach is based on Li and Racine (2003) who employ
  \sQuote{generalized product kernels} that admit a mix of continuous
  and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating
  the cumulative distribution at the point \eqn{x}. Generalized nearest-neighbor
  bandwidths change with the point at which the cumulative distribution is estimated,
  \eqn{x}. Fixed bandwidths are constant over the support of \eqn{x}.

  Data contained in the data frame \code{tdat} (and also \code{edat})
  may be a mix of continuous (default) and ordered discrete (to be
  specified in the data frame \code{tdat} using the
  \code{\link{ordered}} command). Data can be entered in an arbitrary
  order and data types will be detected automatically by the routine
  (see \code{\link{npRmpi}} for details).

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and
  eighth-order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Ordered data types use a variation of the Wang and van Ryzin
  (1981) kernel.

}
\value{

  \code{\link{npudist}} returns a \code{npdistribution} object. The
  generic accessor functions \code{\link{fitted}} and \code{\link{se}}
  extract estimated values and asymptotic standard errors on estimates,
  respectively, from the returned object. Furthermore, the functions
  \code{\link{predict}}, \code{\link{summary}} and \code{\link{plot}}
  support objects of both classes. The returned objects have the
  following components:

  \item{eval}{ the evaluation points. }
  \item{dist}{ estimate of the cumulative distribution at the evaluation points }
  \item{derr}{ standard errors of the cumulative distribution estimates }

}
\section{Book And Method Pointers}{
  The unconditional distribution target is \eqn{F(x)=\Pr(X\le x)}.
  The estimator replaces the density kernel factors by integrated
  distribution-kernel factors where appropriate, while retaining the
  same mixed-data product-kernel logic for continuous and ordered
  components. The returned standard errors are asymptotic standard
  errors for the estimated cumulative distribution values.

  For book-length derivations, see Li and Racine (2007), Chapter 1
  \emph{Density Estimation}, especially Sections 1.4 and 1.5 for CDF
  estimation and bandwidth selection, together with Chapter 4
  \emph{Kernel Estimation with Mixed Data}. The later workflow
  treatment is Racine (2019), Chapter 2 \emph{Continuous Density and
  Cumulative Distribution Functions} and Chapter 3 \emph{Mixed-Data
  Probability Density and Cumulative Distribution Functions}.
}

\references{
  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{ Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2003), \dQuote{Nonparametric estimation of
    distributions with categorical and continuous data,} Journal
  of Multivariate Analysis, 86, 266-292.

  Li, C. and H. Li and J.S. Racine (2017), \dQuote{Cross-Validated Mixed
  Datatype Bandwidth Selection for Nonparametric Cumulative
  Distribution/Survivor Functions,} Econometric Reviews, \bold{36},
  970-987.

  Chen, S.X. (1999), \dQuote{Beta kernel estimators for density functions,}
  Computational Statistics & Data Analysis, 31(2), 131--145.
  doi:10.1016/S0167-9473(99)00010-9.

  Igarashi, G. and Kakizawa, Y. (2020), \dQuote{Higher-order bias corrections
  for asymmetric kernel density estimators,} Journal of Nonparametric
  Statistics, 32(3), 716--752. doi:10.1080/10485252.2020.1770754.

  Mansouri, B., Rastin, A., and Allah Mombeni, H. (2024),
  \dQuote{Beta Kernel Estimator for a Cumulative Distribution Function with
  Bounded Support,} Journal of Sciences, Islamic Republic of Iran,
  \bold{34}(4), 349--361 (volume dated 2023; accepted and DOI registered in
  2024).
  doi:10.22059/jsciences.2024.367933.1007836.

  Ouyang, D. and Q. Li and J.S. Racine (2006), \dQuote{Cross-validation
    and the estimation of probability distributions with categorical
    data,} Journal of Nonparametric Statistics, 18, 69-100.

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

  When reusing a bandwidth object without explicit training or evaluation
  data, the required data are recovered from object-retained state or
  from the stored call environment when available. This is a convenience
  for ordinary same-session use, not a guarantee that a saved object is
  self-contained. For reproducible replay after \code{saveRDS()} and
  \code{readRDS()} in a fresh R session, supply the required data
  arguments explicitly.
}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}, \code{\link{plot.np}} \code{\link{npudistbw}} , \code{\link{density}} }

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

bw <- npudistbw(formula=~ordered(year)+gdp,
                data=Italy)

Fhat <- npudist(bws=bw)

summary(Fhat)

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
\keyword{ nonparametric }% at least one, from doc/KEYWORDS
