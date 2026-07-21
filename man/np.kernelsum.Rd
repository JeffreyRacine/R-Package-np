% 
\name{npksum}
\alias{npksum}
\alias{npksum.default}
\alias{npksum.formula}
\alias{npksum.numeric}
\title{ Kernel Sums with Mixed Data Types }
\description{
  \code{npksum} computes kernel sums on evaluation
  data, given a set of training data, data to be weighted (optional), and a
  bandwidth specification (any bandwidth object).

}
\usage{
npksum(\dots)

\method{npksum}{formula}(formula, 
       data, 
       newdata, 
       subset, 
       na.action, 
       \dots)

\method{npksum}{default}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat = NULL,
       exdat = NULL,
       weights = NULL,
       bandwidth.divide = FALSE,
       compute.ocg = FALSE,
       compute.score = FALSE,
       kernel.pow = 1.0,
       leave.one.out = FALSE,
       operator = names(ALL_OPERATORS),
       permutation.operator = names(PERMUTATION_OPERATORS),
       return.kernel.weights = FALSE,
       \dots)

\method{npksum}{numeric}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat,
       exdat,
       weights,
       bandwidth.divide,
       compute.ocg,
       compute.score,
       kernel.pow,
       leave.one.out,
       operator,
       permutation.operator,
       return.kernel.weights,
       \dots)
}

\arguments{
\subsection{Bandwidth And Data Inputs}{
Core inputs defining the training data, optional weighted response, and
evaluation points for the generalized product kernel sum.
}

  \item{bws}{
    a bandwidth specification. This can be set as any suitable bandwidth
    object returned from a bandwidth-generating function, or a numeric vector.
  }
  
  \item{txdat}{
    a \eqn{p}-variate data frame of sample realizations (training data) used to
    compute the sum.
  }
  
  \item{tydat}{
    a numeric vector of data to be weighted. The \eqn{i}th kernel weight
    is applied to the \eqn{i}th element. Defaults to  \code{1}.
  }
  
  \item{exdat}{
    a \eqn{p}-variate data frame of sum evaluation points (if omitted,
    defaults to the training data itself).
  }
  
  \item{weights}{
    a \eqn{n} by \eqn{q}  matrix of weights which can optionally be
    applied to \code{tydat} in the sum. See details.
  }
  
\subsection{Kernel-Sum Operators And Output}{
Controls for bandwidth normalization, kernel powers, leave-one-out
evaluation, operator variants, and returned kernel weights.
}

  \item{bandwidth.divide}{
    a single, non-missing logical specifying whether or not to divide
    continuous kernel weights by their bandwidths. Use this with
    nearest-neighbor methods. Defaults to \code{FALSE}.
  }

  \item{compute.ocg}{
    a single, non-missing logical specifying whether or not to return a
    separate result for each unordered and ordered dimension, where the
    product kernel term for that dimension is evaluated at an appropriate
    reference category. This is used primarily in \code{npRmpi} to compute
    ordered and unordered categorical gradients. See details.
  }

  \item{compute.score}{
    a single, non-missing logical specifying whether or not to return the
    score (the \sQuote{grad h} terms) for each dimension in addition to the
    kernel sum. Cannot be \code{TRUE} if a permutation operator other
    than \code{"none"} is selected.
  }

  \item{kernel.pow}{
    a single, finite integer specifying the power to which the kernels
    will be raised in the sum. Defaults to \code{1}.
  }

  \item{leave.one.out}{
    a single, non-missing logical value specifying whether or not to
    compute the leave one out sums. Will not work if \code{exdat} is
    specified. Defaults to \code{FALSE}.
  }

  \item{operator}{
    a string specifying whether the \code{normal}, \code{convolution},
    \code{derivative}, or \code{integral} kernels are to be
    used. Operators scale results by factors of \eqn{h} or \eqn{1/h} where
    appropriate. Defaults to \code{normal} and applies to all elements in a
    multivariate object. See details.
  }

  \item{permutation.operator}{
    a string which can have a value of \code{none}, \code{normal},
    \code{derivative}, or \code{integral}. If set to something other
    than \code{none} (the default), then a separate result will be returned for each
    term in the product kernel, with the operator applied to that term.
    Permutation operators scale results by factors of \eqn{h} or \eqn{1/h} where
    appropriate. This is useful for computing gradients, for example.
  }

  \item{return.kernel.weights}{
    a single, non-missing logical specifying whether or not to return the
    matrix of generalized product kernel weights. Defaults to
    \code{FALSE}. See details.
    }

\subsection{Formula Interface}{
Formula-method arguments for symbolic kernel-sum specifications.
}

  \item{formula}{
    a symbolic description of variables on which the sum is
    to be performed. The details of constructing a formula are
    described below.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which the
    function is called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, \code{data} is used.
  }

  \item{subset}{
    an optional vector specifying a subset of observations to be used. 
  }

  \item{na.action}{
    a function which indicates what should happen when the data contain
   \code{NA}s. The default is set by the \code{\link{na.action}} setting of options, and is
    \code{\link{na.fail}} if that is unset. The (recommended) default is
    \code{\link{na.omit}}.  
  }


\subsection{Additional Arguments}{
Further arguments passed to the default kernel-sum method.
}

  \item{\dots}{
    additional arguments supplied to specify the parameters to the
    \code{default} S3 method, which is called during estimation. 
  }
}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, \code{\link{plot}} for plotting options, and \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See \code{\link{npRmpi.init}} details for performance tradeoffs (message passing/startup mode) and the \code{inst/Rprofile} manual-broadcast template.

  \code{npksum}
  exists so that you can create your own kernel objects with
  or without a variable to be weighted (default \eqn{Y=1}). With the options
  available, you could create new nonparametric tests or even new kernel
  estimators. The convolution kernel option would allow you to create,
  say, the least squares cross-validation function for kernel density
  estimation.
  
  \code{npksum} uses highly-optimized C code that strives to minimize
  its \sQuote{memory footprint}, while there is low overhead involved
  when using repeated calls to this function (see, by way of
  illustration, the example below that conducts leave-one-out
  cross-validation for a local constant regression estimator via calls
  to the \code{R} function \code{\link{nlm}}, and compares this to the
  \code{\link{npregbw}} function).

  \code{npksum} implements a variety of methods for computing
  multivariate kernel sums (\eqn{p}-variate) defined over a set of
  possibly continuous and/or discrete (unordered, ordered) data. The
  approach is based on Li and Racine (2003) who employ
  \sQuote{generalized product kernels} that admit a mix of continuous
  and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating
  the kernel sum at the point \eqn{x}. Generalized nearest-neighbor
  bandwidths change with the point at which the sum is computed,
  \eqn{x}. Fixed bandwidths are constant over the support of \eqn{x}.

  \code{npksum} computes \eqn{\sum_{j=1}^{n}{W_j^\prime Y_j
  K(X_j)}}{sum(t(W[j]) \%*\% Y[j] * K(X[j]))}, where \eqn{W_j}{W[j]}
  represents a row vector extracted from \eqn{W}.  That is, it computes
  the kernel weighted sum of the outer product of the rows of \eqn{W}
  and \eqn{Y}. In the examples, the uses of such sums are illustrated.

  \code{npksum} may be invoked \emph{either} with a formula-like
  symbolic 
  description of variables on which the sum is to be
  performed \emph{or} through a simpler interface whereby data is passed
  directly to the function via the \code{txdat} and \code{tydat}
  parameters. Use of these two interfaces is \bold{mutually exclusive}.

  Data contained in the data frame \code{txdat} (and also \code{exdat})
  may be a mix of continuous (default), unordered discrete (to be
  specified in the data frame \code{txdat} using the
  \code{\link{factor}} command), and ordered discrete (to be specified
  in the data frame \code{txdat} using the \code{\link{ordered}}
  command). Data can be entered in an arbitrary order and data types
  will be detected automatically by the routine (see \code{\link{npRmpi}}
  for details).

  Data for which bandwidths are to be estimated may be specified
  symbolically. A typical description has the form \code{dependent data
  ~ explanatory data}, where \code{dependent data} and \code{explanatory
  data} are both series of variables specified by name, separated by the
  separation character '+'. For example, \code{ y1 ~ x1 + x2 } specifies
  that \code{y1} is to be kernel-weighted by \code{x1} and \code{x2}
  throughout the sum. See below for further examples.

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and
  eighth order Gaussian, Epanechnikov, and beta associated kernels, and the
  uniform kernel. The beta family uses finite-support coordinate-aware routes
  described below. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of
  the Wang and van Ryzin (1981) kernel (see \code{\link{npRmpi}} for
  details).

  For \code{ckertype="beta"}, \code{npksum} uses Chen's second-order
  associated PDF kernel as the base contribution on each declared interval.
  If \eqn{L=b-a}, its per-observation contribution is
  \deqn{\frac{1}{L}\,\mathrm{BetaPDF}\left(
    \frac{x_i-a}{L};
    1+\frac{x-a}{L}\left(\frac{L}{h}\right)^2,
    1+\frac{b-x}{L}\left(\frac{L}{h}\right)^2\right).}
  The complete density scaling \eqn{1/L} is part of this contribution.
  Consequently, a univariate beta density estimate is
  \code{npksum(...)$ksum/n}, not \code{npksum(...)$ksum/(n*h)}.
  Local-constant regression uses the ordinary ratio of beta-weighted response
  and unit sums. Do not set \code{bandwidth.divide=TRUE}; that legacy
  normalization is intentionally rejected for beta kernels.

  For order \eqn{2m}, the PDF and CDF operators use
  \eqn{\sum_{s=1}^m(-1)^{s+1}{m\choose s}K_2(h\sqrt{s})}. Thus the coefficient
  vectors for orders 2, 4, 6, and 8 are \eqn{(1)}, \eqn{(2,-1)},
  \eqn{(3,-3,1)}, and \eqn{(4,-6,4,-1)}. Higher-order PDF and CDF
  contributions are signed and are returned without clipping.

  For \code{operator="integral"}, the univariate beta contribution is
  \deqn{I_{(x-a)/L}\left(1+\frac{x_i-a}{L}\left(\frac{L}{h}\right)^2,
    1+\frac{b-x_i}{L}\left(\frac{L}{h}\right)^2\right).}
  This observation-centered analytic CDF is not the integral of the
  evaluation-centered PDF contribution with respect to \eqn{x}. It is exactly
  zero and one at the declared endpoints, is monotone for order 2, and has the
  uniform-CDF limit as \eqn{h} diverges. In several dimensions,
  \code{operator} may mix \code{"normal"}, \code{"integral"}, and
  \code{"convolution"} beta factors coordinate by coordinate.

  For \code{operator="convolution"}, the training and evaluation values are
  interpreted as the centers of two associated beta weights in the common
  observation coordinate. If their shapes are
  \eqn{(\alpha_1,\beta_1)} and \eqn{(\alpha_2,\beta_2)}, the returned factor is
  \deqn{\frac{1}{L}\frac{B(\alpha_1+\alpha_2-1,
                                  \beta_1+\beta_2-1)}
                             {B(\alpha_1,\beta_1)B(\alpha_2,\beta_2)}.}
  This symmetric analytic \eqn{L^2} overlap includes its complete
  original-scale \eqn{1/L} normalization. Its diagonal is the local beta
  roughness. It is neither a translation convolution nor the
  target-coordinate integrated-square term used by bounded density
  \code{cv.ls}; those objectives retain their existing quadrature of the
  actual estimator. Beta PDF, CDF, and overlap factors may be mixed across
  dimensions.

  This beta surface requires continuous \code{txdat}, \code{ckerorder} in
  \code{c(2,4,6,8)}, \code{ckerbound="fixed"} or \code{"range"}, finite
  lower and upper bounds for every dimension, \code{operator="normal"},
  \code{operator="integral"}, and/or \code{operator="convolution"}, and
  \code{kernel.pow=1}. It supports multiple continuous dimensions, response
  and case-weight columns, leave-one-out sums, returned kernel weights, and
  fixed/generalized-NN/adaptive-NN bandwidths. Generalized overlap evaluates
  generalized-NN distances at both centers. Adaptive overlap keeps adaptive
  bandwidths attached to training centers; a genuinely new evaluation center
  uses the established generalized-NN auxiliary bandwidth. Fixed and adaptive
  order-2 CDF factors are monotone in their target. Generalized-NN order-2 CDF
  factors retain exact endpoints and contribution-wise range but can be
  locally nonmonotone. Orders 4, 6, and 8 may additionally leave \eqn{[0,1]}
  or be locally nonmonotone in every bandwidth mode; no clipping or monotone
  repair is applied. One-coordinate analytic first derivatives and their
  continuous-coordinate permutations are available as documented below.
  Higher and mixed derivatives, categorical-product, score, OCG, and other
  unsupported beta operations fail explicitly. See \code{\link{np.kernels}}
  for the bandwidth convention, orientation, and references.

  The option \code{operator=} can be used to \sQuote{mix and match}
  operator strings to create a \sQuote{hybrid} kernel provided they
  match the dimension of the data. For example, for a two-dimensional
  data frame of \link{numeric} datatypes,
  \code{operator=c("normal","derivative")} will use the normal
  (i.e. PDF) kernel for variable one and the derivative of the PDF
  kernel for variable two. Please note that applying operators will scale the
  results by factors of \eqn{h} or \eqn{1/h} where appropriate.

  The option \code{permutation.operator=} can be used to \sQuote{mix and match}
  operator strings to create a \sQuote{hybrid} kernel, in addition to
  the kernel sum with no operators applied, one for each continuous
  dimension in the data. For example, for a two-dimensional
  data frame of \link{numeric} datatypes,
  \code{permutation.operator=c("derivative")} will return the usual
  kernel sum as if \code{operator = c("normal","normal")} in the
  \code{ksum} member, and in the \code{p.ksum} member, it will return
  kernel sums for \code{operator = c("derivative","normal")}, and
  \code{operator = c("normal","derivative")}. This makes the computation
  of gradients much easier.

  For \code{ckertype="beta"}, \code{operator="derivative"} is the analytic
  first derivative of the associated PDF contribution with respect to its
  evaluation coordinate. A direct operator vector may contain one derivative
  coordinate, while \code{permutation.operator="derivative"} returns one
  derivative product per continuous coordinate. Set
  \code{bandwidth.divide=FALSE}; beta contributions already carry their
  original-scale units. Generalized- and adaptive-NN bandwidths are held fixed
  during differentiation. Exact matching support endpoints can produce signed
  infinite scalar derivatives; see \code{\link{np.kernels}} for the endpoint
  and ratio-cancellation conventions. Signed infinities are returned as
  \code{-Inf}/\code{Inf}, not recoded as missing values; \code{NA_real_} is
  reserved for unavailable or indeterminate results.

  The option \code{compute.score=} can be used to compute the gradients
  with respect to \eqn{h} in addition to the normal kernel sum. Like
  permutations, the additional results are returned in the
  \code{p.ksum}. This option does not work in conjunction with
  \code{permutation.operator}.

  The option \code{compute.ocg=} works much like \code{permutation.operator},
  but for discrete variables. The kernel is evaluated at a reference
  category in each dimension: for ordered data, the next lowest category
  is selected, except in the case of the lowest category, where the
  second lowest category is selected; for unordered data, the first
  category is selected. These additional data are returned in the
  \code{p.ksum} member. This option can be set simultaneously with
  \code{permutation.operator}.

  When a permutation, score, or OCG operation selects only a subset of
  variables, the corresponding blocks in \code{p.ksum} and, when requested
  internally, derivative kernel weights are packed without empty coordinates
  and follow the original data-column order among the selected variables.
  
  The option \code{return.kernel.weights=TRUE} returns a matrix of
  dimension \sQuote{number of training observations} by \sQuote{number
  of evaluation observations} and contains only the generalized product
  kernel weights ignoring all other objects and options that may be
  provided to \code{npksum} (e.g. \code{bandwidth.divide=TRUE} will be
  ignored, etc.). Summing the columns of the weight matrix and dividing
  by \sQuote{number of training observations} times the product of the
  bandwidths (i.e. \code{\link{colMeans}(foo$kw)/prod(h)}) would produce
  the kernel estimator of a (multivariate) density
  (\code{operator="normal"}) or multivariate cumulative distribution
  (\code{operator="integral"}).

}
\value{
  \code{npksum} returns a \code{npkernelsum} object
  with the following components:

  \item{eval}{ the evaluation points }
  \item{ksum}{ the sum at the evaluation points }
  \item{kw}{ the kernel weights (when \code{return.kernel.weights=TRUE}
  is specified)}
  \item{p.ksum}{ packed permutation, score, or OCG sums, when requested }
}
\references{
  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{ Multivariate binary
  discrimination by the kernel method,} Biometrika, 63, 413-420.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2003), \dQuote{Nonparametric estimation of
  distributions with categorical and continuous data,} Journal
  of Multivariate Analysis, 86, 266-292.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,} Cambridge
  University Press. 

  Scott, D.W. (1992), \emph{Multivariate Density Estimation. Theory,
  Practice and Visualization,} New York: Wiley.

  Silverman, B.W. (1986), \emph{Density Estimation,} London: Chapman and
  Hall.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
  estimators for discrete distributions,} Biometrika, 68, 301-309.  }

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

\seealso{\code{\link{npRmpi.init}} for MPI startup and workflow guidance.}
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

n <- 100000
x <- rnorm(n)
x.eval <- seq(-4, 4, length=50)

bw <- npudensbw(dat=x, bwmethod="normal-reference")

den.ksum <- npksum(txdat=x, exdat=x.eval, bws=bw$bw,
                   bandwidth.divide=TRUE)$ksum/n

den.ksum

## A beta associated kernel includes its original-scale normalization.
## The target/centre is supplied in exdat and the observation-coordinate
## trace in txdat; return.kernel.weights exposes the individual curve.

z <- seq(0, 1, length.out=101)
beta.trace <- npksum(txdat=z, exdat=0.08, bws=0.12,
                    ckertype="beta", ckerorder=2,
                    ckerbound="fixed", ckerlb=0, ckerub=1,
                    return.kernel.weights=TRUE)$kw[, 1]
plot(z, beta.trace, type="l", xlab="z", ylab="Beta kernel")

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
\keyword{ nonparametric }
