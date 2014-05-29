% $Id: np.kernelsum.Rd,v 1.64 2006/11/03 21:17:20 tristen Exp $
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

\method{npksum}{formula}(formula, data, newdata, subset, na.action, \dots)

\method{npksum}{default}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat = NULL,
       exdat = NULL,
       weights = NULL,
       leave.one.out = FALSE,
       kernel.pow = 1.0,
       bandwidth.divide = FALSE,
       operator = names(ALL_OPERATORS),
       permutation.operator = names(PERMUTATION_OPERATORS),
       compute.score = FALSE,
       compute.ocg = FALSE,
       return.kernel.weights = FALSE,
       \dots)

\method{npksum}{numeric}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat,
       exdat,
       weights,
       leave.one.out,
       kernel.pow,
       bandwidth.divide,
       operator,
       permutation.operator,
       compute.score,
       compute.ocg,
       return.kernel.weights,
       \dots)
}

\arguments{
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

  \item{\dots}{
    additional arguments supplied to specify the parameters to the
    \code{default} S3 method, which is called during estimation. 
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
  
  \item{bws}{
    a bandwidth specification. This can be set as any suitable bandwidth
    object returned from a bandwidth-generating function, or a numeric vector.
  }
  
  \item{weights}{
    a \eqn{n} by \eqn{q}  matrix of weights which can optionally be
    applied to \code{tydat} in the sum. See details.
  }
  
  \item{leave.one.out}{
    a logical value to specify whether or not to compute the leave one
    out sums. Will not work if \code{exdat} is specified. Defaults to
    \code{FALSE}.
  }

  \item{kernel.pow}{
    an integer specifying the power to which the kernels will be raised
    in the sum. Defaults to \code{1}.
  }

  \item{bandwidth.divide}{
    a logical specifying whether or not to divide continuous kernel
    weights by their bandwidths. Use this with nearest-neighbor
    methods. Defaults to \code{FALSE}.
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

  \item{compute.score}{
    a logical specifying whether or not to return the score
    (the \sQuote{grad h} terms) for each dimension in addition to the kernel
    sum. Cannot be \code{TRUE} if a permutation operator other than
    \code{"none"} is selected. 
  }

  \item{compute.ocg}{
    a logical specifying whether or not to return a separate
    result for each unordered and ordered dimension, where the product
    kernel term for that dimension is evaluation at an appropriate
    reference category. This is used primarily in \code{np} to compute ordered and
    unordered categorical gradients. See details.
  }

  \item{return.kernel.weights}{
    a logical specifying whether or not to return the matrix of
    generalized product kernel weights. Defaults to \code{FALSE}. See
    details.
    }
}
\details{
  \code{npksum}
  exists so that you can create your own kernel objects with
  or without a variable to be weighted (default \eqn{Y=1}). With the options
  available, you could create new nonparametric tests or even new kernel
  estimators. The convolution kernel option would allow you to create,
  say, the least squares cross-validation function for kernel density
  estimation.
  
  \code{npksum} uses highly-optimized C code that strives
  to minimize its \sQuote{memory footprint}, while there is low overhead
  involved when using repeated calls to this function (see, by way of
  illustration, the example below that conducts leave-one-out
  cross-validation for a local constant regression estimator via calls
  to the \code{R} function \code{\link{nlm}}, and compares this to the
  \code{\link{npregbw}} function).

  \code{npksum} implements a variety
  of methods for computing
  multivariate kernel sums (\eqn{p}-variate) defined over a set of possibly
  continuous and/or discrete (unordered, ordered) data. The approach is
  based on Li and Racine (2003) who employ \sQuote{generalized product kernels}
  that admit a mix of continuous and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
  kernel sum at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
  with the point at which the sum is computed, \eqn{x}. Fixed bandwidths
  are constant over the support of \eqn{x}.

  \code{npksum} computes
  \eqn{\sum_{j=1}^{n}{W_j^\prime Y_j K(X_j)}}{sum(t(W[j]) \%*\% Y[j] *
  K(X[j]))}, where \eqn{A_j}{A[j]} represents a row vector extracted
  from \eqn{A}.
  That is, it computes the kernel 
  weighted sum of the outer product of the rows of \eqn{W} and \eqn{Y}. In the
  examples, the uses of such sums are illustrated.

  \code{npksum} may be invoked \emph{either} with a formula-like
  symbolic 
  description of variables on which the sum is to be
  performed \emph{or} through a simpler interface whereby data is passed
  directly to the function via the \code{txdat} and \code{tydat}
  parameters. Use of these two interfaces is \bold{mutually exclusive}.

  Data contained in the data frame \code{txdat} (and also \code{exdat})
  may be a mix of continuous (default), unordered discrete (to be
  specified in the data frame \code{txdat} using the \code{\link{factor}} command),
  and ordered discrete (to be specified in the data frame \code{txdat}
  using the \code{\link{ordered}} command). Data can be entered in an arbitrary
  order and data types will be detected automatically by the routine
  (see \code{\link{np}} for details).

  Data for which bandwidths are to be estimated may be specified
  symbolically. A typical description has the form \code{dependent data
  ~ explanatory data},
  where \code{dependent data} and \code{explanatory data} are both
  series of variables specified by name, separated by 
  the separation character '+'. For example, \code{ y1 ~ x1 + x2 }
  specifies that \code{y1} is to be kernel-weighted by \code{x1} and 
  \code{x2} throughout the sum. See below for further examples. 

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and

  Aitken's (1976) kernel, while ordered data types use a variation of
  the Wang and van Ryzin (1981) kernel (see \code{\link{np}} for
  details).

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

  The option \code{compute.score=} can be used compute the gradients
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

\examples{
\dontrun{
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

mpi.spawn.Rslaves(nslaves=1)
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

n <- 100000
x <- rnorm(n)
x.eval <- seq(-4, 4, length=50)

mpi.bcast.Robj2slave(x)
mpi.bcast.Robj2slave(x.eval)

mpi.bcast.cmd(bw <- npudensbw(dat=x, bwmethod="normal-reference"),
              caller.execute=TRUE)

mpi.bcast.cmd(den.ksum <- npksum(txdat=x, exdat=x.eval, bws=bw$bw,
                   bandwidth.divide=TRUE)$ksum/n,
              caller.execute=TRUE)

den.ksum

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

mpi.close.Rslaves()

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun     
}
\keyword{ nonparametric }
