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
\subsection{Data, Bandwidth Inputs And Formula Interface}{
Core inputs defining the symbolic interface, training data, optional weighted
response, and evaluation points for the generalized product kernel sum.
}

  \item{bws}{
    a bandwidth specification. This can be set as any suitable bandwidth
    object returned from a bandwidth-generating function, or a numeric vector.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which the
    function is called.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of sum evaluation points (if omitted,
    defaults to the training data itself).
  }

  \item{formula}{
    a symbolic description of variables on which the sum is
    to be performed. The details of constructing a formula are
    described below.
  }

  \item{na.action}{
    a function which indicates what should happen when the data contain
   \code{NA}s. The default is set by the \code{\link{na.action}} setting of options, and is
    \code{\link{na.fail}} if that is unset. The (recommended) default is
    \code{\link{na.omit}}.
  }

  \item{newdata}{
    an optional data frame in which to look for evaluation data. If
    omitted, \code{data} is used.
  }

  \item{subset}{
    an optional vector specifying a subset of observations to be used.
  }
  
  \item{txdat}{
    a \eqn{p}-variate data frame of sample realizations (training data) used to
    compute the sum.
  }
  
  \item{tydat}{
    a numeric vector of data to be weighted. The \eqn{i}th kernel weight
    is applied to the \eqn{i}th element. Defaults to  \code{1}.
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
    a logical specifying whether or not to divide continuous kernel
    weights by their bandwidths. Use this with nearest-neighbor
    methods. Defaults to \code{FALSE}.
  }

  \item{compute.ocg}{
    a logical specifying whether or not to return a separate result for
    each unordered and ordered dimension, where the product kernel term
    for that dimension is evaluated at an appropriate reference
    category. This is used primarily in \code{np} to compute ordered and
    unordered categorical gradients. See details.
  }

  \item{compute.score}{
    a logical specifying whether or not to return the score
    (the \sQuote{grad h} terms) for each dimension in addition to the kernel
    sum. Cannot be \code{TRUE} if a permutation operator other than
    \code{"none"} is selected. 
  }

  \item{kernel.pow}{
    an integer specifying the power to which the kernels will be raised
    in the sum. Defaults to \code{1}.
  }

  \item{leave.one.out}{
    a logical value to specify whether or not to compute the leave one
    out sums. Will not work if \code{exdat} is specified. Defaults to
    \code{FALSE}.
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
    a logical specifying whether or not to return the matrix of
    generalized product kernel weights. Defaults to \code{FALSE}. See
    details.
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
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, and \code{\link{plot}} for plotting options.

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
  will be detected automatically by the routine (see \code{\link{np}}
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
  the Wang and van Ryzin (1981) kernel (see \code{\link{np}} for
  details).

  For \code{ckertype="beta"}, \code{npksum} uses Chen's second-order
  associated PDF kernel as the base contribution on each declared interval.
  If \eqn{L=b-a}, its
  per-observation contribution is
  \deqn{\frac{1}{L}\,\mathrm{BetaPDF}\left(
    \frac{x_i-a}{L};
    1+\frac{x-a}{L}\left(\frac{L}{h}\right)^2,
    1+\frac{b-x}{L}\left(\frac{L}{h}\right)^2\right).}
  The complete density scaling \eqn{1/L} is part of this contribution.
  Consequently, a univariate beta density estimate is
  \code{npksum(...)$ksum/n}, not \code{npksum(...)$ksum/(n*h)}.
  Local-constant regression uses the ordinary ratio of beta-weighted response
  and unit sums.  Do not set \code{bandwidth.divide=TRUE}; that legacy
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

  This beta surface requires continuous \code{txdat},
  \code{ckerorder} in \code{c(2,4,6,8)},
  \code{ckerbound="fixed"} or \code{"range"}, finite lower and upper bounds for every
  dimension, \code{operator="normal"}, \code{operator="integral"}, and/or
  \code{operator="convolution"}, and \code{kernel.pow=1}. It supports
  multiple continuous dimensions, response and case-weight columns,
  leave-one-out sums, and returned kernel weights. All three operators support
  \code{bwtype="fixed"}. The PDF, CDF, and convolution operators additionally support
  \code{bwtype="generalized_nn"} and \code{bwtype="adaptive_nn"}; their
  evaluation-specific and observation-specific metric bandwidths are the
  distances returned by the existing observation-support nearest-neighbor
  implementation, including its tie and extended-NN policy. Generalized
  overlap evaluates generalized-NN distances at both centers. Adaptive
  overlap keeps adaptive bandwidths attached to training centers; a genuinely
  new evaluation center uses the established generalized-NN auxiliary
  bandwidth. Fixed and adaptive order-2 CDF factors are monotone in their
  target. Generalized-NN CDF factors remain in \eqn{[0,1]} with exact
  endpoints, but their target-varying bandwidth can produce local
  nonmonotonicity. Orders 4, 6, and 8 may additionally leave \eqn{[0,1]} or
  be locally nonmonotone in every bandwidth mode; no clipping or monotone
  repair is applied. Derivative, categorical-product, score, OCG, and
  permutation beta operations are not available and fail with an explicit
  condition. See \code{\link{np.kernels}}
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
  and ratio-cancellation conventions.

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
  estimators for discrete distributions,} Biometrika, 68, 301-309.

  Chen, S.X. (1999), \dQuote{Beta kernel estimators for density functions,}
  Computational Statistics & Data Analysis, 31(2), 131-145.
  doi:10.1016/S0167-9473(99)00010-9.  }

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
# EXAMPLE 1: For this example, we generate 100,000 observations from a
# N(0, 1) distribution, then evaluate the kernel density on a grid of 50
# equally spaced points using the npksum() function, then compare
# results with the (identical) npudens() function output

n <- 100000
x <- rnorm(n)
x.eval <- seq(-4, 4, length=50)

# Compute the bandwidth with the normal-reference rule-of-thumb

bw <- npudensbw(dat=x, bwmethod="normal-reference")

# Compute the univariate kernel density estimate using the 100,000
# training data points evaluated on the 50 evaluation data points, 
# i.e., 1/nh times the sum of the kernel function evaluated at each of
# the 50 points

den.ksum <- npksum(txdat=x, exdat=x.eval, bws=bw$bw)$ksum/(n*bw$bw[1])

# Note that, alternatively (easier perhaps), you could use the
# bandwidth.divide=TRUE argument and drop the *bw$bw[1] term in the
# denominator, as in
# npksum(txdat=x, exdat=x.eval, bws=bw$bw, bandwidth.divide=TRUE)$ksum/n

# Compute the density directly with the npudens() function...

den <- fitted(npudens(tdat=x, edat=x.eval, bws=bw$bw))

# Plot the true DGP, the npksum()/(nh) estimate and (identical)
# npudens() estimate

if (interactive()) plot(x.eval, dnorm(x.eval), xlab="X", ylab="Density", type="l")
points(x.eval, den.ksum, col="blue")
points(x.eval, den, col="red", cex=0.2)
legend(1, .4, 
       c("DGP", "npksum()", 
       "npudens()"), 
       col=c("black", "blue", "red"), 
       lty=c(1, 1, 1))

# Sleep for 5 seconds so that we can examine the output...

if (interactive()) Sys.sleep(5)

# EXAMPLE 2: For this example, we first obtain residuals from a
# parametric regression model, then compute a vector of leave-one-out
# kernel weighted sums of squared residuals where the kernel function is
# raised to the power 2. Note that this returns a vector of kernel
# weighted sums, one for each element of the error vector. Note also
# that you can specify the bandwidth type, kernel function, kernel order
# etc.

data("cps71")
with(cps71, {

error <- residuals(lm(logwage~age+I(age^2)))

bw <- npreg(error~age)

ksum <- npksum(txdat=age, 
               tydat=error^2, 
               bws=bw$bw,
               leave.one.out=TRUE, 
               kernel.pow=2)

ksum

# Obviously, if we wanted the sum of these weighted kernel sums then, 
# trivially, 

sum(ksum$ksum)

})

# Sleep for 5 seconds so that we can examine the output...

if (interactive()) Sys.sleep(5)

# Note that weighted leave-one-out sums of squared residuals are used to
# construct consistent model specification tests. In fact, the
# npcmstest() routine in this package is constructed purely from calls
# to npksum(). You can type npcmstest to see the npcmstest()
# code and also examine some examples of the many uses of
# npksum().

# EXAMPLE 3: For this example, we conduct local-constant (i.e., 
# Nadaraya-Watson) kernel regression. We shall use cross-validated
# bandwidths using npregbw() by way of example. Note we extract
# the kernel sum from npksum() via the `$ksum' argument in both
# the numerator and denominator.

data("cps71")
with(cps71, {

bw <- npregbw(xdat=age, ydat=logwage)

fit.lc <- npksum(txdat=age, tydat=logwage, bws=bw$bw)$ksum/
          npksum(txdat=age, bws=bw$bw)$ksum

if (interactive()) plot(age, logwage, xlab="Age", ylab="log(wage)")
lines(age, fit.lc)

# Sleep for 5 seconds so that we can examine the output...

if (interactive()) Sys.sleep(5)

# If you wished to compute the kernel sum for a set of evaluation points, 
# you first generate the set of points then feed this to npksum, 
# e.g., for the set (20, 30, ..., 60) we would use

age.seq <- seq(20, 60, 10)

fit.lc <- npksum(txdat=age, exdat=age.seq, tydat=logwage, bws=bw$bw)$ksum/
          npksum(txdat=age, exdat=age.seq, bws=bw$bw)$ksum

# Note that now fit.lc contains the 5 values of the local constant
# estimator corresponding to age.seq...

fit.lc

})

# EXAMPLE 4: For this example, we conduct least-squares cross-validation
# for the local-constant regression estimator. We first write an R
# function `ss' that computes the leave-one-out sum of squares using the
# npksum() function, and then feed this function, along with
# random starting values for the bandwidth vector, to the nlm() routine
# in R (nlm = Non-Linear Minimization). Finally, we compare results with
# the function npregbw() that is written solely in C and calls
# a tightly coupled C-level search routine.  Note that one could make
# repeated calls to nlm() using different starting values for h (highly
# recommended in general).

# Increase the number of digits printed out by default in R and avoid
# using scientific notation for this example (we wish to compare
# objective function minima), then restore the user's options on exit.

old <- options(scipen=100, digits=12)
on.exit(options(old), add = TRUE)

# Generate 100 observations from a simple DGP where one explanatory
# variable is irrelevant.

n <- 100

x1 <- runif(n)
x2 <- rnorm(n)
x3 <- runif(n)

txdat <- data.frame(x1, x2, x3)

# Note - x3 is irrelevant

tydat <- x1 + sin(x2) + rnorm(n)

# Write an R function that returns the average leave-one-out sum of
# squared residuals for the local constant estimator based upon
# npksum(). This function accepts one argument and presumes that
# txdat and tydat have been defined already.

ss <- function(h) {

# Test for valid (non-negative) bandwidths - return infinite penalty
# when this occurs

  if(min(h)<=0) {

    return(.Machine$double.xmax)

  } else {

    mean <-  npksum(txdat, 
                    tydat, 
                    leave.one.out=TRUE, 
                    bandwidth.divide=TRUE, 
                    bws=h)$ksum/
             npksum(txdat, 
                    leave.one.out=TRUE, 
                    bandwidth.divide=TRUE, 
                    bws=h)$ksum
  
    return(sum((tydat-mean)^2)/length(tydat))

  }

}

# Now pass this function to R's nlm() routine along with random starting
# values and place results in `nlm.return'.

nlm.return <- nlm(ss, runif(length(txdat)))

bw <- npregbw(xdat=txdat, ydat=tydat)

# Bandwidths from nlm()

nlm.return$estimate

# Bandwidths from npregbw()

bw$bw

# Function value (minimum) from nlm()

nlm.return$minimum

# Function value (minimum) from npregbw()

bw$fval

# Sleep for 5 seconds so that we can examine the output...

if (interactive()) Sys.sleep(5)

# EXAMPLE 5: For this example, we use npksum() to plot the kernel
# function itself. Our `training data' is the singleton, 0, and our
# evaluation data a grid in [-4,4], while we use a bandwidth of 1. By
# way of example we plot a sixth order Gaussian kernel (note that this
# kernel function can assume negative values)

x <- 0
x.eval <- seq(-4,4,length=500)

kernel <- npksum(txdat=x,exdat=x.eval,bws=1,
                 bandwidth.divide=TRUE,
                 ckertype="gaussian",
                 ckerorder=6)$ksum

plot(x.eval,kernel,type="l",xlab="X",ylab="Kernel Function") 
abline(0,0)

# EXAMPLE 6: The beta associated kernel for bounded data.  Its complete
# original-scale normalization is already included, so divide the sum by n,
# not by n*h.  The same weights can be used directly in a local-constant
# numerator/denominator ratio.

x <- c(0, 0.01, 0.04, 0.15, 0.4, 0.75, 1)
x.eval <- seq(0, 1, length.out=101)
h <- 0.12

beta.args <- list(txdat=data.frame(x=x),
                  exdat=data.frame(x=x.eval),
                  bws=h,
                  ckertype="beta",
                  ckerorder=2,
                  ckerbound="fixed",
                  ckerlb=0,
                  ckerub=1)

beta.density <- do.call(npksum, beta.args)$ksum/length(x)

y <- sqrt(x)
beta.mean <- do.call(npksum, c(beta.args, list(tydat=y)))$ksum/
             do.call(npksum, beta.args)$ksum

beta.cdf <- do.call(npksum,
                    c(beta.args,
                      list(operator="integral")))$ksum/length(x)

beta.overlap <- do.call(npksum,
                        c(beta.args,
                          list(operator="convolution")))$ksum
} % enddontrun
}
\keyword{ nonparametric }
