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
       operator = c("normal","convolution","derivative","integral"),
       smooth.coefficient = FALSE,
       \dots)

\method{npksum}{numeric}(bws,
       txdat = stop("training data 'txdat' missing"),
       tydat,
       exdat,
       weights,
       leave.one.out, kernel.pow, bandwidth.divide,
       operator, smooth.coefficient,
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
    used. Defaults to \code{normal}.
  }

  \item{smooth.coefficient}{
    a logical specifying whether or not to use certain optimisations if
    the smooth coefficient estimator is being computed. Currently does
    nothing. Defaults to \code{FALSE}.
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
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel (see \code{\link{np}} for details).

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
  Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
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

plot(x.eval, dnorm(x.eval), xlab="X", ylab="Density", type="l")
points(x.eval, den.ksum, col="blue")
points(x.eval, den, col="red", cex=0.2)
legend(1, .4, 
       c("DGP", "npksum()", 
       "npudens()"), 
       col=c("black", "blue", "red"), 
       lty=c(1, 1, 1))

\dontrun{

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 2: For this example, we first obtain residuals from a parametric
# regression model, then compute a vector of leave-one-out kernel
# weighted sums of squared residuals where the kernel function is raised
# to the power 2. We use a normal-reference rule-of-thumb obtained by
# using the scaling factor 1.06 and bwscaling=TRUE. Note that this
# returns a vector of kernel weighted sums, one for each element of the
# error vector. Note also that you can specify the bandwidth type, 
# kernel function, kernel order etc.

data("cps71")
attach(cps71)

error <- residuals(lm(logwage~age+I(age^2)))

ksum <- npksum(txdat=age, 
               tydat=error^2, 
               bws=1.06, 
               leave.one.out=TRUE, 
               bwscaling=TRUE, 
               kernel.pow=2)

ksum

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Note - npksum also can inherit properties from a bandwidth
# object. To again specify a value for the bandwidth, you might do so in
# the bandwidth object itself then pass this into npksum (e.g., 
# below I use a regression bandwidth object but you could use, say, a
# density bandwidth object instead).

bw <- npregbw(xdat=age, ydat=logwage, bws=c(1.06),
              bwscaling=TRUE, 
              bandwidth.compute=FALSE)

ksum <- npksum(txdat=age, 
               tydat=error^2, 
               bws=bw$bw, 
               leave.one.out=TRUE, 
               kernel.pow=2)

# Obviously, if we wanted the sum of these weighted kernel sums then, 
# trivially, 

sum(ksum$ksum)

detach(cps71)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

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
attach(cps71)

bw <- npregbw(xdat=age, ydat=logwage)

fit.lc <- npksum(txdat=age, tydat=logwage, bws=bw$bw)$ksum/
          npksum(txdat=age, bws=bw$bw)$ksum

plot(age, logwage, xlab="Age", ylab="log(wage)")
lines(age, fit.lc)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# If you wished to compute the kernel sum for a set of evaluation points, 
# you first generate the set of points then feed this to npksum, 
# e.g., for the set (20, 30, ..., 60) we would use

age.seq <- seq(20, 60, 10)

fit.lc <- npksum(txdat=age, exdat=age.seq, tydat=logwage, bws=bw$bw)$ksum/
          npksum(txdat=age, exdat=age.seq, bws=bw$bw)$ksum

# Note that now fit.lc contains the 5 values of the local constant
# estimator corresponding to age.seq...

fit.lc

detach(cps71)

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
# objective function minima)

options(scipen=100, digits=12)

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

Sys.sleep(5)

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
} % enddontrun
}
\keyword{ nonparametric }
