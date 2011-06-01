% $ID$
\name{npplregbw}
\alias{npplregbw}
\alias{npplregbw.NULL}
\alias{npplregbw.default}
\alias{npplregbw.formula}
\alias{npplregbw.plbandwidth}

\title{Partially Linear Kernel Regression Bandwidth Selection with Mixed Data Types}

\description{
  \code{npplregbw} computes a bandwidth object for a partially linear kernel regression estimate of a one (1)
  dimensional dependent variable on
  \eqn{p+q}-variate explanatory data, using the model
  \eqn{Y = X\beta + \Theta (Z) + \epsilon}{Y = XB + theta(Z) + epsilon}
  given a set of estimation points, training
  points (consisting of explanatory data and dependent data), and a
  bandwidth specification, which can be a \code{rbandwidth} object, or a
  bandwidth vector, bandwidth type and kernel type. 
}

\usage{
npplregbw(\dots)

\method{npplregbw}{formula}(formula, data, subset, na.action, call, \dots)

\method{npplregbw}{NULL}(xdat = stop("invoked without data `xdat'"),
          ydat = stop("invoked without data `ydat'"),
          zdat = stop("invoked without data `zdat'"),
          bws,
          \dots)

\method{npplregbw}{default}(xdat = stop("invoked without data `xdat'"),
          ydat = stop("invoked without data `ydat'"),
          zdat = stop("invoked without data `zdat'"),
          bws,
          \dots,
          bandwidth.compute = TRUE,
          nmulti,
          remin,
          itmax,
          ftol,
          tol,
          small)

\method{npplregbw}{plbandwidth}(xdat = stop("invoked without data `xdat'"),
          ydat = stop("invoked without data `ydat'"),
          zdat = stop("invoked without data `zdat'"),
          bws,
          nmulti,
          \dots)

}
\arguments{
  \item{formula}{
    a symbolic description of variables on which bandwidth selection is
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

  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process. 
  }

  \item{na.action}{
    a function which indicates what should happen when the data contain
   \code{NA}s. The default is set by the \code{\link{na.action}} setting of options, and is
    \code{\link{na.fail}} if that is unset. The (recommended) default is
    \code{\link{na.omit}}.  
  }
  
  \item{call}{
    the original function call. This is passed internally by
    \code{\link{np}} when a bandwidth search has been implied by a call to
    another function. It is not recommended that the user set this.
  }

  \item{xdat}{
  a \eqn{p}-variate data frame of explanatory data (training data),
  corresponding to \eqn{X} in the model equation, whose linear
  relationship with the dependent data \eqn{Y} is posited.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
  \code{xdat}. 
  }

  \item{zdat}{
    a \eqn{q}-variate data frame of explanatory data (training data),
    corresponding to \eqn{Z} in the model equation, whose relationship
    to the dependent variable is unspecified (nonparametric)
  }

  \item{bws}{
    a bandwidth specification. This can be set as a \code{plbandwidth}
    object returned from an invocation of \code{npplregbw}, or as
    a matrix of bandwidths, each row is a set of bandwidths for \eqn{Z},
    with a column for each variable \eqn{Z_i}{Z[i]}. In the first row
    are the bandwidths for the regression of \eqn{Y} on \eqn{Z}, the
    following rows contain the bandwidths for the regressions of the
    columns of \eqn{X} on \eqn{Z}.  If specified as a matrix
    additional arguments will need to be supplied as necessary to
    specify the bandwidth type, kernel types, and so on.

    If left unspecified, \code{npplregbw} will search for optimal
    bandwidths using \code{\link{npregbw}} in the course of
    calculations. If specified, \code{npplregbw} will use the given
    bandwidths as the starting point for the numerical search for optimal
    bandwidths, unless you specify bandwidth.compute = FALSE. 
  }

  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, selection methods, and so on. To do
    this, you may specify any of \code{regtype}, 
    \code{bwmethod}, \code{bwscaling}, \code{bwtype}, \code{ckertype},
    \code{ckerorder}, \code{ukertype}, \code{okertype}, as described in
    \code{\link{npregbw}}.
  }
  
  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{plbandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
  }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points. Defaults to \code{min(5,ncol(zdat))}.
  }

  \item{remin}{
    a logical value which when set as \code{TRUE} the search routine
    restarts from located minima for a minor gain in accuracy. Defaults
    to \code{TRUE}
  }

  \item{itmax}{
    integer number of iterations before failure in the numerical
    optimization routine. Defaults to \code{10000}
  }

  \item{ftol}{
    tolerance on the value of the cross-validation function 
    evaluated at located minima. Defaults to \code{1.19e-07
      (FLT_EPSILON)}
  }

  \item{tol}{
    tolerance on the position of located minima of the 
    cross-validation function. Defaults to \code{1.49e-08
      (sqrt(DBL_EPSILON))}
  }

  \item{small}{ a small number, at about the precision of the data type
    used. Defaults to \code{2.22e-16 (DBL_EPSILON)}

  }

}

\details{
 \code{npplregbw} implements a variety of methods for nonparametric regression on
 multivariate (\eqn{q}-variate) explanatory data defined over a set of possibly
 continuous and/or discrete (unordered, ordered) data. The approach is
 based on Li and Racine (2003) who employ \sQuote{generalized product kernels}
 that admit a mix of continuous and discrete data types.

 Three classes of kernel estimators for the continuous data types are
 available: fixed, adaptive nearest-neighbor, and generalized
 nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
 each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
 density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
 with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
 are constant over the support of \eqn{x}.

 \code{npplregbw} may be invoked \emph{either} with a formula-like
 symbolic 
 description of variables on which bandwidth selection is to be
 performed \emph{or} through a simpler interface whereby data is passed
 directly to the function via the \code{xdat}, \code{ydat}, and
 \code{zdat} 
 parameters. Use of these two interfaces is \bold{mutually exclusive}.

 Data contained in the data frame \code{zdat} may be a mix of continuous
 (default), unordered discrete (to be specified in the data frame
 \code{zdat} using \code{\link{factor}}), and ordered discrete (to be
 specified in the data frame \code{zdat} using
 \code{\link{ordered}}). Data can be entered in an arbitrary order and
 data types will be detected automatically by the routine (see
 \code{\link{np}} for details).

 Data for which bandwidths are to be estimated may be specified
 symbolically. A typical description has the form \code{dependent}
 \code{data} \code{~} \code{parametric} \code{explanatory} \code{data}
 \code{|} \code{nonparametric} \code{explanatory} \code{data},
 where \code{dependent} \code{data} is a univariate response, and
 \code{parametric} \code{explanatory} \code{data} and
 \code{nonparametric} \code{explanatory}
 \code{data} are both series of variables specified by name, separated by
 the separation character '+'. For example, \code{ y1 ~ x1 + x2 | z1 }
 specifies that the bandwidth object for the partially linear model with
 response \code{y1}, linear parametric regressors \code{x1} and
 \code{x2}, and 
 nonparametric regressor \code{z1} is to be estimated. See below for
 further examples.  

 A variety of kernels may be specified by the user. Kernels implemented
 for continuous data types include the second, fourth, sixth, and eighth
 order Gaussian and Epanechnikov kernels, and the uniform
 kernel. Unordered discrete data types use a variation on Aitchison and
 Aitken's (1976) kernel, while ordered data types use a variation of the
 Wang and van Ryzin (1981) kernel.
}
\value{
  if \code{bwtype} is set to \code{fixed}, an object containing bandwidths
  (or scale factors if \code{bwscaling = TRUE}) is returned. If it is set to
  \code{generalized_nn} or \code{adaptive_nn}, then instead the \eqn{k}th nearest
  neighbors are returned for the continuous variables while the discrete
  kernel bandwidths are returned for the discrete variables. Bandwidths
  are stored in a list under the component name \code{bw}. Each element
  is an \code{rbandwidth} object. The first
  element of the list corresponds to the regression of \eqn{Y} on \eqn{Z}.  
  Each subsequent element is the bandwidth object corresponding to the
  regression of the \eqn{i}th column of \eqn{X} on \eqn{Z}. See examples
  for more information.
}

\references{

 Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
 discrimination by the kernel method,} Biometrika, 63, 413-420.

 Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
 and Practice,} Princeton University Press.

 Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
   nonparametric regression,} Statistica Sinica, 14, 485-512.

 Racine, J.S. and L. Liu (2006), \dQuote{A partially linear kernel
 estimator for categorical data,} manuscript.

 Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
 Cambridge University Press.
 
 Racine, J.S. and Q. Li (2004), \dQuote{Nonparametric estimation of regression
 functions with both categorical and continuous data,} Journal of
 Econometrics, 119, 99-130.

 Robinson, P.M. (1988), \dQuote{Root-n-consistent semiparametric regression,}
 Econometrica, 56, 931-954. 
 
 Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
 for discrete distributions,}  Biometrika, 68, 301-309.

}
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

  Caution: multivariate data-driven bandwidth selection methods are, by
  their nature, \emph{computationally intensive}. Virtually all methods
  require dropping the \eqn{i}th observation from the data set, computing an
  object, repeating this for all observations in the sample, then
  averaging each of these leave-one-out estimates for a \emph{given}
  value of the bandwidth vector, and only then repeating this a large
  number of times in order to conduct multivariate numerical
  minimization/maximization. Furthermore, due to the potential for local
  minima/maxima, \emph{restarting this procedure a large number of times may
  often be necessary}. This can be frustrating for users possessing
  large datasets. For exploratory purposes, you may wish to override the
  default search tolerances, say, setting ftol=.01 and tol=.01 and
  conduct multistarting (the default is to restart min(5, ncol(zdat))
  times) as is done for a number of examples. Once the procedure
  terminates, you can restart search with default tolerances using those
  bandwidths obtained from the less rigorous search (i.e., set
  \code{bws=bw} on subsequent calls to this routine where \code{bw} is
  the initial bandwidth object).  A version of this package using the
  \code{Rmpi} wrapper is under development that allows one to deploy
  this software in a clustered computing environment to facilitate
  computation involving large datasets.
}

\seealso{
  \code{\link{npregbw}}, \code{\link{npreg}} }
\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we simulate an
# example for a partially linear model and perform bandwidth selection

set.seed(123)

n <- 250
x1 <- rnorm(n)
x2 <- rbinom(n, 5, .3)

z1 <- rbinom(n, 2, .3)
z2 <- rnorm(n)

y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

X <- data.frame(x1, ordered(x2))
Z <- data.frame(ordered(z1), z2)

# Compute data-driven bandwidths... this may take a minute or two
# depending on the speed of your computer...

bw <- npplregbw(formula=y~x1+ordered(x2)|ordered(z1)+z2, tol=.1, ftol=.1)

summary(bw)

\dontrun{

# Note - the default is to use the local constant estimator. If you wish
# to use instead a local linear estimator, this is accomplished via
# npplregbw(xdat=X, zdat=Z, ydat=y, regtype="ll", tol=.1, ftol=.1)

# Note - see the example for npudensbw() for multiple illustrations
# of how to change the kernel function, kernel order, and so forth.

# You may want to manually specify your bandwidths
bw.mat <- matrix(data =  c(0.19, 0.34,  # y on Z
                           0.00, 0.74,  # X[,1] on Z
                           0.99, 0.23), # X[,2] on Z
                ncol = ncol(Z), byrow=TRUE)

bw <- npplregbw(formula=y~x1+ordered(x2)|ordered(z1)+z2, 
                          bws=bw.mat, bandwidth.compute=FALSE)
summary(bw)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# You may want to tweak some of the bandwidths
bw$bw[[1]] # y on Z, alternatively bw$bw$yzbw
bw$bw[[1]]$bw <- c(0.17, 0.30)

bw$bw[[2]] # X[,1] on Z
bw$bw[[2]]$bw[1] <- 0.00054

summary(bw)

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we simulate an
# example for a partially linear model and perform bandwidth selection

set.seed(123)

n <- 250
x1 <- rnorm(n)
x2 <- rbinom(n, 5, .3)

z1 <- rbinom(n, 2, .3)
z2 <- rnorm(n)

y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

X <- data.frame(x1, ordered(x2))
Z <- data.frame(ordered(z1), z2)

# Compute data-driven bandwidths... this may take a minute or two
# depending on the speed of your computer...

bw <- npplregbw(xdat=X, zdat=Z, ydat=y, tol=.1, ftol=.1)

summary(bw)

# Note - the default is to use the local constant estimator. If you wish
# to use instead a local linear estimator, this is accomplished via
# npplregbw(xdat=X, zdat=Z, ydat=y, regtype="ll", tol=.1, ftol=.1)

# Note - see the example for npudensbw() for multiple illustrations
# of how to change the kernel function, kernel order, and so forth.

# You may want to manually specify your bandwidths
bw.mat <- matrix(data =  c(0.19, 0.34,  # y on Z
                           0.00, 0.74,  # X[,1] on Z
                           0.99, 0.23), # X[,2] on Z
                ncol = ncol(Z), byrow=TRUE)

bw <- npplregbw(xdat=X, zdat=Z, ydat=y, 
                          bws=bw.mat, bandwidth.compute=FALSE)
summary(bw)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# You may want to tweak some of the bandwidths
bw$bw[[1]] # y on Z, alternatively bw$bw$yzbw
bw$bw[[1]]$bw <- c(0.17, 0.30)

bw$bw[[2]] # X[,1] on Z
bw$bw[[2]]$bw[1] <- 0.00054

summary(bw)
} % enddontrun
}
\keyword{nonparametric}
