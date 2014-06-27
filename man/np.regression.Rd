% $ID$
\name{npreg}
\alias{npreg}
\alias{npreg.call}
\alias{npreg.default}
\alias{npreg.formula}
\alias{npreg.rbandwidth}

\title{Kernel Regression with Mixed Data Types}

\description{
  \code{npreg} computes a kernel regression estimate of a one
  (1) dimensional dependent variable on \eqn{p}-variate explanatory
  data, given a set of evaluation points, training points (consisting of
  explanatory data and dependent data), and a bandwidth specification
  using the method of Racine and Li (2004) and Li and Racine (2004). A
  bandwidth specification can be a \code{rbandwidth} object, or a
  bandwidth vector, bandwidth type and kernel type.
}

\usage{
npreg(bws, \dots)

\method{npreg}{formula}(bws, data = NULL, newdata = NULL, y.eval =
FALSE, \dots)

\method{npreg}{call}(bws, \dots)

\method{npreg}{default}(bws, txdat, tydat, \dots)

\method{npreg}{rbandwidth}(bws,
      txdat = stop("training data 'txdat' missing"),
      tydat = stop("training data 'tydat' missing"),
      exdat,
      eydat,
      gradients = FALSE,
      residuals = FALSE,
      \dots)
}

\arguments{
  \item{bws}{
    a bandwidth specification. This can be set as a \code{rbandwidth}
    object returned from an invocation of \code{\link{npregbw}}, or
    as a vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{txdat}. If specified as
    a vector, then additional arguments will need to be supplied as
    necessary to specify the bandwidth type, kernel types, and so on.
  }
  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, training data, and so on, detailed
    below. 
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npregbw}} was called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.  
  }
  
  \item{y.eval}{
    If \code{newdata} contains dependent data and \code{y.eval = TRUE},
    \code{\link{np}} will compute goodness of fit statistics on these
    data and return them. Defaults to \code{FALSE}.
  }
  
  \item{txdat}{
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the regression estimators. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{tydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{txdat}. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{eydat}{
    a one (1) dimensional numeric or integer vector of the true values
    of the dependent variable. Optional, and used only to calculate the
    true errors. 
  }

  \item{gradients}{
    a logical value indicating that you want gradients computed and
    returned in the resulting \code{npregression} object. Defaults to
    \code{FALSE}.
  }

  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{npregression} object. Defaults to
    \code{FALSE}.
  }
}

\details{

  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    Usage 1: first compute the bandwidth object via npregbw and then
    compute the conditional mean:
    
    bw <- npregbw(y~x)
    ghat <- npreg(bw)
    
    Usage 2: alternatively, compute the bandwidth object indirectly:
    
    ghat <- npreg(y~x)
    
    Usage 3: modify the default kernel and order:
    
    ghat <- npreg(y~x, ckertype="epanechnikov", ckerorder=4)

    Usage 4: use the data frame interface rather than the formula
    interface:

    ghat <- npreg(tydat=y, txdat=x, ckertype="epanechnikov", ckerorder=4)
  }


  
  \code{npreg} implements a variety of methods for regression on
  multivariate (\eqn{p}-variate) data, the types of which are possibly
  continuous and/or discrete (unordered, ordered). The approach is
  based on Li and Racine (2003) who employ \sQuote{generalized product kernels}
  that admit a mix of continuous and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
  density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
  with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
  are constant over the support of \eqn{x}.

  Data contained in the data frame \code{txdat} may be a mix of
  continuous (default), unordered discrete (to be specified in the data
  frame \code{txdat} using \code{\link{factor}}), and ordered discrete
  (to be specified in the data frame \code{txdat} using
  \code{\link{ordered}}). Data can be entered in an arbitrary order and
  data types will be detected automatically by the routine (see
  \code{\link{np}} for details).

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

}

\value{

  \code{npreg} returns a \code{npregression} object.
  The generic
  functions \code{\link{fitted}}, \code{\link{residuals}},
  \code{\link{se}}, \code{\link{predict}}, and
  \code{\link{gradients}}, extract (or generate) estimated values,
  residuals, asymptotic standard
  errors on estimates, predictions, and gradients, respectively, from
  the returned object. Furthermore, the functions \code{\link{summary}}
  and \code{\link{plot}} support objects of this type. The returned object
  has the following components:


  \item{eval}{ evaluation points }
  \item{mean}{ estimates of the regression function (conditional mean) at the
  evaluation points } 
  \item{merr}{ standard errors of the regression function estimates } 
  \item{grad}{ estimates of the gradients at each evaluation point }
  \item{gerr}{ standard errors of the gradient estimates }
  \item{resid}{ if \code{residuals = TRUE}, in-sample or out-of-sample
    residuals where appropriate (or possible) }
  \item{R2}{ coefficient of determination (Doksum and Samarov (1995))}
  \item{MSE}{ mean squared error }
  \item{MAE}{ mean absolute error }
  \item{MAPE}{ mean absolute percentage error }
  \item{CORR}{ absolute value of Pearson's correlation coefficient }
  \item{SIGN}{ fraction of observations where fitted and observed values
  agree in sign }
  
}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
  discrimination by the kernel method,} Biometrika, 63, 413-420.

  Doksum, K. and A. Samarov (1995), \dQuote{Nonparametric estimation of
  global functionals and a measure of the explanatory power of
  covariates in regression,} The Annals of Statistics, 23 1443-1473.  

  Hall, P. and Q. Li and J.S. Racine (2007), \dQuote{Nonparametric
  estimation of regression functions in the presence of irrelevant
  regressors,} The Review of Economics and Statistics, 89, 784-789.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
  nonparametric regression,} Statistica Sinica, 14, 485-512.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.
  
  Racine, J.S. and Q. Li (2004), \dQuote{Nonparametric estimation of
    regression functions with both categorical and continuous data,}
  Journal of Econometrics, 119, 99-130.
  
  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,} Biometrika, 68, 301-309.

}

\author{
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.
}

\seealso{
  \code{\link{loess}}
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

set.seed(42)

n <- 250

x <- runif(n)
z1 <- rbinom(n,1,.5)
z2 <- rbinom(n,1,.5)
y <- cos(2*pi*x) + z1 + rnorm(n,sd=.25)
z1 <- factor(z1)
z2 <- factor(z2)
mydat <- data.frame(y,x,z1,z2)
rm(x,y,z1,z2)

mpi.bcast.Robj2slave(mydat)

mpi.bcast.cmd(bw <- npregbw(y~x+z1+z2,
                            regtype="lc",
                            bwmethod="cv.ls",
                            data=mydat),
              caller.execute=TRUE)

summary(bw)

mpi.bcast.cmd(model <- npreg(bws=bw,
                             data=mydat),
              caller.execute=TRUE)

summary(model)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

mpi.close.Rslaves()

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun     
}
\keyword{nonparametric}

