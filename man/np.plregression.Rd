% $ID$
\name{npplreg}
\alias{npplreg}
\alias{npplreg.call}
\alias{npplreg.formula}
\alias{npplreg.plbandwidth}

\title{Partially Linear Kernel Regression with Mixed Data Types}

\description{
  \code{npplreg} computes a partially linear kernel regression estimate
  of a one (1) dimensional dependent variable on \eqn{p+q}-variate
  explanatory data, using the model \eqn{Y = X\beta + \Theta (Z) +
  \epsilon}{Y = XB + theta(Z) + epsilon} given a set of estimation
  points, training points (consisting of explanatory data and dependent
  data), and a bandwidth specification, which can be a \code{rbandwidth}
  object, or a bandwidth vector, bandwidth type and kernel type.
}

\usage{
npplreg(bws, \dots)

\method{npplreg}{formula}(bws, data = NULL, newdata = NULL, y.eval =
FALSE, \dots)

\method{npplreg}{call}(bws, \dots)

\method{npplreg}{plbandwidth}(bws,
        txdat = stop("training data txdat missing"),
        tydat = stop("training data tydat missing"),
        tzdat = stop("training data tzdat missing"),
        exdat,
        eydat,
        ezdat,
        residuals = FALSE,
        \dots)
}
\arguments{

  \item{bws}{
  a bandwidth specification. This can be set as a \code{plbandwidth}
  object returned from an invocation of \code{\link{npplregbw}}, or as
  a matrix of bandwidths, where each row is a set of bandwidths for \eqn{Z},
  with a column for each variable \eqn{Z_i}{Z[i]}. In the first row are
  the bandwidths for the regression of \eqn{Y} on \eqn{Z}, the following
  rows contain the bandwidths for the regressions of the columns of
  \eqn{X} on \eqn{Z}.  If specified as a matrix
  additional arguments will need to be supplied as necessary to
  specify the bandwidth type, kernel types, training data, and so on.

  }

  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, selection methods, and so on. To do
    this, you may specify any of \code{regtype}, 
    \code{bwmethod}, \code{bwscaling}, \code{bwtype}, \code{ckertype},
    \code{ckerorder}, \code{ukertype}, \code{okertype}, as described in
    \code{\link{npregbw}}.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
   \code{\link{npplregbw}} was called.
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
  a \eqn{p}-variate data frame of explanatory data (training data),
  corresponding to \eqn{X} in the model equation, whose linear
  relationship with the dependent data \eqn{Y} is posited. Defaults to
  the training data used to 
    compute the bandwidth object.
  }

  \item{tydat}{
  a one (1) dimensional numeric or integer vector of dependent data, each
  element \eqn{i} corresponding to each observation (row) \eqn{i} of
  \code{txdat}. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{tzdat}{
  a \eqn{q}-variate data frame of explanatory data (training data),
  corresponding to \eqn{Z} in the model equation, whose relationship
  to the dependent variable is unspecified (nonparametric). Defaults to
  the training data used to 
    compute the bandwidth object.
  }

  \item{exdat}{
  a \eqn{p}-variate data frame of points on which the regression will be
  estimated (evaluation data).  By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{eydat}{
  a one (1) dimensional numeric or integer vector of the true values
  of the dependent variable. Optional, and used only to calculate the
  true errors.  By default,
    evaluation takes place on the data provided by \code{tydat}.
  }

  \item{ezdat}{
    a \eqn{q}-variate data frame of points on which the regression will
    be estimated (evaluation data).  By default,
    evaluation takes place on the data provided by \code{tzdat}.
  }

  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{plregression} object. Defaults to
    \code{FALSE}.
  }

}

\details{
 \code{npplreg} uses a combination of OLS and nonparametric
 regression to estimate the parameter \eqn{\beta}{B} in the model
 \eqn{Y = X\beta + \Theta (Z) + \epsilon}{Y = XB + theta(Z) + epsilon}.

 \code{npplreg} implements a variety of methods for
 nonparametric regression on multivariate (\eqn{q}-variate) explanatory
 data defined over a set of possibly continuous and/or discrete
 (unordered, ordered) data. The approach is based on Li and Racine
 (2003) who employ \sQuote{generalized product kernels} that admit a mix
 of continuous and discrete data types.

 Three classes of kernel estimators for the continuous data types are
 available: fixed, adaptive nearest-neighbor, and generalized
 nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
 each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
 density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
 with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
 are constant over the support of \eqn{x}.

 Data contained in the data frame \code{tzdat} may be a mix of
 continuous (default), unordered discrete (to be specified in the data
 frame \code{tzdat} using \code{\link{factor}}), and ordered discrete
 (to be specified in the data frame \code{tzdat} using
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
  \code{npplreg} returns a \code{plregression} object. The generic
  accessor functions \code{\link{coef}}, \code{\link{fitted}},
  \code{\link{residuals}}, \code{\link{predict}}, and
  \code{\link{vcov}}, extract (or
  estimate) coefficients, estimated values, residuals,
  predictions, and variance-covariance matrices,
  respectively, from
  the returned object. Furthermore, the functions \code{summary}
  and \code{plot} support objects of this type. The returned object
  has the following components:
  
  \item{evalx}{ evaluation points }
  \item{evalz}{ evaluation points }
  \item{mean}{ estimation of the regression, or conditional mean, at the
    evaluation points } 
  \item{xcoef}{ coefficient(s) corresponding to the components
    \eqn{\beta_i}{B[i]} in the model }
  \item{xcoeferr}{ standard errors of the coefficients }
  \item{xcoefvcov}{ covariance matrix of the coefficients }
  \item{bw}{ the bandwidths, stored as a \code{plbandwidth} object }
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
 global functionals and a measure of the explanatory power of covariates
 in regression,} The Annals of Statistics, 23 1443-1473.

 Gao, Q. and L. Liu and J.S. Racine (forthcoming), \dQuote{A
 partially linear kernel estimator for categorical data,} Econometric
 Reviews.

 Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
 and Practice,} Princeton University Press.

 Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
   nonparametric regression,} Statistica Sinica, 14, 485-512.
 
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
  \code{\link{npregbw}}, \code{\link{npreg}} }
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
x1 <- rnorm(n)
x2 <- rbinom(n, 1, .5)

z1 <- rbinom(n, 1, .5)
z2 <- rnorm(n)

y <- 1 + x1 + x2 + z1 + sin(z2) + rnorm(n)

x2 <- factor(x2)
z1 <- factor(z1)

mpi.bcast.Robj2slave(x1)
mpi.bcast.Robj2slave(x2)
mpi.bcast.Robj2slave(z1)
mpi.bcast.Robj2slave(z2)
mpi.bcast.Robj2slave(y)

mpi.bcast.cmd(bw <- npplregbw(formula=y~x1+x2|z1+z2),
              caller.execute=TRUE)

mpi.bcast.cmd(pl <- npplreg(bws=bw),
              caller.execute=TRUE)

summary(pl)

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
