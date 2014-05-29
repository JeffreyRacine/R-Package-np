% $Id: np.sigtest.Rd,v 1.31 2006/10/24 20:43:29 tristen Exp $

\name{npsigtest}
\alias{npsigtest}
\alias{npsigtest.call}
\alias{npsigtest.default}
\alias{npsigtest.formula}
\alias{npsigtest.rbandwidth}
\alias{npsigtest.npregression}

\title{Kernel Regression Significance Test with Mixed Data Types}
\description{
  \code{npsigtest} implements a consistent test of significance of an
  explanatory variable(s) in a nonparametric regression setting that is
  analogous to a simple \eqn{t}-test (\eqn{F}-test) in a parametric
  regression setting. The test is based on Racine, Hart, and Li (2006)
  and Racine (1997).
}
\usage{
npsigtest(bws, \dots)

\method{npsigtest}{formula}(bws, data = NULL, \dots)

\method{npsigtest}{call}(bws, \dots)

\method{npsigtest}{npregression}(bws, \dots)

\method{npsigtest}{default}(bws, xdat, ydat, \dots)

\method{npsigtest}{rbandwidth}(bws,
          xdat = stop("data xdat missing"),
          ydat = stop("data ydat missing"),
          boot.num = 399,
          boot.method = c("iid","wild","wild-rademacher","pairwise"),
          boot.type = c("I","II"),
          pivot=TRUE,
          joint=FALSE,
          index = seq(1,ncol(xdat)),
          random.seed = 42,
          \dots)

}

\arguments{
  \item{bws}{
    a bandwidth specification. This can be set as a \code{rbandwidth}
    object returned from a previous invocation, or as a vector of
    bandwidths, with each element \eqn{i} corresponding to the bandwidth
    for column \eqn{i} in \code{xdat}. In either case, the bandwidth
    supplied will serve as a starting point in the numerical search for
    optimal bandwidths when using \code{boot.type="II"}. If specified
    as a vector, then additional arguments will need to be supplied as
    necessary to specify the bandwidth type, kernel types, selection
    methods, and so on.
  }

  \item{data}{
    an optional data frame, list or environment (or object coercible to
    a data frame by \code{\link{as.data.frame}}) containing the
    variables in the model. If not found in data, the variables are
    taken from \code{environment(bws)}, typically the environment from
    which \code{\link{npregbw}} was called.
  }

  \item{xdat}{
    a \eqn{p}-variate data frame of explanatory data (training data)
    used to calculate the regression estimators.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data,
    each element \eqn{i} corresponding to each observation (row) \eqn{i}
    of \code{xdat}.
  }

  \item{boot.method}{
    a character string used to specify the bootstrap method for
    determining the null distribution. \code{pairwise} resamples
    pairwise, while the remaining methods use a residual bootstrap
    procedure.  \code{iid} will generate independent identically
    distributed draws. \code{wild} will use a wild
    bootstrap. \code{wild-rademacher} will use a wild bootstrap with
    Rademacher variables. Defaults to \code{iid}.
  }

  \item{boot.num}{
    an integer value specifying the number of bootstrap replications to
    use. Defaults to \code{399}.
  }

  \item{boot.type}{
    a character string specifying whether to use a `Bootstrap I' or
    `Bootstrap II' method (see Racine, Hart, and Li (2006) for
    details). The `Bootstrap II' method re-runs cross-validation for
    each bootstrap replication and uses the new cross-validated
    bandwidth for variable \eqn{i} and the original ones for the
    remaining variables. Defaults to \code{boot.type="I"}.
  }

  \item{pivot}{
    a logical value which specifies whether to bootstrap a pivotal
    statistic or not (pivoting is achieved by dividing gradient
    estimates by their asymptotic standard errors). Defaults to
    \code{TRUE}.
  }

  \item{joint}{
    a logical value which specifies whether to conduct a joint test or
    individual test. This is to be used in conjunction with \code{index}
    where \code{index} contains two or more integers corresponding to
    the variables being tested, where the integers correspond to the
    variables in the order in which they appear among the set of
    explanatory variables in the function call to
    \code{npreg}/\code{npregbw}. Defaults to \code{FALSE}.
  }

  \item{index}{
    a vector of indices for the columns of \code{xdat} for which the
    test of significance is to be conducted. Defaults to
    (1,2,\dots,\eqn{p}) where \eqn{p} is the number of columns in
    \code{xdat}.
  }

  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }

  \item{\dots}{
    additional arguments supplied to specify the bandwidth type,
    kernel types, selection methods, and so on, detailed below.
  }

}
\details{

  \code{npsigtest} implements a variety of methods for computing the
  null distribution of the test statistic and allows the user to
  investigate the impact of a variety of default settings including
  whether or not to pivot the statistic (\code{pivot}), whether pairwise
  or residual resampling is to be used (\code{boot.method}), and whether
  or not to recompute the bandwidths for the variables being tested
  (\code{boot.type}), among others.

  Defaults are chosen so as to provide reasonable behaviour in a broad
  range of settings and this involves a trade-off between computational
  expense and finite-sample performance. However, the default
  \code{boot.type="I"}, though computationally expedient, can deliver a
  test that can be slightly over-sized in small sample settings (e.g.\
  at the 5\% level the test might reject 8\% of the time for samples of
  size \eqn{n=100} for some data generating processes).  If the default
  setting (\code{boot.type="I"}) delivers a P-value that is in the
  neighborhood (i.e.\ slightly smaller) of any classical level
  (e.g. 0.05) and you only have a modest amount of data, it might be
  prudent to re-run the test using the more computationally intensive
  \code{boot.type="II"} setting to confirm the original result. Note
  also that \code{boot.method="pairwise"} is not recommended for the
  multivariate local linear estimator due to substantial size
  distortions that may arise in certain cases.
  
}  

\value{
  \code{npsigtest} returns an object of type
  \code{sigtest}. \code{\link{summary}} supports \code{sigtest} objects. It
  has the
  following components:

  \item{In}{ the vector of statistics \code{In} }
  \item{P}{ the vector of P-values for each statistic in \code{In}}
  \item{In.bootstrap}{ contains a matrix of the bootstrap
    replications of the vector \code{In}, each column corresponding to
    replications associated with explanatory variables in \code{xdat}
    indexed by \code{index} (e.g., if you selected \code{index = c(1,4)}
    then In.bootstrap will have two columns, the first being the
    bootstrap replications of \code{In} associated with variable
    \code{1}, the second with variable \code{4}).}

}
\references{
  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
    Theory and Practice,} Princeton University Press.

  Racine, J.S., J. Hart, and Q. Li (2006), \dQuote{Testing the
    significance of categorical predictor variables in nonparametric
    regression models,} Econometric Reviews, 25, 523-544.

  Racine, J.S. (1997), \dQuote{Consistent significance testing for
    nonparametric regression,} Journal of Business and Economic
  Statistics 15, 369-379.

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
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work
  as intended on mixed data types and will coerce the data to the same
  type.

  Caution: bootstrap methods are, by their nature, \emph{computationally
  intensive}. This can be frustrating for users possessing large
  datasets. For exploratory purposes, you may wish to override the
  default number of bootstrap replications, say, setting them to
  \code{boot.num=99} A version of this package using the \code{Rmpi}
  wrapper is under development that allows one to deploy this software
  in a clustered computing environment to facilitate computation
  involving large datasets.

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

mpi.bcast.cmd(set.seed(42),
              caller.execute=TRUE)

## Significance testing with z irrelevant

n <- 250

z <- factor(rbinom(n,1,.5))
x1 <- rnorm(n)
x2 <- runif(n,-2,2)
y <- x1 + x2 + rnorm(n)
mydat <- data.frame(z,x1,x2,y)
rm(z,x1,x2,y)

mpi.bcast.Robj2slave(mydat)

mpi.bcast.cmd(model <- npreg(y~z+x1+x2,
                             regtype="ll",
                             bwmethod="cv.aic",
                             data=mydat),
              caller.execute=TRUE)

mpi.bcast.cmd(output <- npsigtest(model,boot.num=99),
              caller.execute=TRUE)

summary(output)

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
