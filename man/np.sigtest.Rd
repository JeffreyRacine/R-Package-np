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
  \code{npsigtest} implements a consistent test of significance of
  an explanatory variable in a nonparametric regression setting that is
  analogous to a simple \eqn{t}-test in a parametric regression
  setting. The test is based on Racine, Hart, and Li (2006) and
  Racine (1997).
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
          boot.method = c("iid","wild","wild-rademacher"),
          boot.type = c("I","II"),
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
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npregbw}} was called.
  }

  \item{xdat}{
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the regression estimators.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{xdat}. 
  }


  \item{boot.method}{
    a character string used to specify the bootstrap method.
    \code{iid} will generate independent identically distributed
    draws. \code{wild} will use a wild bootstrap. \code{wild-rademacher}
    will use a wild bootstrap with Rademacher variables. Defaults to
    \code{iid}.
           
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

  \item{index}{
    a vector of indices for the columns of \code{xdat} for which the
    test of significance is to be conducted. Defaults to (1,2,\dots,\eqn{p})
    where \eqn{p} is the number of columns in \code{xdat}.
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
    Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.

  Caution: bootstrap methods are, by their nature, \emph{computationally
    intensive}. This can be frustrating for users possessing large
  datasets. For exploratory purposes, you may wish to override the
  default number of bootstrap replications, say, setting them to
  boot.num=99 A version of this package using the \code{Rmpi} wrapper is
  under development that allows one to deploy this software in a
  clustered computing environment to facilitate computation involving
  large datasets.
}

\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we simulate 100 draws
# from a DGP in which z, the first column of X, is an irrelevant
# discrete variable

set.seed(12345)

n <- 100

z <- rbinom(n,1,.5)
x1 <- rnorm(n)
x2 <- runif(n,-2,2)

y <- x1 + x2 + rnorm(n)

# Next, we must compute bandwidths for our regression model. In this
# case we conduct local linear regression. Note - this may take a few
# minutes depending on the speed of your computer...

bw <- npregbw(formula=y~factor(z)+x1+x2,regtype="ll",bwmethod="cv.aic")

# We then compute a vector of tests corresponding to the columns of
# X. Note - this may take a few minutes depending on the speed of your
# computer... we have to generate the null distribution of the statistic
# for each variable whose significance is being tested using 399
# bootstrap replications for each...

npsigtest(bws=bw)

\dontrun{

# If you wished, you could conduct the test for, say, variables 1 and 3
# only, as in

npsigtest(bws=bw,index=c(1,3))

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we simulate 100
# draws from a DGP in which z, the first column of X, is an irrelevant
# discrete variable

set.seed(12345)

n <- 100

z <- rbinom(n,1,.5)
x1 <- rnorm(n)
x2 <- runif(n,-2,2)

X <- data.frame(factor(z),x1,x2)

y <- x1 + x2 + rnorm(n)

# Next, we must compute bandwidths for our regression model. In this
# case we conduct local linear regression. Note - this may take a few
# minutes depending on the speed of your computer...

bw <- npregbw(xdat=X,ydat=y,regtype="ll",bwmethod="cv.aic")

# We then compute a vector of tests corresponding to the columns of
# X. Note - this may take a few minutes depending on the speed of your
# computer... we have to generate the null distribution of the statistic
# for each variable whose significance is being tested using 399
# bootstrap replications for each...

npsigtest(bws=bw)

# If you wished, you could conduct the test for, say, variables 1 and 3
# only, as in

npsigtest(bws=bw,index=c(1,3))
} % enddontrun
}
\keyword{ nonparametric }
