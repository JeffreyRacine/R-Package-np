%

\name{npsigtest}
\alias{npsigtest}
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
npsigtest(bws,
          \dots,
          B = 399)

\method{npsigtest}{formula}(bws, 
          data = NULL, 
          \dots)

\method{npsigtest}{npregression}(bws, \dots)

\method{npsigtest}{default}(bws, 
          xdat, 
          ydat, 
          \dots)

\method{npsigtest}{rbandwidth}(bws,
          xdat = stop("data xdat missing"),
          ydat = stop("data ydat missing"),
          B = 399,
          boot.method = c("iid","wild","wild-rademacher","pairwise"),
          boot.type = c("I","II"),
          pivot=NULL,
          joint=FALSE,
          index = seq_len(ncol(xdat)),
          random.seed = 42,
          \dots)

}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth object, formula/data interface, and explicit data inputs.
}

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


\subsection{Bootstrap And Test Controls}{
These arguments control bootstrap execution, tested effects, and reproducibility settings.
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


  \item{B}{
    an integer value specifying the number of bootstrap replications to
    use. Defaults to \code{399}.
  }


  \item{boot.type}{
    a character string specifying whether to use a \sQuote{Bootstrap I} or
    \sQuote{Bootstrap II} method (see Racine, Hart, and Li (2006) for
    details). The \sQuote{Bootstrap II} method re-runs cross-validation for
    each bootstrap replication and uses the new cross-validated
    bandwidth for variable \eqn{i} and the original ones for the
    remaining variables. Defaults to \code{boot.type="I"}.
  }


  \item{index}{
    a vector of indices for the columns of \code{xdat} for which the
    test of significance is to be conducted. Defaults to
    (1,2,\dots,\eqn{p}) where \eqn{p} is the number of columns in
    \code{xdat}.
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


  \item{pivot}{
    \code{NULL} or a logical value specifying whether to bootstrap a
    pivotal statistic. With the method-aware default, \code{NULL}, an
    individual test is pivotal for a continuous predictor and
    unstandardized for a categorical predictor. A joint test is pivotal
    only when every tested predictor is continuous; if any tested
    predictor is categorical, the entire joint statistic is
    unstandardized. \code{FALSE} requests the unstandardized statistic
    for every tested predictor. \code{TRUE} requests the pivotal
    statistic and is supported only when every tested predictor is
    continuous; otherwise the function stops with an error. For
    continuous predictors, pivoting divides gradient estimates by their
    asymptotic standard errors.
  }


  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }


\subsection{Additional Arguments}{
Further arguments are passed to the bandwidth-selection routines used by the test.
}

  \item{\dots}{
    additional arguments supplied to specify the bandwidth type,
    kernel types, selection methods, and so on, detailed below.
  }


}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, and \code{\link{plot}} for plotting options.


  \code{npsigtest} implements a variety of methods for computing the
  null distribution of the test statistic and allows the user to
  investigate the impact of a variety of settings including whether or
  not to pivot the statistic (\code{pivot}), whether pairwise
  or residual resampling is to be used (\code{boot.method}), and whether
  or not to recompute the bandwidths for the variables being tested
  (\code{boot.type}), among others.

  Eligible individual tests with fixed bandwidths, IID residual resampling,
  and a supported local-polynomial or categorical-only local-constant engine
  evaluate bootstrap statistics in bounded native tiles. This is an internal
  computational optimization: it preserves bootstrap draws, random-number
  state, statistics, and P-values. Joint tests and configurations outside the
  qualified set continue through the incumbent bootstrap loop.

  The categorical statistic of Racine, Hart, and Li (2006) is the
  unstandardized mean of squared finite-difference estimates. It is not
  divided by the asymptotic standard errors used for continuous
  gradients. Therefore \code{pivot=NULL} selects the statistically
  defined statistic by predictor type, while explicit \code{pivot=TRUE}
  fails for any test containing a categorical predictor. A requested
  continuous pivotal statistic also fails if a required gradient
  standard error is unavailable, non-positive, or non-finite; such a
  value does not define the pivotal statistic.

  Continuous pivotal statistics use the fitted-residual HC0 gradient
  standard errors described in \code{\link{npreg}}. Each bootstrap response
  is fitted and its own training residuals are used in its HC0 quadratic;
  the original-response standard errors are not held fixed across
  replications. The admitted streamed implementation computes this same
  statistic without retaining a full influence matrix. Categorical
  statistics remain unstandardized.

  When the supplied bandwidth object has \code{regtype="lp"}, every
  tested continuous coordinate must have polynomial degree at least one.
  A degree-zero local polynomial does not identify the derivative used by
  this significance test, so the function stops before fitting or
  bootstrap work and names the affected coordinate(s). If the degree was
  selected by NOMAD, refit the bandwidth object with
  \code{degree.min=1}. This restriction does not apply to an ordinary
  \code{regtype="lc"} bandwidth object.

  With nearest-neighbor bandwidths, repeated continuous values can give a
  required radius of zero. In this case the test stops with the shared
  nearest-neighbor diagnostic, identifying the affected variable and the
  failing stage. This can occur during the unrestricted gradient or
  null-model evaluation, before bootstrap replications. It is not a blanket
  restriction on repeated values. Refit using \code{bwtype="fixed"} (the
  default) to avoid this nearest-neighbor radius limitation; see
  \code{\link{np.kernels}}.

  The empirical upper-tail P-value counts bootstrap statistics equal to
  the observed statistic, that is,
  \eqn{P=B^{-1}\sum_{b=1}^B 1\{I_b^*\ge I\}}. Thus equality with the
  bootstrap critical value is non-rejection. In particular, when
  smoothing makes both the observed statistic and every bootstrap
  statistic exactly zero, the returned P-value is one rather than zero.

  Defaults are chosen so as to provide reasonable behaviour in a broad
  range of settings and this involves a trade-off between computational
  expense and finite-sample performance. However, the default
  \code{boot.type="I"}, though computationally expedient, can deliver a
  test that can be slightly over-sized in small sample settings (e.g.
  at the 5\% level the test might reject 8\% of the time for samples of
  size \eqn{n=100} for some data generating processes).  If the default
  setting (\code{boot.type="I"}) delivers a P-value that is in the
  neighborhood (i.e. slightly smaller) of any classical level
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
  \item{P}{ the vector of empirical upper-tail P-values for each statistic
    in \code{In}; bootstrap values equal to the observed statistic are
    included in the upper tail.}
  \item{In.bootstrap}{ contains a matrix of the bootstrap
    replications of the vector \code{In}, each column corresponding to
    replications associated with explanatory variables in \code{xdat}
    indexed by \code{index} (e.g., if you selected \code{index = c(1,4)}
    then In.bootstrap will have two columns, the first being the
    bootstrap replications of \code{In} associated with variable
    \code{1}, the second with variable \code{4}).}
  \item{pivot}{the requested pivot mode: \code{NULL}, \code{TRUE}, or
    \code{FALSE}.}
  \item{pivot.effective}{a named logical vector recording the statistic
    actually used for each tested predictor. For a joint test all entries
    are identical.}

}
\section{Book And Method Pointers}{
  \code{npsigtest} assesses predictor relevance in nonparametric
  regression by comparing the fitted unrestricted regression with
  restrictions implied by excluding or smoothing over selected
  regressors. Bootstrap choices trade computational cost against
  finite-sample behavior.

  For book-length background, see Li and Racine (2007), Chapter 12
  \emph{Model Specification Tests}, especially Sections 12.3 and 12.3.5,
  and Chapter 13 \emph{Nonsmoothing Tests}, especially Section 13.3 where
  relevant. For the underlying conditional-mean setting, see Racine
  (2019), Chapter 6 \emph{Conditional Mean Function Estimation}.
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

  When reusing a bandwidth object without explicit training or evaluation
  data, the required data are recovered from object-retained state or
  from the stored call environment when available. This is a convenience
  for ordinary same-session use, not a guarantee that a saved object is
  self-contained. For reproducible replay after \code{saveRDS()} and
  \code{readRDS()} in a fresh R session, supply the required data
  arguments explicitly.

  Caution: bootstrap methods are, by their nature, \emph{computationally
  intensive}. This can be frustrating for users possessing large
  datasets. For exploratory purposes, you may wish to override the
  default number of bootstrap replications, say, setting them to
  \code{B=99}. A version of this package using the \code{Rmpi}
  wrapper is under development that allows one to deploy this software
  in a clustered computing environment to facilitate computation
  involving large datasets.

}

\examples{
\dontrun{
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
