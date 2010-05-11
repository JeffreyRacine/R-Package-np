
% $Id: np.qcmstest.Rd,v 1.41 2006/11/09 11:43:39 jracine Exp $
\name{npqcmstest}
\alias{npqcmstest}
\title{ Kernel Consistent Quantile Regression Model Specification Test with Mixed Data Types }
\description{
  \code{npqcmstest} implements a consistent test for correct
  specification of parametric quantile regression models (linear or
  nonlinear) as described in Racine (2006) which extends the work of
  Zheng (1998).
}
\usage{
npqcmstest(formula,
           data = NULL,
           subset,
           xdat,
           ydat,
           model = stop(paste(sQuote("model")," has not been provided")),
           tau = 0.5,
           distribution = c("bootstrap", "asymptotic"),
           bwydat = c("y","varepsilon"),
           boot.method = c("iid","wild","wild-rademacher"),
           boot.num = 399,
           pivot = TRUE,
           density.weighted = TRUE,
           random.seed = 42,
           \dots)
}

\arguments{
    \item{formula}{
    a symbolic description of variables on which the test is
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
    an optional vector specifying a subset of observations to be used. 
  }

  \item{model}{
    a model object obtained from a call to \code{\link[quantreg]{rq}}. Important: the
    call to \code{\link[quantreg]{rq}} must have the argument \code{model=TRUE} or
    \code{npqcmstest} will not work.
  }

  \item{xdat}{
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the quantile regression estimators.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{xdat}. 
  }

  \item{tau}{
    a numeric value specifying the \eqn{\tau}{tau}th quantile is desired
  }

  \item{distribution}{
    a character string used to specify the method of estimating the
    distribution of the statistic to be calculated. \code{bootstrap}
    will conduct bootstrapping. \code{asymptotic} will use the normal
    distribution. Defaults to \code{bootstrap}.
  }

  \item{bwydat}{
    a character string used to specify \code{ydat} used in bandwidth
    selection. \code{varepsilon}
    uses \eqn{1-\tau,-\tau}{1-tau,-tau} for \code{ydat} while
    \code{y} will use \eqn{y}. Defaults to \code{y}.
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

  \item{pivot}{
    a logical value specifying whether the statistic should be
    normalised such that it approaches \eqn{N(0,1)} in
    distribution. data. Defaults to \code{TRUE}.
  }

  \item{density.weighted}{
    a logical value specifying whether the statistic should be
    weighted by the density of \code{xdat}. Defaults to \code{TRUE}.
  }

  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }

  \item{\dots}{
    additional arguments supplied to control bandwidth selection on the
    residuals. One can specify the bandwidth type,
    kernel types, and so on. To do this, you may specify any of \code{bwscaling},
    \code{bwtype}, \code{ckertype}, \code{ckerorder}, \code{ukertype},
    \code{okertype}, as described in \code{\link{npregbw}}.
    This is necessary if you specify \code{bws} as a \eqn{p}-vector and not
    a \code{bandwidth} object, and you do not desire the default behaviours.
  }
}

\value{
  \code{npqcmstest} returns an object of type \code{cmstest} with the
  following components, components will contain information
  related to \code{Jn} or \code{In} depending on the value of \code{pivot}:

  \item{Jn}{ the statistic \code{Jn} }
  \item{In}{ the statistic \code{In} }
  \item{Omega.hat}{ as described in Racine, J.S. (2006). }
  \item{q.*}{ the various quantiles of the statistic \code{Jn} (or
    \code{In} if 
    \code{pivot=FALSE})  are in
    components \code{q.90},
    \code{q.95}, \code{q.99} (one-sided 1\%, 5\%, 10\% critical values) }
  \item{P}{ the P-value of the statistic }
  \item{Jn.bootstrap}{ if \code{pivot=TRUE} contains the bootstrap
    replications of \code{Jn} }
  \item{In.bootstrap}{ if \code{pivot=FALSE} contains the bootstrap
    replications of \code{In} }
  
  \code{\link{summary}} supports object of type \code{cmstest}.
}
\references{
  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Koenker, R.W. and G.W. Bassett (1978), \dQuote{Regression quantiles,}
  Econometrica, 46, 33-50.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Murphy, K. M. and F. Welch (1990), \dQuote{Empirical age-earnings
    profiles,} Journal of Labor Economics, 8, 202-229.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.

  Racine, J.S. (2006), \dQuote{Consistent specification testing of
    heteroskedastic parametric regression quantile models with mixed
    data,} manuscript.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
    for discrete distributions,}  Biometrika, 68, 301-309.

  Zheng, J. (1998), \dQuote{A consistent nonparametric test of
    parametric regression models under conditional quantile
    restrictions,} Econometric Theory, 14, 123-138.

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
}

\examples{
# EXAMPLE 1: For this example, we conduct a consistent quantile regression
# model specification test for a parametric wage quantile regression
# model that is quadratic in age. The work of Murphy and Welch (1990)
# would suggest that this parametric quantile regression model is
# misspecified.

library("quantreg")

data("cps71")
attach(cps71)

model <- rq(logwage~age+I(age^2), tau=0.5, model=TRUE)

plot(age, logwage)
lines(age, fitted(model))

X <- data.frame(age)

# Note - this may take a few minutes depending on the speed of your
# computer...

npqcmstest(model = model, xdat = X, ydat = logwage, tau=0.5)

\dontrun{

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Next try Murphy & Welch's (1990) suggested quintic specification.

model <- rq(logwage~age+I(age^2)+I(age^3)+I(age^4)+I(age^5), model=TRUE)

plot(age, logwage)
lines(age, fitted(model))

X <- data.frame(age)

# Note - this may take a few minutes depending on the speed of your
# computer...

npqcmstest(model = model, xdat = X, ydat = logwage, tau=0.5)

detach(cps71)
} % enddontrun
}
\keyword{ nonparametric }
