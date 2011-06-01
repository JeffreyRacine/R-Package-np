% $Id: np.cmstest.Rd,v 1.58 2006/11/03 21:17:20 tristen Exp $
\name{npcmstest}
\alias{npcmstest}
\title{ Kernel Consistent Model Specification Test with Mixed Data Types }
\description{
  \code{npcmstest} implements a consistent test for correct
  specification of parametric regression models (linear or nonlinear) as
  described in Hsiao, Li, and Racine (2007).
}
\usage{
npcmstest(formula,
          data = NULL,
          subset,
          xdat,
          ydat,
          model = stop(paste(sQuote("model")," has not been provided")),
          distribution = c("bootstrap", "asymptotic"),
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
    a model object obtained from a call to \code{\link{lm}} (or
    \code{\link{glm}}). Important: the 
    call to either \code{\link{glm}} or \code{\link{lm}} must have the arguments
    \code{x=TRUE} and 
    \code{y=TRUE} or \code{npcmstest} will not work. Also, the test is
    based on residual bootstrapping hence the outcome must be continuous
    (which rules out Logit, Probit, and Count models).
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
  \item{distribution}{
    a character string used to specify the method of estimating the
    distribution of the statistic to be calculated. \code{bootstrap}
    will conduct bootstrapping. \code{asymptotic} will use the normal
    distribution. Defaults to \code{bootstrap}.
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
    distribution. Defaults to \code{TRUE}.
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
  \code{npcmstest} returns an object of type \code{cmstest} with the
  following components, components will contain information
  related to \code{Jn} or \code{In} depending on the value of \code{pivot}:

  \item{Jn}{ the statistic \code{Jn} }
  \item{In}{ the statistic \code{In} }
  \item{Omega.hat}{ as described in Hsiao, C. and Q. Li and J.S. Racine. }
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

  Hsiao, C. and Q. Li and J.S. Racine (2007), \dQuote{A consistent
    model specification test with mixed categorical and continuous
    data,} Journal of Econometrics, 140, 802-826.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
    and Practice,} Princeton University Press.

  Maasoumi, E. and J.S. Racine and T. Stengos (2007),
  \dQuote{Growth and convergence: a profile of distribution dynamics and
    mobility,} Journal of Econometrics, 136, 483-508.

  Murphy, K. M. and F. Welch (1990), \dQuote{Empirical age-earnings
    profiles,} Journal of Labor Economics, 8, 202-229.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.

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
}

\examples{
# EXAMPLE 1: For this example, we conduct a consistent model
# specification test for a parametric wage regression model that is
# quadratic in age. The work of Murphy and Welch (1990) would suggest
# that this parametric regression model is misspecified.

data("cps71")
attach(cps71)

model <- lm(logwage~age+I(age^2), x=TRUE, y=TRUE)

plot(age, logwage)
lines(age, fitted(model))

# Note - this may take a few minutes depending on the speed of your
# computer...

npcmstest(model = model, xdat = age, ydat = logwage)

\dontrun{

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Next try Murphy & Welch's (1990) suggested quintic specification.

model <- lm(logwage~age+I(age^2)+I(age^3)+I(age^4)+I(age^5), x=TRUE, y=TRUE)

plot(age, logwage)
lines(age, fitted(model))

X <- data.frame(age)

# Note - this may take a few minutes depending on the speed of your
# computer...

npcmstest(model = model, xdat = age, ydat = logwage)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Note - you can pass in multiple arguments to this function. For
# instance, to use local linear rather than local constant regression, 
# you would use npcmstest(model, X, regtype="ll"), while you could also
# change the kernel type (default is second order Gaussian), numerical
# search tolerance, or feed in your own vector of bandwidths and so
# forth.

detach(cps71)

# EXAMPLE 2: For this example, we replicate the application in Maasoumi,
# Racine, and Stengos (2007) (see oecdpanel for details). We
# estimate a parametric model that is used in the literature, then
# subject it to the model specification test.

data("oecdpanel")
attach(oecdpanel)

model <- lm(growth ~ oecd +
            factor(year) +
            initgdp +
            I(initgdp^2) +
            I(initgdp^3) +
            I(initgdp^4) +
            popgro +
            inv +
            humancap +
            I(humancap^2) +
            I(humancap^3) - 1, 
            x=TRUE, 
            y=TRUE)

X <- data.frame(factor(oecd), factor(year), initgdp, popgro, inv, humancap)

# Note - we override the default tolerances for the sake of this example
# (don't of course do this in general). This example may take a few
# minutes depending on the speed of your computer (data-driven bandwidth
# selection is, by its nature, time consuming, while the bootstrapping
# also takes some time).

npcmstest(model = model, xdat = X, ydat = growth, tol=.1, ftol=.1)

detach(oecdpanel)
} % enddontrun
}
\keyword{ nonparametric }
