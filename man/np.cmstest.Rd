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
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  \code{npcmstest} supports regression objects generated by
  \code{\link{lm}} and uses features specific to objects of type
  \code{\link{lm}} hence if you attempt to pass objects of a different
  type the function cannot be expected to work.

  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.
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

mpi.bcast.cmd(data(cps71),
              caller.execute=TRUE)

mpi.bcast.cmd(attach(cps71),
              caller.execute=TRUE)

mpi.bcast.cmd(model <- lm(logwage~age+I(age^2), x=TRUE, y=TRUE),
              caller.execute=TRUE)

mpi.bcast.cmd(npcmstest(model = model, xdat = age, ydat = logwage,
              boot.num = 99),
              caller.execute=TRUE)

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
