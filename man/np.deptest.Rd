% $Id: np.cmstest.Rd,v 1.58 2006/11/03 21:17:20 tristen Exp $
\name{npdeptest}
\alias{npdeptest}
\title{ Kernel Consistent Pairwise Nonlinear Dependence Test for Univariate Processes }
\description{
  
  \code{npdeptest} implements the consistent metric entropy test of
  pairwise independence as described in Maasoumi and Racine (2002).
  
}
\usage{

npdeptest(data.x = NULL,
          data.y = NULL,
          method = c("integration","summation"),
          bootstrap = TRUE,
          boot.num = 399,
          random.seed = 42)
  
}

\arguments{
  \item{data.x, data.y}{
    
    two univariate vectors containing two variables that are of type
    \code{\link{numeric}}.
    
  }
  
  \item{method}{

    a character string used to specify whether to compute the integral
    version or the summation version of the statistic. Can be set as
    \code{integration} or \code{summation} (see below for
    details). Defaults to \code{integration}.

  }
  
  \item{bootstrap}{
    a logical value which specifies whether to conduct
    the bootstrap test or not. If set to \code{FALSE}, only the
    statistic will be computed. Defaults to \code{TRUE}.    
  }
  
  \item{boot.num}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }
  
  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }

}
            
\value{
  \code{npdeptest} returns an object of type \code{deptest} with the
  following components
  
  \item{Srho}{ the statistic \code{Srho} }
  \item{Srho.bootstrap.vec}{ contains the bootstrap replications of
    \code{Srho} }
  \item{P}{ the P-value of the Srho statistic }
  \item{bootstrap}{ a logical value indicating whether bootstrapping was
  performed }
  \item{boot.num}{ number of bootstrap replications }
  \item{bw.data.x}{ the numeric bandwidth for \code{data.x} marginal
    density}
  \item{bw.data.y}{ the numeric bandwidth for
    \code{data.y} marginal density}  
  \item{bw.joint}{ the numeric matrix of bandwidths for \code{data}
  and lagged \code{data} joint density  at lag \code{num.lag}}      

\code{\link{summary}} supports object of type \code{deptest}.
}
\references{
  
  Granger, C.W. and E. Maasoumi and J.S. Racine (2004), \dQuote{A
  dependence metric for possibly nonlinear processes}, Journal of Time
  Series Analysis, 25, 649-669.
  
  Maasoumi, E. and J.S. Racine (2002), \dQuote{Entropy and
    Predictability of Stock Market Returns,} Journal of Econometrics,
    107, 2, pp 291-312.

}
\author{
    Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{

  \code{npsdeptest} computes the nonparametric metric entropy
   (normalized Hellinger of Granger, Maasoumi and Racine (2004)) for
   testing pairwise nonlinear dependence between the densities of two
   data series. See Maasoumi and Racine (2002) for details. Default
   bandwidths are of the Kullback-Leibler variety obtained via
   likelihood cross-validation. The null distribution is obtained via
   bootstrap resampling under the null of pairwise independence.

  \code{npdeptest} computes the distance between the joint distribution
  and the product of marginals (i.e. the joint distribution under the
  null), \eqn{D[f(y, \hat y), f(y)\times f(\hat y)]}{D[f(y, y-hat),
  f(y)*f(y-hat)]}. Examples include, (a) a measure/test of \dQuote{fit},
  for in-sample values of a variable \eqn{y}{y} and its fitted values,
  \eqn{\hat y}{y-hat}, and (b) a measure of \dQuote{predictability} for
  a variable \eqn{y} and its predicted values \eqn{\hat y}{y-hat} (from
  a user implemented model).

  The summation version of this statistic will be numerically unstable
  when \code{data.x} and \code{data.y} lack common support or are sparse
  (the summation version involves division of densities while the
  integration version involves differences). Warning messages are
  produced should this occur (`integration recommended') and should be
  heeded.
 
}

\section{Usage Issues}{
  
  The \code{integration} version of the statistic uses multdimensional
  numerical methods from the \code{\link{cubature}} package. See
  \code{\link{adaptIntegrate}} for details. The \code{integration}
  version of the statistic will be substantially slower than the
  \code{summation} version, however, it will likely be both more
  accurate and powerful.
  
}

\seealso{
  \code{\link{npdeneqtest},\link{npsdeptest},\link{npsymtest},\link{npunitest}}
}

\examples{
set.seed(1234)

## Test/measure lack of fit between y and its fitted value from a
## regression model when x is relevant using the `summation' version.

n <- 100
x <- rnorm(n)
y <- 1 + x + rnorm(n)
model <- lm(y~x)
y.fit <- fitted(model)

npdeptest(y,y.fit,boot.num=99,method="summation")

Sys.sleep(5)

## Test/measure lack of fit between y and its fitted value from a
## regression model when x is irrelevant using the `summation' version.

n <- 100
x <- runif(n,-2,2)
y <- 1 + rnorm(n)
model <- lm(y~x)
y.fit <- fitted(model)

npdeptest(x,y,boot.num=99,method="summation")

\dontrun{
## Test/measure lack of fit between y and its fitted value from a
## regression model when x is relevant using the `integration'
## version (default, slower than summation version).

n <- 100
x <- rnorm(n)
y <- 1 + x + rnorm(n)
model <- lm(y~x)
y.fit <- fitted(model)

npdeptest(y,y.fit,boot.num=99)

Sys.sleep(5)

## Test/measure lack of fit between y and its fitted value from a
## regression model when x is irrelevant using the `integration'
## version (default, slower than summation version).

n <- 100
x <- runif(n,-2,2)
y <- 1 + rnorm(n)
model <- lm(y~x)
y.fit <- fitted(model)

npdeptest(x,y,boot.num=99)
} % enddontrun
}
\keyword{ univar }
\keyword{ nonparametric }
