% $Id: np.cmstest.Rd,v 1.58 2006/11/03 21:17:20 tristen Exp $
\name{npdeneqtest}
\alias{npdeneqtest}
\title{ Kernel Consistent Density Equality Test with Mixed Data Types }
\description{
  \code{npdeneqtest} implements a consistent integrated squared
  difference test for equality of densities as described in Li, Maasoumi,
  and Racine (2009).
}
\usage{
npdeneqtest(x = NULL,
            y = NULL,
            bw.x = NULL,
            bw.y = NULL,
            boot.num = 399,
            random.seed = 42,
            \dots)
}

\arguments{
  \item{x,y}{
    data frames for the two samples for which one wishes to
    test equality of densities. The variables in each data
    frame must be the same (i.e. have identical names).
  }
  
  \item{bw.x,bw.y}{
    optional bandwidth objects for \code{x,y}
  }
  
  \item{boot.num}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }
  
  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }

  \item{\dots}{ additional arguments supplied to specify the bandwidth
    type, kernel types, and so on.  This is used if you do not pass in
    bandwidth objects and you do not desire the default behaviours. To
    do this, you may specify any of \code{bwscaling}, \code{bwtype},
    \code{ckertype}, \code{ckerorder}, \code{ukertype},
    \code{okertype}.}
  
}

\value{
  \code{npdeneqtest} returns an object of type \code{deneqtest} with the
  following components
  
  \item{Tn}{ the (standardized) statistic \code{Tn} }
  \item{In}{ the (unstandardized) statistic \code{In} }
  \item{Tn.bootstrap}{ contains the bootstrap replications of \code{Tn} }
  \item{In.bootstrap}{ contains the bootstrap replications of \code{In} }
  \item{Tn.P}{ the P-value of the \code{Tn} statistic }
  \item{In.P}{ the P-value of the \code{In} statistic }
  \item{boot.num}{ number of bootstrap replications }
    
  \code{\link{summary}} supports object of type \code{deneqtest}.
}
\references{
  Li, Q. and E. Maasoumi and J.S. Racine (2009), \dQuote{A Nonparametric
    Test for Equality of Distributions with Mixed Categorical and
    Continuous Data,} Journal of Econometrics, 148, pp 186-200.
}
\author{
    Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{
  
  \code{npdeneqtest} computes the integrated squared density difference
  between the estimated densities/probabilities of two samples having
  identical variables/datatypes. See Li, Maasoumi, and Racine (2009) for
  details.
 
}

\section{Usage Issues}{
  
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.

  It is crucial that both data frames have the same variable names.

}

\seealso{
  \code{\link{npdeptest},\link{npsdeptest},\link{npsymtest},\link{npunitest}}
}

\examples{
set.seed(1234)

## Distributions are equal

n <- 250

sample.A <- data.frame(x=rnorm(n))
sample.B <- data.frame(x=rnorm(n))

npdeneqtest(sample.A,sample.B,boot.num=99)

Sys.sleep(5)

## Distributions are unequal

sample.A <- data.frame(x=rnorm(n))
sample.B <- data.frame(x=rchisq(n,df=5))

npdeneqtest(sample.A,sample.B,boot.num=99)

\dontrun{
## Mixed datatypes, distributions are equal

sample.A <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))
sample.B <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))

npdeneqtest(sample.A,sample.B,boot.num=99)

Sys.sleep(5)

## Mixed datatypes, distributions are unequal

sample.A <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))
sample.B <- data.frame(a=rnorm(n,sd=10),b=factor(rbinom(n,2,.25)))

npdeneqtest(sample.A,sample.B,boot.num=99)
} % enddontrun
}

\keyword{ nonparametric }
