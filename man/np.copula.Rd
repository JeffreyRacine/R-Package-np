\name{npcopula}
\alias{npcopula}
\title{ Kernel Copula Estimation with Mixed Data Types }

\description{
  \code{npcopula} implements the nonparametric mixed data kernel copula
  approach of Racine (2012) for an arbitrary number of dimensions
}

\usage{
npcopula(bws,
         bws.pdf,
         data,
         u = NULL,
         density = FALSE,
         n.quasi.inv = 1000,
         er.quasi.inv = 1) 
}

\arguments{
  \item{bws}{
    an unconditional joint distribution bandwidth object
  }
  \item{bws.pdf}{
    an unconditional joint density bandwidth object (only needed when \code{density=TRUE})
  }
  \item{data}{
    a data frame containing variables used to construct \code{bws}
  }
  \item{u}{
    an optional matrix of real numbers lying in [0,1], each column of
    which corresponds to the vector of uth quantile values desired for
    each variable in the copula (otherwise the u values returned are
    those corresponding to the sample realizations)
  }
  \item{density}{
    a logical value indicating whether to compute and return the copula
    or copula density (default is \code{FALSE} [copula])
  }
  \item{n.quasi.inv}{
    number of grid points generated when \code{u} is provided in order to
    compute the quasi-inverse of each marginal distribution (see
    details)
  }
  \item{er.quasi.inv}{
    number passed to \code{\link{extendrange}} when \code{u} is provided
    specifying the fraction by which the data range should be extended when
    constructing the grid used to compute the quasi-inverse (see details)
  }
}

\value{

  \code{npcopula} returns an object of type \code{\link{data.frame}}
  with the following components

  \item{copula}{ the copula (default) or copula density
  (\code{density=TRUE}) }

  \item{u}{ the matrix of marginal u values associated with the sample
  realizations (\code{u=NULL}) or those created via
  \code{\link{expand.grid}} when \code{u} is provided}

  \item{data}{ the matrix of marginal quantiles constructed when
  \code{u} is provided (\code{data} returned has the same names as
  \code{data} inputted)}

}
\references{

  Nelsen, R. B. (2006), \emph{An Introduction to Copulas,} Second
  Edition, Springer-Verlag.

  Racine, J.S. (2012), \dQuote{Mixed Data Kernel Copulas,}
  manuscript.

}
\author{
    Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\details{

  \code{npcopula} computes the nonparametric copula or copula density
  using inversion (Nelsen (2006), page 51). For the inversion approach,
  we exploit Sklar's theorem (Corollary 2.3.7, Nelsen (2006)) to produce
  copulas directly from the joint distribution function using
  \eqn{C(u,v) = H(F^{-1}(u),G^{-1}(v))} rather than the typical approach
  that instead uses \eqn{H(x,y) = C(F(x),G(y))}. Whereas the latter
  requires kernel density estimation on a d-dimensional unit hypercube
  which necessitates the use of boundary correction methods, the former
  does not.

  Note that if \code{u} is provided then \code{\link{expand.grid}} is
  called on \code{u}. As the dimension increases this can become
  unwieldy and potentially consume an enormous amount of memory unless
  the number of grid points is kept very small. As the reason for
  computing the copula on a grid is typically for graphical purposes,
  then providing \code{u} is typically done for two-dimensional problems
  only. Even here, however, providing a grid of length 100 will expand
  into a matrix of dimension 10000 by 2 which, though not memory
  intensive, may be computationally burdensome.

  The \sQuote{quasi-inverse} is computed via Definition 2.3.6 from
  Nelsen (2006). We compute an equi-quantile grid on the range of the
  data of length \code{n.quasi.inv/2}. We then extend the range of the
  data by the factor \code{er.quasi.inv} and compute an equi-spaced grid
  of points of length \code{n.quasi.inv/2} (e.g. using the default
  \code{er.quasi.inv=1} we go from the minimum data value minus
  \eqn{1\times}{1x} the range to the maximum data value plus
  \eqn{1\times}{1x} the range for each marginal). We then take these two
  grids, concatenate and sort, and these form the final grid of length
  \code{n.quasi.inv} for computing the quasi-inverse.

  Note that \code{npcopula}, when called with the option
  \code{density=TRUE}, allows for two levels of smoothing if desired,
  one for the distribution, \eqn{F(x,y)}, and one for the density,
  \eqn{f(x,y)}. These can be different or the same, for instance, that
  for the density only (we include this flexibility for extensibility
  purposes). If the user wishes to use the same bandwidths throughout
  then simply copy the bandwidths from, e.g. the \code{npudensbw}
  bandwidth object, \code{bw.pdf}, into those from, e.g the
  \code{npudistbw} bandwidth object, \code{bw.cdf}, as follows:

  \preformatted{

    bw.cdf$bws$bandwidth$x <- bw.pdf$bws$bandwidth$x
    bw.cdf$bw <- bw.pdf$bw

  }

  When these two bandwidth objects are passed to \code{npcopula} the
  same bandwidths will be used everywhere (i.e when inverting
  \eqn{F(x,y)} and when computing
  \eqn{f(F_x^{-1}(u_x),F_y^{-1}(u_y))/(f_x(F_x^{-1}(u_x))f_x(F_y^{-1}(u_y)))}).

  Note that copula are only defined for data of type
  \code{\link{numeric}} or \code{\link{ordered}}.

}
\section{Usage Issues}{

  See the example below for proper usage.

}

\seealso{
  \link{npudensbw},\link{npudens},\link{npudist}
}

\examples{
\dontrun{
require(MASS)

set.seed(42)

n <- 100
n.eval <- 25
rho <- 0.95
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
mydat <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=mydat[,1],y=mydat[,2])

q.min <- 0.0
q.max <- 1.0
grid.seq <- seq(q.min,q.max,length=n.eval)
grid.dat <- cbind(grid.seq,grid.seq)

## Estimate the copula

bw <- npudistbw(~x+y,data=mydat)
copula <- npcopula(bws=bw,data=mydat,u=grid.dat)

## Plot the copula

contour(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),xlab="u1",ylab="u2")
persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),ticktype="detailed",xlab="u1",ylab="u2",
      zlab="Copula",zlim=c(0,1))

## Estimate the copula density

bw.pdf <- npudensbw(~x+y,data=mydat)
copula <- npcopula(bws=bw,bws.pdf=bw.pdf,data=mydat,u=grid.dat,density=TRUE)

## Plot the copula density

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),ticktype="detailed",xlab="u1",
      ylab="u2",zlab="Copula Density")

## Plot the empirical copula

copula.emp <- npcopula(bws=bw,data=mydat)
plot(copula.emp$u1,copula.emp$u2,xlab="u1",ylab="u2",cex=.25)
} % enddontrun
}

\keyword{ nonparametric }
