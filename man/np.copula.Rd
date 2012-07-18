\name{npcopula}
\alias{npcopula}
\title{ Kernel Copula Estimation }
\description{
  \code{npcopula} implements the kernel copula approach 
  of Racine (2012) for an arbitrary number of dimensions
}
\usage{
npcopula(bws,
         data,
         u = NULL,
         density = FALSE,
         n.pseudo.inv = 500) 
}

\arguments{
  \item{bws}{
    an unconditional joint distribution/density bandwidth object
  }
  \item{data}{
    a dataframe containing variables used to construct \code{bws}
  }
  \item{u}{
    a matrix of real numbers lying in [0,1], each column of which
    corresponds to the vector of uth quantile values desired (otherwise
    the u values returned are those corresponding to the sample realizations)
  }
  \item{density}{
    a logical value indicating whether to compute and return the copula
    or copula density (default is \code{FALSE})
  }
  \item{n.pseudo.inv}{
    a value used to contruct a grid of length \eqn{2\times
    n.pseudo.inv}{2 x n.pseudo.inv} that is used to compute the
    quantiles given \code{u} via pseudo-inverse of the marginal
    distributions
  }
}

\value{

  \code{npcopula} returns an object of type \code{dataframe} with the
  following components

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

  Racine, J.S. (2012), \dQuote{Kernel Copula Estimation Via Inversion,}
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
  does not (unless of course the data itself is bounded requiring the
  use of boundary corrections).

  Note that if \code{u} is provided then \code{\link{expand.grid}} is
  called on \code{u}. As the dimension increases this can become
  unwieldy and potentially consume an enormous amount of memory unless
  the number of grid points is kept very small. As the reason for
  computing the copula on a grid is typically for graphical purposes,
  then providing \code{u} is typically done for two dimension problems
  only. Even here, however, providing a grid of length 100 will expand
  into a matrix of dimension 10000 by 2 which, though not memory
  intensive, may be computationally burdensome.

  Note that copula are only defined for data of type
  \code{\link{numeric}} or \code{\link{ordered}}.

  The pseudo-inverse (\sQuote{quasi-inverse}) is computed via Definition
  2.3.6 from Nelson (2006).

}

\section{Usage Issues}{

  See the example below for proper usage.

}

\seealso{
  \link{npudensbw},\link{npudens},\link{npudist}
}

\examples{
require(MASS)

set.seed(42)

n <- 100
n.eval <- 25
rho <- 0.95
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
mydat <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=mydat[,1],y=mydat[,2])
bw <- npudensbw(~x+y,data=mydat)

q.min <- 0.0
q.max <- 1.0
grid.seq <- seq(q.min,q.max,length=n.eval)
grid.dat <- cbind(grid.seq,grid.seq)

## Here we use the same bandwidths for the copula and copula density
## (temporary awaiting CDF bandwidth selector)

copula <- npcopula(bws=bw,data=mydat,u=grid.dat)
contour(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),xlab="u1",ylab="u2")
persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),ticktype="detailed",xlab="u1",ylab="u2",
      zlab="Copula",zlim=c(0,1))
copula <- npcopula(bws=bw,data=mydat,u=grid.dat,density=TRUE)
persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),ticktype="detailed",xlab="u1",
      ylab="u2",zlab="Copula Density")
copula.emp <- npcopula(bws=bw,data=mydat)
plot(copula.emp$u1,copula.emp$u2,xlab="u1",ylab="u2",cex=.25)
}

\keyword{ nonparametric }
