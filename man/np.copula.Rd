\name{npcopula}
\alias{npcopula}
\title{ Kernel Copula Estimation }
\description{
  \code{npcopula} implements the kernel copula approach 
  of Racine (2012) for an arbitrary number of dimensions
}
\usage{
npcopula(bws.joint,
         bws.univariate=NULL,
         u=NULL) 
}

\arguments{
  \item{bws.joint}{
    an unconditional joint bandwidth object
  }
  \item{bws.univariate}{
    an optional list of unconditional univariate bandwidth objects
  }
  \item{u}{
    a matrix of real numbers lying in [0,1], each column of which
    corresponds to the vector of uth quantile values desired (otherwise
    the u values returned are those corresponding to the sample realizations)
  }
}

\value{
  \code{npcopula} returns an object of type \code{dataframe} with the
  following components

  \item{copula}{ the copula }
  \item{copula density}{ the copula density}
  \item{u}{ the matrix of quantiles associated with the sample
  realizations (\code{u=NULL}) or those created via \code{expand.grid}
  when \code{u=u} is specified}

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

  \code{npcopula} computes the nonparametric copula and copula density
  using inversion (Nelsen (2006), page 51). For the inversion approach,
  we exploit Sklar's theorem (Corollary 2.3.7, Nelsen (2006)) to produce
  copulas directly from the joint distribution function using
  \eqn{C(u,v) = H(F^{-1}(u),G^{-1}(v))} rather than the typical approach
  that instead uses \eqn{H(x,y) = C(F(x),G(y))}. Whereas the latter
  requires kernel density estimation on a d-dimensional unit hypercube
  which necessitates the use of boundary correction methods, the former
  does not (unless of course the data itself is bounded requiring the
  use of boundary corrections).

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
data <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=data[,2],y=data[,1])
attach(mydat)

bw <- npudensbw(~x+y,nmulti=1,tol=.1,ftol=.1)

q.min <- 0.0
q.max <- 1.0
grid.seq <- seq(q.min,q.max,length=n.eval)
grid.dat <- cbind(grid.seq,grid.seq)

copula <- npcopula(bws.joint=bw,u=grid.dat)
contour(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),xlab="u1",ylab="u2")
persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),ticktype="detailed",xlab="u1",ylab="u2",zlab="Copula",zlim=c(0,1))
persp(grid.seq,grid.seq,matrix(copula$copula.density,n.eval,n.eval),ticktype="detailed",xlab="u1",ylab="u2",zlab="Copula Density")
copula.emp <- npcopula(bws.joint=bw)
plot(copula.emp$u1,copula.emp$u2,xlab="u1",ylab="u2",cex=.25)
}

\keyword{ nonparametric }
