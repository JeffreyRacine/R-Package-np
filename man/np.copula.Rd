\name{npcopula}
\alias{npcopula}
\title{ Kernel Copula Estimation with Mixed Data Types }

\description{
  \code{npcopula} implements the nonparametric mixed data kernel copula
  approach of Racine (2012) for an arbitrary number of dimensions
}

\usage{
npcopula(bws,
         data,
         u = NULL,
         n.quasi.inv = 1000,
         er.quasi.inv = 1) 
}

\arguments{
  \item{bws}{
    an unconditional joint distribution (\code{npudistbw}) or joint
    density (\code{npudensbw}) bandwidth object (if \code{bws} is
    delivered via \code{npudistbw} the copula distribution is estimated,
    while if \code{bws} is delivered via \code{npudensbw} the copula
    density is estimated)
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

  \item{copula}{
    the copula (bandwidth object obtained from \code{npudistbw}) or
    copula density (bandwidth object obtained from \code{npudensbw})
  }

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
  the number of grid points is kept very small. Given that computing the
  copula on a grid is typically done for graphical purposes, providing
  \code{u} is typically done for two-dimensional problems only. Even
  here, however, providing a grid of length 100 will expand into a
  matrix of dimension 10000 by 2 which, though not memory intensive, may
  be computationally burdensome.

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

  Note that if \code{u} is provided and any elements of (the columns of)
  \code{u} are such that they lie beyond the respective values of
  \code{F} for the evaluation data for the respective marginal, such
  values are reset to the minimum/maximum values of \code{F} for the
  respective marginal. It is therefore prudent to inspect the values of
  \code{u} returned by \code{\link{npcopula}} when \code{u} is provided.

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
## Example 1: Bivariate Mixed Data

## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

mpi.spawn.Rslaves(nslaves=1)
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

mpi.bcast.cmd(library(MASS),caller.execute=TRUE)

set.seed(42)

## Simulate correlated Gaussian data (rho(x,y)=0.99)

n <- 1000
n.eval <- 100
rho <- 0.99
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
mydat <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=mydat[,1],
                    y=ordered(as.integer(cut(mydat[,2],
                      quantile(mydat[,2],seq(0,1,by=.1)),
                      include.lowest=TRUE))-1))
q.min <- 0.0
q.max <- 1.0
grid.seq <- seq(q.min,q.max,length=n.eval)
grid.dat <- cbind(grid.seq,grid.seq)

mpi.bcast.Robj2slave(mydat)
mpi.bcast.Robj2slave(grid.dat)

## Estimate the copula (bw object obtained from npudistbw())

mpi.bcast.cmd(bw.cdf <- npudistbw(~x+y,data=mydat),
              caller.execute=TRUE)
mpi.bcast.cmd(copula <- npcopula(bws=bw.cdf,data=mydat,u=grid.dat),
              caller.execute=TRUE)

## Plot the copula


contour(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
        xlab="u1",
        ylab="u2",
        main="Copula Contour")

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
      ticktype="detailed",
      xlab="u1",
      ylab="u2",
      zlab="Copula",zlim=c(0,1))

## Plot the empirical copula

mpi.bcast.cmd(copula.emp <- npcopula(bws=bw.cdf,data=mydat),
              caller.execute=TRUE)
plot(copula.emp$u1,copula.emp$u2,xlab="u1",ylab="u2",cex=.25,main="Empirical Copula")

## Estimate the copula density (bw object obtained from npudensbw())

mpi.bcast.cmd(bw.pdf <- npudensbw(~x+y,data=mydat),
              caller.execute=TRUE)
mpi.bcast.cmd(copula <- npcopula(bws=bw.pdf,data=mydat,u=grid.dat),
              caller.execute=TRUE)

## Plot the copula density

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
      ticktype="detailed",
      xlab="u1",
      ylab="u2",
      zlab="Copula Density")

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

mpi.close.Rslaves()

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)

## Example 2: Bivariate Continuous Data

## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

mpi.spawn.Rslaves(nslaves=1)
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

mpi.bcast.cmd(library(MASS),caller.execute=TRUE)

set.seed(42)

## Simulate correlated Gaussian data (rho(x,y)=0.99)

n <- 1000
n.eval <- 100
rho <- 0.99
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
mydat <- mvrnorm(n=n, mu, Sigma)
mydat <- data.frame(x=mydat[,1],y=mydat[,2])

q.min <- 0.0
q.max <- 1.0
grid.seq <- seq(q.min,q.max,length=n.eval)
grid.dat <- cbind(grid.seq,grid.seq)

mpi.bcast.Robj2slave(mydat)
mpi.bcast.Robj2slave(grid.dat)

## Estimate the copula (bw object obtained from npudistbw())

mpi.bcast.cmd(bw.cdf <- npudistbw(~x+y,data=mydat),
              caller.execute=TRUE)
mpi.bcast.cmd(copula <- npcopula(bws=bw.cdf,data=mydat,u=grid.dat),
              caller.execute=TRUE)

## Plot the copula

contour(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
        xlab="u1",
        ylab="u2",
        main="Copula Contour")

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
      ticktype="detailed",
      xlab="u1",
      ylab="u2",
      zlab="Copula",
      zlim=c(0,1))

## Plot the empirical copula

mpi.bcast.cmd(copula.emp <- npcopula(bws=bw.cdf,data=mydat),
              caller.execute=TRUE)
plot(copula.emp$u1,copula.emp$u2,xlab="u1",ylab="u2",cex=.25,main="Empirical Copula")

## Estimate the copula density (bw object obtained from npudensbw())

mpi.bcast.cmd(bw.pdf <- npudensbw(~x+y,data=mydat),
              caller.execute=TRUE)
mpi.bcast.cmd(copula <- npcopula(bws=bw.pdf,data=mydat,u=grid.dat),
              caller.execute=TRUE)

## Plot the copula density

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),ticktype="detailed",xlab="u1",
      ylab="u2",zlab="Copula Density")

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
