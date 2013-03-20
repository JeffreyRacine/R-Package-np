\name{mpi.cart.coords}
\title{MPI\_Cart\_coords}
\alias{mpi.cart.coords}
\description{
 \code{mpi.cart.coords} translates a rank to its Cartesian topology coordinate.
}
\usage{mpi.cart.coords(comm=3, rank, maxdims)
}
\arguments{
 \item{comm}{Communicator with Cartesian structure}
 \item{rank}{rank of a process within group}
 \item{maxdims}{length of vector coord in the calling program}
}
\value{
 \code{mpi.cart.coords} returns an integer array containing the Cartesian
 coordinates of specified process. 
}
\details{
 This function is the rank-to-coordinates translator.  It is the inverse map of
 \code{mpi.cart.rank}.  maxdims is at least as big as ndims as
 returned by \code{mpi.cartdim.get}.
}
\author{
Alek Hunchak and Hao Yu
}
\seealso{
\code{\link{mpi.cart.rank}}
}
\references{
  \url{http://www.lam-mpi.org/},
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\examples{
\dontrun{
#Need at least 9 slaves
mpi.bcast.cmd(mpi.cart.create(1,c(3,3),c(F,T)))
mpi.cart.create(1,c(3,3),c(F,T))
mpi.cart.coords(3,4,2)
}

} 
\keyword{interface}
