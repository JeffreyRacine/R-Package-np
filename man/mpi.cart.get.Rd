\name{mpi.cart.get}
\title{MPI_Cart_get}
\alias{mpi.cart.get}
\description{
 \code{mpi.cart.get} provides the user with information on the Cartesian topology
 associated with a comm.
}
\usage{
 mpi.cart.get(comm=3, maxdims)
}
\arguments{
 \item{comm}{Communicator with Cartesian structure}
 \item{maxdims}{length of vectors dims, periods, and coords in the calling program}
}
\value{
 \code{mpi.cart.get} returns a vector containing information on the Cartesian topology 
 associated with comm.  maxdims must be at least ndims as returned by \code{mpi.cartdim.get}.  
}
\details{
 The coords are as given for the rank of the calling process as shown.
}
\author{
Alek Hunchak and Hao Yu
}
\seealso{
\code{\link{mpi.cart.create},\link{mpi.cartdim.get}}
}
\references{
  \url{https://www.mpich.org/},
  \url{https://www.mpich.org/static/docs/latest/www3/}
}
\examples{
\dontrun{
# Not run in checks: requires a Cartesian communicator built from spawned slaves.
#Need at least 9 slaves
mpi.bcast.cmd(mpi.cart.create(1,c(3,3),c(F,T)))
mpi.cart.create(1,c(3,3),c(F,T))
mpi.remote.exec(mpi.cart.get(3,2))
}
}
 
\keyword{interface}
