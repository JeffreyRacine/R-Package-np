\name{mpi.cartdim.get}
\title{MPI\_Cartdim\_get}
\alias{mpi.cartdim.get}
\description{
 \code{mpi.cartdim.get} gets dim information about a Cartesian topology. 
}
\usage{
 mpi.cartdim.get(comm=3)
}
\arguments{
 \item{comm}{Communicator with Cartesian structure}
}
\value{
 \code{mpi.cartdim.get} returns the number of dimensions of the Cartesian structure
}
\details{
 Can be used to provide other functions with the correct size of arrays.
}
\author{
Alek Hunchak and Hao Yu
}
\seealso{
\code{\link{mpi.cart.get}}
}
\references{
  \url{http://www.lam-mpi.org/},
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\examples{
\donttest{
#Need at least 9 slaves
mpi.bcast.cmd(mpi.cart.create(1,c(3,3),c(F,T)))
mpi.cart.create(1,c(3,3),c(F,T))
mpi.cartdim.get(comm=3)
}
}
\keyword{interface} 
