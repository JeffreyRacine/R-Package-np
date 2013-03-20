\name{mpi.cart.create}
\title{MPI\_Cart\_create}
\alias{mpi.cart.create}
\description{
 \code{mpi.cart.create} creates a Cartesian structure of arbitrary dimension.
}
\usage{
 mpi.cart.create(commold=1, dims, periods, reorder=FALSE, commcart=3)
}
\arguments{
 \item{commold}{Input communicator}
 \item{dims}{Integery array of size ndims specifying the number of processes in 
 each dimension}
 \item{periods}{Logical array of size ndims specifying whether the grid is periodic
 or not in each dimension}
 \item{reorder}{ranks may be reordered or not}
 \item{commcart}{The new communicator to which the Cartesian topology information is attached}
}
\value{
 \code{mpi.cart.create} returns 1 if success and 0 otherwise.
}
\details{
 If reorder = false, then the rank of each process in the new group is the same
 as its rank in the old group.  If the total size of the Cartesian grid is smaller
 than the size of the group of commold, then some processes are returned 
 mpi.comm.null.  The call is erroneous if it specifies a grid that is larger than
 the group size.
}
\author{
Alek Hunchak and Hao Yu
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
}
}
 
\keyword{interface}
