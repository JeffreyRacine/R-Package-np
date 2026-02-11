\name{mpi.cart.rank}
\title{MPI_Cart_rank}
\alias{mpi.cart.rank}
\description{
 \code{mpi.cart.rank} translates a Cartesian topology coordinate to its rank.
}
\usage{
 mpi.cart.rank(comm=3, coords)
}
\arguments{
 \item{comm}{Communicator with Cartesian structure}
 \item{coords}{Specifies the Cartesian coordinates of a process}
}
\value{
 \code{mpi.cart.rank} returns the rank of the specified process.
}
\details{
 For a process group with a Cartesian topology, this function translates the logical
 process coordinates to process ranks as they are used by the point-to-point routines.
 It is the inverse map of \code{mpi.cart.coords}. }
\author{
Alek Hunchak and Hao Yu
}
\seealso{
\code{\link{mpi.cart.coords}}
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
mpi.cart.rank(3,c(1,0))
}
} 
\keyword{interface}
