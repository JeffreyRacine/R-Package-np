\name{mpi.cart.shift}
\title{MPI\_Cart\_shift}
\alias{mpi.cart.shift}
\description{
 \code{mpi.cart.shift} shifts the Cartesian topology in both manners, displacement
 and direction.
}
\usage{
 mpi.cart.shift(comm=3, direction, disp)
}
\arguments{
 \item{comm}{Communicator with Cartesian structure}
 \item{direction}{Coordinate dimension of the shift}
 \item{disp}{displacement (>0 for upwards or left shift, <0 for downwards or right shift)}
}
\value{
 \code{mpi.cart.shift} returns a vector containing information regarding the 
 rank of the source process and rank of the destination process.
}
\details{
 \code{mpi.cart.shift} provides neighbor ranks from given direction and displacement.
 The direction argument indicates the dimension of the shift. direction=1 means the first dim, 
 direction=2 means the second dim, etc. disp=1 or -1 provides immediate neighbor ranks and disp=2 
 or -2 provides neighbor's neighbor ranks. Negative ranks mean out of boundary. They correspond to 
 \code{mpi.proc.null}.  
}
\author{
Alek Hunchak and Hao Yu
}
\seealso{
\code{\link{mpi.cart.create}},\code{\link{mpi.proc.null}}
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
mpi.remote.exec(mpi.cart.shift(3,2,1))#get neighbor ranks
mpi.remote.exec(mpi.cart.shift(3,1,1))
}
} 
\keyword{interface}
