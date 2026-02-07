\name{mpi.comm.disconnect}
\title{MPI_Comm_disconnect API}
\usage{
mpi.comm.disconnect(comm=1)
}
\alias{mpi.comm.disconnect}
\arguments{
  \item{comm}{a communicator number}
}
\description{
  \code{mpi.comm.disconnect} disconnects itself from a communicator and then 
  deallocates the communicator so it points to MPI_COMM_NULL.
}
\details{
  When members associated with a communicator finish jobs or exit, they have to 
  call \code{mpi.comm.disconnect} to release resource if the communicator was 
  created from an intercommunicator by \code{\link{mpi.intercomm.merge}}. If 
  \code{\link{mpi.comm.free}} is used instead, \code{\link{mpi.finalize}} called 
  by slaves may cause undefined impacts on master who wishes to stay.
}
\value{
  1 if success. Otherwise 0. 
}
\author{
Hao Yu
}
\references{
  \url{https://www.mpich.org/}, 
  \url{https://www.mpich.org/static/docs/latest/www3/}
}
\seealso{
 \code{\link{mpi.comm.free}}
}
%\examples{
%}
\keyword{interface}
