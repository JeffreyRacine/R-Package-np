\name{mpi.comm.free}
\title{MPI_Comm_free API}
\usage{
  mpi.comm.free(comm=1)
}
\alias{mpi.comm.free}
\arguments{
  \item{comm}{a communicator number}
}
\description{
  \code{mpi.comm.free}  deallocates a communicator so it 
  points to MPI_COMM_NULL.
}
\details{
  When members associated with a communicator finish jobs or exit, they have to 
  call \code{mpi.comm.free} to release resource so \code{\link{mpi.comm.size}} 
  will return 0. If the comm  was created from an intercommunicator by 
  \code{\link{mpi.intercomm.merge}}, use \code{\link{mpi.comm.disconnect}} instead.
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
 \code{\link{mpi.comm.disconnect}}
}
%\examples{
%}
\keyword{interface}
 
