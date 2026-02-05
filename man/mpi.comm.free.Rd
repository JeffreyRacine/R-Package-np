\name{mpi.comm.free}
\title{MPI\_Comm\_free API}
\usage{
  mpi.comm.free(comm=1)
}
\alias{mpi.comm.free}
\arguments{
  \item{comm}{a communicator number}
}
\description{
  \code{mpi.comm.free}  deallocates a communicator so it 
  points to MPI\_COMM\_NULL.
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
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
 \code{\link{mpi.comm.disconnect}}
}
%\examples{
%}
\keyword{interface}
 
