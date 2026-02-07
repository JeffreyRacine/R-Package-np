\name{mpi.abort}
\title{MPI_Abort API}
\usage{
  mpi.abort(comm = 1)
}
\alias{mpi.abort}
\arguments{
  \item{comm}{a communicator number}
}
\description{
  \code{mpi.abort} makes a ``best attempt" to abort all tasks in a comm.
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
 \code{\link{mpi.finalize}}
}
%\examples{
%}
\keyword{interface}
 
