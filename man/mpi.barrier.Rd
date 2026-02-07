\name{mpi.barrier}
\title{MPI_Barrier API}
\usage{
  mpi.barrier(comm = 1)
}
\alias{mpi.barrier}
\arguments{
  \item{comm}{a communicator number}
}
\description{
  \code{mpi.barrier} blocks the caller until all members have called it.
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
%\examples{
%}
\keyword{interface}
 
