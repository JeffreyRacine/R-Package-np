\name{mpi.barrier}
\title{MPI\_Barrier API}
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
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
%\examples{
%}
\keyword{interface}
 
