\name{mpi.universe.size}
\title{MPI_Universe_size API}
\usage{
  mpi.universe.size()
}
\alias{mpi.universe.size}
\arguments{
  None.
}
\description{
  \code{mpi.universe.size} returns the total number of CPUs available in a 
cluster. Some MPI implements may not have this MPI call available.
}
\value{
An integer giving the total number of CPUs available in the MPI universe
for the current configuration.
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
 
