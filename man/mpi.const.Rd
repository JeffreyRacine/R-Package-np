\name{mpi.any.source}
\alias{mpi.any.source}
\alias{mpi.any.tag}
\alias{mpi.proc.null}
\title{MPI Constants}
\description{
  Find MPI constants: MPI_ANY_SOURCE, MPI_ANY_TAG, or MPI_PROC_NULL
}
\usage{
mpi.any.source()
mpi.any.tag()
mpi.proc.null()
}
\arguments{
  None
}
\value{
  Each function returns an integer value.
}
\details{
  These constants are mainly used by 
\code{\link{mpi.send}}, \code{\link{mpi.recv}}, and 
\code{\link{mpi.probe}}. 
Different implementation of MPI may use different 
  integers for MPI_ANY_SOURCE, MPI_ANY_TAG, and MPI_PROC_NULL. Hence one 
  should use these functions instead of real integers for MPI communications.
} 
\seealso{
  \code{\link{mpi.send}}, \code{\link{mpi.recv}}.
}
\references{
  \url{http://www.lam-mpi.org}, \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\keyword{utilities}
