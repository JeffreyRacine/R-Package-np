\name{mpi.any.source}
\alias{mpi.any.source}
\alias{mpi.any.tag}
\alias{mpi.proc.null}
\title{MPI Constants}
\description{
  Find MPI constants: MPI\_ANY\_SOURCE, MPI\_ANY\_TAG, or MPI\_PROC\_NULL
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
  integers for MPI\_ANY\_SOURCE, MPI\_ANY\_TAG, and MPI\_PROC\_NULL. Hence one 
  should use these functions instead  real integers for MPI communications.
} 
\seealso{
  \code{\link{mpi.send}}, \code{\link{mpi.recv}}.
}
\references{
  \url{http://www.lam-mpi.org}, \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\keyword{utilities}
