\name{mpi.comm.size}
\title{MPI_Comm_c2f, MPI_Comm_dup, MPI_Comm_rank, and MPI_Comm_size APIs}
\usage{
  mpi.comm.c2f(comm=1)
  mpi.comm.dup(comm, newcomm)
  mpi.comm.rank(comm = 1)
  mpi.comm.size(comm = 1)
}
\alias{mpi.comm.c2f}
\alias{mpi.comm.dup}
\alias{mpi.comm.rank}
\alias{mpi.comm.size}
\arguments{
  \item{comm}{a communicator number}
  \item{newcomm}{a new communicator number}
}
\description{
\code{mpi.comm.c2f} converts the comm (a C communicator) and returns an integer that can be
used as the communicator in external FORTRAN code. \code{mpi.comm.dup} duplicates 
(copies) a comm to a new comm. \code{mpi.comm.rank} 
returns its rank in a comm. \code{mpi.comm.size} returns 
the total number of members in a comm. 
}
\value{
\itemize{
  \item \code{mpi.comm.c2f}: integer communicator for use in FORTRAN code.
  \item \code{mpi.comm.dup}: integer identifier of the duplicated communicator.
  \item \code{mpi.comm.rank}: integer rank within the communicator.
  \item \code{mpi.comm.size}: integer size of the communicator.
}
}
\author{
Hao Yu
}
\references{
  \url{https://www.mpich.org/}, 
  \url{https://www.mpich.org/static/docs/latest/www3/}
}
\examples{
\dontrun{
## Not run in checks when toggled to dontrun: communicator examples are
## documented for manual MPI sessions.
mpi.comm.rank(comm=0)
mpi.comm.size(comm=0)
mpi.comm.dup(comm=0, newcomm=5)
} % enddontrun
}
\keyword{interface}
 
