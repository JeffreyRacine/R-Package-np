\name{mpi.comm.set.errhandler}
\title{MPI_Comm_set_errhandler API}
\usage{
mpi.comm.set.errhandler(comm = 1)
}
\alias{mpi.comm.set.errhandler}
\arguments{
  \item{comm}{a communicator number}
}
\description{
  \code{mpi.comm.set.errhandler} sets a communicator to MPI_ERRORS_RETURN 
instead of  
MPI_ERRORS_ARE_FATAL (default) which crashes R on any type of MPI errors.  
Almost all MPI API calls return errcodes which can map to specific MPI error 
messages. All MPI related error messages come from predefined 
MPI_Error_string. 
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
