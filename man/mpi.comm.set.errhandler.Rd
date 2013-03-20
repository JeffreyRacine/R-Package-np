\name{mpi.comm.set.errhandler}
\title{MPI\_Comm\_set\_errhandler API}
\usage{
mpi.comm.set.errhandler(comm = 1)
}
\alias{mpi.comm.set.errhandler}
\arguments{
  \item{comm}{a communicator number}
}
\description{
  \code{mpi.comm.set.errhandler} sets a communicator to MPI\_ERRORS\_RETURN 
instead of  
MPI\_ERRORS\_ARE\_FATAL (default) which crashes R on any type of MPI errors.  
Almost all MPI API calls return errcodes which can map to specific MPI error 
messages. All MPI related error messages come from predefined 
MPI\_Error\_string. 
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
