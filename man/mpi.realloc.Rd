\name{mpi.realloc}
\title{Find and increase the lengthes of MPI opaques comm, request, and 
status}
\alias{mpi.comm.maxsize}
\alias{mpi.request.maxsize}
\alias{mpi.status.maxsize}
\alias{mpi.realloc.comm}
\alias{mpi.realloc.request}
\alias{mpi.realloc.status}
\description{
  \code{mpi.comm.maxsize}, \code{mpi.request.maxsize}, and 
\code{mpi.status.maxsize} find the lengthes of comm, request, and status 
arrayes respectively.

  \code{mpi.realloc.comm}, \code{mpi.realloc.request} and 
\code{mpi.realloc.status} increase the lengthes of comm, request and 
status arrayes to \code{newmaxsize} respectively if \code{newmaxsize} is 
bigger than the original maximum size. 
}
\usage{
mpi.realloc.comm(newmaxsize)
mpi.realloc.request(newmaxsize)
mpi.realloc.status(newmaxsize)
mpi.comm.maxsize()
mpi.request.maxsize()
mpi.status.maxsize()
}
\arguments{
  \item{newmaxsize}{an integer.}
}
\details{
When \pkg{Rmpi} is loaded, Rmpi allocs comm array with size 10, request 
array with 10,000 and status array with 5,000. They should be enough in 
most cases. They use less than 150KB system memory. In rare case, one can 
use \code{mpi.realloc.comm}, \code{mpi.realloc.request} and
\code{mpi.realloc.status} to increase them to bigger arrayes.} 
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
