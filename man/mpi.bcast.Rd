\name{mpi.bcast}
\title{MPI\_Bcast API}
\alias{mpi.bcast}
\description{
  \code{mpi.bcast} is a collective call among all members in a comm. It 
  broadcasts a message from the specified rank to all members. 
}
\usage{
mpi.bcast(x, type, rank = 0, comm = 1, buffunit=100)
}
\arguments{
  \item{x}{data to be sent or received. Must be the same 
  type among all members.}
  \item{type}{1 for integer, 2 for double, and 3 for 
  character. Others are not supported.}
  \item{rank}{the sender.}
  \item{comm}{a communicator number.}
  \item{buffunit}{a buffer unit number.}
}
\value{
  \code{mpi.bcast} returns the message broadcasted by the sender 
  (specified by the rank).
} 
\details{
  \code{mpi.bcast} is a blocking call among all members in a comm, i.e, 
  all members have to wait until everyone calls it. All members have to 
  prepare the same type of messages (buffers). Hence it is relatively 
  difficult to use in R environment since the receivers may not know what 
  types of data to receive, not mention the length of data. Users should 
  use various extensions of \code{mpi.bcast} in R. They are 
  \code{\link{mpi.bcast.Robj}}, \code{\link{mpi.bcast.cmd}}, and 
  \code{\link{mpi.bcast.Robj2slave}}.
  
  When type=5, MPI continuous datatype (double) is defined with unit given by 
  \code{buffunit}. It is used to transfer huge data where a double vector or matrix
  is divided into many chunks with unit \code{buffunit}. Total 
  ceiling(length(obj)/buffunit) units are transferred. Due to MPI specification, both
  \code{buffunit} and total units transferred cannot be over 2^31-1. Notice that the last 
  chunk may not have full length of data due to rounding. Special care is needed.
  }
\references{
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
  \code{\link{mpi.bcast.Robj}},
  \code{\link{mpi.bcast.cmd}},
  \code{\link{mpi.bcast.Robj2slave}}.
}
%\examples{
%}
\keyword{interface}
