\name{mpi.bcast}
\title{MPI\_Bcast API}
\alias{mpi.bcast}
\description{
  \code{mpi.bcast} is a collective call among all members in a comm. It 
  broadcasts a message from the specified rank to all members. 
}
\usage{
mpi.bcast(x, type, rank = 0, comm = 1)
}
\arguments{
  \item{x}{data to be sent or received. Must be the same 
  type among all members.}
  \item{type}{1 for integer, 2 for double, and 3 for 
  character. Others are not supported.}
  \item{rank}{the sender.}
  \item{comm}{a communicator number.}
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
