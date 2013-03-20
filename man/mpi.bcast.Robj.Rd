\name{mpi.bcast.Robj}
\title{Extensions of MPI\_Bcast API}
\alias{mpi.bcast.Robj}
\alias{mpi.bcast.Robj2slave}
\alias{mpi.bcast.Rfun2slave}
\description{
  \code{mpi.bcast.Robj} and \code{mpi.bcast.Robj2slave} are used to move 
a general R object around among master and all slaves. 
}
\usage{
mpi.bcast.Robj(obj = NULL, rank = 0, comm = 1)
mpi.bcast.Robj2slave(obj, comm = 1, all = FALSE)
mpi.bcast.Rfun2slave(comm = 1)
}
\arguments{
  \item{obj}{an R object to be transmitted from the sender}
  \item{rank}{the sender.}
  \item{comm}{a communicator number.}
  \item{all}{a logical. If TRUE, all R objects on master are transmitted to slaves.}
}
\value{
  \code{mpi.bcast.Robj} returns no value for the sender and the 
transmitted one for others. \code{mpi.bcast.Robj2slave} returns no value for the master 
and the transmitted R object along its name on slaves. \code{mpi.bcast.Rfun2slave} 
transmits all master's functions to slaves and returns no value.
} 
\details{
  \code{mpi.bcast.Robj} is an extension of \code{\link{mpi.bcast}} for 
  moving a general R object around from a sender to everyone. 
  \code{mpi.bcast.Robj2slave} does an R object transmission from 
  master to all slaves unless all=TRUE in which case, all master's objects with 
  the global enviroment are transmitted to all slavers.  
}
\author{
Hao Yu
}
\seealso{
  \code{\link{mpi.send.Robj}},
  \code{\link{mpi.recv.Robj}},
}
%\examples{
%}
\keyword{utilities}
