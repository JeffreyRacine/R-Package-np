\name{mpi.gather}
\title{MPI\_Gather, MPI\_Gatherv, MPI\_Allgather, and MPI\_Allgatherv APIs}
\alias{mpi.gather}
\alias{mpi.gatherv}
\alias{mpi.allgather}
\alias{mpi.allgatherv}
\description{
  \code{mpi.gather} and \code{mpi.gatherv} (vector variant) gather each 
  member's message to the member specified by the argument \code{root}.
  The root member receives the messages and stores them in rank 
  order. \code{mpi.allgather} and \code{mpi.allgatherv} are the same as
  \code{mpi.gather} and \code{mpi.gatherv} except that all members receive
  the result instead of just the root.
}
\usage{
mpi.gather(x, type, rdata, root = 0, comm = 1) 
mpi.gatherv(x, type, rdata, rcounts, root = 0, comm = 1) 

mpi.allgather(x, type, rdata, comm = 1) 
mpi.allgatherv(x, type, rdata, rcounts, comm = 1) 
}
\arguments{
  \item{x}{data to be gathered. Must be the same type.}
  \item{type}{1 for integer, 2 for double, and 3 for character. Others are not 
supported.}
  \item{rdata}{the receive buffer. Must be the same type as the sender and big 
  enough to include all message gathered.}
  \item{rcounts}{int vector specifying the length of each message.}
  \item{root}{rank of the receiver}
  \item{comm}{a communicator number}
}
\value{
  For \code{mpi.gather} or \code{mpi.gatherv}, it returns the gathered 
  message for the root member. For other members, it returns what is in rdata, 
  i.e., rdata (or rcounts) is ignored. For \code{mpi.allgather} or 
  \code{mpi.allgatherv}, it returns the gathered message for all members. 
} 
\details{
  For \code{mpi.gather} and \code{mpi.allgather}, the message to be gathered 
must be the same dim and the same type. The receive buffer can be prepared as 
either integer(size * dim) or double(size * dim), where size is the total 
number of members in a comm. For \code{mpi.gatherv} and \code{mpi.allgatherv}, 
the message to be gathered can have different dims but must be the same type. 
The argument \code{rcounts} records these different dims into an integer vector 
in rank order. Then the receive buffer can be prepared as either 
integer(sum(rcounts)) or double(sum(rcounts)).  
}
\author{
Hao Yu
}
\references{
  \url{http://www.lam-mpi.org/},
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
  \code{\link{mpi.scatter}}, \code{\link{mpi.scatterv}}.
}
\examples{
\donttest{
#Need 3 slaves to run properly
#Or use mpi.spawn.Rslaves(nslaves=3)
 mpi.bcast.cmd(id <-mpi.comm.rank(.comm), comm=1)
mpi.bcast.cmd(mpi.gather(letters[id],type=3,rdata=string(1)))
mpi.gather(letters[10],type=3,rdata=string(4))

 mpi.bcast.cmd(x<-rnorm(id))
 mpi.bcast.cmd(mpi.gatherv(x,type=2,rdata=double(1),rcounts=1))
 mpi.gatherv(double(1),type=2,rdata=double(sum(1:3)+1),rcounts=c(1,1:3))

mpi.bcast.cmd(out1<-mpi.allgatherv(x,type=2,rdata=double(sum(1:3)+1),
		rcounts=c(1,1:3)))
mpi.allgatherv(double(1),type=2,rdata=double(sum(1:3)+1),rcounts=c(1,1:3))
}
}
\keyword{interface}
