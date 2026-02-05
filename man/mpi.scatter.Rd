\name{mpi.scatter}
\title{MPI\_Scatter and MPI\_Scatterv APIs}
\alias{mpi.scatter}
\alias{mpi.scatterv}
\description{
  \code{mpi.scatter} and \code{mpi.scatterv} are the inverse operations of 
  \code{\link{mpi.gather}} and \code{\link{mpi.gatherv}} respectively.
}
\usage{
mpi.scatter(x, type, rdata, root = 0,  comm = 1) 
mpi.scatterv(x, scounts, type, rdata, root = 0, comm = 1)
}
\arguments{
  \item{x}{data to be scattered.}
  \item{type}{1 for integer, 2 for double, and 3 for character. Others are not 
supported.}
  \item{rdata}{the receive buffer. Must be the same type as the sender} 
  \item{scounts}{int vector specifying the block length inside a message to be 
  scattered to other members.}
  \item{root}{rank of the receiver}
  \item{comm}{a communicator number}
}
\details{
  \code{mpi.scatter} scatters the message x to all members. Each member receives 
  a portion of x with dim as length(x)/size in rank order, where size is the 
  total number of members in a comm. So the receive buffer can be prepared as 
  either integer(length(x)/size) or double(length(x)/size). For 
\code{mpi.scatterv}, scounts counts the portions (different dims) of x sent to 
each member. Each member needs to prepare the receive buffer as either 
integer(scounts[i]) or double(scounts[i]). 
}
\value{
For non-root members, \code{mpi.scatter} or \code{scatterv} returns the
scattered message and ignores whatever is in x (or scounts). For the root 
  member, it returns the portion belonging to itself.
} 
\author{
Hao Yu
}
\references{
  \url{http://www.lam-mpi.org/},
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
  \code{\link{mpi.gather}}, \code{\link{mpi.gatherv}}.
}
\examples{
\donttest{
#Need 3 slaves to run properly
#Or run  mpi.spawn.Rslaves(nslaves=3)
  num="123456789abcd"
  scounts<-c(2,3,1,7)
  mpi.bcast.cmd(strnum<-mpi.scatter(integer(1),type=1,rdata=integer(1),root=0))
  strnum<-mpi.scatter(scounts,type=1,rdata=integer(1),root=0)
  mpi.bcast.cmd(ans <- mpi.scatterv(string(1),scounts=0,type=3,rdata=string(strnum),
					root=0))
  mpi.scatterv(as.character(num),scounts=scounts,type=3,rdata=string(strnum),root=0)
  mpi.remote.exec(ans)
}
}
\keyword{interface}
