\name{mpi.gather.Robj}
\title{Extentions of MPI\_Gather and MPI\_Allgather APIs}
\alias{mpi.gather.Robj}
\alias{mpi.allgather.Robj}
\description{
  \code{mpi.gather.Robj} gathers each member's object to the member 
specified by the argument \code{root}.
  The root member receives the objects as a list.
  \code{mpi.allgather.Robj} is the same as \code{mpi.gather.Robj}
 except that all members receive the result instead of just the root.
}
\usage{
mpi.gather.Robj(obj=NULL, root = 0, comm = 1, ...)

mpi.allgather.Robj(obj=NULL, comm = 1)
}
\arguments{
  \item{obj}{data to be gathered. Could be different type.}
  \item{root}{rank of the gather}
  \item{comm}{a communicator number}
  \item{...}{optional arugments to \code{sapply}.}
}
\details{
Since sapply is used to gather all results, its default option 
"simplify=TRUE" is to simplify outputs. In some situations, this option 
is not desirable. Using "simplify=FALSE" as in the place of ... will tell
sapply not to simplify and a list of outputs will be returned.
}
\value{
  For \code{mpi.gather.Robj}, it returns a list, the gathered message 
for the root member. For
  \code{mpi.allgatherv.Robj}, it returns a list, the gathered message 
for all members.
}
\author{
Hao Yu and Wei Xia
}
\references{
  \url{http://www.lam-mpi.org/},
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
  \code{\link{mpi.gather}}, \code{\link{mpi.allgatherv}}.
}
\examples{
\dontrun{
#Assume that there are some slaves running
mpi.bcast.cmd(id<-mpi.comm.rank())
mpi.bcast.cmd(x<-rnorm(id))
mpi.bcast.cmd(mpi.gather.Robj(x))
x<-"test mpi.gather.Robj"
mpi.gather.Robj(x)

mpi.bcast.cmd(obj<-rnorm(id+10))
mpi.bcast.cmd(nn<-mpi.allgather.Robj(obj))
obj<-rnorm(5)
mpi.allgather.Robj(obj)
mpi.remote.exec(nn)
}
}
\keyword{interface}

