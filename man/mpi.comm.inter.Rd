\name{mpi.comm.get.parent}
\title{MPI\_Comm\_get\_parent, MPI\_Comm\_remote\_size, MPI\_Comm\_test\_inter 
APIs}
\usage{
  mpi.comm.get.parent(comm = 2)
  mpi.comm.remote.size(comm = 2)
  mpi.comm.test.inter(comm = 2)
}
\alias{mpi.comm.get.parent}
\alias{mpi.comm.remote.size}
\alias{mpi.comm.test.inter}
\arguments{
  \item{comm}{an intercommunicator number.}
}
\description{
  \code{mpi.comm.get.parent} is mainly used by slaves to find the 
intercommunicator or the parent who spawns them. The intercommunicator is saved 
in the specified comm number.

\code{mpi.comm.remote.size} is mainly used by master to find the total number of 
slaves spawned.

\code{mpi.comm.test.inter} tests if a comm is an intercomm or not.
}

\value{
\code{mpi.comm.get.parent} and \code{mpi.comm.test.inter} return  1 if success 
and 0 otherwise. 

\code{mpi.comm.remote.size} returns the total number of members in the remote 
group in an intercomm. 
}
\author{
Hao Yu
}
\references{
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
 \code{\link{mpi.intercomm.merge}}
}
%\examples{
%}
\keyword{interface}
 
