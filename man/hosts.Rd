\name{lamhosts}
\alias{mpi.is.master}
\alias{lamhosts}
\title{Hosts Information}
\description{
\code{lamhosts} finds the host name associated with its node number. It can be
used with \code{\link{npRmpi.init}} (\code{mode="spawn"}), which internally
uses \code{mpi.spawn.Rslaves}, to start R slaves on selected hosts. This is a
MPI implementation specific function.

\code{mpi.is.master} checks if it is running on master or slaves.
}
\usage{
lamhosts()
mpi.is.master()
}
\value{
  \code{lamhosts} returns CPUs nodes numbers with their host names.

  \code{mpi.is.master} returns TRUE if it is on master and FALSE otherwise.
} 
\seealso{
  \code{\link{npRmpi.init}},
  \code{\link{mpi.hostinfo}},
  \code{\link{slave.hostinfo}}
}
\author{
Hao Yu (minor modifications by Jeffrey S. Racine \email{racinej@mcmaster.ca})
}
\keyword{utilities}
