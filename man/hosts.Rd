\name{lamhosts}
\alias{mpi.is.master}
\alias{lamhosts}
\alias{mpi.hostinfo}
\alias{slave.hostinfo}
\alias{mpichhosts}
\title{Hosts Information}
\description{
\code{lamhosts} finds the host name associated with its node number. Can be used 
by \code{\link{mpi.spawn.Rslaves}} to spawn R slaves on selected hosts. This is 
a LAM-MPI specific function.

\code{mpi.is.master} checks if it is running on master or slaves.

\code{mpi.hostinfo} finds an individual host information including rank and 
size in a comm.

\code{slave.hostinfo} is executed only by master and find all master and slaves 
host information in a comm.

\code{mpichhosts} finds host names from master Windows registery database.
}
\usage{
lamhosts()
mpi.is.master()
mpi.hostinfo(comm = 1)
slave.hostinfo(comm = 1, short=TRUE) 
mpichhosts()
}
\arguments{
 \item{comm}{a communicator number}
 \item{short}{if true, a short form is printed}
}
\value{
  \code{lamhosts} returns CPUs nodes numbers with their host names.

  \code{mpi.is.master} returns TRUE if it is on master and FALSE otherwise.

  \code{mpi.hostinfo} sends to stdio a host name, rank, size and comm.

  \code{slave.hostname} sends to stdio a list of host, rank, size, and comm 
information for all master and slaves. With short=TRUE and 8 slaves or more, 
the first 3 and last 2 slaves are shown.

  \code{mpichhosts} returns hostnames including multiple same hostname if a host 
has multiple CPUs. The first host is labeled as "master", the rest is labeled as 
"slave1", "slave2", and so on. \code{mpi.universe.size} will return the length of 
\code{mpichhosts}-1 on Windows platform. } 
\seealso{
  \code{\link{mpi.spawn.Rslaves}}
}
\author{
Hao Yu (minor modifications by Jeffrey S. Racine \email{racinej@mcmaster.ca})
}
\keyword{utilities}
