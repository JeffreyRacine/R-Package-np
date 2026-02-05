\name{mpi.comm.spawn}
\alias{mpi.comm.spawn}
\title{MPI\_Comm\_spawn API }
\description{
\code{mpi.comm.spawn} tries to start \code{nslaves} identical copies of 
\code{slaves}, establishing communication with them and returning an 
intercommunicator. The spawned slaves are referred to as children, and the 
process that spawned them is called the parent (master). The children have 
their own MPI\_COMM\_WORLD represented by comm 0. To make communication 
possible among master and slaves, all slaves should use 
\code{\link{mpi.comm.get.parent}} to find their parent and use 
\code{\link{mpi.intercomm.merge}} to merger an intercomm to a comm.
 }
\usage{
 mpi.comm.spawn(slave, slavearg = character(0),
                nslaves = mpi.universe.size(), info = 0,
                root = 0, intercomm = 2, quiet = FALSE)
}
\arguments{
  \item{slave}{a file name to an executable program.}
  \item{slavearg}{an argument list (a char vector) to slave.}
  \item{nslaves}{number of slaves to be spawned.}
  \item{info}{an info number.}
  \item{root}{the root member who spawns slaves.}
  \item{intercomm}{an intercomm number.}
  \item{quiet}{a logical. If TRUE, do not print anything unless an error occurs.}
}
\value{
Unless \code{quiet = TRUE}, a message is printed to indicate how many slaves are successfully 
spawned and how many failed.
}
\seealso{
  \code{\link{mpi.comm.get.parent}},
  \code{\link{mpi.intercomm.merge}}.
}
\references{
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}

\author{
Hao Yu
}
\keyword{interface}
