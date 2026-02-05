\name{mpi.intercomm.merge}
\title{MPI\_Intercomm\_merge API}
\usage{
mpi.intercomm.merge(intercomm=2, high=0, comm=1)
}
\alias{mpi.intercomm.merge}
\arguments{
  \item{intercomm}{an intercommunicator number}
  \item{high}{Used to order the groups of the two intracommunicators within comm 
when creating the new communicator}
  \item{comm}{a (intra)communicator number}
}
\description{
  Creates an intracommunicator from an intercommunicator 
}
\details{
  When master spawns slaves, an intercommunicator is created. To make 
  communications (point-to-point or groupwise) among master and slaves, an 
  intracommunicator must be created. \code{mpi.intercomm.merge} is used for 
  that purpose. This is a collective call so all master and slaves call 
  together. R slaves spawned by \code{\link{mpi.spawn.Rslaves}} should use
\code{\link{mpi.comm.get.parent}} to get (set) an intercomm to a number followed 
by merging antercomm to an intracomm. One 
can use \code{\link{mpi.comm.test.inter}} to test if a 
  communicator is an intercommunicator or not.
}
\value{
  1 if success. Otherwise 0. 
}
\author{
Hao Yu
}
\references{
  \url{http://www.lam-mpi.org}, \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
 \code{\link{mpi.comm.test.inter}}
}
%\examples{
%}
\keyword{interface}
 
