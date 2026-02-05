\name{mpi.get.count}
\title{MPI\_Get\_count API}
\usage{
mpi.get.count(type, status = 0)
}
\alias{mpi.get.count}
\arguments{
  \item{type}{1 for integer, 2 for double, 3 for char.}
  \item{status}{a status number}
}
\description{
  \code{mpi.get.count} finds the length of a received message.
}
\details{
  When \code{\link{mpi.recv}} is used to receive a message, the receiver 
  buffer can be set to be bigger than the incoming message. To find the 
  exact length of the received message, \code{mpi.get.count} is used to 
  find its exact length. \code{mpi.get.count} must be called 
  immediately after calling \code{mpi.recv} otherwise the status may be 
  changed.
}
\value{
  the length of a received message.
}
\author{
Hao Yu
}
\references{
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
 \code{\link{mpi.send}}, \code{\link{mpi.recv}},
 \code{\link{mpi.get.sourcetag}}, \code{\link{mpi.probe}}.
}
%\examples{
%}
\keyword{interface}
 
