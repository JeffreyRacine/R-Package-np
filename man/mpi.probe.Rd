\name{mpi.probe}
\title{MPI\_Probe and MPI\_Iprobe APIs}
\usage{
mpi.probe(source, tag, comm = 1, status = 0)
mpi.iprobe(source, tag, comm = 1, status = 0)
}
\alias{mpi.probe}
\alias{mpi.iprobe}
\arguments{
  \item{source}{the source of incoming message or mpi.any.source() for any 
source.}
  \item{tag}{a tag number or mpi.any.tag() for any tag.}
  \item{comm}{a communicator number}
  \item{status}{a status number}
}
\description{
  \code{mpi.probe} uses the source and tag of incoming message to set a 
status. \code{mpi.iprobe} does the same except it is a nonblocking call, 
i.e., returns immediately.}
\details{
  When \code{\link{mpi.send}} or other nonblocking sends are used to send 
a message, the receiver may not know the exact length before receiving 
it. \code{mpi.probe} is used to probe the incoming message and put some 
information into a status. Then the exact length can be found by using 
\code{\link{mpi.get.count}} to such a status. If the wild card 
\code{mpi.any.source} or \code{mpi.any.tag} are used, then one 
can use \code{\link{mpi.get.sourcetag}} to find the exact source or tag of 
a sender. }
\value{
 \code{mpi.probe} returns 1 only after a matching message has been found.

 \code{mpi.iproble} returns TRUE if there is a message that can be 
received; FALSE otherwise.
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
 \code{\link{mpi.get.count}}
}
%\examples{
%}
\keyword{interface}
 
