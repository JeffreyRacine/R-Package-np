\name{mpi.sendrecv}
\title{MPI\_Sendrecv and MPI\_Sendrecv\_replace APIs}
\alias{mpi.sendrecv}
\alias{mpi.sendrecv.replace}
\description{
\code{mpi.sendrecv} and \code{mpi.sendrecv.replace} execute blocking send 
and receive operations. Both of them combine the sending of one message to 
a destination and the receiving of another message from a source in one 
call. The source and destination are possibly the same. The send buffer and 
receive buffer are disjoint for \code{mpi.sendrecv}, while the buffers are 
not disjoint for \code{mpi.sendrecv.replace}.  
}
\usage{
mpi.sendrecv(senddata, sendtype, dest, sendtag, recvdata, recvtype, 
source, recvtag, comm = 1, status = 0)

mpi.sendrecv.replace(x, type, dest, sendtag, source, recvtag, 
comm = 1, status = 0)
}
\arguments{
  \item{x}{data to be sent or recieved. Must be the same 
  type for source and destination.}
  \item{senddata}{data to be sent. May have different datatypes and 
  lengths} 
  \item{recvdata}{data to be recieved. May have different datatypes and 
  lengths} 
  \item{type}{type of the data to be sent or recieved. 
   1 for integer, 2 for double, and 3 for character. Others are not supported.}
  \item{sendtype}{type of the data to be sent. 
   1 for integer, 2 for double, and 3 for character. Others are not supported.}
  \item{recvtype}{type of the data to be recieved. 
   1 for integer, 2 for double, and 3 for character. Others are not supported.}
  \item{dest}{the destination rank. Use \code{mpi.proc.null} for a 
  fake destination.}
  \item{source}{the source rank. Use \code{mpi.any.source} for any source. 
   Use \code{mpi.proc.null} for a fake source. }
  \item{sendtag}{non-negative integer. Use \code{mpi.any.tag} 
   for any tag flag.}
  \item{recvtag}{non-negative integer. Use \code{mpi.any.tag} 
   for any tag flag.}
  \item{comm}{a communicator number.}
  \item{status}{a status number.}
}
\value{
  Returns the int, double or char vector sent from the send buffers. 
} 
\details{
  The receive buffer must be large enough to contain an incoming message 
  otherwise programs will be crashed. There is compatibility between 
  send-receive and normal sends and receives. A message sent by a 
send-receive can be received  by a regular receive and a send-receive 
can receive a message sent by a regular send.
 }
\author{
Kris Chen
}
\references{
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
  \code{\link{mpi.send.Robj}},
  \code{\link{mpi.recv.Robj}},
  \code{\link{mpi.probe}}.  
  \code{\link{mpi.get.sourcetag}}. 
} 
\examples{
mpi.sendrecv(as.integer(11:20),1,0,33,integer(10),1,0,33,comm=0)
mpi.sendrecv.replace(seq(1,2,by=0.1),2,0,99,0,99,comm=0)
} 
\keyword{interface}

