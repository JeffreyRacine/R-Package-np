\name{mpi.send}
\title{MPI\_Send, MPI\_Isend, MPI\_Recv, and MPI\_Irecv APIs}
\alias{mpi.send}
\alias{mpi.isend}
\alias{mpi.recv}
\alias{mpi.irecv}
\description{
  The pair \code{mpi.send} and \code{mpi.recv} are two most used blocking 
  calls for point-to-point communications. An int, double or char vector 
  can be transmitted from any source to any destination.  

  The pair \code{mpi.isend} and \code{mpi.irecv} are the same except that 
  they are nonblocking calls.

  Blocking and nonblocking calls are interchangeable, e.g., nonblocking 
sends can be matched with blocking receives, and vice-versa.
}
\usage{
mpi.send(x, type, dest, tag,  comm = 1)
mpi.isend(x, type, dest, tag,  comm = 1, request=0)
mpi.recv(x, type, source, tag,  comm = 1, status = 0)
mpi.irecv(x, type, source, tag,  comm = 1, request = 0)
}
\arguments{
  \item{x}{data to be sent or received. Must be the same 
  type for source and destination. The receive buffer must be as large as 
  the send buffer.}
  \item{type}{1 for integer, 2 for double, and 3 for 
  character. Others are not supported.}
  \item{dest}{the destination rank. Use \code{mpi.proc.null} for a 
  fake destination.}
  \item{source}{the source rank. Use \code{mpi.any.source} for any source. 
   Use \code{mpi.proc.null} for a fake source. }
  \item{tag}{non-negative integer. Use \code{mpi.any.tag} for any tag flag.}
  \item{comm}{a communicator number.}
  \item{request}{a request number.}
  \item{status}{a status number.}
}
\value{
  \code{mpi.send} and \code{mpi.isend} return no value. \code{mpi.recv} 
returns the int, double or char vector sent from \code{source}. However, 
\code{mpi.irecv} returns no value. See details for explanation.
} 
\details{
  The pair \code{mpi.send} (or \code{mpi.isend}) and \code{mpi.recv} 
(or \code{mpi.irecv}) must be used together, i.e., if there is a sender, 
then there must be a receiver. Any mismatch will result a deadlock 
situation, i.e., programs stop responding. The receive buffer must be 
large enough to contain an incoming message otherwise programs will be 
crashed. One can use \code{\link{mpi.probe}} (or \link{mpi.iprobe}) and 
\code{\link{mpi.get.count}} to find the length of an incoming message 
before calling \code{mpi.recv}. If \code{\link{mpi.any.source}} or  
\code{\link{mpi.any.tag}} is used in \code{mpi.recv}, one can use  
\code{\link{mpi.get.sourcetag}} to find out the source or tag of the  
received message. To send/receive an R object rather than an int, double  
or char vector, please use the pair \code{\link{mpi.send.Robj}} and  
\code{\link{mpi.recv.Robj}}.

Since \code{mpi.irecv} is a nonblocking call, \code{x} with enough buffer 
must be created before using it. Then use nonblocking completion calls 
such as \code{\link{mpi.wait}} or \code{\link{mpi.test}} to test if 
\code{x} contains data from sender.

If multiple nonblocking sends or receives are used, please use request 
number consecutively from 0. For example, to receive two messages from two 
slaves, try 
mpi.irecv(x,1,source=1,tag=0,comm=1,request=0)
mpi.irecv(y,1,source=2,tag=0,comm=1,request=1)
Then \code{mpi.waitany}, \code{mpi.waitsome} or \code{mpi.waitall} can be 
used to complete the operations.
 }
\author{
Hao Yu
}
\references{
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
  \code{\link{mpi.send.Robj}},
  \code{\link{mpi.recv.Robj}},
  \code{\link{mpi.probe}},  
  \code{\link{mpi.wait}},  
  \code{\link{mpi.get.count}}, 
  \code{\link{mpi.get.sourcetag}}. 
} 
\examples{ 
\donttest{
#on a slave
mpi.send(1:10,1,0,0)

#on master
x <- integer(10)
mpi.irecv(x,1,1,0)
x	
mpi.wait()
x
}
} 
\keyword{interface}

