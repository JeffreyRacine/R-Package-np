\name{mpi.send.Robj}
\title{Extensions of MPI\_Send and MPI\_Recv APIs}
\alias{mpi.send.Robj}
\alias{mpi.isend.Robj}
\alias{mpi.recv.Robj}
\description{
  \code{mpi.send.Robj} and \code{mpi.recv.Robj} are two 
  extensions of \code{mpi.send} and \code{mpi.recv}. They are used to 
  transmit a general R object from any source to any destination.  

  \code{mpi.isend.Robj} is a nonblocking version of \code{mpi.send.Robj}.
}
\usage{
mpi.send.Robj(obj, dest, tag, comm = 1)
mpi.isend.Robj(obj, dest, tag, comm = 1, request=0)
mpi.recv.Robj(source, tag, comm = 1, status = 0)
}
\arguments{
  \item{obj}{an R object. Can be any R object. }
  \item{dest}{the destination rank.}
  \item{source}{the source rank or mpi.any.source() for any source.} 
  \item{tag}{non-negative integer or mpi.any.tag() for any tag.}
  \item{comm}{a communicator number.}
  \item{request}{a request number.}
  \item{status}{a status number.}
}
\value{
  \code{mpi.send.Robj} or \code{mpi.isend.Robj} return no value. 
  \code{mpi.recv.Robj} returns the the transmitted R object. 
} 
\details{
   \code{mpi.send.Robj} and \code{mpi.isend.Robj} use 
\code{serialize} to encode  an R object into a binary 
char vector. It sends the message to the destination. The receiver 
decode the message back into an R object by using  
\code{unserialize}. 

If \code{mpi.isend.Robj} is used, \code{mpi.wait} or \code{mpi.test} must 
be used to check the object has been sent.
 }
\author{
Hao Yu
}
\references{
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
  \code{\link{mpi.send}},
  \code{\link{mpi.recv}},
  \code{\link{mpi.wait}},
  \code{\link{serialize}},
  \code{\link{unserialize}}, 
}
%\examples{
%}
\keyword{utilities}

