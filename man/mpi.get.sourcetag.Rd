\name{mpi.get.sourcetag}
\title{Utility for finding the source and tag of a received message}
\usage{
  mpi.get.sourcetag(status = 0)
}
\alias{mpi.get.sourcetag}
\arguments{
  \item{status}{a status number}
}
\description{
  \code{mpi.get.sourcetag} finds the source and tag of a received message.
}
\details{
  When \code{\link{mpi.any.source}} and/or \code{\link{mpi.any.tag}} are 
  used by \code{\link{mpi.recv}} or \code{\link{mpi.probe}}, one can use 
\code{mpi.get.sourcetag} 
  to find who sends the message or with what a tag number. 
  \code{\link{mpi.get.sourcetag}} must be called immediately after calling 
  \code{mpi.recv} or \code{mpi.probe} otherwise the obtained information may not 
be right.
}
\value{
  2 dim int vector. The first integer is the source and the second is the 
tag.
}
\author{
Hao Yu
}
\references{
  \url{http://www.lam-mpi.org/}, 
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
 \code{\link{mpi.send}}, \code{\link{mpi.recv}}, \code{\link{mpi.probe}},
 \code{\link{mpi.get.count}}
}
%\examples{
%}
\keyword{utilities}
 
