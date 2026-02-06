\name{mpi.apply}
\alias{mpi.apply}
\alias{mpi.iapply}
\title{Scatter an array to slaves and then apply a FUN}
\description{
An array (length <= total number of slaves) is scattered to slaves so that the first 
slave calls \code{FUN} with arguments \code{x[[1]]} and \code{...}, the second one 
calls with arguments \code{x[[2]]} and \code{...}, and so on. \code{mpi.iapply} is a 
nonblocking version of \code{mpi.apply} so that it will not consume CPU on master node.} 
\usage{
mpi.apply(X, FUN, ..., comm=1)  
mpi.iapply(X, FUN, ..., comm=1, sleep=0.01)  
}
\arguments{
  \item{X}{an array}
  \item{FUN}{a function}
  \item{...}{optional arguments to \code{FUN}}
  \item{comm}{a communicator number}
  \item{sleep}{a sleep interval on master node (in sec)}
}
\value{
A list of the results is returned. Its length is the same as that of \code{x}. In 
case the call \code{FUN} with arguments \code{x[[i]]} and \code{...} fails on ith 
slave, corresponding error message will be returned in the returning list.  }

%\seealso{
%  \code{\link{mpi.spawn.Rslaves}},
%  \code{\link{mpi.bcast.cmd}}
%  \code{\link{mpi.bcast.Robj2slave}}
%}
\author{
Hao Yu
}
\examples{
\dontrun{
#Assume that there are at least 5 slaves running
#Otherwise run mpi.spawn.Rslaves(nslaves=5)
x=c(10,20)
mpi.apply(x,runif)
meanx=1:5
mpi.apply(meanx,rnorm,n=2,sd=4)
}
}
\keyword{utilities}
