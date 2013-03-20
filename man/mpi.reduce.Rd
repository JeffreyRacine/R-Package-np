\name{mpi.reduce}
\title{MPI\_Reduce and MPI\_Allreduce APIs}
\alias{mpi.allreduce}
\alias{mpi.reduce}
\description{
  \code{mpi.reduce} and \code{mpi.allreduce} are global reduction operations. 
  \code{mpi.reduce} combines each member's result, using the operation 
  \code{op}, and returns the combined value(s) to the member specified by 
  the argument \code{dest}. \code{mpi.allreduce} is the same as 
  \code{mpi.reduce} except that all members receive the combined value(s).
}
\usage{
mpi.reduce(x, type=2, op=c("sum","prod","max","min","maxloc","minloc"), 
	dest = 0, comm = 1) 

mpi.allreduce(x, type=2, op=c("sum","prod","max","min","maxloc","minloc"), 
	comm = 1) 
}
\arguments{
  \item{x}{data to be reduced. Must be the same dim and the same type for all 
  members.}
  \item{type}{1 for integer and 2 for double. Others are not supported.}
  \item{op}{one of "sum", "prod", "max", "min", "maxloc", or "minloc".} 
  \item{dest}{rank of destination}
  \item{comm}{a communicator number}
}
\value{
  \code{mpi.reduce} returns the combined value(s) to the member specified  
  by \code{dest}. \code{mpi.allreduce} returns the combined values(s) to 
  every member in a comm. The combined value(s) may be the summation, 
  production, maximum, or minimum specified by the argument \code{op}. If 
  the \code{op} is either "maxloc" or "minloc", then the maximum (minimum) 
  value(s) along the maximum (minimum) rank(s) will be returned.
}
\details{
  It is important that all members in a comm call either all \code{mpi.reduce} 
  or all \code{mpi.allreduce} even though the master may not be in 
  computation. They must provide exactly the same type and dim 
  vectors to be reduced. If the operation "maxloc" or "minloc" is used, 
  the combined vector is twice as long as the original one since the 
  maximum or minimum ranks are included. 
} 
\author{
  Hao Yu
}
\references{
  \url{http://www.lam-mpi.org/},  
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\seealso{
  \code{\link{mpi.gather}}.
}
%\examples{
%}
\keyword{interface}
