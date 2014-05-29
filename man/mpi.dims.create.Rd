\name{mpi.dims.create}
\title{MPI\_Dims\_create}
\alias{mpi.dims.create}
\description{
 \code{mpi.dims.create} Create a Cartesian dimension used by \code{mpi.cart.create}.
}
\usage{
 mpi.dims.create(nnodes, ndims, dims=integer(ndims))
}
\arguments{
 \item{nnodes}{Number of nodes in a cluster}
 \item{ndims}{Number of dimension in a Cartesian topology}
 \item{dims}{Initial dimension numbers}
}
\value{
 \code{mpi.dims.create} returns the dimension vector used by 
that in \code{mpi.cart.create}.
}
\details{
The entries in the return value are set to describe a Cartesian grid with 
\code{ndims} dimensions and a total of \code{nnodes} nodes. The dimensions are set 
to be as close to each other as possible, using an appropriate divisibility 
algorithm. The return value can be constrained by specifying positive number(s) in 
\code{dims}. Only those 0 values in \code{dims} are modified by 
\code{mpi.dims.create}.} 
\author{
Hao Yu
}
\seealso{
\code{\link{mpi.cart.create}}
}
\references{
  \url{http://www.lam-mpi.org/},
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
\examples{
\donttest{
#What is the dim numbers of 2 dim Cartersian topology under a grid of 36 nodes
mpi.dims.create(36,2)	#return c(6,6)
#Constrained dim numbers
mpi.dims.create(12,2,c(0,4)) #return c(9,4)
}
} 
\keyword{interface}
