\name{mpi.scatter.Robj}
\title{Extensions of MPI\_ SCATTER and MPI\_SCATTERV }
\alias{mpi.scatter.Robj}
\alias{mpi.scatter.Robj2slave}
\description{
  \code{mpi.scatter.Robj} and \code{mpi.scatter.Robj2slave} are used to scatter a list 
   to all members. They are more efficient than using any parallel apply functions.
}
\usage{
mpi.scatter.Robj(obj = NULL, root = 0, comm = 1)
mpi.scatter.Robj2slave(obj, comm = 1)
}
\arguments{
  \item{obj}{a list object to be scattered from the root or master}
  \item{root}{rank of the scatter.}
  \item{comm}{a communicator number.}
}
\value{
  \code{mpi.scatter.Robj} for non-root members,  returns the
     scattered R object. For the root  member, it returns the
     portion belonging to itself. \code{mpi.scatter.Robj2slave} returns no value
     for the master and all slaves get their corresponding components in the list,
     i.e., the first slave gets the first component in the list. 
}
\details{
  \code{mpi.scatter.Robj} is an extension of \code{\link{mpi.scatter}} for
  scattering a list object from a sender (root) to everyone. \code{mpi.scatter.Robj2slave}
  scatters a list to all slaves.
}
\author{
Hao Yu  and   Wei Xia
}
\seealso{
  \code{\link{mpi.scatter}},
  \code{\link{mpi.gather.Robj}},
}
\examples{
\donttest{
#assume that there are three slaves running
mpi.bcast.cmd(x<-mpi.scatter.Robj())

xx <- list("master",rnorm(3),letters[2],1:10)
mpi.scatter.Robj(obj=xx)

mpi.remote.exec(x)

#scatter a matrix to slaves
dat=matrix(1:24,ncol=3)
splitmatrix = function(x, ncl) lapply(.splitIndices(nrow(x), ncl), function(i) x[i,])
dat2=splitmatrix(dat,3)
mpi.scatter.Robj2slave(dat2)
mpi.remote.exec(dat2)
}
}
\keyword{utilities}

