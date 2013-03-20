\name{mpi.finalize}
\alias{mpi.finalize}
\title{MPI\_Finalize API}
\description{
  Terminates MPI execution environment. 
}
\usage{
  mpi.finalize()
}
\arguments{
  None
}
\value{
  Always return 1
}
\details{
  This routines must be called by each slave (master) before it exits. This 
  call cleans all MPI state. Once \code{mpi.finalize} has been called, no MPI 
  routine may be called. To be more safe leaving MPI, please use 
  \code{\link{mpi.exit}} which not only calls \code{mpi.finalize} but also 
  detaches the library Rmpi. This will make reload the library Rmpi impossible. 
} 
\author{
Hao Yu
}
\seealso{
  \code{\link{mpi.exit}}
}
\references{
  \url{http://www.lam-mpi.org/},
  \url{http://www-unix.mcs.anl.gov/mpi/www/www3/}
}
%\examples{ 
%}
\keyword{interface}
