\name{mpi.exit}
\alias{mpi.exit}
\alias{mpi.quit}
\title{Exit MPI Environment }
\description{
\code{mpi.exit} terminates MPI execution environment and detaches the 
library Rmpi. After that, you can still work on R.

\code{mpi.quit} terminates MPI execution environment and quits R.
}
\usage{
mpi.exit()
mpi.quit(save = "no")
}
\arguments{
\item{save}{the same argument as \code{quit} but default to "no".}
}
\value{
 \code{mpi.exit} always returns 1
}
\details{
  Normally, \code{\link{mpi.finalize}} is used to clean all MPI states. 
  However, it will not detach the library Rmpi. To be more safe leaving MPI,  
  \code{mpi.exit} not only calls \code{mpi.finalize} but also detaches the 
  library Rmpi. This will make reload the library Rmpi impossible.

 If leaving MPI and R altogether, one simply uses \code{mpi.quit}. 
} 
\seealso{
  \code{\link{mpi.finalize}}
}
\author{
Hao Yu
}
\keyword{utilities}
