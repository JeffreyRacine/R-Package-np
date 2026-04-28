\name{mpi.exit}
\alias{mpi.exit}
\alias{mpi.quit}
\title{Exit MPI Environment }
\description{
\code{mpi.exit} terminates MPI execution environment and detaches the 
package \code{npRmpi}. After that, you can still work in R.

\code{mpi.quit} terminates MPI execution environment and quits R.
}
\usage{
mpi.exit()
mpi.quit(save = "no")
}
\arguments{
\subsection{Exit Controls}{
R-session save behavior used by \code{mpi.quit()}.
}

\item{save}{the same argument as \code{quit} but default to "no".}
}
\value{
 \code{mpi.exit} always returns 1
}
\details{
  Normally, MPI finalization is used to clean all MPI states.
  However, it will not detach the loaded \code{npRmpi} package. To be safer when leaving MPI,
  \code{mpi.exit} not only calls \code{mpi.finalize} but also detaches the 
  \code{npRmpi} package. This will make reloading \code{npRmpi} impossible
  in the same R session.

 If leaving MPI and R altogether, one simply uses \code{mpi.quit}. 
} 
\author{
Hao Yu
}
\keyword{utilities}
