\name{mpi.remote.exec}
\alias{mpi.remote.exec}
\title{Remote Executions on R slaves}
\description{
 Remotely execute a command on R slaves spawned by using slavedaemon.R script 
 and return all executed results back to master. 
} 
\usage{
mpi.remote.exec(cmd, ..., simplify = TRUE, comm = 1, ret = TRUE)  
}
\arguments{
  \item{cmd}{the command to be executed on R slaves}
  \item{...}{used as arguments to cmd (function command) for passing their 
  (master) values to R slaves, i.e., if `myfun(x)' will be executed on R slaves 
with `x' as master variable, use mpi.remote.exec(cmd=myfun, x).}
  \item{simplify}{logical; should the result be simplified to a data.frame if possible?}
  \item{comm}{a communicator number.}
  \item{ret}{return executed results from R slaves if TRUE.}
}
\value{
  return executed results from R slaves if the argument \code{ret} is 
  set to be TRUE. The value could be a data.frame if values 
  (integer or double) from each slave have the same dimension. 
  Otherwise a list is returned.}
\details{
  Once R slaves are spawned by \code{\link{mpi.spawn.Rslaves}} with the 
  slavedaemon.R script, they are waiting for instructions from master. One can 
  use \code{\link{mpi.bcast.cmd}} to send a command to R slaves. However it 
  will not return executed results. Hence \code{mpi.remote.exec} can be 
  considered an extension to \code{\link{mpi.bcast.cmd}}.    
}

\section{Warning}{
  \code{mpi.remote.exec} may have difficult guessing invisible results 
  on R slaves. Use \code{ret = FALSE} instead. 
} 
\seealso{
  \code{\link{mpi.spawn.Rslaves}},
  \code{\link{mpi.bcast.cmd}}
}
\author{
Hao Yu
}
\examples{
\donttest{
mpi.remote.exec(mpi.comm.rank())
 x=5
mpi.remote.exec(rnorm,x)
}
}

\keyword{utilities}
