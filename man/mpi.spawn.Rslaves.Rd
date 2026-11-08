\name{mpi.close.Rslaves}
\alias{mpi.close.Rslaves}
\alias{tailslave.log}
\title{Close and Inspect R Slaves}
\description{
\code{mpi.close.Rslaves} shuts down (or soft-closes) R slave daemons managed by
\code{npRmpi}. \code{tailslave.log} shows tail output from slave log files.
}
\usage{
mpi.close.Rslaves(dellog = TRUE, comm = 1, force = FALSE)
tailslave.log(nlines = 3, comm = 1)
}
\arguments{
\subsection{Slave Shutdown Controls}{
Communicator, log-deletion, and hard-shutdown controls for slave daemons.
}

  \item{dellog}{a logical specifying if R slave log files are deleted.}
  \item{comm}{a communicator number.}
  \item{force}{a logical. If \code{TRUE}, force a hard shutdown of slave daemons.
  When \code{options(npRmpi.reuse.slaves=TRUE)} and \code{force=FALSE},
  \code{mpi.close.Rslaves()} performs a soft-close (keeps daemons alive for reuse).}

\subsection{Slave Log Inspection}{
Tail length used when inspecting slave log files.
}

  \item{nlines}{number of lines shown from the tail of each slave log file.}
}
\value{
\code{mpi.close.Rslaves} returns a status code (in soft-close mode this may be
an invisible no-op code). \code{tailslave.log} returns the tail output from
slave log files.
}
\details{
The package-owned startup routes compare the loaded R and native source
identities on participating ranks. A mismatch or a worker from an older build
without an identity record refuses admission. Close the pool and restart all
ranks with matching package sources; changing files underneath a running R
session does not update its loaded code. Different library paths and supported
platform builds of the same sources are allowed. Identity checking runs at
startup/reuse, not for every numerical call. It is an accidental source-mismatch
check, not authentication of modified binaries or subsequent namespace changes.
Low-level MPI worlds that bypass these startup helpers and the shipped profile
remain caller-managed and are not certified by this check.


In normal user workflows, call \code{npRmpi.quit()} rather than using
\code{mpi.close.Rslaves()} directly.

\code{tailslave.log()} is useful for debugging worker startup or teardown issues.
}
\seealso{
  \code{\link{npRmpi.init}},
  \code{\link{npRmpi.quit}}.
}
\author{
Hao Yu
}
\examples{
\dontrun{
# Inspect slave logs from the current communicator.
tailslave.log()

# Close active slave daemons.
mpi.close.Rslaves()
}
}
\keyword{utilities}
