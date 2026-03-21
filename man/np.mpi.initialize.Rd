\name{np.mpi.initialize}
\alias{np.mpi.initialize}
\title{Initialize Ranks for Manual-Broadcast npRmpi Workflows}
\description{
  \code{np.mpi.initialize} initializes the caller and worker ranks for the
  profile/manual-broadcast \pkg{npRmpi} workflow.
}
\details{
  \code{np.mpi.initialize()} is the helper used after ranks have already
  been started with the profile/manual-broadcast route. The usual pattern
  is:
  \code{mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)}.

  This helper is \emph{not} the ordinary entry point for session or attach
  workflows. For those routes, use \code{\link{npRmpi.init}} and
  \code{\link{npRmpi.quit}} instead.

  For MPI startup/performance guidance (including message-passing tradeoffs and the manual-broadcast template), see \code{\link{npRmpi.init}} details and \code{inst/Rprofile}.
  Documentation guide: see \code{\link{np.kernels}} for kernels,
  \code{\link{np.options}} for global options, and
  \code{\link{plot}} for plotting options.
}
\usage{
np.mpi.initialize()
}

\value{
  \code{np.mpi.initialize} returns no value for the sender and an expression
  of the transmitted command for others.
} 

\author{
 Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}},
  \code{\link{plot}}, \code{\link{npRmpi.init}}.
}
