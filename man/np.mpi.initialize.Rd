\name{np.mpi.initialize}
\alias{np.mpi.initialize}
\title{Initialize Master and Slave Nodes for the np Package}
\description{
  \code{np.mpi.initialize} is used to initialize master and slave nodes.
}
\details{
  For MPI startup/performance guidance (including message-passing tradeoffs and the manual-broadcast template), see \code{\link{npRmpi.init}} details and \code{inst/.Rprofile}.
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
