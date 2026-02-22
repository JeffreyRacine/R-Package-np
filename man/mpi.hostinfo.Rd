\name{mpi.hostinfo}
\alias{mpi.hostinfo}
\alias{slave.hostinfo}
\title{Host Information Utilities}
\description{
\code{mpi.hostinfo} prints host and rank information for the calling rank.
\code{slave.hostinfo} prints host/rank summaries for active slave ranks.
}
\usage{
mpi.hostinfo(comm = 1)
slave.hostinfo(comm = 1, short = TRUE)
}
\arguments{
  \item{comm}{a communicator number.}
  \item{short}{a logical; if \code{TRUE}, abbreviate output when there are many slaves.}
}
\details{
\code{slave.hostinfo()} must be called on rank 0 of the target communicator.
}
\value{
Both functions print informational output and return invisibly.
}
\seealso{
  \code{\link{mpi.get.processor.name}},
  \code{\link{mpi.comm.size}},
  \code{\link{mpi.comm.rank}}.
}
\author{
Hao Yu
}
\keyword{utilities}
