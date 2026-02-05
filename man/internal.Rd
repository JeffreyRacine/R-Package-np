\name{string}
\alias{string}
\alias{mpi.comm.is.null}
\title{Internal functions}
\description{ Internal functions used by other MPI functions.

\code{mpi.comm.is.null} is used to test if a comm is MPI\_COMM\_NULL (empty 
members).

\code{string} create a string (empty space character) buffer.
}

\usage{
mpi.comm.is.null(comm)
string(length)
}
\arguments{
  \item{comm}{a communicator number.}
  \item{length}{length of a string.}
}
\value{
\code{string} returns an empty character string.
}
\seealso{
  \code{\link{mpi.spawn.Rslaves}}
}
\author{
Hao Yu
}
\keyword{utilities}

