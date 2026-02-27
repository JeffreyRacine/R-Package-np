\name{mpi.get.version}
\title{MPI_Get_version API}
\usage{
  mpi.get.version()
}
\alias{mpi.get.version}
\description{
  \code{mpi.get.version} returns the runtime MPI API version as reported by
  \code{MPI_Get_version}.
}
\value{
An integer vector of length two named \code{major} and \code{minor}.
}
\author{
Hao Yu
}
\references{
  \url{https://www.mpich.org/},
  \url{https://www.mpich.org/static/docs/latest/www3/}
}
\keyword{interface}
