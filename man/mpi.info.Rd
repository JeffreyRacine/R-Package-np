\name{mpi.info.create}
\alias{mpi.info.create}
\alias{mpi.info.free}
\alias{mpi.info.get}
\alias{mpi.info.set}
\title{MPI\_Info\_create, MPI\_Info\_free, MPI\_Info\_get, MPI\_Info\_set APIs}
\description{
Many MPI APIs take an info argument for additional information passing. An info 
is an object which consists of many (key,value) pairs. Rmpi uses an internal 
memory to store an info object.

\code{mpi.info.create} creates a new info object.

\code{mpi.info.free} frees an info object and sets it to MPI\_INFO\_NULL.

\code{mpi.info.get}  retrieves the value associated with key in an info.

\code{mpi.info.set} adds the key and value pair to info.
}
\usage{
  mpi.info.create(info = 0)
  mpi.info.free(info = 0)
  mpi.info.get(info = 0, key, valuelen)
  mpi.info.set(info = 0, key, value)
}
\arguments{
  \item{info}{an info number.}
  \item{key}{a char (length 1).}
  \item{valuelen}{the length (nchar) of a key}
  \item{value}{a char (length 1).}
}
\value{
  \code{mpi.info.create}, \code{mpi.info.free}, and \code{mpi.info.set}  return 
1 if success and 0 otherwise.

\code{mpi.info.get} returns the value (a char) for a given info and valuelen.
}
\seealso{
  \code{\link{mpi.spawn.Rslaves}}
}
\author{
Hao Yu
}
\keyword{interface}
