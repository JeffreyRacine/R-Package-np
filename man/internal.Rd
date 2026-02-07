\name{npRmpi-internal}
\alias{string}
\alias{mpi.comm.is.null}
\alias{.docall}
\alias{.force.type}
\alias{.mpi.undefined}
\alias{.mpi.worker.apply}
\alias{.mpi.worker.applyLB}
\alias{.mpi.worker.exec}
\alias{.mpi.worker.sim}
\alias{.simplify}
\alias{.splitIndices}
\alias{.typeindex}
\title{Internal npRmpi functions}
\description{
  Internal functions used by other MPI functions. These are not intended to be called directly by the user.
}
\usage{
mpi.comm.is.null(comm)
string(length)
.docall(fun, args)
.force.type(x, type)
.mpi.undefined()
.mpi.worker.apply(n, tag)
.mpi.worker.applyLB(n)
.mpi.worker.exec(tag, ret, simplify)
.mpi.worker.sim(n, nsim, run)
.simplify(n, answer, simplify, len = 1, recursive = FALSE)
.splitIndices(nx, ncl)
.typeindex(x)
}
\arguments{
  \item{comm}{a communicator number.}
  \item{length}{length of a string.}
  \item{fun}{a function or name of a function.}
  \item{args}{a list of arguments.}
  \item{x}{an object.}
  \item{type}{a type indicator.}
  \item{n}{number of tasks.}
  \item{tag}{an MPI tag.}
  \item{ret}{logical; whether to return a value.}
  \item{simplify}{logical; whether to simplify the result.}
  \item{nsim}{number of simulations.}
  \item{run}{run indicator.}
  \item{answer}{a result list.}
  \item{len}{expected length.}
  \item{recursive}{logical; whether to unlist recursively.}
  \item{nx}{number of elements.}
  \item{ncl}{number of clusters.}
}
\details{
  These functions are required for internal MPI communication and slave execution.
}
\author{
Hao Yu and Jeffrey Racine
}
\keyword{internal}