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
\value{
Internal helpers; return values vary by function:
\itemize{
  \item \code{mpi.comm.is.null}: logical indicator.
  \item \code{string}: character string of requested length.
  \item \code{.docall}: result of calling \code{fun} with \code{args}.
  \item \code{.force.type}: coerced object of the requested type.
  \item \code{.mpi.undefined}: integer constant used by MPI.
  \item \code{.mpi.worker.apply}, \code{.mpi.worker.applyLB},
        \code{.mpi.worker.exec}, \code{.mpi.worker.sim}: internal worker
        results (typically lists or vectors).
  \item \code{.simplify}: simplified result (vector, matrix, or list).
  \item \code{.splitIndices}: list of index vectors.
  \item \code{.typeindex}: integer type code.
}
}
\details{
  These functions are required for internal MPI communication and slave execution.
}
\author{
Hao Yu and Jeffrey Racine
}
\keyword{internal}
