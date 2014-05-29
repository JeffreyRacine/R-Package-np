\name{se}
\alias{se}
\title{Extract Standard Errors}
\usage{
se(x)
}
\arguments{
  \item{x}{an object for which the extraction of standard errors is
    meaningful.}
  \item{\dots}{other arguments.}
}
\value{
  Standard errors extracted from the model object \code{x}.
}
\description{
  \code{se} is a generic function which extracts standard errors 
  from objects. 
}
\details{
  This function provides a generic interface for extraction of standard
  errors from objects.
}
\seealso{
  \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{coef}},
  and \code{\link{gradients}}, for related methods;
  \code{\link{np}} for supported objects. 
}

\author{
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Note}{
  This method currently only supports objects from the \code{\link{np}} library.
}
\examples{
\dontrun{
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

mpi.spawn.Rslaves(nslaves=1)
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

x <- rnorm(10)
mpi.bcast.Robj2slave(x)
mpi.bcast.cmd(se(npudens(npudensbw(~x))),
              caller.execute=TRUE)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

mpi.close.Rslaves()

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun     
}
\keyword{nonparametric}
