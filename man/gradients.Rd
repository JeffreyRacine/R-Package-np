\name{gradients}
\title{Extract Gradients}
\alias{gradients}
\alias{gradients.condensity}
\alias{gradients.condistribution}
\alias{gradients.npregression}
\alias{gradients.qregression}  
\alias{gradients.singleindex}
\usage{
gradients(x, \dots)

\method{gradients}{condensity}(x, errors = FALSE, \dots)

\method{gradients}{condistribution}(x, errors = FALSE, \dots)

\method{gradients}{npregression}(x, errors = FALSE, \dots)

\method{gradients}{qregression}(x, errors = FALSE, \dots)

\method{gradients}{singleindex}(x, errors = FALSE, \dots)
}

\arguments{
  \item{x}{an object for which the extraction of gradients is
    meaningful.}
  \item{\dots}{other arguments.}
  \item{errors}{ a logical value specifying whether or not standard
    errors of gradients are desired. Defaults to \code{FALSE}.}
}
\value{
  Gradients extracted from the model object \code{x}.
}
\description{
  \code{gradients} is a generic function which extracts gradients 
  from objects. 
}
\details{
  This function provides a generic interface for extraction of gradients
  from objects.
}
\seealso{
  \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{coef}},
  and \code{\link{se}}, for related methods;
  \code{\link{np}} for supported objects. 
}
\references{
  See the references for the method being interrogated via
  \code{\link{gradients}} in the appropriate help file. For example, for
  the particulars of the gradients for nonparametric regression see the
  references in \code{\link{npreg}}
}
\author{
  Tristen Hayfield \email{hayfield@mpia.de}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Note}{ This method currently only supports objects from the
  \code{\link{np}} library.  }

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

x <- runif(10)
y <- x + rnorm(10, sd = 0.1)

mpi.bcast.Robj2slave(x)
mpi.bcast.Robj2slave(y)
mpi.bcast.cmd(gradients(npreg(y~x, gradients=TRUE)),
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
