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
  Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Note}{
  This method currently only supports objects from the \code{\link{np}} library.
}
\examples{
x <- runif(10)
y <- x + rnorm(10, sd = 0.1)
gradients(npreg(y~x, gradients=TRUE))
}
\keyword{nonparametric}
