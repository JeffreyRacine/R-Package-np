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
  Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Note}{
  This method currently only supports objects from the \code{\link{np}} library.
}
\examples{
x <- rnorm(10)
se(npudens(npudensbw(~x)))
}
\keyword{nonparametric}
