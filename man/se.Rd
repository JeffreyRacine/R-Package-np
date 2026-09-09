\name{se}
\alias{se}
\title{Extract Standard Errors}
\usage{
se(x)
}
\arguments{
\subsection{Object Input}{
Object to interrogate for standard errors.
}

  \item{x}{an object for which the extraction of standard errors is
    meaningful.}
}
\value{
  Standard errors extracted from the model object \code{x}.
}
\description{
  \code{se} is a generic function which extracts standard errors 
  from objects. 
}
\details{
  This function extracts stored standard errors; it does not select
  bandwidths, refit a model, or run a bootstrap. For estimators with a
  \code{se} fitting control, request \code{se = TRUE} when those errors
  will be needed. With \code{se = FALSE}, errors are deliberately absent
  and this accessor explains how to refit with \code{se = TRUE}, reusing
  the bandwidth object without repeating bandwidth selection.
  Requested but unavailable standard errors may remain \code{NA};
  extraction does not replace or recompute them.
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
x <- rnorm(10)
se(npudens(npudensbw(~x), se = TRUE))
}
\keyword{nonparametric}
