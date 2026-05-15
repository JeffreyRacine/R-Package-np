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

\method{gradients}{condensity}(x, errors = FALSE, gradient.order = NULL, \dots)

\method{gradients}{condistribution}(x, errors = FALSE, gradient.order = NULL, \dots)

\method{gradients}{npregression}(x, errors = FALSE, gradient.order = NULL, \dots)

\method{gradients}{qregression}(x, errors = FALSE, \dots)

\method{gradients}{singleindex}(x, errors = FALSE, \dots)
}

\arguments{
\subsection{Object And Output Controls}{
Object to interrogate and whether gradient standard errors are requested.
}

  \item{x}{an object for which the extraction of gradients is
    meaningful.}
  \item{errors}{ a logical value specifying whether or not standard
    errors of gradients are desired. Defaults to \code{FALSE}.}

\subsection{Derivative Order Controls}{
Optional local-polynomial derivative order controls.
}

  \item{gradient.order}{
    for \code{npregression}, \code{condensity}, and
    \code{condistribution} objects fitted with \code{regtype="lp"},
    optional derivative order request (scalar or one entry per
    continuous predictor). Orders exceeding the fitted polynomial degree
    are returned as \code{NA}. If a different admissible derivative order
    is desired, compute the fit, prediction, or evaluation with
    \code{gradients=TRUE} and the desired \code{gradient.order};
    \code{gradients()} extracts stored derivatives and does not
    recompute another derivative order after fitting.
  }

\subsection{Additional Arguments}{
Further method-specific arguments.
}

  \item{\dots}{other arguments.}
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
  from objects. For \code{npregression}, \code{condensity}, and
  \code{condistribution} objects fitted with \code{regtype="lp"},
  \code{gradient.order} identifies the stored continuous-predictor
  derivative order to extract; it is not a post-fit recomputation
  control. For \code{qregression} objects,
  \code{errors=TRUE} returns
  asymptotic standard errors for the quantile gradients when the object was
  fitted with \code{gradients=TRUE}.
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
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
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
