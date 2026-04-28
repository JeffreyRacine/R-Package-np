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
    for \code{npregression} objects fitted with \code{regtype="lp"},
    optional derivative order request (scalar or one entry per
    continuous predictor). Orders exceeding the fitted polynomial degree
    (or greater than one, pending future C-level support) are returned
    as \code{NA}.
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
  from objects.
}
\seealso{
  \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{coef}},
  and \code{\link{se}}, for related methods;
  \code{\link{npRmpi}} for supported objects;
  \code{\link{npRmpi.init}} for MPI session startup.
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
  This method currently only supports objects from the \code{\link{npRmpi}} library.
}
\examples{
## Not run in checks: excluded to keep MPI examples stable and check times short.
\dontrun{

force.run <- nzchar(Sys.getenv("NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK"))

in.check <- !force.run && (
            nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
            nzchar(Sys.getenv("_R_CHECK_INTERNALS2_")) ||
            nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_")) ||
            nzchar(Sys.getenv("_R_CHECK_EXAMPLE_TIMING_THRESHOLD_"))
          )

if (!in.check) {
npRmpi.init(nslaves=1)

x <- runif(10)
y <- x + rnorm(10, sd = 0.1)
model <- npreg(y~x, gradients=TRUE)
gradients(model)

npRmpi.quit()
} else {
  message("Skipping MPI spawn in check context (set NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1).")
}
} % enddontrun     
}
\keyword{nonparametric}
