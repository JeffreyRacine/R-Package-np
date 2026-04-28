\name{dimBS}
\alias{dimBS}
\title{Local-Polynomial Basis Dimension Helper}
\description{
  \code{dimBS} returns the number of columns implied by an additive,
  generalized local-polynomial, or tensor-product basis specification.
  It is a compatibility wrapper around the internal \code{dim_basis()}
  helper used by \pkg{npRmpi}.
}
\usage{
dimBS(basis = "additive",
      kernel = TRUE,
      degree = NULL,
      segments = NULL,
      include = NULL,
      categories = NULL)
}
\arguments{
\subsection{Basis Specification}{
Basis family, continuous-kernel counting mode, polynomial degree, and
segment controls.
}

  \item{basis}{
    basis family. One of \code{"additive"}, \code{"glp"}, or
    \code{"tensor"}.
  }
  \item{kernel}{
    logical indicating whether only the continuous-kernel basis should
    be counted. When \code{FALSE}, optional categorical augmentation
    controlled by \code{include} and \code{categories} is included.
  }
  \item{degree}{
    non-negative integer vector of local-polynomial degrees.
  }
  \item{segments}{
    positive integer vector giving the number of segments for each
    continuous predictor. Defaults to one segment per degree entry.
  }

\subsection{Categorical Augmentation}{
Optional categorical-component controls used when \code{kernel = FALSE}.
}

  \item{include}{
    non-negative integer vector indicating which categorical
    components are included when \code{kernel = FALSE}.
  }
  \item{categories}{
    non-negative integer vector giving category counts for included
    categorical components when \code{kernel = FALSE}.
  }
}
\value{
  A numeric scalar giving the implied basis dimension.
}
\details{
  \code{dimBS()} is provided for compatibility with \pkg{crs}. In
  \pkg{npRmpi}, the underlying implementation lives in
  \code{dim_basis()}, which is used internally for LP basis-dimension
  checks and safe NOMAD restart initialization.
}
\examples{
dimBS(basis = "tensor", degree = c(2, 2))
dimBS(basis = "glp", degree = c(3, 1, 0))
}
