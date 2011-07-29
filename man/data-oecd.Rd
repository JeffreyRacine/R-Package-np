% $ID$
\name{oecdpanel}
\docType{data}
\alias{oecdpanel}
\alias{bw}
\title{ Cross Country Growth Panel }
\description{
  Cross country GDP growth panel covering the period 1960-1995 used by
  Liu and Stengos (2000) and Maasoumi, Racine, and Stengos (2007). There
  are 616 observations in total. \code{data("oecdpanel")} makes available
  the dataset \code{"oecdpanel"} plus an additional object \code{"bw"}.
}
\usage{data("oecdpanel")}
\format{

  A data frame with 7 columns, and 616 rows.  This panel covers 7
  5-year periods: 1960-1964, 1965-1969, 1970-1974, 1975-1979, 1980-1984,
  1985-1989 and 1990-1994.

  A separate local-linear \code{rbandwidth} object (`\code{bw}') has
  been computed for the user's convenience which can be used to
  visualize this dataset using \code{\link{plot}(bw)}.
  
  \describe{

    \item{growth}{ the first column, of type \code{numeric}: growth rate
    of real GDP per capita for each 5-year period}

    \item{oecd}{ the second column, of type \code{factor}: equal to
    1 for OECD members, 0 otherwise}

    \item{year}{ the third column, of type \code{integer}}

    \item{initgdp}{ the fourth column, of type \code{numeric}: per capita
    real GDP at the beginning of each 5-year period}

    \item{popgro}{ the fifth column, of type \code{numeric}: average
    annual population growth rate for each 5-year period}

    \item{inv}{ the sixth column, of type \code{numeric}: average
    investment/GDP ratio for each 5-year period}

    \item{humancap}{ the seventh column, of type \code{numeric}:
      average secondary school enrollment rate for each 5-year period}

    }

}
\source{ Thanasis Stengos }
\examples{
data("oecdpanel")
attach(oecdpanel)
summary(oecdpanel)
detach(oecdpanel)
}
\references{

  Liu, Z. and T. Stengos (1999), \dQuote{Non-linearities in cross
  country growth regressions: a semiparametric approach,}
  Journal of Applied Econometrics, 14, 527-538.

  Maasoumi, E. and J.S. Racine and T. Stengos (2007),
  \dQuote{Growth and convergence: a profile of distribution dynamics and
  mobility,} Journal of Econometrics, 136, 483-508 
  
}
\keyword{datasets}
