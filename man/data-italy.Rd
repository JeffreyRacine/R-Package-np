% $ID$
\name{Italy}
\docType{data}
\alias{Italy}
\title{ Italian GDP Panel }
\description{
  Italian GDP growth panel for 21 regions covering the period 1951-1998
  (millions of Lire, 1990=base). There are 1008 observations in total.
}
\usage{data("Italy")}
\format{ A data frame with 2 columns, and 1008 rows.
  \describe{
    \item{year}{ the first column, of type \code{ordered}}

    \item{gdp}{ the second column, of type \code{numeric}: millions of
    Lire, 1990=base}
  }
}
\source{Giovanni Baiocchi}
\examples{
data("Italy")
attach(Italy)

plot(ordered(year), gdp, xlab="Year (ordered factor)",
     ylab="GDP (millions of Lire, 1990=base)")

detach(Italy)
}
\references{
  Baiocchi, G. (2006), \dQuote{Economic Applications of Nonparametric
    Methods,} Ph.D. Thesis, University of York.
}
\keyword{datasets}
