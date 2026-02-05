\name{Engel95}
\docType{data}
\alias{Engel95}
\title{ 1995 British Family Expenditure Survey  }
\description{
British cross-section data consisting of a random sample taken from
the British Family Expenditure Survey for 1995. The households consist
of married couples with an employed head-of-household between the ages
of 25 and 55 years. There are 1655 household-level observations in
total.
}
\usage{data("Engel95")}
\format{ A data frame with 10 columns, and 1655 rows.
  \describe{
 \item{food}{ expenditure share on food, of type \code{numeric}}
 \item{catering}{ expenditure share on catering, of type \code{numeric}}
 \item{alcohol}{ expenditure share on alcohol, of type \code{numeric}}
 \item{fuel}{ expenditure share on fuel, of type \code{numeric}}
 \item{motor}{ expenditure share on motor, of type \code{numeric}}
 \item{fares}{ expenditure share on fares, of type \code{numeric}}
 \item{leisure}{ expenditure share on leisure, of type \code{numeric}}
 \item{logexp}{ logarithm of total expenditure, of type \code{numeric}}
 \item{logwages}{ logarithm of total earnings, of type \code{numeric}}
 \item{nkids}{ number of children, of type \code{numeric}}  
  }
}
\source{ Richard Blundell and Dennis Kristensen }
\references{

  Blundell, R. and X. Chen and D. Kristensen (2007),
  \dQuote{Semi-Nonparametric IV Estimation of Shape-Invariant Engel
  Curves,} Econometrica, 75, 1613-1669.
    
  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.
  
}

\keyword{datasets}

\examples{

data(Engel95)
summary(Engel95)

}
