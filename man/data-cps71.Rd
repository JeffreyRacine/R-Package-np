% $ID$
\name{cps71}
\docType{data}
\alias{cps71}
\title{ Canadian High School Graduate Earnings  }
\description{
	Canadian cross-section wage data consisting of a random sample taken
  from the 1971 Canadian Census Public Use Tapes for male individuals
  having common education (grade 13). There are 205 observations in total.
}
\usage{data("cps71")}
\format{ A data frame with 2 columns, and 205 rows.
  \describe{
    \item{logwage}{ the first column, of type \code{numeric}}
    \item{age}{ the second column, of type \code{integer}}
  }
}
\source{ Aman Ullah }
\examples{
data("cps71")
attach(cps71)

plot(age, logwage, xlab="Age", ylab="log(wage)")

detach(cps71)
}
\references{
  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge	University Press.
}
\keyword{datasets}
