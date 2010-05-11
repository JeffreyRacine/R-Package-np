% $ID$
\name{wage1}
\docType{data}
\alias{wage1}
\alias{bw.subset}
\alias{bw.all}
\title{ Cross-Sectional Data on Wages }
\description{ Cross-section wage data consisting of a random sample
taken from the U.S. Current Population Survey for the year 1976. There
are 526 observations in total. \code{data("wage1")} makes available the
dataset \code{"wage"} plus additional objects \code{"bw.all"} and
\code{"bw.subset"}.}

\usage{data("wage1")}

\format{ A data frame with 24 columns, and 526 rows.

    Two local-linear \code{rbandwidth} objects (`\code{bw.all}' and
    `\code{bw.subset}') have been computed for the user's convenience
    which can be used to visualize this dataset using
    \code{\link{npplot}(bws=bw.all)}
  
\describe{
    \item{wage}{column 1, of type \code{numeric}, average hourly earnings}
    \item{educ}{column 2, of type \code{numeric}, years of education}
    \item{exper}{column 3, of type \code{numeric}, years potential experience}
    \item{tenure}{column 4, of type \code{numeric}, years with current employer}
    \item{nonwhite}{column 5, of type \code{factor}, =\dQuote{Nonwhite} if nonwhite, \dQuote{White} otherwise}
    \item{female}{column 6, of type \code{factor}, =\dQuote{Female} if female, \dQuote{Male} otherwise}
    \item{married}{column 7, of type \code{factor}, =\dQuote{Married} if Married, \dQuote{Nonmarried} otherwise}
    \item{numdep}{column 8, of type \code{numeric}, number of dependents}
    \item{smsa}{column 9, of type \code{numeric}, =1 if live in SMSA}
    \item{northcen}{column 10, of type \code{numeric}, =1 if live in north central U.S}
    \item{south}{column 11, of type \code{numeric}, =1 if live in southern region}
    \item{west}{column 12, of type \code{numeric}, =1 if live in western region}
    \item{construc}{column 13, of type \code{numeric}, =1 if work in construc. indus.}
    \item{ndurman}{column 14, of type \code{numeric}, =1 if in nondur. manuf. indus.}
    \item{trcommpu}{column 15, of type \code{numeric}, =1 if in trans, commun, pub ut}
    \item{trade}{column 16, of type \code{numeric}, =1 if in wholesale or retail}
    \item{services}{column 17, of type \code{numeric}, =1 if in services indus.}
    \item{profserv}{column 18, of type \code{numeric}, =1 if in prof. serv. indus.}
    \item{profocc}{column 19, of type \code{numeric}, =1 if in profess. occupation}
    \item{clerocc}{column 20, of type \code{numeric}, =1 if in clerical occupation}
    \item{servocc}{column 21, of type \code{numeric}, =1 if in service occupation}
    \item{lwage}{column 22, of type \code{numeric}, log(wage)}
    \item{expersq}{column 23, of type \code{numeric}, exper\eqn{^2}{**2}}
    \item{tenursq}{column 24, of type \code{numeric}, tenure\eqn{^2}{**2}}
  }
}
\source{ Jeffrey M. Wooldridge }
\examples{
data("wage1")
attach(wage1)
summary(wage1)
detach(wage1)
}


\references{

  Wooldridge, J.M. (2000), \emph{Introductory Econometrics: A Modern
  Approach}, South-Western College Publishing.

}

\keyword{datasets}
