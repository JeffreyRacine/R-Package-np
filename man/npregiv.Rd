\name{npregiv}
\alias{npregiv}

\title{
Nonparametric Instrumental Regression
}
\description{
\code{npregiv} computes nonparametric estimation of an instrumental
regression function \eqn{\varphi}{phi} defined by conditional moment
restrictions stemming from a structural econometric model: \eqn{E [Y -
\varphi (Z) | W ] = 0}{E [Y - phi (Z) | W ] = 0}, and involving
endogenous variables \eqn{Y} and \eqn{Z} and instruments \eqn{W}. The
function \eqn{\varphi}{phi} is the solution of an ill-posed inverse
problem and the estimation procedure is based on Tikhonov
regularization.

\code{npregiv} is an implementation of the approach of S. Darolles,
Y. Fan, J.P. Florens and E. Renault (forthcoming, Econometrica).
}
\usage{
npregiv(y,
        z,
        w,
        yeval = NULL,
        zeval = NULL,
        weval = NULL,
        alpha.min = 1e-10,
        alpha.max = 1,
        p = 1)
}

\arguments{
  \item{y}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{z}.
  }

  \item{z}{
    a one (1) dimensional numeric or integer vector for the endogenous
    regressor.
  }
  \item{w}{
    a one (1) dimensional numeric or integer vector for the instrument.
}
  \item{yeval}{
    a one (1) dimensional numeric or integer vector on which the
    regression will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{y}.
}
  \item{zeval}{
    a one (1) dimensional numeric or integer vector on which the
    regression will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{z}.
}
  \item{weval}{
    a one (1) dimensional numeric or integer vector on which the
    regression will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{w}.
}
  \item{alpha.min}{
    mimimum of search range for \eqn{\alpha}{alpha}, the Tikhonov
    regularization parameter.
}
  \item{alpha.max}{
    maximum of search range for \eqn{\alpha}{alpha}, the Tikhonov
    regularization parameter.
}
  \item{p}{
    the order of the local polynomial regression (defaults to
    \code{p=1}, i.e. local linear)
}
}
\details{
  TBA
}
\value{
  \code{npregiv} returns a list with components \code{phihat} and
    \code{alpha}.
}
\references{
  
  Darolles, S. and Y. Fan and J.P. Florens and E. Renault (forthcoming),
  \dQuote{Nonparametric Instrumental Regression,} Econometrica.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
  nonparametric regression,} Statistica Sinica, 14, 485-512.

}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}, Samuele Centorrino
 \email{samuele.centorrino@univ-tlse1.fr}
}
\note{
This function should be considered to be in `beta test' status until further notice.
}

\seealso{
\code{\link{npreg}}
}
\examples{
## This illustration was made possible by Samuele Centorrino
## <samuele.centorrino@univ-tlse1.fr>

set.seed(42)
n <- 1500

## The DGP is as follows:

## 1) y = phi(z) + u

## 2) E(z|u) !=0 hence we have endogeneity

## 3) Suppose there exists an instrument w such that z = beta w + v

## 4) We generate v, w, and generate u such that u and z are
## correlated. To achieve this we express u as a function of v (i.e. u =
## gamma v + eps)

v <- rnorm(n,mean=0,sd=0.27)
eps <- rnorm(n,mean=0,sd=0.05)
u <- -0.5*v + eps
w <- rnorm(n,mean=0,sd=1)

# In Darolles et al (2011) there exist two DGPs. The first is
# phi(z)=z^2 and the second is phi(z)=exp(-abs(z)) (which is
# discontinuous and has a kink at zero).

fun1 <- function(z) { z^2 }
fun2 <- function(z) { exp(-abs(z)) }

z <- 0.4*w + v

## Generate two y vectors for each function.

y1 <- fun1(z) + u
y2 <- fun2(z) + u

## You set y to be either y1 or y2 depending on which DGP you are
## considering:

y <- y1
phi <- fun1

## sort on z

ivdata <- data.frame(y,z,w)
ivdata <- ivdata[order(ivdata$z),]
rm(y,z,w)
attach(ivdata)

phihat.iv <- npregiv(y=y,z=z,w=w)$phihat

## Now the non-iv local linear estimator of E(y|z)

ll.mean <- fitted(npreg(y~z,regtype="ll"))

## For the plots, restrict focal attention to the bulk of the data
## (i.e. for the plotting area trim out 1/4 of one percent from each
## tail of y and z)

trim <- 0.0025

curve(phi,min(z),max(z),
     xlim=quantile(z,c(trim,1-trim)),
     ylim=quantile(y,c(trim,1-trim)),
      ylab="Y",
      xlab="Z",
      main="Nonparametric IV (Darolles, Fan, Florens and Renault (2011))",
      lwd=2,lty=1)

points(z,y,type="p",cex=.25,col="grey")

lines(z,phihat.iv,col="blue",lwd=2,lty=2)

lines(z,ll.mean,col="red",lwd=2,lty=4)

legend(quantile(z,trim),quantile(y,1-trim),
       c(expression(paste(varphi(z))),
         expression(paste("Nonparametric ",hat(varphi)(z))),
         "Nonparametric E(y|z)"),
       lty=c(1,2,4),
       col=c("black","blue","red"),
       lwd=c(2,2,2))

}
\keyword{ instrument }
