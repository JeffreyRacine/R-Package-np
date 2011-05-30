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
problem.

When \code{method="Tikhonov"}, \code{npregiv} uses the approach of
Darolles, Fan, Florens and Renault (forthcoming) modified for local
polynomial kernel regression of any order (Darolles et al use local constant
kernel weighting). When \code{method="Landweber-Fridman"},
\code{npregiv} uses the approach of Horowitz (2011) again using local
polynomial kernel regresion.
}
\usage{
npregiv(y,
        z,
        w,
        yeval = NULL,
        zeval = NULL,
        weval = NULL,
        p = 1,
        alpha.min = 1e-10,
        alpha.max = 1,
        tol = .Machine$double.eps^0.25,
        num.iterations = 10,
        constant = 0.5,
        method = c("Landweber-Fridman","Tikhonov"),
        ...)
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
  \item{p}{
    the order of the local polynomial regression (defaults to
    \code{p=1}, i.e. local linear)
  }
  \item{alpha.min}{
    minimum of search range for \eqn{\alpha}{alpha}, the Tikhonov
    regularization parameter, when using
    Tikhonov regularization.
  }
  \item{alpha.max}{
    maximum of search range for \eqn{\alpha}{alpha}, the Tikhonov
    regularization parameter, when using
    Tikhonov regularization.
  }
  \item{tol}{
    the search tolerance for \code{optimize} when solving for
    \eqn{\alpha}{alpha}, the Tikhonov regularization parameter, when using
    Tikhonov regularization.
  }
  \item{num.iterations}{
    an integer indicating the number of iterations to conduct when using
    Landweber-Fridman regularization.
  }
  \item{constant}{
    the constant to use when using Landweber-Fridman regularization.
  }
  \item{method}{
    the regularization method employed (defaults to \code{"Tikhonov"}, see
    Darolles, Fan, Florens and Renault (forthcoming); see Horowitz
    (2011) for details for \code{"Landweber-Fridman"}).
  }
  \item{...}{
    additional arguments supplied to \code{crs}.
  }
}
\details{
  Tikhonov regularization requires computation of weight matrices of
  dimension \eqn{n\times n}{n x n} which can be computationally costly
  in terms of memory requirements and may be unsuitable for large
  datasets. Landweber-Fridman will be preferred in such settings as it
  does not require construction and storage of these weight matrices
  while it also avoids the need for numerical optimization methods to
  determine \eqn{\alpha}{alpha}.
}
\value{
  \code{npregiv} returns a list with components \code{phihat} and
    \code{alpha}.
}
\references{
  
  Carrasco, M. and J.P. Florens and E. Renault (2007), \dQuote{Linear
  Inverse Problems in Structural Econometrics Estimation Based on
  Spectral Decomposition and Regularization,} In: James J. Heckman and
  Edward E. Leamer, Editor(s), Handbook of Econometrics, Elsevier, 2007,
  Volume 6, Part 2, Chapter 77, Pages 5633-5751

  Darolles, S. and Y. Fan and J.P. Florens and E. Renault (forthcoming),
  \dQuote{Nonparametric Instrumental Regression,} Econometrica.

  Feve, F. and J.P. Florens (2010), \dQuote{The practice of
  non-parametric estimation by solving inverse problems: the example of
  transformation models,} Econometrics Journal, 13, S1-S27.

  Fridman, V. M. (1956), \dQuote{A Method of Successive Approximations
  for Fredholm Integral Equations of the First Kind,} Uspeskhi,
  Math. Nauk., 11, 233-334, in Russian.

  Horowitz, J.L. (2011), \dQuote{Applied Nonparametric Instrumental
  Variables Estimation,}, Econometrica, 79, 347-394.
  
  Landweber, L. (1951), \dQuote{An iterative formula for Fredholm
  integral equations of the first kind,} American Journal of
  Mathematics, 73, 615-24.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

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

## In Darolles et al (2011) there exist two DGPs. The first is
## phi(z)=z^2 and the second is phi(z)=exp(-abs(z)) (which is
## discontinuous and has a kink at zero).

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
