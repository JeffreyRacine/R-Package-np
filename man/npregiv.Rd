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
polynomial kernel regression of any order (Darolles et al use local
constant kernel weighting which corresponds to setting \code{p=0}; see
below for details). When \code{method="Landweber-Fridman"},
\code{npregiv} uses the approach of Horowitz (2011) again using local
polynomial kernel regression (Horowitz uses B-spline weighting).
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
        start.iterations = 10,
        max.iterations = 100,
        iterate.smoothing = TRUE,
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
    a \eqn{q}-variate data frame of instruments. The data types may be
    continuous, discrete (unordered and ordered factors), or some
    combination thereof.
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
    a \eqn{q}-variate data frame of instruments on which the regression
    will be estimated (evaluation data). By default, evaluation
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
  \item{start.iterations}{
    an integer indicating the minimum number of iterations to conduct
    when using Landweber-Fridman regularization (here at each iteration
    the smoothing parameters are refined and optimized).
  }
  \item{max.iterations}{
    an integer indicating the maximum number of iterations permitted
    before termination occurs.
  }
  \item{iterate.smoothing}{
    a logical value indicating whether, after the initial
    \code{start.iterations} have been conducted, to continue to optimize the
    smoothing parameters for each successive iteration.
  }
  \item{constant}{
    the constant to use when using Landweber-Fridman regularization.
  }
  \item{method}{
  the regularization method employed (defaults to
  \code{"Landweber-Fridman"}, see Horowitz (2011) for details; see
  Darolles, Fan, Florens and Renault (forthcoming) for details on
  \code{"Tikhonov"}).
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
  
  \code{method="Landweber-Fridman"} uses an optimal stopping rule
  whereby if \eqn{||E(y|w)-E(\varphi_k(z)|w)||^2
  }{||E(y|w)-E(phi(z)|w)||^2 } is unchanged for at least
  \code{num.iterations} then iteration terminates. However, if
  insufficient training is conducted the estimates can be overly noisy
  and iteration may fail to terminate when \code{num.iterations} is
  large. To best guard against this eventuality set \code{nmulti} to a
  larger number than the default \code{nmulti=0} for
  \code{\link{npreg}}.

  When using \code{method="Landweber-Fridman"}, iteration will terminate
  when either \code{start.iterations} have passed with no change in the
  value of \eqn{||E(y|w)-E(\varphi_k(z)|w)||^2
  }{||E(y|w)-E(phi(z)|w)||^2 } or we hit \code{max.iterations} or
  \eqn{||E(y|w)-E(\varphi_k(z)|w)||^2 }{||E(y|w)-E(phi(z)|w)||^2 } stops
  falling in value and starts rising.
}
\value{
  \code{crsiv} returns a list with components \code{phihat} and either
  \code{alpha} when \code{method="Tikhonov"} or \code{num.iterations}
  when \code{method="Landweber-Fridman"} and \code{norm.stop}.
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

## 2) E(u|z) != 0 (endogeneity present)

## 3) Suppose there exists an instrument w such that z = f(w) + v and
## E(u|w) = 0

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

## You set y to be either y1 or y2 (ditto for phi) depending on which
## DGP you are considering:

y <- y1
phi <- fun1

## Sort on z (for plotting)

ivdata <- data.frame(y,z,w)
ivdata <- ivdata[order(ivdata$z),]
rm(y,z,w)
attach(ivdata)

model.iv <- npregiv(y=y,z=z,w=w)
phihat.iv <- model.iv$phihat

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
      main="Nonparametric Instrumental Kernel Regression",
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
