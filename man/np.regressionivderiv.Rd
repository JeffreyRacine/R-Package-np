\name{npregivderiv}
\alias{npregivderiv}

\title{
Nonparametric Instrumental Derivatives
}
\description{
\code{npregivderiv} uses the approach of Florens and Racine (2012) to
compute the partial derivative of a nonparametric estimation of an
instrumental regression function \eqn{\varphi}{phi} defined by
conditional moment restrictions stemming from a structural econometric
model: \eqn{E [Y - \varphi (Z,X) | W ] = 0}{E [Y - phi (Z,X) | W ] = 0},
and involving endogenous variables \eqn{Y} and \eqn{Z} and exogenous
variables \eqn{X} and instruments \eqn{W}. The derivative function
\eqn{\varphi'}{phi'} is the solution of an ill-posed inverse problem,
and is computed using Landweber-Fridman regularization.
}
\usage{
npregivderiv(y,
             z,
             w,
             x = NULL,
             zeval = NULL,
             weval = NULL,
             xeval = NULL,
             p = 1,
             nmulti = 1,
             random.seed = 42,
             optim.maxattempts = 10,
             optim.method = c("Nelder-Mead", "BFGS", "CG"),
             optim.reltol = sqrt(.Machine$double.eps),
             optim.abstol = .Machine$double.eps,
             optim.maxit = 500,
             iterate.max = 1000,
             iterate.diff.tol = 1.0e-08,
             constant = 0.5,
             penalize.iteration = TRUE,
             start.from = c("Eyz","EEywz"),
             starting.values = NULL,
             stop.on.increase = TRUE,
             smooth.residuals = TRUE,
             ...)
}
\arguments{
  \item{y}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{z}.
  }

  \item{z}{
    a \eqn{p}-variate data frame of endogenous regressors. The data
    types may be continuous, discrete (unordered and ordered factors),
    or some combination thereof.
  }
  \item{w}{
    a \eqn{q}-variate data frame of instruments. The data types may be
    continuous, discrete (unordered and ordered factors), or some
    combination thereof.
  }
  \item{x}{
    an \eqn{r}-variate data frame of exogenous regressors. The data
    types may be continuous, discrete (unordered and ordered factors),
    or some combination thereof.
  }
  \item{zeval}{
    a \eqn{p}-variate data frame of endogenous regressors on which the
    regression will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{z}.
  }
  \item{weval}{
    a \eqn{q}-variate data frame of instruments on which the regression
    will be estimated (evaluation data). By default, evaluation
    takes place on the data provided by \code{w}.
  }
  \item{xeval}{
    an \eqn{r}-variate data frame of exogenous regressors on which the
    regression will be estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{x}.
  }
  \item{p}{
    the order of the local polynomial regression (defaults to
    \code{p=1}, i.e. local linear).
  }
  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points.
  }

  \item{random.seed}{
    an integer used to seed R's random number generator. This ensures
    replicability of the numerical search. Defaults to 42.
  }

  \item{optim.method}{ method used by \code{\link{optim}} for minimization of
    the objective function. See \code{?optim} for references. Defaults
    to \code{"Nelder-Mead"}.

    the default method is an implementation of that of Nelder and Mead
    (1965), that uses only function values and is robust but relatively
    slow.  It will work reasonably well for non-differentiable
    functions.

    method \code{"BFGS"} is a quasi-Newton method (also known as a
    variable metric algorithm), specifically that published
    simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno.
    This uses function values and gradients to build up a picture of the
    surface to be optimized.

    method \code{"CG"} is a conjugate gradients method based
    on that by Fletcher and Reeves (1964) (but with the option of
    Polak-Ribiere or Beale-Sorenson updates).  Conjugate gradient
    methods will generally be more fragile than the BFGS method, but as
    they do not store a matrix they may be successful in much larger
    optimization problems.
  }

  \item{optim.maxattempts}{
    maximum number of attempts taken trying to achieve successful
    convergence in \code{\link{optim}}. Defaults to \code{100}.
  }

  \item{optim.abstol}{
    the absolute convergence tolerance used by \code{\link{optim}}. Only useful
    for non-negative functions, as a tolerance for reaching
    zero. Defaults to \code{.Machine$double.eps}.
  }

  \item{optim.reltol}{
    relative convergence tolerance used by \code{\link{optim}}.  The algorithm
    stops if it is unable to reduce the value by a factor of 'reltol *
    (abs(val) + reltol)' at a step.  Defaults to
    \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}.
  }

  \item{optim.maxit}{
    maximum number of iterations used by \code{\link{optim}}. Defaults
     to \code{500}.
   }
\item{iterate.max}{
  an integer indicating the maximum number of iterations permitted
  before termination occurs for Landweber-Fridman iteration.
}
\item{iterate.diff.tol}{
  the search tolerance for the difference in the stopping rule from
  iteration to iteration when using Landweber-Fridman iteration
  (disable by setting to zero).
}
\item{constant}{
  the constant to use for Landweber-Fridman iteration.
}
\item{penalize.iteration}{ a logical value indicating whether to
    penalize the norm by the number of iterations or not (default
    \code{TRUE})
}    
\item{start.from}{ a character string indicating whether to start from
  \eqn{E(Y|z)} (default, \code{"Eyz"}) or from \eqn{E(E(Y|z)|z)} (this can
  be overridden by providing \code{starting.values} below)
}  
\item{starting.values}{ a value indicating whether to commence
  Landweber-Fridman assuming
  \eqn{\varphi'_{-1}=starting.values}{phi'[-1]=starting.values} (proper
  Landweber-Fridman) or instead begin from \eqn{E(y|z)} (defaults to
  \code{NULL}, see details below)
}
\item{stop.on.increase}{
  a logical value (defaults to \code{TRUE}) indicating whether to halt
  iteration if the stopping criterion (see below) increases over the
  course of one iteration (i.e. it may be above the iteration tolerance
  but increased).
}
\item{smooth.residuals}{
  a logical value (defaults to \code{TRUE}) indicating whether to
  optimize bandwidths for the regression of \eqn{y-\varphi(z)}{y-phi(z)}
  on \eqn{w} or for the regression of \eqn{\varphi(z)}{phi(z)} on
  \eqn{w} during Landweber-Fridman iteration.
 }
\item{...}{ 
  additional arguments supplied to \code{\link{npksum}}.
}

}

\details{

  Note that Landweber-Fridman iteration presumes that
  \eqn{\varphi_{-1}=0}{phi[-1]=0}, and so for derivative estimation we
  commence iterating from a model having derivatives all equal to
  zero. Given this starting point it may require a fairly large number
  of iterations in order to converge. Other perhaps more reasonable
  starting values might present themselves. When \code{start.phi.zero}
  is set to \code{FALSE} iteration will commence instead using
  derivatives from the conditional mean model \eqn{E(y|z)}. Should the
  default iteration terminate quickly or you are concerned about your
  results, it would be prudent to verify that this alternative starting
  value produces the same result. Also, check the norm.stop vector for
  any anomalies (such as the error criterion increasing immediately).
  
  Landweber-Fridman iteration uses an optimal stopping rule based upon
  \eqn{||E(y|w)-E(\varphi_k(z,x)|w)||^2 }{||E(y|w)-E(phi(z,x)|w)||^2
  }. However, if insufficient training is conducted the estimates can be
  overly noisy. To best guard against this eventuality set \code{nmulti}
  to a larger number than the default \code{nmulti=0} for
  \code{\link{npreg}}.

  Iteration will terminate when either the change in the value of
  \eqn{||(E(y|w)-E(\varphi_k(z,x)|w))/E(y|w)||^2
  }{||(E(y|w)-E(phi(z,x)|w))/E(y|w)||^2 } from iteration to iteration is
  less than \code{iterate.diff.tol} or we hit \code{iterate.max} or
  \eqn{||(E(y|w)-E(\varphi_k(z,x)|w))/E(y|w)||^2
  }{||(E(y|w)-E(phi(z,x)|w))/E(y|w)||^2 } stops falling in value and
  starts rising.  }

\value{

  \code{npregivderiv} returns a list with components \code{phi.prime},
  \code{phi}, \code{num.iterations}, \code{norm.stop} and
  \code{convergence}.
  
}

\references{
  
  Carrasco, M. and J.P. Florens and E. Renault (2007), \dQuote{Linear
  Inverse Problems in Structural Econometrics Estimation Based on
  Spectral Decomposition and Regularization,} In: James J. Heckman and
  Edward E. Leamer, Editor(s), Handbook of Econometrics, Elsevier, 2007,
  Volume 6, Part 2, Chapter 77, Pages 5633-5751

  Darolles, S. and Y. Fan and J.P. Florens and E. Renault (2011),
  \dQuote{Nonparametric Instrumental Regression,} Econometrica, 79,
  1541-1565.
  
  Feve, F. and J.P. Florens (2010), \dQuote{The Practice of
  Non-parametric Estimation by Solving Inverse Problems: The Example of
  Transformation Models,} Econometrics Journal, 13, S1-S27.

  Florens, J.P. and J.S. Racine (2012), \dQuote{Nonparametric
  Instrumental Derivatives,} Working Paper.

  Fridman, V. M. (1956), \dQuote{A Method of Successive Approximations
  for Fredholm Integral Equations of the First Kind,} Uspeskhi,
  Math. Nauk., 11, 233-334, in Russian.

  Horowitz, J.L. (2011), \dQuote{Applied Nonparametric Instrumental
  Variables Estimation,} Econometrica, 79, 347-394.
  
  Landweber, L. (1951), \dQuote{An Iterative Formula for Fredholm
  Integral Equations of the First Kind,} American Journal of
  Mathematics, 73, 615-24.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated Local Linear
  Nonparametric Regression,} Statistica Sinica, 14, 485-512.

}
\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\note{
  This function currently supports univariate \code{z} only.  This
  function is not MPI-aware. This function should be considered to be in
  `beta test' status until further notice.
}

\seealso{
\code{\link{npregiv},\link{npreg}}
}
\examples{
\dontrun{
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

mpi.spawn.Rslaves(nslaves=1)
mpi.bcast.cmd(np.mpi.initialize(),caller.execute=TRUE)

# ## This illustration was made possible by Samuele Centorrino
# ## <samuele.centorrino@univ-tlse1.fr>
# 
# set.seed(42)
# n <- 1500
# 
# ## For trimming the plot (trim .5\% from each tail)
# 
# trim <- 0.005
# 
# ## The DGP is as follows:
# 
# ## 1) y = phi(z) + u
# 
# ## 2) E(u|z) != 0 (endogeneity present)
# 
# ## 3) Suppose there exists an instrument w such that z = f(w) + v and
# ## E(u|w) = 0
# 
# ## 4) We generate v, w, and generate u such that u and z are
# ## correlated. To achieve this we express u as a function of v (i.e. u =
# ## gamma v + eps)
# 
# v <- rnorm(n,mean=0,sd=0.27)
# eps <- rnorm(n,mean=0,sd=0.05)
# u <- -0.5*v + eps
# w <- rnorm(n,mean=0,sd=1)
# 
# ## In Darolles et al (2011) there exist two DGPs. The first is
# ## phi(z)=z^2 and the second is phi(z)=exp(-abs(z)) (which is
# ## discontinuous and has a kink at zero).
# 
# fun1 <- function(z) { z^2 }
# fun2 <- function(z) { exp(-abs(z)) }
# 
# z <- 0.2*w + v
# 
# ## Generate two y vectors for each function.
# 
# y1 <- fun1(z) + u
# y2 <- fun2(z) + u
# 
# ## You set y to be either y1 or y2 (ditto for phi) depending on which
# ## DGP you are considering:
# 
# y <- y1
# phi <- fun1
# 
# ## Sort on z (for plotting)
# 
# ivdata <- data.frame(y,z,w,u,v)
# ivdata <- ivdata[order(ivdata$z),]
# rm(y,z,w,u,v)
# 
# mpi.bcast.Robj2slave(ivdata)
# mpi.bcast.cmd(attach(ivdata),
#               caller.execute=TRUE)
# 
# mpi.bcast.cmd(model.ivderiv <- npregivderiv(y=y,z=z,w=w),
#               caller.execute=TRUE)
# 
# ylim <-c(quantile(model.ivderiv$phi.prime,trim),
#          quantile(model.ivderiv$phi.prime,1-trim))
# 
# plot(z,model.ivderiv$phi.prime,
#      xlim=quantile(z,c(trim,1-trim)),
#      main="",
#      ylim=ylim,
#      xlab="Z",
#      ylab="Derivative",
#      type="l",
#      lwd=2)
# rug(z)     
# 
## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

mpi.close.Rslaves()

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun
}
\keyword{ instrument }
