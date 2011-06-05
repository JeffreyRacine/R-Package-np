% $ID$
\name{npreg}
\alias{npreg}
\alias{npreg.call}
\alias{npreg.default}
\alias{npreg.formula}
\alias{npreg.rbandwidth}

\title{Kernel Regression with Mixed Data Types}

\description{
  \code{npreg} computes a kernel regression estimate of a one
  (1) dimensional dependent variable on \eqn{p}-variate explanatory
  data, given a set of evaluation points, training points (consisting of
  explanatory data and dependent data), and a bandwidth specification
  using the method of Racine and Li (2004) and Li and Racine (2004). A
  bandwidth specification can be a \code{rbandwidth} object, or a
  bandwidth vector, bandwidth type and kernel type.
}

\usage{
npreg(bws, \dots)

\method{npreg}{formula}(bws, data = NULL, newdata = NULL, \dots)

\method{npreg}{call}(bws, \dots)

\method{npreg}{default}(bws, txdat, tydat, \dots)

\method{npreg}{rbandwidth}(bws,
      txdat = stop("training data 'txdat' missing"),
      tydat = stop("training data 'tydat' missing"),
      exdat,
      eydat,
      gradients = FALSE,
      residuals = FALSE,
      \dots)
}

\arguments{
  \item{bws}{
    a bandwidth specification. This can be set as a \code{rbandwidth}
    object returned from an invocation of \code{\link{npregbw}}, or
    as a vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{txdat}. If specified as
    a vector, then additional arguments will need to be supplied as
    necessary to specify the bandwidth type, kernel types, and so on.
  }
  \item{\dots}{
    additional arguments supplied to specify the regression type,
    bandwidth type, kernel types, training data, and so on, detailed
    below. 
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npregbw}} was called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.  
  }
  
  \item{txdat}{
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the regression estimators. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{tydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{txdat}. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{eydat}{
    a one (1) dimensional numeric or integer vector of the true values
    of the dependent variable. Optional, and used only to calculate the
    true errors. 
  }

  \item{gradients}{
    a logical value indicating that you want gradients computed and
    returned in the resulting \code{npregression} object. Defaults to
    \code{FALSE}.
  }

  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{npregression} object. Defaults to
    \code{FALSE}.
  }
}

\details{

  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    Usage 1: first compute the bandwidth object via npregbw and then
    compute the conditional mean:
    
    bw <- npregbw(y~x)
    ghat <- npreg(bw)
    
    Usage 2: alternatively, compute the bandwidth object indirectly:
    
    ghat <- npreg(y~x)
    
    Usage 3: modify the default kernel and order:
    
    ghat <- npreg(y~x, ckertype="epanechnikov", ckerorder=4)

    Usage 4: use the data frame interface rather than the formula
    interface:

    ghat <- npreg(tydat=y, txdat=x, ckertype="epanechnikov", ckerorder=4)
  }


  
  \code{npreg} implements a variety of methods for regression on
  multivariate (\eqn{p}-variate) data, the types of which are possibly
  continuous and/or discrete (unordered, ordered). The approach is
  based on Li and Racine (2003) who employ \sQuote{generalized product kernels}
  that admit a mix of continuous and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
  density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
  with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
  are constant over the support of \eqn{x}.

  Data contained in the data frame \code{txdat} may be a mix of
  continuous (default), unordered discrete (to be specified in the data
  frame \code{txdat} using \code{\link{factor}}), and ordered discrete
  (to be specified in the data frame \code{txdat} using
  \code{\link{ordered}}). Data can be entered in an arbitrary order and
  data types will be detected automatically by the routine (see
  \code{\link{np}} for details).

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

}

\value{

  \code{npreg} returns a \code{npregression} object.
  The generic
  functions \code{\link{fitted}}, \code{\link{residuals}},
  \code{\link{se}}, \code{\link{predict}}, and
  \code{\link{gradients}}, extract (or generate) estimated values,
  residuals, asymptotic standard
  errors on estimates, predictions, and gradients, respectively, from
  the returned object. Furthermore, the functions \code{\link{summary}}
  and \code{\link{plot}} support objects of this type. The returned object
  has the following components:


  \item{eval}{ evaluation points }
  \item{mean}{ estimates of the regression function (conditional mean) at the
  evaluation points } 
  \item{merr}{ standard errors of the regression function estimates } 
  \item{grad}{ estimates of the gradients at each evaluation point }
  \item{gerr}{ standard errors of the gradient estimates }
  \item{resid}{ if \code{residuals = TRUE}, in-sample or out-of-sample
    residuals where appropriate (or possible) }
  \item{R2}{ coefficient of determination }
  \item{MSE}{ mean squared error }
  \item{MAE}{ mean absolute error }
  \item{MAPE}{ mean absolute percentage error }
  \item{CORR}{ absolute value of Pearson's correlation coefficient }
  \item{SIGN}{ fraction of observations where fitted and observed values
  agree in sign }
  
}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
  discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and Q. Li and J.S. Racine (2007), \dQuote{Nonparametric
  estimation of regression functions in the presence of irrelevant
  regressors,} The Review of Economics and Statistics, 89, 784-789.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
  nonparametric regression,} Statistica Sinica, 14, 485-512.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.
  
  Racine, J.S. and Q. Li (2004), \dQuote{Nonparametric estimation of
    regression functions with both categorical and continuous Data,}
  Journal of Econometrics, 119, 99-130.
  
  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,} Biometrika, 68, 301-309.

}

\author{
  Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.
}

\seealso{
  \code{\link{loess}}
}
\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we compute a
# bivariate nonparametric regression estimate for Giovanni Baiocchi's
# Italian income panel (see Italy for details)

data("Italy")
attach(Italy)

# First, compute the least-squares cross-validated bandwidths for the
# local constant estimator (default). We override the default tolerances
# for the search method as the objective function is well-behaved (don't
# of course do this in general).

bw <- npregbw(formula=gdp~ordered(year), tol=.1, ftol=.1)

# Now take these bandwidths and fit the model and gradients

model <- npreg(bws = bw, gradients = TRUE)

summary(model)

\dontrun{

# Use npplot() to visualize the regression function, add bootstrap
# error bars, and overlay the data on the same plot.

# Note - this may take a minute or two depending on the speed of your
# computer due to bootstrapping being conducted (<ctrl>-C will
# interrupt). Note - nothing will appear in the graphics window until
# all computations are completed (if you use
# plot.errors.method="asymptotic" the figure will instantly appear).

npplot(bws=bw, plot.errors.method="bootstrap")
points(ordered(year), gdp, cex=.2, col="red")

detach(Italy)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we compute a
# bivariate nonparametric regression estimate for Giovanni Baiocchi's
# Italian income panel (see Italy for details)

data("Italy")
attach(Italy)

# First, compute the least-squares cross-validated bandwidths for the
# local constant estimator (default). We override the default tolerances
# for the search method as the objective function is well-behaved (don't
# of course do this in general).

bw <- npregbw(xdat=ordered(year), ydat=gdp, tol=.1, ftol=.1)

# Now take these bandwidths and fit the model and gradients

model <- npreg(bws = bw, gradients = TRUE)

summary(model)

# Use npplot() to visualize the regression function, add bootstrap
# error bars, and overlay the data on the same plot.

# Note - this may take a minute or two depending on the speed of your
# computer due to bootstrapping being conducted (<ctrl>-C will
# interrupt). Note - nothing will appear in the graphics window until
# all computations are completed (if you use
# plot.errors.method="asymptotic" the figure will instantly appear).

npplot(bws=bw, plot.errors.method="bootstrap")
points(ordered(year), gdp, cex=.2, col="red")

detach(Italy)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 2 (INTERFACE=FORMULA): For this example, we compute a local
# linear fit using the AIC_c bandwidth selection criterion. We then plot
# the estimator and its gradient using asymptotic standard errors.

data("cps71")
attach(cps71)

bw <- npregbw(logwage~age, regtype="ll", bwmethod="cv.aic")

# Next, plot the regression function...

npplot(bws=bw, plot.errors.method="asymptotic")
points(age, logwage, cex=.2, col="red")

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Next, plot the derivative...

npplot(bws=bw, plot.errors.method="asymptotic", gradient=TRUE)

detach(cps71)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 2 (INTERFACE=DATA FRAME): For this example, we compute a local
# linear fit using the AIC_c bandwidth selection criterion. We then plot
# the estimator and its gradient using asymptotic standard errors.

data("cps71")
attach(cps71)

bw <- npregbw(xdat=age, ydat=logwage, regtype="ll", bwmethod="cv.aic")

# Next, plot the regression function...

npplot(bws=bw, plot.errors.method="asymptotic")
points(age, logwage, cex=.2, col="red")

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Next, plot the derivative...

npplot(bws=bw, plot.errors.method="asymptotic", gradient=TRUE)

detach(cps71)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 3 (INTERFACE=FORMULA): For this example, we replicate the
# nonparametric regression in Maasoumi, Racine, and Stengos
# (2007) (see oecdpanel for details). Note that X is multivariate
# containing a mix of unordered, ordered, and continuous data types. Note
# - this may take a few minutes depending on the speed of your computer.

data("oecdpanel")
attach(oecdpanel)

# We override the default tolerances for the search method as the
# objective function is well-behaved (don't of course do this in
# general).

bw <- npregbw(formula=growth~
              factor(oecd)+
              factor(year)+
              initgdp+
              popgro+
              inv+
              humancap,
              tol=.1, ftol=.1)

npplot(bws=bw, plot.errors.method="asymptotic")

detach(oecdpanel)

# EXAMPLE 3 (INTERFACE=DATA FRAME): For this example, we replicate the
# nonparametric regression in Maasoumi, Racine, and Stengos
# (2007) (see oecdpanel for details). Note that X is multivariate
# containing a mix of unordered, ordered, and continuous data types. Note
# - this may take a few minutes depending on the speed of your computer.

data("oecdpanel")
attach(oecdpanel)

y <- growth
X <- data.frame(factor(oecd), factor(year), initgdp, popgro, inv, humancap)

# We override the default tolerances for the search method as the
# objective function is well-behaved (don't of course do this in
# general).

bw <- npregbw(xdat=X, ydat=y, tol=.1, ftol=.1)

npplot(bws=bw, plot.errors.method="asymptotic")

detach(oecdpanel)


# EXAMPLE 4 (INTERFACE=FORMULA): Experimental data - the effect of
# vitamin C on tooth growth in guinea pigs
#
# Description:
#
#     The response is the length of odontoblasts (teeth) in each of 10
#     guinea pigs at each of three dose levels of Vitamin C (0.5, 1, and
#     2 mg) with each of two delivery methods (orange juice or ascorbic
#     acid).
#
# Usage:
#
#     ToothGrowth
#
# Format:
#
#     A data frame with 60 observations on 3 variables.
#
#       [,1]  len   numeric  Tooth length
#       [,2]  supp  factor   Supplement type (VC or OJ).
#       [,3]  dose  numeric  Dose in milligrams.

library("datasets")
attach(ToothGrowth)

# Note - in this example, there are six cells. 

bw <- npregbw(formula=len~factor(supp)+ordered(dose))

# Now plot the partial regression surfaces with bootstrapped
# nonparametric confidence bounds

npplot(bws=bw, plot.errors.method="bootstrap", plot.errors.type="quantile")

detach(ToothGrowth)

# EXAMPLE 4 (INTERFACE=DATA FRAME): Experimental data - the effect of
# vitamin C on tooth growth in guinea pigs
#
# Description:
#
#     The response is the length of odontoblasts (teeth) in each of 10
#     guinea pigs at each of three dose levels of Vitamin C (0.5, 1, and
#     2 mg) with each of two delivery methods (orange juice or ascorbic
#     acid).
#
# Usage:
#
#     ToothGrowth
#
# Format:
#
#     A data frame with 60 observations on 3 variables.
#
#       [,1]  len   numeric  Tooth length
#       [,2]  supp  factor   Supplement type (VC or OJ).
#       [,3]  dose  numeric  Dose in milligrams.

library("datasets")
attach(ToothGrowth)

# Note - in this example, there are six cells. 

y <- len
X <- data.frame(supp=factor(supp), dose=ordered(dose))

bw <- npregbw(X, y)

# Now plot the partial regression surfaces with bootstrapped
# nonparametric confidence bounds

npplot(bws=bw, plot.errors.method="bootstrap", plot.errors.type="quantile")

detach(ToothGrowth)

# EXAMPLE 5 (INTERFACE=FORMULA): a pretty 2-d smoothing example adapted
# from the R mgcv library which was written by Simon N. Wood.

set.seed(12345)

# This function generates a smooth nonlinear DGP

dgp.func <- function(x, z, sx=0.3, sz=0.4)
  { (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
                 0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
  }

# Generate 500 observations, compute the true DGP (i.e., no noise), 
# then a noisy sample

n <- 500

x <- runif(n)
z <- runif(n)

xs <- seq(0, 1, length=30)
zs <- seq(0, 1, length=30)

X.eval <- data.frame(x=rep(xs, 30), z=rep(zs, rep(30, 30)))

dgp <- matrix(dgp.func(X.eval$x, X.eval$z), 30, 30)

y <- dgp.func(x, z)+rnorm(n)*0.1

# Prepare the screen for output... first, plot the true DGP

split.screen(c(2, 1))

screen(1)

persp(xs, zs, dgp, xlab="x1", ylab="x2", zlab="y", main="True DGP")

# Next, compute a local linear fit and plot that

bw <- npregbw(formula=y~x+z, regtype="ll", bwmethod="cv.aic",
                       tol=.01, ftol=.01)
fit <- fitted(npreg(bws=bw, newdata=X.eval))
fit.mat <- matrix(fit, 30, 30)

screen(2)

persp(xs, zs, fit.mat, xlab="x1", ylab="x2", zlab="y",
      main="Local linear estimate")

# EXAMPLE 5 (INTERFACE=DATA FRAME): a pretty 2-d smoothing example
# adapted from the R mgcv library which was written by Simon N. Wood.

set.seed(12345)

# This function generates a smooth nonlinear DGP

dgp.func <- function(x, z, sx=0.3, sz=0.4)
  { (pi**sx*sz)*(1.2*exp(-(x-0.2)^2/sx^2-(z-0.3)^2/sz^2)+
                 0.8*exp(-(x-0.7)^2/sx^2-(z-0.8)^2/sz^2))
  }

# Generate 500 observations, compute the true DGP (i.e., no noise), 
# then a noisy sample

n <- 500

x <- runif(n)
z <- runif(n)

xs <- seq(0, 1, length=30)
zs <- seq(0, 1, length=30)

X <- data.frame(x, z)
X.eval <- data.frame(x=rep(xs, 30), z=rep(zs, rep(30, 30)))

dgp <- matrix(dgp.func(X.eval$x, X.eval$z), 30, 30)

y <- dgp.func(x, z)+rnorm(n)*0.1

# Prepare the screen for output... first, plot the true DGP

split.screen(c(2, 1))

screen(1)

persp(xs, zs, dgp, xlab="x1", ylab="x2", zlab="y", main="True DGP")

# Next, compute a local linear fit and plot that

bw <- npregbw(xdat=X, ydat=y, regtype="ll", bwmethod="cv.aic",
                       tol=.01, ftol=.01)
fit <- fitted(npreg(exdat=X.eval, bws=bw))
fit.mat <- matrix(fit, 30, 30)

screen(2)

persp(xs, zs, fit.mat, xlab="x1", ylab="x2", zlab="y",
      main="Local linear estimate")
} % enddontrun
}
\keyword{nonparametric}

