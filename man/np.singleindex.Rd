% $ID$
\name{npindex}
\alias{npindex}
\alias{npindex.call}
\alias{npindex.default}
\alias{npindex.formula}
\alias{npindex.sibandwidth}

\title{Semiparametric Single Index Model}

\description{
  \code{npindex} computes a semiparametric single index model
  for a dependent variable and \eqn{p}-variate explanatory data using
  the model \eqn{Y = G(X\beta) + \epsilon}{Y = G(XB) + epsilon}, given a
  set of evaluation points, training points (consisting of explanatory
  data and dependent data), and a \code{npindexbw} bandwidth
  specification. Note that for this semiparametric estimator, the
  bandwidth object contains parameters for the single index model and
  the (scalar) bandwidth for the index function.
}

\usage{
npindex(bws, \dots)

\method{npindex}{formula}(bws, data = NULL, newdata = NULL, \dots)

\method{npindex}{call}(bws, \dots)

\method{npindex}{default}(bws, txdat, tydat, \dots) 

\method{npindex}{sibandwidth}(bws,
        txdat = stop("training data 'txdat' missing"),
        tydat = stop("training data 'tydat' missing"),
        exdat,
        eydat,
        gradients = FALSE,
        residuals = FALSE,
        errors = FALSE,
        boot.num = 399,
        \dots)
}

\arguments{

  \item{bws}{
    a bandwidth specification. This can be set as a
    \code{sibandwidth}
    object returned from an invocation of \code{npindexbw}, or
    as a vector of parameters (beta) with each element \eqn{i}
    corresponding to the coefficient for column \eqn{i} in \code{txdat}
    where the first element is normalized to 1, and a scalar bandwidth
    (h). 
  }

  \item{gradients}{
    a logical value indicating that you want gradients computed and
    returned in the resulting \code{singleindex} object. Defaults to
    \code{FALSE}.
  }

  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{singleindex} object. Defaults to
    \code{FALSE}.
  }

  \item{errors}{
    a logical value indicating that you want (bootstrapped)
    standard errors for the conditional mean, gradients (when
    \code{gradients=TRUE} is set), and average gradients (when
    \code{gradients=TRUE} is set), computed and returned in the
    resulting \code{singleindex} object. Defaults to \code{FALSE}.
  }

  \item{boot.num}{
    an integer specifying the number of bootstrap replications to use
    when performing standard error calculations. Defaults to
    \code{399}. 
  }

  \item{\dots}{
    additional arguments supplied to specify the parameters to the
    \code{sibandwidth} S3 method, which is called during estimation. 
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npindexbw}} was called.
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
  
}

\details{

  A matrix of gradients along with average derivatives are computed and
  returned if \code{gradients=TRUE} is used.

}

\value{

  \code{npindex} returns a \code{npsingleindex} object.   The generic
  functions \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{coef}},
  \code{\link{se}}, \code{\link{predict}}, and
  \code{\link{gradients}}, extract (or generate) estimated values,
  residuals, coefficients, bootstrapped standard
  errors on estimates, predictions, and gradients, respectively, from
  the returned object. Furthermore, the functions \code{\link{summary}}
  and \code{\link{plot}} support objects of this type. The returned object
  has the following components:

  \item{eval}{ evaluation points }
  \item{mean}{ estimates of the regression function (conditional mean) at the
  evaluation points } 
  \item{beta}{ the model coefficients }
  \item{merr}{ standard errors of the regression function estimates } 
  \item{grad}{ estimates of the gradients at each evaluation point }
  \item{gerr}{ standard errors of the gradient estimates }
  \item{mean.grad}{ mean (average) gradient over the evaluation points}
  \item{mean.gerr}{ bootstrapped standard error of the mean gradient estimates }
  \item{R2}{ if \code{method="ichimura"}, coefficient of determination }
  \item{MSE}{if \code{method="ichimura"}, mean squared error }
  \item{MAE}{if \code{method="ichimura"}, mean absolute error }
  \item{MAPE}{if \code{method="ichimura"}, mean absolute percentage error }
  \item{CORR}{if \code{method="ichimura"}, absolute value of Pearson's correlation coefficient }
  \item{SIGN}{if \code{method="ichimura"}, fraction of observations where fitted and observed values
    agree in sign }
  \item{confusion.matrix}{if \code{method="kleinspady"}, the confusion matrix or \code{NA} if outcomes
    are not available }
  \item{CCR.overall}{if \code{method="kleinspady"},  the overall correct
    classification ratio, or \code{NA} if outcomes are not available }
  \item{CCR.byoutcome}{if \code{method="kleinspady"},  a numeric vector containing the correct
    classification ratio by outcome, or \code{NA} if outcomes are not
    available }
  \item{fit.mcfadden}{if \code{method="kleinspady"},  the McFadden-Puig-Kerschner performance measure
    or \code{NA} if outcomes are not available }
}



\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Ichimura, H., (1993), \dQuote{Semiparametric least squares (SLS) and
  weighted SLS estimation of single-index models,} Journal of
  Econometrics, 58, 71-120.

  Klein, R. W. and R. H. Spady (1993), \dQuote{An efficient semiparametric
  estimator for binary response models,} Econometrica, 61, 387-421.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice,} Princeton University Press.

  McFadden, D. and C. Puig and D. Kerschner (1977), \dQuote{Determinants
    of the long-run demand for electricity,} Proceedings of the
  American Statistical Association (Business and Economics Section),
  109-117. 

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

\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): Generate a simple linear model then
# estimate it using a semiparametric single index specification and
# Ichimura's nonlinear least squares coefficients and bandwidth
# (default). Also compute the matrix of gradients and average derivative
# estimates.

set.seed(12345)

n <- 100

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

y <- x1 - x2 + rnorm(n)

# Note - this may take a minute or two depending on the speed of your
# computer. Note also that the first element of the vector beta is
# normalized to one for identification purposes, and that X must contain
# at least one continuous variable.

bw <- npindexbw(formula=y~x1+x2)

summary(bw)

model <- npindex(bws=bw, gradients=TRUE)

summary(model)

\dontrun{

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Or you can visualize the input with npplot.

npplot(bws=bw)

Sys.sleep(5)

# EXAMPLE 1 (INTERFACE=DATA FRAME): Generate a simple linear model then
# estimate it using a semiparametric single index specification and
# Ichimura's nonlinear least squares coefficients and bandwidth
# (default). Also compute the matrix of gradients and average derivative
# estimates.

set.seed(12345)

n <- 100

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

y <- x1 - x2 + rnorm(n)

X <- cbind(x1, x2)

# Note - this may take a minute or two depending on the speed of your
# computer. Note also that the first element of the vector beta is
# normalized to one for identification purposes, and that X must contain
# at least one continuous variable.

bw <- npindexbw(xdat=X, ydat=y)

summary(bw)

model <- npindex(bws=bw, gradients=TRUE)

summary(model)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Or you can visualize the input with npplot.

npplot(bws=bw)

Sys.sleep(5)

# EXAMPLE 2 (INTERFACE=FORMULA): Generate a simple binary outcome linear
# model then estimate it using a semiparametric single index
# specification and Klein and Spady's likelihood-based coefficients and
# bandwidth (default). Also compute the matrix of gradients and average
# derivative estimates.

n <- 100

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)

# Note that the first element of the vector beta is normalized to one
# for identification purposes, and that X must contain at least one
# continuous variable.

bw <- npindexbw(formula=y~x1+x2, method="kleinspady")

summary(bw)

model <- npindex(bws=bw, gradients=TRUE)

# Note that, since the outcome is binary, we can assess model
# performance using methods appropriate for binary outcomes. We look at
# the confusion matrix, various classification ratios, and McFadden et
# al's measure of predictive performance.

summary(model)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 2 (INTERFACE=DATA FRAME): Generate a simple binary outcome
# linear model then estimate it using a semiparametric single index
# specification and Klein and Spady's likelihood-based coefficients and
# bandwidth (default). Also compute the matrix of gradients and average
# derivative estimates.

n <- 100

x1 <- runif(n, min=-1, max=1)
x2 <- runif(n, min=-1, max=1)

y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)

X <- cbind(x1, x2)

# Note that the first element of the vector beta is normalized to one
# for identification purposes, and that X must contain at least one
# continuous variable.

bw <- npindexbw(xdat=X, ydat=y, method="kleinspady")

summary(bw)

model <- npindex(bws=bw, gradients=TRUE)

# Note that, since the outcome is binary, we can assess model
# performance using methods appropriate for binary outcomes. We look at
# the confusion matrix, various classification ratios, and McFadden et
# al's measure of predictive performance.

summary(model)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# EXAMPLE 3 (INTERFACE=FORMULA): Replicate the DGP of Klein & Spady
# (1993) (see their description on page 405, pay careful attention to
# footnote 6 on page 405).

set.seed(123)

n <- 1000

# x1 is chi-squared having 3 df truncated at 6 standardized by
# subtracting 2.348 and dividing by 1.511

x <- rchisq(n, df=3)
x1 <- (ifelse(x < 6, x, 6) - 2.348)/1.511

# x2 is normal (0, 1) truncated at +- 2 divided by 0.8796

x <- rnorm(n)
x2 <- ifelse(abs(x) < 2 , x, 2) / 0.8796

# y is 1 if y* > 0, 0 otherwise.

y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)

# Compute the parameter vector and bandwidth. Note that the first
# element of the vector beta is normalized to one for identification
# purposes, and that X must contain at least one continuous variable.


bw <- npindexbw(formula=y~x1+x2, method="kleinspady")

# Next, create the evaluation data in order to generate a perspective
# plot

# Create an evaluation data matrix

x1.seq <- seq(min(x1), max(x1), length=50)
x2.seq <- seq(min(x2), max(x2), length=50)
X.eval <- expand.grid(x1=x1.seq, x2=x2.seq)

# Now evaluate the single index model on the evaluation data

fit <- fitted(npindex(exdat=X.eval, 
               eydat=rep(1, nrow(X.eval)), 
               bws=bw))

# Finally, coerce the fitted model into a matrix suitable for 3D
# plotting via persp()

fit.mat <- matrix(fit, 50, 50)

# Generate a perspective plot similar to Figure 2 b of Klein and Spady
# (1993)

persp(x1.seq, 
      x2.seq, 
      fit.mat, 
      col="white", 
      ticktype="detailed", 
      expand=0.5, 
      axes=FALSE, 
      box=FALSE, 
      main="Estimated Semiparametric Probability Perspective", 
      theta=310, 
      phi=25)

# EXAMPLE 3 (INTERFACE=DATA FRAME): Replicate the DGP of Klein & Spady
# (1993) (see their description on page 405, pay careful attention to
# footnote 6 on page 405).

set.seed(123)

n <- 1000

# x1 is chi-squared having 3 df truncated at 6 standardized by
# subtracting 2.348 and dividing by 1.511

x <- rchisq(n, df=3)
x1 <- (ifelse(x < 6, x, 6) - 2.348)/1.511

# x2 is normal (0, 1) truncated at +- 2 divided by 0.8796

x <- rnorm(n)
x2 <- ifelse(abs(x) < 2 , x, 2) / 0.8796

# y is 1 if y* > 0, 0 otherwise.

y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)

# Create the X matrix

X <- cbind(x1, x2)

# Compute the parameter vector and bandwidth. Note that the first
# element of the vector beta is normalized to one for identification
# purposes, and that X must contain at least one continuous variable.


bw <- npindexbw(xdat=X, ydat=y, method="kleinspady")

# Next, create the evaluation data in order to generate a perspective
# plot

# Create an evaluation data matrix

x1.seq <- seq(min(x1), max(x1), length=50)
x2.seq <- seq(min(x2), max(x2), length=50)
X.eval <- expand.grid(x1=x1.seq, x2=x2.seq)

# Now evaluate the single index model on the evaluation data

fit <- fitted(npindex(exdat=X.eval, 
               eydat=rep(1, nrow(X.eval)), 
               bws=bw))

# Finally, coerce the fitted model into a matrix suitable for 3D
# plotting via persp()

fit.mat <- matrix(fit, 50, 50)

# Generate a perspective plot similar to Figure 2 b of Klein and Spady
# (1993)

persp(x1.seq, 
      x2.seq, 
      fit.mat, 
      col="white", 
      ticktype="detailed", 
      expand=0.5, 
      axes=FALSE, 
      box=FALSE, 
      main="Estimated Semiparametric Probability Perspective", 
      theta=310, 
      phi=25)
} % enddontrun
}
\keyword{nonparametric}
