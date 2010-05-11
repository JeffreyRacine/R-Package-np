% $Id: np.condensity.Rd,v 1.67 2006/11/03 21:19:03 tristen Exp $
\name{npcdens}
\alias{npcdens}
\alias{npcdens.call}
\alias{npcdens.conbandwidth}
\alias{npcdens.default}
\alias{npcdens.formula}
\alias{npcdist}
\title{Kernel Conditional Density and Distribution Estimates with Mixed Data Types}

\description{
  \code{npcdens} computes kernel conditional density estimates on
  \eqn{p+q}-variate evaluation data, given a set of training data (both
  explanatory and dependent) and a bandwidth specification (a
  \code{conbandwidth} object or a bandwidth vector, bandwidth type, and
  kernel type) using the method of Hall, Racine, and Li (2004).
  Similarly \code{npcdist} computes kernel
  conditional cumulative distribution estimates.
  The data may be continuous, discrete (unordered and ordered
  factors), or some combination thereof. 
}

\usage{
npcdens(bws, \dots)

\method{npcdens}{formula}(bws, data = NULL, newdata = NULL, \dots)

\method{npcdens}{call}(bws, \dots)

\method{npcdens}{conbandwidth}(bws,
        txdat = stop("invoked without training data 'txdat'"),
        tydat = stop("invoked without training data 'tydat'"),
        exdat,
        eydat,
        gradients = FALSE,
        \dots)

\method{npcdens}{default}(bws, txdat, tydat, \dots)

}

\arguments{
  \item{bws}{
    a bandwidth specification. This can be set as a \code{conbandwidth}
    object returned from a previous invocation of \code{\link{npcdensbw}}, or
    as a \eqn{p+q}-vector 
    of bandwidths, with each element \eqn{i} up to \eqn{i=p}
    corresponding to the bandwidth for column \eqn{i} in \code{txdat}, 
    and each element \eqn{i} from \eqn{i=p+1} to \eqn{i=p+q}
    corresponding to the bandwidth for column \eqn{i-p} in
    \code{tydat}. If specified as a vector, then
    additional arguments will need to be supplied as necessary to
    specify the bandwidth type, kernel types, training data, and so
    on. 
  }

  \item{gradients}{
    a logical value specifying whether to return estimates of the
    gradients at the evaluation points. Defaults to \code{FALSE}.
  }

  \item{\dots}{
    additional arguments supplied to specify the bandwidth type,
    kernel types, and so on.  This is necessary if you
    specify bws as a \eqn{p+q}-vector and not a \code{conbandwidth}
    object, and you do not desire the default behaviours. To do this,
    you may specify any of \code{bwmethod}, \code{bwscaling},
    \code{bwtype}, \code{cxkertype}, \code{cxkerorder},
    \code{cykertype}, \code{cykerorder}, \code{uxkertype},
    \code{uykertype}, \code{oxkertype}, \code{oykertype}, as described
    in \code{\link{npcdensbw}}. 

  }
  
  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
   \code{\link{npcdensbw}} was called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.  
  }


  \item{txdat}{
    a \eqn{p}-variate data frame of sample realizations of explanatory
    data (training data). Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{tydat}{
    a \eqn{q}-variate data frame of sample realizations of dependent
    data (training data). Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of explanatory data on which
    conditional densities will be evaluated. By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{eydat}{
    a \eqn{q}-variate data frame of dependent data on which conditional
    densities will be evaluated. By default,
    evaluation takes place on the data provided by \code{tydat}.
  }

}
\details{
  \code{npcdens} and \code{\link{npcdist}} implement a
  variety of methods for estimating multivariate conditional
  distributions (\eqn{p+q}-variate) defined over a set of possibly
  continuous and/or discrete (unordered, ordered) data. The approach is
  based on Li and Racine (2004) who employ \sQuote{generalized product
    kernels} that admit a mix of continuous and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating
  the density at the point \eqn{x}. Generalized nearest-neighbor
  bandwidths change with the point at which the density is estimated,
  \eqn{x}. Fixed bandwidths are constant over the support of \eqn{x}.

  Training and evaluation input data  may be a
  mix of continuous (default), unordered discrete (to be specified in
  the data frames using \code{\link{factor}}), and ordered discrete (to be
  specified in the data frames using \code{\link{ordered}}). Data can be
  entered in an arbitrary order and data types will be detected
  automatically by the routine (see \code{\link{np}} for details).

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

}
\value{
  \code{npcdens} returns a \code{condensity} object, similarly
  \code{\link{npcdist}} returns a \code{condistribution} object. The generic
  accessor functions \code{\link{fitted}}, \code{\link{se}}, and
  \code{\link{gradients}}, extract estimated values, asymptotic standard
  errors on estimates, and gradients, respectively, from
  the returned object. Furthermore, the functions \code{\link{predict}},
  \code{\link{summary}}
  and \code{\link{plot}} support objects of both classes. The returned objects
  have the following components:

  \item{xbw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    explanatory data, \code{txdat} }
  \item{ybw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    dependent data, \code{tydat} }
  \item{xeval}{ the evaluation points of the explanatory data }
  \item{yeval}{ the evaluation points of the dependent data }
  \item{condens or condist}{ estimates of the conditional density (cumulative
    distribution) at the evaluation points }
  \item{conderr}{ standard errors of the conditional density (cumulative distribution)
    estimates }
  \item{congrad}{ if invoked with \code{gradients = TRUE}, estimates of
    the gradients at the evaluation points }
  \item{congerr}{ if invoked with \code{gradients = TRUE}, standard
    errors of the gradients at the evaluation points }
  \item{log_likelihood}{ log likelihood of the conditional density estimate }

}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and J.S. Racine and Q. Li (2004), \dQuote{Cross-validation and the
    estimation of conditional probability densities,} Journal of the
  American Statistical Association, 99, 1015-1026.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,} Cambridge
  University Press. 

  Scott, D.W. (1992), \emph{Multivariate Density Estimation. Theory,
  Practice and Visualization,} New York: Wiley.

  Silverman, B.W. (1986), \emph{Density Estimation,} London: Chapman and
  Hall.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
    for discrete distributions,}  Biometrika, 68, 301-309.

}
\author{
  Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
  \email{racinej@mcmaster.ca}    
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.
}

\seealso{
  \code{\link{npudens}}
}
\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we load Giovanni
# Baiocchi's Italian GDP panel (see Italy for details), and compute the
# likelihood cross-validated bandwidths (default) using a second-order
# Gaussian kernel (default). Note - this may take a minute or two
# depending on the speed of your computer.

data("Italy")
attach(Italy)

# First, compute the bandwidths... note that this may take a minute or
# two depending on the speed of your computer. We override the default
# tolerances for the search method as the objective function is
# well-behaved (don't of course do this in general).

bw <- npcdensbw(formula=gdp~ordered(year), tol=.1, ftol=.1)

# Next, compute the condensity object...

fhat <- npcdens(bws=bw)

# The object fhat now contains results such as the estimated conditional
# density function (fhat$condens) and so on...

summary(fhat)

\dontrun{

# Call the npplot() function to visualize the results (<ctrl>-C will
# interrupt on *NIX systems, <esc> will interrupt on MS Windows
# systems).

npplot(bws=bw)

# To plot the conditional distribution we use cdf=TRUE in npplot
# (<ctrl>-C will interrupt on *NIX systems, <esc> will interrupt on MS
# Windows systems)

npplot(bws=bw, cdf=TRUE)

detach(Italy)

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we load Giovanni
# Baiocchi's Italian GDP panel (see Italy for details), and compute the
# likelihood cross-validated bandwidths (default) using a second-order
# Gaussian kernel (default). Note - this may take a minute or two
# depending on the speed of your computer.

data("Italy")
attach(Italy)

# First, compute the bandwidths... note that this may take a minute or
# two depending on the speed of your computer. We override the default
# tolerances for the search method as the objective function is
# well-behaved (don't of course do this in general).

# Note - we cast `X' and `y' as data frames so that npplot() can
# automatically grab names (this looks like overkill, but in
# multivariate settings you would do this anyway, so may as well get in
# the habit).

X <- data.frame(year=ordered(year))
y <- data.frame(gdp)

bw <- npcdensbw(xdat=X, ydat=y, tol=.1, ftol=.1)

# Next, compute the condensity object...

fhat <- npcdens(bws=bw)

# The object fhat now contains results such as the estimated conditional
# density function (fhat$condens) and so on...

summary(fhat)

# Call the npplot() function to visualize the results (<ctrl>-C will
# interrupt on *NIX systems, <esc> will interrupt on MS Windows systems).

npplot(bws=bw)

# To plot the conditional distribution we use cdf=TRUE in npplot
# (<ctrl>-C will interrupt on *NIX systems, <esc> will interrupt on MS
# Windows systems)

npplot(bws=bw, cdf=TRUE)

detach(Italy)

# EXAMPLE 2 (INTERFACE=FORMULA): For this example, we load the old
# faithful geyser data from the R `datasets' library and compute the
# conditional density and conditional distribution functions.

library("datasets")
data("faithful")
attach(faithful)

# Note - this may take a few minutes depending on the speed of your
# computer...

bw <- npcdensbw(formula=eruptions~waiting)

summary(bw)

# Plot the density function (<ctrl>-C will interrupt on *NIX systems,
# <esc> will interrupt on MS Windows systems).

npplot(bws=bw)

# Plot the distribution function (cdf=TRUE) (<ctrl>-C will interrupt on
# *NIX systems, <esc> will interrupt on MS Windows systems)

npplot(bws=bw, cdf=TRUE)

detach(faithful)

# EXAMPLE 2 (INTERFACE=DATA FRAME): For this example, we load the old
# faithful geyser data from the R `datasets' library and compute the
# conditional density and conditional distribution functions.

library("datasets")
data("faithful")
attach(faithful)

# Note - this may take a few minutes depending on the speed of your
# computer...

# Note - we cast `X' and `y' as data frames so that npplot() can
# automatically grab names (this looks like overkill, but in
# multivariate settings you would do this anyway, so may as well get in
# the habit).

X <- data.frame(waiting)
y <- data.frame(eruptions)

bw <- npcdensbw(xdat=X, ydat=y)

summary(bw)

# Plot the density function (<ctrl>-C will interrupt on *NIX systems, 
# <esc> will interrupt on MS Windows systems)

npplot(bws=bw)

# Plot the distribution function (cdf=TRUE) (<ctrl>-C will interrupt on
# *NIX systems, <esc> will interrupt on MS Windows systems)

npplot(bws=bw, cdf=TRUE)

detach(faithful)

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

# Generate data-driven bandwidths (likelihood cross-validation).  We
# override the default tolerances for the search method as the objective
# function is well-behaved (don't of course do this in general).  Note -
# this may take a few minutes depending on the speed of your computer...

bw <- npcdensbw(formula=factor(y)~x1+x2, tol=.1, ftol=.1)

# Next, create the evaluation data in order to generate a perspective
# plot

x1.seq <- seq(min(x1), max(x1), length=50)
x2.seq <- seq(min(x2), max(x2), length=50)
X.eval <- expand.grid(x1=x1.seq,x2=x2.seq)

data.eval <- data.frame(y=factor(rep(1, nrow(X.eval))),x1=X.eval[,1],x2=X.eval[,2])

# Now evaluate the conditional probability for y=1 and for the
# evaluation Xs

fit <- fitted(npcdens(bws=bw,newdata=data.eval))

# Finally, coerce the data into a matrix for plotting with persp()

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
      main="Estimated Nonparametric Probability Perspective", 
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

# Generate data-driven bandwidths (likelihood cross-validation).  We
# override the default tolerances for the search method as the objective
# function is well-behaved (don't of course do this in general).  Note -
# this may take a few minutes depending on the speed of your computer...

bw <- npcdensbw(xdat=X, ydat=factor(y), tol=.1, ftol=.1)

# Next, create the evaluation data in order to generate a perspective
# plot

x1.seq <- seq(min(x1), max(x1), length=50)
x2.seq <- seq(min(x2), max(x2), length=50)
X.eval <- expand.grid(x1=x1.seq,x2=x2.seq)

# Now evaluate the conditional probability for y=1 and for the
# evaluation Xs

fit <- fitted(npcdens(exdat=X.eval, 
               eydat=factor(rep(1, nrow(X.eval))), 
               bws=bw))

# Finally, coerce the data into a matrix for plotting with persp()

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
      main="Estimated Nonparametric Probability Perspective", 
      theta=310, 
      phi=25)
} % enddontrun
}
\keyword{nonparametric}

