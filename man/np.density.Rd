% $Id: np.density.Rd,v 1.75 2006/11/03 21:19:03 tristen Exp $
\name{npudens}
\alias{npudens}
\alias{npudens.bandwidth}
\alias{npudens.call}
\alias{npudens.default}
\alias{npudens.formula}

\alias{npudist}
\title{ Kernel Density and Distribution Estimation with Mixed Data Types }
\description{

  \code{npudens} computes  kernel unconditional density estimates on
  evaluation data, given a set of training data and a bandwidth
  specification (a \code{bandwidth} object or a bandwidth vector,
  bandwidth type, and kernel type) using the method of Li and Racine
  (2003).  Similarly \code{\link{npudist}} computes  kernel
  unconditional cumulative distribution estimates.

}
\usage{
npudens(bws, \dots)

\method{npudens}{formula}(bws, data = NULL, newdata = NULL, \dots)

\method{npudens}{bandwidth}(bws,
        tdat = stop("invoked without training data 'tdat'"),
        edat,
        \dots)

\method{npudens}{call}(bws, \dots)

\method{npudens}{default}(bws, tdat, \dots)
}

\arguments{
  \item{bws}{
    a bandwidth specification. This can be set as a \code{bandwidth}
    object returned from an invocation of \code{\link{npudensbw}}, or as a
    \eqn{p}-vector of bandwidths, with an element for each variable in the
    training data. If specified as a vector, then additional arguments
    will need to be supplied as necessary to change them from the
    defaults to specify the bandwidth type, kernel types, training data,
    and so on.
  }

  \item{\dots}{
    additional arguments supplied to specify, the training
    data, the
    bandwidth type, kernel types, and so on.
    This is necessary if you specify bws as a \eqn{p}-vector and not
    a \code{bandwidth} object, and you do not desire the default
    behaviours. To do this, you may specify any of \code{bwscaling},
    \code{bwtype}, \code{ckertype}, \code{ckerorder}, \code{ukertype},
    \code{okertype}, as described in \code{\link{npudensbw}}.
  }

  \item{tdat}{
    a \eqn{p}-variate data frame of sample realizations (training data)
    used to estimate the density. Defaults to the training data used to
    compute the bandwidth object.
  }
  \item{edat}{
    a \eqn{p}-variate data frame of density evaluation points. By default,
    evaluation takes place on the data provided by \code{tdat}.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npudensbw}} was called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.  
  }
  
}

\details{

  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    Usage 1: first compute the bandwidth object via npudensbw and then
    compute the density:
    
    bw <- npudensbw(~y)
    fhat <- npudens(bw)
    
    Usage 2: alternatively, compute the bandwidth object indirectly:
    
    fhat <- npudens(~y)
    
    Usage 3: modify the default kernel and order:
    
    fhat <- npudens(~y, ckertype="epanechnikov", ckerorder=4)

    Usage 4: use the data frame interface rather than the formula
    interface:

    fhat <- npudens(tdat = y, ckertype="epanechnikov", ckerorder=4)
  }

  \code{npudens} and \code{\link{npudist}} implement a variety
  of methods for estimating multivariate distributions (\eqn{p}-variate)
  defined over a set of possibly continuous and/or discrete (unordered,
  ordered) data. The approach is based on Li and Racine (2003) who
  employ \sQuote{generalized product kernels} that admit a mix of
  continuous and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating
  the density at the point \eqn{x}. Generalized nearest-neighbor
  bandwidths change with the point at which the density is estimated,
  \eqn{x}. Fixed bandwidths are constant over the support of \eqn{x}.

  Data contained in the data frame \code{tdat} (and also \code{edat})
  may be a mix of continuous (default), unordered discrete (to be
  specified in the data frame \code{tdat} using the \code{\link{factor}}
  command), and ordered discrete (to be specified in the data frame
  \code{tdat} using the \code{\link{ordered}} command). Data can be
  entered in an arbitrary order and data types will be detected
  automatically by the routine (see \code{\link{np}} for details).

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and
  eighth order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of
  the Wang and van Ryzin (1981) kernel.

}
\value{
  
  \code{npudens} returns a \code{npdensity} object, similarly
  \code{\link{npudist}} returns a \code{npdistribution} object. The
  generic accessor functions \code{\link{fitted}}, and \code{\link{se}},
  extract estimated values and asymptotic standard errors on estimates,
  respectively, from the returned object. Furthermore, the functions
  \code{\link{predict}}, \code{\link{summary}} and \code{\link{plot}}
  support objects of both classes. The returned objects have the
  following components:

  \item{eval}{ the evaluation points. }
  \item{dens or dist}{ estimation of the density (cumulative
   distribution) at the evaluation points }
  \item{derr}{ standard errors of the density (cumulative distribution)
   estimates }
  \item{log_likelihood}{ log likelihood of the density estimates }
  
}
\references{
  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{ Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2003), \dQuote{Nonparametric estimation of
    distributions with categorical and continuous data,} Journal
  of Multivariate Analysis, 86, 266-292.

  Ouyang, D. and Q. Li and J.S. Racine (2006), \dQuote{Cross-validation
    and the estimation of probability distributions with categorical
    data,} Journal of Nonparametric Statistics, 18, 69-100.

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

\seealso{ \code{\link{npudensbw}} , \code{\link{density}} }

\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we load Giovanni
# Baiocchi's Italian GDP panel (see Italy for details), then create a
# data frame in which year is an ordered factor, GDP is continuous,
# compute bandwidths using likelihood cross-validation, then create a
# grid of data on which the density will be evaluated for plotting
# purposes.

data("Italy")
attach(Italy)

# Compute bandwidths using likelihood cross-validation (default), but we
# override the default search tolerances to speed things up as this is a
# well-behaved dataset (don't of course do this in general).

bw <- npudensbw(formula=~ordered(year)+gdp, tol=.01, ftol=.01)

# At this stage you could use npudens() to do a variety of
# things. Here we compute the npudens() object and place it in fhat.

fhat <- npudens(bws=bw)

# Note that simply typing the name of the object returns some useful
# information. For more info, one can call summary: 

summary(fhat)

\dontrun{

# Next, we illustrate how to create a grid of `evaluation data' and feed
# it to the perspective plotting routines in R, among others.

# Create an evaluation data matrix

year.seq <- sort(unique(year))
gdp.seq <- seq(1,36,length=50)
data.eval <- expand.grid(year=year.seq,gdp=gdp.seq)

# Generate the estimated density computed for the evaluation data

fhat <- fitted(npudens(bws=bw, newdata=data.eval))

# Coerce the data into a matrix for plotting with persp()

f <- matrix(fhat, length(unique(year)), 50)

# Next, create a 3D perspective plot of the PDF f, and a 2D
# contour plot.

persp(as.integer(levels(year.seq)), gdp.seq, f, col="lightblue",
      ticktype="detailed", ylab="GDP", xlab="Year", zlab="Density",
      theta=300, phi=50)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

contour(as.integer(levels(year.seq)), 
        gdp.seq, 
        f, 
        xlab="Year", 
        ylab="GDP", 
        main = "Density Contour Plot", 
        col=topo.colors(100))

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Alternatively, you could use the npplot() command (<ctrl>-C will
# interrupt on *NIX systems, <esc> will interrupt on MS Windows
# systems).

npplot(bws=bw)

detach(Italy)

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we load Giovanni
# Baiocchi's Italian GDP panel (see Italy for details), then create a
# data frame in which year is an ordered factor, GDP is continuous,
# compute bandwidths using likelihood cross-validation, then create a
# grid of data on which the density will be evaluated for plotting
# purposes.

data("Italy")
attach(Italy)

data <- data.frame(year=ordered(year), gdp)

# Compute bandwidths using likelihood cross-validation (default), but we
# override the default search tolerances to speed things up as this is a
# well-behaved dataset (don't of course do this in general).

bw <- npudensbw(dat=data, tol=.01, ftol=.01)

# At this stage you could use npudens() to do a variety of
# things. Here we compute the npudens() object and place it in fhat.

fhat <- npudens(bws=bw)

# Note that simply typing the name of the object returns some useful
# information. For more info, one can call summary: 

summary(fhat)

# Next, we illustrate how to create a grid of `evaluation data' and feed
# it to the perspective plotting routines in R, among others.

# Create an evaluation data matrix

year.seq <- sort(unique(year))
gdp.seq <- seq(1,36,length=50)
data.eval <- expand.grid(year=year.seq,gdp=gdp.seq)

# Generate the estimated density computed for the evaluation data

fhat <- fitted(npudens(edat = data.eval, bws=bw))

# Coerce the data into a matrix for plotting with persp()

f <- matrix(fhat, length(unique(year)), 50)

# Next, create a 3D perspective plot of the PDF f, and a 2D
# contour plot.

persp(as.integer(levels(year.seq)), gdp.seq, f, col="lightblue",
      ticktype="detailed", ylab="GDP", xlab="Year", zlab="Density",
      theta=300, phi=50)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

contour(as.integer(levels(year.seq)),
        gdp.seq, 
        f, 
        xlab="Year", 
        ylab="GDP", 
        main = "Density Contour Plot", 
        col=topo.colors(100))

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Alternatively, you could use the npplot() command (<ctrl>-C will
# interrupt on *NIX systems, <esc> will interrupt on MS Windows
# systems).

npplot(bws=bw)

detach(Italy)

# EXAMPLE 2 (INTERFACE=FORMULA): For this example, we load the old
# faithful geyser data and compute the density and distribution
# functions.

library("datasets")
data("faithful")
attach(faithful)

# Note - this may take a few minutes depending on the speed of your
# computer...

bw <- npudensbw(formula=~eruptions+waiting)

summary(bw)

# Plot the density function (<ctrl>-C will interrupt on *NIX systems, 
# <esc> will interrupt on MS Windows systems). Note that we use xtrim =
# -0.2 to extend the plot outside the support of the data (i.e., extend
# the tails of the estimate to meet the horizontal axis).

npplot(bws=bw, xtrim=-0.2)

# Plot the distribution function (cdf=TRUE) (<ctrl>-C will interrupt)

npplot(bws=bw, cdf=TRUE)

detach(faithful)

# EXAMPLE 2 (INTERFACE=DATA FRAME): For this example, we load the old
# faithful geyser data and compute the density and distribution
# functions.

library("datasets")
data("faithful")
attach(faithful)

# Note - this may take a few minutes depending on the speed of your
# computer...

bw <- npudensbw(dat=faithful)

summary(bw)

# Plot the density function (<ctrl>-C will interrupt on *NIX systems, 
# <esc> will interrupt on MS Windows systems). Note that we use xtrim =
# -0.2 to extend the plot outside the support of the data (i.e., extend
# the tails of the estimate to meet the horizontal axis).

npplot(bws=bw, xtrim=-0.2)

# Plot the distribution function (cdf=TRUE) (<ctrl>-C will interrupt)

npplot(bws=bw, cdf=TRUE)

detach(faithful)
} % enddontrun
}
\keyword{ nonparametric }% at least one, from doc/KEYWORDS

