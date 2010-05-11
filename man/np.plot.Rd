% $Id: np.plot.Rd,v 1.62 2006/11/03 21:17:20 tristen Exp $
\name{npplot}
\alias{npplot}
\alias{npplot.bandwidth}
\alias{npplot.conbandwidth}
\alias{npplot.plbandwidth}
\alias{npplot.rbandwidth}
\alias{npplot.scbandwidth}
\alias{npplot.sibandwidth}

\title{General Purpose Plotting of Nonparametric Objects}

\description{
  \code{npplot} generates plots of nonparametric statistical objects
  such as regressions, quantile regressions, partially linear regressions,
  single-index models, densities and distributions, given
  training data and a bandwidth object.
}

\usage{
npplot(bws = stop("'bws' has not been set"), \dots, random.seed = 42)

\method{npplot}{bandwidth}(bws,
       xdat,
       data = NULL,
       xq = 0.5,
       xtrim = 0.0,
       neval = 50,
       common.scale = TRUE,
       perspective = TRUE,
       main = "",
       theta = 0.0,
       phi = 10.0,
       view = c("rotate","fixed"),
       type = "l",
       ylim = NULL,
       cdf = FALSE,
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.boot.num = 399,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("bar","band"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       \dots,
       random.seed)

\method{npplot}{conbandwidth}(bws,
       xdat,
       ydat,
       data = NULL,
       xq = 0.5,
       yq = 0.5,
       xtrim = 0.0,
       ytrim = 0.0,
       neval = 50,
       quantreg = FALSE,
       gradients = FALSE,
       cdf = FALSE,
       common.scale = TRUE,
       perspective = TRUE,
       main = "",
       theta = 0.0,
       phi = 10.0,
       tau = 0.5,
       view = c("rotate","fixed"),
       type = "l",
       ylim = NULL,
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.boot.num = 399,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("bar","band"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       \dots,
       random.seed)

\method{npplot}{plbandwidth}(bws,
       xdat,
       ydat,
       zdat,
       data = NULL,
       xq = 0.5,
       zq = 0.5,
       xtrim = 0.0,
        ztrim = 0.0,
       neval = 50,
       common.scale = TRUE,
       perspective = TRUE,
       gradients = FALSE,
       main = "",
       theta = 0.0,
       phi = 10.0,
       view = c("rotate","fixed"),
       type = "l",
       ylim = NULL,
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.boot.num = 399,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("bar","band"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       \dots,
       random.seed)

\method{npplot}{rbandwidth}(bws,
       xdat,
       ydat,
       data = NULL,
       xq = 0.5,
       xtrim = 0.0,
       neval = 50,
       common.scale = TRUE,
       perspective = TRUE,
       gradients = FALSE,
       main = "",
       theta = 0.0,
       phi = 10.0,
       view = c("rotate","fixed"),
       type = "l",
       ylim = NULL,
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.num = 399,
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("bar","band"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       \dots,
       random.seed)

\method{npplot}{scbandwidth}(bws,
       xdat,
       ydat,
       zdat = NULL,
       data = NULL,
       xq = 0.5,
       zq = 0.5,
       xtrim = 0.0,
       ztrim = 0.0,
       neval = 50,
       common.scale = TRUE,
       perspective = TRUE,
       gradients = FALSE,
       main = "",
       theta = 0.0,
       phi = 10.0,
       view = c("rotate","fixed"),
       type = "l",
       ylim = NULL,
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.num = 399,
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("bar","band"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       \dots,
       random.seed)

\method{npplot}{sibandwidth}(bws,
       xdat,
       ydat,
       data = NULL,
       common.scale = TRUE,
       gradients = FALSE,
       main = "",
       type = "l",
       ylim = NULL,
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.num = 399,
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("bar","band"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = NULL,
       \dots,
       random.seed)

}
\arguments{

  \item{bws}{
    a bandwidth specification. This should be a bandwidth object
    returned from an invocation of \code{\link{npudensbw}}, \link{npcdensbw},
    \code{\link{npregbw}}, \link{npplregbw}, \link{npindexbw}, or
    \code{\link{npscoefbw}}.
  }

  \item{\dots}{
    additional arguments supplied to control various aspects of
    plotting, depending on the type of object to be plotted, detailed
    below. 
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment where the
    bandwidth object was generated.
  }

  
  \item{xdat}{
    a \eqn{p}-variate data frame of sample realizations (training
    data). 
  }

  \item{ydat}{
    a \eqn{q}-variate data frame of sample realizations (training
    data). In a regression or conditional density context, this is the
    dependent data.  
  }

  \item{zdat}{
    a \eqn{p}-variate data frame of sample realizations (training
    data).
  }

  \item{xq}{
    a numeric \eqn{p}-vector of quantiles. Each element \eqn{i} of
    \code{xq} corresponds to the \eqn{i}th column of
    \code{txdat}. Defaults to the median (0.5). See
    details. 
  }

  \item{yq}{
    a numeric \eqn{q}-vector of quantiles. Each element \eqn{i} of
    \code{yq} corresponds to the \eqn{i}th column of \code{tydat}. Only
    to be specified in a conditional density context. Defaults to the
    median (0.5). See details. 
  }

  \item{zq}{
    a numeric \eqn{q}-vector of quantiles. Each element \eqn{i} of
    \code{zq} corresponds to the \eqn{i}th column of \code{tzdat}. Only
    to be specified in a semiparametric model context. Defaults to the
    median (0.5). See details. 
  }

  \item{xtrim}{
    a numeric \eqn{p}-vector of quantiles. Each element \eqn{i} of
    \code{xtrim} corresponds to the \eqn{i}th column of
    \code{txdat}. Defaults to \code{0.0}. See details.
  }

  \item{ytrim}{
    a numeric \eqn{q}-vector of quantiles. Each element \eqn{i} of
    \code{ytrim} corresponds to the \eqn{i}th column of
    \code{tydat}. Defaults to \code{0.0}. See details.
  }

  \item{ztrim}{
    a numeric \eqn{q}-vector of quantiles. Each element \eqn{i} of
    \code{ztrim} corresponds to the \eqn{i}th column of
    \code{tzdat}. Defaults to \code{0.0}. See details.
  }


  \item{neval}{
    an integer specifying the number of evaluation points. Only applies
    to continuous variables however, as discrete variables will be
    evaluated once at each category. Defaults to \code{50}. 
  }

  \item{common.scale}{
    a logical value specifying whether or not all graphs are to be
    plotted on a common scale. Defaults to \code{TRUE}.
  }

  \item{perspective}{
    a logical value specifying whether a perspective plot should be
    displayed (if possible). Defaults to \code{TRUE}.
  }

  \item{quantreg}{
    a logical value specifying whether a quantile regression should be
    plotted (\code{conbandwidth} objects only). See related parameter
    'tau'. Defaults to \code{FALSE}.
  }

  \item{gradients}{
    a logical value specifying whether gradients should be plotted
    (if possible). Defaults to \code{FALSE}.
  }

  \item{cdf}{
    a logical value specifying whether cumulative distributions should
    be plotted instead of densities (when appropriate). Defaults to
    \code{FALSE}.
  }

  \item{main}{ the title, see \code{\link{title}} }

  \item{type}{
    character indicating the type of plotting; actually any of
    the types as in \code{\link{plot.default}}
  }
  \item{ylim}{
    a two-element numeric vector of the minimum and maximum y plotting
    limits. Defaults to \code{NULL}.
  }

  \item{theta}{
    a numeric value specifying the starting azimuthal angle of the
    perspective plot. Defaults to \code{0.0}.
  }

  \item{phi}{
    a numeric value specifying the starting zenith angle of the
    perspective plot. Defaults to \code{10.0}.
  }

  \item{tau}{
    a numeric value specifying the \eqn{\tau}{tau}th quantile is
    desired when plotting quantile regressions.
  }

  \item{view}{
    a character string used to specify the viewing mode of the
    perspective plot. Can be set as \code{rotate} or
    \code{fixed}. Defaults to \code{rotate}.
  }

  \item{plot.behavior}{
    a character string used to specify the net behavior of \code{npplot}. Can
    be set as \code{plot}, \code{plot-data} or \code{data}. Defaults to
    \code{plot}. See value.
  }

  \item{plot.errors.method}{
    a character string used to specify the method to calculate
    errors. Can be set as \code{none}, \code{bootstrap}, or
    \code{asymptotic}. Defaults to \code{none}.
  }

  \item{plot.errors.boot.method}{
    a character string used to specify the bootstrap method. Can be set
    as \code{inid}, \code{fixed}, or \code{geom} (see
    below for details). Defaults to \code{inid}.
  }

  \item{plot.errors.boot.blocklen}{
    an integer used to specify the block length \eqn{b} for the
    \code{fixed} or \code{geom} bootstrap (see below for
    details).
  }

  \item{plot.errors.boot.num}{
    an integer used to specify the number of bootstrap samples to use
    for the calculation of errors. Defaults to \code{399}.
  }

  \item{plot.errors.center}{
    a character string used to specify where to center the errors on the
    plot(s). Can be set as \code{estimate} or
    \code{bias-corrected}. Defaults to \code{estimate}.
  }

  \item{plot.errors.type}{
    a character string used to specify the type of error to
    calculate. Can be set as \code{standard} or
    \code{quantiles}. Defaults to \code{standard}.
  }

  \item{plot.errors.quantiles}{
    a numeric vector specifying the quantiles of the statistic to
    calculate for the purpose of error plotting. Defaults to
    \code{c(0.025,0.975)}.
  }

  \item{plot.errors.style}{
    a character string used to specify the style of error plotting. Can
    be set as \code{bar} or \code{band}. Defaults to \code{bar}. Bands
    are not drawn for discrete variables.
  }

  \item{plot.errors.bar}{
    a character string used to specify the error bar shape. Can be set
    as \code{|} (vertical bar character) for a dashed vertical bar, or
    as \code{I} for an \sQuote{I} shaped error bar with horizontal
    bounding bars. Defaults to \code{|}.
  }

  \item{plot.errors.bar.num}{
    an integer specifying the number of error bars to plot. Defaults to
    \code{min(neval,25)}. 
  }

  \item{plot.bxp}{
    a logical value specifying whether boxplots should be produced when
    appropriate. Defaults to \code{FALSE}.
  }

  \item{plot.bxp.out}{
    a logical value specifying whether outliers should be plotted on
    boxplots. Defaults to \code{TRUE}.
  }

  \item{random.seed}{
    an integer used to seed R's random number generator. This ensures
    replicability of the bootstrapped errors. Defaults to 42.
  }

}

\details{
  \code{npplot} is a general purpose plotting routine for visually
  exploring objects generated by the \code{np} library, such as
  regressions, quantile regressions, partially linear regressions,
  single-index models, densities and distributions.

  Visualizing one and two dimensional datasets is a straightforward
  process. The default behavior of \code{npplot} is to generate a
  standard 2D plot to visualize univariate data, and a perspective plot
  for bivariate data. When visualizing higher dimensional data,
  \code{npplot} resorts to plotting a series of 1D slices of the
  data. For a slice along dimension \eqn{i}, all other variables at
  indices \eqn{j \ne i}{j != i} are held constant at the quantiles
  specified in the \eqn{j}th element of \code{xq}. The default is the
  median. 

  The slice itself is evaluated on a uniformly spaced sequence of
  \eqn{neval} points. The interval of evaluation is determined by the
  training data. The default behavior is to evaluate from
  \code{min(txdat[,i])} to \code{max(txdat[,i])}. The \code{xtrim}
  variable allows for control over this behavior. When \code{xtrim} is
  set, data is evaluated from the \code{xtrim[i]}th quantile of
  \code{txdat[,i]} to the \code{1.0-xtrim[i]}th quantile of
  \code{txdat[,i]}.

  Furthermore, \code{xtrim} can be set to a negative
  value in which case it will expand the limits of the evaluation
  interval beyond the support of the training data, by measuring the
  distance between \code{min(txdat[,i])} and the \code{xtrim[i]}th
  quantile of \code{txdat[,i]}, and extending the support by that
  distance on the lower limit of the interval. \code{npplot} uses an
  analogous procedure to extend the upper limit of the interval.

  Bootstrap resampling is conducted pairwise on \eqn{(y,X,Z)} (i.e., by
  resampling from rows of the \eqn{(y,X)} data or \eqn{(y,X,Z)} data
  where appropriate). \code{inid} admits general
  heteroskedasticity of unknown form, though it does not allow for
  dependence. \code{fixed} conducts Kunsch's (1988) block bootstrap
  for dependent data, while \code{geom} conducts Politis and
  Romano's (1994) stationary bootstrap.

  For consistency of the block and stationary bootstrap, the (mean)
  block length \eqn{b} should grow with the sample size \eqn{n} at an
  appropriate rate. If \eqn{b} is not given, then a default growth rate
  of \eqn{const \times n^{1/3}}{const n^{1/3}} is used. This rate is
  \dQuote{optimal} under certain conditions (see Politis and Romano
  (1994) for more details). However, in general the growth rate depends on
  the specific properties of the DGP. A default value for \eqn{const}
  (\eqn{3.15}) has been determined by a Monte Carlo simulation using a
  Gaussian AR(1) process (AR(1)-parameter of 0.5, 500
  observations). \eqn{const} has been chosen such that the mean square
  error for the bootstrap estimate of the variance of the empirical mean
  is minimized.

}

\value{
  Setting \code{plot.behavior} will instruct \code{npplot} what data
  to return. Option summary:\cr
  \code{plot}: instruct \code{npplot} to just plot the data and
  return \code{NULL} \cr
  \code{plot-data}: instruct \code{npplot} to plot the data and return
  the data used to generate the plots. The data will be a \code{list} of
  objects of the appropriate type, with one object per plot. For
  example, invoking \code{npplot} on 3D density data will have it
  return a list of three npdensity objects. If biases were calculated,
  they are stored in a component named \code{bias}\cr
  \code{data}: instruct \code{npplot} to generate data only and no plots
}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and J.S. Racine and Q. Li (2004), \dQuote{Cross-validation and the
    estimation of conditional probability densities,} Journal of the
  American Statistical Association, 99, 1015-1026.

  Kunsch, H.R. (1989), \dQuote{The jackknife and the bootstrap for
    general stationary observations,} The Annals of Statistics, 17,
  1217-1241.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,} Cambridge
  University Press. 

  Politis, D.N. and J.P. Romano (1994), \dQuote{The stationary
    bootstrap,} Journal of the American Statistical Association, 89,
  1303-1313.

  Scott, D.W. (1992), \emph{Multivariate Density Estimation. Theory,
  Practice and Visualization,} New York: Wiley.

  Silverman, B.W. (1986), \emph{Density Estimation,} London: Chapman and
  Hall.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,} Biometrika, 68, 301-309.

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

\examples{
\dontrun{
# EXAMPLE 1: For this example, we load Giovanni Baiocchi's Italian GDP
# panel (see Italy for details), then create a data frame in which year
# is an ordered factor, GDP is continuous, compute bandwidths using
# likelihood cross-validation, then create a grid of data on which the
# density will be evaluated for plotting purposes

data("Italy")
attach(Italy)

data <- data.frame(ordered(year), gdp)

# Compute bandwidths using likelihood cross-validation (default). We
# override the default tolerances for the search method as the objective
# function is well-behaved (don't of course do this in general). Note
# that this may take a minute or two depending on the speed of your
# computer...

bw <- npudensbw(dat=data, tol=.1, ftol=.1)

# You can always do things manually, as the following example demonstrates

# Create an evaluation data matrix

year.seq <- sort(unique(year))
gdp.seq <- seq(1,36,length=50)
data.eval <- expand.grid(year=year.seq,gdp=gdp.seq)

# Generate the estimated density computed for the evaluation data

fhat <- fitted(npudens(tdat = data, edat = data.eval, bws=bw))

# Coerce the data into a matrix for plotting with persp()

f <- matrix(fhat, length(unique(year)), 50)

# Next, create a 3D perspective plot of the PDF f

persp(as.integer(levels(year.seq)), gdp.seq, f, col="lightblue",
      ticktype="detailed", ylab="GDP", xlab="Year", zlab="Density",
      theta=300, phi=50)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# However, npplot simply streamlines this process and aids in the
# visualization process (<ctrl>-C will interrupt on *NIX systems, <esc>
# will interrupt on MS Windows systems).

npplot(bws=bw)

# npplot also streamlines construction of variability bounds (<ctrl>-C
# will interrupt on *NIX systems, <esc> will interrupt on MS Windows
# systems)

npplot(bws = bw, plot.errors.method = "asymptotic")

# EXAMPLE 2: For this example, we simulate multivariate data, and plot the
# partial regression surfaces for a locally linear estimator and its
# derivatives.

set.seed(123)

n <- 100

x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
x4 <- rbinom(n, 2, .3)

y <- 1 + x1 + x2 + x3 + x4 + rnorm(n)

X <- data.frame(x1, x2, x3, ordered(x4))

# We override the default tolerances for the search method as the
# objective function is well-behaved (don't of course do this in
# general).

bw <- npregbw(xdat=X, ydat=y, regtype="ll",
              bwmethod="cv.aic", tol=.1, ftol=.1)

npplot(bws=bw)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Now plot the gradients...

npplot(bws=bw, gradients=TRUE)

# Plot the partial regression surfaces with bias-corrected bootstrapped
# nonparametric confidence intervals... this may take a minute or two
# depending on the speed of your computer as the bootstrapping must be
# completed prior to results being displayed...

npplot(bws=bw,
       plot.errors.method="bootstrap", 
       plot.errors.center="bias-corrected",
       plot.errors.type="quantiles")

# Note - if you wished to create, say, a postscript graph for inclusion
# in, say, a latex document, use R's `postscript' command to switch to
# the postscript device (turn off the device once completed). The
# following will create a disk file `graph.ps' that can be pulled into,
# say, a latex document via \includegraphics[width=5in, height=5in,
# angle=270]{graph.ps}

# Note - make sure to include the graphicx package in your latex
# document via adding \usepackage{graphicx} in your latex file. Also, 
# you might was to use larger fonts, which can be achieved by adding the
# pointsize= argument, e.g., postscript(file="graph.ps", pointsize=20)

postscript(file="graph.ps")
npplot(bws=bw)
dev.off()

# The following latex file compiled in the same directory as graph.ps
# ought to work (remove the #s and place in a file names, e.g., 
# test.tex).
# \documentclass[]{article}
# \usepackage{graphicx}
# \begin{document}
# \begin{figure}[!ht]
# \includegraphics[width=5in, height=5in, angle=270]{graph.ps}
# \caption{Local linear partial regression surfaces.}
# \end{figure}
# \end{document}


# EXAMPLE 3: This example demonstrates how to retrieve plotting data from
# npplot(). When npplot() is called with the arguments
# `plot.behavior="plot-data"' (or "data"), it returns plotting objects
# named r1, r2, and so on (rg1, rg2, and so on when `gradients=TRUE' is
# set).  Each plotting object's index (1,2,...) corresponds to the index
# of the explanatory data data frame xdat (and zdat if appropriate). 

# Take the cps71 data by way of example. In this case, there is only one
# object returned by default, `r1', since xdat is univariate.

data("cps71")
attach(cps71)

# Compute bandwidths for local linear regression using cv.aic...

bw <- npregbw(xdat=age,ydat=logwage,regtype="ll",bwmethod="cv.aic")

# Generate the plot and return plotting data, and store output in
# `plot.out' (NOTE: the call to `plot.behavior' is necessary).

plot.out <- npplot(bws=bw,
                   perspective=FALSE,
                   plot.errors.method="bootstrap",
                   plot.errors.boot.num=25,
                   plot.behavior="plot-data")

# Now grab the r1 object that npplot plotted on the screen, and take
# what you need. First, take the output, lower error bound and upper
# error bound...

logwage.eval <- fitted(plot.out$r1)
logwage.se <- se(plot.out$r1)
logwage.lower.ci <- logwage.eval + logwage.se[,1]
logwage.upper.ci <- logwage.eval + logwage.se[,2]

# Next grab the x data evaluation data. xdat is a data.frame(), so we
# need to coerce it into a vector (take the `first column' of data frame
# even though there is only one column)

age.eval <- plot.out$r1$eval[,1]

# Now we could plot this if we wished, or direct it to whatever end use
# we envisioned. We plot the results using R's plot() routines...

plot(age,logwage,cex=0.2,xlab="Age",ylab="log(Wage)")
lines(age.eval,logwage.eval)
lines(age.eval,logwage.lower.ci,lty=3)
lines(age.eval,logwage.upper.ci,lty=3)

# If you wanted npplot() data for gradients, you would use the argument
# `gradients=TRUE' in the call to npplot() as the following
# demonstrates...

plot.out <- npplot(bws=bw,
                   perspective=FALSE,
                   plot.errors.method="bootstrap",
                   plot.errors.boot.num=25,
                   plot.behavior="plot-data",
                   gradients=TRUE)

# Now grab object that npplot() plotted on the screen. First, take the
# output, lower error bound and upper error bound... note that gradients
# are stored in objects rg1, rg2 etc.

grad.eval <- gradients(plot.out$rg1)
grad.se <- gradients(plot.out$rg1, errors = TRUE)
grad.lower.ci <- grad.eval + grad.se[,1]
grad.upper.ci <- grad.eval + grad.se[,2]

# Next grab the x evaluation data. xdat is a data.frame(), so we need to
# coerce it into a vector (take `first column' of data frame even though
# there is only one column)

age.eval <- plot.out$rg1$eval[,1]

# We plot the results using R's plot() routines...

plot(age.eval,grad.eval,cex=0.2,
     ylim=c(min(grad.lower.ci),max(grad.upper.ci)),
     xlab="Age",ylab="d log(Wage)/d Age",type="l")
lines(age.eval,grad.lower.ci,lty=3)
lines(age.eval,grad.upper.ci,lty=3)

detach(cps71)
} % enddontrun
}

\keyword{nonparametric}

