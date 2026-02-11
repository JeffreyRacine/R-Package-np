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
  \code{npplot} is invoked by \code{\link{plot}} and generates plots of
  nonparametric statistical objects such as regressions, quantile
  regressions, partially linear regressions, single-index models,
  densities and distributions, given training data and a bandwidth
  object.
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
       main = NULL,
       type = NULL,
       border = NULL,
       cex.axis = NULL,
       cex.lab = NULL,
       cex.main = NULL,
       cex.sub = NULL,
       col = NULL,
       ylab = NULL,
       xlab = NULL,
       zlab = NULL,
       sub = NULL,
       ylim = NULL,
       xlim = NULL,
       zlim = NULL,
       lty = NULL,
       lwd = NULL,
       theta = 0.0,
       phi = 10.0,
       view = c("rotate","fixed"),
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.boot.num = 399,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("band","bar"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       plot.par.mfrow = TRUE,
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
       gradients = FALSE,
       common.scale = TRUE,
       perspective = TRUE,
       main = NULL,
       type = NULL,
       border = NULL,
       cex.axis = NULL,
       cex.lab = NULL,
       cex.main = NULL,
       cex.sub = NULL,
       col = NULL,
       ylab = NULL,
       xlab = NULL,
       zlab = NULL,
       sub = NULL,
       ylim = NULL,
       xlim = NULL,
       zlim = NULL,
       lty = NULL,
       lwd = NULL,
       theta = 0.0,
       phi = 10.0,
       tau = 0.5,
       view = c("rotate","fixed"),
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.boot.num = 399,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("band","bar"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       plot.par.mfrow = TRUE,
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
       main = NULL,
       type = NULL,
       border = NULL,
       cex.axis = NULL,
       cex.lab = NULL,
       cex.main = NULL,
       cex.sub = NULL,
       col = NULL,
       ylab = NULL,
       xlab = NULL,
       zlab = NULL,
       sub = NULL,
       ylim = NULL,
       xlim = NULL,
       zlim = NULL,
       lty = NULL,
       lwd = NULL,
       theta = 0.0,
       phi = 10.0,
       view = c("rotate","fixed"),
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.boot.num = 399,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("band","bar"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       plot.par.mfrow = TRUE,
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
       main = NULL,
       type = NULL,
       border = NULL,
       cex.axis = NULL,
       cex.lab = NULL,
       cex.main = NULL,
       cex.sub = NULL,
       col = NULL,
       ylab = NULL,
       xlab = NULL,
       zlab = NULL,
       sub = NULL,
       ylim = NULL,
       xlim = NULL,
       zlim = NULL,
       lty = NULL,
       lwd = NULL,
       theta = 0.0,
       phi = 10.0,
       view = c("rotate","fixed"),
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.num = 399,
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("band","bar"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       plot.par.mfrow = TRUE,
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
       main = NULL,
       type = NULL,
       border = NULL,
       cex.axis = NULL,
       cex.lab = NULL,
       cex.main = NULL,
       cex.sub = NULL,
       col = NULL,
       ylab = NULL,
       xlab = NULL,
       zlab = NULL,
       sub = NULL,
       ylim = NULL,
       xlim = NULL,
       zlim = NULL,
       lty = NULL,
       lwd = NULL,
       theta = 0.0,
       phi = 10.0,
       view = c("rotate","fixed"),
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.num = 399,
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("band","bar"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = min(neval,25),
       plot.bxp = FALSE,
       plot.bxp.out = TRUE,
       plot.par.mfrow = TRUE,
       \dots,
       random.seed)

\method{npplot}{sibandwidth}(bws,
       xdat,
       ydat,
       data = NULL,
       common.scale = TRUE,
       gradients = FALSE,
       main = NULL,
       type = NULL,
       cex.axis = NULL,
       cex.lab = NULL,
       cex.main = NULL,
       cex.sub = NULL,
       col = NULL,
       ylab = NULL,
       xlab = NULL,
       sub = NULL,
       ylim = NULL,
       xlim = NULL,
       lty = NULL,
       lwd = NULL,
       plot.behavior = c("plot","plot-data","data"),
       plot.errors.method = c("none","bootstrap","asymptotic"),
       plot.errors.boot.num = 399,
       plot.errors.boot.method = c("inid", "fixed", "geom"),
       plot.errors.boot.blocklen = NULL,
       plot.errors.center = c("estimate","bias-corrected"),
       plot.errors.type = c("standard","quantiles"),
       plot.errors.quantiles = c(0.025,0.975),
       plot.errors.style = c("band","bar"),
       plot.errors.bar = c("|","I"),
       plot.errors.bar.num = NULL,
       plot.par.mfrow = TRUE,
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

  \item{gradients}{
    a logical value specifying whether gradients should be plotted
    (if possible). Defaults to \code{FALSE}.
  }

  \item{main}{
    optional title, see \code{\link{title}}. Defaults to \code{NULL}.
  }

  \item{sub}{
    optional subtitle, see \code{\link{sub}}. Defaults to \code{NULL}.
  }

  \item{type}{
    optional character indicating the type of plotting; actually any of
    the types as in \code{\link{plot.default}}. Defaults to \code{NULL}.
  }

  \item{border}{
    optional character indicating the border of plotting; actually any of
    the borders as in \code{\link{plot.default}}. Defaults to \code{NULL}.
  }

  \item{cex.axis}{
    The magnification to be used for axis annotation relative to the
    current setting of cex.
  }

  \item{cex.lab}{
    The magnification to be used for x and y labels relative to the
    current setting of cex.
  }

  \item{cex.main}{
    The magnification to be used for main titles relative to the
    current setting of cex.
  }

  \item{cex.sub}{
    The magnification to be used for sub-titles relative to the
    current setting of cex.
  }

  \item{col}{
    optional character indicating the color of plotting; actually any of
    the colours as in \code{\link{plot.default}}. Defaults to \code{NULL}.
  }

  \item{ylab}{
    optional character indicating the y axis label of plotting; actually any of
    the ylabs as in \code{\link{plot.default}}. Defaults to \code{NULL}.
  }

  \item{xlab}{
    optional character indicating the x axis label of plotting; actually any of
    the xlabs as in \code{\link{plot.default}}. Defaults to \code{NULL}.
  }

  \item{zlab}{
    optional character indicating the z axis label of plotting; actually any of
    the zlabs as in \code{\link{plot.default}}. Defaults to \code{NULL}.
  }

  \item{ylim}{
    optional a two-element numeric vector of the minimum and maximum y plotting
    limits. Defaults to \code{NULL}.
  }

  \item{xlim}{
    a two-element numeric vector of the minimum and maximum x plotting
    limits. Defaults to \code{NULL}.
  }

  \item{zlim}{
    a two-element numeric vector of the minimum and maximum z plotting
    limits. Defaults to \code{NULL}.
  }

  \item{lty}{
    a numeric value indicating the line type of plotting; actually any of
    the ltys as in \code{\link{plot.default}}. Defaults to \code{1}.
  }

  \item{lwd}{
    a numeric value indicating the width of the line of plotting;
    actually any of the lwds as in \code{\link{plot.default}}. Defaults
    to \code{1}.
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
    desired when plotting quantile regressions. Defaults to \code{0.5}.
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
    be set as \code{band} or \code{bar}. Defaults to \code{band}. Bands
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

  \item{plot.par.mfrow}{
    a logical value specifying whether \code{par(mfrow=c(,))} should be called
    before plotting. Defaults to \code{TRUE}.
  }

  \item{random.seed}{
    an integer used to seed R's random number generator. This ensures
    replicability of the bootstrapped errors. Defaults to 42.
  }

}

\details{
  \code{npplot} is a general purpose plotting routine for visually
  exploring objects generated by the \code{npRmpi} library, such as
  regressions, quantile regressions, partially linear regressions,
  single-index models, densities and distributions. There is no need to
  call \code{npplot} directly as it is automatically invoked when
  \code{\link{plot}} is used with an object generated by the \pkg{npRmpi}
  package.

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
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
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
## Not run in checks: excluded to keep MPI examples stable and check times short.
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi",package="npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.stop(force=TRUE)` then restart.
npRmpi.start(nslaves=1)

x <- rnorm(100)
mpi.bcast.Robj2slave(x)
mpi.bcast.cmd(fhat <- npudens(~x),
              caller.execute=TRUE)

mpi.bcast.cmd(plot(fhat),
              caller.execute=TRUE)

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

## Note: on some systems (notably macOS+MPICH), repeatedly spawning and
## tearing down slaves in the same R session can lead to hangs/crashes.
## npRmpi may therefore keep slave daemons alive by default and
## `npRmpi.stop()` performs a "soft close". Use `force=TRUE` to
## actually shut down the slaves.
##
## You can disable reuse via `options(npRmpi.reuse.slaves=FALSE)` or by
## setting the environment variable `NP_RMPI_NO_REUSE_SLAVES=1` before
## loading the package.

npRmpi.stop()               ## soft close (may keep slaves alive)
## npRmpi.stop(force=TRUE)  ## hard close

## Note that in order to exit npRmpi properly avoid quit(), and instead
## use mpi.quit() as follows.

## mpi.bcast.cmd(mpi.quit(),
##               caller.execute=TRUE)
} % enddontrun     
}

\keyword{nonparametric}
