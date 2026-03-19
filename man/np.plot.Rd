% 
\name{plot.np}
\alias{plot.np}

\alias{plot.bandwidth}
\alias{plot.conbandwidth}
\alias{plot.plbandwidth}
\alias{plot.rbandwidth}
\alias{plot.scbandwidth}
\alias{plot.sibandwidth}

\title{General Purpose Plotting of Nonparametric Objects}

\description{
  Plotting is provided via \code{\link{plot}} S3 methods, which generate
  plots of
  nonparametric statistical objects such as regressions, quantile
  regressions, partially linear regressions, single-index models,
  densities and distributions, given training data and a bandwidth
  object. \code{plot(...)} is the supported public interface.
}

\usage{
\method{plot}{bandwidth}(x, \dots)

\method{plot}{conbandwidth}(x, \dots)

\method{plot}{plbandwidth}(x, \dots)

\method{plot}{rbandwidth}(x, \dots)

\method{plot}{scbandwidth}(x, \dots)

\method{plot}{sibandwidth}(x, \dots)

}
\arguments{

  \item{x}{
    a bandwidth specification. This should be a bandwidth object
    returned from an invocation of \code{\link{npudensbw}}, \link{npcdensbw},
    \code{\link{npregbw}}, \link{npplregbw}, \link{npindexbw}, or
    \code{\link{npscoefbw}}.
  }

  \item{\dots}{
    additional arguments supplied to control plotting behavior or passed
    through to underlying plotting helpers where supported. Named options
    passed via \code{\dots} include:

    \describe{
      \item{\code{xdat}, \code{ydat}, \code{zdat}}{
        training data used when plotting bandwidth objects. \code{xdat}
        is the regressor/conditioning data, \code{ydat} is the response
        or conditional outcome data, and \code{zdat} is the auxiliary
        semiparametric covariate block when relevant.
      }
      \item{\code{data}}{
        an optional \code{\link{data.frame}} used to resolve stored
        formula calls when explicit training data arguments are not
        supplied.
      }
      \item{\code{xq}, \code{yq}, \code{zq}}{
        quantiles used to hold non-plotted variables fixed when
        constructing slices.
      }
      \item{\code{xtrim}, \code{ytrim}, \code{ztrim}}{
        trimming/expansion controls for the evaluation support.
      }
      \item{\code{neval}}{
        the number of evaluation points used for each plotted continuous
        margin.
      }
      \item{\code{common.scale}}{
        a logical value specifying whether multiple panels should share a
        common plotting scale when feasible.
      }
      \item{\code{perspective}}{
        a logical value specifying whether bivariate continuous surfaces
        are displayed with \code{\link{persp}}-style plots when relevant.
      }
      \item{\code{gradients}, \code{gradient.order}}{
        derivative plotting controls for families that support
        gradients/derivatives.
      }
      \item{\code{coef}, \code{coef.index}}{
        coefficient-function plotting controls for semiparametric
        families where supported.
      }
      \item{\code{tau}}{
        the quantile index used for quantile-regression-related plotting
        routes when relevant.
      }
      \item{\code{view}, \code{theta}, \code{phi}}{
        view controls for bivariate continuous surface plots.
      }
      \item{\code{plot.behavior}}{
        one of \code{"plot"}, \code{"plot-data"}, or \code{"data"}.
      }
      \item{\code{plot.data.overlay}}{
        logical; when \code{TRUE}, overlay raw data points on regression
        surface plots (currently \code{npreg}, \code{npscoef}, and
        \code{npplreg}) and expand plot limits to include the data. The
        option is ignored for derivative/gradient plots and coefficient
        plots. Overlay points default to small, light markers and can be
        customized with standard \code{points()} arguments passed via
        \code{...} (for example \code{pch}, \code{cex}, \code{col}, and
        \code{bg}). Default is \code{TRUE}.
      }
      \item{\code{plot.errors.method}}{
        one of \code{"none"}, \code{"bootstrap"}, or
        \code{"asymptotic"}.
      }
      \item{\code{plot.errors.boot.method}}{
        the bootstrap resampling scheme. Regression-family plots support
        \code{"wild"}, \code{"inid"}, \code{"fixed"}, and
        \code{"geom"}; density/distribution-family plots support
        \code{"inid"}, \code{"fixed"}, and \code{"geom"}.
      }
      \item{\code{plot.errors.boot.nonfixed}}{
        for non-fixed density/distribution bootstrap routes, one of
        \code{"exact"} or \code{"frozen"}.
      }
      \item{\code{plot.errors.boot.wild}}{
        the wild-bootstrap multiplier distribution, one of
        \code{"rademacher"} or \code{"mammen"}, when relevant.
      }
      \item{\code{plot.errors.boot.blocklen}, \code{plot.errors.boot.num}}{
        the dependent-bootstrap block length and the bootstrap
        replication count.
      }
      \item{\code{plot.errors.center}}{
        one of \code{"estimate"} or \code{"bias-corrected"}.
      }
      \item{\code{plot.errors.type}}{
        one of \code{"pmzsd"}, \code{"pointwise"},
        \code{"bonferroni"}, \code{"simultaneous"}, or \code{"all"}.
      }
      \item{\code{plot.errors.alpha}}{
        the nominal significance level used for interval/band
        construction.
      }
      \item{\code{plot.errors.style}, \code{plot.errors.bar},
      \code{plot.errors.bar.num}}{
        controls for band- versus bar-style uncertainty displays.
      }
      \item{\code{plot.bxp}, \code{plot.bxp.out}}{
        boxplot-style bootstrap summary controls for factor/categorical
        displays when available.
      }
      \item{\code{plot.par.mfrow}}{
        a logical value specifying whether
        \code{par(mfrow=...)} should be managed automatically for
        multi-panel layouts.
      }
      \item{\code{proper}, \code{proper.method}, \code{proper.control}}{
        controls for proper conditional density/distribution handling
        where implemented.
      }
      \item{\code{main}, \code{sub}, \code{xlab}, \code{ylab},
      \code{zlab}, \code{col}, \code{border}, \code{type}, \code{lty},
      \code{lwd}, \code{xlim}, \code{ylim}, \code{zlim},
      \code{cex.axis}, \code{cex.lab}, \code{cex.main},
      \code{cex.sub}}{
        standard graphical controls forwarded to the underlying plotting
        primitives when relevant.
      }
      \item{\code{random.seed}}{
        an optional random seed used to make bootstrap plotting
        reproducible.
      }
    }
  }

}

\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels and
  \code{\link{np.options}} for global options.

  The preferred public interface is \code{\link{plot}} on fitted or
  bandwidth objects (e.g., \code{plot(fit)} or \code{plot(bw)}).

  \code{plot} is a general purpose plotting routine for visually
  exploring objects generated by the \code{np} library, such as
  regressions, quantile regressions, partially linear regressions,
  single-index models, densities and distributions. There is no need to
  call \code{plot} directly: plotting is handled by class-specific S3
  \code{\link{plot}} methods for objects generated by the \pkg{np}
  package.

  Visualizing one and two dimensional datasets is a straightforward
  process. The default behavior of \code{plot} is to generate a
  standard 2D plot to visualize univariate data, and a perspective plot
  for bivariate data. When visualizing higher dimensional data,
  \code{plot} resorts to plotting a series of 1D slices of the
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
  distance on the lower limit of the interval. \code{plot} uses an
  analogous procedure to extend the upper limit of the interval.

  Plot interval/error types are:
  \preformatted{
  "pmzsd"        : point estimate +/- z_(1-alpha/2) * standard error
  "pointwise"    : per-point two-sided interval from N(0,1) quantiles
  "bonferroni"   : pointwise interval with alpha replaced by alpha/m
  "simultaneous" : joint-band interval (bootstrap route)
  }
  where \eqn{m} is the number of evaluation points used in the plotted
  curve/surface (\eqn{m=\texttt{neval}}{m=neval} for univariate curves,
  typically \eqn{m=\texttt{neval}^2}{m=neval^2} for full 2D perspective
  surfaces).

  For asymptotic intervals, let \eqn{T(x)} denote the plotted functional
  (mean, gradient, density, distribution, etc.) and \eqn{\widehat{se}(x)}
  its asymptotic standard error:
  \eqn{T(x)\pm z_{1-\alpha/2}\widehat{se}(x)} for \code{"pmzsd"} and
  \eqn{[T(x)+z_{\alpha/2}\widehat{se}(x),\ T(x)+z_{1-\alpha/2}\widehat{se}(x)]}
  for \code{"pointwise"}.
  \code{"bonferroni"} applies the same pointwise construction with
  \eqn{\alpha/m} in place of \eqn{\alpha}. For the kernel estimators in
  this package, asymptotic simultaneous bands are not generally
  available, so \code{"simultaneous"} with
  \code{plot.errors.method="asymptotic"} returns \code{NA} bands.
  Asymptotic standard errors are taken from fitted-object components such as
  \code{merr}, \code{gerr}, \code{derr}, \code{conderr}, and
  \code{congerr} where implemented.

  Bootstrap resampling is conducted pairwise on \eqn{(y,X,Z)} (i.e., by
  resampling rows of \eqn{(y,X)} or \eqn{(y,X,Z)} as appropriate).
  Bootstrap method support differs by estimator family:
  \preformatted{
  Regression-family (npreg/npindex/npscoef/npplreg): wild, inid, fixed, geom
  Density/distribution-family (npudens/npudist/npcdens/npcdist): inid, fixed, geom
  }
  hence \code{"wild"} is only available for regression-family plotting.

  Implementation notes for speed:
  \preformatted{
  wild            : fast np*hat linear-operator bootstrap path
  inid/fixed/geom : fast direct helper path (no internal bandwidth search)
  }
  For non-fixed density/distribution bootstrap, an explicit experimental
  approximation is available via
  \code{plot.errors.boot.nonfixed=c("exact","frozen")}. The default
  \code{"exact"} route recomputes the non-fixed geometry for each resample.
  The experimental \code{"frozen"} route reuses the original-sample
  non-fixed geometry throughout the bootstrap run. This option is currently
  implemented only for unconditional and conditional density/distribution
  bootstrap routes and remains off by default. For generalized/adaptive
  nearest-neighbor runs, \code{"frozen"} is an approximation that can alter
  interval/band width by holding the original-sample nearest-neighbor
  geometry fixed across bootstrap resamples; \code{"exact"} remains the
  recommended setting for production inference. This approximation can be
  more noticeable for conditional density/distribution plotting than for the
  regression-style plot families because the conditional bootstrap paths
  freeze both numerator and denominator nearest-neighbor geometry before
  recombining them. In practice, conditional distribution bands are often
  closer, while conditional density bands can differ more materially from
  \code{"exact"} under generalized/adaptive nearest-neighbor bandwidths.
  For smooth coefficient plots (\code{npscoef}) under non-fixed bandwidths,
  \code{"exact"} can also be much more expensive than \code{"frozen"} on
  large jobs, because the coefficient field must be recomputed for each
  bootstrap resample rather than reusing the original-sample geometry.
  \code{inid} admits general heteroskedasticity of unknown form, though
  it does not allow for dependence. \code{fixed} conducts Kunsch's (1988)
  block bootstrap for dependent data, while \code{geom} conducts Politis
  and Romano's (1994) stationary bootstrap.

  For local polynomial conditional density/distribution plotting
  (\code{npcdens}/\code{npcdist} with \code{regtype="ll"} or
  \code{regtype="lp"}) and \code{proper=TRUE}, the plotted estimate is
  rendered proper slice-by-slice on the fixed evaluation grid: each
  conditional density slice is projected to be nonnegative and to integrate
  to one using trapezoidal quadrature weights from the evaluation
  \eqn{y}-grid, while each conditional distribution slice is projected to
  be monotone and bounded in \eqn{[0,1]}. When
  \code{plot.errors.method="bootstrap"}, the bootstrap resample surfaces
  are computed first on that same fixed grid and then properized
  resample-by-resample using the same grid geometry before
  \code{"pointwise"}, \code{"bonferroni"}, \code{"simultaneous"}, and
  \code{"all"} bands are constructed. Thus the bootstrap distribution used
  to form these bands is built from properized resample surfaces. The final
  lower/upper band surfaces are interval envelopes and are not themselves
  separately re-projected to satisfy the density/distribution shape
  constraints.

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

  The default bootstrap replication count is
  \code{plot.errors.boot.num=1999}. For pointwise tails, ensure
  \eqn{B \ge \lceil 2/\alpha - 1 \rceil}{B >= ceiling(2/alpha - 1)} so
  \eqn{\alpha(B+1)} is feasible on the bootstrap rank grid. For interval
  types \code{"bonferroni"}, \code{"simultaneous"}, and \code{"all"},
  the minimum recommended count is
  \eqn{B_{\min}=\lceil 2m/\alpha-1 \rceil}{Bmin=ceiling(2*m/alpha-1)},
  where \eqn{m} is the number of evaluation points used by the plotted
  curve/surface. For full 2D perspective grids this is typically
  \eqn{m=\texttt{neval}^2}{m=neval^2}. When \code{B} is below these
  thresholds, plotting proceeds but warning guidance is reported.

  Typical plotting calls:
  \preformatted{
  ## Asymptotic pointwise/bonferroni intervals
  plot(bw, plot.errors.method="asymptotic", plot.errors.type="pointwise")
  plot(bw, plot.errors.method="asymptotic", plot.errors.type="bonferroni")

  ## Regression-family bootstrap (wild available)
  plot(bw, plot.errors.method="bootstrap", plot.errors.boot.method="wild")

  ## Density/distribution-family bootstrap (use inid/fixed/geom)
  plot(bw, plot.errors.method="bootstrap", plot.errors.boot.method="inid")
  }

}

\value{
  Setting \code{plot.behavior} will instruct \code{plot} what data
  to return. Option summary:\cr
  \code{plot}: instruct \code{plot} to just plot the data and
  return \code{NULL} \cr
  \code{plot-data}: instruct \code{plot} to plot the data and return
  the data used to generate the plots. The data will be a \code{list} of
  objects of the appropriate type, with one object per plot. For
  example, invoking \code{plot} on 3D density data will have it
  return a list of three npdensity objects. If biases were calculated,
  they are stored in a component named \code{bias}\cr
  \code{data}: instruct \code{plot} to generate data only and no plots
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
# EXAMPLE 1: For this example, we load Giovanni Baiocchi's Italian GDP
# panel (see Italy for details), then create a data frame in which year
# is an ordered factor, GDP is continuous, compute bandwidths using
# likelihood cross-validation, then create a grid of data on which the
# density will be evaluated for plotting purposes

data("Italy")
attach(Italy)

data <- data.frame(ordered(year), gdp)

# Compute bandwidths using likelihood cross-validation (default). Note
# that this may take a minute or two depending on the speed of your
# computer...

bw <- npudensbw(dat=data)

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

# However, plot simply streamlines this process and aids in the
# visualization process (<ctrl>-C will interrupt on *NIX systems, <esc>
# will interrupt on MS Windows systems).

plot(bw)

# plot also streamlines construction of variability bounds (<ctrl>-C
# will interrupt on *NIX systems, <esc> will interrupt on MS Windows
# systems)

plot(bw, plot.errors.method = "asymptotic")

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

bw <- npregbw(xdat=X, ydat=y, regtype="ll", bwmethod="cv.aic")

plot(bw)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Now plot the gradients...

plot(bw, gradients=TRUE)

# Plot the partial regression surfaces with bias-corrected bootstrapped
# nonparametric confidence intervals... this may take a minute or two
# depending on the speed of your computer as the bootstrapping must be
# completed prior to results being displayed...

plot(bw,
     plot.errors.method="bootstrap", 
     plot.errors.center="bias-corrected",
     plot.errors.type="simultaneous")

# EXAMPLE 3: This example demonstrates how to retrieve plotting data from
# plot(). When plot() is called with the arguments
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

plot.out <- plot(bw,
                 perspective=FALSE,
                 plot.errors.method="bootstrap",
                 plot.errors.boot.num=25,
                 plot.behavior="plot-data")

# Now grab the r1 object that plot plotted on the screen, and take
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

# If you wanted plot() data for gradients, you would use the argument
# `gradients=TRUE' in the call to plot() as the following
# demonstrates...

plot.out <- plot(bw,
                 perspective=FALSE,
                 plot.errors.method="bootstrap",
                 plot.errors.boot.num=25,
                 plot.behavior="plot-data",
                 gradients=TRUE)

# Now grab object that plot() plotted on the screen. First, take the
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

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}
  \code{\link{np.options}}
}

\keyword{nonparametric}
