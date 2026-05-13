\name{npcopula}
\alias{npcopula}
\title{ Kernel Copula Estimation with Mixed Data Types }

\description{
  \code{npcopula} implements the nonparametric mixed data kernel copula
  approach of Racine (2015) for an arbitrary number of dimensions
}

\usage{
npcopula(bws, ...)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification and source data.
}

  \item{bws}{
    an unconditional joint distribution (\code{npudistbw}) bandwidth
    object, an unconditional joint density (\code{npudensbw}) bandwidth
    object, or a one-sided formula. If a formula is supplied,
    \code{npcopula} first computes the required bandwidth object using
    \code{npudistbw} when \code{target="distribution"} and
    \code{npudensbw} when \code{target="density"}.
  }

  \item{data}{
    a data frame containing variables used to construct \code{bws}; when
    \code{bws} is a formula, \code{data} is passed through to the
    bandwidth selector.
  }

  \item{...}{
    additional arguments passed to \code{npudistbw} or \code{npudensbw}
    when \code{bws} is a formula.
  }


\subsection{Copula Evaluation Grid}{
These arguments control the marginal probability grid and numerical inversion used for copula evaluation.
}

  \item{er.quasi.inv}{
    number passed to \code{\link{extendrange}} when \code{u} is provided
    specifying the fraction by which the data range should be extended when
    constructing the grid used to compute the quasi-inverse (see details)
  }

  \item{evaluation}{
    evaluation route used when \code{bws} is a formula. The default
    \code{"grid"} constructs a plot-ready probability grid for
    two-dimensional copulas. Use \code{"sample"} to evaluate the copula
    at the sample realizations.
  }

  \item{n.quasi.inv}{
    number of grid points generated when \code{u} is provided in order to
    compute the quasi-inverse of each marginal distribution (see
    details)
  }

  \item{neval}{
    number of probability values per margin used for the automatically
    generated two-dimensional grid when \code{bws} is a formula and
    \code{evaluation="grid"}.
  }

  \item{target}{
    target used when \code{bws} is a formula. The default
    \code{"distribution"} estimates the copula distribution using
    \code{npudistbw}; \code{"density"} estimates the copula density using
    \code{npudensbw}. When \code{bws} is already a bandwidth object, the
    target is inferred from that object and a conflicting explicit
    \code{target} is rejected.
  }

  \item{u}{
    an optional matrix of real numbers lying in [0,1], each column of
    which corresponds to the vector of uth quantile values desired for
    each variable in the copula (otherwise the u values returned are
    those corresponding to the sample realizations)
  }


\subsection{Copula Plotting}{
These arguments are used by the \code{plot.npcopula} method for two-dimensional grid output.
}

  \item{alpha}{
    nominal size used for asymptotic or bootstrap intervals when
    \code{errors!="none"}.
  }

  \item{B}{
    number of bootstrap replications when \code{errors="bootstrap"}.
  }

  \item{band}{
    interval type for plotted surfaces. Supported values are
    \code{"pointwise"}, \code{"pmzsd"}, \code{"bonferroni"},
    \code{"simultaneous"}, and \code{"all"}. The \code{"all"} option
    overlays pointwise, simultaneous, and Bonferroni wireframes where
    available.
  }

  \item{border}{
    border color for the surface facets when \code{view="surface"}.
  }

  \item{bootstrap}{
    bootstrap resampling method used when \code{errors="bootstrap"};
    supported values are \code{"inid"}, \code{"fixed"}, and
    \code{"geom"}. Wild bootstrap is not defined for copula surfaces.
  }

  \item{boot_control}{
    optional \code{\link{np_boot_control}} object. For copula surfaces
    the block length is used by the \code{"fixed"} and \code{"geom"}
    block bootstrap routes.
  }

  \item{center}{
    centering convention for bootstrap intervals, either
    \code{"estimate"} or \code{"bias-corrected"}.
  }

  \item{col}{
    optional surface or image colors. If omitted, \code{plot.npcopula}
    uses the same viridis \code{\link{hcl.colors}} palette used by the
    other modern surface plot methods.
  }

  \item{errors}{
    interval route for \code{plot.npcopula}: \code{"none"},
    \code{"asymptotic"}, or \code{"bootstrap"}. Intervals are available
    for two-dimensional grid evaluation output and are drawn as
    transparent wireframes over the copula surface.
  }

  \item{legend}{
    logical or legend-control value used for \code{band="all"} interval
    overlays.
  }

  \item{output}{
    plot return mode: \code{"plot"} draws the plot, \code{"data"}
    returns the plotted data with interval columns when requested, and
    \code{"plot-data"} draws and returns the data.
  }

  \item{renderer}{
    plotting renderer for surface displays. \code{"base"} uses
    \code{\link{persp}}; \code{"rgl"} uses the shared package
    \code{rgl} surface renderer when the suggested package \pkg{rgl} is
    installed.
  }

  \item{zlim}{
    optional z-axis limits for surface displays.
  }

}
\value{

  \code{npcopula} returns an object of class \code{"npcopula"} and
  \code{\link{data.frame}} with the following components

  \item{copula}{
    the copula (bandwidth object obtained from \code{npudistbw}) or
    copula density (bandwidth object obtained from \code{npudensbw})
  }

  \item{u}{ the matrix of marginal u values associated with the sample
  realizations (\code{u=NULL}) or those created via
  \code{\link{expand.grid}} when \code{u} is provided}

  \item{data}{ the matrix of marginal quantiles constructed when
  \code{u} is provided (\code{data} returned has the same names as
  \code{data} inputted)}

  The selected bandwidth object and display metadata are retained as
  attributes. Use \code{attr(object, "bws")} to retrieve the bandwidth
  object selected by the formula route.

}
\references{

  Nelsen, R. B. (2006), \emph{An Introduction to Copulas,} Second
  Edition, Springer-Verlag.

  Racine, J.S. (2015), \dQuote{Mixed Data Kernel Copulas,} Empirical
  Economics, 48, 37-59.

}
\author{
    Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, \code{\link{plot}} for plotting options, and \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See \code{\link{npRmpi.init}} details for performance tradeoffs (message passing/startup mode) and the \code{inst/Rprofile} manual-broadcast template.


  \code{npcopula} computes the nonparametric copula or copula density
  using inversion (Nelsen (2006), page 51). For the inversion approach,
  we exploit Sklar's theorem (Corollary 2.3.7, Nelsen (2006)) to produce
  copulas directly from the joint distribution function using
  \eqn{C(u,v) = H(F^{-1}(u),G^{-1}(v))} rather than the typical approach
  that instead uses \eqn{H(x,y) = C(F(x),G(y))}. Whereas the latter
  requires kernel density estimation on a d-dimensional unit hypercube
  which necessitates the use of boundary correction methods, the former
  does not.

  Note that if \code{u} is provided then \code{\link{expand.grid}} is
  called on \code{u}. As the dimension increases this can become
  unwieldy and potentially consume an enormous amount of memory unless
  the number of grid points is kept very small. Given that computing the
  copula on a grid is typically done for graphical purposes, providing
  \code{u} is typically done for two-dimensional problems only. Even
  here, however, providing a grid of length 100 will expand into a
  matrix of dimension 10000 by 2 which, though not memory intensive, may
  be computationally burdensome.

  If \code{bws} is a formula and \code{u} is omitted, \code{npcopula}
  automatically creates a two-dimensional probability grid when
  \code{evaluation="grid"}. For higher-dimensional copulas, either
  supply \code{u} explicitly or use \code{evaluation="sample"} to avoid
  accidental grid expansion.

  The \code{plot.npcopula} method supports base \code{\link{persp}},
  \code{\link{contour}}, and \code{\link{image}} displays for
  two-dimensional grid output. Surface plots use the package-wide
  viridis default palette; \code{renderer="rgl"} requests the shared
  interactive \pkg{rgl} surface renderer.

  For grid surfaces, \code{plot.npcopula} can add asymptotic or
  bootstrap intervals via \code{errors}. Distribution-copula
  asymptotic intervals use the joint distribution standard error
  evaluated at the marginal quasi-inverse grid. Density-copula
  asymptotic intervals use the plug-in delta-method denominator
  correction corresponding to the estimator
  \eqn{c(u)=f(x_u)/\prod_j f_j(x_{u,j})}. Bootstrap intervals resample
  rows and recompute the plotted copula surface on the same probability
  grid; \code{band="all"} overlays transparent pointwise,
  simultaneous, and Bonferroni wireframes.

  The \sQuote{quasi-inverse} is computed via Definition 2.3.6 from
  Nelsen (2006). We compute an equi-quantile grid on the range of the
  data of length \code{n.quasi.inv/2}. We then extend the range of the
  data by the factor \code{er.quasi.inv} and compute an equi-spaced grid
  of points of length \code{n.quasi.inv/2} (e.g. using the default
  \code{er.quasi.inv=1} we go from the minimum data value minus
  \eqn{1\times}{1x} the range to the maximum data value plus
  \eqn{1\times}{1x} the range for each marginal). We then take these two
  grids, concatenate and sort, and these form the final grid of length
  \code{n.quasi.inv} for computing the quasi-inverse.

  Note that if \code{u} is provided and any elements of (the columns of)
  \code{u} are such that they lie beyond the respective values of
  \code{F} for the evaluation data for the respective marginal, such
  values are reset to the minimum/maximum values of \code{F} for the
  respective marginal. It is therefore prudent to inspect the values of
  \code{u} returned by \code{\link{npcopula}} when \code{u} is provided.

  Note that copula are only defined for data of type
  \code{\link{numeric}} or \code{\link{ordered}}.

}
\section{Usage Issues}{

  See the example below for proper usage.

}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}
  \link{npudensbw},\link{npudens},\link{npudist}
}

\examples{
\dontrun{
## Not run in checks: excluded to keep MPI examples stable and check times short.
## Example 1: Bivariate Mixed Data

## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi_getting_started", package = "npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.quit(force=TRUE)` then restart.

force.run <- nzchar(Sys.getenv("NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK"))

in.check <- !force.run && (
            nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
            nzchar(Sys.getenv("_R_CHECK_INTERNALS2_")) ||
            nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_")) ||
            nzchar(Sys.getenv("_R_CHECK_EXAMPLE_TIMING_THRESHOLD_"))
          )

if (!in.check) {
npRmpi.init(nslaves=1)

library(MASS)

## Minimal bivariate continuous example using the formula interface

set.seed(42)
n <- 1000
rho <- 0.95
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
dat <- as.data.frame(mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma))
names(dat) <- c("x", "y")

copula_fit <- npcopula(~ x + y, data = dat, neval = 20)
summary(copula_fit)
plot(copula_fit)

## Reuse the selected bandwidth object for a custom probability grid
u_grid <- data.frame(x = seq(0, 1, length.out = 10),
                     y = seq(0, 1, length.out = 10))
copula_custom <- npcopula(bws = attr(copula_fit, "bws"),
                          data = dat,
                          u = u_grid)

set.seed(42)

## Simulate correlated Gaussian data (rho(x,y)=0.99)

n <- 300
n.eval <- 30
rho <- 0.99
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
xy <- mvrnorm(n=n, mu, Sigma)
x <- xy[,1]
y <- ordered(as.integer(cut(xy[,2],
                      quantile(xy[,2],seq(0,1,by=.1)),
                      include.lowest=TRUE))-1)
mydat <- data.frame(x=x, y=y)
q.min <- 0.0
q.max <- 1.0
grid.seq <- seq(q.min,q.max,length=n.eval)
grid.dat <- cbind(grid.seq,grid.seq)

## Estimate the copula (bw object obtained from npudistbw())

bw.cdf <- npudistbw(~x+y, data=mydat)
copula <- npcopula(bws=bw.cdf, data=mydat, u=grid.dat)

## Plot the copula


contour(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
        xlab="u1",
        ylab="u2",
        main="Copula Contour")

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
      ticktype="detailed",
      xlab="u1",
      ylab="u2",
      zlab="Copula",zlim=c(0,1))

## Plot the empirical copula

copula.emp <- npcopula(bws=bw.cdf, data=mydat)
if (interactive()) plot(copula.emp$u1,
                        copula.emp$u2,
                        xlab="u1",
                        ylab="u2",
                        cex=.25,
                        main="Empirical Copula")

## Estimate the copula density (bw object obtained from npudensbw())

bw.pdf <- npudensbw(~x+y, data=mydat)
copula <- npcopula(bws=bw.pdf, data=mydat, u=grid.dat)

## Plot the copula density

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
      ticktype="detailed",
      xlab="u1",
      ylab="u2",
      zlab="Copula Density")

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

## Note: on some systems (notably macOS+MPICH), repeatedly spawning and
## tearing down slaves in the same R session can lead to hangs/crashes.
## npRmpi may therefore keep slave daemons alive by default and
## `npRmpi.quit()` performs a "soft close". Use `force=TRUE` to
## actually shut down the slaves.
##
## You can disable reuse via `options(npRmpi.reuse.slaves=FALSE)` or by
## setting the environment variable `NP_RMPI_NO_REUSE_SLAVES=1` before
## loading the package.

npRmpi.quit()               ## soft close (may keep slaves alive)
## npRmpi.quit(force=TRUE)  ## hard close

## Example 2: Bivariate Continuous Data

## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi_getting_started", package = "npRmpi").

## Start npRmpi for interactive execution. If slaves are already running and
## `options(npRmpi.reuse.slaves=TRUE)` (default on some systems), this will
## reuse the existing pool instead of respawning. To change the number of
## slaves, call `npRmpi.quit(force=TRUE)` then restart.
npRmpi.init(nslaves=1)

library(MASS)

set.seed(42)

## Simulate correlated Gaussian data (rho(x,y)=0.99)

n <- 300
n.eval <- 30
rho <- 0.99
mu <- c(0,0)
Sigma <- matrix(c(1,rho,rho,1),2,2)
xy <- mvrnorm(n=n, mu, Sigma)
x <- xy[,1]
y <- xy[,2]
mydat <- data.frame(x=x, y=y)

q.min <- 0.0
q.max <- 1.0
grid.seq <- seq(q.min,q.max,length=n.eval)
grid.dat <- cbind(grid.seq,grid.seq)

## Estimate the copula (bw object obtained from npudistbw())

bw.cdf <- npudistbw(~x+y, data=mydat)
copula <- npcopula(bws=bw.cdf, data=mydat, u=grid.dat)

## Plot the copula

contour(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
        xlab="u1",
        ylab="u2",
        main="Copula Contour")

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),
      ticktype="detailed",
      xlab="u1",
      ylab="u2",
      zlab="Copula",
      zlim=c(0,1))

## Plot the empirical copula

copula.emp <- npcopula(bws=bw.cdf, data=mydat)
if (interactive()) plot(copula.emp$u1,
                        copula.emp$u2,
                        xlab="u1",
                        ylab="u2",
                        cex=.25,
                        main="Empirical Copula")

## Estimate the copula density (bw object obtained from npudensbw())

bw.pdf <- npudensbw(~x+y, data=mydat)
copula <- npcopula(bws=bw.pdf, data=mydat, u=grid.dat)

## Plot the copula density

persp(grid.seq,grid.seq,matrix(copula$copula,n.eval,n.eval),ticktype="detailed",xlab="u1",
      ylab="u2",zlab="Copula Density")

## For the interactive run only we close the slaves perhaps to proceed
## with other examples and so forth. This is redundant in batch mode.

## Note: on some systems (notably macOS+MPICH), repeatedly spawning and
## tearing down slaves in the same R session can lead to hangs/crashes.
## npRmpi may therefore keep slave daemons alive by default and
## `npRmpi.quit()` performs a "soft close". Use `force=TRUE` to
## actually shut down the slaves.
##
## You can disable reuse via `options(npRmpi.reuse.slaves=FALSE)` or by
## setting the environment variable `NP_RMPI_NO_REUSE_SLAVES=1` before
## loading the package.

npRmpi.quit()               ## soft close (may keep slaves alive)
## npRmpi.quit(force=TRUE)  ## hard close
} else {
  message("Skipping MPI spawn in check context (set NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1).")
}
} % enddontrun
}

\keyword{ nonparametric }
