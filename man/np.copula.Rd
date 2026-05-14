\name{npcopula}
\alias{npcopula}
\alias{npcopula.default}
\alias{npcopula.formula}
\alias{plot.npcopula}
\alias{print.npcopula}
\alias{summary.npcopula}
\alias{fitted.npcopula}
\alias{predict.npcopula}
\alias{se.npcopula}
\alias{as.data.frame.npcopula}
\title{Kernel Copula Estimation with Mixed Data Types}

\description{
  \code{npcopula} estimates a mixed-data kernel copula distribution or
  copula density. It can be called with an existing unconditional
  distribution/density bandwidth object, or with a one-sided formula in
  which case the appropriate bandwidth object is selected internally.
}

\usage{
npcopula(bws, \dots)

\method{npcopula}{formula}(bws,
         data = NULL,
         u = NULL,
         target = c("distribution", "density"),
         evaluation = c("grid", "sample"),
         neval = 30,
         n.quasi.inv = 1000,
         er.quasi.inv = 1,
         \dots)

\method{npcopula}{default}(bws,
         data,
         u = NULL,
         target = NULL,
         evaluation = c("sample", "grid"),
         neval = 30,
         n.quasi.inv = 1000,
         er.quasi.inv = 1,
         \dots)

\method{predict}{npcopula}(object,
        newdata = NULL,
        u = NULL,
        se.fit = FALSE,
        output = c("vector", "object", "data"),
        \dots)

\method{plot}{npcopula}(x,
     perspective = TRUE,
     view = c("rotate", "fixed", "surface", "contour", "image",
              "empirical", "all"),
     renderer = c("base", "rgl"),
     errors = c("none", "bootstrap", "asymptotic"),
     band = c("pointwise", "pmzsd", "bonferroni",
              "simultaneous", "all"),
     alpha = 0.05,
     bootstrap = c("inid", "fixed", "geom"),
     B = 1999,
     center = c("estimate", "bias-corrected"),
     boot_control = np_boot_control(),
     output = c("plot", "data", "plot-data", "both"),
     legend = TRUE,
     theta = 0.0,
     phi = 20.0,
     xlab = "u1",
     ylab = "u2",
     zlab = NULL,
     main = NULL,
     col = NULL,
     border = "black",
     zlim = NULL,
     \dots)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, source data, and
formula route.
}

  \item{bws}{
    bandwidth specification or one-sided formula. A bandwidth
    specification is either an unconditional distribution bandwidth
    object returned by \code{\link{npudistbw}}, for a copula
    distribution, or an unconditional density bandwidth object returned
    by \code{\link{npudensbw}}, for a copula density. If \code{bws} is a
    formula such as \code{~ x + y}, \code{npcopula} first calls
    \code{\link{npudistbw}} when \code{target="distribution"} and
    \code{\link{npudensbw}} when \code{target="density"}; additional
    bandwidth-selection arguments in \code{\dots} are forwarded to that
    bandwidth selector.
  }

  \item{data}{
    data frame containing the variables used to construct \code{bws};
    when \code{bws} is a formula, \code{data} is passed to the
    bandwidth selector. Copulas are defined here for numeric and ordered
    variables; unordered factors are rejected.
  }

\subsection{Copula Target And Evaluation Grid}{
These arguments control whether a copula distribution or density is
estimated and where it is evaluated.
}

  \item{er.quasi.inv}{
    fraction passed to \code{\link{extendrange}} when constructing the
    marginal quasi-inverse grid used for supplied probability values
    \code{u}. See Details.
  }

  \item{evaluation}{
    evaluation route used when \code{bws} is a formula. The default
    \code{"grid"} constructs a plot-ready two-dimensional probability
    grid when \code{u} is omitted. Use \code{"sample"} to evaluate at
    the sample realizations.
  }

  \item{n.quasi.inv}{
    number of grid points used to compute each marginal quasi-inverse
    when \code{u} is supplied or automatically generated.
  }

  \item{newdata}{
    optional prediction data for \code{predict.npcopula}. This is a
    compatibility alias for \code{u}: it should contain marginal
    probability values in \eqn{[0,1]}, either in columns named as the
    original variables or in columns named \code{u1}, \code{u2}, ....
    If both \code{u} and \code{newdata} are supplied, \code{u} takes
    precedence.
  }

  \item{neval}{
    number of probability values per margin in the automatically
    generated two-dimensional grid used when \code{bws} is a formula,
    \code{evaluation="grid"}, and \code{u} is omitted.
  }

  \item{target}{
    target used when \code{bws} is a formula. The default
    \code{"distribution"} estimates a copula distribution using
    \code{\link{npudistbw}}; \code{"density"} estimates a copula
    density using \code{\link{npudensbw}}. When \code{bws} is already a
    bandwidth object, the target is inferred from that object and a
    conflicting explicit \code{target} is rejected.
  }

  \item{u}{
    optional matrix or data frame of marginal probability values in
    \eqn{[0,1]}. Each column corresponds to one variable in the copula.
    If supplied, \code{npcopula} evaluates on the Cartesian product
    \code{expand.grid(u)}. For two-dimensional displays, the clearest
    spelling is often \code{data.frame(u1 = u1.seq, u2 = u2.seq)};
    names matching the original variables are also accepted for
    compatibility. If omitted with a formula route and
    \code{evaluation="grid"}, a two-dimensional grid is generated
    automatically.
  }

  \item{object}{
    an object of class \code{"npcopula"} returned by \code{npcopula}.
  }

  \item{se.fit}{
    logical value. If \code{TRUE}, \code{predict.npcopula} returns a
    list with fitted values and the stored standard-error slot.
  }

\subsection{Plot Display Controls}{
These arguments control how a two-dimensional grid copula object is
displayed.
}

  \item{x}{
    an object of class \code{"npcopula"} returned by \code{npcopula}.
  }

  \item{border}{
    border color for surface facets and interval wireframes when
    \code{view="surface"}.
  }

  \item{col}{
    optional surface or image colors. If omitted, \code{plot.npcopula}
    uses the same viridis \code{\link{hcl.colors}} palette used by the
    other modern surface plot methods.
  }

  \item{legend}{
    logical or legend-control value used for \code{band="all"} interval
    overlays. Use \code{FALSE}, \code{NULL}, or \code{NA} to suppress
    the interval legend.
  }

  \item{main, xlab, ylab, zlab}{
    plot titles and axis labels.
  }

  \item{output}{
    return mode. For \code{plot.npcopula}, \code{"plot"} draws the
    plot, \code{"data"} returns the plotted data with interval columns
    when requested, and \code{"plot-data"} draws and returns the
    plotted data. \code{"both"} is accepted as an alias for
    \code{"plot-data"}. For \code{predict.npcopula}, \code{"vector"}
    returns fitted copula values, \code{"object"} returns the evaluated
    \code{"npcopula"} object, and \code{"data"} returns
    \code{as.data.frame()} on that object.
  }

  \item{perspective}{
    logical value. If \code{TRUE}, draw a surface display; otherwise
    use \code{view} to choose \code{"contour"} or \code{"image"}.
    The \code{"empirical"} and \code{"all"} views ignore this argument.
  }

  \item{phi, theta}{
    viewing angles passed to \code{\link{persp}} or to the shared
    \pkg{rgl} surface renderer. The defaults match the package-wide
    surface-plot defaults used by the other perspective plot methods.
    As with those methods, the exact default pair \code{theta = 0} and
    \code{phi = 20} is remapped internally for \code{renderer="rgl"} to
    account for the different viewing-angle convention used by \pkg{rgl};
    explicitly supplied non-default angles are passed through.
  }

  \item{renderer}{
    plotting renderer for surface displays. \code{"base"} uses
    \code{\link{persp}}. \code{"rgl"} uses the shared interactive
    \pkg{rgl} surface renderer when the suggested package \pkg{rgl} is
    installed.
  }

  \item{view}{
    display type for grid output. The default \code{"rotate"} draws a
    rotating base \code{\link{persp}} surface using the same frame step
    and delay as the other package perspective plots. Use
    \code{"fixed"} for a single fixed surface, \code{"surface"} as a
    backward-compatible fixed-surface alias, or \code{"contour"} and
    \code{"image"} with \code{perspective=FALSE}. Use
    \code{"empirical"} to plot the empirical copula coordinates. Use
    \code{"all"} for a base-graphics four-panel display containing the
    copula contour, copula surface, empirical copula coordinates, and a
    copula-density surface. The \code{"all"} view is not currently
    supported with \code{renderer="rgl"}.
  }

  \item{zlim}{
    optional z-axis limits for surface displays.
  }

\subsection{Plot Interval Controls}{
These arguments add asymptotic or bootstrap intervals to two-dimensional
surface plots.
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

  \item{errors}{
    interval route for \code{plot.npcopula}: \code{"none"},
    \code{"asymptotic"}, or \code{"bootstrap"}. Intervals are available
    for two-dimensional grid evaluation output and are drawn as
    transparent wireframes over the copula surface.
  }

\subsection{Additional Arguments}{
Further arguments are passed to the bandwidth-selection counterpart,
prediction/evaluation route, or graphics renderer as appropriate.
}

  \item{\dots}{
    additional arguments supplied to \code{\link{npudistbw}} or
    \code{\link{npudensbw}} when \code{npcopula} computes bandwidths
    internally, or arguments needed to interpret a numeric \code{bws}
    vector. This is where bandwidth-selection controls such as
    \code{bwmethod}, \code{bwtype}, and \code{bwscaling},
    kernel/support controls such as \code{ckertype}, \code{ckerorder},
    and \code{ckerbound}, categorical kernel controls such as
    \code{ukertype} and \code{okertype}, and search controls such as
    \code{nmulti} and \code{scale.factor.search.lower} are supplied.
    In \code{predict.npcopula}, additional arguments are passed to
    \code{npcopula} for evaluation with the stored bandwidth object and
    training data. In \code{plot.npcopula}, additional arguments are
    passed to the selected graphics routine, such as
    \code{\link{persp}}, \code{\link{contour}}, \code{\link{image}}, or
    the shared \pkg{rgl} renderer.
  }
}

\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels,
  \code{\link{np.options}} for global options, and \code{\link{plot}}
  for plotting options.

  \code{npcopula} computes the nonparametric copula distribution or
  copula density using marginal quasi-inversion. For the distribution
  target, Sklar's theorem gives
  \deqn{
    C(u_1,\ldots,u_d)
      =
    H(F_1^{-1}(u_1),\ldots,F_d^{-1}(u_d)),
  }
  where \eqn{H} is the joint distribution and \eqn{F_j^{-1}} is the
  quasi-inverse of marginal distribution \eqn{F_j}. For the density
  target, the estimated copula density is
  \deqn{
    c(u)
      =
    \frac{f(x_u)}
         {\prod_{j=1}^d f_j(x_{u,j})},
    \quad x_{u,j}=F_j^{-1}(u_j),
  }
  with numerator and marginal denominators estimated using the selected
  mixed-data kernel bandwidths.

  If \code{u} is provided, \code{\link{expand.grid}} is called on
  \code{u}. As the dimension increases this can become unwieldy because
  a grid with \eqn{m} points in each of \eqn{d} margins has
  \eqn{m^d} rows. Therefore the formula route automatically generates a
  probability grid only for two-dimensional copulas. For
  higher-dimensional copulas, supply \code{u} explicitly or use
  \code{evaluation="sample"}.

  The \sQuote{quasi-inverse} is computed via Definition 2.3.6 from
  Nelsen (2006). An equi-quantile grid on the data range of length
  \code{n.quasi.inv/2} is combined with an equi-spaced grid on the data
  range extended by \code{er.quasi.inv}; the sorted union forms the
  grid used for marginal inversion. If requested probability values lie
  outside the attainable estimated marginal distribution range, they
  are reset to the nearest attainable endpoint. Inspect the returned
  \code{u} columns when endpoint behavior matters.

  The \code{plot.npcopula} method supports base \code{\link{persp}},
  \code{\link{contour}}, and \code{\link{image}} displays for
  two-dimensional grid output. Surface plots use the package-wide
  viridis default palette, detailed perspective ticks, and the same
  default viewing angles and base-graphics rotation cadence as the
  other surface plot methods. \code{renderer="rgl"} requests the shared
  interactive \pkg{rgl} surface renderer, using the same default-angle
  remapping used by the other package surface plots. For mixed ordered
  margins, grid displays are drawn against the requested probability
  grid, while the returned \code{u} columns retain the attainable
  marginal probability values produced by quasi-inversion. The
  \code{"empirical"} view plots empirical copula coordinates, and
  \code{"all"} gives a base-graphics four-panel diagnostic display.

  For grid surfaces, \code{plot.npcopula} can add asymptotic or
  bootstrap intervals. Distribution-copula asymptotic intervals use the
  joint distribution standard error evaluated at the marginal
  quasi-inverse grid. Density-copula asymptotic intervals use the
  plug-in delta-method denominator correction corresponding to
  \eqn{c(u)=f(x_u)/\prod_j f_j(x_{u,j})}. Bootstrap intervals resample
  rows and recompute the plotted copula surface on the same probability
  grid; \code{band="all"} overlays transparent pointwise,
  simultaneous, and Bonferroni wireframes.
}

\value{
  \code{npcopula} returns an object of class \code{"npcopula"}. The main
  components are:

  \item{copula}{
    estimated copula distribution value or copula density value.
  }

  \item{u1, u2, ...}{
    marginal probability coordinates associated with the sample
    realizations or evaluation grid.
  }

  \item{x, y, ...}{
    marginal quasi-inverse coordinates corresponding to the requested
    probability grid when grid evaluation is used.
  }

  \item{bws}{
    selected unconditional distribution or density bandwidth object.
  }

  \item{eval}{
    data frame containing the copula values, probability coordinates, and
    quasi-inverse coordinates. \code{as.data.frame(object)} returns this
    component for data-frame workflows.
  }

  \item{copulaerr}{
    asymptotic or bootstrap standard-error slot. The fitted object stores
    \code{NA} unless an interval-producing plotting route constructs
    evaluation-specific intervals.
  }

  The source data, target, evaluation route, grid dimensions, and timing
  metadata are retained as list components. The functions
  \code{\link{fitted}}, \code{\link{predict}}, \code{\link{se}},
  \code{\link{summary}}, \code{\link{as.data.frame}}, and
  \code{\link{plot}} support \code{"npcopula"} objects.
}

\section{Book And Method Pointers}{
  The copula distribution target is
  \eqn{C(u)=H(F_1^{-1}(u_1),\ldots,F_d^{-1}(u_d))}; the copula density
  target is \eqn{c(u)=f(x_u)/\prod_j f_j(x_{u,j})}. The mixed-data
  kernel implementation follows Racine (2015), with quasi-inversion in
  the sense of Nelsen (2006). For the underlying mixed-data density and
  distribution estimators, see Li and Racine (2007), Chapter 1
  \emph{Density Estimation}, Chapter 3 \emph{Kernel Estimation with
  Mixed Data}, and Racine (2019), Chapter 2 \emph{Continuous Density
  and Cumulative Distribution Functions}.
}

\references{
  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
  Theory and Practice}, Princeton University Press.

  Nelsen, R.B. (2006), \emph{An Introduction to Copulas}, Second
  Edition, Springer.

  Racine, J.S. (2015), \dQuote{Mixed Data Kernel Copulas}, \emph{Empirical
  Economics}, 48, 37--59.

  Racine, J.S. (2019), \emph{An Introduction to the Advanced Theory and
  Practice of Nonparametric Econometrics: A Replicable Approach Using R}.
}

\author{
  Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  Use a \code{\link{data.frame}} rather than \code{\link{cbind}} for
  mixed data so that ordered variables remain ordered. Unordered
  factors are not valid for copula estimation in this implementation.
}

\seealso{
  \code{\link{npudistbw}}, \code{\link{npudist}},
  \code{\link{npudensbw}}, \code{\link{npudens}},
  \code{\link{np.kernels}}, \code{\link{np.options}},
  \code{\link{plot}}
}

\examples{
\dontrun{
library("MASS")

## Example 1: bivariate mixed data, continuous x and ordered y.

set.seed(42)
n <- 1000
n.eval <- 30
rho <- 0.99
mu <- c(0, 0)
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
xy <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
mydat <- data.frame(
  x = xy[, 1],
  y = ordered(as.integer(cut(xy[, 2],
    quantile(xy[, 2], seq(0, 1, by = .1)),
    include.lowest = TRUE)) - 1)
)

grid.seq <- seq(0, 1, length.out = n.eval)
grid.dat <- data.frame(u1 = grid.seq, u2 = grid.seq)

## Estimate the copula distribution from an npudistbw() object.

bw.cdf <- npudistbw(~ x + y, data = mydat, nmulti = 1)
copula <- npcopula(bws = bw.cdf, data = mydat, u = grid.dat)
summary(copula)

## Native plotting replaces the older manual contour(), persp(), and
## empirical scatterplot calls.

plot(copula, perspective = FALSE, view = "contour")
plot(copula, perspective = FALSE, view = "image")
plot(copula, view = "fixed", zlim = c(0, 1))
if (requireNamespace("rgl", quietly = TRUE))
  plot(copula, view = "fixed", renderer = "rgl", zlim = c(0, 1))
plot(copula)

## Plot empirical copula coordinates from the retained sample data.

plot(copula, view = "empirical")

## Or request the four-panel base-graphics diagnostic display.

plot(copula, view = "all")

## Estimate and plot the copula density from an npudensbw() object.

bw.pdf <- npudensbw(~ x + y, data = mydat, nmulti = 1)
copula.dens <- npcopula(bws = bw.pdf, data = mydat, u = grid.dat)
summary(copula.dens)
plot(copula.dens, view = "fixed")
if (requireNamespace("rgl", quietly = TRUE))
  plot(copula.dens, view = "fixed", renderer = "rgl")
plot(copula.dens)

## Intervals are available for two-dimensional grid surfaces.

plot(copula, errors = "asymptotic", band = "pointwise")
plot(copula, errors = "bootstrap", bootstrap = "inid", B = 399,
     band = "pointwise")

## Prediction evaluates the retained bandwidth object on a supplied
## probability grid.

predict(copula, u = data.frame(x = c(0.25, 0.75),
                               y = c(0.25, 0.75)))
predict(copula, newdata = data.frame(u1 = c(0.25, 0.75),
                                     u2 = c(0.25, 0.75)))

## Example 2: bivariate continuous data.

set.seed(42)
n <- 1000
n.eval <- 30
rho <- 0.99
mu <- c(0, 0)
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
xy <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
mydat <- data.frame(x = xy[, 1], y = xy[, 2])

grid.seq <- seq(0, 1, length.out = n.eval)
grid.dat <- data.frame(u1 = grid.seq, u2 = grid.seq)

bw.cdf <- npudistbw(~ x + y, data = mydat, nmulti = 1)
copula <- npcopula(bws = bw.cdf, data = mydat, u = grid.dat)
summary(copula)

plot(copula, perspective = FALSE, view = "contour")
plot(copula, perspective = FALSE, view = "image")
plot(copula, view = "fixed", zlim = c(0, 1))
if (requireNamespace("rgl", quietly = TRUE))
  plot(copula, view = "fixed", renderer = "rgl", zlim = c(0, 1))
plot(copula)
plot(copula, view = "empirical")
plot(copula, view = "all")

bw.pdf <- npudensbw(~ x + y, data = mydat, nmulti = 1)
copula.dens <- npcopula(bws = bw.pdf, data = mydat, u = grid.dat)
summary(copula.dens)
plot(copula.dens, view = "fixed", zlim = c(0, 40))
if (requireNamespace("rgl", quietly = TRUE))
  plot(copula.dens, view = "fixed", renderer = "rgl",
       zlim = c(0, 40))
plot(copula.dens, zlim = c(0, 40))

## The formula interface is a shorter route when bandwidths do not need
## to be reused explicitly.

copula.short <- npcopula(~ x + y, data = mydat, neval = n.eval,
                         nmulti = 1)
plot(copula.short, view = "all")
} % enddontrun
}

\keyword{nonparametric}
