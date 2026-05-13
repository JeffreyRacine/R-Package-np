\name{npcopula}
\alias{npcopula}
\alias{npcopula.default}
\alias{npcopula.formula}
\alias{plot.npcopula}
\alias{print.npcopula}
\alias{summary.npcopula}
\alias{fitted.npcopula}
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
        target = c("distribution", "density"),
        evaluation = c("grid", "sample"),
        u = NULL,
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

\method{plot}{npcopula}(x,
        perspective = TRUE,
        view = c("surface", "contour", "image"),
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
        theta = 45,
        phi = 30,
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

  \item{\dots}{
    additional arguments. In \code{npcopula.formula}, these are passed
    to \code{\link{npudistbw}} or \code{\link{npudensbw}}. In
    \code{plot.npcopula}, these are passed to the selected graphics
    routine, such as \code{\link{persp}}, \code{\link{contour}},
    \code{\link{image}}, or the shared \pkg{rgl} renderer.
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
    \code{expand.grid(u)}. If omitted with a formula route and
    \code{evaluation="grid"}, a two-dimensional grid is generated
    automatically.
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
    plot return mode. \code{"plot"} draws the plot, \code{"data"}
    returns the plotted data with interval columns when requested, and
    \code{"plot-data"} draws and returns the plotted data. \code{"both"}
    is accepted as an alias for \code{"plot-data"}.
  }

  \item{perspective}{
    logical value. If \code{TRUE}, draw a surface display; otherwise
    use \code{view} to choose \code{"contour"} or \code{"image"}.
  }

  \item{phi, theta}{
    viewing angles passed to \code{\link{persp}} or to the shared
    \pkg{rgl} surface renderer.
  }

  \item{renderer}{
    plotting renderer for surface displays. \code{"base"} uses
    \code{\link{persp}}. \code{"rgl"} uses the shared interactive
    \pkg{rgl} surface renderer when the suggested package \pkg{rgl} is
    installed.
  }

  \item{view}{
    display type for grid output: surface, contour, or image. The
    \code{"contour"} and \code{"image"} views are used with
    \code{perspective=FALSE}.
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
  viridis default palette. \code{renderer="rgl"} requests the shared
  interactive \pkg{rgl} surface renderer.

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
  \code{npcopula} returns an object of class \code{"npcopula"} and
  \code{\link{data.frame}}. The main columns are:

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

  The selected bandwidth object, source data, target, evaluation route,
  grid dimensions, and timing metadata are retained as attributes. Use
  \code{attr(object, "bws")} to retrieve the selected bandwidth object.
  The functions \code{\link{fitted}}, \code{\link{summary}}, and
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
## Minimal bivariate continuous example using the formula interface.

library("MASS")

set.seed(42)
n <- 500
rho <- 0.95
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
dat <- as.data.frame(mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma))
names(dat) <- c("x", "y")

fit <- npcopula(~ x + y, data = dat, neval = 20, nmulti = 1)
summary(fit)
plot(fit)

## Intervals are available for two-dimensional grid surfaces.

plot(fit, errors = "asymptotic", band = "pointwise")
plot(fit, errors = "bootstrap", bootstrap = "inid", B = 399,
     band = "pointwise")

## Reuse the selected bandwidth object for a custom probability grid.

u_grid <- data.frame(x = seq(0.05, 0.95, length.out = 15),
                     y = seq(0.05, 0.95, length.out = 15))
custom <- npcopula(bws = attr(fit, "bws"), data = dat, u = u_grid)
plot(custom, perspective = FALSE, view = "contour")

## Copula density uses npudensbw internally when target = "density".

fit_dens <- npcopula(~ x + y, data = dat, target = "density",
                     neval = 20, nmulti = 1)
plot(fit_dens)
} % enddontrun
}

\keyword{nonparametric}
