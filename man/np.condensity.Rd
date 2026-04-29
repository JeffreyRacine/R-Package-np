%
\name{npcdens}
\alias{npcdens}
\alias{npcdens.conbandwidth}
\alias{npcdens.default}
\alias{npcdens.formula}
\title{Kernel Conditional Density Estimation with Mixed Data Types}

\description{
  \code{npcdens} computes kernel conditional density estimates on
  \eqn{p+q}-variate evaluation data, given a set of training data (both
  explanatory and dependent) and a bandwidth specification (a
  \code{conbandwidth} object or a bandwidth vector, bandwidth type, and
  kernel type) using the method of Hall, Racine, and Li (2004).
  The data may be continuous, discrete (unordered and ordered
  factors), or some combination thereof.
}

\usage{
npcdens(bws, \dots)

\method{npcdens}{formula}(bws, data = NULL, newdata = NULL, \dots)


\method{npcdens}{conbandwidth}(bws,
        txdat = stop("invoked without training data 'txdat'"),
        tydat = stop("invoked without training data 'tydat'"),
        exdat,
        eydat,
        gradients = FALSE,
        proper = FALSE,
        proper.method = c("project"),
        proper.control = list(),
        \dots)

\method{npcdens}{default}(bws, txdat, tydat, nomad = FALSE, \dots)

}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface, and training data.
}

  \item{bws}{
    a bandwidth specification. This can be set as a \code{conbandwidth}
    object returned from a previous invocation of
    \code{\link{npcdensbw}}, or as a \eqn{p+q}-vector of bandwidths,
    with each element \eqn{i} up to \eqn{i=q} corresponding to the
    bandwidth for column \eqn{i} in \code{tydat}, and each element
    \eqn{i} from \eqn{i=q+1} to \eqn{i=p+q} corresponding to the
    bandwidth for column \eqn{i-q} in \code{txdat}. If specified as a
    vector, then additional arguments will need to be supplied as
    necessary to specify the bandwidth type, kernel types, training
    data, and so on.
  }


  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
   \code{\link{npcdensbw}} was called.
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


\subsection{Bandwidth Search Shortcut}{
This argument passes the recommended automatic local-polynomial NOMAD preset to \code{\link{npcdensbw}} when bandwidths are computed inside \code{npcdens}.
}

  \item{nomad}{
    logical shortcut passed through to \code{\link{npcdensbw}} when
    bandwidths are computed inside \code{npcdens}. When \code{TRUE},
    the bandwidth route fills any missing values among \code{regtype},
    \code{search.engine}, \code{degree.select},
    \code{bernstein.basis}, \code{degree.min}, \code{degree.max},
    \code{degree.verify}, and \code{bwtype} with the recommended
    automatic LP NOMAD preset documented in \code{\link{npcdensbw}}.
    Additional NOMAD tuning arguments such as \code{nomad.nmulti} may
    also be supplied through \code{\dots}; \code{nmulti} remains the
    outer restart count while \code{nomad.nmulti} controls inner
    \code{crs::snomadr()} multistarts within each outer restart.
    After fitting, inspect \code{fit$bws$nomad.shortcut} on the
    returned object \code{fit} to see the normalized shortcut metadata.
  }


\subsection{Evaluation Data And Returned Quantities}{
These arguments control where the fitted conditional density is evaluated and which estimates are returned.
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


  \item{gradients}{
    a logical value specifying whether to return estimates of the
    gradients at the evaluation points. Defaults to \code{FALSE}.
  }


  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.
  }



\subsection{Fit Properization Controls}{
These arguments control optional post-estimation properization of the fitted conditional density.
}

  \item{proper}{
    a logical value specifying whether to apply post-estimation
    properization to the conditional density estimate. Defaults to
    \code{FALSE}.
  }


  \item{proper.control}{
    a list of controls for properization. Supported entries are
    \code{tol}, \code{grid.check}, \code{store.raw}, and
    \code{fail.on.unsupported}.
  }



  \item{proper.method}{
    a character string specifying the properization method. Currently
    \code{"project"} is supported.
  }


\subsection{Additional Arguments}{
Further arguments are passed to \code{\link{npcdensbw}} when bandwidths are computed internally, or used to interpret a numeric \code{bws} vector.
}

  \item{\dots}{
    additional arguments supplied to \code{\link{npcdensbw}} when
    \code{npcdens} computes bandwidths internally, or arguments needed
    to interpret a numeric \code{bws} vector. This is where bandwidth
    selection controls such as \code{bwmethod}, \code{bwtype},
    kernel/support controls such as \code{cxkertype},
    \code{cykertype}, \code{cxkerbound}, and \code{cykerbound},
    search controls such as \code{nmulti},
    \code{scale.factor.search.lower}, and \code{nomad.nmulti}, and
    local-polynomial controls such as \code{regtype}, \code{degree},
    \code{basis}, and \code{bernstein.basis} are supplied. See
    \code{\link{npcdensbw}} for the complete bandwidth-selection
    argument surface.

  }

}
\details{
  Documentation guide: see \code{\link{npcdensbw}} for bandwidth selection and search controls, \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, and \code{\link{plot}}, \code{\link{plot.np}} for plotting options.

  When \code{bws} is omitted, the formula and default methods call
  \code{\link{npcdensbw}} first and pass bandwidth-selection arguments
  from \code{\dots} to that call. When \code{bws} is already a
  \code{conbandwidth} object, \code{npcdens} estimates with the stored
  bandwidth metadata in that object.

  Argument groups for bandwidth selection are documented on
  \code{\link{npcdensbw}}. The most common workflow is to choose data
  and bandwidth inputs first, then bandwidth criterion and
  representation, then kernel/support controls, numerical search
  controls, bounded \code{cv.ls} quadrature controls if relevant, and
  finally local-polynomial/NOMAD controls for polynomial-adaptive fits.

  For S3 plotting help, use \code{methods("plot")} and query
  class-specific help topics such as \code{?plot.npregression} and
  \code{?plot.rbandwidth}. You can inspect implementations with
  \code{getS3method("plot","npregression")}.

  \code{npcdens} implements a variety of methods for estimating
  multivariate conditional distributions (\eqn{p+q}-variate) defined
  over a set of possibly continuous and/or discrete (unordered, ordered)
  data. The approach is based on Li and Racine (2004) who employ
  \sQuote{generalized product kernels} that admit a mix of continuous
  and discrete data types.

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

  For practitioners who want the recommended automatic LP NOMAD route
  without spelling out all LP tuning arguments,
  \code{npcdens(..., nomad=TRUE)} and \code{npcdensbw(..., nomad=TRUE)}
  expand missing settings to the same documented preset. Explicit
  incompatible settings fail fast rather than being silently rewritten.

}
\value{

  \code{npcdens} returns a \code{condensity} object. The generic
  accessor functions \code{\link{fitted}}, \code{\link{se}}, and
  \code{\link{gradients}}, extract estimated values, asymptotic standard
  errors on estimates, and gradients, respectively, from the returned
  object. Furthermore, the functions \code{\link{predict}},
  \code{\link{summary}} and \code{\link{plot}} support objects of both
  classes. The returned objects have the following components:

  \item{xbw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    explanatory data, \code{txdat} }
  \item{ybw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    dependent data, \code{tydat} }
  \item{xeval}{ the evaluation points of the explanatory data }
  \item{yeval}{ the evaluation points of the dependent data }
  \item{condens}{ estimates of the conditional density at the evaluation
    points }
  \item{conderr}{ standard errors of the conditional density
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

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}, \code{\link{plot.np}}
  \code{\link{npudens}}
}
\examples{
\dontrun{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we load Giovanni
# Baiocchi's Italian GDP panel (see Italy for details), and compute the
# likelihood cross-validated bandwidths (default) using a second-order
# Gaussian kernel (default). Note - this may take a minute or two
# depending on the speed of your computer.

data("Italy")
Italy <- Italy[seq_len(min(300, nrow(Italy))), ]
attach(Italy)

# First, compute the bandwidths... note that this may take a minute or
# two depending on the speed of your computer.

bw <- npcdensbw(formula=gdp~ordered(year), nmulti=1)

# Next, compute the condensity object...

fhat <- npcdens(bws=bw)

# The object fhat now contains results such as the estimated conditional
# density function (fhat$condens) and so on...

summary(fhat)

# Call the plot() function to visualize the results (<ctrl>-C will
# interrupt on *NIX systems, <esc> will interrupt on MS Windows
# systems).

if (interactive()) plot(bw)

detach(Italy)

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we load Giovanni
# Baiocchi's Italian GDP panel (see Italy for details), and compute the
# likelihood cross-validated bandwidths (default) using a second-order
# Gaussian kernel (default). Note - this may take a minute or two
# depending on the speed of your computer.

data("Italy")
Italy <- Italy[seq_len(min(300, nrow(Italy))), ]
attach(Italy)

# First, compute the bandwidths... note that this may take a minute or
# two depending on the speed of your computer.

# Note - we cast `X' and `y' as data frames so that plot() can
# automatically grab names (this looks like overkill, but in
# multivariate settings you would do this anyway, so may as well get in
# the habit).

X <- data.frame(year=ordered(year))
y <- data.frame(gdp)

bw <- npcdensbw(xdat=X, ydat=y, nmulti=1)

# Next, compute the condensity object...

fhat <- npcdens(bws=bw)

# The object fhat now contains results such as the estimated conditional
# density function (fhat$condens) and so on...

summary(fhat)

# Call the plot() function to visualize the results (<ctrl>-C will
# interrupt on *NIX systems, <esc> will interrupt on MS Windows systems).

if (interactive()) plot(bw)

detach(Italy)

# EXAMPLE 2 (INTERFACE=FORMULA): For this example, we load the old
# faithful geyser data from the R `datasets' library and compute the
# conditional density function.

library("datasets")
data("faithful")
attach(faithful)

# Note - this may take a few minutes depending on the speed of your
# computer...

bw <- npcdensbw(formula=eruptions~waiting, nmulti=1)

summary(bw)

# Plot the density function (<ctrl>-C will interrupt on *NIX systems,
# <esc> will interrupt on MS Windows systems).

if (interactive()) plot(bw)

detach(faithful)

# EXAMPLE 2 (INTERFACE=DATA FRAME): For this example, we load the old
# faithful geyser data from the R `datasets' library and compute the
# conditional density function.

library("datasets")
data("faithful")
attach(faithful)

# Note - this may take a few minutes depending on the speed of your
# computer...

# Note - we cast `X' and `y' as data frames so that plot() can
# automatically grab names (this looks like overkill, but in
# multivariate settings you would do this anyway, so may as well get in
# the habit).

X <- data.frame(waiting)
y <- data.frame(eruptions)

bw <- npcdensbw(xdat=X, ydat=y, nmulti=1)

summary(bw)

# Plot the density function (<ctrl>-C will interrupt on *NIX systems,
# <esc> will interrupt on MS Windows systems)

if (interactive()) plot(bw)

detach(faithful)

# EXAMPLE 3 (INTERFACE=FORMULA): Replicate the DGP of Klein & Spady
# (1993) (see their description on page 405, pay careful attention to
# footnote 6 on page 405).

set.seed(123)

n <- 300

# x1 is chi-squared having 3 df truncated at 6 standardized by
# subtracting 2.348 and dividing by 1.511

x <- rchisq(n, df=3)
x1 <- (ifelse(x < 6, x, 6) - 2.348)/1.511

# x2 is normal (0, 1) truncated at +- 2 divided by 0.8796

x <- rnorm(n)
x2 <- ifelse(abs(x) < 2 , x, 2) / 0.8796

# y is 1 if y* > 0, 0 otherwise.

y <- ifelse(x1 + x2 + rnorm(n) > 0, 1, 0)

# Generate data-driven bandwidths (likelihood cross-validation). Note -
# this may take a few minutes depending on the speed of your computer...

bw <- npcdensbw(formula=factor(y)~x1+x2, nmulti=1)

# Next, create the evaluation data in order to generate a perspective
# plot

x1.seq <- seq(min(x1), max(x1), length=30)
x2.seq <- seq(min(x2), max(x2), length=30)
X.eval <- expand.grid(x1=x1.seq,x2=x2.seq)

data.eval <- data.frame(y=factor(rep(1, nrow(X.eval))),x1=X.eval[,1],x2=X.eval[,2])

# Now evaluate the conditional probability for y=1 and for the
# evaluation Xs

fit <- fitted(npcdens(bws=bw,newdata=data.eval))

# Finally, coerce the data into a matrix for plotting with persp()

fit.mat <- matrix(fit, 30, 30)

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

n <- 300

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

# Generate data-driven bandwidths (likelihood cross-validation). Note -
# this may take a few minutes depending on the speed of your computer...

bw <- npcdensbw(xdat=X, ydat=factor(y), nmulti=1)

# Next, create the evaluation data in order to generate a perspective
# plot

x1.seq <- seq(min(x1), max(x1), length=30)
x2.seq <- seq(min(x2), max(x2), length=30)
X.eval <- expand.grid(x1=x1.seq,x2=x2.seq)

# Now evaluate the conditional probability for y=1 and for the
# evaluation Xs

fit <- fitted(npcdens(exdat=X.eval,
               eydat=factor(rep(1, nrow(X.eval))),
               bws=bw))

# Finally, coerce the data into a matrix for plotting with persp()

fit.mat <- matrix(fit, 30, 30)

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

# EXAMPLE 4: Variations on local polynomial conditional density
# estimation with proper = TRUE.

data("Italy")

Italy2 <- within(Italy, {
  year <- as.numeric(as.character(year))
})

# Plot only: make the plotted surface proper on the plot evaluation grid.

fhat <- npcdens(gdp ~ year, data = Italy2,
                regtype = "lp", degree = 3, nmulti = 1)

plot(fhat, proper = TRUE)

# Fit an object whose fitted values are themselves proper.

ctrl_fit <- list(
  mode = "slice",
  apply = "fitted",
  slice.grid.size = 101L,
  slice.extend.factor = 0.1
)

fhat_fit <- npcdens(
  gdp ~ year,
  data = Italy2,
  regtype = "lp",
  degree = 3,
  nmulti = 1,
  proper = TRUE,
  proper.control = ctrl_fit
)

fit_proper <- fitted(fhat_fit)
fit_raw <- fhat_fit$condens.raw

# Display the repaired and raw fitted values for cases where the raw
# fitted density is negative.

head(cbind(fit_proper, fit_raw)[which(fit_raw < 0), ])

# Predict on a common explicit y-grid for several years, and render
# those predictions proper.

g.grid <- seq(min(Italy2$gdp), max(Italy2$gdp), length.out = 200)

nd_grid <- expand.grid(
  gdp = g.grid,
  year = c(1955, 1975, 1995)
)

pred_grid <- predict(fhat, newdata = nd_grid, proper = TRUE)

# Predict on paired rows with different gdp grids by year, and still
# make the predictions proper via slice mode.

g1 <- seq(quantile(Italy2$gdp, 0.10),
          quantile(Italy2$gdp, 0.60), length.out = 60)
g2 <- seq(quantile(Italy2$gdp, 0.30),
          quantile(Italy2$gdp, 0.90), length.out = 35)

nd_slice <- rbind(
  data.frame(gdp = g1, year = rep(1960, length(g1))),
  data.frame(gdp = g2, year = rep(1985, length(g2)))
)

pred_slice <- predict(
  fhat,
  newdata = nd_slice,
  proper = TRUE,
  proper.control = list(mode = "slice")
)

# One object that carries properization for fitted values and for later
# predict() calls.

ctrl_both <- list(
  mode = "slice",
  apply = "both",
  slice.grid.size = 101L,
  slice.extend.factor = 0.1
)

fhat_both <- npcdens(
  gdp ~ year,
  data = Italy2,
  regtype = "lp",
  degree = 3,
  nmulti = 1,
  proper = TRUE,
  proper.control = ctrl_both
)

fit_both <- fitted(fhat_both)
pred_both <- predict(
  fhat_both,
  newdata = nd_slice,
  proper.control = ctrl_both
)
} % enddontrun
}
\keyword{nonparametric}
