%
\name{npreg}
\alias{npreg}
\alias{npreg.default}
\alias{npreg.formula}
\alias{npreg.rbandwidth}

\title{Kernel Regression with Mixed Data Types}

\description{
  \code{npreg} computes a kernel regression estimate of a one
  (1) dimensional dependent variable on \eqn{p}-variate explanatory
  data, given a set of evaluation points, training points (consisting of
  explanatory data and dependent data), and a bandwidth specification
  using the method of Racine and Li (2004) and Li and Racine (2004). A
  bandwidth specification can be a \code{rbandwidth} object, or a
  bandwidth vector, bandwidth type and kernel type.
}

\usage{
npreg(bws, \dots)

\method{npreg}{formula}(bws, 
      data = NULL, 
      newdata = NULL, 
      y.eval = FALSE, 
      \dots)

\method{npreg}{default}(bws, 
      txdat, 
      tydat, 
      nomad = FALSE, 
      \dots)

\method{npreg}{rbandwidth}(bws,
      txdat = stop("training data 'txdat' missing"),
      tydat = stop("training data 'tydat' missing"),
      exdat,
      eydat,
      gradient.order = 1L,
      gradients = FALSE,
      residuals = FALSE,
      \dots)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface, and training data.
}

  \item{bws}{
    a bandwidth specification. This can be set as a \code{rbandwidth}
    object returned from an invocation of \code{\link{npregbw}}, or
    as a vector of bandwidths, with each element \eqn{i} corresponding
    to the bandwidth for column \eqn{i} in \code{txdat}. If specified as
    a vector, then additional arguments will need to be supplied as
    necessary to specify the bandwidth type, kernel types, and so on.
  }


  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npregbw}} was called.
  }


  \item{txdat}{
    a \eqn{p}-variate data frame of explanatory data (training data) used to
    calculate the regression estimators. Defaults to the corresponding
    training data recoverable from the bandwidth object in the current
    session.
  }


  \item{tydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{txdat}. Defaults to the corresponding
    training data recoverable from the bandwidth object in the current
    session.
  }


\subsection{Local-Polynomial Degree And Bandwidth Search}{
This argument controls the recommended automatic local-polynomial NOMAD route, which jointly selects continuous polynomial degree and bandwidths when these are computed inside \code{npreg}.
}

  \item{nomad}{
    logical shortcut passed through to \code{\link{npregbw}} when
    bandwidths are computed inside \code{npreg}. When \code{TRUE}, the
    regression bandwidth route fills any missing values among
    \code{regtype}, \code{search.engine}, \code{degree.select},
    \code{bernstein.basis}, \code{degree.min}, \code{degree.max},
    \code{degree.verify}, and \code{bwtype} with the recommended
    automatic local-polynomial degree-and-bandwidth NOMAD preset documented in \code{\link{npregbw}}.
    Additional NOMAD tuning arguments such as \code{nomad.nmulti} may
    also be supplied through \code{\dots}; \code{nmulti} remains the
    outer restart count while \code{nomad.nmulti} controls inner
    native \pkg{crs} NOMAD API multistarts within each outer restart.
    After fitting, inspect \code{fit$bws$nomad.shortcut} on the
    returned object \code{fit} to see the normalized shortcut metadata.
  }


\subsection{Evaluation Data And Returned Quantities}{
These arguments control where the regression is evaluated and which fitted quantities are returned.
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


  \item{gradient.order}{
    for \code{regtype="lp"} with \code{gradients=TRUE}, a positive
    integer (or integer vector with one entry per continuous predictor)
    specifying derivative order(s). Defaults to \code{1L}. Orders
    exceeding \code{degree} for a variable are returned as
    \code{NA}, apart from the established all-degree-zero first-derivative
    compatibility route. Higher-order derivatives are available for continuous
    predictors when the requested order does not exceed the
    corresponding local-polynomial degree.
  }


  \item{gradients}{
    a logical value indicating that you want gradients computed and
    returned in the resulting \code{npregression} object. Defaults to
    \code{FALSE}. For \code{regtype="lp"}, derivative components that
    are not defined (requested order exceeds the variable-specific
    polynomial degree) are returned as \code{NA}. For ordered and
    unordered predictors the gradient columns are first-order
    effects/contrasts where supported; \code{gradient.order} controls
    derivative order for continuous predictors.
    Plotting retains one panel per requested predictor: an unavailable
    continuous derivative is shown as an empty, correctly labelled panel when
    another continuous derivative or categorical effect is available. A
    request for which no component is available is rejected with the affected
    orders and degrees reported. Gradients are not currently available with
    \code{ckertype="beta"}.
  }


  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.
  }


  \item{residuals}{
    a logical value indicating that you want residuals computed and
    returned in the resulting \code{npregression} object. Defaults to
    \code{FALSE}.
  }


  \item{y.eval}{
    If \code{newdata} contains dependent data and \code{y.eval = TRUE},
    \code{\link{npRmpi}} will compute goodness of fit statistics on these
    data and return them. Defaults to \code{FALSE}.
  }


\subsection{Additional Arguments}{
Further arguments are passed to \code{\link{npregbw}} when bandwidths are computed internally, or used to interpret a numeric \code{bws} vector.
}

  \item{\dots}{
    additional arguments supplied to \code{\link{npregbw}} when
    \code{npreg} computes bandwidths internally, or arguments needed to
    interpret a numeric \code{bws} vector. This is where bandwidth
    selection controls such as \code{bwmethod}, \code{bwtype},
    \code{bwscaling}, kernel/support controls such as \code{ckertype},
    \code{ckerorder}, and \code{ckerbound}, categorical kernel controls
    such as \code{ukertype} and \code{okertype}, search controls such
    as \code{nmulti} and \code{scale.factor.search.lower}, and
    local-polynomial/NOMAD controls such as \code{regtype},
    \code{degree}, \code{basis}, \code{bernstein.basis},
    \code{degree.select}, and \code{nomad.nmulti} are supplied. See
    \code{\link{npregbw}} for the complete bandwidth-selection argument
    surface.
  }
}
\details{
  Documentation guide: see \code{\link{npregbw}} for bandwidth
  selection and search controls, \code{\link{np.kernels}} for kernels,
  \code{\link{np.options}} for global options, \code{\link{plot}},
  \code{\link{plot.np}} for plotting options, and
  \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See
  \code{\link{npRmpi.init}} details for performance tradeoffs (message
  passing/startup mode) and the \code{inst/Rprofile} manual-broadcast
  template.

  When \code{bws} is omitted, the formula and default methods call
  \code{\link{npregbw}} first and pass bandwidth-selection arguments
  from \code{\dots} to that call. When \code{bws} is already an
  \code{rbandwidth} object, \code{npreg} estimates with the stored
  bandwidth metadata in that object.

  Argument groups for bandwidth selection are documented on
  \code{\link{npregbw}}. The most common workflow is to choose data and
  bandwidth inputs first, then bandwidth criterion and representation,
  then kernel/support controls, and finally local-polynomial/NOMAD
  controls when using polynomial-adaptive fits.

  For S3 plotting help, see \code{\link{plot.np}}. You can list
  available plot methods with \code{methods("plot")}.


  Typical usages are (see below for a complete list of options and also
  the examples at the end of this help file)
  \preformatted{
    Usage 1: first compute the bandwidth object via npregbw and then
    compute the conditional mean:

    bw <- npregbw(y~x)
    ghat <- npreg(bw)

    Usage 2: alternatively, compute the bandwidth object indirectly:

    ghat <- npreg(y~x)

    Usage 3: modify the default kernel and order:

    ghat <- npreg(y~x, ckertype="epanechnikov", ckerorder=4)

    Usage 4: use the data frame interface rather than the formula
    interface:

    ghat <- npreg(tydat=y, txdat=x, ckertype="epanechnikov", ckerorder=4)
  }



  \code{npreg} implements a variety of methods for regression on
  multivariate (\eqn{p}-variate) data, the types of which are possibly
  continuous and/or discrete (unordered, ordered). The approach is
  based on Li and Racine (2003) who employ \sQuote{generalized product kernels}
  that admit a mix of continuous and discrete data types.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
  density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
  with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
  are constant over the support of \eqn{x}.

  Data contained in the data frame \code{txdat} may be a mix of
  continuous (default), unordered discrete (to be specified in the data
  frame \code{txdat} using \code{\link{factor}}), and ordered discrete
  (to be specified in the data frame \code{txdat} using
  \code{\link{ordered}}). Data can be entered in an arbitrary order and
  data types will be detected automatically by the routine (see
  \code{\link{npRmpi}} for details).

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. A second-order beta associated kernel is also available for the
  manual, fixed-bandwidth local-constant route described below. Unordered
  discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

  For \code{ckertype="beta"}, suppose continuous predictor \eqn{j} has
  declared finite support \eqn{[a_j,b_j]}, width
  \eqn{L_j=b_j-a_j}, public metric bandwidth \eqn{h_j}, evaluation coordinate
  \eqn{x_j}, and observation \eqn{X_{ij}}. Define
  \deqn{t_j=(x_j-a_j)/L_j,\quad z_{ij}=(X_{ij}-a_j)/L_j,\quad
        \tau_j=(L_j/h_j)^2.}
  The local-constant product weight is
  \deqn{K_i(x)=\prod_j \left\{\frac{1}{L_j}
    \frac{z_{ij}^{t_j\tau_j}(1-z_{ij})^{(1-t_j)\tau_j}}
    {B(1+t_j\tau_j,1+(1-t_j)\tau_j)}\right\},}
  and the fitted conditional mean is the Nadaraya--Watson ratio
  \deqn{\widehat m(x)=\frac{\sum_i K_i(x)Y_i}{\sum_i K_i(x)}.}
  This uses Chen's evaluation-centred beta associated-kernel contribution
  inside the established \pkg{np} local-constant ratio. It should not be
  confused with the Gasser--Muller beta regression smoother studied by Chen
  (2000). The support scale is already present in each contribution, and the
  implementation applies neither the legacy \eqn{1/h_j} factor nor the generic
  bounded-kernel normalizer. Although that normalizer cancels from ordinary
  local-constant numerator/denominator ratios, the beta kernel is not merely
  such a correction: its shape parameters change with the evaluation point.

  Products and ratios are evaluated in the log domain, so representable beta
  regression results are retained even when every unscaled product would
  underflow to zero. If \eqn{\ell_i=K_i/\sum_r K_r}, the reported beta-route
  standard error is
  \deqn{\left\{\left[\sum_i \ell_i(Y_i-\widehat m)^2\right]
                    \left[\sum_i \ell_i^2\right]\right\}^{1/2},}
  preserving the established uncorrected local-variance convention while
  using the exact effective-sample factor for these associated weights.

  The present beta regression route requires \code{regtype="lc"},
  \code{ckerorder=2}, \code{bwtype="fixed"},
  \code{ckerbound="fixed"}, explicit finite \code{ckerlb}/\code{ckerub},
  continuous predictors only, manual bandwidths, and
  \code{gradients=FALSE}. Automatic bandwidth selection,
  nearest-neighbour bandwidths, sample-range support, categorical predictors,
  and local-linear/local-polynomial fits fail explicitly. An evaluation point
  with no supported observation weight also produces an explicit error.

  When bandwidths are obtained with \code{regtype="lp"}, C-level
  \code{npreg} supports heterogeneous continuous polynomial degrees via
  \code{degree}. The \code{basis} selector currently supports
  \code{basis="glp"}, \code{"additive"}, and \code{"tensor"}.
  For continuous predictors with degree vector
  \eqn{d}{d}, additive basis size is
  \eqn{1+\sum_j d_j}{1+sum(d[j])}, tensor basis size is
  \eqn{\prod_j (d_j+1)}{prod(d[j]+1)}, and GLP uses admissible
  multi-indices \eqn{\alpha}{alpha} with
  \eqn{\alpha_j \le d_j}{alpha[j] <= d[j]} and
  \eqn{0<\sum_j \alpha_j \le \max_j d_j}{0<sum(alpha)<=max(d)} plus
  an intercept. Thus, for \eqn{d=(2,1)}{d=(2,1)}, GLP contains
  \eqn{1,x_1,x_1^2,x_2,x_1x_2}{1, x1, x1^2, x2, x1*x2}; the additive
  basis omits \eqn{x_1x_2}{x1*x2}, while the tensor basis additionally
  contains \eqn{x_1^2x_2}{x1^2*x2}.

  The optional flag \code{bernstein.basis} controls basis construction:
  \code{FALSE} (default) uses raw local-polynomial powers. For
  \code{basis="glp"}, the historical selector \code{TRUE} uses a
  degree-graded, orthonormal shifted-Legendre representation. On the unit-scaled support
  \eqn{u\in[0,1]}{u in [0,1]}, its order-\eqn{k}{k} factor is
  \eqn{\sqrt{2k+1}P_k(2u-1)}{sqrt(2*k+1) P[k](2*u-1)}, where
  \eqn{P_k}{P[k]} is a Legendre polynomial. This is an exact
  degree-\eqn{k}{k} linear combination of Bernstein polynomials and, because
  the GLP multi-index set is downward closed, spans exactly the same
  heterogeneous polynomial space as the raw representation. It is evaluated
  by a deterministic Legendre recurrence, not by sample-dependent
  Gram--Schmidt or QR orthogonalization. Additive and
  tensor bases retain their fixed-coordinate-degree Bernstein
  representation. The
  homogeneous degree-0 and degree-1 cases remain
  equivalent to \code{lc} and \code{ll}, respectively. For continuous
  predictors, \code{gradients=TRUE} returns the derivative order(s)
  requested by \code{gradient.order}; the associated \code{gerr}
  entries are standard errors for those same derivative orders. Both
  raw and shifted-Legendre representations report fitted values and derivatives
  on the original response/predictor scale. In mixed-data GLP settings,
  unordered and ordered predictors contribute first-order
  effects/contrasts rather than higher-order derivatives.
  When \code{npregbw(..., regtype="lp")} is used with
  \code{degree.select="manual"}, the degree vector remains fixed user
  input. When \code{degree.select != "manual"}, \code{npregbw} can
  jointly select polynomial degree and bandwidth using either the
  cached cell-search backend or the direct
  \code{search.engine="nomad"}/\code{"nomad+powell"} route described in
  \code{\link{npregbw}}; the latter follows Hall and Racine (2015). For
  practitioners who want that recommended route without spelling out all
  LP tuning arguments, \code{npreg(..., nomad=TRUE)} and
  \code{npregbw(..., nomad=TRUE)} expand missing settings to the same
  documented automatic-LP NOMAD preset. Explicit incompatible settings
  fail fast rather than being silently rewritten. The direct NOMAD
  backend is provided by the suggested package \code{crs}, so install
  \code{crs} before using \code{search.engine="nomad"},
  \code{"nomad+powell"}, or \code{nomad=TRUE}. For
  \code{bernstein.basis=TRUE}, evaluation points for continuous predictors
  must lie within training support; use \code{bernstein.basis=FALSE} for
  extrapolation. For \code{regtype="ll"} and \code{regtype="lp"}, the
  training continuous design is checked for rank deficiency and extreme
  condition number before estimation proceeds.

  The use of compactly supported kernels or the occurrence of small
  bandwidths can lead to numerical problems for the local linear
  estimator when computing the locally weighted least squares
  solution. To overcome this problem we rely on a form or
  \sQuote{ridging} proposed by Cheng, Hall, and Titterington (1997),
  modified so that we solve the problem pointwise rather than globally
  (i.e. only when it is needed).

}

\value{

  \code{npreg} returns a \code{npregression} object.
  The generic
  functions \code{\link{fitted}}, \code{\link{residuals}},
  \code{\link{se}}, \code{\link{predict}}, and
  \code{\link{gradients}}, extract (or generate) estimated values,
  residuals, asymptotic standard
  errors on estimates, predictions, and gradients, respectively, from
  the returned object. Furthermore, the functions \code{\link{summary}}
  and \code{\link{plot}} support objects of this type. The returned object
  has the following components:


  \item{eval}{ evaluation points }
  \item{mean}{ estimates of the regression function (conditional mean) at the
  evaluation points }
  \item{merr}{ standard errors of the regression function estimates }
  \item{grad}{ estimates of the gradients at each evaluation point }
  \item{gerr}{ standard errors of the gradient estimates }
  \item{resid}{ if \code{residuals = TRUE}, in-sample or out-of-sample
    residuals where appropriate (or possible) }
  \item{R2}{ coefficient of determination (Doksum and Samarov (1995))}
  \item{MSE}{ mean squared error }
  \item{MAE}{ mean absolute error }
  \item{MAPE}{ mean absolute percentage error }
  \item{CORR}{ absolute value of Pearson's correlation coefficient }
  \item{SIGN}{ fraction of observations where fitted and observed values
  agree in sign }

}

\section{Book And Method Pointers}{
  The regression target is the conditional mean
  \eqn{m(x)=E[Y\mid X=x]}. Local-constant estimation can be written
  schematically as \eqn{\hat m(x)=\sum_i W_i(x)Y_i}, where the weights
  are normalized generalized product-kernel weights over the continuous,
  unordered, and ordered components of \eqn{X}. Local-polynomial
  regression replaces this scalar weighted average by a weighted least
  squares fit in powers \eqn{b(X_i)}{b(X[i])} of the original continuous
  coordinates, with kernel weights centered at \eqn{x}. If
  \eqn{\widehat\beta_x}{beta-hat[x]} denotes that local weighted fit, the
  reported level and derivative are respectively
  \eqn{b(x)^T\widehat\beta_x}{b(x)' beta-hat[x]} and
  \eqn{D^\gamma b(x)^T\widehat\beta_x}{D^gamma b(x)' beta-hat[x]}.
  For the downward-closed GLP space this uncentered representation is
  numerically equivalent, up to floating-point solution paths, to fitting
  powers of \eqn{X_i-x}{X[i]-x} and reading the corresponding intercept and
  derivative coefficients.

  For book-length derivations, see Li and Racine (2007), Chapter 2
  \emph{Regression}, especially Sections 2.1, 2.2, 2.4, and 2.5, and
  Chapter 4 \emph{Kernel Estimation with Mixed Data}, especially Sections
  4.2 and 4.4. The later workflow treatment is Racine (2019), Chapter 6
  \emph{Conditional Mean Function Estimation}, including the local
  polynomial and mixed-data marginal-effect discussion.
}

\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
  discrimination by the kernel method,} Biometrika, 63, 413-420.

  Cheng, M.-Y. and P. Hall and D.M. Titterington (1997), \dQuote{On the
  shrinkage of local linear curve estimators,} Statistics and Computing,
  7, 11-17.

  Chen, S.X. (1999), \dQuote{Beta kernel estimators for density functions,}
  Computational Statistics and Data Analysis, 31, 131-145.

  Chen, S.X. (2000), \dQuote{Beta kernel smoothers for regression curves,}
  Statistica Sinica, 10, 73-91.

  Cohen, A. and G. Migliorati (2018), \dQuote{Multivariate approximation
  in downward closed polynomial spaces,} in R. Cools and D. Nuyens (eds.),
  \emph{Contemporary Computational Mathematics}, Springer, 233-282.

  Dung, D. and T. Ullrich (2011), \dQuote{Whitney type inequalities for
  local anisotropic polynomial approximation,} Journal of Approximation
  Theory, 163, 1590-1605.

  Fan, J. and I. Gijbels (1996), \emph{Local Polynomial Modelling and
  Its Applications,} Chapman and Hall.

  Doksum, K. and A. Samarov (1995), \dQuote{Nonparametric estimation of
  global functionals and a measure of the explanatory power of
  covariates in regression,} The Annals of Statistics, 23 1443-1473.

  Hall, P. and Q. Li and J.S. Racine (2007), \dQuote{Nonparametric
  estimation of regression functions in the presence of irrelevant
  regressors,} The Review of Economics and Statistics, 89, 784-789.

  Hall, P. and J.S. Racine (2015), \dQuote{Infinite Order Cross-Validated
  Local Polynomial Regression,} Journal of Econometrics, 185, 510-525.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
  nonparametric regression,} Statistica Sinica, 14, 485-512.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.

  Racine, J.S. and Q. Li (2004), \dQuote{Nonparametric estimation of
    regression functions with both categorical and continuous data,}
  Journal of Econometrics, 119, 99-130.

  Ruppert, D. and M.P. Wand (1994), \dQuote{Multivariate locally weighted
  least squares regression,} The Annals of Statistics, 22, 1346-1370.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth
    estimators for discrete distributions,} Biometrika, 68, 301-309.

}

\author{
  Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.

  When reusing a bandwidth object without explicit training or evaluation
  data, the required data are recovered from object-retained state or
  from the stored call environment when available. This is a convenience
  for ordinary same-session use, not a guarantee that a saved object is
  self-contained. For reproducible replay after \code{saveRDS()} and
  \code{readRDS()} in a fresh R session, supply the required data
  arguments explicitly.
}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}, \code{\link{plot.np}}
  \code{\link{loess}}
}
\examples{
\dontrun{
## The following example is adapted for interactive parallel execution
## in R. Here we spawn 1 slave so that there will be two compute nodes
## (master and slave).  Kindly see the batch examples in the demos
## directory (npRmpi/demos) and study them carefully. Also kindly see
## the more extensive examples in the np package itself. See the npRmpi
## vignette for further details on running parallel np programs via
## vignette("npRmpi_getting_started", package = "npRmpi").
##
## R CMD check can set up a process environment where spawn-mode MPI
## teardown may terminate the check parent process. Keep this example
## fully runnable for users while skipping MPI spawn only in check.
force.run <- nzchar(Sys.getenv("NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK"))

in.check <- !force.run && (
            nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
            nzchar(Sys.getenv("_R_CHECK_INTERNALS2_")) ||
            nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_")) ||
            nzchar(Sys.getenv("_R_CHECK_EXAMPLE_TIMING_THRESHOLD_"))
          )

if (!in.check) {
  npRmpi.init(nslaves=1)

  set.seed(42)

  n <- 250

  x <- runif(n)
  z1 <- rbinom(n,1,.5)
  z2 <- rbinom(n,1,.5)
  y <- cos(2*pi*x) + z1 + rnorm(n,sd=.25)
  z1 <- factor(z1)
  z2 <- factor(z2)

  bw <- npregbw(y~x+z1+z2,
                regtype="lc",
                bwmethod="cv.ls",
                nmulti=1)

  summary(bw)

  model <- npreg(bws=bw,
                 gradients=FALSE)

  summary(model)

  npRmpi.quit()
  ## npRmpi.quit(force=TRUE)
} else {
  message("Skipping MPI spawn in check context (set NP_RMPI_RUN_MPI_EXAMPLES_IN_CHECK=1).")
}
} % enddontrun
}
\keyword{nonparametric}
