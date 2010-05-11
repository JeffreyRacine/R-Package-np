% $ID$
\name{npregbw}
\alias{npregbw}
\alias{npregbw.NULL}
\alias{npregbw.default}
\alias{npregbw.formula}
\alias{npregbw.rbandwidth}

\title{Kernel Regression Bandwidth Selection with Mixed Data Types}

\description{
  \code{npregbw} computes a bandwidth object for a
  \eqn{p}-variate kernel regression estimator defined over mixed
  continuous and discrete (unordered, ordered) data using expected
  Kullback-Leibler cross-validation, or least-squares cross validation
  using the method of Racine and Li (2004) and Li and Racine (2004).
}

\usage{
npregbw(\dots)

\method{npregbw}{formula}(formula, data, subset, na.action, call, \dots)

\method{npregbw}{NULL}(xdat = stop("invoked without data 'xdat'"),
        ydat = stop("invoked without data 'ydat'"),
        bws,
        \dots)

\method{npregbw}{default}(xdat = stop("invoked without data 'xdat'"),
        ydat = stop("invoked without data 'ydat'"),
        bws,
        bandwidth.compute = TRUE,
        nmulti,
        remin,
        itmax,
        ftol,
        tol,
        small,
        regtype,
        bwmethod,
        bwscaling,
        bwtype,
        ckertype,
        ckerorder,
        ukertype,
        okertype,
        \dots)

\method{npregbw}{rbandwidth}(xdat = stop("invoked without data 'xdat'"),
        ydat = stop("invoked without data 'ydat'"),
        bws,
        bandwidth.compute = TRUE,
        nmulti,
        remin = TRUE,
        itmax = 10000,
        ftol = 1.19209e-07,
        tol = 1.49012e-08,
        small = 2.22045e-16,
        \dots)
}
\arguments{
  \item{formula}{
    a symbolic description of variables on which bandwidth selection is
    to be performed. The details of constructing a formula are
    described below.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which the
    function is called.
  }

  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process. 
  }

  \item{na.action}{
    a function which indicates what should happen when the data contain
   \code{NA}s. The default is set by the \code{\link{na.action}} setting of options, and is
    \code{\link{na.fail}} if that is unset. The (recommended) default is
    \code{\link{na.omit}}.  
  }
  
  \item{call}{
    the original function call. This is passed internally by
    \code{\link{np}} when a bandwidth search has been implied by a call to
    another function. It is not recommended that the user set this.
  }

  \item{xdat}{
    a \eqn{p}-variate data frame of regressors which bandwidth selection will
    be performed. The data types may be continuous, discrete (unordered
    and ordered factors), or some combination thereof.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{xdat}. 
  }

  \item{bws}{
    a bandwidth specification. This can be set as a \code{rbandwidth}
    object returned from a previous invocation, or as a vector of
    bandwidths, with each element \eqn{i} corresponding to the bandwidth
    for column \eqn{i} in \code{xdat}. In either case, the bandwidth
    supplied will serve as a starting point in the numerical search for
    optimal bandwidths. If specified as a vector, then additional
    arguments will need to be supplied as necessary to specify the
    bandwidth type, kernel types, selection methods, and so on. This can
    be left unset.
  }

  \item{\dots}{
    additional arguments supplied to specify the bandwidth type,
    kernel types, selection methods, and so on, detailed below.
  }
  
  \item{regtype}{
    a character string specifying which type of kernel regression
    estimator to use. \code{lc} specifies a local-constant estimator
    (Nadaraya-Watson) and \code{ll} specifies a local-linear
    estimator. Defaults to \code{lc}.
  }
  
  \item{bwmethod}{
    which method to use to select bandwidths. \code{cv.aic} specifies
    expected Kullback-Leibler cross-validation (Hurvich, Simonoff, and
    Tsai (1998)), and \code{cv.ls} specifies least-squares
    cross-validation. Defaults to \code{cv.ls}.
  }

  \item{bwscaling}{ a logical value that when set to \code{TRUE} the
    supplied bandwidths are interpreted as `scale factors'
    (\eqn{c_j}{c[j]}), otherwise when the value is \code{FALSE} they are
    interpreted as `raw bandwidths' (\eqn{h_j}{h[j]} for continuous data
    types, \eqn{\lambda_j}{lambda[j]} for discrete data types). For
    continuous data types, \eqn{c_j}{c[j]} and \eqn{h_j}{h[j]} are
    related by the formula \eqn{h_j = c_j \sigma_j n^{-1/(2P+l)}}{h[j] =
    c[j]*sigma[j]*n^(-1/(2*P+l))}, where \eqn{\sigma_j}{sigma[j]} is an
    adaptive measure of spread of continuous variable \eqn{j} defined as
    min(standard deviation, interquartile range/1.349), \eqn{n} the
    number of observations, \eqn{P}{P} the order of the kernel, and
    \eqn{l}{l} the number of continuous variables. For discrete data
    types, \eqn{c_j}{c[j]} and \eqn{h_j}{h[j]} are related by the
    formula \eqn{h_j = c_jn^{-2/(2P+l)}}{h[j] = c[j]*n^(-2/(2*P+l))},
    where here \eqn{j}{[j]} denotes discrete variable \eqn{j}{j}.
    Defaults to \code{FALSE}.  }

  \item{bwtype}{
    character string used for the continuous variable bandwidth type,
    specifying the type of bandwidth to compute and return in the
    \code{bandwidth} object. Defaults to \code{fixed}. Option
    summary:\cr
    \code{fixed}: compute fixed bandwidths \cr
    \code{generalized_nn}: compute generalized nearest neighbors \cr
    \code{adaptive_nn}: compute adaptive nearest neighbors
    
  }

  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{rbandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
  }


  \item{ckertype}{
    character string used to specify the continuous kernel type.
    Can be set as \code{gaussian}, \code{epanechnikov}, or
    \code{uniform}. Defaults to \code{gaussian}.
  }

  \item{ckerorder}{
    numeric value specifying kernel order (one of
    \code{(2,4,6,8)}). Kernel order specified along with a
  \code{uniform} continuous kernel type will be ignored. Defaults to
  \code{2}.
  }

  \item{ukertype}{
    character string used to specify the unordered categorical kernel type.
    Can be set as \code{aitchisonaitken} or \code{liracine}. Defaults to
    \code{aitchisonaitken}.
  }

  \item{okertype}{
    character string used to specify the ordered categorical kernel type.
    Can be set as \code{wangvanryzin} or \code{liracine}. Defaults to
    \code{wangvanryzin}.
  }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points. Defaults to \code{min(5,ncol(xdat))}.
  }

  \item{remin}{
    a logical value which when set as \code{TRUE} the search routine
    restarts from located minima for a minor gain in accuracy. Defaults
    to \code{TRUE}.
  }

  \item{itmax}{
    integer number of iterations before failure in the numerical
    optimization routine. Defaults to \code{10000}.
  }

  \item{ftol}{
    tolerance on the value of the cross-validation function 
    evaluated at located minima. Defaults to \code{1.19e-07
      (FLT_EPSILON)}.
  }

  \item{tol}{
    tolerance on the position of located minima of the 
    cross-validation function. Defaults to \code{1.49e-08
      (sqrt(DBL_EPSILON))}.
  }

  \item{small}{ a small number, at about the precision of the data type
    used. Defaults to \code{2.22e-16 (DBL_EPSILON)}.

  }

}
\details{
  \code{npregbw} implements a variety of methods for choosing
  bandwidths for multivariate (\eqn{p}-variate) regression data defined
  over a set of possibly continuous and/or discrete (unordered, ordered)
  data. The approach is based on Li and Racine (2003) who employ
  \sQuote{generalized product kernels} that admit a mix of continuous
  and discrete data types.

  The cross-validation methods employ multivariate numerical search
  algorithms (direction set (Powell's) methods in multidimensions).

  Bandwidths can (and will) differ for each variable which is, of
  course, desirable.

  Three classes of kernel estimators for the continuous data types are
  available: fixed, adaptive nearest-neighbor, and generalized
  nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
  each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
  density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
  with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
  are constant over the support of \eqn{x}.

  \code{npregbw} may be invoked \emph{either} with a formula-like
  symbolic 
  description of variables on which bandwidth selection is to be
  performed \emph{or} through a simpler interface whereby data is passed
  directly to the function via the \code{xdat} and \code{ydat}
  parameters. Use of these two interfaces is \bold{mutually exclusive}.

  Data contained in the data frame \code{xdat} may be a mix of
  continuous (default), unordered discrete (to be specified in the data
  frame \code{xdat} using \code{\link{factor}}), and ordered discrete
  (to be specified in the data frame \code{xdat} using
  \code{\link{ordered}}). Data can be entered in an arbitrary order and
  data types will be detected automatically by the routine (see
  \code{\link{np}} for details).

  Data for which bandwidths are to be estimated may be specified
  symbolically. A typical description has the form \code{dependent data
    ~ explanatory data},
  where \code{dependent data} is a univariate response, and
  \code{explanatory data} is a
  series of variables specified by name, separated by 
  the separation character '+'. For example, \code{ y1 ~ x1 + x2 }
  specifies that the bandwidths for the regression of response \code{y1}
  and 
  nonparametric regressors \code{x1} and \code{x2} are to be estimated.
  See below for further examples. 


  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

}
\value{
  \code{npregbw} returns a \code{rbandwidth} object, with the
  following components:

  \item{bw}{ bandwidth(s), scale factor(s) or nearest neighbours for the
    data, \code{xdat} }
  \item{fval}{ objective function value at minimum }

  if \code{bwtype} is set to \code{fixed}, an object containing bandwidths
  (or scale factors if \code{bwscaling = TRUE}) is returned. If it is set to
  \code{generalized_nn} or \code{adaptive_nn}, then instead the \eqn{k}th nearest
  neighbors are returned for the continuous variables while the discrete
  kernel bandwidths are returned for the discrete variables. Bandwidths
  are stored under the component name \code{bw}, with each
  element \eqn{i} corresponding to column \eqn{i} of input data
  \code{xdat}.

  The functions \code{\link{predict}}, \code{\link{summary}}, and \code{\link{plot}} support
  objects of this class.
}

\references{
  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and Q. Li and J.S. Racine (2007), \dQuote{Nonparametric
  estimation of regression functions in the presence of irrelevant
  regressors,} The Review of Economics and Statistics, 89, 784-789.

  Hurvich, C.M. and J.S. Simonoff and C.L. Tsai (1998),
  \dQuote{Smoothing parameter selection in nonparametric regression
    using an improved Akaike information criterion,} Journal of the
  Royal Statistical Society B, 60, 271-293.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
    nonparametric regression,} Statistica Sinica, 14, 485-512.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.
  
  Racine, J.S. and Q. Li (2004), \dQuote{Nonparametric estimation of regression
    functions with both categorical and continuous data,} Journal of
  Econometrics, 119, 99-130.
  
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

  Caution: multivariate data-driven bandwidth selection methods are, by
  their nature, \emph{computationally intensive}. Virtually all methods
  require dropping the \eqn{i}th observation from the data set, computing an
  object, repeating this for all observations in the sample, then
  averaging each of these leave-one-out estimates for a \emph{given}
  value of the bandwidth vector, and only then repeating this a large
  number of times in order to conduct multivariate numerical
  minimization/maximization. Furthermore, due to the potential for local
  minima/maxima, \emph{restarting this procedure a large number of times may
  often be necessary}. This can be frustrating for users possessing
  large datasets. For exploratory purposes, you may wish to override the
  default search tolerances, say, setting ftol=.01 and tol=.01 and
  conduct multistarting (the default is to restart min(5, ncol(xdat))
  times) as is done for a number of examples. Once the procedure
  terminates, you can restart search with default tolerances using those
  bandwidths obtained from the less rigorous search (i.e., set
  \code{bws=bw} on subsequent calls to this routine where \code{bw} is
  the initial bandwidth object).  A version of this package using the
  \code{Rmpi} wrapper is under development that allows one to deploy
  this software in a clustered computing environment to facilitate
  computation involving large datasets.
}

\seealso{
  \code{\link{npreg}}
}
\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we compute a
# Bivariate nonparametric regression estimate for Giovanni Baiocchi's
# Italian income panel (see Italy for details)

data("Italy")
attach(Italy)

# Compute the least-squares cross-validated bandwidths for the local
# constant estimator (default)

bw <- npregbw(formula=gdp~ordered(year), tol=.1, ftol=.1)

summary(bw)

\dontrun{

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Supply your own bandwidth...

bw <- npregbw(formula=gdp~ordered(year), bws=c(0.75),
              bandwidth.compute=FALSE)

summary(bw)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Treat year as continuous and supply your own scaling factor c in
# c sigma n^{-1/(2p+q)}

bw <- npregbw(formula=gdp~year, bws=c(1.06),
              bandwidth.compute=FALSE, 
              bwscaling=TRUE)

summary(bw)

# Note - see also the example for npudensbw() for more extensive
# multiple illustrations of how to change the kernel function, kernel
# order, bandwidth type and so forth.

detach(Italy)

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we compute a
# Bivariate nonparametric regression estimate for Giovanni Baiocchi's
# Italian income panel (see Italy for details)

data("Italy")
attach(Italy)

# Compute the least-squares cross-validated bandwidths for the local
# constant estimator (default)

bw <- npregbw(xdat=ordered(year), ydat=gdp, tol=.1, ftol=.1)

summary(bw)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Supply your own bandwidth...

bw <- npregbw(xdat=ordered(year), ydat=gdp, bws=c(0.75),
              bandwidth.compute=FALSE)

summary(bw)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Treat year as continuous and supply your own scaling factor c in
# c sigma n^{-1/(2p+q)}

bw <- npregbw(xdat=year, ydat=gdp, bws=c(1.06),
              bandwidth.compute=FALSE, 
              bwscaling=TRUE)

summary(bw)

# Note - see also the example for npudensbw() for more extensive
# multiple illustrations of how to change the kernel function, kernel
# order, bandwidth type and so forth.

detach(Italy)
} % enddontrun
}
\keyword{nonparametric}
