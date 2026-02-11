% $Id: np.condensity.Rd,v 1.67 2006/11/03 21:19:03 tristen Exp $
\name{npcdens}
\alias{npcdens}
\alias{npcdens.call}
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

\method{npcdens}{call}(bws, \dots)

\method{npcdens}{conbandwidth}(bws,
        txdat = stop("invoked without training data 'txdat'"),
        tydat = stop("invoked without training data 'tydat'"),
        exdat,
        eydat,
        gradients = FALSE,
        \dots)

\method{npcdens}{default}(bws, txdat, tydat, \dots)

}

\arguments{
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

  \item{gradients}{
    a logical value specifying whether to return estimates of the
    gradients at the evaluation points. Defaults to \code{FALSE}.
  }

  \item{\dots}{
    additional arguments supplied to specify the bandwidth type,
    kernel types, and so on.  This is necessary if you
    specify bws as a \eqn{p+q}-vector and not a \code{conbandwidth}
    object, and you do not desire the default behaviours. To do this,
    you may specify any of \code{bwmethod}, \code{bwscaling},
    \code{bwtype}, \code{cxkertype}, \code{cxkerorder},
    \code{cykertype}, \code{cykerorder}, \code{uxkertype},
    \code{uykertype}, \code{oxkertype}, \code{oykertype}, as described
    in \code{\link{npcdensbw}}. 

  }
  
  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
   \code{\link{npcdensbw}} was called.
  }

  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.  
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

}
\details{
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
  automatically by the routine (see \code{\link{npRmpi}} for details).

  A variety of kernels may be specified by the user. Kernels implemented
  for continuous data types include the second, fourth, sixth, and eighth
  order Gaussian and Epanechnikov kernels, and the uniform
  kernel. Unordered discrete data types use a variation on Aitchison and
  Aitken's (1976) kernel, while ordered data types use a variation of the
  Wang and van Ryzin (1981) kernel.

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
  \code{\link{npudens}}
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

mpi.bcast.cmd(data("Italy"),
              caller.execute=TRUE)
mpi.bcast.cmd(attach(Italy),
              caller.execute=TRUE)

mpi.bcast.cmd(bw <- npcdensbw(formula=gdp~ordered(year)),
              caller.execute=TRUE)

mpi.bcast.cmd(fhat <- npcdens(bws=bw),
              caller.execute=TRUE)

summary(fhat)

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
