%
\name{npdeneqtest}
\alias{npdeneqtest}
\title{ Kernel Consistent Density Equality Test with Mixed Data Types }
\description{
  \code{npdeneqtest} implements a consistent integrated squared
  difference test for equality of densities as described in Li, Maasoumi,
  and Racine (2009).
}
\usage{
npdeneqtest(x = NULL,
            y = NULL,
            bw.x = NULL,
            bw.y = NULL,
            B = 399,
            random.seed = 42,
            \dots)
}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the two samples and any supplied bandwidths.
}

  \item{bw.x,bw.y}{
    optional common bandwidth specification. If neither is supplied, one
    bandwidth is selected on the pooled sample by least-squares cross-validation
    (\code{bwmethod = "cv.ls"}, overridable through \code{\dots}), then
    rescaled to the harmonic-mean sample size as described below.
    If only one is supplied it is reused for both samples without a search.
    If both are supplied their effective bandwidths, kernels and support
    specifications must agree; incompatible specifications are rejected.
    This test requires fixed bandwidths and symmetric kernels without
    boundary normalization. Nearest-neighbour bandwidths, beta kernels and
    finite continuous normalization bounds are rejected: the test's directed
    cross sum and variance formula assume symmetry. Unnormalized Gaussian,
    Epanechnikov and uniform kernels remain supported, including mixed
    continuous/categorical data. An active ordered Racine--Li--Yan kernel is
    also rejected because it is asymmetric and has a nonconstant self-diagonal.
    Other estimators' kernel support is unchanged.
  }


  \item{x,y}{
    data frames for the two samples for which one wishes to
    test equality of densities. The variables in each data
    frame must be the same (i.e. have identical names).
  }


\subsection{Bootstrap Controls}{
These arguments control bootstrap replication and reproducibility settings.
}

  \item{B}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }


  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }


\subsection{Additional Arguments}{
Further arguments are passed to the bandwidth-selection routines used by the test.
}

  \item{\dots}{ additional arguments supplied to specify the bandwidth
    type, kernel types, and so on. These selection arguments are used only
    when both \code{bw.x} and \code{bw.y} are absent; they do not modify a
    supplied bandwidth specification. To
    do this, you may specify any of \code{bwscaling}, \code{bwtype},
    \code{ckertype}, \code{ckerorder}, \code{ukertype},
    \code{okertype}.}


}
\section{Bootstrap P-Values}{
  The bootstrap upper-tail p-value is the proportion of bootstrap statistics
  greater than or equal to the observed statistic. Ties count as at least as
  extreme, giving the conservative tie convention for discrete statistics.
  This is the ordinary bootstrap proportion, without an add-one correction;
  a reported zero is possible with a finite number of replications.
}
\value{
  \code{npdeneqtest} returns an object of type \code{deneqtest} with the
  following components

  \item{Tn}{ the (standardized) statistic \code{Tn} }
  \item{In}{ the (unstandardized) statistic \code{In} }
  \item{Tn.bootstrap}{ contains the bootstrap replications of \code{Tn} }
  \item{In.bootstrap}{ contains the bootstrap replications of \code{In} }
  \item{Tn.P}{ the P-value of the \code{Tn} statistic }
  \item{In.P}{ the P-value of the \code{In} statistic }
  \item{boot.num}{ number of bootstrap replications }

  \code{\link{summary}} supports object of type \code{deneqtest}.
}
\section{Book And Method Pointers}{
  \code{npdeneqtest} compares two distributions by estimating the
  integrated squared difference between their density or probability
  functions. The two samples must have matching variables and data
  types because the statistic compares like with like.

  For book-length background, see Li and Racine (2007), Chapter 12
  \emph{Model Specification Tests}, especially Section 12.2, and Chapter
  13 \emph{Nonsmoothing Tests}, especially Section 13.2. For density and
  mixed-data density context, see Racine (2019), Chapters 2 and 3.
}

\references{
  Li, Q. and E. Maasoumi and J.S. Racine (2009), \dQuote{A Nonparametric
    Test for Equality of Distributions with Mixed Categorical and
    Continuous Data,} Journal of Econometrics, 148, pp 186-200.
}
\author{
    Tristen Hayfield \email{tristen.hayfield@gmail.com}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, \code{\link{plot}} for plotting options, and \code{\link{npRmpi.init}} for interactive/cluster MPI startup. See \code{\link{npRmpi.init}} details for performance tradeoffs (message passing/startup mode) and the \code{inst/Rprofile} manual-broadcast template.


  \code{npdeneqtest} computes the integrated squared density difference
  between the estimated densities/probabilities of two samples having
  identical variables/datatypes. See Li, Maasoumi, and Racine (2009) for
  details.

  The same bandwidth is held fixed in both samples and every bootstrap
  replication. Automatic selection uses one pooled-sample search, followed by
  a \emph{harmonic-mean reference-size adjustment}. Write
  \eqn{N=n_1+n_2} and \eqn{n_H=2/(1/n_1+1/n_2)}, where the sample sizes count
  complete observations. Retaining the pooled spread estimates, continuous
  bandwidths are multiplied by \eqn{(N/n_H)^{1/(2P+q)}} and categorical
  bandwidths by \eqn{(N/n_H)^{2/(2P+q)}}, following the scale-factor convention
  of \code{\link{npudensbw}}. Here \eqn{P} is the continuous kernel order
  (two by default), and \eqn{q} is the number of continuous variables.
  For the default second-order kernel the continuous exponent is
  \eqn{1/(4+q)}. Categorical parameters cannot exceed their kernel's permitted
  maximum: an adjusted value above that maximum is set to the maximum and a
  warning names the affected variables.

  Under a common density, the leading variance coefficient of the equally
  weighted average of the two density-estimation risks is
  \eqn{(1/n_1+1/n_2)/2=1/n_H}. This motivates the common reference size;
  it is not a claim of optimal test power or exact finite-sample calibration.
  The mixed-data and higher-order adjustments use the existing bandwidth
  scaling convention. This rule differs from the archived first-sample
  simulation implementation and from unadjusted pooled selection.
  To use an unadjusted pooled bandwidth, construct it with
  \code{npudensbw(dat = rbind(x, y), bwmethod = "cv.ls")} and supply it as
  \code{bw.x}. Explicit \code{bw.x}/\code{bw.y} are never rescaled.
  Likewise, \code{bws = ..., bandwidth.compute = FALSE} in \code{\dots}
  remains a manual specification, with scale factors interpreted using
  \code{x}'s complete observations rather than the pool.

  Categorical coordinates use the union of the declared levels in both samples
  and supplied bandwidth objects, including unused levels. The same domain
  and effective kernel metadata are retained in every contraction and bootstrap
  draw. Ordered domains must have an unambiguous common order; otherwise declare
  the same ordered levels in both samples. Automatic selection uses both
  samples with this common metadata.
  A raw numeric bandwidth retains the default raw kernel-sum conventions;
  a density bandwidth object retains its selected density-kernel conventions,
  including normalized ordered Li--Racine masses. Selection arguments in
  \code{\dots} do not override an explicitly supplied bandwidth specification.

  With fixed, unbounded bandwidths, pooled bootstrap samples are represented
  by multiplicities on their common support and processed in bounded groups
  on each MPI rank. Multi-response kernel sums then reuse the fixed geometry
  for the first- and second-power contractions without constructing a dense
  kernel matrix. Boundary-normalized and nearest-neighbour specifications are
  not supported by this density-equality statistic.

}

\section{Usage Issues}{

  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.

  It is crucial that both data frames have the same variable names.

  Incomplete rows are omitted independently from each sample before bandwidth
  selection, bootstrap resampling and computation of the statistic. Factor
  levels are retained. Each sample must contain at least two complete
  observations. Supplying one bandwidth reuses it for both samples; it is not
  reselected when the other bandwidth argument is missing.

}

\seealso{
  \code{\link{np.kernels}}, \code{\link{np.options}}, \code{\link{plot}}
  \code{\link{npdeptest},\link{npsdeptest},\link{npsymtest},\link{npunitest}}
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
set.seed(42)

n <- 100

sample.A <- data.frame(x=rnorm(n))
sample.B <- data.frame(x=rnorm(n))

output <- npdeneqtest(sample.A,sample.B,B=29)

output


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
