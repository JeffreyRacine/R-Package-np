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
    bandwidth is selected on \code{x} by least-squares cross-validation
    (\code{bwmethod = "cv.ls"}, overridable through \code{\dots}).
    If only one is supplied it is reused for both samples without a search.
    If both are supplied their effective bandwidths, kernels and support
    specifications must agree; incompatible specifications are rejected.
    This test requires fixed bandwidths and symmetric kernels without
    boundary normalization. Nearest-neighbour bandwidths, beta kernels and
    finite continuous normalization bounds are rejected: the test's directed
    cross sum and variance formula assume symmetry. Unnormalized Gaussian,
    Epanechnikov and uniform kernels remain supported, including mixed
    continuous/categorical data. Other estimators' kernel support is unchanged.
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
  Documentation guide: see \code{\link{np.kernels}} for kernels, \code{\link{np.options}} for global options, and \code{\link{plot}} for plotting options.


  \code{npdeneqtest} computes the integrated squared density difference
  between the estimated densities/probabilities of two samples having
  identical variables/datatypes. See Li, Maasoumi, and Racine (2009) for
  details.

  The same bandwidth is held fixed in both samples and every bootstrap
  replication. Automatic selection follows the archived simulation drivers
  used for the paper's tables: select on the first sample. The methodological
  description in the paper instead discusses pooled-sample selection.
  First-sample selection makes sample ordering relevant when a bandwidth is
  selected automatically. To choose on the pool explicitly, construct a
  bandwidth with \code{npudensbw(dat = rbind(x, y), bwmethod = "cv.ls")}
  and supply it as \code{bw.x}.

  With fixed, unbounded bandwidths, pooled bootstrap samples are represented
  by multiplicities on their common support and processed in bounded groups.
  Multi-response kernel sums then reuse the fixed geometry for the first- and
  second-power contractions without constructing a dense kernel matrix.
  Bounded-kernel and nearest-neighbour bandwidth routes retain literal
  duplicate resampling because their boundary corrections or sample-owned
  radii cannot in general be represented by fixed-support counts.

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
set.seed(1234)

## Distributions are equal

n <- 250

sample.A <- data.frame(x=rnorm(n))
sample.B <- data.frame(x=rnorm(n))

npdeneqtest(sample.A,sample.B,B=99)

if (interactive()) Sys.sleep(5)

## Distributions are unequal

sample.A <- data.frame(x=rnorm(n))
sample.B <- data.frame(x=rchisq(n,df=5))

npdeneqtest(sample.A,sample.B,B=99)

## Mixed datatypes, distributions are equal

sample.A <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))
sample.B <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))

npdeneqtest(sample.A,sample.B,B=99)

if (interactive()) Sys.sleep(5)

## Mixed datatypes, distributions are unequal

sample.A <- data.frame(a=rnorm(n),b=factor(rbinom(n,2,.5)))
sample.B <- data.frame(a=rnorm(n,sd=10),b=factor(rbinom(n,2,.25)))

npdeneqtest(sample.A,sample.B,B=99)
} % enddontrun
}

\keyword{ nonparametric }
