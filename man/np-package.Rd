\name{np}
\alias{np}
\docType{package}
\title{Nonparametric Kernel Smoothing Methods for Mixed Data Types}
\description{

	This package provides a variety of nonparametric and semiparametric
	kernel methods that seamlessly handle a mix of continuous, unordered,
	and ordered factor data types (unordered and ordered factors are often
	referred to as \sQuote{nominal} and \sQuote{ordinal} categorical
	variables respectively). A vignette containing many of the examples
	found in the help files accompanying the \pkg{np} package that is
	intended to serve as a gentle introduction to this package can be
	accessed via \code{vignette("np", package="np")}.

  For a listing of all routines in the \pkg{np} package type:
  \sQuote{library(help="np")}.  

  Bandwidth selection is a key aspect of sound nonparametric and
  semiparametric kernel estimation. \code{np} is designed from the
  ground up to make bandwidth selection the focus of attention. To this
  end, one typically begins by creating a \sQuote{bandwidth object}
  which embodies all aspects of the method, including specific kernel
  functions, data names, data types, and the like. One then passes these
  bandwidth objects to other functions, and those functions can grab the
  specifics from the bandwidth object thereby removing potential
  inconsistencies and unnecessary repetition. Furthermore, many
  functions such as \code{npplot} can work with the bandwidth object
  directly without having to do the subsequent companion function
  evaluation.

  As of \code{np} version 0.20-0, we allow the user to combine these
  steps. When using \code{np} versions 0.20-0 and higher, if the first
  step (bandwidth selection) is not performed explicitly then the second
  step will automatically call the omitted first step bandwidth selector
  using defaults unless otherwise specified, and the bandwidth object
  could then be retrieved retroactively if so desired via
  \code{objectname$bws}. Furthermore, options for bandwidth selection
  will be passed directly to the bandwidth selector function. Note that
  the combined approach would not be a wise choice for certain
  applications such as when bootstrapping (as it would involve
  unnecessary computation since the bandwidths would properly be those
  for the original sample and not the bootstrap resamples) or when
  conducting quantile regression (as it would involve unnecessary
  computation when different quantiles are computed from the same
  conditional cumulative distribution estimate).
  
  There are two ways in which you can interact with functions in
  \code{np}, either i) using data frames, or ii) using a formula
  interface, where appropriate.

  To some, it may be natural to use the data frame interface.  The R
  \code{\link{data.frame}} function preserves a variable's type once it
  has been cast (unlike \code{\link{cbind}}, which we avoid for this
  reason).  If you find this most natural for your project, you first
  create a data frame casting data according to their type (i.e., one of
  continuous (default, \code{\link{numeric}}), \code{\link{factor}},
  \code{\link{ordered}}). Then you would simply pass this data frame to
  the appropriate \code{np} function, for example
  \code{npudensbw(dat=data)}.

  To others, however, it may be natural to use the formula interface
  that is used for the regression examples, among others. For
  nonparametric regression functions such as \code{\link{npreg}}, you
  would proceed as you would using \code{\link{lm}} (e.g., \code{bw <-
  npregbw(y~factor(x1)+x2))} except that you would of course not need to
  specify, e.g., polynomials in variables, interaction terms, or create
  a number of dummy variables for a factor.  Every function in np
  supports both interfaces, where appropriate.

  Note that if your factor is in fact a character string such as, say,
	\code{X} being either \code{"MALE"} or \code{"FEMALE"}, np will handle
	this directly, i.e., there is no need to map the string values into
	unique integers such as (0,1). Once the user casts a variable as a
	particular data type (i.e., \code{\link{factor}},
	\code{\link{ordered}}, or continuous (default,
	\code{\link{numeric}})), all subsequent methods automatically detect
	the type and use the appropriate kernel function and method where
	appropriate.

	All estimation methods are fully multivariate, i.e., there are no
	limitations on the number of variables one can model (or number of
	observations for that matter). Execution time for most routines is,
	however, exponentially increasing in the number of observations and
	increases with the number of variables involved.

	Nonparametric methods include unconditional density (distribution),
	conditional density (distribution), regression, mode, and quantile
	estimators along with gradients where appropriate, while
	semiparametric methods include single index, partially linear, and
	smooth (i.e., varying) coefficient models.
	
	A number of tests are included such as consistent specification tests
	for parametric regression and quantile regression models along with
	tests of significance for nonparametric regression.
	
	A variety of bootstrap methods for computing standard errors,
	nonparametric confidence bounds, and bias-corrected bounds are
	implemented.

	A variety of bandwidth methods are implemented including fixed,
	nearest-neighbor, and adaptive nearest-neighbor.

	A variety of data-driven methods of bandwidth selection are
	implemented, while the user can specify their own bandwidths should
	they so choose (either a raw bandwidth or scaling factor).

	A flexible plotting utility, \code{\link{npplot}}, facilitates
	graphing of multivariate objects. An example for creating postscript
	graphs using the \code{\link{npplot}} utility and pulling this into a
	LaTeX document is provided.

	The function \code{\link{npksum}} allows users to create or implement
	their own kernel estimators or tests should they so desire.

	The underlying functions are written in C for computational
	efficiency. Despite this, due to their nature, data-driven bandwidth
	selection methods involving multivariate numerical search can be
	time-consuming, particularly for large datasets. A version of this
	package using the \code{Rmpi} wrapper is under development that allows
	one to deploy this software in a clustered computing environment to
	facilitate computation involving large datasets.

  To cite the \code{np} package, type \code{citation("np")} from within
  \code{R} for details.

}

\details{
	The kernel methods in \code{np} employ the so-called `generalized
	product kernels' found in Hall, Racine, and Li (2004), Li and Racine
  (2003), Li and Racine (2004), Li and Racine (2007), Ouyang, Li,
  and Racine (2006), and Racine and Li (2004), among others. For
  details on a particular method, kindly refer to the original
  references listed above.

  We briefly describe the particulars of various univariate kernels used
  to generate the generalized product kernels that underlie the kernel
  estimators implemented in the \code{np} package. In a nutshell, the
  generalized kernel functions that underlie the kernel estimators in
  \code{np} are formed by taking the product of univariate kernels such
  as those listed below. When you cast your data as a particular type
  (continuous, factor, or ordered factor) in a data frame or formula,
  the routines will automatically recognize the type of variable being
  modelled and use the appropriate kernel type for each variable in the
  resulting estimator.

	\describe{

		\item{Second Order Gaussian (\eqn{x} is continuous)}{
			\eqn{k(z) = \exp(-z^2/2)/\sqrt{2\pi}}{k(z) =
				exp(-z^2/2)/sqrt{2*pi}}, where \eqn{z=(x_i-x)/h}{z=(x[i]-x)/h},
      and \eqn{h>0}.
    }

    \item{Second Order Epanechnikov (\eqn{x} is continuous)}{
      \eqn{k(z) =  3\left(1 - z^2/5\right)/(4\sqrt{5})}{
        k(z) =  3(1 - z^2/5)/(4*sqrt(5))}
      if \eqn{z^2<5}, \eqn{0} otherwise, where
      \eqn{z=(x_i-x)/h}{z=(x[i]-x)/h}, and \eqn{h>0}.
    }

    \item{Uniform (\eqn{x} is continuous)}{
      \eqn{k(z) = 1/2} if \eqn{|z|<1}, \eqn{0} otherwise, where
      \eqn{z=(x_i-x)/h}{z=(x[i]-x)/h}, and \eqn{h>0}.
    }

    \item{Aitchison and Aitken (\eqn{x} is a (discrete) factor)}{
      \eqn{l(x_i,x,\lambda) = 1 - \lambda}{l(x[i],x,lambda) = 1 -
        lambda} if \eqn{x_i=x}{x[i] = x}, and
      \eqn{\lambda/(c-1)}{lambda/(c-1)} if \eqn{x_i \neq x}{x[i] != x},
      where \eqn{c} is the number of (discrete) outcomes assumed by the
      factor \eqn{x}. Note that \eqn{\lambda}{lambda} must lie between
      \eqn{0} and \eqn{(c-1)/c}.
    }

    \item{Wang and van Ryzin (\eqn{x} is a (discrete) ordered factor)}{
      \eqn{l(x_i,x,\lambda) = 1 - \lambda}{l(x[i],x,lambda) = 1 -
        lambda} if \eqn{|x_i-x|=0}{|x[i] - x| = 0}, and
      \eqn{(1-\lambda)\lambda^{|x_i-x|}/2}{(1-lambda)lambda^|x_i-x|/2}
      if \eqn{|x_i - x|\ge1}{|x[i] - x|>=1}. Note that
      \eqn{\lambda}{lambda} must lie between \eqn{0} and \eqn{1}.
    }
      

    \item{Li and Racine (\eqn{x} is a (discrete) factor)}{
      \eqn{l(x_i,x,\lambda) = 1 }{l(x[i],x,lambda) =
        1} if \eqn{x_i=x}{x[i] = x}, and
      \eqn{\lambda}{lambda} if \eqn{x_i \neq x}{x[i] != x}. Note that
      \eqn{\lambda}{lambda} must lie between \eqn{0} and \eqn{1}.
      }

    \item{Li and Racine (\eqn{x} is a (discrete) ordered factor)}{
      \eqn{l(x_i,x,\lambda) = 1}{l(x[i],x,lambda) = 1} if
      \eqn{|x_i-x|=0}{|x[i] - x| = 0}, and
      \eqn{\lambda^{|x_i-x|}}{lambda^|x_i-x|} if \eqn{|x_i -
        x|\ge1}{|x[i] - x|>=1}. Note that \eqn{\lambda}{lambda} must lie
      between \eqn{0} and \eqn{1}.
    }

  }

  So, if you had two variables, \eqn{x_{i1}}{x1[i]} and
  \eqn{x_{i2}}{x2[i]}, and \eqn{x_{i1}}{x1[i]} was continuous while
  \eqn{x_{i2}}{x2[i]} was, say, binary (0/1), and you created a data
  frame of the form \code{X <- data.frame(x1,factor(x2))}, then the
  kernel function used by \code{np} would be
  \eqn{K(\cdot)=k(\cdot)\times l(\cdot)}{K(.)=k(.)*l(.)} where the
  particular kernel functions \eqn{k(\cdot)}{k(.)} and
  \eqn{l(\cdot)}{l(.)} would be, say, the second order Gaussian
  (\code{ckertype="gaussian"}) and Aitchison and Aitken
  (\code{ukertype="aitchisonaitken"}) kernels by default, respectively.

  Note that higher order continuous kernels (i.e., fourth, sixth, and
  eighth order) are derived from the second order kernels given above
  (see Li and Racine (2007) for details).

  For particulars on any given method, kindly see the references listed
  for the method in question.

}

\author{
  Tristen Hayfield <hayfield@phys.ethz.ch>, Jeffrey S. Racine
  <racinej@mcmaster.ca>

  Maintainer: Jeffrey S. Racine <racinej@mcmaster.ca>

  We are grateful to John Fox and Achim Zeleis for their valuable input
  and encouragement. We would like to gratefully acknowledge support
  from the Natural Sciences and Engineering Research Council of Canada
  (NSERC:www.nserc.ca), the Social Sciences and Humanities Research
  Council of Canada (SSHRC:www.sshrc.ca), and the Shared Hierarchical
  Academic Research Computing Network (SHARCNET:www.sharcnet.ca)

  
}
\references{

  Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and J.S. Racine and Q. Li (2004), \dQuote{Cross-validation and the
    estimation of conditional probability densities,} Journal of the
  American Statistical Association, 99, 1015-1026.

  Li, Q. and J.S. Racine (2003), \dQuote{Nonparametric estimation of
    distributions with categorical and continuous data,} Journal
  of Multivariate Analysis, 86, 266-292.

  Li, Q. and J.S. Racine (2004), \dQuote{Cross-validated local linear
  nonparametric regression,} Statistica Sinica, 14, 485-512.

 Ouyang, D. and Q. Li and J.S. Racine (2006), \dQuote{Cross-validation
   and the estimation of probability distributions with categorical data,}
 Journal of Nonparametric Statistics, 18, 69-100.

 Racine, J.S. and Q. Li (2004), \dQuote{Nonparametric estimation of
	regression functions with both categorical and continuous data,}
	Journal of Econometrics, 119, 99-130.
  
  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics:
		Theory and Practice,} Princeton University Press.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,} Cambridge
  University Press. 

  Scott, D.W. (1992), \emph{Multivariate Density Estimation. Theory,
  Practice and Visualization,} New York: Wiley.

  Silverman, B.W. (1986), \emph{Density Estimation,} London: Chapman and
  Hall.

  Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
    for discrete distributions,}  Biometrika, 68, 301-309.

}
\keyword{package}
