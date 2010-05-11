% $Id: np.conmode.Rd,v 1.38 2006/11/03 21:19:03 tristen Exp $
\name{npconmode}
\alias{npconmode}
\alias{npconmode.call}
\alias{npconmode.default}
\alias{npconmode.formula}
\alias{npconmode.conbandwidth}
\title{ Kernel Modal Regression with Mixed Data Types }
\description{
  \code{npconmode} performs kernel modal regression on mixed data,
  and finds the
  conditional mode given a set of training data, consisting of
  explanatory data and dependent data, and possibly evaluation data.
  Automatically computes various in sample and out of sample measures of
  accuracy.
}

\usage{
npconmode(bws, \dots)

\method{npconmode}{formula}(bws, data = NULL, newdata = NULL, \dots)

\method{npconmode}{call}(bws, \dots)

\method{npconmode}{default}(bws, txdat, tydat, \dots)

\method{npconmode}{conbandwidth}(bws,
          txdat = stop("invoked without training data 'txdat'"),
          tydat = stop("invoked without training data 'tydat'"),
          exdat,
          eydat,
          \dots)

}

\arguments{
  \item{bws}{
    a bandwidth specification. This can be set as a \code{conbandwidth}
    object returned from an invocation of \code{\link{npcdensbw}}
  }

  \item{\dots}{
    additional arguments supplied to specify the bandwidth type,
    kernel types, and so on, detailed below.
    This is necessary if you specify bws as a \eqn{p+q}-vector and not
    a \code{conbandwidth} object, and you do not desire the default behaviours.
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
    a \eqn{p}-variate data frame of explanatory data (conditioning data) used to
    calculate the regression estimators. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{tydat}{
    a one (1) dimensional vector of unordered or ordered factors,
    containing the dependent data. Defaults to the training data used to
    compute the bandwidth object.
  }

  \item{exdat}{
    a \eqn{p}-variate data frame of points on which the regression will be
    estimated (evaluation data). By default,
    evaluation takes place on the data provided by \code{txdat}.
  }

  \item{eydat}{
    a one (1) dimensional numeric or integer vector of the true values
    (outcomes) of the dependent variable. By default,
    evaluation takes place on the data provided by \code{tydat}.
  }



}

\value{
  \code{npconmode} returns a \code{conmode} object with the following
  components:

  \item{conmode}{ a vector of type \code{factor} (or \code{ordered
      factor}) containing the conditional mode at each evaluation
    point}
  \item{condens}{ a vector of numeric type containing the density
    estimates at each evaluation point}
  \item{xeval}{ a data frame of evaluation points }
  \item{yeval}{ a vector of type \code{factor} (or \code{ordered
      factor}) containing the actual outcomes, or \code{NA} if not
    provided }
  \item{confusion.matrix}{ the confusion matrix or \code{NA} if outcomes
    are not available }
  \item{CCR.overall}{ the overall correct
    classification ratio, or \code{NA} if outcomes are not available }
  \item{CCR.byoutcome}{ a numeric vector containing the correct
    classification ratio by outcome, or \code{NA} if outcomes are not
    available }
  \item{fit.mcfadden}{ the McFadden-Puig-Kerschner performance measure
    or \code{NA} if outcomes are not available }

  The functions \code{\link{mode}}, and \code{\link{fitted}} may be used to
  extract the conditional mode estimates, and the conditional density
  estimates at the conditional mode, respectively,
  from the resulting object. Also, \code{\link{summary}} supports
  \code{conmode} objects.
}
\references{
  Aitchison, J. and  C.G.G. Aitken (1976), \dQuote{Multivariate binary
    discrimination by the kernel method,} Biometrika, 63, 413-420.

  Hall, P. and J.S. Racine and Q. Li (2004), \dQuote{Cross-validation and the
    estimation of conditional probability densities,} Journal of the
  American Statistical Association, 99, 1015-1026.

  Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
  and Practice,} Princeton University Press.

  McFadden, D. and C. Puig and D. Kerschner (1977), \dQuote{Determinants
    of the long-run demand for electricity,} Proceedings of the
  American Statistical Association (Business and Economics Section),
  109-117. 

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
  Cambridge University Press.

  Scott, D.W. (1992), \emph{Multivariate Density Estimation. Theory,
  Practice and Visualization,} New York: Wiley.

  Silverman, B.W. (1986), \emph{Density Estimation,} London: Chapman and
  Hall.

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
}

\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we use the
# birthweight data taken from the MASS library, and compute a parametric
# logit model and a nonparametric conditional mode model. We then
# compare their confusion matrices and summary measures of
# classification ability.

library("MASS")
data("birthwt")
attach(birthwt)

# Fit a parametric logit model with low (0/1) as the dependent
# variable and age, lwt, and smoke (0/1) as the covariates

# From ?birthwt
# 'low' indicator of birth weight less than 2.5kg
# 'smoke' smoking status during pregnancy
# 'race' mother's race ('1' = white, '2' = black, '3' = other)
# 'ht' history of hypertension
# 'ui' presence of uterine irritability
# 'ftv' number of physician visits during the first trimester
# 'age' mother's age in years
# 'lwt' mother's weight in pounds at last menstrual period

model.logit <- glm(low~factor(smoke)+
                   factor(race)+
                   factor(ht)+
                   factor(ui)+
                   ordered(ftv)+
                   age+
                   lwt, 
                   family=binomial(link=logit))

# Generate the confusion matrix and correct classification ratio

cm <- table(low, ifelse(fitted(model.logit)>0.5, 1, 0))
ccr <- sum(diag(cm))/sum(cm)

# Now do the same with a nonparametric model.  Note - this may take a
# few minutes depending on the speed of your computer... we override the
# default search tolerances for speed considerations...

bw <- npcdensbw(formula=factor(low)~factor(smoke)+
                factor(race)+
                factor(ht)+
                factor(ui)+
                ordered(ftv)+
                age+
                lwt,
                tol=.1, ftol=.1)

model.np <- npconmode(bws=bw)

# Compare confusion matrices from the logit and nonparametric model

# Logit

cm
ccr

# Nonparametric
summary(model.np)

detach(birthwt)

\dontrun{

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we use the
# birthweight data taken from the MASS library, and compute a parametric
# logit model and a nonparametric conditional mode model. We then
# compare their confusion matrices and summary measures of
# classification ability.

library("MASS")
data("birthwt")
attach(birthwt)

# Fit a parametric logit model with low (0/1) as the dependent
# variable and age, lwt, and smoke (0/1) as the covariates

# From ?birthwt
# 'low' indicator of birth weight less than 2.5kg
# 'smoke' smoking status during pregnancy
# 'race' mother's race ('1' = white, '2' = black, '3' = other)
# 'ht' history of hypertension
# 'ui' presence of uterine irritability
# 'ftv' number of physician visits during the first trimester
# 'age' mother's age in years
# 'lwt' mother's weight in pounds at last menstrual period

model.logit <- glm(low~factor(smoke)+
                   factor(race)+
                   factor(ht)+
                   factor(ui)+
                   ordered(ftv)+
                   age+
                   lwt, 
                   family=binomial(link=logit))

# Generate the confusion matrix and correct classification ratio

cm <- table(low, ifelse(fitted(model.logit)>0.5, 1, 0))
ccr <- sum(diag(cm))/sum(cm)

# Now do the same with a nonparametric model...

X <- data.frame(factor(smoke), 
                factor(race), 
                factor(ht), 
                factor(ui), 
                ordered(ftv), 
                age, 
                lwt)

y <- factor(low)

# Note - this may take a few minutes depending on the speed of your
# computer... we override the default search tolerances for speed
# considerations...

bw <- npcdensbw(xdat=X, ydat=y, tol=.1, ftol=.1)

model.np <- npconmode(bws=bw)

# Compare confusion matrices from the logit and nonparametric model

# Logit

cm
ccr

# Nonparametric
summary(model.np)

detach(birthwt)
} % enddontrun
}
\keyword{ nonparametric }% at least one, from doc/KEYWORDS

