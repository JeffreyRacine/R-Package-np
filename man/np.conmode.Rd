%
\name{npconmode}
\alias{npconmode}
\alias{npconmode.default}
\alias{npconmode.formula}
\alias{npconmode.conbandwidth}
\alias{predict.conmode}
\alias{plot.conmode}
\title{ Kernel Modal Regression with Mixed Data Types }
\description{
  \code{npconmode} performs kernel modal regression on mixed data,
  and finds the
  conditional mode given a set of training data, consisting of
  explanatory data and dependent data, and possibly evaluation data.
  Automatically computes various in sample and out of sample measures of
  accuracy.
}
\details{
  Documentation guide: see \code{\link{np.kernels}} for kernels,
  \code{\link{np.options}} for global options, and
  \code{\link{plot}} for plotting options.

  If \code{bws} is not an explicit \code{conbandwidth} object,
  \code{npconmode} computes conditional density bandwidths by forwarding the
  call to \code{\link{npcdensbw}}. The resulting \code{conbandwidth} object is
  stored unchanged in the returned object's \code{bws} component. Regression
  type, local-polynomial degree, and NOMAD/search metadata are also mirrored on
  the returned \code{conmode} object for convenient inspection and summary
  reporting; the bandwidth object remains the canonical source.

  For non-local-constant conditional-mode fits, \code{npconmode} constructs
  the full set of fitted probabilities over the discrete response support and
  projects them onto the probability simplex before selecting and reporting
  the modal outcome. Thus the probabilities used for modal selection are
  non-negative and sum to one over the discrete support. For binary outcomes
  this is the complement contract
  \eqn{\Pr(Y=\ell_2\mid X=x)=1-\Pr(Y=\ell_1\mid X=x)}. Local-constant fits
  are already proper by construction, so omitted \code{proper} resolves to
  \code{FALSE} for \code{regtype="lc"} and to \code{TRUE} otherwise.

  The \code{\link{predict}} method follows the usual S3
  \code{newdata} convention. With no evaluation arguments,
  \code{predict(fit)} extracts the stored modal class and
  \code{predict(fit, type="prob")} extracts stored class probabilities
  when the object was fitted with \code{probabilities=TRUE}. With
  \code{newdata}, \code{predict} evaluates the conditional mode at the
  supplied rows using the stored bandwidth object. Use
  \code{type="prob"} to return the full matrix of fitted class
  probabilities at the evaluation rows. Use \code{se.fit=TRUE} with
  \code{type="prob"} to return a list containing the class-probability matrix
  and the matching asymptotic standard-error matrix. Native evaluator
  arguments \code{exdat} and \code{eydat} remain available for advanced
  workflows and take precedence if supplied.
}

\usage{
npconmode(bws, 
          \dots)

\method{npconmode}{formula}(bws, 
          data = NULL,
          newdata = NULL, 
          \dots)

\method{npconmode}{default}(bws, 
          txdat, 
          tydat, 
          nomad = FALSE,
          proper = NULL,
          proper.control = list(),
          probabilities = FALSE,
          gradients = FALSE,
          level = NULL,
          \dots)

\method{npconmode}{conbandwidth}(bws,
          txdat = stop("invoked without training data 'txdat'"),
          tydat = stop("invoked without training data 'tydat'"),
          exdat,
          eydat,
          proper = NULL,
          proper.control = list(),
          probabilities = FALSE,
          gradients = FALSE,
          level = NULL,
          \dots)

\method{predict}{conmode}(object,
          newdata = NULL,
          type = c("class", "prob"),
          se.fit = FALSE,
          \dots)

\method{plot}{conmode}(x, \dots)

}

\arguments{
\subsection{Data, Bandwidth Inputs And Formula Interface}{
These arguments identify the bandwidth specification, formula/data interface, and training data.
}

  \item{bws}{
    a bandwidth specification. This can be set as a \code{conbandwidth}
    object returned from an invocation of \code{\link{npcdensbw}}
  }



  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(bws)}, typically the environment from which
    \code{\link{npcdensbw}} was called.
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

  \item{object}{
    an object of class \code{"conmode"} returned by \code{npconmode}.
  }

  \item{x}{
    an object of class \code{"conmode"} returned by \code{npconmode}.
  }


\subsection{Evaluation Data}{
These arguments control where the conditional mode is evaluated.
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



  \item{newdata}{
    An optional data frame in which to look for evaluation data. If
    omitted, the training data are used.
  }

  \item{type}{
    prediction type for \code{predict.conmode}. Use \code{"class"} to
    return modal class predictions and \code{"prob"} to return the full
    class-probability matrix. Probability prediction requires either a
    stored probability matrix or an evaluation path with
    \code{probabilities=TRUE}.
  }

  \item{se.fit}{
    a logical value for \code{predict.conmode}. If \code{TRUE} and
    \code{type="prob"}, return a list with components \code{fit} and
    \code{se.fit}, containing the class-probability matrix and matching
    asymptotic standard-error matrix. Standard errors are for the smooth
    class-probability target, not for the discrete modal class label.
  }


\subsection{Probability And Gradient Controls}{
These arguments control stored class probabilities, proper probability
normalization, and level-specific class-probability effects.
}

  \item{probabilities}{
    a logical value. If \code{TRUE}, store the full matrix of fitted
    probabilities over the discrete response support in the returned object.
    The default \code{FALSE} keeps returned objects compact. If
    \code{gradients=TRUE}, probabilities are stored automatically because
    they are the object-fed target for \code{\link{plot.conmode}}.
  }

  \item{gradients}{
    a logical value. If \code{TRUE}, compute and store class-probability
    gradients/effects for one response level with respect to the conditioning
    variables. These are derivatives for continuous conditioning variables and
    the corresponding smooth finite-change effects for discrete conditioning
    variables as produced by the underlying \code{\link{npcdens}} gradient
    route. The target is the smooth class probability
    \eqn{\Pr(Y=\ell\mid X=x)}, not the derivative of the discrete modal class
    label.
  }

  \item{level}{
    response level for class-probability gradients/effects and
    \code{plot.conmode} class-probability displays. If omitted, the default is
    the base/reference response level, \code{levels(y)[1]}, for both factor and
    ordered responses. Users should supply \code{level} explicitly when a
    particular binary or multinomial outcome probability is of substantive
    interest.
  }

  \item{proper}{
    a logical value or \code{NULL}. If \code{TRUE}, fitted probabilities over
    the discrete response support are projected onto the probability simplex
    before modal selection. If \code{FALSE}, legacy raw per-level values are
    used. If \code{NULL}, the default is estimator-aware:
    \code{FALSE} for \code{regtype="lc"} and \code{TRUE} otherwise.
  }

  \item{proper.control}{
    a list of controls for discrete probability properization. Currently
    \code{tol} may be supplied to set the numerical tolerance for
    non-negativity and unit-mass checks.
  }


\subsection{Local-Polynomial Degree And Bandwidth Search}{
This argument controls the recommended automatic local-polynomial NOMAD route, which jointly selects continuous polynomial degree and bandwidths when conditional-density bandwidths are computed inside \code{npconmode}.
}

  \item{nomad}{
    logical shortcut passed through to \code{\link{npcdensbw}} when
    bandwidths are computed inside \code{npconmode}. When \code{TRUE},
    the conditional-density bandwidth route fills any missing values
    among \code{regtype}, \code{search.engine}, \code{degree.select},
    \code{bernstein.basis}, \code{degree.min}, \code{degree.max},
    \code{degree.verify}, and \code{bwtype} with the recommended
    automatic local-polynomial degree-and-bandwidth NOMAD preset documented in \code{\link{npcdensbw}}.
    Additional NOMAD tuning arguments such as \code{nomad.nmulti} may
    also be supplied through \code{\dots}; \code{nmulti} remains the
    outer restart count while \code{nomad.nmulti} controls inner
    \code{crs::snomadr()} multistarts within each outer restart.
    Automatic NOMAD degree search requires at least one continuous
    explanatory variable. After fitting, inspect
    \code{fit$bws$nomad.shortcut} on the returned object \code{fit} to
    see the normalized shortcut metadata.
  }

\subsection{Additional Arguments}{
Further arguments are passed to the bandwidth-selection counterpart,
prediction/evaluation route, or plot route as appropriate.
}

  \item{\dots}{
    additional arguments supplied to \code{\link{npcdensbw}} when bandwidths
    are computed internally, or arguments needed to interpret a numeric
    \code{bws} vector. This is where bandwidth-selection controls such
    as \code{bwmethod}, \code{bwtype}, and \code{bwscaling},
    kernel/support controls such as \code{cxkertype}, \code{cykertype},
    \code{cxkerorder}, \code{cykerorder}, \code{cxkerbound}, and
    \code{cykerbound}, categorical kernel controls such as
    \code{uxkertype}, \code{uykertype}, \code{oxkertype}, and
    \code{oykertype}, search controls such as \code{nmulti} and
    \code{scale.factor.search.lower}, and local-polynomial/NOMAD
    controls such as \code{regtype}, \code{degree},
    \code{bernstein.basis}, \code{degree.select}, and
    \code{nomad.nmulti} are supplied. In \code{predict.conmode},
    additional arguments are passed to \code{npconmode} for evaluation
    with the stored bandwidth object; native \code{exdat} and
    \code{eydat} take precedence over \code{newdata}. In
    \code{plot.conmode}, additional arguments control the object-fed
    probability/effect display; common controls include
    \code{gradients}, \code{level}, \code{output}, \code{data_rug},
    \code{layout}, \code{legend}, \code{view}, \code{neval},
    \code{grid_control}, \code{perspective}, \code{plot.vars},
    \code{renderer}, and graphics arguments. The default
    \code{view="sample"} draws the stored fitted probabilities/effects at the
    fitted evaluation points. The \code{view="fixed"} route constructs
    one-dimensional object-fed slices over each conditioning variable, holding
    the remaining variables at their \code{uocquantile} median/mode values;
    \code{np_grid_control(xtrim=..., xq=...)} may be used to adjust the slice
    range and holdout quantiles. With \code{perspective=TRUE}, the same
    fixed-grid route draws a base-graphics probability surface for one selected
    response level over two continuous conditioning variables; use
    \code{plot.vars} to name the displayed pair when more than two continuous
    conditioning variables are present; \code{renderer="rgl"} requests the
    corresponding interactive surface when \pkg{rgl} is available. The public
    \code{output="both"} spelling is
    accepted as an alias for the historical \code{output="plot-data"}
    behavior. For class-probability plots, \code{errors="asymptotic"}
    adds probability-level asymptotic standard errors and intervals for the
    selected response level. These intervals target the smooth probability
    \eqn{\Pr(Y=\ell\mid X=x)}, not the discrete modal class label. Rows whose
    probabilities were materially repaired by simplex projection report
    \code{NA} standard errors rather than treating raw conditional-density
    standard errors as projected-probability standard errors. Bootstrap
    intervals are available for one-dimensional \code{view="sample"} and
    \code{view="fixed"} probability plots and use the stored bandwidth object
    without recomputing bandwidths. Surface asymptotic intervals are available
    as \code{output="data"} payload columns; surface band rendering and
    surface bootstrap intervals are intentionally deferred.
    See \code{\link{npcdensbw}} and \code{\link{plot.np}} for the complete
    bandwidth-selection and plotting argument surfaces.
  }

}
\value{
  \code{npconmode} returns a \code{conmode} object with the following
  components:

  \item{bws}{ the \code{conbandwidth} object used to compute the conditional
    mode }
  \item{conmode}{ a vector of type \code{factor} (or \code{ordered
      factor}) containing the conditional mode at each evaluation
    point}
  \item{condens}{ a vector of numeric type containing the modal density
    estimates at each evaluation point}
  \item{conderr}{ a vector of numeric type containing asymptotic
    standard errors for the modal density estimates at each evaluation point.
    If a row is materially properized by probability projection, this value is
    set to \code{NA} for that row rather than reporting a standard error for
    the unrepaired raw probability.}
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
  \item{probabilities}{ if requested, a matrix of fitted probabilities over
    the discrete response support }
  \item{probability.levels}{ if \code{probabilities} are stored, the response
    support corresponding to the probability matrix columns }
  \item{probability.errors}{ if \code{probabilities} are stored, a matrix of
    asymptotic standard errors for the raw selected-level conditional
    probability estimates, with rows set to \code{NA} where simplex projection
    materially repaired the fitted probability vector }
  \item{probability.repaired.rows}{ if \code{probabilities} are stored, a
    logical vector identifying rows whose fitted probability vector was
    materially repaired by proper probability projection }
  \item{probability.gradients}{ if \code{gradients=TRUE}, a matrix whose
    entries are class-probability gradients/effects for one selected response
    level, indexed by evaluation row and conditioning variable }
  \item{probability.gradient.level, probability.gradient.names,
    probability.gradient.info}{ metadata describing the response support,
    conditioning variables, and interpretation of
    \code{probability.gradients} }
  \item{proper.requested, proper.applied, proper.info}{ metadata describing
    whether discrete probability properization was requested, whether any row
    was materially projected, and diagnostics for non-negativity and unit-mass
    checks }
  \item{regtype, degree, nomad, search.engine}{ metadata mirrored from
    \code{bws} describing the conditional density regression type,
    local-polynomial degree, whether NOMAD automatic degree search was used,
    and the search engine when available }
  \item{degree.search, nomad.shortcut, nomad.time, powell.time}{ detailed
    search metadata mirrored from \code{bws} when automatic degree search or
    NOMAD/Powell refinement was used }

  The function \code{\link{predict}} may be used to extract conditional mode
  class predictions, while \code{\link{fitted}} extracts the conditional
  density estimates at the conditional mode from the resulting object. The
  function \code{\link{gradients}} extracts
  class-probability gradients/effects when \code{gradients=TRUE}. Also,
  \code{\link{summary}} and \code{\link{plot}} support \code{conmode}
  objects. For \code{plot.conmode}, first fit with
  \code{probabilities=TRUE} so the plot can remain object-fed. The default
  plot displays the fitted probability for the base/reference response level
  \code{levels(y)[1]}; use \code{plot(fit, level=...)} to select another
  outcome and \code{plot(fit, gradients=TRUE)} to display stored
  class-probability effects. The default \code{view="sample"} draws stored
  object-fed probability/effect payloads at the fitted evaluation points. Use
  \code{view="fixed"} and \code{neval} to draw object-fed one-dimensional
  slices over each conditioning variable, with other variables held at their
  median/mode values. Use \code{perspective=TRUE} to draw a base-graphics
  probability surface for one selected response level over two continuous
  conditioning variables, and \code{renderer="rgl"} for the corresponding
  interactive surface. Use \code{errors="asymptotic"} for probability-level
  standard errors and intervals in one-dimensional plots, or with
  \code{output="data"} for surface interval payloads. Bootstrap intervals and
  surface band rendering are not yet implemented for \code{plot.conmode}.
}
\section{Book And Method Pointers}{
  The conditional-mode target is
  \eqn{\arg\max_y \Pr(Y=y\mid X=x)} for a discrete response. In practice
  \code{npconmode} estimates the conditional probability of each
  response support point using \code{\link{npcdensbw}} /
  \code{\link{npcdens}}, optionally projects non-local-constant fitted
  probabilities onto the probability simplex, and then reports the
  modal support point.

  Setting \code{gradients=TRUE} stores class-probability effects for one
  response level. If \code{level} is omitted, the base/reference response level
  \code{levels(y)[1]} is used. These effects are useful for asking how the
  fitted probability of a selected class changes with each covariate; they are
  not gradients of the \code{arg max} classification rule itself.

  For book-length background, see Racine (2019), Chapter 4
  \emph{Conditional Probability Density and Cumulative Distribution
  Functions}, especially the binary and multinomial choice material,
  and Li and Racine (2007), Chapter 5 \emph{Conditional Density
  Estimation} together with Chapter 4 \emph{Kernel Estimation with Mixed
  Data}.
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
  \code{\link{np.kernels}}, \code{\link{np.options}},
  \code{\link{plot}}, \code{\link{npcdensbw}}.
}

\examples{
\dontrun{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we use the
# birthweight data taken from the MASS library, and compute a parametric
# logit model and a nonparametric conditional mode model. We then
# compare their confusion matrices and summary measures of
# classification ability.

library("MASS")
data("birthwt")

birthwt$low <- factor(birthwt$low)
birthwt$smoke <- factor(birthwt$smoke)
birthwt$race <- factor(birthwt$race)
birthwt$ht <- factor(birthwt$ht)
birthwt$ui <- factor(birthwt$ui)
birthwt$ftv <- ordered(birthwt$ftv)

with(birthwt, {

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

model.logit <- glm(low~smoke+
                   race+
                   ht+
                   ui+
                   ftv+
                   age+
                   lwt,
                   family=binomial(link=logit))

# Generate the confusion matrix and correct classification ratio

cm <- table(low, ifelse(fitted(model.logit)>0.5, 1, 0))
ccr <- sum(diag(cm))/sum(cm)

# Now do the same with a nonparametric model.  Note - this may take a
# few minutes depending on the speed of your computer...

bw <- npcdensbw(formula=low~smoke+
                race+
                ht+
                ui+
                ftv+
                age+
                lwt,
                data=birthwt)

model.np <- npconmode(bws=bw)

# Compare confusion matrices from the logit and nonparametric model

# Logit

cm
ccr

# Nonparametric
summary(model.np)

# Predict modal classes and fitted class probabilities at selected rows
new.birthwt <- birthwt[1:5, c("smoke", "race", "ht", "ui", "ftv", "age", "lwt")]
predict(model.np, newdata=new.birthwt)
predict(npconmode(bws=bw, probabilities=TRUE),
        newdata=new.birthwt, type="prob")

})

# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we use the
# birthweight data taken from the MASS library, and compute a parametric
# logit model and a nonparametric conditional mode model. We then
# compare their confusion matrices and summary measures of
# classification ability.

library("MASS")
data("birthwt")
with(birthwt, {

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
# computer...

bw <- npcdensbw(xdat=X, ydat=y)

model.np <- npconmode(bws=bw)

# Compare confusion matrices from the logit and nonparametric model

# Logit

cm
ccr

# Nonparametric
summary(model.np)

})

# EXAMPLE 3 (CLASS PROBABILITY EFFECTS): compute and plot
# class-probability gradients/effects for one selected response level.
# This example uses a small artificial sample so that it runs quickly.

set.seed(42)
n <- 100
x <- seq(-1, 1, length.out=n)
y <- factor(rbinom(n, 1, plogis(1.5*x)), levels=0:1)

model.effects <- npconmode(y~x,
                           regtype="ll",
                           bwmethod="cv.ls",
                           nmulti=1,
                           probabilities=TRUE,
                           gradients=TRUE)

gradients(model.effects)
plot(model.effects)
plot(model.effects, view="fixed", neval=25)
plot(model.effects, gradients=TRUE)
plot(model.effects, gradients=TRUE, view="fixed", neval=25)

## A two-continuous-predictor fit can be displayed as a probability surface.
## The default level is the base/reference response level.

z <- runif(n, -1, 1)
y2 <- factor(rbinom(n, 1, plogis(1.5*x - z)), levels=0:1)
model.surface <- npconmode(y2~x+z,
                           regtype="ll",
                           bwmethod="cv.ls",
                           nmulti=1,
                           probabilities=TRUE)
plot(model.surface, perspective=TRUE, neval=15)
} % enddontrun
}
\keyword{ nonparametric }% at least one, from doc/KEYWORDS
