\name{np.kernels}
\alias{np.kernels}
\title{Kernel Functions Used In \pkg{npRmpi}}
\description{
  Summary of continuous, unordered-categorical, and ordered-categorical
  kernels used by \pkg{npRmpi} (including higher-order continuous kernels and
  compact-support variants used in C-level code paths).
}
\details{
  Documentation guide: see \code{\link{np.options}} for global options
  and \code{\link{npplot}} for plotting options.

  Kernel option names used in \pkg{npRmpi}:

  \itemize{
    \item Continuous kernels:
    \code{ckertype} (and \code{ckerorder}, \code{ckerbound} where applicable).
    \item Unordered kernels:
    \code{ukertype}.
    \item Ordered kernels:
    \code{okertype}.
    \item Conditional density/distribution bandwidth objects split kernel
    choices by response and regressor blocks:
    \code{cykertype}/\code{cxkertype},
    \code{uykertype}/\code{uxkertype},
    \code{oykertype}/\code{oxkertype}
    (with matching order/bound options for continuous kernels).
  }

  Let \eqn{u = (x_i-x)/h} for continuous variables.

  Continuous kernels (called via \code{ckertype}):
  \deqn{K_{G,2}(u)=\phi(u)}
  \deqn{K_{G,4}(u)=\left(\frac{3}{2}-\frac{1}{2}u^2\right)\phi(u)}
  \deqn{K_{G,6}(u)=\left(\frac{15}{8}-\frac{5}{4}u^2+\frac{1}{8}u^4\right)\phi(u)}
  \deqn{K_{G,8}(u)=\left(\frac{35}{16}-\frac{35}{16}u^2+\frac{7}{16}u^4-\frac{1}{48}u^6\right)\phi(u)}
  where \eqn{\phi(u)} is the standard normal density.

  \code{ckertype="gaussian"} with \code{ckerorder=2,4,6,8}.

  The compact-support Epanechnikov-family kernels implemented in C use
  support \eqn{|u|<\sqrt{5}}:
  \deqn{K_{E,2}(u)=\frac{3}{4\sqrt{5}}\left(1-\frac{u^2}{5}\right)\mathbf{1}(|u|<\sqrt{5})}
  \deqn{K_{E,4}(u)=0.008385254916(-15+7u^2)(-5+u^2)\mathbf{1}(u^2<5)}
  \deqn{K_{E,6}(u)=0.33541019662496845446\left(2.734375-3.28125u^2+0.721875u^4\right)\left(1-0.2u^2\right)\mathbf{1}(u^2<5)}
  \deqn{K_{E,8}(u)=0.33541019662496845446\left(3.5888671875-7.8955078125u^2+4.1056640625u^4-0.5865234375u^6\right)\left(1-0.2u^2\right)\mathbf{1}(u^2<5)}

  \code{ckertype="epanechnikov"} with \code{ckerorder=2,4,6,8}.

  Uniform (rectangular) kernel:
  \deqn{K_U(u)=\frac{1}{2}\mathbf{1}(|u|<1)}
  via \code{ckertype="uniform"} (order ignored).

  Truncated-Gaussian (second-order) kernel via
  \code{ckertype="truncated gaussian"}:
  \deqn{K_{TG,2}(u)=\left[\alpha\phi(u)-c_0\right]\mathbf{1}(|u|<b)}
  with defaults \eqn{b=3} and internal constants calibrated in C.

  Unordered-categorical kernels (called via \code{ukertype}; for category
  count \eqn{c}):
  \deqn{L_{AA}(x_i,x;\lambda)=\mathbf{1}(x_i=x)(1-\lambda)+\mathbf{1}(x_i\neq x)\frac{\lambda}{c-1}}
  (Aitchison-Aitken)
  via \code{ukertype="aitchisonaitken"}.

  \deqn{L_{LR,u}(x_i,x;\lambda)=\mathbf{1}(x_i=x)+\mathbf{1}(x_i\neq x)\lambda}
  (Li-Racine unordered kernel)
  via \code{ukertype="liracine"}.

  Ordered-categorical kernels (called via \code{okertype}):
  \deqn{L_{WvR}(x_i,x;\lambda)=
    \begin{cases}
      1-\lambda, & x_i=x\\
      \frac{1-\lambda}{2}\lambda^{|x_i-x|}, & x_i\neq x
    \end{cases}}
  (Wang-van Ryzin)
  via \code{okertype="wangvanryzin"}.

  \deqn{L_{LR,o}(x_i,x;\lambda)=\lambda^{|x_i-x|}}
  (Li-Racine ordered kernel)
  via \code{okertype="liracine"}.

  \deqn{L_{NLR,o}(x_i,x;\lambda)=\lambda^{|x_i-x|}\frac{1-\lambda}{1+\lambda}}
  (normalized Li-Racine ordered kernel; used internally)

  \deqn{L_{RLY}(x_i,x;\lambda)=\frac{\lambda^{|x_i-x|}}{\sum_{z\in\mathcal{S}(x)}\lambda^{|x_i-z|}}}
  (Racine-Li-Yan ordered kernel, normalized on support \eqn{\mathcal{S}(x)}).
  exposed as \code{okertype="racineliyan"}.

  These univariate kernels are combined as generalized product kernels over
  mixed data types in the estimators and cross-validation criteria.
}
\references{
  Aitchison, J. and Aitken, C. G. G. (1976),
  \dQuote{Multivariate binary discrimination by the kernel method,}
  \emph{Biometrika}, \bold{63}, 413--420.

  Wang, M. C. and Van Ryzin, J. (1981),
  \dQuote{A class of smooth estimators for discrete distributions,}
  \emph{Biometrika}, \bold{68}, 301--309.

  Li, Q. and Racine, J. S. (2007),
  \emph{Nonparametric Econometrics: Theory and Practice},
  Princeton University Press.

  Racine, J. S. and Li, Q. (2004),
  \dQuote{Nonparametric estimation of regression functions with both
  categorical and continuous data,}
  \emph{Journal of Econometrics}, \bold{119}, 99--130.

  Racine, J. S., Li, Q., and Yan, K. X. (2020),
  \dQuote{Kernel Smoothed Probability Mass Functions for Ordered Datatypes,}
  \emph{Journal of Nonparametric Statistics}, \bold{32}(3), 563--586.
  doi:10.1080/10485252.2020.1759595

  Hall, P., Racine, J. S., and Li, Q. (2004),
  \dQuote{Cross-validation and the estimation of conditional probability
  densities,}
  \emph{Journal of the American Statistical Association}, \bold{99}, 1015--1026.
}
\seealso{
  \code{\link{np.options}}, \code{\link{npplot}}
  \code{\link{npregbw}},
  \code{\link{npudensbw}},
  \code{\link{npudistbw}},
  \code{\link{npcdensbw}},
  \code{\link{npcdistbw}},
  \code{\link{npksum}},
  \code{\link{np.options}}.
}
\keyword{nonparametric}
\keyword{smooth}
