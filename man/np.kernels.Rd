\name{np.kernels}
\alias{np.kernels}
\title{Kernel Functions Used In \pkg{npRmpi}}
\description{
  Summary of continuous, unordered-categorical, and ordered-categorical
  kernels used by \pkg{npRmpi} (including higher-order continuous kernels and
  compact-support variants used in C-level code paths).
}
\details{
  For MPI startup/performance guidance (including message-passing tradeoffs and the manual-broadcast template), see \code{\link{npRmpi.init}} details and \code{inst/Rprofile}.
  Documentation guide: see \code{\link{np.options}} for global options For interactive and cluster batch workflows, see \code{\link{npRmpi.init}}.
  and \code{\link{plot}} for plotting options.

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

  The beta associated-kernel family is available for continuous data with
  finite support. For an
  evaluation point \eqn{x}, observation \eqn{x_i}, fixed support
  \eqn{[a,b]}, \eqn{L=b-a}, and public metric bandwidth \eqn{h}, define
  \deqn{t=\frac{x-a}{L},\quad z_i=\frac{x_i-a}{L},\quad
        \tau=\left(\frac{L}{h}\right)^2.}
  The implemented PDF contribution is Chen's simple beta kernel
  \deqn{K_{\mathit{Beta},2}(x_i;x,h,a,b)=
    \frac{1}{L}\frac{z_i^{t\tau}(1-z_i)^{(1-t)\tau}}
    {B(1+t\tau,1+(1-t)\tau)}.}
  Thus \code{ckertype="beta"} uses the evaluation-centered associated-kernel
  orientation: the beta shape parameters depend on \eqn{x}, and the density
  is evaluated at \eqn{x_i}.  The bandwidth convention is deliberately
  metric-compatible with other continuous kernels: Chen's smoothing
  parameter is \eqn{\rho=(h/L)^2}.  The factor \eqn{1/L} is already included,
  so callers must not apply either an additional \eqn{1/h} factor or the
  generic bounded-kernel normalization below.

  Orders 4, 6, and 8 use the additive higher-order construction of Igarashi
  and Kakizawa (2020). For public order \eqn{2m},
  \deqn{K_{\mathit{Beta},2m}(x_i;x,h)=
    \sum_{s=1}^m c_{m,s}K_{\mathit{Beta},2}(x_i;x,h\sqrt{s}),
    \qquad c_{m,s}=(-1)^{s+1}{m\choose s}.}
  The coefficient vectors for orders 2, 4, 6, and 8 are respectively
  \eqn{(1)}, \eqn{(2,-1)}, \eqn{(3,-3,1)}, and
  \eqn{(4,-6,4,-1)}. The construction preserves unit mass and the original
  support while deliberately producing signed kernels above order 2.

  In the present implementation orders 2, 4, 6, and 8 are available through
  \code{npksum}, \code{npudens}, \code{npudist}, local-constant
  \code{npreg}, and the local-constant conditional density/distribution
  families. Manual and cross-validated bandwidth objects support
  \code{bwtype="fixed"}, \code{"generalized_nn"}, and
  \code{"adaptive_nn"}. Bounds may be explicit
  (\code{ckerbound="fixed"}) or the exact training-sample range
  (\code{ckerbound="range"}); all active beta bounds must be finite.
  The beta routes currently require continuous variables only. Regression and
  conditional estimators additionally require \code{regtype="lc"}; if either
  conditional side uses beta, all X and Y coordinates must be continuous.
  \code{npksum} requires \code{bandwidth.divide=FALSE} and
  \code{kernel.pow=1}. Normal, integral, convolution, and one-coordinate
  PDF-derivative operators are supported;
  \code{permutation.operator="derivative"} returns every continuous first
  derivative. Score, OCG, categorical-product, higher/mixed derivative, and
  local-polynomial combinations fail explicitly rather than falling through
  a legacy kernel code.

  For scale component \eqn{s}, put \eqn{\tau_s=\tau/s},
  \eqn{\alpha_s=1+t\tau_s}, and \eqn{\beta_s=1+(1-t)\tau_s}. At an interior
  observation \eqn{0<z_i<1}, its evaluation-coordinate derivative is
  \deqn{\frac{\partial K_s}{\partial x}=\frac{\tau_s}{L}K_s
    \left\{\log z_i-\log(1-z_i)-\psi(\alpha_s)+\psi(\beta_s)\right\},}
  where \eqn{\psi} is the digamma function. Higher-order derivatives use the
  same signed multiscale coefficients as the PDF. Bandwidths are held fixed
  while differentiating: this includes the evaluation bandwidth selected by a
  generalized-NN rule and each observation bandwidth selected by an
  adaptive-NN rule. These are derivatives of the fitted smoother conditional
  on its realized bandwidths, not derivatives of the NN selection map.

  An observation exactly equal to the lower (upper) support endpoint makes
  its matching scalar PDF contribution jump downward (upward) when the target
  crosses that endpoint. The scalar derivative is consequently reported as
  \eqn{-\infty} at a matching lower endpoint and \eqn{+\infty} at a matching
  upper endpoint; nonmatching endpoint contributions have derivative zero.
  Local-constant regression and conditional estimators cancel numerator and
  denominator jumps before classifying their gradient, so a normalized ratio
  may still have a finite derivative (and a constant response or constant Y
  factor gives exactly zero). A genuine ratio jump is returned as a signed
  infinity. Numerically indeterminate cancellation is returned as \code{NA}
  with a warning. Derivative hat matrices cannot encode response-specific
  cancellation and therefore mark an affected endpoint row \code{NA}.

  For nearest-neighbor modes, the
  existing observation-support nearest-neighbor engine supplies the
  evaluation-specific metric bandwidth \eqn{h(x)} or observation-specific
  metric bandwidth \eqn{h_i}, after which the same beta concentration
  mapping is used.
  For overlap, generalized bandwidths are evaluated at both centers;
  adaptive bandwidths remain attached to training centers, while a genuinely
  new evaluation center receives the established generalized-NN auxiliary
  bandwidth. PDF, CDF, and overlap operators support both NN modes. See
  \code{\link{npksum}} for the exact kernel-sum normalization and
  \code{\link{npudens}} and \code{\link{npudist}} for estimator-specific
  standard errors and examples.

  The implemented beta CDF factor intentionally uses the complementary
  observation-centered orientation.  With the notation above,
  \deqn{H_{\mathit{Beta},2}(x;x_i,h,a,b)=
    I_t(1+z_i\tau,1+(1-z_i)\tau),}
  where \eqn{I_t(\alpha,\beta)} is the regularized incomplete beta function.
  It is evaluated as an analytic incomplete-beta probability, with a
  symmetry transformation in the upper tail. Thus every order-2
  contribution is in \eqn{[0,1]}, and with fixed or adaptive bandwidth is
  nondecreasing in \eqn{x}; it is exactly
  zero at \eqn{a}, exactly one at \eqn{b}, and tends to
  \eqn{(x-a)/(b-a)} as \eqn{h\to\infty}. Multivariate CDFs use products of
  these factors; \code{npksum} may also mix beta PDF and CDF factors across
  dimensions.

  This observation-centered beta CDF is not the integral, with respect to
  \eqn{x}, of the evaluation-centered Chen PDF contribution. It is also not
  the target-centered bounded-CDF estimator of Mansouri, Rastin, and Mombeni
  (2024). The orientation was selected because it preserves the essential
  empirical-CDF contribution limit as \eqn{h\to0}, exact endpoint behavior,
  monotonicity for order 2, and a direct analytic implementation suitable for
  the existing kernel-sum architecture. Orders 4, 6, and 8 apply the same
  signed multiscale coefficients to these incomplete-beta factors. The
  Igarashi--Kakizawa density theorem does not establish a corresponding CDF
  bias order; the public order identifies the construction. Signed
  higher-order beta CDFs may overshoot \eqn{[0,1]} or be locally nonmonotone,
  just as existing higher-order distribution kernels may. Values are not
  silently clipped or monotonized.

  With \code{bwtype="generalized_nn"}, the CDF bandwidth changes with the
  target. Every returned order-2 contribution remains an incomplete-beta
  probability in \eqn{[0,1]} with exact support endpoints, but the resulting
  curve need not be monotone. A deterministic boundary-heavy comparison was
  monotone across the recorded chi-square(1), beta-boundary, U-shaped, and
  uniform designs, while a seeded 300-design adversarial search found local
  decreases as large as approximately \eqn{8.6\times 10^{-5}}. The
  implementation neither clips nor monotonizes those values. Inversion-based
  uses therefore remain outside this generalized-NN CDF contract.

  The beta derivative operator differentiates the PDF contribution with
  respect to its evaluation coordinate. It does not differentiate the
  observation-centered CDF factor with respect to Y, and it is not a
  bandwidth score. As with existing higher-order kernels and local-polynomial
  derivatives, signed kernels can yield non-proper density or distribution
  estimates; derivative output is not clipped, monotonized, or otherwise
  repaired.

  The beta \code{"convolution"} operator is the analytic overlap of two
  associated beta density weights in their common observation coordinate. If
  the two centers and bandwidths give shapes
  \eqn{(\alpha_1,\beta_1)} and \eqn{(\alpha_2,\beta_2)}, then
  \deqn{\mathcal O=
    \frac{1}{L}\frac{B(\alpha_1+\alpha_2-1,
                             \beta_1+\beta_2-1)}
                        {B(\alpha_1,\beta_1)B(\alpha_2,\beta_2)}.}
  The diagonal is the exact location- and bandwidth-dependent roughness.
  This is an \eqn{L^2} overlap, is symmetric, has units \eqn{L^{-1}}, and is
  evaluated analytically in C with stable log-gamma ratios. It is not a
  translation convolution and is not Chen's target-coordinate integrated
  square. In particular, bounded \code{cv.ls} density and conditional-density
  objectives continue to integrate the actual bounded estimator by their
  existing quadrature route; they do not substitute this public overlap.

  Bounded continuous-kernel normalization (\code{ckerbound} and, for
  conditional objects, \code{cxkerbound}/\code{cykerbound}) reuses the
  selected continuous kernel and renormalizes it on the declared support.
  For a base kernel \eqn{K} and support \eqn{[a,b]}, the bounded kernel is
  \deqn{K_{[a,b]}(u;x,h)=\frac{K(u)}{\int_{(a-x)/h}^{(b-x)/h}K(t)dt}}
  with \eqn{u=(x_i-x)/h}. Option \code{ckerbound="range"} uses
  sample bounds for \eqn{a,b}; \code{ckerbound="fixed"} uses user-supplied
  bounds via \code{ckerlb}/\code{ckerub} (or the corresponding
  \code{cx*}/\code{cy*} arguments). Infinite bounds recover the unbounded
  kernel. This support-normalization strategy follows the same Racine-Li-Yan
  finite-support normalization principle and is useful when data exhibit
  non-negligible probability mass near boundaries.

  Computationally, let
  \eqn{H_K(r)=\int_0^r K(t)dt} and
  \eqn{S_K(r)=H_K(r)/r}, with the continuous value
  \eqn{S_K(0)=K(0)}. For finite bounds the implementation evaluates the scaled
  normalization mass directly as
  \deqn{hD(x,h)=(x-a)S_K((x-a)/h)+(b-x)S_K((b-x)/h).}
  Each supported kernel uses its analytic centered primitive. This avoids
  subtracting nearly equal whole-kernel CDF values and avoids constructing a
  normalization mass that vanishes numerically as \eqn{h} grows. Consequently,
  on finite \eqn{[a,b]} the bounded density kernel converges to the uniform
  density \eqn{1/(b-a)} as \eqn{h\to\infty}. One-sided infinite bounds use the
  corresponding exact half-mass term; two infinite bounds retain the ordinary
  unbounded kernel path.

  For distribution estimation, finite bounds use the distinct
  observation-centered truncated-kernel CDF
  \deqn{B_i(y;h)=\frac{G((y-Y_i)/h)-G((a-Y_i)/h)}
  {G((b-Y_i)/h)-G((a-Y_i)/h)},\quad y\in[a,b].}
  Thus \eqn{B_i(a;h)=0}, \eqn{B_i(b;h)=1}, and
  \eqn{B_i(y;h)\to(y-a)/(b-a)} as \eqn{h\to\infty}. The implementation
  evaluates centered analytic interval masses rather than subtracting
  nearly equal whole-kernel CDF values. This distribution operator is not the
  integral of the evaluation-centered bounded density kernel above. For
  nonnegative second-order kernels each contribution is monotone and lies in
  \eqn{[0,1]}; higher-order signed kernels retain their usual possible interior
  nonmonotonicity or overshoot.

  Typical bounded-kernel calls:
  \preformatted{
  ## Unconditional density on [0,1]
  bw <- npudensbw(dat=data.frame(x),
                  ckertype="gaussian",
                  ckerbound="fixed", ckerlb=0, ckerub=1)

  ## Regression with automatic sample-range bounds
  bw <- npregbw(xdat=data.frame(x), ydat=y, ckerbound="range")

  ## Conditional density with separate x/y support controls
  bw <- npcdensbw(xdat=data.frame(x), ydat=data.frame(y),
                  cxkerbound="fixed", cxkerlb=0, cxkerub=1,
                  cykerbound="range")
  }

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

  Chen, S. X. (1999),
  \dQuote{Beta kernel estimators for density functions,}
  \emph{Computational Statistics & Data Analysis}, \bold{31}(2), 131--145.
  doi:10.1016/S0167-9473(99)00010-9

  Jones, M. C. and Henderson, D. A. (2007),
  \dQuote{Kernel-type density estimation on the unit interval,}
  \emph{Biometrika}, \bold{94}(4), 977--984.
  doi:10.1093/biomet/asm068

  Igarashi, G. and Kakizawa, Y. (2020),
  \dQuote{Higher-order bias corrections for asymmetric kernel density
  estimators,}
  \emph{Journal of Nonparametric Statistics}, \bold{32}(3), 716--752.
  doi:10.1080/10485252.2020.1770754

  Mansouri, B., Rastin, A., and Allah Mombeni, H. (2024),
  \dQuote{Beta Kernel Estimator for a Cumulative Distribution Function with
  Bounded Support,}
  \emph{Journal of Sciences, Islamic Republic of Iran}, \bold{34}(4), 349--361
  (volume dated 2023; accepted and DOI registered in 2024).
  doi:10.22059/jsciences.2024.367933.1007836

  Hall, P., Racine, J. S., and Li, Q. (2004),
  \dQuote{Cross-validation and the estimation of conditional probability
  densities,}
  \emph{Journal of the American Statistical Association}, \bold{99}, 1015--1026.
}
\seealso{
  \code{\link{np.options}}, \code{\link{plot}}
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
