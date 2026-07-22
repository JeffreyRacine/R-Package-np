## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## weval: optional evaluation data for the instruments
## xeval: optional evaluation data for the exogenous predictors

## ... optional arguments for npreg()

## This function returns a list with the following elements:

## phi: the IV estimator of phi(z) corresponding to the estimated
## derivative phi(z)
## phi.prime: the IV derivative estimator
## phi.mat: the matrix with colums phi_1, phi_2 etc. over all iterations
## phi.prime.mat: the matrix with colums phi'_1, phi'_2 etc. over all iterations
## num.iterations: number of iterations taken by Landweber-Fridman
## norm.stop: the stopping rule for each Landweber-Fridman iteration
## convergence: a character string indicating whether/why iteration terminated

## First, a series of functions for local polynomial kernel regression
## (regtype = "lp", with lc/ll as special degree-restricted cases)

## This function returns the weight matrix for a local polynomial

## supports mixed data types. It presumes that Y is in column 1. Basic
## error checking is undertaken. j.reg= strips off weights for mean
## (1), partials up to order p, and cross-partials. All partials and
## cross partials are wrt continuous regressors, and cross-partials
## require k > 1 and p > 1. Shrinking towards the local constant mean,
## first, and second partial derivatives is implemented for regions
## where the local polynomial estimator is ill-conditioned (sparse
## data, small h etc.).

.npregivderiv_select_stop_index <- function(norm.stop) {
  if(is.monotone.increasing(norm.stop)) {
    return(list(index=1L, monotone.failure=TRUE))
  }

  N <- 1L
  while(N < length(norm.stop) && norm.stop[N+1L] >= norm.stop[N]) N <- N+1L
  while(N < length(norm.stop) && norm.stop[N+1L] < norm.stop[N]) N <- N+1L

  list(index=N, monotone.failure=FALSE)
}

.npregivderiv_progress_stage <- function(state,
                                         label,
                                         iteration=NULL,
                                         show.now=FALSE,
                                         force=FALSE) {
  state$current_detail <- label
  if(isTRUE(show.now)) {
    return(.np_progress_show_now(state, done=iteration, detail=label))
  }

  .np_progress_step_at(
    state=state,
    now=.np_progress_now(),
    done=iteration,
    detail=label,
    force=force)
}

npregivderiv <- function(y, ...) UseMethod("npregivderiv")

.np_iv_deriv_adjoint_dots <- function(dots) {
  dot.names <- names(dots)
  if(is.null(dot.names)) return(dots)
  owned <- c("bandwidth.divide", "ukertype", "okertype")
  dots[!(!is.na(dot.names) & dot.names %in% owned)]
}

npregivderiv.default <- function(y,
                         z,
                         w,
                         x=NULL,
                         zeval=NULL,
                         weval=NULL,
                         xeval=NULL,
                         constant=0.5,
                         iterate.break=TRUE,
                         iterate.max=1000,
                         nmulti=NULL,
                         random.seed=42,
                         smooth.residuals=TRUE,
                         start.from=c("Eyz","EEywz"),
                         starting.values=NULL,
                         stop.on.increase=TRUE,
                         regtype=c("ll","lc","lp"),
                         degree=NULL,
                         nomad=FALSE,
                         ...) {

  ptm.start <- proc.time()
  regtype.missing <- missing(regtype)
  degree.missing <- missing(degree)
  nomad.missing <- missing(nomad)
  cl <- match.call()
  cl[[1L]] <- quote(npregivderiv)
  if(any(c("data", "newdata", "subset", "na.action") %in% names(cl)))
    stop("data, newdata, subset, and na.action require the formula interface",
         call. = FALSE)
  smoothing.spec <- .np_iv_resolve_deriv_smoothing(
    regtype = regtype,
    degree = degree,
    nomad = nomad,
    regtype.missing = regtype.missing,
    degree.missing = degree.missing,
    nomad.missing = nomad.missing
  )
  iterate.break <- .np_iv_validate_flag(iterate.break, "iterate.break")
  smooth.residuals <- .np_iv_validate_flag(smooth.residuals, "smooth.residuals")
  stop.on.increase <- .np_iv_validate_flag(stop.on.increase, "stop.on.increase")
  constant <- .np_iv_validate_number(constant, "constant", 0, 1, TRUE, TRUE)
  iterate.max <- .np_iv_validate_integer(iterate.max, "iterate.max", 2)
  random.seed <- .np_iv_validate_integer(random.seed, "random.seed", 0)
  iv.npreg <- function(...) {
    args <- list(...)
    if(is.null(args$bws)) args$random.seed <- random.seed
    do.call(npreg, c(args, .np_iv_deriv_stage_args(smoothing.spec, args$txdat)))
  }

  ## Basic error checking

  start.from <- match.arg(start.from)

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(!is.null(x) && NROW(y) != NROW(x)) stop("y and x have differing numbers of rows")
  if(!is.null(x))
    stop("npregivderiv currently supports one univariate endogenous z and does not support a separate exogenous x",
         call. = FALSE)

  ## The inverse problem is always solved on training rows. zeval is an
  ## output grid for the structural derivative only; instrument evaluation
  ## does not define a separate estimand.

  evaluation.requested <- !is.null(zeval)
  trainiseval <- !evaluation.requested
  if(is.null(x) && !is.null(xeval))
    stop("xeval requires exogenous training data x", call. = FALSE)
  if(!is.null(weval) &&
     !identical(as.data.frame(weval), as.data.frame(w))) {
    stop("weval cannot differ from the training instruments; use zeval to evaluate the structural derivative",
         call.=FALSE)
  }
  if(is.null(zeval)) zeval <- z
  weval <- w
  if(is.null(xeval)) xeval <- x

  if(NROW(zeval) < 1L)
    stop("zeval must contain at least one evaluation row", call.=FALSE)

  if(!is.null(starting.values) &&
     (!is.numeric(starting.values) || !is.null(dim(starting.values)) ||
      length(starting.values) != NROW(y) || anyNA(starting.values) ||
      any(!is.finite(starting.values)))) {
    stop("starting.values must be a finite numeric vector with one value per training row",
         call.=FALSE)
  }

  ## Need to determine how many x, w, z are numeric

  z <- data.frame(z)
  w <- data.frame(w)
  zeval <- data.frame(zeval)
  if(!is.null(x)) z <- data.frame(z,x)
  if(!is.null(xeval)) zeval <- data.frame(zeval,xeval)

  if(!is.null(nmulti)) nmulti <- npValidateNmulti(nmulti)
  nmulti.loop <- if(!is.null(nmulti)) nmulti else 1
  nmulti <- if(!is.null(nmulti)) nmulti else npDefaultNmulti(max(NCOL(z), NCOL(w)))

  z.numeric <- sapply(seq_len(NCOL(z)),function(i){is.numeric(z[,i])})
  num.z.numeric <- NCOL(as.data.frame(z[,z.numeric]))

  w.numeric <- sapply(seq_len(NCOL(w)),function(i){is.numeric(w[,i])})
  num.w.numeric <- NCOL(as.data.frame(w[,w.numeric]))

  ## Landweber-Fridman

  ## We begin the iteration computing phi.prime.0

  ## Note - here I am only treating the univariate case, so let's
  ## throw a stop with warning for now...

  if(NCOL(z) > 1)
    stop("npregivderiv currently supports one univariate endogenous z",
         call. = FALSE)
  if(NCOL(zeval) != 1L || !is.numeric(zeval[[1L]]))
    stop("zeval must contain one continuous numeric endogenous variable",
         call.=FALSE)
  if(!is.numeric(z[[1L]]))
    stop("npregivderiv requires one continuous numeric endogenous z",
         call.=FALSE)
  if(anyNA(z[[1L]]) || any(!is.finite(z[[1L]])) ||
     length(unique(z[[1L]])) < 2L) {
    stop("npregivderiv requires finite z values with at least two distinct support points",
         call.=FALSE)
  }
  if(anyNA(zeval[[1L]]) || any(!is.finite(zeval[[1L]])))
    stop("zeval must contain only finite values", call.=FALSE)
  if(!is.null(starting.values) && evaluation.requested &&
     !identical(as.numeric(zeval[[1L]]), as.numeric(z[[1L]]))) {
    stop("starting.values with a separate zeval grid is not supported; fit the training state first",
         call.=FALSE)
  }

  ## For all results we need the density function for Z and the
  ## survivor function for Z (1-CDF of Z)

  .np_progress_note("Preparing f_Z(z), S_Z(z)")

  ## Let's compute the bandwidth object for the unconditional
  ## density for the moment. Use the normal-reference rule for speed
  ## considerations, same smoothing for PDF and CDF.

  bw <- .np_progress_with_legacy_suppressed(npudensbw(dat=z, bwmethod="normal-reference"))
  model.fz <- .np_progress_with_legacy_suppressed(npudens(tdat=z, bws=bw$bw))
  f.z <- predict(model.fz, newdata=z)
  model.Sz <- .np_progress_with_legacy_suppressed(npudist(tdat=z, bws=bw$bw))
  S.z <- 1-predict(model.Sz, newdata=z)
  f.z.eval <- if(evaluation.requested) predict(model.fz, newdata=zeval) else NULL
  S.z.eval <- if(evaluation.requested) {
    1-predict(model.Sz, newdata=zeval)
  } else {
    NULL
  }

  progress <- .np_progress_begin("IV derivative", surface = "iv_solve")
  progress.finished <- FALSE
  on.exit({
    if(!isTRUE(progress.finished)) .np_progress_abort(progress)
  }, add = TRUE)

  progress <- .npregivderiv_progress_stage(
    progress,
    label="E[y|w]",
    show.now=TRUE)

  ## For stopping rule...

  model.E.y.w <- .np_progress_with_legacy_suppressed(iv.npreg(tydat=y,
                                                           txdat=w,
                                                           exdat=w,
                                                           nmulti=nmulti,
                                                           ...))
  E.y.w <- model.E.y.w$mean
  bw.E.y.w <- model.E.y.w$bws

  ## Potential alternative starting rule (consistent with
  ## npregiv). Here we start with E(Y|Z) rather than zero

  if(is.null(starting.values)) {

    progress <- .npregivderiv_progress_stage(
      progress,
      label=if(start.from=="Eyz") "d/dz E[y|z]" else "d/dz E[E[y|w]|z]",
      force=TRUE)

    model.phi.prime <- .np_progress_with_legacy_suppressed(iv.npreg(tydat=if(start.from=="Eyz") y else E.y.w,
                                                                 txdat=z,
                                                                 exdat=z,
                                                                 gradients=TRUE,
                                                                 nmulti=nmulti,
                                                                 ...))
    phi.prime <- model.phi.prime$grad[,1]
    bw.E.y.z <- model.phi.prime$bws

    phi.prime.eval <- if(evaluation.requested) {
      .np_progress_with_legacy_suppressed(iv.npreg(
        tydat=if(start.from=="Eyz") y else E.y.w,
        txdat=z,
        exdat=zeval,
        bws=bw.E.y.z,
        gradients=TRUE,
        nmulti=nmulti.loop,
        ...))$grad[,1]
    } else {
      NULL
    }

  } else {

    phi.prime <- starting.values
    phi.prime.eval <- if(evaluation.requested) starting.values else NULL
    bw.E.y.z <- NULL

  }

  ## Now we can compute phi.0 by integrating phi.prime.0 up to each
  ## sample realization (here we use the trapezoidal rule)

  ## NOTE - this presumes univariate z case... in general this would
  ## be a continuous variable's index

  # Pre-calculate components for integrate.trapezoidal
  z.val <- z[,1]
  n.z <- length(z.val)
  order.z <- order(z.val)
  z.sorted <- z.val[order.z]
  dz.sorted <- diff(z.sorted)
  cz.z <- dz.sorted[1]
  inv.order.z <- order(order.z)

  integrate.trapezoidal.internal <- function(y.val) {
      y.sorted <- y.val[order.z]
      dy.sorted <- diff(y.sorted)
      ca.z <- dy.sorted[1] / cz.z
      cb.z <- dy.sorted[n.z - 1] / cz.z
      cf.z <- cz.z^2 / 12 * (cb.z - ca.z)
      if (!is.finite(cf.z)) cf.z <- 0
      int.vec <- c(0, cumsum(dz.sorted * (y.sorted[-n.z] + y.sorted[-1]) / 2))
      int.vec <- int.vec - cf.z
      int.vec[inv.order.z]
  }

  phi <- integrate.trapezoidal.internal(phi.prime)

  ## In the definition of phi we have the integral minus the mean of
  ## the integral with respect to z, so subtract the mean here

  phi <- phi - mean(phi) + mean(y)

  starting.values.phi <- phi
  starting.values.phi.prime <- phi.prime

  ## Now we compute mu.0 (a residual of sorts)

  mu <- y - phi

  ## Now we repeat this entire process using mu = y = phi.0 rather than y

  if(smooth.residuals) {
    progress <- .npregivderiv_progress_stage(
      progress,
      label="E[y-phi_0(z)|w]",
      force=TRUE)

    ## Additional smoothing on top of the stopping rule required, but
    ## we have computed the stopping rule so reuse the bandwidth
    ## vector to be passed below. Here we compute the bandwidth
    ## optimal for the regression of mu on w.

    ## Next, we regress require \mu_{0,i} W using bws optimal for phi on w

    model.mu.w <- .np_progress_with_legacy_suppressed(iv.npreg(tydat=mu,
                                                            txdat=w,
                                                            exdat=w,
                                                            ...))
    predicted.E.mu.w <- model.mu.w$mean
    bw.mu.w <- model.mu.w$bws

  } else {

    progress <- .npregivderiv_progress_stage(
      progress,
      label="E[phi_0(z)|w]",
      force=TRUE)

    model.phi.w <- .np_progress_with_legacy_suppressed(iv.npreg(tydat=phi,
                                                             txdat=w,
                                                             exdat=w,
                                                             ...))
    E.phi.w <- model.phi.w$mean
    bw.mu.w <- model.phi.w$bws

    predicted.E.mu.w <- E.y.w - E.phi.w

  }

  ## Equation (14) requires an ordinary CDF-weighted average. Keep
  ## bandwidth.divide in ... available to the regression calls above, but do
  ## not allow it to change the normalization of this private adjoint.

  npksum.dots <- .np_iv_deriv_adjoint_dots(list(...))

  cdf.weighted.average.apply <- function(rhs, evaluation) {
    do.call(npksum,
            c(list(txdat=z,
                   exdat=evaluation,
                   tydat=as.matrix(rhs),
                   operator="integral",
                   ukertype="liracine",
                   okertype="liracine",
                   bws=bw$bw,
                   bandwidth.divide=TRUE),
              npksum.dots))$ksum/length(y)
  }

  adjoint.apply <- function(rhs, evaluation, density, survivor) {
    rhs.mean <- mean(rhs)
    cdf.average <- cdf.weighted.average.apply(rhs, evaluation)
    survivor.average <- rhs.mean - cdf.average
    (survivor.average-survivor*rhs.mean)/NZD(density)
  }

  ## Now we compute T^* applied to mu

  progress <- .npregivderiv_progress_stage(
    progress,
    label="T*{E[y-phi_0(z)|w]}",
    force=TRUE)
  T.star.mu <- adjoint.apply(predicted.E.mu.w, z, f.z, S.z)
  T.star.mu.eval <- if(evaluation.requested) {
    adjoint.apply(predicted.E.mu.w, zeval, f.z.eval, S.z.eval)
  } else {
    NULL
  }

  phi.mat <- matrix(NA, length(phi), iterate.max)
  phi.prime.mat <- matrix(NA, length(phi.prime), iterate.max)
  phi.prime.eval.mat <- if(evaluation.requested) {
    matrix(NA, length(phi.prime.eval), iterate.max)
  } else {
    NULL
  }
  norm.stop <- numeric(iterate.max)
  convergence <- "ITERATE_MAX"
  N.evaluated <- 0L

  ## Column N records the complete state after N derivative updates:
  ## phi_N, phi'_N, and the stopping rule evaluated at phi_N.

  for(N in seq_len(iterate.max)) {

    phi.prime <- phi.prime + constant*T.star.mu
    if(evaluation.requested)
      phi.prime.eval <- phi.prime.eval + constant*T.star.mu.eval

    ## NOTE - this presumes univariate z case... in general this would
    ## be a continuous variable's index

    phi <- integrate.trapezoidal.internal(phi.prime)

    ## In the definition of phi we have the integral minus the mean of
    ## the integral with respect to z, so subtract the mean here

    phi <- phi - mean(phi) + mean(y)

    ## Now we compute mu.0 (a residual of sorts)

    mu <- y - phi

    ## Now we repeat this entire process using mu = y = phi.0 rather than y

    if(smooth.residuals) {
      progress <- .npregivderiv_progress_stage(
        progress,
        label="E[y-phi(z)|w]",
        iteration=N)


      ## Additional smoothing on top of the stopping rule required, but
      ## we have computed the stopping rule so reuse the bandwidth
      ## vector to be passed below. Here we compute the bandwidth
      ## optimal for the regression of mu on w.

      ## Next, we regress require \mu_{0,i} W using bws optimal for phi on w

      model.mu.w <- .np_progress_with_legacy_suppressed(iv.npreg(tydat=mu,
                                                              txdat=w,
                                                              eydat=mu,
                                                              exdat=w,
                                                              bws=bw.mu.w,
                                                              nmulti=nmulti.loop,
                                                              ...))
      predicted.E.mu.w <- model.mu.w$mean
      bw.mu.w <- model.mu.w$bws

    } else {

      progress <- .npregivderiv_progress_stage(
        progress,
        label="E[phi(z)|w]",
        iteration=N)

      model.phi.w <- .np_progress_with_legacy_suppressed(iv.npreg(tydat=phi,
                                                               txdat=w,
                                                               eydat=phi,
                                                               exdat=w,
                                                               bws=bw.mu.w,
                                                               nmulti=nmulti.loop,
                                                               ...))
      E.phi.w <- model.phi.w$mean
      bw.mu.w <- model.phi.w$bws

      predicted.E.mu.w <- E.y.w - E.phi.w

    }

    norm.stop[N] <- N*sum(predicted.E.mu.w^2)/NZD_pos(sum(E.y.w^2))
    phi.mat[,N] <- phi
    phi.prime.mat[,N] <- phi.prime
    if(evaluation.requested) phi.prime.eval.mat[,N] <- phi.prime.eval
    N.evaluated <- N

    ## The number of iterations in LF is asymptotically equivalent to
    ## 1/alpha (where alpha is the regularization parameter in
    ## Tikhonov).  Plus the criterion function we use is increasing
    ## for very small number of iterations. So we need a threshold
    ## after which we can pretty much confidently say that the
    ## stopping criterion is decreasing.  In Darolles et al. (2011)
    ## \alpha ~ O(N^(-1/(min(beta,2)+2)), where beta is the so called
    ## qualification of your regularization method. Take the worst
    ## case in which beta = 0 and then the number of iterations is ~
    ## N^0.5. Note that derivative estimation seems to require more
    ## iterations hence the heuristic sqrt(N)

    should.break <- FALSE
    if(N > round(sqrt(nrow(z))) &&
       !is.monotone.increasing(norm.stop[seq_len(N)])) {

      ## If stopping rule criterion increases or we are below stopping
      ## tolerance then break

      if(stop.on.increase && norm.stop[N] > norm.stop[N-1L]) {
        convergence <- "STOP_ON_INCREASE"
        should.break <- isTRUE(iterate.break)
      }

    }

    if(should.break) break

    convergence <- "ITERATE_MAX"

    ## The current state's residual supplies the adjoint for update N+1.
    ## Do not compute it when state N is the terminal iterate.max state.

    if(N < iterate.max) {
      progress <- .npregivderiv_progress_stage(
        progress,
        label="T*{E[y-phi(z)|w]}",
        iteration=N)
      ## Equation (14) uses the same fitted conditional-residual vector in
      ## both empirical-adjoint terms.

      T.star.mu <- adjoint.apply(predicted.E.mu.w, z, f.z, S.z)
      if(evaluation.requested) {
        T.star.mu.eval <- adjoint.apply(predicted.E.mu.w,
                                        zeval,
                                        f.z.eval,
                                        S.z.eval)
      }
    }

  }

  ## Trim matrices and norm.stop to the actual number of iterations performed
  phi.mat <- phi.mat[, seq_len(N.evaluated), drop=FALSE]
  phi.prime.mat <- phi.prime.mat[, seq_len(N.evaluated), drop=FALSE]
  if(evaluation.requested)
    phi.prime.eval.mat <- phi.prime.eval.mat[, seq_len(N.evaluated), drop=FALSE]
  norm.stop <- norm.stop[seq_len(N.evaluated)]

  ## Extract minimum, and check for monotone increasing function and
  ## issue warning in that case. Otherwise allow for an increasing
  ## then decreasing (and potentially increasing thereafter) portion
  ## of the stopping function, ignore the initial increasing portion,
  ## and take the min from where the initial inflection point occurs
  ## to the length of norm.stop.

  stop.pick <- .npregivderiv_select_stop_index(norm.stop)
  N.selected <- stop.pick$index

  if(stop.pick$monotone.failure) {
    .np_warning("Stopping rule increases monotonically (consult model$norm.stop):\nThis could be the result of an inspired initial value (unlikely)\nNote: we suggest manually choosing phi.0 and restarting (e.g., instead set `start.from' to EEywz or provide a vector of starting values")
    convergence <- "FAILURE_MONOTONE_INCREASING"
  }

  phi <- phi.mat[,N.selected]
  phi.prime <- phi.prime.mat[,N.selected]
  phi.prime.eval <- if(evaluation.requested) {
    phi.prime.eval.mat[,N.selected]
  } else {
    NULL
  }
  
  progress <- .np_progress_end(progress, detail=progress$current_detail)
  progress.finished <- TRUE

  if(N.evaluated == iterate.max) .np_warning(" iterate.max reached: increase iterate.max or inspect norm.stop vector")

  ret <- list(phi=phi,
              phi.prime=phi.prime,
              phi.prime.eval=phi.prime.eval,
              phi.mat=phi.mat,
              phi.prime.mat=phi.prime.mat,
              phi.prime.eval.mat=phi.prime.eval.mat,
              num.iterations=N.selected,
              norm.stop=norm.stop,
              convergence=convergence,
              starting.values.phi=starting.values.phi,
              starting.values.phi.prime=starting.values.phi.prime,
              bw.E.y.w=bw.E.y.w,
              bw.E.y.z=bw.E.y.z,
              call=cl,
              y=y,
              z=z,
              w=w,
              x=x,
              zeval=zeval,
              weval=weval,
              xeval=xeval,
              bws=list(E.y.w=.np_iv_bw_payload(bw.E.y.w),
                       E.y.z=.np_iv_bw_payload(bw.E.y.z),
                       residual.w=.np_iv_bw_payload(bw.mu.w)),
              smoothing.spec=smoothing.spec,
              stage.specs=list(
                E.y.w=.np_iv_stage_spec("E.y.w", smoothing.spec, w),
                E.y.z=.np_iv_stage_spec("E.y.z", smoothing.spec, z),
                residual.w=.np_iv_stage_spec("residual.w", smoothing.spec, w)),
              trainiseval=trainiseval,
              nmulti=nmulti,
              ptm=proc.time() - ptm.start)
  class(ret) <- "npregivderiv"
  return(ret)

}

print.npregivderiv <- function(x, ...) {
  print(summary.npregivderiv(x, ...))
  invisible(x)
}

summary.npregivderiv <- function(object, ...) {
  ans <- .np_iv_summary_payload(object, derivative = TRUE)
  class(ans) <- "summary.npregivderiv"
  ans
}

plot.npregivderiv <- function(x,
                              gradients = TRUE,
                              data_overlay = TRUE,
                              data_rug = FALSE,
                              ...) {
  matched.call <- match.call(expand.dots = FALSE)
  .np_plot_validate_npregiv_call(
    sys.call(),
    method_args = c("gradients", "data_overlay", "data_rug"),
    context = "plot.npregivderiv",
    dots_call = matched.call[["..."]]
  )
  spec <- .np_iv_plot_spec(
    object = x,
    family = "npregivderiv",
    gradients = gradients,
    data_overlay = data_overlay,
    data_rug = data_rug,
    dots = list(...)
  )
  .np_iv_plot_render(spec)
}
