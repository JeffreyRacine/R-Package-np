nplsqreg <-
  function(bws, ...) {
    args <- list(...)

    if (!missing(bws)) {
      if (isa(bws, "lsqregressionbandwidth"))
        UseMethod("nplsqreg", bws)
      if (inherits(bws, "formula") && is.null(args$txdat))
        UseMethod("nplsqreg", bws)
      if (is.recursive(bws)) {
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("nplsqreg", bws$formula)
        else if (!is.call(bws))
          UseMethod("nplsqreg", bws)
        else
          UseMethod("nplsqreg", NULL)
      } else {
        UseMethod("nplsqreg", NULL)
      }
    } else {
      UseMethod("nplsqreg", NULL)
    }
  }

nplsqregbw <-
  function(bws, ...) {
    args <- list(...)

    if (!missing(bws)) {
      if (isa(bws, "lsqregressionbandwidth"))
        UseMethod("nplsqregbw", bws)
      if (inherits(bws, "formula") && is.null(args$xdat))
        UseMethod("nplsqregbw", bws)
      if (is.recursive(bws) && !is.call(bws))
        UseMethod("nplsqregbw", bws)
      else
        UseMethod("nplsqregbw", NULL)
    } else {
      UseMethod("nplsqregbw", NULL)
    }
  }

.nplsqreg_validate_tau <- function(tau) {
  if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) ||
      !is.finite(tau) || tau <= 0 || tau >= 1)
    stop("'tau' must be a single finite numeric value in (0, 1)",
         call. = FALSE)
  as.numeric(tau)
}

.nplsqreg_check_loss <- function(u, tau) {
  mean(u * (tau - (u < 0)), na.rm = TRUE)
}

.nplsqreg_validate_scale <- function(scale, n) {
  if (!is.numeric(scale) || length(scale) != n ||
      any(!is.finite(scale)) || any(scale <= 0))
    stop("'scale' must be a finite positive numeric vector matching the training data length",
         call. = FALSE)
  as.numeric(scale)
}

.nplsqreg_scale_pilot <- function(xdat, ydat, dots) {
  mean.args <- c(list(txdat = xdat, tydat = ydat), dots)
  mean.fit <- do.call(npreg, mean.args)
  res <- as.numeric(ydat) - as.numeric(fitted(mean.fit))
  scale.args <- c(list(txdat = xdat, tydat = res^2), dots)
  scale.fit <- do.call(npreg, scale.args)
  scale <- sqrt(abs(as.numeric(fitted(scale.fit))))
  floor <- sqrt(.Machine$double.eps)
  scale[!is.finite(scale) | scale <= floor] <- floor
  list(scale = scale, mean.fit = mean.fit, scale.fit = scale.fit)
}

.nplsqreg_rebuild_rbandwidth <- function(template, bw, ydat,
                                         fval = NA, num.feval = NA,
                                         bandwidth.compute = FALSE) {
  out <- rbandwidth(
    bw = as.numeric(bw),
    regtype = template$regtype,
    basis = template$basis,
    degree = template$degree,
    bernstein.basis = template$bernstein.basis,
    bwmethod = if (is.null(template$method)) "cv.ls" else template$method,
    bwscaling = template$scaling,
    bwtype = template$type,
    ckertype = template$ckertype,
    ckerorder = template$ckerorder,
    ckerbound = template$ckerbound,
    ckerlb = template$ckerlb,
    ckerub = template$ckerub,
    ukertype = template$ukertype,
    okertype = template$okertype,
    fval = fval,
    num.feval = num.feval,
    nobs = template$nobs,
    xdati = template$xdati,
    ydati = untangle(data.frame(ydat)),
    xnames = template$xnames,
    ynames = template$ynames,
    sfactor = as.numeric(bw),
    bandwidth = as.numeric(bw),
    rows.omit = template$rows.omit,
    nconfac = template$nconfac,
    ncatfac = template$ncatfac,
    sdev = template$sdev,
    bandwidth.compute = bandwidth.compute,
    timing = template$timing,
    total.time = template$total.time)
  out
}

.nplsqreg_bw_bounds <- function(template, xdat, start) {
  lower <- rep.int(0, length(start))
  upper <- rep.int(Inf, length(start))

  if (any(template$icon)) {
    con.idx <- which(template$icon)
    for (j in con.idx) {
      xj <- as.numeric(xdat[[j]])
      sx <- stats::sd(xj, na.rm = TRUE)
      rx <- diff(range(xj, na.rm = TRUE))
      scale <- max(abs(start[j]), sx, rx, 1)
      lower[j] <- max(.Machine$double.eps, scale * 1e-5)
      upper[j] <- scale * 20
    }
  }

  if (any(template$iuno | template$iord)) {
    cat.idx <- which(template$iuno | template$iord)
    cat.max <- template$sumNum$x
    for (j in cat.idx) {
      uj <- if (length(cat.max) >= j && is.finite(cat.max[j])) cat.max[j] else 1
      upper[j] <- max(uj, abs(start[j]), .Machine$double.eps)
    }
  }

  lower <- pmin(lower, upper - sqrt(.Machine$double.eps))
  start <- pmin(pmax(start, lower + sqrt(.Machine$double.eps)),
                upper - sqrt(.Machine$double.eps))
  list(start = start, lower = lower, upper = upper)
}

.nplsqreg_objective_factory <- function(xdat, ydat, scale, tau, template,
                                        delta.bounds) {
  invalid.count <- 0L
  f.history <- numeric()
  p.history <- list()

  objective <- function(par) {
    bw <- par[-length(par)]
    delta <- par[length(par)]
    qdat <- as.numeric(ydat) + scale * stats::qnorm(delta)
    rbw <- try(.nplsqreg_rebuild_rbandwidth(template, bw, qdat),
               silent = TRUE)
    if (inherits(rbw, "try-error")) {
      invalid.count <<- invalid.count + 1L
      return(.Machine$double.xmax^0.5)
    }
    qloo <- try(npreghat(bws = rbw, txdat = xdat, y = qdat,
                         output = "apply", leave.one.out = TRUE),
                silent = TRUE)
    if (inherits(qloo, "try-error") || any(!is.finite(qloo))) {
      invalid.count <<- invalid.count + 1L
      return(.Machine$double.xmax^0.5)
    }
    value <- .nplsqreg_check_loss(as.numeric(ydat) - as.numeric(qloo), tau)
    f.history <<- c(f.history, value)
    p.history[[length(p.history) + 1L]] <<- par
    value
  }

  environment(objective)$state <- function() {
    list(invalid.count = invalid.count,
         f.history = f.history,
         p.history = p.history)
  }
  objective
}

nplsqregbw.formula <-
  function(bws, data = NULL, tau = 0.5, ...) {
    tt <- terms(bws)
    mc <- match.call(expand.dots = FALSE)
    m <- match(c("bws", "data", "subset", "na.action"),
               names(mc), nomatch = 0)
    tmf <- mc[c(1, m)]
    if ("bws" %in% names(tmf))
      names(tmf)[names(tmf) == "bws"] <- "formula"
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    mf <- do.call(stats::model.frame, mf.args, envir = environment(tt))
    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"), "term.labels"), drop = FALSE]
    out <- nplsqregbw(xdat = xdat, ydat = ydat, tau = tau, ...)
    out$formula <- bws
    out$call <- match.call(expand.dots = FALSE)
    environment(out$call) <- parent.frame()
    out
  }

nplsqregbw.lsqregressionbandwidth <- function(bws, tau = bws$tau, ...) {
  tau <- .nplsqreg_validate_tau(tau)
  if (!isTRUE(all.equal(tau, bws$tau)))
    stop("cross-tau nplsqreg bandwidth-object reuse is not supported",
         call. = FALSE)
  bws
}

nplsqregbw.NULL <- function(...) {
  nplsqregbw.default(...)
}

nplsqregbw.default <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           bws,
           tau = 0.5,
           delta = NULL,
           scale = NULL,
           bandwidth.compute = TRUE,
           delta.bounds = c(1e-4, 1 - 1e-4),
           optim.control = list(maxit = 50L),
           ...) {

    elapsed.start <- proc.time()[3]
    tau <- .nplsqreg_validate_tau(tau)
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute,
                                                 "bandwidth.compute")
    xdat <- toFrame(xdat)
    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")
    if (nrow(xdat) != length(ydat))
      stop("number of explanatory data and response data do not match")
    if (is.factor(ydat))
      stop("nplsqreg requires a numeric dependent variable", call. = FALSE)
    ydat <- as.numeric(ydat)
    if (any(!is.finite(ydat)))
      stop("'ydat' must be finite")
    if (!is.numeric(delta.bounds) || length(delta.bounds) != 2L ||
        any(!is.finite(delta.bounds)) || delta.bounds[1] <= 0 ||
        delta.bounds[2] >= 1 || delta.bounds[1] >= delta.bounds[2])
      stop("'delta.bounds' must be a two-value numeric interval inside (0, 1)",
           call. = FALSE)

    dots <- list(...)
    if (!is.null(dots$nomad) && isTRUE(dots$nomad))
      stop("nplsqregbw does not yet support nomad=TRUE degree search",
           call. = FALSE)

    if (is.null(scale)) {
      pilot <- .nplsqreg_scale_pilot(xdat, ydat, dots)
      scale <- pilot$scale
      scale.type <- "fan-yao"
      mean.fit <- pilot$mean.fit
      scale.fit <- pilot$scale.fit
    } else {
      scale <- .nplsqreg_validate_scale(scale, length(ydat))
      scale.type <- "supplied"
      mean.fit <- NULL
      scale.fit <- NULL
    }

    if (missing(bws)) {
      start.bws <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat), dots))
    } else if (isa(bws, "lsqregressionbandwidth")) {
      return(nplsqregbw.lsqregressionbandwidth(bws, tau = tau))
    } else if (isa(bws, "rbandwidth")) {
      start.bws <- bws
    } else {
      start.bws <- do.call(npregbw, c(list(xdat = xdat, ydat = ydat,
                                           bws = bws,
                                           bandwidth.compute = FALSE), dots))
    }

    if (!identical(start.bws$type, "fixed"))
      stop("nplsqregbw currently supports bwtype='fixed' only",
           call. = FALSE)

    if (is.null(delta))
      delta <- tau
    if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
        delta <= delta.bounds[1] || delta >= delta.bounds[2])
      stop("'delta' must be a single finite numeric value inside delta.bounds",
           call. = FALSE)

    start.bw <- as.numeric(start.bws$bw)
    bounds <- .nplsqreg_bw_bounds(start.bws, xdat, start.bw)
    start.par <- c(bounds$start, delta)
    lower <- c(bounds$lower, delta.bounds[1])
    upper <- c(bounds$upper, delta.bounds[2])

    objective <- .nplsqreg_objective_factory(
      xdat = xdat,
      ydat = ydat,
      scale = scale,
      tau = tau,
      template = start.bws,
      delta.bounds = delta.bounds)

    if (bandwidth.compute) {
      opt <- stats::optim(start.par, objective, method = "L-BFGS-B",
                          lower = lower, upper = upper,
                          control = optim.control)
      best.par <- opt$par
      best.value <- opt$value
    } else {
      best.par <- start.par
      best.value <- objective(best.par)
      opt <- list(par = best.par, value = best.value, counts = c("function" = 1L),
                  convergence = 0L, message = "manual bandwidth/delta evaluation")
    }

    best.bw <- best.par[-length(best.par)]
    best.delta <- best.par[length(best.par)]
    qdat <- ydat + scale * stats::qnorm(best.delta)
    state <- environment(objective)$state()
    reg.bws <- .nplsqreg_rebuild_rbandwidth(
      template = start.bws,
      bw = best.bw,
      ydat = qdat,
      fval = best.value,
      num.feval = length(state$f.history),
      bandwidth.compute = bandwidth.compute)
    reg.bws$method <- "cv.check"
    reg.bws$pmethod <- "Check-Loss Cross-Validation"
    reg.bws$total.time <- proc.time()[3] - elapsed.start

    out <- lsqregressionbandwidth(
      reg.bws = reg.bws,
      tau = tau,
      delta = best.delta,
      scale = scale,
      scale.type = scale.type,
      xdat = xdat,
      ydat = ydat,
      qdat = qdat,
      objective = best.value,
      optim = opt,
      mean.fit = mean.fit,
      scale.fit = scale.fit,
      formula = NULL,
      call = match.call(expand.dots = FALSE))
    environment(out$call) <- parent.frame()
    out
  }

nplsqreg.formula <-
  function(bws, data = NULL, newdata = NULL, tau = 0.5,
           gradients = FALSE, residuals = FALSE, ...) {

    tt <- terms(bws)
    mc <- match.call(expand.dots = FALSE)
    m <- match(c("bws", "data", "subset", "na.action"),
               names(mc), nomatch = 0)
    tmf <- mc[c(1, m)]
    if ("bws" %in% names(tmf))
      names(tmf)[names(tmf) == "bws"] <- "formula"
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    mf <- do.call(stats::model.frame, mf.args, envir = environment(tt))
    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"), "term.labels"), drop = FALSE]

    has.eval <- !is.null(newdata)
    if (has.eval) {
      npValidateNewdataFormula(newdata, delete.response(tt),
                               include.response = FALSE)
      emf <- do.call(stats::model.frame,
                     list(formula = delete.response(tt), data = newdata),
                     envir = parent.frame())
      exdat <- emf[, attr(attr(emf, "terms"), "term.labels"), drop = FALSE]
    }

    bw <- nplsqregbw(xdat = xdat, ydat = ydat, tau = tau, ...)
    bw$formula <- bws
    fit.args <- list(bws = bw, txdat = xdat, tydat = ydat,
                     gradients = gradients, residuals = residuals)
    if (has.eval)
      fit.args$exdat <- exdat
    out <- do.call(nplsqreg, fit.args)
    out$call <- match.call(expand.dots = FALSE)
    environment(out$call) <- parent.frame()
    out$bws$formula <- bws
    out
  }

nplsqreg.lsqregressionbandwidth <-
  function(bws, txdat = NULL, tydat = NULL, tau = bws$tau, ...) {
    tau <- .nplsqreg_validate_tau(tau)
    if (!isTRUE(all.equal(tau, bws$tau)))
      stop("cross-tau nplsqreg bandwidth-object reuse is not supported",
           call. = FALSE)
    if (is.null(txdat))
      txdat <- bws$xdat
    if (is.null(tydat))
      tydat <- bws$ydat
    nplsqreg.default(bws = bws, txdat = txdat, tydat = tydat, tau = tau, ...)
  }

nplsqreg.NULL <- function(...) {
  nplsqreg.default(...)
}

nplsqreg.default <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           tau = 0.5,
           exdat,
           gradients = FALSE,
           residuals = FALSE,
           ...) {

    tau <- .nplsqreg_validate_tau(tau)
    gradients <- npValidateScalarLogical(gradients, "gradients")
    residuals <- npValidateScalarLogical(residuals, "residuals")

    if (missing(bws) || !isa(bws, "lsqregressionbandwidth")) {
      bw.args <- list(xdat = txdat, ydat = tydat, tau = tau)
      if (!missing(bws))
        bw.args$bws <- bws
      bw <- do.call(nplsqregbw, c(bw.args, list(...)))
      fit.args <- list(bws = bw, txdat = txdat, tydat = tydat,
                       gradients = gradients, residuals = residuals)
      if (!missing(exdat))
        fit.args$exdat <- exdat
      return(do.call(nplsqreg, fit.args))
    }

    if (!isTRUE(all.equal(tau, bws$tau)))
      stop("cross-tau nplsqreg bandwidth-object reuse is not supported",
           call. = FALSE)

    txdat <- toFrame(txdat)
    tydat <- as.numeric(tydat)
    if (nrow(txdat) != length(tydat))
      stop("number of explanatory data and response data do not match")

    fit.args <- list(bws = bws$reg.bws, txdat = txdat, tydat = bws$qdat,
                     gradients = gradients, residuals = residuals)
    if (!missing(exdat))
      fit.args$exdat <- exdat
    fit <- do.call(npreg, c(fit.args, list(...)))

    quant <- fitted(fit)
    qerr <- se(fit)
    qgrad <- if (gradients) gradients(fit) else NA
    qgerr <- if (gradients) gradients(fit, errors = TRUE) else NA
    resid.out <- if (residuals) tydat - fitted(npreg(bws = bws$reg.bws,
                                                    txdat = txdat,
                                                    tydat = bws$qdat)) else NA

    out <- lsqregression(
      bws = bws,
      fit = fit,
      xeval = fit$eval,
      tau = tau,
      delta = bws$delta,
      quantile = quant,
      quanterr = qerr,
      quantgrad = qgrad,
      quantgerr = qgerr,
      ntrain = nrow(txdat),
      trainiseval = missing(exdat),
      gradients = gradients,
      residuals = residuals,
      resid = resid.out,
      call = match.call(expand.dots = FALSE))
    environment(out$call) <- parent.frame()
    out
  }
