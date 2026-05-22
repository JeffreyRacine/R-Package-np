lsqregressionbandwidth <-
  function(reg.bws,
           tau,
           delta,
           scale,
           scale.type,
           xdat,
           ydat,
           qdat,
           objective,
           optim = NULL,
           mean.fit = NULL,
           scale.fit = NULL,
           formula = NULL,
           call = NULL) {

    if (missing(reg.bws) || missing(tau) || missing(delta) ||
        missing(scale) || missing(xdat) || missing(ydat) || missing(qdat))
      stop("improper invocation of lsqregressionbandwidth constructor")

    d <- list(
      bw = reg.bws$bw,
      bws = reg.bws,
      reg.bws = reg.bws,
      tau = tau,
      delta = delta,
      scale = scale,
      scale.type = scale.type,
      xdat = xdat,
      ydat = ydat,
      qdat = qdat,
      objective = objective,
      optim = optim,
      mean.fit = mean.fit,
      scale.fit = scale.fit,
      formula = formula,
      call = call,
      xnames = reg.bws$xnames,
      ynames = reg.bws$ynames,
      nobs = reg.bws$nobs,
      ndim = reg.bws$ndim,
      nord = reg.bws$nord,
      nuno = reg.bws$nuno,
      ncon = reg.bws$ncon,
      pscaling = reg.bws$pscaling,
      ptype = reg.bws$ptype,
      pckertype = reg.bws$pckertype,
      pukertype = reg.bws$pukertype,
      pokertype = reg.bws$pokertype,
      regtype = reg.bws$regtype,
      pregtype = reg.bws$pregtype,
      basis = reg.bws$basis,
      degree = reg.bws$degree,
      bernstein.basis = reg.bws$bernstein.basis,
      method = "cv.check",
      pmethod = "Check-Loss Cross-Validation")

    class(d) <- "lsqregressionbandwidth"
    d
  }

lsqregression <-
  function(bws,
           fit,
           xeval,
           tau,
           delta,
           quantile,
           quanterr = NA,
           quantgrad = NA,
           quantgerr = NA,
           ntrain,
           trainiseval = FALSE,
           gradients = FALSE,
           residuals = FALSE,
           resid = NA,
           call = NULL) {

    if (missing(bws) || missing(fit) || missing(xeval) ||
        missing(tau) || missing(delta) || missing(quantile) ||
        missing(ntrain))
      stop("improper invocation of lsqregression constructor")

    d <- list(
      bw = bws$bw,
      bws = bws,
      reg.bws = bws$reg.bws,
      fit = fit,
      xnames = bws$xnames,
      ynames = bws$ynames,
      nobs = nrow(xeval),
      ndim = bws$ndim,
      nord = bws$nord,
      nuno = bws$nuno,
      ncon = bws$ncon,
      pscaling = bws$pscaling,
      ptype = bws$ptype,
      pckertype = bws$pckertype,
      pukertype = bws$pukertype,
      pokertype = bws$pokertype,
      xeval = xeval,
      tau = tau,
      delta = delta,
      scale = bws$scale,
      scale.type = bws$scale.type,
      quantile = quantile,
      quanterr = quanterr,
      quantgrad = quantgrad,
      quantgerr = quantgerr,
      ntrain = ntrain,
      trainiseval = trainiseval,
      gradients = gradients,
      residuals = residuals,
      resid = resid,
      objective = bws$objective,
      optim = bws$optim,
      call = call)

    class(d) <- "lsqregression"
    d
  }

.nplsqreg_format_compact <- function(x, digits = NULL) {
  if (is.null(x) || !length(x))
    return("")
  x <- as.numeric(x)
  if (is.null(digits))
    paste(format(x, trim = TRUE, scientific = FALSE), collapse = ",")
  else
    paste(format(x, digits = digits, trim = TRUE, scientific = FALSE),
          collapse = ",")
}

.nplsqreg_tau_state_table <- function(x, digits = NULL) {
  bw.list <- if (!is.null(x$tau.bws)) {
    x$tau.bws
  } else if (!is.null(x$bws) && !is.null(x$bws$tau.bws)) {
    x$bws$tau.bws
  } else {
    NULL
  }
  if (is.null(bw.list))
    return(data.frame(tau = x$tau,
                      delta = as.numeric(x$delta),
                      objective = as.numeric(x$objective)))

  rows <- lapply(seq_along(bw.list), function(j) {
    one <- bw.list[[j]]
    rbw <- one$reg.bws
    dsearch <- rbw$degree.search
    data.frame(
      tau = as.numeric(x$tau[[j]]),
      delta = as.numeric(one$delta),
      objective = as.numeric(one$objective),
      regtype = if (is.null(rbw$regtype)) "" else rbw$regtype,
      degree = .nplsqreg_format_compact(rbw$degree, digits = digits),
      bwtype = if (is.null(rbw$type)) "" else rbw$type,
      bandwidth = .nplsqreg_format_compact(rbw$bw, digits = digits),
      search = if (!is.null(dsearch$mode)) dsearch$mode else "powell",
      nomad = isTRUE(rbw$nomad.shortcut$enabled),
      feval = if (is.null(rbw$num.feval)) NA_real_ else as.numeric(rbw$num.feval),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

print.lsqregressionbandwidth <- function(x, digits = NULL, ...) {
  cat("\nLocation-scale quantile regression bandwidth object\n", sep = "")
  if (length(x$tau) > 1L) {
    cat("Tau values: ", paste(format(x$tau, trim = TRUE), collapse = ", "),
        "\n", sep = "")
    cat("Tau search: ", if (is.null(x$tau.search)) "full" else x$tau.search,
        "\n\n", sep = "")
    print(.nplsqreg_tau_state_table(x, digits = digits), row.names = FALSE)
  } else {
    cat("Tau: ", format(x$tau, trim = TRUE),
        "  Delta: ", format(x$delta, trim = TRUE),
        "  Objective: ", format(x$objective, trim = TRUE), "\n\n", sep = "")
    print(x$reg.bws, digits = digits)
  }
  invisible(x)
}

summary.lsqregressionbandwidth <- function(object, ...) {
  print.lsqregressionbandwidth(object, ...)
}

print.lsqregression <- function(x, digits = NULL, ...) {
  cat("\nLocation-scale quantile regression data: ", x$ntrain,
      " training points,",
      if (x$trainiseval) "" else paste(" and ", x$nobs,
                                       " evaluation points,", sep = ""),
      " in ", x$ndim, " variable(s)\n", sep = "")
  if (length(x$tau) > 1L) {
    cat("Tau values: ", paste(format(x$tau, trim = TRUE), collapse = ", "),
        "\n", sep = "")
    cat("Tau search: ", if (is.null(x$tau.search)) "full" else x$tau.search,
        "\n\n", sep = "")
    print(.nplsqreg_tau_state_table(x, digits = digits), row.names = FALSE)
  } else {
    cat("Tau: ", format(x$tau, trim = TRUE),
        "  Delta: ", format(x$delta, trim = TRUE),
        "  Objective: ", format(x$objective, trim = TRUE), "\n\n", sep = "")
    print(matrix(x$bw, ncol = x$ndim,
                 dimnames = list(paste(x$pscaling, ":", sep = ""),
                                 x$xnames)))
    cat(genRegEstStr(x$fit))
    cat(genBwKerStrs(x$reg.bws))
  }
  cat("\n\n")
  if (!missing(...))
    print(..., digits = digits)
  invisible(x)
}

summary.lsqregression <- function(object, ...) {
  print.lsqregression(object, ...)
}

fitted.lsqregression <- function(object, ...) {
  object$quantile
}

quantile.lsqregression <- function(x, ...) {
  x$quantile
}

se.lsqregression <- function(x) {
  x$quanterr
}

gradients.lsqregression <- function(x, errors = FALSE, ...) {
  errors <- npValidateScalarLogical(errors, "errors")
  gout <- if (!errors) x$quantgrad else x$quantgerr
  if (is.null(gout) || (length(gout) == 1L && is.logical(gout) && is.na(gout)))
    stop(if (!errors)
      "gradients are not available: fit the model with gradients=TRUE"
    else
      "gradient standard errors are not available: fit the model with gradients=TRUE")
  gout
}

.nplsqreg_predict_formula_newdata_to_exdat <- function(object, newdata) {
  tt <- stats::terms(object$bws$formula)
  rhs <- stats::delete.response(tt)
  npValidateNewdataFormula(newdata, rhs, include.response = FALSE)
  mf <- stats::model.frame(formula = rhs, data = newdata)
  mf[, attr(attr(mf, "terms"), "term.labels"), drop = FALSE]
}

predict.lsqregression <- function(object, se.fit = FALSE, ...) {
  se.fit <- npValidateScalarLogical(se.fit, "se.fit")
  if (length(object$tau) > 1L) {
    if (is.null(object$tau.fits) || length(object$tau.fits) != length(object$tau))
      stop("vector nplsqreg object lacks per-tau fit state", call. = FALSE)
    labels <- .nplsqreg_tau_labels(object$tau)
    pred <- lapply(object$tau.fits, predict.lsqregression,
                   se.fit = se.fit, ...)
    if (se.fit) {
      fit <- do.call(cbind, lapply(pred, `[[`, "fit"))
      se.out <- do.call(cbind, lapply(pred, `[[`, "se.fit"))
      colnames(fit) <- labels
      colnames(se.out) <- labels
      return(list(fit = fit,
                  se.fit = se.out,
                  df = pred[[1L]]$df,
                  residual.scale = pred[[1L]]$residual.scale))
    }
    out <- do.call(cbind, pred)
    colnames(out) <- labels
    return(out)
  }
  dots <- list(...)
  has.formula.route <- !is.null(object$bws$formula)

  if (!is.null(dots$exdat) && !is.null(dots$newdata))
    dots$newdata <- NULL

  if (has.formula.route) {
    if (is.null(dots$exdat) && !is.null(dots$newdata)) {
      dots$exdat <- .nplsqreg_predict_formula_newdata_to_exdat(object, dots$newdata)
      dots$newdata <- NULL
    }
  } else {
    if (is.null(dots$exdat) && !is.null(dots$newdata)) {
      dots$exdat <- dots$newdata
      dots$newdata <- NULL
    }
  }

  fit.args <- list(bws = object$reg.bws,
                   txdat = object$bws$xdat,
                   tydat = object$bws$qdat)
  if (!is.null(dots$exdat))
    fit.args$exdat <- dots$exdat
  fit.args <- c(fit.args, dots[setdiff(names(dots), c("newdata", "exdat"))])

  tr <- do.call(npreg, fit.args)
  if (se.fit) {
    list(fit = fitted(tr), se.fit = se(tr),
         df = tr$nobs, residual.scale = tr$MSE)
  } else {
    fitted(tr)
  }
}

plot.lsqregression <- function(x, tau = NULL, gradient = FALSE,
                               gradients = gradient,
                               plot.data.overlay = NULL,
                               xlab = NULL, ylab = NULL, ylim = NULL, ...) {
  dots.call <- match.call(expand.dots = FALSE)$...
  .np_plot_validate_public_dots(
    dots.call,
    method = .np_plot_rbandwidth_engine,
    bws = x$reg.bws,
    context = "plot.lsqregression"
  )
  dots <- .np_plot_normalize_public_dots(list(...),
                                         context = "plot.lsqregression")
  if (!is.null(dots$gradient)) {
    gradients <- dots$gradient
    dots$gradient <- NULL
  }
  if (!is.null(dots$gradients)) {
    gradients <- dots$gradients
    dots$gradients <- NULL
  }
  gradients <- npValidateScalarLogical(gradients, "gradients")
  if (is.null(plot.data.overlay)) {
    plot.data.overlay <- if (!is.null(dots$plot.data.overlay)) {
      dots$plot.data.overlay
    } else {
      !isTRUE(gradients)
    }
  }
  dots$plot.data.overlay <- NULL
  plot.data.overlay <- npValidateScalarLogical(plot.data.overlay,
                                               "plot.data.overlay")
  plot.behavior <- if (is.null(dots$plot.behavior)) {
    "plot"
  } else {
    match.arg(dots$plot.behavior, c("plot", "plot-data", "data"))
  }
  dots$plot.behavior <- NULL
  plot.errors.method <- if (is.null(dots$plot.errors.method)) {
    "none"
  } else {
    match.arg(dots$plot.errors.method,
              c("none", "bootstrap", "asymptotic"))
  }
  dots$plot.errors.method <- plot.errors.method
  plot.errors <- !identical(plot.errors.method, "none")
  dots$plot.rug <- NULL
  plot.legend <- if (is.null(dots$legend)) TRUE else dots$legend
  dots$legend <- NULL

  refit_for_plot <- function(object, tau.values, want.gradients) {
    nplsqreg(
      bws = object$bws,
      txdat = object$bws$xdat,
      tydat = object$bws$ydat,
      exdat = object$xeval,
      tau = tau.values,
      gradients = want.gradients
    )
  }
  if (!is.null(tau)) {
    tau <- .nplsqreg_validate_tau_values(tau)
    if (length(tau) != length(x$tau) || !isTRUE(all.equal(tau, x$tau))) {
      if (length(x$tau) != 1L)
        stop("plot tau expansion from an existing vector-tau nplsqreg object is not supported; refit with the requested tau values",
             call. = FALSE)
      x <- refit_for_plot(x, tau, gradients)
    }
  }
  if (isTRUE(gradients) && !isTRUE(x$gradients))
    x <- refit_for_plot(x, x$tau, TRUE)
  if (x$ndim != 1L)
    stop("plot.lsqregression currently supports one explanatory variable",
         call. = FALSE)

  plot_data_one <- function(object) {
    args <- c(
      list(
        bws = object$reg.bws,
        xdat = object$bws$xdat,
        ydat = object$bws$qdat,
        gradients = gradients,
        plot.behavior = "data",
        plot.par.mfrow = FALSE
      ),
      dots
    )
    do.call(.np_plot_rbandwidth_engine, args)[[1L]]
  }

  multi.tau <- length(x$tau) > 1L
  tau.fits <- if (multi.tau) x$tau.fits else list(x)
  tau.labels <- .nplsqreg_tau_labels(x$tau)
  plot.out <- stats::setNames(lapply(tau.fits, plot_data_one), tau.labels)
  get.value <- function(object)
    as.numeric(if (gradients) gradients(object) else fitted(object))
  get.se <- function(object)
    if (gradients) gradients(object, errors = TRUE) else se(object)
  y <- do.call(cbind, lapply(plot.out, get.value))
  colnames(y) <- tau.labels
  yerr <- NULL
  if (plot.errors) {
    yerr <- array(NA_real_, dim = c(nrow(y), 3L, ncol(y)),
                  dimnames = list(NULL, NULL, tau.labels))
    for (j in seq_along(plot.out)) {
      se.j <- get.se(plot.out[[j]])
      if (is.null(se.j))
        next
      se.j <- as.matrix(se.j)
      yerr[, 1L, j] <- -se.j[, 1L]
      yerr[, 2L, j] <- se.j[, 2L]
      bias.j <- if (gradients) plot.out[[j]]$gbias else plot.out[[j]]$bias
      if (!is.null(bias.j) && length(bias.j) == nrow(y))
        yerr[, 3L, j] <- y[, j] - bias.j
    }
  }

  if (identical(plot.behavior, "data"))
    return(plot.out)

  if (is.null(xlab))
    xlab <- x$xnames[1L]
  if (is.null(ylab))
    ylab <- if (gradients) {
      "Gradient"
    } else if (multi.tau) {
      "Conditional quantile"
    } else {
      paste(x$tau, "quantile")
    }

  xeval <- plot.out[[1L]]$eval[[1L]]
  xfactor <- is.factor(xeval)

  plot.user.args <- .np_plot_user_args(dots, "plot")
  .np_plot_quantile_overlay_1d(
    ei = xeval,
    value = y,
    xi.factor = xfactor,
    xlab.value = xlab,
    ylab.value = ylab,
    err = yerr,
    all.err = NULL,
    tau.labels = tau.labels,
    overlay.x = x$bws$xdat[[1L]],
    overlay.y = x$bws$ydat,
    plot.data.overlay = plot.data.overlay,
    gradients = gradients,
    plotOnEstimate = identical(plot.errors.method, "none") ||
      identical(dots$plot.errors.center, "estimate"),
    plot.errors.type = if (is.null(dots$plot.errors.type)) "pmzsd" else dots$plot.errors.type,
    plot.errors.style = if (is.null(dots$plot.errors.style)) "band" else dots$plot.errors.style,
    plot.errors.bar = if (is.null(dots$plot.errors.bar)) "|" else dots$plot.errors.bar,
    plot.errors.bar.num = if (is.null(dots$plot.errors.bar.num)) min(nrow(y), 25) else dots$plot.errors.bar.num,
    plot.user.args = plot.user.args,
    legend = plot.legend,
    col = dots$col,
    lty = dots$lty,
    lwd = dots$lwd,
    type = dots$type,
    xlim = dots$xlim,
    ylim = ylim,
    main = dots$main,
    sub = dots$sub,
    cex.axis = dots$cex.axis,
    cex.lab = dots$cex.lab,
    cex.main = dots$cex.main,
    cex.sub = dots$cex.sub
  )
  if (identical(plot.behavior, "plot-data"))
    return(plot.out)
  invisible(x)
}
