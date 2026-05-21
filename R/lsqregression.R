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

print.lsqregressionbandwidth <- function(x, digits = NULL, ...) {
  cat("\nLocation-scale quantile regression bandwidth object\n", sep = "")
  if (length(x$tau) > 1L) {
    cat("Tau values: ", paste(format(x$tau, trim = TRUE), collapse = ", "),
        "\n", sep = "")
    cat("Tau search: ", if (is.null(x$tau.search)) "full" else x$tau.search,
        "\n\n", sep = "")
    print(data.frame(tau = x$tau,
                     delta = as.numeric(x$delta),
                     objective = as.numeric(x$objective)))
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
    print(data.frame(tau = x$tau,
                     delta = as.numeric(x$delta),
                     objective = as.numeric(x$objective)))
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

.nplsqreg_predict_newdata_to_exdat <- function(object, newdata) {
  if (is.data.frame(newdata) && !is.null(names(newdata))) {
    xnames <- object$xnames
    if (length(xnames)) {
      missing.names <- setdiff(xnames, names(newdata))
      if (length(missing.names))
        stop(sprintf(
          "newdata must contain columns: %s",
          paste(shQuote(xnames), collapse = ", ")
        ), call. = FALSE)
      return(newdata[, xnames, drop = FALSE])
    }
  }
  newdata
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
    if (!is.null(dots$exdat)) {
      dots$exdat <- .nplsqreg_predict_newdata_to_exdat(object, dots$exdat)
    } else if (!is.null(dots$newdata)) {
      dots$exdat <- .nplsqreg_predict_newdata_to_exdat(object, dots$newdata)
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

plot.lsqregression <- function(x, gradient = FALSE,
                               plot.data.overlay = !gradient,
                               xlab = NULL, ylab = NULL, ylim = NULL, ...) {
  gradient <- npValidateScalarLogical(gradient, "gradient")
  plot.data.overlay <- npValidateScalarLogical(plot.data.overlay,
                                               "plot.data.overlay")
  if (x$ndim != 1L)
    stop("plot.lsqregression currently supports one explanatory variable",
         call. = FALSE)

  xraw <- x$xeval[[1L]]
  xpos <- if (is.factor(xraw)) as.numeric(xraw) else as.numeric(xraw)
  y <- if (gradient) gradients.lsqregression(x) else fitted.lsqregression(x)
  multi.tau <- length(x$tau) > 1L
  if (gradient && multi.tau) {
    if (length(dim(y)) != 3L || dim(y)[2L] != 1L)
      stop("vector-tau gradient plots require one explanatory variable",
           call. = FALSE)
    y <- y[, 1L, , drop = FALSE]
    dim(y) <- c(dim(y)[1L], dim(y)[3L])
    colnames(y) <- .nplsqreg_tau_labels(x$tau)
  }
  if (!gradient && multi.tau && is.null(colnames(y)))
    colnames(y) <- .nplsqreg_tau_labels(x$tau)
  ord <- order(xpos)

  if (is.null(xlab))
    xlab <- x$xnames[1L]
  if (is.null(ylab))
    ylab <- if (gradient) {
      "Gradient"
    } else if (multi.tau) {
      "Conditional quantile"
    } else {
      paste("Quantile tau =", format(x$tau, trim = TRUE))
    }

  if (is.null(ylim)) {
    yr <- range(y, finite = TRUE)
    if (plot.data.overlay && !gradient)
      yr <- range(yr, x$bws$ydat, finite = TRUE)
    ylim <- grDevices::extendrange(yr)
  }

  graphics::plot(xpos[ord], if (multi.tau) y[ord, 1L] else y[ord],
                 type = "n", xlab = xlab, ylab = ylab,
                 ylim = ylim, xaxt = if (is.factor(xraw)) "n" else "s", ...)
  if (is.factor(xraw)) {
    graphics::axis(1, at = seq_along(levels(xraw)), labels = levels(xraw))
  }
  if (plot.data.overlay && !gradient) {
    xtrain <- x$bws$xdat[[1L]]
    if (is.factor(xtrain)) {
      graphics::boxplot(x$bws$ydat ~ xtrain, add = TRUE, axes = FALSE,
                        outline = FALSE, col = "grey90", border = "grey60")
    } else {
      xtrain <- as.numeric(xtrain)
      graphics::points(xtrain, x$bws$ydat)
    }
  }
  if (multi.tau) {
    cols <- seq_len(length(x$tau))
    mat <- as.matrix(y)
    for (j in seq_along(x$tau))
      graphics::lines(xpos[ord], mat[ord, j], col = cols[j], lty = j)
    graphics::legend("topright", legend = .nplsqreg_tau_labels(x$tau),
                     col = cols, lty = seq_along(x$tau), bty = "n")
  } else {
    graphics::lines(xpos[ord], y[ord])
  }
  invisible(x)
}
