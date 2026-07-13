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
      xbw = reg.bws$bw,
      ybw = rep.int(NA_real_, 1L),
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
      xndim = reg.bws$ndim,
      yndim = 1L,
      nord = reg.bws$nord,
      nuno = reg.bws$nuno,
      ncon = reg.bws$ncon,
      xnord = reg.bws$nord,
      xnuno = reg.bws$nuno,
      xncon = reg.bws$ncon,
      ynord = 0L,
      ynuno = 0L,
      yncon = 1L,
      ixcon = reg.bws$icon,
      ixuno = reg.bws$iuno,
      ixord = reg.bws$iord,
      iycon = TRUE,
      iyuno = FALSE,
      iyord = FALSE,
      pscaling = reg.bws$pscaling,
      ptype = reg.bws$ptype,
      pckertype = reg.bws$pckertype,
      pukertype = reg.bws$pukertype,
      pokertype = reg.bws$pokertype,
      pcxkertype = reg.bws$pckertype,
      puxkertype = reg.bws$pukertype,
      poxkertype = reg.bws$pokertype,
      pcykertype = "unused",
      puykertype = "unused",
      poykertype = "unused",
      xdati = reg.bws$xdati,
      ydati = reg.bws$ydati,
      regtype = reg.bws$regtype,
      pregtype = reg.bws$pregtype,
      basis = reg.bws$basis,
      degree = reg.bws$degree,
      bernstein.basis = reg.bws$bernstein.basis,
      method = "cv.check",
      pmethod = "Check-Loss Cross-Validation",
      fval = objective,
      num.feval = reg.bws$num.feval,
      num.feval.fast = reg.bws$num.feval.fast,
      rows.omit = reg.bws$rows.omit,
      nobs.omit = reg.bws$nobs.omit,
      total.time = reg.bws$total.time,
      optim.time = reg.bws$total.time,
      fit.time = NA)

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
           call = NULL,
           gradient.order = fit$gradient.order) {

    if (missing(bws) || missing(fit) || missing(xeval) ||
        missing(tau) || missing(delta) || missing(quantile) ||
        missing(ntrain))
      stop("improper invocation of lsqregression constructor")

    d <- list(
      bw = bws$bw,
      xbw = bws$xbw,
      ybw = bws$ybw,
      bws = bws,
      reg.bws = bws$reg.bws,
      fit = fit,
      xnames = bws$xnames,
      ynames = bws$ynames,
      nobs = nrow(xeval),
      ndim = bws$ndim,
      xndim = bws$xndim,
      yndim = bws$yndim,
      nord = bws$nord,
      nuno = bws$nuno,
      ncon = bws$ncon,
      xnord = bws$xnord,
      xnuno = bws$xnuno,
      xncon = bws$xncon,
      ynord = bws$ynord,
      ynuno = bws$ynuno,
      yncon = bws$yncon,
      pscaling = bws$pscaling,
      ptype = bws$ptype,
      pckertype = bws$pckertype,
      pukertype = bws$pukertype,
      pokertype = bws$pokertype,
      pcxkertype = bws$pcxkertype,
      puxkertype = bws$puxkertype,
      poxkertype = bws$poxkertype,
      pcykertype = bws$pcykertype,
      puykertype = bws$puykertype,
      poykertype = bws$poykertype,
      regtype = bws$regtype,
      pregtype = bws$pregtype,
      basis = bws$basis,
      degree = bws$degree,
      bernstein.basis = bws$bernstein.basis,
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
      gradient.order = gradient.order,
      residuals = residuals,
      resid = resid,
      objective = bws$objective,
      optim = bws$optim,
      total.time = .nplsqreg_sum_times(bws$total.time, fit$total.time),
      optim.time = bws$total.time,
      fit.time = fit$total.time,
      call = call)

    class(d) <- "lsqregression"

    reset <- get0(".npRmpi_nplsqreg_reset_worker_comm_state",
                  envir = asNamespace("npRmpi"),
                  inherits = FALSE)
    if (is.function(reset) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
      reset(comm = 1L)

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

.nplsqreg_sum_times <- function(...) {
  x <- unlist(list(...), recursive = TRUE, use.names = FALSE)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x))
    return(NA_real_)
  sum(x)
}

.nplsqreg_primary_bws <- function(x) {
  if (inherits(x, "lsqregression"))
    x <- x$bws
  if (length(x$tau) > 1L && !is.null(x$tau.bws))
    return(x$tau.bws[[1L]])
  x
}

.nplsqreg_primary_reg_bws <- function(x) {
  bws <- .nplsqreg_primary_bws(x)
  if (!is.null(bws$reg.bws))
    return(bws$reg.bws)
  bws
}

.nplsqreg_search_label <- function(x) {
  methods <- if (length(x$tau) > 1L && !is.null(x$tau.bws)) {
    vapply(x$tau.bws, function(z) {
      ds <- z$reg.bws$degree.search
      if (!is.null(ds$mode))
        return(as.character(ds$mode)[1L])
      if (!is.null(ds$engine))
        return(as.character(ds$engine)[1L])
      if (isTRUE(z$reg.bws$nomad.shortcut$enabled))
        return("nomad")
      "powell"
    }, character(1L))
  } else {
    rbw <- .nplsqreg_primary_reg_bws(x)
    ds <- rbw$degree.search
    if (!is.null(ds$mode))
      as.character(ds$mode)[1L]
    else if (!is.null(ds$engine))
      as.character(ds$engine)[1L]
    else if (isTRUE(rbw$nomad.shortcut$enabled))
      "nomad"
    else
      "powell"
  }
  if (length(unique(methods)) == 1L) methods[[1L]] else "tau-specific"
}

.nplsqreg_bandwidth_display <- function(x) {
  bws <- if (inherits(x, "lsqregression")) x$bws else x
  label <- paste("Exp. Var. ", bws$pscaling, ":", sep = "")
  if (length(bws$tau) > 1L) {
    bw <- as.matrix(bws$bw)
    rownames(bw) <- bws$xnames
    colnames(bw) <- .nplsqreg_tau_labels(bws$tau)
    cat(label, "\n", sep = "")
    print(t(bw))
  } else {
    print(matrix(bws$bw, ncol = bws$xndim,
                 dimnames = list(label, bws$xnames)))
  }
}

.nplsqreg_delta_objective_display <- function(x) {
  bws <- if (inherits(x, "lsqregression")) x$bws else x
  cat("\nDelta:", paste(format(bws$delta, trim = TRUE), collapse = ", "))
  cat("\nCheck-Loss Objective:",
      paste(format(bws$objective, trim = TRUE), collapse = ", "))
  cat("\nSearch Method:", .nplsqreg_search_label(bws))
  if (length(bws$tau) > 1L)
    cat("\nTau Search:", if (is.null(bws$tau.search)) "full" else bws$tau.search)
  cat("\n")
}

.nplsqreg_degree_matrix <- function(x) {
  bws <- if (inherits(x, "lsqregression")) x$bws else x
  if (length(bws$tau) <= 1L || is.null(bws$tau.bws))
    return(NULL)

  first.rbw <- bws$tau.bws[[1L]]$reg.bws
  if (is.null(first.rbw) || !identical(first.rbw$regtype, "lp") ||
      isTRUE(first.rbw$ncon <= 0L))
    return(NULL)

  has.degree.search <- any(vapply(
    bws$tau.bws,
    function(z) {
      rbw <- z$reg.bws
      isTRUE(rbw$nomad.shortcut$enabled) || !is.null(rbw$degree.search)
    },
    logical(1L)))
  if (!isTRUE(has.degree.search) && isTRUE(bws$child.degree.common))
    return(NULL)

  cont.names <- first.rbw$xnames[first.rbw$icon]
  if (!length(cont.names))
    cont.names <- paste0("x", seq_len(first.rbw$ncon))

  degree.list <- lapply(bws$tau.bws, function(z) {
    d <- z$reg.bws$degree
    if (is.null(d) || !length(d))
      return(rep.int(NA_integer_, length(cont.names)))
    d <- as.integer(d)
    if (length(d) == length(cont.names))
      return(d)
    rep_len(d, length(cont.names))
  })

  out <- do.call(rbind, degree.list)
  rownames(out) <- .nplsqreg_tau_labels(bws$tau)
  colnames(out) <- cont.names
  out
}

.nplsqreg_degree_display <- function(x) {
  degree.mat <- .nplsqreg_degree_matrix(x)
  if (is.null(degree.mat))
    return(invisible(FALSE))

  cat("\nContinuous LP Degree(s):\n")
  print(degree.mat)
  invisible(TRUE)
}

.nplsqreg_summary_kernel_source <- function(x) {
  .nplsqreg_primary_reg_bws(x)
}

.nplsqreg_regression_type_label <- function(x) {
  bws <- if (inherits(x, "lsqregression")) x$bws else x
  rbw <- .nplsqreg_summary_kernel_source(x)
  if (length(bws$tau) > 1L && identical(rbw$regtype, "lp")) {
    rbw <- rbw
    rbw$child.degree.common <- isTRUE(bws$child.degree.common)
  }
  npFormatRegressionType(rbw)
}

.nplsqreg_summary_header <- function(x, bandwidth = FALSE) {
  bws <- if (inherits(x, "lsqregression")) x$bws else x
  if (bandwidth) {
    cat("\nLocation-Scale Quantile Regression Bandwidth Data: ",
        bws$nobs, " observations, ",
        bws$xndim + bws$yndim, " variable(s)",
        "\n(", bws$yndim, " dependent variable(s), and ", bws$xndim,
        " explanatory variable(s))\n\n", sep = "")
  } else {
    cat("\nLocation-Scale Quantile Regression Data: ", x$ntrain,
        " training points,",
        if (x$trainiseval) "" else paste(" and ", x$nobs,
                                         " evaluation points,\n", sep = ""),
        " in ", bws$xndim + bws$yndim, " variable(s)",
        "\n(", bws$yndim, " dependent variable(s), and ", bws$xndim,
        " explanatory variable(s))\n\n", sep = "")
  }
}

.nplsqreg_summary_common <- function(x, bandwidth = FALSE) {
  bws <- if (inherits(x, "lsqregression")) x$bws else x
  rbw <- .nplsqreg_summary_kernel_source(x)
  .nplsqreg_summary_header(x, bandwidth = bandwidth)
  cat(genOmitStr(bws))
  if (length(bws$tau) > 1L)
    cat("Tau values:", paste(format(bws$tau, trim = TRUE), collapse = ", "),
        "\n\n")

  cat("Dep. Var. Name:", paste(bws$ynames, collapse = ", "), "\n\n")
  .nplsqreg_bandwidth_display(bws)
  .nplsqreg_delta_objective_display(bws)
  .nplsqreg_degree_display(bws)

  if (bandwidth) {
    cat("\nRegression Type:", .nplsqreg_regression_type_label(bws))
    cat("\nBandwidth Selection Method:", bws$pmethod)
    if (!identical(bws$formula, NULL))
      cat("\nFormula:", paste(deparse(bws$formula), collapse = "\n"))
    cat("\nBandwidth Type:", bws$ptype)
    if (npRegressionHasExtendedNn(bws))
      cat("\nExtended NN: K above n-1 scales the saturated nearest-neighbor bandwidth")
  } else {
    cat(genRegEstStr(x))
  }
  cat(genBwKerStrs(rbw))
  cat(genTimingStr(x))
  cat("\n\n")
  invisible(x)
}

.nplsqreg_formula_response_name <- function(formula) {
  paste(deparse(stats::formula(formula)[[2L]]), collapse = "")
}

.nplsqreg_set_response_name <- function(x, response.name) {
  x$ynames <- response.name
  if (!is.null(x$reg.bws))
    x$reg.bws$ynames <- response.name
  if (!is.null(x$bws))
    x$bws <- .nplsqreg_set_response_name(x$bws, response.name)
  if (!is.null(x$tau.bws))
    x$tau.bws <- lapply(x$tau.bws, .nplsqreg_set_response_name,
                        response.name = response.name)
  if (!is.null(x$tau.fits))
    x$tau.fits <- lapply(x$tau.fits, .nplsqreg_set_response_name,
                         response.name = response.name)
  x
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
  .nplsqreg_summary_common(object, bandwidth = TRUE)
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
  .nplsqreg_summary_common(object, bandwidth = FALSE)
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

.nplsqreg_residuals_unavailable <- function(x) {
  is.null(x) ||
    !length(x) ||
    (length(x) == 1L && (is.logical(x) || is.numeric(x)) && is.na(x))
}

residuals.lsqregression <- function(object, ...) {
  if (length(object$tau) > 1L) {
    if (is.null(object$tau.fits) || length(object$tau.fits) != length(object$tau))
      stop("residuals are not available: vector nplsqreg object lacks per-tau fit state",
           call. = FALSE)
    resid.list <- lapply(object$tau.fits, function(z) z$resid)
    if (any(vapply(resid.list, .nplsqreg_residuals_unavailable, logical(1L))))
      stop("residuals are not available: refit with residuals=TRUE",
           call. = FALSE)
    out <- do.call(cbind, resid.list)
    colnames(out) <- .nplsqreg_tau_labels(object$tau)
    return(out)
  }

  if (.nplsqreg_residuals_unavailable(object$resid))
    stop("residuals are not available: refit with residuals=TRUE",
         call. = FALSE)
  object$resid
}

gradients.lsqregression <- function(x, errors = FALSE,
                                    gradient.order = NULL, ...) {
  errors <- npValidateScalarLogical(errors, "errors")
  gout <- if (!errors) x$quantgrad else x$quantgerr
  if (is.null(gout) || (length(gout) == 1L && is.logical(gout) && is.na(gout)))
    stop(if (!errors)
      "gradients are not available: fit the model with gradients=TRUE"
    else
      "gradient standard errors are not available: fit the model with gradients=TRUE")
  if (!is.null(gradient.order)) {
    fit.list <- if (inherits(x$fit, "npregression")) list(x$fit) else x$fit
    if (!is.list(fit.list) || !length(fit.list) ||
        !all(vapply(fit.list, inherits, logical(1L), "npregression"))) {
      stop("cannot validate requested gradient.order: fitted npregression state is unavailable",
           call. = FALSE)
    }
    invisible(lapply(fit.list, function(fit)
      gradients(fit, errors = errors, gradient.order = gradient.order)))
  }
  gout
}

.nplsqreg_match_plot_tau <- function(requested, stored) {
  requested <- .nplsqreg_validate_tau_values(requested,
                                             allow.duplicates = TRUE)
  stored <- .nplsqreg_validate_tau_values(stored,
                                          allow.duplicates = TRUE)
  idx <- integer(length(requested))
  tol <- sqrt(.Machine$double.eps)
  for (j in seq_along(requested)) {
    hit <- which(abs(stored - requested[[j]]) <= tol)
    if (!length(hit))
      stop("plot tau values for nplsqreg must have been fitted; refit nplsqreg with the requested tau values before plotting",
           call. = FALSE)
    idx[[j]] <- hit[[1L]]
  }
  idx
}

.nplsqreg_plot_progress_target <- function(progress.target, tau.label,
                                           include.tau = TRUE) {
  if (!isTRUE(include.tau))
    return(progress.target)
  tau.label <- if (is.null(tau.label)) NULL else as.character(tau.label)[1L]
  if (is.null(tau.label) || !nzchar(tau.label) || is.na(tau.label))
    return(progress.target)
  progress.target <- if (is.null(progress.target)) NULL else as.character(progress.target)[1L]
  if (is.null(progress.target) || !nzchar(progress.target) || is.na(progress.target))
    return(tau.label)
  sprintf("%s, %s", tau.label, progress.target)
}

.np_plot_lsqregression_eval <- function(bws,
                                        txdat,
                                        tydat,
                                        exdat,
                                        tau = bws$tau,
                                        gradients = FALSE,
                                        need.errors = TRUE,
                                        gradient.order = 1L,
                                        ...) {
  .np_plot_activity_run(
    label = if (length(tau) == 1L) {
      "Computing location-scale quantile plot fit"
    } else {
      sprintf("Computing location-scale quantile plot fit (%d tau values)",
              length(tau))
    },
    expr = {
      tau <- .nplsqreg_validate_tau_values(tau, allow.duplicates = TRUE)
      gradients <- npValidateScalarLogical(gradients, "gradients")
      need.errors <- npValidateScalarLogical(need.errors, "need.errors")
      txdat <- toFrame(txdat)
      tydat <- as.numeric(toFrame(tydat)[[1L]])
      exdat <- toFrame(exdat)

      fit.one <- function(one.bws, tau.i) {
        fit <- .np_plot_regression_eval(
          bws = one.bws$reg.bws,
          xdat = txdat,
          ydat = one.bws$qdat,
          exdat = exdat,
          gradients = gradients,
          gradient.order = gradient.order,
          need.asymptotic = need.errors
        )
        lsqregression(
          bws = one.bws,
          fit = fit,
          xeval = exdat,
          tau = tau.i,
          delta = one.bws$delta,
          quantile = fit$mean,
          quanterr = if (need.errors) fit$merr else rep.int(NA_real_, length(fit$mean)),
          quantgrad = if (gradients) fit$grad else NA,
          quantgerr = if (gradients && need.errors) fit$gerr else NA,
          ntrain = nrow(txdat),
          trainiseval = FALSE,
          gradients = gradients,
          residuals = FALSE,
          resid = NA,
          call = NULL
        )
      }

      if (length(tau) == 1L) {
        idx.one <- .nplsqreg_match_plot_tau(tau, bws$tau)
        one.bws <- if (length(bws$tau) == 1L) bws else bws$tau.bws[[idx.one]]
        if (is.null(one.bws))
          stop("vector nplsqreg bandwidth object lacks per-tau bandwidth state",
               call. = FALSE)
        return(fit.one(one.bws, tau))
      }

      idx <- .nplsqreg_match_plot_tau(tau, bws$tau)
      if (is.null(bws$tau.bws) || length(bws$tau.bws) < max(idx))
        stop("vector nplsqreg bandwidth object lacks per-tau bandwidth state",
             call. = FALSE)
      fit.list <- lapply(seq_along(idx), function(j)
        fit.one(bws$tau.bws[[idx[[j]]]], tau[[j]]))
      sub.bws <- .nplsqreg_combine_bandwidths(
        bw.list = lapply(idx, function(i) bws$tau.bws[[i]]),
        tau = tau,
        tau.search = if (is.null(bws$tau.search)) "full" else bws$tau.search
      )
      .nplsqreg_combine_fits(
        fit.list = fit.list,
        tau = tau,
        bws = sub.bws,
        tau.search = if (is.null(bws$tau.search)) "full" else bws$tau.search,
        call = NULL
      )
    }
  )
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
    method = .np_plot_condbandwidth_engine,
    bws = x$bws,
    context = "plot.lsqregression"
  )
  dots <- .np_plot_normalize_public_dots(list(...),
                                         context = "plot.lsqregression")
  random.seed <- if (!is.null(dots$random.seed)) dots$random.seed else 42L
  dots$random.seed <- NULL
  if (!is.null(dots$gradient)) {
    gradients <- dots$gradient
    dots$gradient <- NULL
  }
  if (!is.null(dots$gradients)) {
    gradients <- dots$gradients
    dots$gradients <- NULL
  }
  gradients <- npValidateScalarLogical(gradients, "gradients")
  tau <- if (is.null(tau)) {
    x$tau
  } else {
    .nplsqreg_validate_tau_values(tau, allow.duplicates = TRUE)
  }

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

  if (!is.null(xlab)) dots$xlab <- xlab
  if (!is.null(ylab)) dots$ylab <- ylab
  if (!is.null(ylim)) dots$ylim <- ylim

  args <- c(
    list(
      bws = x$bws,
      xdat = x$bws$xdat,
      ydat = toFrame(x$bws$ydat),
      quantreg = TRUE,
      tau = tau,
      gradients = gradients,
      plot.data.overlay = plot.data.overlay
    ),
    dots
  )
  .np_with_seed(random.seed, do.call(.np_plot_condbandwidth_engine, args))
}

compute.bootstrap.errors.lsqregressionbandwidth <-
  function(xdat, ydat,
           exdat, eydat,
           cdf,
           quantreg,
           tau,
           gradients,
           gradient.index,
           gradient.order = 1L,
           slice.index,
           plot.errors.boot.method,
           plot.errors.boot.nonfixed = c("exact", "frozen"),
           plot.errors.boot.blocklen,
           plot.errors.boot.num,
           plot.errors.center,
           plot.errors.type,
           plot.errors.alpha,
           progress.target = NULL,
           ...,
           bws) {
    if (!isTRUE(quantreg))
      stop("location-scale quantile bootstrap requires quantreg=TRUE",
           call. = FALSE)
    tau <- .nplsqreg_validate_tau_values(tau, allow.duplicates = TRUE)
    idx <- .nplsqreg_match_plot_tau(tau, bws$tau)
    if (length(idx) == 1L && length(bws$tau) == 1L) {
      one <- bws
      return(compute.bootstrap.errors.rbandwidth(
        xdat = xdat,
        ydat = one$qdat,
        exdat = exdat,
        gradients = gradients,
        gradient.order = gradient.order,
        slice.index = slice.index,
        plot.errors.boot.method = plot.errors.boot.method,
        plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
        plot.errors.boot.blocklen = plot.errors.boot.blocklen,
        plot.errors.boot.num = plot.errors.boot.num,
        plot.errors.center = plot.errors.center,
        plot.errors.type = plot.errors.type,
        plot.errors.alpha = plot.errors.alpha,
        progress.target = .nplsqreg_plot_progress_target(
          progress.target = progress.target,
          tau.label = .npqreg_tau_labels(tau)[[1L]],
          include.tau = length(tau) > 1L || length(bws$tau) > 1L
        ),
        bws = one$reg.bws
      ))
    }

    if (is.null(bws$tau.bws) || length(bws$tau.bws) < max(idx))
      stop("vector nplsqreg bandwidth object lacks per-tau bandwidth state",
           call. = FALSE)
    tau.labels <- .npqreg_tau_labels(tau)
    one.out <- lapply(seq_along(idx), function(j) {
      one <- bws$tau.bws[[idx[[j]]]]
      compute.bootstrap.errors.rbandwidth(
        xdat = xdat,
        ydat = one$qdat,
        exdat = exdat,
        gradients = gradients,
        gradient.order = gradient.order,
        slice.index = slice.index,
        plot.errors.boot.method = plot.errors.boot.method,
        plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
        plot.errors.boot.blocklen = plot.errors.boot.blocklen,
        plot.errors.boot.num = plot.errors.boot.num,
        plot.errors.center = plot.errors.center,
        plot.errors.type = plot.errors.type,
        plot.errors.alpha = plot.errors.alpha,
        progress.target = .nplsqreg_plot_progress_target(
          progress.target = progress.target,
          tau.label = tau.labels[[j]]
        ),
        bws = one$reg.bws
      )
    })

    boot.err <- array(
      NA_real_,
      dim = c(nrow(exdat), 3L, length(tau)),
      dimnames = list(NULL, c("lower", "upper", "center"), tau.labels)
    )
    boot.all.err <- vector("list", length(tau))
    names(boot.all.err) <- tau.labels
    for (j in seq_along(one.out)) {
      boot.err[, , j] <- one.out[[j]]$boot.err
      boot.all.err[[j]] <- one.out[[j]]$boot.all.err
    }
    list(boot.err = boot.err, bxp = list(), boot.all.err = boot.all.err)
  }
