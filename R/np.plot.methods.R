.np_plot_call_method <- function(method, bws, ..., where = "plot()") {
  dots <- list(...)
  random.seed <- if (!is.null(dots$random.seed)) dots$random.seed else 42L
  dots$random.seed <- NULL

  # Keep backward compatibility with legacy plot argument aliases.
  if (!is.null(dots$gradient) && is.null(dots$gradients)) {
    dots$gradients <- dots$gradient
  }
  dots$gradient <- NULL

  if (!is.null(dots$persp) && is.null(dots$perspective)) {
    dots$perspective <- dots$persp
  }
  dots$persp <- NULL

  if (!is.null(dots$plot.rug)) {
    dots$plot.rug <- .np_plot_match_flag(dots$plot.rug, "plot.rug")
    if (isTRUE(dots$plot.rug) &&
        !identical(method, .np_plot_rbandwidth_engine) &&
        !identical(method, .np_plot_bandwidth_engine) &&
        !identical(method, .np_plot_dbandwidth_engine) &&
        !identical(method, .np_plot_conbandwidth_engine) &&
        !identical(method, .np_plot_condbandwidth_engine) &&
        !identical(method, .np_plot_plbandwidth_engine) &&
        !identical(method, .np_plot_scbandwidth_engine) &&
        !identical(method, .np_plot_sibandwidth_engine) &&
        !identical(method, .np_plot_compat_dispatch)) {
      stop("plot.rug=TRUE is not yet implemented for this plot route.",
           call. = FALSE)
    }
  }

  if (!is.null(dots$renderer)) {
    dots$renderer <- .np_plot_match_renderer(dots$renderer)
    if (identical(dots$renderer, "rgl") &&
        !identical(method, .np_plot_rbandwidth_engine) &&
        !identical(method, .np_plot_bandwidth_engine) &&
        !identical(method, .np_plot_dbandwidth_engine) &&
        !identical(method, .np_plot_conbandwidth_engine) &&
        !identical(method, .np_plot_condbandwidth_engine) &&
        !identical(method, .np_plot_scbandwidth_engine) &&
        !identical(method, .np_plot_plbandwidth_engine) &&
        !identical(method, .np_plot_compat_dispatch)) {
      stop("renderer='rgl' is not yet implemented for this plot route. Use renderer='base'.",
           call. = FALSE)
    }
  }

  .npRmpi_require_active_slave_pool(where = where)
  .npRmpi_guard_no_auto_object_in_manual_bcast(bws, where = where)

  if (isTRUE(.npRmpi_autodispatch_called_from_bcast())) {
    req.method <- if (!is.null(dots$plot.errors.method) &&
                      length(dots$plot.errors.method) >= 1L &&
                      !is.na(dots$plot.errors.method[1L])) {
      as.character(dots$plot.errors.method[1L])
    } else {
      NA_character_
    }
    if (identical(req.method, "bootstrap")) {
      stop(
        sprintf(
          "%s cannot run plot.errors.method='bootstrap' inside mpi.bcast.cmd(...); call plot(...) from master context with npRmpi.autodispatch=TRUE",
          where
        ),
        call. = FALSE
      )
    }
  }

  .np_with_seed(
    random.seed,
    do.call(method, c(list(bws = .npRmpi_autodispatch_untag(bws)), dots))
  )
}

.np_plot_compat_dispatch <- function(bws, ...) {
  cls <- class(bws)
  dots <- list(...)

  if (!is.null(dots$renderer)) {
    dots$renderer <- .np_plot_match_renderer(dots$renderer)
    if (identical(dots$renderer, "rgl") &&
        !any(c("rbandwidth", "bandwidth", "dbandwidth",
               "conbandwidth", "condbandwidth", "scbandwidth", "plbandwidth") %in% cls)) {
      stop("renderer='rgl' is not yet implemented for this plot route. Use renderer='base'.",
           call. = FALSE)
    }
  }

  if (!is.null(dots$plot.rug)) {
    dots$plot.rug <- .np_plot_match_flag(dots$plot.rug, "plot.rug")
    if (isTRUE(dots$plot.rug) &&
        !any(c("rbandwidth", "bandwidth", "dbandwidth",
               "conbandwidth", "condbandwidth",
               "plbandwidth", "scbandwidth", "sibandwidth") %in% cls)) {
      stop("plot.rug=TRUE is not yet implemented for this plot route.",
           call. = FALSE)
    }
  }

  if ("rbandwidth" %in% cls)
    return(do.call(.np_plot_rbandwidth_engine, c(list(bws = bws), dots)))
  if ("conbandwidth" %in% cls)
    return(do.call(.np_plot_conbandwidth_engine, c(list(bws = bws), dots)))
  if ("condbandwidth" %in% cls)
    return(do.call(.np_plot_condbandwidth_engine, c(list(bws = bws), dots)))
  if ("plbandwidth" %in% cls)
    return(do.call(.np_plot_plbandwidth_engine, c(list(bws = bws), dots)))
  if ("sibandwidth" %in% cls)
    return(do.call(.np_plot_sibandwidth_engine, c(list(bws = bws), dots)))
  if ("scbandwidth" %in% cls)
    return(do.call(.np_plot_scbandwidth_engine, c(list(bws = bws), dots)))
  if ("dbandwidth" %in% cls)
    return(do.call(.np_plot_dbandwidth_engine, c(list(bws = bws), dots)))
  if ("bandwidth" %in% cls)
    return(do.call(.np_plot_bandwidth_engine, c(list(bws = bws), dots)))

  stop("unsupported bandwidth class for plotting")
}

.np_plot_from_slot <- function(object, slot = "bws", ...) {
  bws <- object[[slot]]
  if (is.null(bws))
    stop("plot object does not contain expected bandwidth slot")
  .np_plot_call_method(.np_plot_compat_dispatch, bws = bws, ...,
                       where = "plot()")
}

.np_plot_restore_bandwidth_from_call <- function(object, bws, caller_env = parent.frame()) {
  if (!is.null(bws$formula) || is.null(object$call))
    return(bws)

  bws.orig <- tryCatch(
    .np_eval_call_arg(object$call, "bws", caller_env = caller_env),
    error = function(e) NULL
  )

  if (!is.null(bws.orig) && any(grepl("bandwidth$", class(bws.orig))))
    return(bws.orig)

  bws
}

.np_plot_plreg_training_data <- function(bws) {
  out <- list(xdat = NULL, ydat = NULL, zdat = NULL)

  if (is.null(bws))
    return(out)

  if (!is.null(bws$formula) && !is.null(bws$call)) {
    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)

    tmf.xf <- tmf <- bws$call[c(1, m)]
    tmf[[1]] <- tmf.xf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt

    mf.args <- as.list(tmf)[-1L]
    tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    bronze <- lapply(bws$chromoly, paste, collapse = " + ")
    tmf.xf[["formula"]] <- as.formula(paste(" ~ ", bronze[[2]]),
                                      env = environment(tt))
    mf.xf.args <- as.list(tmf.xf)[-1L]
    tmf.xf <- do.call(stats::model.frame, mf.xf.args, envir = environment(tt))

    out$ydat <- model.response(tmf)
    out$xdat <- tmf.xf
    out$zdat <- tmf[, bws$chromoly[[3]], drop = FALSE]
    return(out)
  }

  if (!is.null(bws$call)) {
    out$xdat <- tryCatch(data.frame(.np_eval_bws_call_arg(bws, "xdat")),
                         error = function(e) NULL)
    out$ydat <- tryCatch(.np_eval_bws_call_arg(bws, "ydat"),
                         error = function(e) NULL)
    out$zdat <- tryCatch(data.frame(.np_eval_bws_call_arg(bws, "zdat")),
                         error = function(e) NULL)
  }

  out
}

.np_plot_npregression <- function(object, ...) {
  dots <- list(...)
  if (is.null(dots$xdat) && is.null(dots$ydat) &&
      is.null(object$bws$formula) &&
      !is.null(object$call)) {
    bws.orig <- tryCatch(.np_eval_call_arg(object$call, "bws", caller_env = parent.frame(2L)),
                         error = function(e) NULL)
    if (!is.null(bws.orig) && any(grepl("bandwidth$", class(bws.orig))))
      object$bws <- bws.orig
  }

  if (is.null(dots$xdat) && is.null(dots$ydat) &&
      isTRUE(object$trainiseval) &&
      !is.null(object$eval) &&
      !is.null(object$mean) &&
      !is.null(object$resid) &&
      NROW(object$eval) == length(object$mean) &&
      length(object$mean) == length(object$resid) &&
      all(is.finite(object$mean)) &&
      all(is.finite(object$resid))) {
    dots$xdat <- object$eval
    dots$ydat <- object$mean + object$resid
    dots$fit.mean.train <- object$mean
  }

  plot.args <- c(list(object = object, slot = "bws"), dots)
  use.local.plot <- !.npRmpi_has_active_slave_pool(comm = 1L) &&
    identical(as.character(dots$plot.behavior)[1L], "data")

  if (use.local.plot)
    return(.npRmpi_with_local_regression(do.call(.np_plot_from_slot, plot.args)))

  do.call(.np_plot_from_slot, plot.args)
}
.np_plot_npdensity <- function(object, ...) {
  dots <- list(...)
  plot.rug <- FALSE
  if (!is.null(dots$plot.rug)) {
    plot.rug <- .np_plot_match_flag(dots$plot.rug, "plot.rug")
    dots$plot.rug <- NULL
  }

  direct.args <- c("plot.behavior", "plot.errors.method", "plot.errors.type",
                   "plot.errors.boot.num", "plot.errors.boot.method",
                   "plot.errors.boot.nonfixed",
                   "plot.errors.alpha", "perspective", "gradients",
                   "xdat", "data", "neval", "xtrim", "xq")
  use.direct <- isTRUE(object$ndim == 1) &&
    isTRUE(object$trainiseval) &&
    !any(direct.args %in% names(dots))

  if (use.direct) {
    ex <- object$eval[[1]]
    if (!is.factor(ex)) {
      ord <- order(ex)
      xlab.val <- if (!is.null(dots$xlab)) dots$xlab else gen.label(object$xnames[1], "X1")
      ylab.val <- if (!is.null(dots$ylab)) dots$ylab else "Density"
      type.val <- if (!is.null(dots$type)) dots$type else "l"

      dots$xlab <- NULL
      dots$ylab <- NULL
      dots$type <- NULL

      do.call(plot, c(list(x = as.numeric(ex[ord]),
                           y = as.numeric(object$dens[ord]),
                           type = type.val,
                           xlab = xlab.val,
                           ylab = ylab.val),
                      dots))
      if (isTRUE(plot.rug)) {
        .np_plot_validate_rug_request(
          plot.rug = TRUE,
          route = "plot.npdensity()",
          supported.route = TRUE,
          renderer = "base"
        )
        .np_plot_draw_rug_1d(as.numeric(ex))
      }
      return(invisible(object))
    }
  }

  .np_plot_from_slot(object, "bws", ...)
}
.np_plot_condensity <- function(object, ...) {
  dots <- list(...)
  if (is.null(dots$proper) && isTRUE(object$proper.requested))
    dots$proper <- TRUE
  do.call(.np_plot_from_slot, c(list(object = object, slot = "bws"), dots))
}
.np_plot_condistribution <- function(object, ...) {
  dots <- list(...)
  if (is.null(dots$proper) && isTRUE(object$proper.requested))
    dots$proper <- TRUE
  do.call(.np_plot_from_slot, c(list(object = object, slot = "bws"), dots))
}
.np_plot_npdistribution <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_qregression <- function(object, ...) {
  dots <- list(...)
  if (is.null(dots$quantreg))
    dots$quantreg <- TRUE
  if (is.null(dots$tau) && !is.null(object$tau))
    dots$tau <- object$tau
  do.call(.np_plot_from_slot, c(list(object = object, slot = "bws"), dots))
}
.np_plot_singleindex <- function(object, ...) .np_plot_from_slot(object, "bws", ...)
.np_plot_smoothcoefficient <- function(object, ...) {
  dots <- list(...)
  obj.bws <- .np_plot_restore_bandwidth_from_call(
    object = object,
    bws = object$bws,
    caller_env = parent.frame(2L)
  )

  if (is.null(dots$xdat) && is.null(dots$ydat) && is.null(dots$zdat) &&
      isTRUE(object$trainiseval) &&
      !is.null(object$eval) &&
      !is.null(object$mean) &&
      !is.null(object$resid) &&
      length(object$mean) == length(object$resid) &&
      all(is.finite(object$mean)) &&
      all(is.finite(object$resid))) {
    if (is.list(object$eval) && !is.null(object$eval$exdat)) {
      dots$xdat <- object$eval$exdat
      if (!is.null(object$eval$ezdat))
        dots$zdat <- object$eval$ezdat
    } else {
      dots$xdat <- object$eval
    }
    dots$ydat <- object$mean + object$resid
  }

  do.call(.np_plot_call_method,
          c(list(method = .np_plot_compat_dispatch, bws = obj.bws), dots))
}
.np_plot_plregression <- function(object, ...) {
  dots <- list(...)
  obj.bws <- .np_plot_restore_bandwidth_from_call(
    object = object,
    bws = .np_plreg_bws(object, where = "plot.plregression"),
    caller_env = parent.frame(2L)
  )

  if (is.null(dots$xdat) && is.null(dots$ydat) && is.null(dots$zdat) &&
      isTRUE(object$trainiseval) &&
      !is.null(object$evalx) &&
      !is.null(object$evalz) &&
      !is.null(object$mean) &&
      !is.null(object$resid) &&
      NROW(object$evalx) == NROW(object$evalz) &&
      NROW(object$evalx) == length(object$mean) &&
      length(object$mean) == length(object$resid) &&
      all(is.finite(object$mean)) &&
      all(is.finite(object$resid))) {
    dots$xdat <- object$evalx
    dots$ydat <- object$mean + object$resid
    dots$zdat <- object$evalz
  }

  if (is.null(dots$xdat) && is.null(dots$ydat) && is.null(dots$zdat)) {
    training <- .np_plot_plreg_training_data(obj.bws)
    if (!is.null(training$xdat) &&
        !is.null(training$ydat) &&
        !is.null(training$zdat)) {
      dots$xdat <- training$xdat
      dots$ydat <- training$ydat
      dots$zdat <- training$zdat
    }
  }

  do.call(.np_plot_call_method,
          c(list(method = .np_plot_compat_dispatch,
                 bws = obj.bws,
                 where = "plot.plregression()"),
            dots))
}

.np_plot_bandwidth <- function(bws, ...) {
  .np_plot_call_method(.np_plot_bandwidth_engine, bws = bws, ..., where = "plot.bandwidth()")
}
.np_plot_rbandwidth <- function(bws, ...) {
  .np_plot_call_method(.np_plot_rbandwidth_engine, bws = bws, ..., where = "plot.rbandwidth()")
}
.np_plot_dbandwidth <- function(bws, ...) {
  .np_plot_call_method(.np_plot_dbandwidth_engine, bws = bws, ..., where = "plot.dbandwidth()")
}
.np_plot_conbandwidth <- function(bws, ...) {
  .np_plot_call_method(.np_plot_conbandwidth_engine, bws = bws, ..., where = "plot.conbandwidth()")
}
.np_plot_condbandwidth <- function(bws, ...) {
  .np_plot_call_method(.np_plot_condbandwidth_engine, bws = bws, ..., where = "plot.condbandwidth()")
}
.np_plot_plbandwidth <- function(bws, ...) {
  .np_plot_call_method(.np_plot_plbandwidth_engine, bws = bws, ..., where = "plot.plbandwidth()")
}
.np_plot_sibandwidth <- function(bws, ...) {
  .np_plot_call_method(.np_plot_sibandwidth_engine, bws = bws, ..., where = "plot.sibandwidth()")
}
.np_plot_scbandwidth <- function(bws, ...) {
  .np_plot_call_method(.np_plot_scbandwidth_engine, bws = bws, ..., where = "plot.scbandwidth()")
}

plot.bandwidth <- function(x, ...) .np_plot_bandwidth(x, ...)
plot.rbandwidth <- function(x, ...) .np_plot_rbandwidth(x, ...)
plot.dbandwidth <- function(x, ...) .np_plot_dbandwidth(x, ...)
plot.conbandwidth <- function(x, ...) .np_plot_conbandwidth(x, ...)
plot.condbandwidth <- function(x, ...) .np_plot_condbandwidth(x, ...)
plot.plbandwidth <- function(x, ...) .np_plot_plbandwidth(x, ...)
plot.sibandwidth <- function(x, ...) .np_plot_sibandwidth(x, ...)
plot.scbandwidth <- function(x, ...) .np_plot_scbandwidth(x, ...)

plot.npregression <- function(x, ...) .np_plot_npregression(x, ...)
plot.npdensity <- function(x, ...) .np_plot_npdensity(x, ...)
plot.condensity <- function(x, ...) .np_plot_condensity(x, ...)
plot.condistribution <- function(x, ...) .np_plot_condistribution(x, ...)
plot.npdistribution <- function(x, ...) .np_plot_npdistribution(x, ...)
plot.qregression <- function(x, ...) .np_plot_qregression(x, ...)
plot.singleindex <- function(x, ...) .np_plot_singleindex(x, ...)
plot.smoothcoefficient <- function(x, ...) .np_plot_smoothcoefficient(x, ...)
plot.plregression <- function(x, ...) .np_plot_plregression(x, ...)
