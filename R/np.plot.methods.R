.np_plot_scalar_match <- function(value, choices, argname) {
  if (is.null(value))
    return(NULL)
  if (length(value) != 1L || is.na(value))
    stop(sprintf("%s must be one of %s",
                 argname,
                 paste(sprintf("\"%s\"", choices), collapse = ", ")),
         call. = FALSE)
  value <- as.character(value)
  if (!(value %in% choices))
    stop(sprintf("%s must be one of %s",
                 argname,
                 paste(sprintf("\"%s\"", choices), collapse = ", ")),
         call. = FALSE)
  value
}

np_boot_control <- function(nonfixed = c("exact", "frozen"),
                            wild = c("rademacher", "mammen"),
                            blocklen = NULL) {
  nonfixed <- match.arg(nonfixed)
  wild <- match.arg(wild)
  if (!is.null(blocklen) &&
      (!is.numeric(blocklen) || length(blocklen) != 1L ||
       is.na(blocklen) || blocklen <= 0))
    stop("blocklen must be a positive numeric scalar", call. = FALSE)
  structure(
    list(nonfixed = nonfixed, wild = wild, blocklen = blocklen),
    class = "np_boot_control"
  )
}

np_grid_control <- function(xtrim = NULL, xq = NULL, slices = NULL) {
  if (!is.null(xtrim) &&
      (!is.numeric(xtrim) || length(xtrim) != 2L ||
       any(is.na(xtrim)) || any(xtrim < 0) || any(xtrim > 1) ||
       xtrim[1L] >= xtrim[2L]))
    stop("xtrim must be a numeric length-two vector with 0 <= xtrim[1] < xtrim[2] <= 1",
         call. = FALSE)
  structure(
    list(xtrim = xtrim, xq = xq, slices = slices),
    class = "np_grid_control"
  )
}

np_render_control <- function(style = c("band", "bar"),
                              bar = c("|", "I"),
                              bar.num = NULL,
                              rug = FALSE,
                              overlay = FALSE,
                              rotate = FALSE) {
  style <- match.arg(style)
  bar <- match.arg(bar)
  if (!is.null(bar.num) &&
      (!is.numeric(bar.num) || length(bar.num) != 1L ||
       is.na(bar.num) || bar.num < 1))
    stop("bar.num must be a positive numeric scalar", call. = FALSE)
  structure(
    list(style = style, bar = bar, bar.num = bar.num,
         rug = isTRUE(rug), overlay = isTRUE(overlay), rotate = isTRUE(rotate)),
    class = "np_render_control"
  )
}

.np_plot_dot_names <- function(dots_call) {
  if (is.null(dots_call) || length(dots_call) == 0L)
    return(character())
  nms <- names(dots_call)
  if (is.null(nms))
    rep.int("", length(dots_call))
  else
    nms
}

.np_plot_stop_unused_args <- function(bad, allowed) {
  bad <- unique(bad[nzchar(bad)])
  if (!length(bad))
    return(invisible(NULL))
  msg <- if (length(bad) == 1L)
    sprintf("unused plot argument: %s", bad)
  else
    sprintf("unused plot arguments: %s", paste(bad, collapse = ", "))
  close <- vapply(bad, function(x) {
    d <- utils::adist(x, allowed, partial = FALSE, ignore.case = FALSE)
    if (!length(d))
      return(NA_character_)
    i <- which.min(d)
    if (is.finite(d[i]) && d[i] <= max(2L, floor(nchar(x) / 3L)))
      allowed[i]
    else
      NA_character_
  }, character(1L))
  close <- unique(stats::na.omit(close))
  if (length(close))
    msg <- paste0(msg, "; did you mean ",
                  paste(close, collapse = " or "), "?")
  stop(msg, call. = FALSE)
}

.np_plot_graphics_arg_names <- function() {
  unique(c(
    setdiff(names(formals(graphics::plot.default)), c("x", "y", "...")),
    names(graphics::par(no.readonly = TRUE)),
    "panel.first", "panel.last", "zlab", "zlim", "theta", "phi", "border",
    "view", "type", "lty", "lwd", "col", "pch", "cex", "main", "sub",
    "xlab", "ylab", "xlim", "ylim"
  ))
}

.np_plot_validate_qregression_dots <- function(dots_call) {
  dot.names <- .np_plot_dot_names(dots_call)
  if (any(!nzchar(dot.names)))
    stop("unnamed plot arguments are not supported for plot.qregression",
         call. = FALSE)
  if ("intervals" %in% dot.names)
    stop("unused plot argument: intervals; did you mean errors?",
         call. = FALSE)
  if ("boot" %in% dot.names)
    stop("unused plot argument: boot; did you mean bootstrap?",
         call. = FALSE)
  if ("bands" %in% dot.names)
    stop("unused plot argument: bands; did you mean band?",
         call. = FALSE)
  canonical <- c("errors", "band", "alpha", "bootstrap", "B", "center",
                 "behavior", "renderer", "neval", "perspective",
                 "boot_control", "grid_control", "render_control")
  legacy <- c("plot.errors.method", "plot.errors.type", "plot.errors.alpha",
              "plot.errors.boot.method", "plot.errors.boot.num",
              "plot.errors.boot.nonfixed", "plot.errors.boot.wild",
              "plot.errors.boot.blocklen", "plot.errors.center",
              "plot.errors.style", "plot.errors.bar",
              "plot.errors.bar.num", "plot.behavior", "gradient", "persp",
              "plot.par.mfrow")
  engine <- setdiff(names(formals(.np_plot_condbandwidth_engine)), "...")
  allowed <- unique(c(canonical, legacy, engine, .np_plot_graphics_arg_names()))
  bad <- setdiff(dot.names[nzchar(dot.names)], allowed)
  .np_plot_stop_unused_args(bad, allowed)
  invisible(NULL)
}

.np_plot_set_normalized_arg <- function(dots, public, internal, value) {
  same <- !is.null(dots[[internal]]) &&
    isTRUE(all.equal(dots[[internal]], value, check.attributes = FALSE))
  if (!is.null(dots[[internal]]) && !same)
    stop(sprintf("conflicting plot arguments: %s and %s specify different values",
                 public, internal),
         call. = FALSE)
  dots[[internal]] <- value
  dots
}

.np_plot_normalize_qregression_dots <- function(dots) {
  supplied <- names(dots)
  has <- function(x) x %in% supplied

  if (has("intervals"))
    stop("unused plot argument: intervals; did you mean errors?", call. = FALSE)
  if (has("boot"))
    stop("unused plot argument: boot; did you mean bootstrap?", call. = FALSE)
  if (has("bands"))
    stop("unused plot argument: bands; did you mean band?", call. = FALSE)

  if (has("errors")) {
    errors <- .np_plot_scalar_match(dots$errors,
                                    c("none", "bootstrap", "asymptotic"),
                                    "errors")
    dots$errors <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "errors",
                                        "plot.errors.method", errors)
  }
  if (has("band")) {
    band <- .np_plot_scalar_match(dots$band,
                                  c("pmzsd", "pointwise", "bonferroni",
                                    "simultaneous", "all"),
                                  "band")
    dots$band <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "band",
                                        "plot.errors.type", band)
  }
  if (has("alpha")) {
    alpha <- dots$alpha
    if (!is.numeric(alpha) || length(alpha) != 1L ||
        is.na(alpha) || alpha <= 0 || alpha >= 0.5)
      stop("alpha must lie in (0, 0.5)", call. = FALSE)
    dots$alpha <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "alpha",
                                        "plot.errors.alpha", alpha)
  }
  if (has("bootstrap")) {
    if (is.list(dots$bootstrap))
      stop("unused plot argument: bootstrap list; use scalar bootstrap, B, and boot_control",
           call. = FALSE)
    bootstrap <- .np_plot_scalar_match(dots$bootstrap,
                                       c("wild", "inid", "fixed", "geom"),
                                       "bootstrap")
    dots$bootstrap <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "bootstrap",
                                        "plot.errors.boot.method", bootstrap)
  }
  if (has("B")) {
    B <- dots$B
    if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1)
      stop("B must be a positive numeric scalar", call. = FALSE)
    dots$B <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "B",
                                        "plot.errors.boot.num", as.integer(B))
  }
  if (has("center")) {
    center <- .np_plot_scalar_match(dots$center,
                                    c("estimate", "bias-corrected"),
                                    "center")
    dots$center <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "center",
                                        "plot.errors.center", center)
  }
  if (has("behavior")) {
    behavior <- .np_plot_scalar_match(dots$behavior,
                                      c("plot", "plot-data", "data"),
                                      "behavior")
    dots$behavior <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "behavior",
                                        "plot.behavior", behavior)
  }
  if (has("perspective")) {
    perspective <- isTRUE(dots$perspective)
    dots$perspective <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "perspective",
                                        "perspective", perspective)
  }
  if (has("boot_control")) {
    if (!inherits(dots$boot_control, "np_boot_control"))
      stop("boot_control must be created by np_boot_control()", call. = FALSE)
    ctrl <- dots$boot_control
    dots$boot_control <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "boot_control$nonfixed",
                                        "plot.errors.boot.nonfixed",
                                        ctrl$nonfixed)
    dots <- .np_plot_set_normalized_arg(dots, "boot_control$wild",
                                        "plot.errors.boot.wild",
                                        ctrl$wild)
    if (!is.null(ctrl$blocklen))
      dots <- .np_plot_set_normalized_arg(dots, "boot_control$blocklen",
                                          "plot.errors.boot.blocklen",
                                          ctrl$blocklen)
  }
  if (has("grid_control")) {
    if (!inherits(dots$grid_control, "np_grid_control"))
      stop("grid_control must be created by np_grid_control()", call. = FALSE)
    ctrl <- dots$grid_control
    dots$grid_control <- NULL
    if (!is.null(ctrl$xtrim))
      dots <- .np_plot_set_normalized_arg(dots, "grid_control$xtrim",
                                          "xtrim", ctrl$xtrim)
    if (!is.null(ctrl$xq))
      dots <- .np_plot_set_normalized_arg(dots, "grid_control$xq",
                                          "xq", ctrl$xq)
    if (!is.null(ctrl$slices))
      stop("grid_control$slices is not yet supported for plot.qregression",
           call. = FALSE)
  }
  if (has("render_control")) {
    if (!inherits(dots$render_control, "np_render_control"))
      stop("render_control must be created by np_render_control()", call. = FALSE)
    ctrl <- dots$render_control
    dots$render_control <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "render_control$style",
                                        "plot.errors.style", ctrl$style)
    dots <- .np_plot_set_normalized_arg(dots, "render_control$bar",
                                        "plot.errors.bar", ctrl$bar)
    if (!is.null(ctrl$bar.num))
      dots <- .np_plot_set_normalized_arg(dots, "render_control$bar.num",
                                          "plot.errors.bar.num", ctrl$bar.num)
    if (isTRUE(ctrl$rug))
      dots <- .np_plot_set_normalized_arg(dots, "render_control$rug",
                                          "plot.rug", TRUE)
    if (isTRUE(ctrl$overlay))
      stop("render_control$overlay is not yet supported for plot.qregression",
           call. = FALSE)
    if (isTRUE(ctrl$rotate))
      stop("render_control$rotate is not yet supported for plot.qregression",
           call. = FALSE)
  }

  method <- if (!is.null(dots$plot.errors.method))
    as.character(dots$plot.errors.method)[1L]
  else
    "none"
  boot.only <- c("plot.errors.boot.method", "plot.errors.boot.num",
                 "plot.errors.boot.nonfixed", "plot.errors.boot.wild",
                 "plot.errors.boot.blocklen")
  boot.supplied <- any(c("bootstrap", "B", "boot_control", "center",
                         "plot.errors.center", boot.only) %in% supplied)
  if (!identical(method, "bootstrap") && boot.supplied)
    stop("bootstrap controls require errors = \"bootstrap\"",
         call. = FALSE)
  error.only <- c("plot.errors.type", "plot.errors.alpha",
                  "plot.errors.style", "plot.errors.bar",
                  "plot.errors.bar.num")
  error.supplied <- any(c("band", "alpha", "render_control",
                          error.only) %in% supplied)
  if (identical(method, "none") && error.supplied)
    stop("band, alpha, and interval rendering controls require errors != \"none\"",
         call. = FALSE)

  dots
}

.np_plot_call_method <- function(method, bws, ...) {
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

  .np_with_seed(random.seed, do.call(method, c(list(bws = bws), dots)))
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
  .np_plot_call_method(.np_plot_compat_dispatch, bws = bws, ...)
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
      length(object$mean) == length(object$resid)) {
    dots$xdat <- object$eval
    dots$ydat <- object$mean + object$resid
  }

  do.call(.np_plot_from_slot, c(list(object = object, slot = "bws"), dots))
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
.np_plot_qregression <- function(object, ..., .plot_dots_call = NULL) {
  .np_plot_validate_qregression_dots(.plot_dots_call)
  dots <- list(...)
  dots <- .np_plot_normalize_qregression_dots(dots)
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
          c(list(method = .np_plot_compat_dispatch, bws = obj.bws), dots))
}

plot.bandwidth <- function(x, ...) .np_plot_call_method(.np_plot_bandwidth_engine, bws = x, ...)
plot.rbandwidth <- function(x, ...) .np_plot_call_method(.np_plot_rbandwidth_engine, bws = x, ...)
plot.dbandwidth <- function(x, ...) .np_plot_call_method(.np_plot_dbandwidth_engine, bws = x, ...)
plot.conbandwidth <- function(x, ...) .np_plot_call_method(.np_plot_conbandwidth_engine, bws = x, ...)
plot.condbandwidth <- function(x, ...) .np_plot_call_method(.np_plot_condbandwidth_engine, bws = x, ...)
plot.plbandwidth <- function(x, ...) .np_plot_call_method(.np_plot_plbandwidth_engine, bws = x, ...)
plot.sibandwidth <- function(x, ...) .np_plot_call_method(.np_plot_sibandwidth_engine, bws = x, ...)
plot.scbandwidth <- function(x, ...) .np_plot_call_method(.np_plot_scbandwidth_engine, bws = x, ...)

plot.npregression <- function(x, ...) .np_plot_npregression(x, ...)
plot.npdensity <- function(x, ...) .np_plot_npdensity(x, ...)
plot.condensity <- function(x, ...) .np_plot_condensity(x, ...)
plot.condistribution <- function(x, ...) .np_plot_condistribution(x, ...)
plot.npdistribution <- function(x, ...) .np_plot_npdistribution(x, ...)
plot.qregression <- function(x, ...)
  .np_plot_qregression(x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.singleindex <- function(x, ...) .np_plot_singleindex(x, ...)
plot.smoothcoefficient <- function(x, ...) .np_plot_smoothcoefficient(x, ...)
plot.plregression <- function(x, ...) .np_plot_plregression(x, ...)
