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
                              bar_num = NULL) {
  style <- match.arg(style)
  bar <- match.arg(bar)
  if (!is.null(bar_num) &&
      (!is.numeric(bar_num) || length(bar_num) != 1L ||
       is.na(bar_num) || bar_num < 1))
    stop("bar_num must be a positive numeric scalar", call. = FALSE)
  structure(
    list(style = style, bar = bar, bar_num = bar_num),
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

.np_plot_canonical_arg_names <- function() {
  c("errors", "band", "alpha", "bootstrap", "B", "center",
    "output", "data_overlay", "data_rug", "layout", "legend",
    "factor_boxplot", "boxplot_outliers", "coef_index",
    "gradient_order", "common_scale", "proper_method", "proper_control",
    "renderer", "neval", "perspective",
    "boot_control", "grid_control", "render_control")
}

.np_plot_legacy_arg_names <- function() {
  c("plot.errors.method", "plot.errors.type", "plot.errors.alpha",
    "plot.errors.boot.method", "plot.errors.boot.num",
    "plot.errors.boot.nonfixed", "plot.errors.boot.wild",
    "plot.errors.boot.blocklen", "plot.errors.center",
    "plot.errors.style", "plot.errors.bar", "plot.errors.bar.num",
    "plot.behavior", "gradient", "persp", "plot.par.mfrow",
    "plot.data.overlay", "plot.rug", "plot.bxp", "plot.bxp.out",
    "behavior")
}

.np_plot_engine_for_bws <- function(bws) {
  cls <- class(bws)
  if ("rbandwidth" %in% cls)
    return(.np_plot_rbandwidth_engine)
  if ("conbandwidth" %in% cls)
    return(.np_plot_conbandwidth_engine)
  if ("condbandwidth" %in% cls)
    return(.np_plot_condbandwidth_engine)
  if ("plbandwidth" %in% cls)
    return(.np_plot_plbandwidth_engine)
  if ("sibandwidth" %in% cls)
    return(.np_plot_sibandwidth_engine)
  if ("scbandwidth" %in% cls)
    return(.np_plot_scbandwidth_engine)
  if ("dbandwidth" %in% cls)
    return(.np_plot_dbandwidth_engine)
  if ("bandwidth" %in% cls)
    return(.np_plot_bandwidth_engine)
  NULL
}

.np_plot_allowed_engine_args <- function(method = NULL, bws = NULL) {
  if (!is.null(method) && identical(method, .np_plot_compat_dispatch))
    method <- .np_plot_engine_for_bws(bws)
  if (is.null(method))
    return(character())
  setdiff(names(formals(method)), c("bws", "..."))
}

.np_plot_validate_public_dots <- function(dots_call,
                                          method = NULL,
                                          bws = NULL,
                                          context = "plot") {
  dot.names <- .np_plot_dot_names(dots_call)
  if (any(!nzchar(dot.names)))
    stop(sprintf("unnamed plot arguments are not supported for %s", context),
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
  canonical <- .np_plot_canonical_arg_names()
  legacy <- .np_plot_legacy_arg_names()
  engine <- .np_plot_allowed_engine_args(method = method, bws = bws)
  dispatcher <- "random.seed"
  allowed <- unique(c(canonical, legacy, dispatcher, engine,
                      .np_plot_graphics_arg_names()))
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

.np_plot_match_layout <- function(value) {
  if (is.logical(value)) {
    if (length(value) != 1L || is.na(value))
      stop("layout must be TRUE/FALSE or one of \"auto\", \"current\"",
           call. = FALSE)
    return(isTRUE(value))
  }
  layout <- .np_plot_scalar_match(value, c("auto", "current"), "layout")
  identical(layout, "auto")
}

.np_plot_normalize_public_dots <- function(dots, context = "plot") {
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
  if (has("output")) {
    output <- .np_plot_scalar_match(dots$output,
                                    c("plot", "data", "both", "plot-data"),
                                    "output")
    if (identical(output, "both"))
      output <- "plot-data"
    dots$output <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "output",
                                        "plot.behavior", output)
  }
  if (has("behavior")) {
    behavior <- .np_plot_scalar_match(dots$behavior,
                                      c("plot", "plot-data", "data"),
                                      "behavior")
    dots$behavior <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "behavior",
                                        "plot.behavior", behavior)
  }
  if (has("data_overlay")) {
    data_overlay <- .np_plot_match_flag(dots$data_overlay, "data_overlay")
    dots$data_overlay <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "data_overlay",
                                        "plot.data.overlay", data_overlay)
  }
  if (has("data_rug")) {
    data_rug <- .np_plot_match_flag(dots$data_rug, "data_rug")
    dots$data_rug <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "data_rug",
                                        "plot.rug", data_rug)
  }
  if (has("layout")) {
    layout <- .np_plot_match_layout(dots$layout)
    dots$layout <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "layout",
                                        "plot.par.mfrow", layout)
  }
  if (has("factor_boxplot")) {
    factor_boxplot <- .np_plot_match_flag(dots$factor_boxplot,
                                          "factor_boxplot")
    dots$factor_boxplot <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "factor_boxplot",
                                        "plot.bxp", factor_boxplot)
  }
  if (has("boxplot_outliers")) {
    boxplot_outliers <- .np_plot_match_flag(dots$boxplot_outliers,
                                            "boxplot_outliers")
    dots$boxplot_outliers <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "boxplot_outliers",
                                        "plot.bxp.out", boxplot_outliers)
  }
  if (has("coef_index")) {
    coef_index <- dots$coef_index
    if (!is.numeric(coef_index) || length(coef_index) != 1L ||
        is.na(coef_index) || coef_index < 1L)
      stop("coef_index must be a positive numeric scalar", call. = FALSE)
    dots$coef_index <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "coef_index",
                                        "coef.index", as.integer(coef_index))
  }
  if (has("gradient_order")) {
    gradient_order <- dots$gradient_order
    if (!is.numeric(gradient_order) || any(is.na(gradient_order)) ||
        any(gradient_order < 1L))
      stop("gradient_order must contain positive numeric values",
           call. = FALSE)
    dots$gradient_order <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "gradient_order",
                                        "gradient.order",
                                        as.integer(gradient_order))
  }
  if (has("common_scale")) {
    common_scale <- .np_plot_match_flag(dots$common_scale, "common_scale")
    dots$common_scale <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "common_scale",
                                        "common.scale", common_scale)
  }
  if (has("proper_method")) {
    proper_method <- dots$proper_method
    if (!is.character(proper_method) || length(proper_method) != 1L ||
        is.na(proper_method))
      stop("proper_method must be a character scalar", call. = FALSE)
    dots$proper_method <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "proper_method",
                                        "proper.method", proper_method)
  }
  if (has("proper_control")) {
    proper_control <- dots$proper_control
    if (!is.list(proper_control))
      stop("proper_control must be a list", call. = FALSE)
    dots$proper_control <- NULL
    dots <- .np_plot_set_normalized_arg(dots, "proper_control",
                                        "proper.control", proper_control)
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
      stop(sprintf("grid_control$slices is not yet supported for %s", context),
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
    if (!is.null(ctrl$bar_num))
      dots <- .np_plot_set_normalized_arg(dots, "render_control$bar_num",
                                          "plot.errors.bar.num", ctrl$bar_num)
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

.np_plot_call_method <- function(method, bws, ..., .plot_dots_call = NULL,
                                 .plot_context = "plot") {
  .np_plot_validate_public_dots(
    .plot_dots_call,
    method = method,
    bws = bws,
    context = .plot_context
  )
  dots <- list(...)
  dots <- .np_plot_normalize_public_dots(dots, context = .plot_context)
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

  if (!is.null(dots$plot.data.overlay) &&
      !any(c("rbandwidth", "plbandwidth", "scbandwidth") %in% cls) &&
      !("condbandwidth" %in% cls && isTRUE(dots$quantreg))) {
    dots$plot.data.overlay <- .np_plot_match_flag(dots$plot.data.overlay,
                                                  "plot.data.overlay")
    if (isTRUE(dots$plot.data.overlay)) {
      stop("plot.data.overlay=TRUE is available only for regression, quantile regression, partially linear, and smooth coefficient plot surfaces.",
           call. = FALSE)
    }
    dots$plot.data.overlay <- NULL
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

.np_plot_from_slot <- function(object, slot = "bws", ...,
                               .plot_dots_call = NULL,
                               .plot_context = "plot") {
  bws <- object[[slot]]
  if (is.null(bws))
    stop("plot object does not contain expected bandwidth slot")
  .np_plot_call_method(.np_plot_compat_dispatch, bws = bws, ...,
                       .plot_dots_call = .plot_dots_call,
                       .plot_context = .plot_context)
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

.np_plot_npregression <- function(object, ..., .plot_dots_call = NULL) {
  .np_plot_validate_public_dots(
    .plot_dots_call,
    method = .np_plot_compat_dispatch,
    bws = object$bws,
    context = "plot.npregression"
  )
  dots <- list(...)
  dots <- .np_plot_normalize_public_dots(dots, context = "plot.npregression")
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
.np_plot_npdensity <- function(object, ..., .plot_dots_call = NULL) {
  .np_plot_validate_public_dots(
    .plot_dots_call,
    method = .np_plot_compat_dispatch,
    bws = object$bws,
    context = "plot.npdensity"
  )
  dots <- list(...)
  dots <- .np_plot_normalize_public_dots(dots, context = "plot.npdensity")
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

  do.call(.np_plot_from_slot, c(list(object = object, slot = "bws"), dots))
}
.np_plot_condensity <- function(object, ..., .plot_dots_call = NULL) {
  .np_plot_validate_public_dots(
    .plot_dots_call,
    method = .np_plot_compat_dispatch,
    bws = object$bws,
    context = "plot.condensity"
  )
  dots <- list(...)
  dots <- .np_plot_normalize_public_dots(dots, context = "plot.condensity")
  if (is.null(dots$proper) && isTRUE(object$proper.requested))
    dots$proper <- TRUE
  do.call(.np_plot_from_slot, c(list(object = object, slot = "bws"), dots))
}
.np_plot_condistribution <- function(object, ..., .plot_dots_call = NULL) {
  .np_plot_validate_public_dots(
    .plot_dots_call,
    method = .np_plot_compat_dispatch,
    bws = object$bws,
    context = "plot.condistribution"
  )
  dots <- list(...)
  dots <- .np_plot_normalize_public_dots(dots, context = "plot.condistribution")
  if (is.null(dots$proper) && isTRUE(object$proper.requested))
    dots$proper <- TRUE
  do.call(.np_plot_from_slot, c(list(object = object, slot = "bws"), dots))
}
.np_plot_npdistribution <- function(object, ..., .plot_dots_call = NULL)
  .np_plot_from_slot(object, "bws", ...,
                     .plot_dots_call = .plot_dots_call,
                     .plot_context = "plot.npdistribution")
.np_plot_qregression <- function(object, ..., .plot_dots_call = NULL) {
  .np_plot_validate_public_dots(
    .plot_dots_call,
    method = .np_plot_compat_dispatch,
    bws = object$bws,
    context = "plot.qregression"
  )
  dots <- list(...)
  dots <- .np_plot_normalize_public_dots(dots, context = "plot.qregression")
  if (is.null(dots$quantreg))
    dots$quantreg <- TRUE
  if (is.null(dots$tau) && !is.null(object$tau))
    dots$tau <- object$tau
  if (is.null(dots$plot.data.overlay) && !isTRUE(dots$gradients))
    dots$plot.data.overlay <- TRUE
  do.call(.np_plot_from_slot, c(list(object = object, slot = "bws"), dots))
}
.np_plot_conmode_data <- function(object, gradients = FALSE, level = NULL) {
  lev <- if (!is.null(object$probability.levels)) {
    as.character(object$probability.levels)
  } else if (!is.null(object$probability.gradient.level)) {
    as.character(object$probability.gradient.level)
  } else {
    levels(object$conmode)
  }
  lev <- unique(lev)
  if (!length(lev))
    stop("conmode object does not contain response-level metadata", call. = FALSE)

  if (is.null(level)) {
    if (isTRUE(gradients) && !is.null(object$probability.gradient.level)) {
      level <- as.character(object$probability.gradient.level)
    } else {
      level <- lev[1L]
    }
  }
  if (length(level) != 1L || is.na(level) || !(as.character(level) %in% lev))
    stop("'level' must identify one response level in the fitted conmode object",
         call. = FALSE)
  level <- as.character(level)
  level.idx <- match(level, lev)

  if (isTRUE(gradients)) {
    values <- object$probability.gradients
    if (is.null(values))
      stop("class-probability gradients/effects are not available: fit with gradients=TRUE",
           call. = FALSE)
    stored.level <- object$probability.gradient.level
    if (!identical(as.character(level), as.character(stored.level)))
      stop(sprintf("stored class-probability gradients are for level %s; refit with level=%s to plot that level",
                   sQuote(as.character(stored.level)),
                   sQuote(as.character(level))),
           call. = FALSE)
    p <- ncol(values)
    ymat <- as.matrix(values)
    value.name <- "effect"
  } else {
    values <- object$probabilities
    if (is.null(values))
      stop("class probabilities are not available: fit with probabilities=TRUE or gradients=TRUE",
           call. = FALSE)
    p <- object$xndim
    ymat <- matrix(values[, level.idx], nrow = nrow(values), ncol = 1L)
    if (p > 1L)
      ymat <- ymat[, rep(1L, p), drop = FALSE]
    value.name <- "probability"
  }

  xeval <- as.data.frame(object$xeval)
  out <- vector("list", p)
  xnames <- object$xnames
  if (is.null(xnames) || length(xnames) != p)
    xnames <- names(xeval)[seq_len(p)]
  for (j in seq_len(p)) {
    out[[j]] <- data.frame(
      variable = xnames[j],
      x = xeval[[j]],
      value = ymat[, j],
      level = level,
      gradients = isTRUE(gradients),
      stringsAsFactors = FALSE
    )
    names(out[[j]])[3L] <- value.name
  }
  names(out) <- xnames
  out
}

.np_plot_conmode_panel <- function(dat,
                                   gradients = FALSE,
                                   level,
                                   plot.user.args = list(),
                                   line.user.args = list(),
                                   data_rug = FALSE,
                                   main = NULL,
                                   xlab = NULL,
                                   ylab = NULL) {
  x <- dat$x
  y <- if (isTRUE(gradients)) dat$effect else dat$probability
  vname <- dat$variable[1L]
  if (is.null(main))
    main <- if (isTRUE(gradients)) {
      paste0("Effect on Pr(Y=", level, "|X=x)")
    } else {
      paste0("Pr(Y=", level, "|X=x)")
    }
  if (is.null(xlab))
    xlab <- vname
  if (is.null(ylab))
    ylab <- if (isTRUE(gradients)) "Effect" else "Probability"

  if (is.factor(x) && !is.ordered(x)) {
    args <- .np_plot_merge_user_args(
      list(formula = y ~ x, xlab = xlab, ylab = ylab, main = main),
      plot.user.args
    )
    do.call(graphics::boxplot, args)
    return(invisible(NULL))
  }

  xplot <- if (is.ordered(x)) as.numeric(x) else x
  ord <- order(xplot)
  args <- .np_plot_merge_user_args(
    list(x = xplot[ord], y = y[ord],
         type = if (is.ordered(x) || is.numeric(xplot)) "l" else "p",
         xlab = xlab, ylab = ylab, main = main),
    line.user.args
  )
  args <- .np_plot_merge_user_args(
    args,
    plot.user.args
  )
  if (is.ordered(x)) {
    args$xaxt <- "n"
    do.call(graphics::plot, args)
    graphics::axis(1L, at = seq_along(levels(x)), labels = levels(x))
  } else {
    do.call(graphics::plot, args)
  }
  if (isTRUE(data_rug) && is.numeric(xplot))
    graphics::rug(xplot)
  invisible(NULL)
}

.np_plot_conmode <- function(object, ...,
                             gradients = FALSE,
                             level = NULL,
                             output = c("plot", "data", "plot-data", "both"),
                             data_rug = FALSE,
                             layout = TRUE,
                             legend = TRUE,
                             .plot_dots_call = NULL) {
  dots <- list(...)
  dot.names <- .np_plot_dot_names(.plot_dots_call)
  if (any(!nzchar(dot.names)))
    stop("unnamed plot arguments are not supported for plot.conmode",
         call. = FALSE)
  interval.args <- c("errors", "band", "alpha", "bootstrap", "B", "center",
                     "boot_control", "plot.errors.method", "plot.errors.type",
                     "plot.errors.alpha", "plot.errors.boot.method",
                     "plot.errors.boot.num", "plot.errors.boot.nonfixed",
                     "plot.errors.boot.wild", "plot.errors.boot.blocklen",
                     "plot.errors.center", "plot.errors.style",
                     "plot.errors.bar", "plot.errors.bar.num")
  interval.supplied <- intersect(dot.names[nzchar(dot.names)], interval.args)
  if (length(interval.supplied))
    stop("class-probability/effect intervals are not yet implemented for plot.conmode",
         call. = FALSE)
  grid.args <- c("renderer", "perspective", "persp", "neval",
                 "grid_control", "view")
  grid.supplied <- intersect(dot.names[nzchar(dot.names)], grid.args)
  if (length(grid.supplied))
    stop("grid/surface plotting is not yet implemented for plot.conmode",
         call. = FALSE)
  allowed <- unique(c("gradients", "level", "output", "data_rug",
                      "layout", "legend",
                      .np_plot_graphics_arg_names()))
  .np_plot_stop_unused_args(setdiff(dot.names[nzchar(dot.names)], allowed),
                            allowed)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  data_rug <- npValidateScalarLogical(data_rug, "data_rug")
  plot.par.mfrow <- .np_plot_match_layout(layout)
  output <- match.arg(output)
  if (identical(output, "both"))
    output <- "plot-data"
  output <- .np_plot_scalar_match(output, c("plot", "data", "plot-data"), "output")
  plot.user.args <- .np_plot_user_args(dots, "plot")
  line.user.args <- .np_plot_user_args(dots, "lines")
  if (!is.null(dots$type))
    line.user.args$type <- dots$type

  plot.data <- .np_plot_conmode_data(object, gradients = gradients, level = level)
  level <- plot.data[[1L]]$level[1L]

  if (identical(output, "data"))
    return(plot.data)

  oldpar <- NULL
  if (isTRUE(plot.par.mfrow) && length(plot.data) > 1L) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = grDevices::n2mfrow(length(plot.data)))
  }

  for (i in seq_along(plot.data)) {
    panel.dots <- plot.user.args
    if (!is.null(dots$main) && length(plot.data) > 1L)
      panel.dots$main <- paste(dots$main, names(plot.data)[i], sep = ": ")
    .np_plot_conmode_panel(
      plot.data[[i]],
      gradients = gradients,
      level = level,
      plot.user.args = panel.dots,
      line.user.args = line.user.args,
      data_rug = data_rug
    )
  }

  if (identical(output, "plot-data"))
    return(plot.data)
  invisible(object)
}
.np_plot_singleindex <- function(object, ..., .plot_dots_call = NULL)
  .np_plot_from_slot(object, "bws", ...,
                     .plot_dots_call = .plot_dots_call,
                     .plot_context = "plot.singleindex")
.np_plot_smoothcoefficient <- function(object, ..., .plot_dots_call = NULL) {
  obj.bws <- .np_plot_restore_bandwidth_from_call(
    object = object,
    bws = object$bws,
    caller_env = parent.frame(2L)
  )
  .np_plot_validate_public_dots(
    .plot_dots_call,
    method = .np_plot_compat_dispatch,
    bws = obj.bws,
    context = "plot.smoothcoefficient"
  )
  dots <- list(...)
  dots <- .np_plot_normalize_public_dots(dots, context = "plot.smoothcoefficient")

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
.np_plot_plregression <- function(object, ..., .plot_dots_call = NULL) {
  obj.bws <- .np_plot_restore_bandwidth_from_call(
    object = object,
    bws = .np_plreg_bws(object, where = "plot.plregression"),
    caller_env = parent.frame(2L)
  )
  .np_plot_validate_public_dots(
    .plot_dots_call,
    method = .np_plot_compat_dispatch,
    bws = obj.bws,
    context = "plot.plregression"
  )
  dots <- list(...)
  dots <- .np_plot_normalize_public_dots(dots, context = "plot.plregression")

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

plot.bandwidth <- function(x, ...)
  .np_plot_call_method(.np_plot_bandwidth_engine, bws = x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...,
                       .plot_context = "plot.bandwidth")
plot.rbandwidth <- function(x, ...)
  .np_plot_call_method(.np_plot_rbandwidth_engine, bws = x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...,
                       .plot_context = "plot.rbandwidth")
plot.dbandwidth <- function(x, ...)
  .np_plot_call_method(.np_plot_dbandwidth_engine, bws = x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...,
                       .plot_context = "plot.dbandwidth")
plot.conbandwidth <- function(x, ...)
  .np_plot_call_method(.np_plot_conbandwidth_engine, bws = x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...,
                       .plot_context = "plot.conbandwidth")
plot.condbandwidth <- function(x, ...)
  .np_plot_call_method(.np_plot_condbandwidth_engine, bws = x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...,
                       .plot_context = "plot.condbandwidth")
plot.plbandwidth <- function(x, ...)
  .np_plot_call_method(.np_plot_plbandwidth_engine, bws = x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...,
                       .plot_context = "plot.plbandwidth")
plot.sibandwidth <- function(x, ...)
  .np_plot_call_method(.np_plot_sibandwidth_engine, bws = x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...,
                       .plot_context = "plot.sibandwidth")
plot.scbandwidth <- function(x, ...)
  .np_plot_call_method(.np_plot_scbandwidth_engine, bws = x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...,
                       .plot_context = "plot.scbandwidth")

plot.npregression <- function(x, ...)
  .np_plot_npregression(x, ...,
                        .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.npdensity <- function(x, ...)
  .np_plot_npdensity(x, ...,
                     .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.condensity <- function(x, ...)
  .np_plot_condensity(x, ...,
                      .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.condistribution <- function(x, ...)
  .np_plot_condistribution(x, ...,
                           .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.npdistribution <- function(x, ...)
  .np_plot_npdistribution(x, ...,
                          .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.qregression <- function(x, ...)
  .np_plot_qregression(x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.conmode <- function(x, ...)
  .np_plot_conmode(x, ...,
                   .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.singleindex <- function(x, ...)
  .np_plot_singleindex(x, ...,
                       .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.smoothcoefficient <- function(x, ...)
  .np_plot_smoothcoefficient(x, ...,
                             .plot_dots_call = match.call(expand.dots = FALSE)$...)
plot.plregression <- function(x, ...)
  .np_plot_plregression(x, ...,
                        .plot_dots_call = match.call(expand.dots = FALSE)$...)
