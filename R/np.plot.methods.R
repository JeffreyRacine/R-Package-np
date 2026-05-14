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
.np_plot_conmode_data <- function(object,
                                  gradients = FALSE,
                                  level = NULL,
                                  errors = "none",
                                  alpha = 0.05,
                                  band = "pointwise",
                                  center = "estimate",
                                  bootstrap = "inid",
                                  B = 399L,
                                  blocklen = NULL) {
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
  boot.raw <- NULL
  boot.se <- NULL
  if (!isTRUE(gradients) && identical(errors, "bootstrap")) {
    training <- .np_plot_conmode_training_data(object)
    boot.raw <- .np_plot_conmode_bootstrap_raw(
      object = object,
      xtrain = training$xtrain,
      ytrain = training$ytrain,
      exdat = xeval,
      level = level,
      levels = lev,
      t0 = ymat[, 1L],
      alpha = alpha,
      band = band,
      center = center,
      method = bootstrap,
      B = B,
      blocklen = blocklen
    )
    boot.se <- rep(NA_real_, nrow(xeval))
  }
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
    if (!isTRUE(gradients) && !identical(errors, "none")) {
      se <- rep(NA_real_, nrow(out[[j]]))
      if (identical(errors, "asymptotic")) {
        se <- object$probability.errors
        if (is.null(se))
          stop("class-probability standard errors are not available: refit with probabilities=TRUE",
               call. = FALSE)
        se <- as.matrix(se)[, level.idx]
        repaired <- object$probability.repaired.rows
        if (!is.null(repaired)) {
          repaired <- as.logical(repaired)
          se[is.na(repaired) | repaired] <- NA_real_
        }
      }
      out[[j]] <- .np_plot_conmode_add_interval_columns(
        out[[j]],
        se = if (identical(errors, "bootstrap")) boot.se else se,
        errors = errors,
        alpha = alpha,
        band = band,
        center = center,
        bootstrap_raw = boot.raw
      )
    }
  }
  names(out) <- xnames
  out
}

.np_plot_conmode_cast_like <- function(x, template) {
  if (is.factor(template))
    return(factor(as.character(x),
                  levels = levels(template),
                  ordered = is.ordered(template)))
  if (is.integer(template))
    return(as.integer(x))
  if (is.numeric(template))
    return(as.numeric(x))
  x
}

.np_plot_conmode_grid_values <- function(x, neval, xtrim) {
  if (is.factor(x)) {
    return(factor(levels(x), levels = levels(x), ordered = is.ordered(x)))
  }
  rng <- stats::quantile(x, probs = xtrim, names = FALSE, na.rm = TRUE)
  seq(rng[1L], rng[2L], length.out = neval)
}

.np_plot_conmode_base_row <- function(xtrain, xq) {
  out <- xtrain[1L, , drop = FALSE]
  for (j in seq_along(xtrain)) {
    val <- uocquantile(xtrain[[j]], prob = xq[j])
    out[[j]] <- .np_plot_conmode_cast_like(val, xtrain[[j]])
  }
  out
}

.np_plot_conmode_level_factor <- function(ytrain, level, n) {
  y <- ytrain[[1L]]
  factor(rep(level, n),
         levels = levels(y),
         ordered = is.ordered(y))
}

.np_plot_conmode_training_data <- function(object) {
  xtrain <- object$xtrain
  ytrain <- object$ytrain
  if ((is.null(xtrain) || is.null(ytrain)) &&
      isTRUE(object$trainiseval) &&
      !is.null(object$xeval) &&
      !is.null(object$yeval)) {
    xtrain <- object$xeval
    ytrain <- object$yeval
  }
  if (is.null(xtrain) || is.null(ytrain))
    stop("fixed-grid plot.conmode displays require an object fitted with probabilities=TRUE or gradients=TRUE",
         call. = FALSE)
  list(xtrain = as.data.frame(xtrain), ytrain = as.data.frame(ytrain))
}

.np_plot_conmode_levels <- function(object, ytrain) {
  lev <- if (!is.null(object$probability.levels)) {
    as.character(object$probability.levels)
  } else {
    levels(ytrain[[1L]])
  }
  lev <- unique(lev)
  if (!length(lev))
    stop("conmode object does not contain response-level metadata", call. = FALSE)
  lev
}

.np_plot_conmode_resolve_level <- function(object, ytrain, gradients, level) {
  lev <- .np_plot_conmode_levels(object, ytrain)
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
  if (isTRUE(gradients) &&
      !identical(level, as.character(object$probability.gradient.level)))
    stop(sprintf("stored class-probability gradients are for level %s; refit with level=%s to plot that level",
                 sQuote(as.character(object$probability.gradient.level)),
                 sQuote(level)),
         call. = FALSE)
  list(level = level, levels = lev)
}

.np_plot_conmode_resolve_neval <- function(neval) {
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("neval must be a numeric scalar at least 2 for fixed-grid plot.conmode displays",
         call. = FALSE)
  as.integer(neval)
}

.np_plot_conmode_resolve_xtrim <- function(xtrim) {
  if (!is.numeric(xtrim) || length(xtrim) != 2L ||
      any(is.na(xtrim)) || any(xtrim < 0) || any(xtrim > 1) ||
      xtrim[1L] >= xtrim[2L])
    stop("xtrim must be a numeric length-two vector with 0 <= xtrim[1] < xtrim[2] <= 1",
         call. = FALSE)
  xtrim
}

.np_plot_conmode_resolve_xq <- function(xq, p) {
  if (is.null(xq))
    xq <- rep(0.5, p)
  if (!is.numeric(xq) || any(is.na(xq)) || any(xq < 0) || any(xq > 1))
    stop("xq must contain probabilities in [0,1]", call. = FALSE)
  if (length(xq) == 1L)
    xq <- rep(xq, p)
  if (length(xq) != p)
    stop("xq must be a scalar or have one entry for each conditioning variable",
         call. = FALSE)
  xq
}

.np_plot_conmode_probability_se <- function(object,
                                            xtrain,
                                            ytrain,
                                            exdat,
                                            level) {
  dens.obj <- npcdens(
    txdat = xtrain,
    tydat = ytrain,
    exdat = exdat,
    eydat = .np_plot_conmode_level_factor(ytrain, level, nrow(exdat)),
    bws = object$bws
  )
  as.vector(dens.obj$conderr)
}

.np_plot_conmode_proper_control <- function(object) {
  out <- list()
  if (!is.null(object$proper.info$tol))
    out$tol <- object$proper.info$tol
  out
}

.np_plot_conmode_probability_matrix <- function(object,
                                                xtrain,
                                                ytrain,
                                                exdat,
                                                levels) {
  rhs <- rep.int(1.0, nrow(xtrain))
  pmat <- matrix(NA_real_, nrow(exdat), length(levels),
                 dimnames = list(NULL, levels))
  for (k in seq_along(levels)) {
    pmat[, k] <- npcdenshat(
      bws = object$bws,
      txdat = xtrain,
      tydat = ytrain,
      exdat = exdat,
      eydat = .np_plot_conmode_level_factor(ytrain, levels[k], nrow(exdat)),
      y = rhs,
      output = "apply"
    )
  }
  .npConmodeProperProbabilities(
    pmat,
    levels = levels,
    proper = isTRUE(object$proper.requested),
    proper.control = .np_plot_conmode_proper_control(object)
  )
}

.np_plot_conmode_boot_counts <- function(n, B, method, blocklen) {
  if (identical(method, "inid"))
    return(.np_inid_counts_matrix(n = n, B = B))
  drawer <- .np_block_counts_drawer(
    n = n,
    sim = method,
    blocklen = blocklen,
    B = B
  )
  .np_inid_counts_matrix(n = n, B = B, counts = drawer(1L, B))
}

.np_plot_conmode_bootstrap_values <- function(object,
                                              xtrain,
                                              ytrain,
                                              exdat,
                                              level,
                                              levels,
                                              t0,
                                              B,
                                              method,
                                              blocklen) {
  if (identical(method, "wild"))
    stop("wild bootstrap is not available for class-probability plot.conmode intervals; use bootstrap=\"inid\", \"fixed\", or \"geom\"",
         call. = FALSE)
  B <- as.integer(B)
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1L)
    stop("plot.errors.boot.num must be a positive integer", call. = FALSE)
  n <- nrow(xtrain)
  counts <- .np_plot_conmode_boot_counts(
    n = n,
    B = B,
    method = method,
    blocklen = blocklen
  )
  level.idx <- match(level, levels)
  tmat <- matrix(NA_real_, nrow = B, ncol = length(t0))
  progress <- .np_plot_bootstrap_progress_begin(
    total = B,
    label = "Conmode probability bootstrap"
  )
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (bb in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, bb])
    proper.out <- .np_plot_conmode_probability_matrix(
      object = object,
      xtrain = xtrain[idx, , drop = FALSE],
      ytrain = ytrain[idx, , drop = FALSE],
      exdat = exdat,
      levels = levels
    )
    tmat[bb, ] <- proper.out$probabilities[, level.idx]
    progress <- .np_plot_progress_tick(progress, bb)
  }
  list(t = tmat, t0 = t0)
}

.np_plot_conmode_bootstrap_raw <- function(object,
                                           xtrain,
                                           ytrain,
                                           exdat,
                                           level,
                                           levels,
                                           t0,
                                           alpha,
                                           band,
                                           center,
                                           method,
                                           B,
                                           blocklen) {
  boot.out <- .np_plot_conmode_bootstrap_values(
    object = object,
    xtrain = xtrain,
    ytrain = ytrain,
    exdat = exdat,
    level = level,
    levels = levels,
    t0 = t0,
    B = B,
    method = method,
    blocklen = blocklen
  )
  if (identical(band, "pmzsd")) {
    boot.err <- matrix(NA_real_, nrow = length(t0), ncol = 3L)
    boot.sd <- .np_plot_bootstrap_col_sds(boot.out$t)
    boot.err[, 1:2] <- stats::qnorm(alpha / 2, lower.tail = FALSE) * boot.sd
    boot.all.err <- NULL
  } else {
    interval.summary <- .np_plot_bootstrap_interval_summary(
      boot.t = boot.out$t,
      t0 = boot.out$t0,
      alpha = alpha,
      band.type = band
    )
    boot.err <- matrix(NA_real_, nrow = length(t0), ncol = 3L)
    boot.err[, 1:2] <- interval.summary$err
    boot.all.err <- interval.summary$all.err
  }
  if (identical(center, "bias-corrected"))
    boot.err[, 3L] <- 2 * boot.out$t0 - colMeans(boot.out$t)
  list(
    boot.err = boot.err,
    boot.all.err = boot.all.err,
    bxp = list(),
    boot.t = boot.out$t
  )
}

.np_plot_conmode_add_interval_columns <- function(dat,
                                                  se,
                                                  errors,
                                                  alpha,
                                                  band,
                                                  center = "estimate",
                                                  bootstrap_raw = NULL) {
  payload <- .np_plot_interval_payload(
    estimate = dat$probability,
    se = se,
    plot.errors.method = errors,
    plot.errors.alpha = alpha,
    plot.errors.type = band,
    plot.errors.center = center,
    bootstrap_raw = bootstrap_raw
  )
  dat$stderr <- if (identical(errors, "bootstrap") &&
                    !is.null(bootstrap_raw$boot.t)) {
    .np_plot_bootstrap_col_sds(bootstrap_raw$boot.t)
  } else {
    as.numeric(se)
  }
  dat$center <- payload$center
  dat$lower <- payload$center - payload$err[, 1L]
  dat$upper <- payload$center + payload$err[, 2L]
  if (!is.null(payload$all.err)) {
    for (nm in names(payload$all.err)) {
      dat[[paste0(nm, ".lower")]] <- payload$center - payload$all.err[[nm]][, 1L]
      dat[[paste0(nm, ".upper")]] <- payload$center + payload$all.err[[nm]][, 2L]
    }
  }
  dat
}

.np_plot_conmode_grid_data <- function(object,
                                       gradients = FALSE,
                                       level = NULL,
                                       neval = 50L,
                                       xtrim = c(0, 1),
                                       xq = NULL,
                                       errors = "none",
                                       alpha = 0.05,
                                       band = "pointwise",
                                       center = "estimate",
                                       bootstrap = "inid",
                                       B = 399L,
                                       blocklen = NULL) {
  if (!isTRUE(gradients) && is.null(object$probabilities))
    stop("class probabilities are not available: fit with probabilities=TRUE or gradients=TRUE",
         call. = FALSE)
  if (isTRUE(gradients) && is.null(object$probability.gradients))
    stop("class-probability gradients/effects are not available: fit with gradients=TRUE",
         call. = FALSE)

  training <- .np_plot_conmode_training_data(object)
  xtrain <- training$xtrain
  ytrain <- training$ytrain
  resolved <- .np_plot_conmode_resolve_level(object, ytrain, gradients, level)
  level <- resolved$level
  lev <- resolved$levels
  neval <- .np_plot_conmode_resolve_neval(neval)
  xtrim <- .np_plot_conmode_resolve_xtrim(xtrim)

  p <- ncol(xtrain)
  xnames <- object$xnames
  if (is.null(xnames) || length(xnames) != p)
    xnames <- names(xtrain)
  xq <- .np_plot_conmode_resolve_xq(xq, p)

  rhs <- rep.int(1.0, nrow(xtrain))
  base <- .np_plot_conmode_base_row(xtrain, xq)
  proper.control <- list()
  if (!is.null(object$proper.info$tol))
    proper.control$tol <- object$proper.info$tol
  out <- vector("list", p)
  names(out) <- xnames
  ixcon <- object$bws$ixcon
  if (is.null(ixcon) || length(ixcon) != p)
    ixcon <- vapply(xtrain, is.numeric, logical(1L))

  for (j in seq_len(p)) {
    grid <- .np_plot_conmode_grid_values(xtrain[[j]], neval = neval,
                                         xtrim = xtrim)
    exdat <- base[rep(1L, length(grid)), , drop = FALSE]
    exdat[[j]] <- .np_plot_conmode_cast_like(grid, xtrain[[j]])

    if (isTRUE(gradients)) {
      if (isTRUE(ixcon[j])) {
        s <- integer(sum(ixcon))
        names(s) <- xnames[ixcon]
        s[xnames[j]] <- 1L
        vals <- npcdenshat(
          bws = object$bws,
          txdat = xtrain,
          tydat = ytrain,
          exdat = exdat,
          eydat = .np_plot_conmode_level_factor(ytrain, level, nrow(exdat)),
          y = rhs,
          output = "apply",
          s = s
        )
      } else {
        vals <- rep(0, nrow(exdat))
      }
      out[[j]] <- data.frame(
        variable = xnames[j],
        x = exdat[[j]],
        effect = as.vector(vals),
        level = level,
        gradients = TRUE,
        view = "fixed",
        stringsAsFactors = FALSE
      )
    } else {
      pmat <- matrix(NA_real_, nrow(exdat), length(lev),
                     dimnames = list(NULL, lev))
      for (k in seq_along(lev)) {
        pmat[, k] <- npcdenshat(
          bws = object$bws,
          txdat = xtrain,
          tydat = ytrain,
          exdat = exdat,
          eydat = .np_plot_conmode_level_factor(ytrain, lev[k], nrow(exdat)),
          y = rhs,
          output = "apply"
        )
      }
      proper.out <- .npConmodeProperProbabilities(
        pmat,
        levels = lev,
        proper = isTRUE(object$proper.requested),
        proper.control = proper.control
      )
      dat <- data.frame(
        variable = xnames[j],
        x = exdat[[j]],
        probability = proper.out$probabilities[, match(level, lev)],
        level = level,
        gradients = FALSE,
        view = "fixed",
        stringsAsFactors = FALSE
      )
      if (!identical(errors, "none")) {
        boot.raw <- NULL
        if (identical(errors, "bootstrap")) {
          boot.raw <- .np_plot_conmode_bootstrap_raw(
            object = object,
            xtrain = xtrain,
            ytrain = ytrain,
            exdat = exdat,
            level = level,
            levels = lev,
            t0 = dat$probability,
            alpha = alpha,
            band = band,
            center = center,
            method = bootstrap,
            B = B,
            blocklen = blocklen
          )
          se <- rep(NA_real_, nrow(exdat))
        } else {
          se <- .np_plot_conmode_probability_se(
            object = object,
            xtrain = xtrain,
            ytrain = ytrain,
            exdat = exdat,
            level = level
          )
          se[proper.out$repaired.rows] <- NA_real_
        }
        dat <- .np_plot_conmode_add_interval_columns(
          dat,
          se = se,
          errors = errors,
          alpha = alpha,
          band = band,
          center = center,
          bootstrap_raw = boot.raw
        )
      }
      out[[j]] <- dat
    }
  }
  out
}

.np_plot_conmode_resolve_surface_vars <- function(xtrain, object, plot.vars = NULL) {
  p <- ncol(xtrain)
  xnames <- object$xnames
  if (is.null(xnames) || length(xnames) != p)
    xnames <- names(xtrain)
  ixcon <- object$bws$ixcon
  if (is.null(ixcon) || length(ixcon) != p)
    ixcon <- vapply(xtrain, is.numeric, logical(1L))
  con.names <- xnames[ixcon]
  if (is.null(plot.vars)) {
    if (length(con.names) != 2L)
      stop("surface plot.conmode payloads require exactly two continuous conditioning variables, or a length-two 'plot.vars' selection",
           call. = FALSE)
    return(con.names)
  }
  if (!is.character(plot.vars) || length(plot.vars) != 2L ||
      any(is.na(plot.vars)) || any(!nzchar(plot.vars)))
    stop("'plot.vars' must be a character vector naming exactly two conditioning variables",
         call. = FALSE)
  if (any(!(plot.vars %in% xnames)))
    stop("'plot.vars' must name conditioning variables in the fitted conmode object",
         call. = FALSE)
  if (any(!(plot.vars %in% con.names)))
    stop("surface plot.conmode payloads currently require continuous variables in 'plot.vars'",
         call. = FALSE)
  plot.vars
}

.np_plot_conmode_surface_data <- function(object,
                                          level = NULL,
                                          neval = 50L,
                                          xtrim = c(0, 1),
                                          xq = NULL,
                                          plot.vars = NULL,
                                          errors = "none",
                                          alpha = 0.05,
                                          band = "pointwise") {
  if (is.null(object$probabilities))
    stop("class probabilities are not available: fit with probabilities=TRUE or gradients=TRUE",
         call. = FALSE)
  training <- .np_plot_conmode_training_data(object)
  xtrain <- training$xtrain
  ytrain <- training$ytrain
  resolved <- .np_plot_conmode_resolve_level(object, ytrain, FALSE, level)
  level <- resolved$level
  lev <- resolved$levels
  neval <- .np_plot_conmode_resolve_neval(neval)
  xtrim <- .np_plot_conmode_resolve_xtrim(xtrim)

  p <- ncol(xtrain)
  xnames <- object$xnames
  if (is.null(xnames) || length(xnames) != p)
    xnames <- names(xtrain)
  xq <- .np_plot_conmode_resolve_xq(xq, p)
  plot.vars <- .np_plot_conmode_resolve_surface_vars(
    xtrain = xtrain,
    object = object,
    plot.vars = plot.vars
  )
  ivars <- match(plot.vars, xnames)

  grid1 <- .np_plot_conmode_grid_values(xtrain[[ivars[1L]]],
                                        neval = neval,
                                        xtrim = xtrim)
  grid2 <- .np_plot_conmode_grid_values(xtrain[[ivars[2L]]],
                                        neval = neval,
                                        xtrim = xtrim)
  base <- .np_plot_conmode_base_row(xtrain, xq)
  grid.index <- expand.grid(
    i1 = seq_along(grid1),
    i2 = seq_along(grid2),
    KEEP.OUT.ATTRS = FALSE
  )
  exdat <- base[rep(1L, nrow(grid.index)), , drop = FALSE]
  exdat[[ivars[1L]]] <- .np_plot_conmode_cast_like(grid1[grid.index$i1],
                                                   xtrain[[ivars[1L]]])
  exdat[[ivars[2L]]] <- .np_plot_conmode_cast_like(grid2[grid.index$i2],
                                                   xtrain[[ivars[2L]]])

  rhs <- rep.int(1.0, nrow(xtrain))
  pmat <- matrix(NA_real_, nrow(exdat), length(lev),
                 dimnames = list(NULL, lev))
  for (k in seq_along(lev)) {
    pmat[, k] <- npcdenshat(
      bws = object$bws,
      txdat = xtrain,
      tydat = ytrain,
      exdat = exdat,
      eydat = .np_plot_conmode_level_factor(ytrain, lev[k], nrow(exdat)),
      y = rhs,
      output = "apply"
    )
  }
  proper.control <- list()
  if (!is.null(object$proper.info$tol))
    proper.control$tol <- object$proper.info$tol
  proper.out <- .npConmodeProperProbabilities(
    pmat,
    levels = lev,
    proper = isTRUE(object$proper.requested),
    proper.control = proper.control
  )

  surface <- data.frame(
    variable1 = plot.vars[1L],
    variable2 = plot.vars[2L],
    grid_index1 = grid.index$i1,
    grid_index2 = grid.index$i2,
    x1 = exdat[[ivars[1L]]],
    x2 = exdat[[ivars[2L]]],
    probability = proper.out$probabilities[, match(level, lev)],
    level = level,
    view = "fixed",
    perspective = TRUE,
    stringsAsFactors = FALSE
  )
  if (!identical(errors, "none")) {
    se <- .np_plot_conmode_probability_se(
      object = object,
      xtrain = xtrain,
      ytrain = ytrain,
      exdat = exdat,
      level = level
    )
    se[proper.out$repaired.rows] <- NA_real_
    surface <- .np_plot_conmode_add_interval_columns(
      surface,
      se = se,
      errors = errors,
      alpha = alpha,
      band = band
    )
  }
  held.vars <- setdiff(xnames, plot.vars)
  held <- if (length(held.vars)) base[, held.vars, drop = FALSE] else data.frame()
  structure(
    list(
      surface = surface,
      grid = list(
        variables = plot.vars,
        values = setNames(list(grid1, grid2), plot.vars),
        dims = c(length(grid1), length(grid2))
      ),
      held = held,
      level = level,
      levels = lev,
      proper = proper.out$proper.info
    ),
    class = "np_conmode_surface_data"
  )
}

.np_plot_conmode_surface_matrix <- function(payload) {
  dims <- payload$grid$dims
  matrix(payload$surface$probability, nrow = dims[1L], ncol = dims[2L])
}

.np_plot_conmode_validate_zlim <- function(zlim) {
  if (is.null(zlim))
    return(NULL)
  if (!is.numeric(zlim) || length(zlim) != 2L || any(is.na(zlim)) ||
      zlim[1L] >= zlim[2L])
    stop("zlim must be a numeric length-two vector with zlim[1] < zlim[2]",
         call. = FALSE)
  zlim
}

.np_plot_conmode_surface_panel <- function(payload, dots, renderer = "base") {
  vars <- payload$grid$variables
  x1 <- payload$grid$values[[1L]]
  x2 <- payload$grid$values[[2L]]
  z <- .np_plot_conmode_surface_matrix(payload)
  zlim <- .np_plot_conmode_validate_zlim(dots$zlim)
  if (is.null(zlim))
    zlim <- c(0, 1)

  theta <- .np_plot_scalar_default(dots$theta, 0)
  phi <- .np_plot_scalar_default(dots$phi, 20)
  main <- .np_plot_scalar_default(
    dots$main,
    paste0("Pr(Y=", payload$level, "|X=x)")
  )
  xlab <- .np_plot_scalar_default(dots$xlab, vars[1L])
  ylab <- .np_plot_scalar_default(dots$ylab, vars[2L])
  zlab <- .np_plot_scalar_default(dots$zlab, "Probability")
  border <- .np_plot_scalar_default(dots$border, "black")

  if (identical(renderer, "rgl")) {
    rgl.view <- .np_plot_rgl_view_angles(theta = theta, phi = phi)
    return(.np_plot_render_surface_rgl(
      x = x1,
      y = x2,
      z = z,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      main = main,
      theta = rgl.view$theta,
      phi = rgl.view$phi,
      col = dots$col,
      border = border,
      zlim = zlim,
      par3d.args = .np_plot_collect_rgl_args(dots, "rgl.par3d", "rgl.par3d."),
      view3d.args = .np_plot_collect_rgl_args(dots, "rgl.view3d", "rgl.view3d."),
      persp3d.args = .np_plot_collect_rgl_args(dots, "rgl.persp3d", "rgl.persp3d."),
      grid3d.args = .np_plot_collect_rgl_args(dots, "rgl.grid3d", "rgl.grid3d."),
      widget.args = .np_plot_collect_rgl_args(dots, "rgl.widget", "rgl.widget.")
    ))
  }

  persp.args <- list(
    x = x1,
    y = x2,
    z = z,
    theta = theta,
    phi = phi,
    ticktype = "detailed",
    xlab = xlab,
    ylab = ylab,
    zlab = zlab,
    zlim = zlim,
    main = main,
    col = grDevices::adjustcolor(
      .np_plot_persp_surface_colors(z = z, col = dots$col),
      alpha.f = 0.5
    ),
    border = border
  )
  persp.args <- .np_plot_merge_user_args(
    persp.args,
    .np_plot_user_args(dots, "persp")
  )
  do.call(graphics::persp, persp.args)
  invisible(NULL)
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

  fixed.grid <- identical(dat$view[1L], "fixed")
  if (is.factor(x) && !is.ordered(x) && !isTRUE(fixed.grid)) {
    args <- .np_plot_merge_user_args(
      list(formula = y ~ x, xlab = xlab, ylab = ylab, main = main),
      plot.user.args
    )
    do.call(graphics::boxplot, args)
    return(invisible(NULL))
  }

  xplot <- if (is.factor(x)) as.numeric(x) else x
  ord <- order(xplot)
  interval.cols <- !isTRUE(gradients) && all(c("lower", "upper") %in% names(dat))
  ylim.interval <- NULL
  if (interval.cols) {
    ylim.interval <- range(c(y, dat$lower, dat$upper), finite = TRUE)
    if (!all(is.finite(ylim.interval)))
      ylim.interval <- NULL
  }
  args <- .np_plot_merge_user_args(
    list(x = xplot[ord], y = y[ord],
         type = if (is.factor(x) && !is.ordered(x)) "p" else "l",
         xlab = xlab, ylab = ylab, main = main),
    line.user.args
  )
  if (is.null(args$ylim) && !is.null(ylim.interval))
    args$ylim <- ylim.interval
  args <- .np_plot_merge_user_args(
    args,
    plot.user.args
  )
  if (is.factor(x)) {
    args$xaxt <- "n"
    do.call(graphics::plot, args)
    graphics::axis(1L, at = seq_along(levels(x)), labels = levels(x))
  } else {
    do.call(graphics::plot, args)
  }
  if (isTRUE(data_rug) && is.numeric(xplot))
    graphics::rug(xplot)
  if (interval.cols) {
    lower <- dat$lower[ord]
    upper <- dat$upper[ord]
    xord <- xplot[ord]
    good <- is.finite(xord) & is.finite(lower) & is.finite(upper)
    if (any(good)) {
      if (is.factor(x) && !is.ordered(x)) {
        graphics::segments(xord[good], lower[good], xord[good], upper[good],
                           col = "gray45", lty = 2)
      } else {
        graphics::lines(xord[good], lower[good], col = "gray45", lty = 2)
        graphics::lines(xord[good], upper[good], col = "gray45", lty = 2)
      }
    }
  }
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
                     "render_control",
                     "boot_control", "plot.errors.method", "plot.errors.type",
                     "plot.errors.alpha", "plot.errors.boot.method",
                     "plot.errors.boot.num", "plot.errors.boot.nonfixed",
                     "plot.errors.boot.wild", "plot.errors.boot.blocklen",
                     "plot.errors.center", "plot.errors.style",
                     "plot.errors.bar", "plot.errors.bar.num")
  surface.args <- c("renderer", "perspective", "persp", "plot.vars")
  grid.args <- c("neval", "grid_control", "view", "xtrim", "xq")
  grid.supplied <- intersect(dot.names[nzchar(dot.names)], grid.args)
  rgl.prefixed.args <- dot.names[startsWith(dot.names, "rgl.")]
  allowed <- unique(c("gradients", "level", "output", "data_rug",
                      "layout", "legend", grid.args, surface.args,
                      interval.args, rgl.prefixed.args,
                      .np_plot_graphics_arg_names()))
  .np_plot_stop_unused_args(setdiff(dot.names[nzchar(dot.names)], allowed),
                            allowed)
  perspective.raw <- dots$perspective
  persp.raw <- dots$persp
  dots$perspective <- NULL
  dots$persp <- NULL
  dots <- .np_plot_normalize_public_dots(dots, context = "plot.conmode")
  if (!is.null(perspective.raw))
    dots$perspective <- perspective.raw
  if (!is.null(persp.raw))
    dots$persp <- persp.raw
  if (!is.null(dots$plot.behavior)) {
    output <- dots$plot.behavior
    dots$plot.behavior <- NULL
  }
  if (!is.null(dots$plot.rug)) {
    data_rug <- dots$plot.rug
    dots$plot.rug <- NULL
  }
  if (!is.null(dots$plot.par.mfrow)) {
    layout <- dots$plot.par.mfrow
    dots$plot.par.mfrow <- NULL
  }
  errors <- if (is.null(dots$plot.errors.method)) "none" else
    .np_plot_scalar_match(dots$plot.errors.method,
                          c("none", "bootstrap", "asymptotic"),
                          "plot.errors.method")
  band <- if (is.null(dots$plot.errors.type)) "pointwise" else
    .np_plot_scalar_match(dots$plot.errors.type,
                          c("pmzsd", "pointwise", "bonferroni",
                            "simultaneous", "all"),
                          "plot.errors.type")
  alpha <- if (is.null(dots$plot.errors.alpha)) 0.05 else
    dots$plot.errors.alpha
  center <- if (is.null(dots$plot.errors.center)) "estimate" else
    .np_plot_scalar_match(dots$plot.errors.center,
                          c("estimate", "bias-corrected"),
                          "plot.errors.center")
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha <= 0 || alpha >= 0.5)
    stop("plot.errors.alpha must lie in (0, 0.5)", call. = FALSE)
  bootstrap <- if (is.null(dots$plot.errors.boot.method)) "inid" else
    .np_plot_scalar_match(dots$plot.errors.boot.method,
                          c("wild", "inid", "fixed", "geom"),
                          "plot.errors.boot.method")
  B <- if (is.null(dots$plot.errors.boot.num)) 399L else
    dots$plot.errors.boot.num
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1)
    stop("plot.errors.boot.num must be a positive numeric scalar",
         call. = FALSE)
  B <- as.integer(B)
  boot.blocklen <- dots$plot.errors.boot.blocklen
  if (!identical(center, "estimate") && !identical(errors, "bootstrap"))
    stop("bias-corrected interval centering requires errors=\"bootstrap\"",
         call. = FALSE)
  unsupported.interval <- intersect(
    names(dots),
    c("plot.errors.style", "plot.errors.bar", "plot.errors.bar.num")
  )
  if (length(unsupported.interval))
    stop("plot.conmode probability intervals currently use simple asymptotic line/error-bar rendering; render_control/style arguments are not yet implemented",
         call. = FALSE)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  if (!identical(errors, "none") && isTRUE(gradients))
    stop("plot.conmode intervals are available for class probabilities, not probability gradients/effects",
         call. = FALSE)
  data_rug <- npValidateScalarLogical(data_rug, "data_rug")
  renderer <- "base"
  if (!is.null(dots$renderer))
    renderer <- .np_plot_match_renderer(dots$renderer)
  perspective <- FALSE
  perspective.supplied <- !is.null(dots$perspective)
  persp.supplied <- !is.null(dots$persp)
  if (perspective.supplied)
    perspective <- npValidateScalarLogical(dots$perspective, "perspective")
  if (persp.supplied) {
    persp.value <- npValidateScalarLogical(dots$persp, "persp")
    if (perspective.supplied && !identical(perspective, persp.value))
      stop("conflicting plot.conmode arguments: 'perspective' and 'persp'",
           call. = FALSE)
    perspective <- persp.value
  }
  if (!is.null(dots$plot.vars) && !isTRUE(perspective))
    stop("'plot.vars' is only supported for fixed surface payloads with perspective=TRUE",
         call. = FALSE)
  if (!identical(renderer, "base") && !isTRUE(perspective))
    stop("renderer='rgl' is supported only for fixed surface payloads with perspective=TRUE",
         call. = FALSE)
  if (isTRUE(perspective) && isTRUE(gradients))
    stop("surface plot.conmode payloads are currently available for probabilities, not gradients",
         call. = FALSE)
  view <- dots$view
  if (is.null(view))
    view <- if (isTRUE(perspective) || length(grid.supplied)) "fixed" else "sample"
  view <- .np_plot_scalar_match(view, c("sample", "fixed"), "view")
  if (isTRUE(perspective) && !identical(view, "fixed"))
    stop("surface plot.conmode payloads require view='fixed'",
         call. = FALSE)
  neval <- dots$neval
  if (is.null(neval))
    neval <- 50L
  xtrim <- dots$xtrim
  xq <- dots$xq
  if (!is.null(dots$grid_control)) {
    if (!inherits(dots$grid_control, "np_grid_control"))
      stop("grid_control must be created by np_grid_control()", call. = FALSE)
    if (!is.null(dots$grid_control$slices))
      stop("grid_control$slices is not yet supported for plot.conmode",
           call. = FALSE)
    if (!is.null(dots$grid_control$xtrim))
      xtrim <- dots$grid_control$xtrim
    if (!is.null(dots$grid_control$xq))
      xq <- dots$grid_control$xq
  }
  if (is.null(xtrim))
    xtrim <- c(0, 1)
  plot.par.mfrow <- .np_plot_match_layout(layout)
  output <- match.arg(output)
  if (identical(output, "both"))
    output <- "plot-data"
  output <- .np_plot_scalar_match(output, c("plot", "data", "plot-data"), "output")
  if (!identical(errors, "none") && isTRUE(perspective) &&
      !identical(output, "data"))
    stop("surface probability intervals are available as output=\"data\" only; surface band rendering is not yet implemented for plot.conmode",
         call. = FALSE)
  if (identical(errors, "bootstrap") && isTRUE(perspective))
    stop("surface bootstrap intervals are not yet implemented for plot.conmode; use one-dimensional views or errors=\"asymptotic\" with output=\"data\"",
         call. = FALSE)
  plot.user.args <- .np_plot_user_args(dots, "plot")
  line.user.args <- .np_plot_user_args(dots, "lines")
  if (!is.null(dots$type))
    line.user.args$type <- dots$type

  plot.data <- if (isTRUE(perspective)) {
    .np_plot_conmode_surface_data(
      object,
      level = level,
      neval = neval,
      xtrim = xtrim,
      xq = xq,
      plot.vars = dots$plot.vars,
      errors = errors,
      alpha = alpha,
      band = band
    )
  } else if (identical(view, "fixed")) {
    .np_plot_conmode_grid_data(
      object,
      gradients = gradients,
      level = level,
      neval = neval,
      xtrim = xtrim,
      xq = xq,
      errors = errors,
      alpha = alpha,
      band = band,
      center = center,
      bootstrap = bootstrap,
      B = B,
      blocklen = boot.blocklen
    )
  } else {
    .np_plot_conmode_data(
      object,
      gradients = gradients,
      level = level,
      errors = errors,
      alpha = alpha,
      band = band,
      center = center,
      bootstrap = bootstrap,
      B = B,
      blocklen = boot.blocklen
    )
  }

  if (isTRUE(perspective)) {
    if (identical(output, "data"))
      return(plot.data)
    render.out <- .np_plot_conmode_surface_panel(plot.data, dots,
                                                 renderer = renderer)
    if (identical(renderer, "rgl"))
      return(.np_plot_rgl_finalize(
        rgl.out = render.out,
        plot.behavior = output,
        plot.data = plot.data
      ))
    if (identical(output, "plot-data"))
      return(plot.data)
    return(invisible(object))
  }
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
