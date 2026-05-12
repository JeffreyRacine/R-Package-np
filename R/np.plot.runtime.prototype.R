.np_plot_proto_check_npcdens_fixed_surface <- function(bws,
                                                       xdat,
                                                       ydat,
                                                       neval,
                                                       xtrim,
                                                       ytrim,
                                                       cdf = FALSE) {
  expected.class <- if (isTRUE(cdf)) "condbandwidth" else "conbandwidth"
  expected.label <- if (isTRUE(cdf)) "conditional distribution" else "conditional density"
  if (!inherits(bws, expected.class))
    stop(sprintf("prototype route requires a %s bandwidth object", expected.label), call. = FALSE)
  regtype.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  if (!is.element(regtype.engine, c("lc", "lp")))
    stop("prototype route currently supports regtype='lc', 'll', or 'lp' only", call. = FALSE)
  if (!is.element(as.character(bws$type), c("fixed", "generalized_nn", "adaptive_nn")))
    stop("prototype route currently supports fixed, generalized-NN, or adaptive-NN bandwidths only", call. = FALSE)
  if (bws$xndim != 1L || bws$yndim != 1L)
    stop("prototype route currently supports one x variable and one y variable", call. = FALSE)
  if (bws$xnuno + bws$ynuno != 0L)
    stop("prototype route currently supports continuous/ordered surface variables only", call. = FALSE)
  if (bws$xncon + bws$xnord + bws$yncon + bws$ynord != 2L)
    stop(sprintf("prototype route requires a two-dimensional %s surface", expected.label), call. = FALSE)
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("prototype route requires scalar neval >= 2", call. = FALSE)
  invisible(TRUE)
}

.np_plot_proto_check_npcdens_lc_fixed_none <- .np_plot_proto_check_npcdens_fixed_surface

.np_plot_proto_clean_conditional_data <- function(xdat, ydat) {
  ## Contract: align explicit training data for the first prototype route. This
  ## stage intentionally does not recover data from formula/call objects; that
  ## wider state-resolution contract belongs to a later slice.
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE
  if (!any(keep.rows))
    stop("Data has no rows without NAs")
  xdat <- xdat[keep.rows, , drop = FALSE]
  ydat <- ydat[keep.rows, , drop = FALSE]

  list(xdat = xdat, ydat = ydat)
}

.np_plot_proto_conditional_surface_grid <- function(bws,
                                                    xdat,
                                                    ydat,
                                                    neval,
                                                    xtrim,
                                                    ytrim) {
  ## Contract: build the two-dimensional conditional surface grid. This helper
  ## owns target construction only; it must not evaluate the estimator or create
  ## plot objects.
  xtrim <- double(bws$xndim) + xtrim
  ytrim <- double(bws$yndim) + ytrim

  if (is.ordered(xdat[, 1L])) {
    x1.eval <- bws$xdati$all.ulev[[1L]]
    x1.neval <- length(x1.eval)
  } else {
    x1.neval <- as.integer(neval)
    qi <- trim.quantiles(xdat[, 1L], xtrim[1L])
    x1.eval <- seq(qi[1L], qi[2L], length.out = x1.neval)
  }

  tx2 <- ydat[, 1L]
  txdati <- bws$ydati
  txtrim <- ytrim

  if (txdati$iord[1L]) {
    x2.eval <- txdati$all.ulev[[1L]]
    x2.neval <- length(x2.eval)
  } else {
    x2.neval <- as.integer(neval)
    qi <- trim.quantiles(tx2, txtrim[1L])
    x2.eval <- seq(qi[1L], qi[2L], length.out = x2.neval)
  }

  x.eval <- expand.grid(x1.eval, x2.eval)
  if (bws$xdati$iord[1L])
    x1.eval <- bws$xdati$all.dlev[[1L]][as.integer(x1.eval)]
  if (txdati$iord[1L])
    x2.eval <- txdati$all.dlev[[1L]][as.integer(x2.eval)]

  list(
    x.eval = x.eval,
    x1.eval = x1.eval,
    x2.eval = x2.eval,
    x1.neval = x1.neval,
    x2.neval = x2.neval
  )
}

.np_plot_proto_npcdens_fixed_data <- function(bws,
                                              xdat,
                                              ydat,
                                              neval = 50,
                                              xtrim = 0.0,
                                              ytrim = 0.0,
                                              plot.errors.method = c("none", "asymptotic", "bootstrap"),
                                              plot.errors.boot.method = c("inid", "fixed", "geom"),
                                              plot.errors.boot.nonfixed = c("exact", "frozen"),
                                              plot.errors.boot.blocklen = NULL,
                                              plot.errors.boot.num = 1999,
                                              plot.errors.center = c("estimate", "bias-corrected"),
                                              plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                   "simultaneous", "all"),
                                              plot.errors.alpha = 0.05,
                                              proper = FALSE,
                                              proper.method = c("project"),
                                              proper.control = list(),
                                              return.stages = FALSE,
                                              cdf = FALSE) {
  ## Contract: private conditional density/distribution data-only prototype.
  ## This owns explicit data cleanup, target construction, evaluator invocation,
  ## optional asymptotic/bootstrap interval construction, and old-compatible
  ## plot-data assembly.
  ## It must not draw graphics, bootstrap, change RNG state, or recover formula
  ## data until those stages receive their own slice.
  if (missing(xdat) || missing(ydat))
    stop("prototype route requires explicit xdat and ydat", call. = FALSE)
  plot.errors.method <- match.arg(plot.errors.method)
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)
  if (missing(proper.method))
    proper.method <- if (isTRUE(cdf)) "isotonic" else "project"
  proper.args <- if (isTRUE(cdf)) {
    .np_condist_validate_proper_args(
      proper = proper,
      proper.method = proper.method,
      proper.control = proper.control
    )
  } else {
    .np_condens_validate_proper_args(
      proper = proper,
      proper.method = proper.method,
      proper.control = proper.control
    )
  }
  dat <- .np_plot_proto_clean_conditional_data(xdat = xdat, ydat = ydat)
  xdat <- dat$xdat
  ydat <- dat$ydat
.np_plot_proto_check_npcdens_fixed_surface(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    cdf = cdf
  )

  grid <- .np_plot_proto_conditional_surface_grid(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim
  )
  fit <- .np_plot_conditional_eval(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = grid$x.eval[, 1L, drop = FALSE],
    eydat = grid$x.eval[, 2L, drop = FALSE],
    cdf = isTRUE(cdf),
    gradients = FALSE,
    proper = isTRUE(proper.args$proper.requested),
    proper.method = proper.args$proper.method,
    proper.control = proper.args$proper.control
  )

  tcomp <- if (isTRUE(cdf)) fit$condist else fit$condens
  terr <- matrix(fit$conderr, nrow = length(tcomp), ncol = 3L)
  terr[, 3L] <- NA_real_
  interval <- NULL
  bootstrap <- NULL
  if (identical(plot.errors.method, "asymptotic")) {
    interval <- .np_plot_asymptotic_error_from_se(
      se = fit$conderr,
      alpha = plot.errors.alpha,
      band.type = plot.errors.type,
      m = nrow(grid$x.eval[, 1L])
    )
    terr[, 1:2] <- interval$err
  } else if (identical(plot.errors.method, "bootstrap")) {
    bootstrap <- compute.bootstrap.errors(
      xdat = xdat,
      ydat = ydat,
      exdat = grid$x.eval[, 1L],
      eydat = grid$x.eval[, 2L],
      cdf = isTRUE(cdf),
      quantreg = FALSE,
      tau = 0.5,
      gradients = FALSE,
      gradient.index = 0L,
      slice.index = 0L,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      progress.target = NULL,
      proper = isTRUE(proper.args$proper.requested),
      proper.method = proper.args$proper.method,
      proper.control = proper.args$proper.control,
      bws = bws
    )
    terr <- bootstrap$boot.err
    interval <- list(
      err = bootstrap$boot.err[, 1:2, drop = FALSE],
      all.err = bootstrap$boot.all.err
    )
  }

  cd.args <- list(
    bws = bws,
    xeval = grid$x.eval[, 1L],
    yeval = grid$x.eval[, 2L],
    ntrain = nrow(xdat),
    conderr = terr[, 1:2, drop = FALSE],
    proper.requested = fit$proper.requested,
    proper.applied = fit$proper.applied,
    proper.method = fit$proper.method,
    proper.info = fit$proper.info
  )
  if (isTRUE(cdf)) {
    cd.args$condist <- fit$condist
    cd.args$condist.raw <- fit$condist.raw
    cd1 <- do.call(condistribution, cd.args)
  } else {
    cd.args$condens <- fit$condens
    cd.args$condens.raw <- fit$condens.raw
    cd1 <- do.call(condensity, cd.args)
  }
  cd1$bias <- NA

  plot.data <- list(cd1 = cd1)
  if (!isTRUE(return.stages))
    return(plot.data)

  list(
    state = list(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      ntrain = nrow(xdat),
      family = if (isTRUE(cdf)) "npcdist" else "npcdens",
      cdf = isTRUE(cdf),
      gradients = FALSE,
      proper = proper.args
    ),
    target_grid = grid,
    evaluator = fit,
    intervals = if (is.null(interval)) NULL else list(
      method = plot.errors.method,
      type = plot.errors.type,
      alpha = plot.errors.alpha,
      err = interval$err,
      all.err = interval$all.err
    ),
    bootstrap = if (is.null(bootstrap)) NULL else list(
      method = plot.errors.boot.method,
      nonfixed = plot.errors.boot.nonfixed,
      blocklen = plot.errors.boot.blocklen,
      B = plot.errors.boot.num,
      center = plot.errors.center,
      boot.err = bootstrap$boot.err,
      boot.all.err = bootstrap$boot.all.err,
      bxp = bootstrap$bxp
    ),
    proper_projection = if (!isTRUE(proper.args$proper.requested)) NULL else list(
      requested = fit$proper.requested,
      applied = fit$proper.applied,
      method = fit$proper.method,
      info = fit$proper.info
    ),
    plot_data = plot.data
  )
}

.np_plot_proto_npcdens_lc_fixed_data <- .np_plot_proto_npcdens_fixed_data

.np_plot_proto_npcdens_lc_fixed_none_data <- function(bws,
                                                      xdat,
                                                      ydat,
                                                      neval = 50,
                                                      xtrim = 0.0,
                                                      ytrim = 0.0,
                                                      proper = FALSE,
                                                      proper.method = c("project"),
                                                      proper.control = list(),
                                                      return.stages = FALSE) {
  .np_plot_proto_npcdens_lc_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    plot.errors.method = "none",
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control,
    return.stages = return.stages
  )
}

.np_plot_proto_npcdens_fixed_none_data <- function(bws,
                                                   xdat,
                                                   ydat,
                                                   neval = 50,
                                                   xtrim = 0.0,
                                                   ytrim = 0.0,
                                                   proper = FALSE,
                                                   proper.method = c("project"),
                                                   proper.control = list(),
                                                   return.stages = FALSE) {
  .np_plot_proto_npcdens_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    plot.errors.method = "none",
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control,
    return.stages = return.stages
  )
}

.np_plot_proto_npcdist_fixed_none_data <- function(bws,
                                                   xdat,
                                                   ydat,
                                                   neval = 50,
                                                   xtrim = 0.0,
                                                   ytrim = 0.0,
                                                   proper = FALSE,
                                                   proper.method = c("isotonic"),
                                                   proper.control = list(),
                                                   return.stages = FALSE) {
  .np_plot_proto_npcdens_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    plot.errors.method = "none",
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control,
    return.stages = return.stages,
    cdf = TRUE
  )
}

.np_plot_proto_npcdist_fixed_bootstrap_inid_data <- function(bws,
                                                            xdat,
                                                            ydat,
                                                            neval = 50,
                                                            xtrim = 0.0,
                                                            ytrim = 0.0,
                                                            plot.errors.boot.num = 1999,
                                                            plot.errors.center = c("estimate", "bias-corrected"),
                                                            plot.errors.type = c("pmzsd", "pointwise",
                                                                                 "bonferroni", "simultaneous",
                                                                                 "all"),
                                                            plot.errors.alpha = 0.05,
                                                            proper = FALSE,
                                                            proper.method = c("isotonic"),
                                                            proper.control = list(),
                                                            return.stages = FALSE) {
  .np_plot_proto_npcdens_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control,
    return.stages = return.stages,
    cdf = TRUE
  )
}

.np_plot_proto_npcdens_lc_fixed_asymptotic_data <- function(bws,
                                                           xdat,
                                                           ydat,
                                                           neval = 50,
                                                           xtrim = 0.0,
                                                           ytrim = 0.0,
                                                           plot.errors.type = c("pmzsd", "pointwise",
                                                                                "bonferroni", "simultaneous",
                                                                                "all"),
                                                           plot.errors.alpha = 0.05,
                                                           proper = FALSE,
                                                           proper.method = c("project"),
                                                           proper.control = list(),
                                                           return.stages = FALSE) {
  .np_plot_proto_npcdens_lc_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    plot.errors.method = "asymptotic",
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control,
    return.stages = return.stages
  )
}

.np_plot_proto_npcdens_fixed_bootstrap_inid_data <- function(bws,
                                                            xdat,
                                                            ydat,
                                                            neval = 50,
                                                            xtrim = 0.0,
                                                            ytrim = 0.0,
                                                            plot.errors.boot.num = 1999,
                                                            plot.errors.center = c("estimate", "bias-corrected"),
                                                            plot.errors.type = c("pmzsd", "pointwise",
                                                                                 "bonferroni", "simultaneous",
                                                                                 "all"),
                                                            plot.errors.alpha = 0.05,
                                                            proper = FALSE,
                                                            proper.method = c("project"),
                                                            proper.control = list(),
                                                            return.stages = FALSE) {
  .np_plot_proto_npcdens_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control,
    return.stages = return.stages
  )
}

.np_plot_proto_npcdens_lc_fixed_bootstrap_inid_data <- function(bws,
                                                               xdat,
                                                               ydat,
                                                               neval = 50,
                                                               xtrim = 0.0,
                                                               ytrim = 0.0,
                                                               plot.errors.boot.num = 1999,
                                                               plot.errors.center = c("estimate", "bias-corrected"),
                                                               plot.errors.type = c("pmzsd", "pointwise",
                                                                                    "bonferroni", "simultaneous",
                                                                                    "all"),
                                                               plot.errors.alpha = 0.05,
                                                               proper = FALSE,
                                                               proper.method = c("project"),
                                                               proper.control = list(),
                                                               return.stages = FALSE) {
  .np_plot_proto_npcdens_lc_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control,
    return.stages = return.stages
  )
}

.np_plot_proto_npcdens_lc_fixed_bootstrap_block_data <- function(bws,
                                                                xdat,
                                                                ydat,
                                                                neval = 50,
                                                                xtrim = 0.0,
                                                                ytrim = 0.0,
                                                                plot.errors.boot.method = c("fixed", "geom"),
                                                                plot.errors.boot.blocklen,
                                                                plot.errors.boot.num = 1999,
                                                                plot.errors.center = c("estimate", "bias-corrected"),
                                                                plot.errors.type = c("pmzsd", "pointwise",
                                                                                     "bonferroni", "simultaneous",
                                                                                     "all"),
                                                                plot.errors.alpha = 0.05,
                                                                proper = FALSE,
                                                                proper.method = c("project"),
                                                                proper.control = list(),
                                                                return.stages = FALSE) {
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  .np_plot_proto_npcdens_lc_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control,
    return.stages = return.stages
  )
}

.np_plot_proto_npcdens_surface_base_render <- function(plot.data,
                                                       perspective = TRUE,
                                                       main = NULL,
                                                       xlab = NULL,
                                                       ylab = NULL,
                                                       zlab = NULL,
                                                       col = "lightblue",
                                                       theta = 0.0,
                                                       phi = 20.0,
                                                       ...) {
  ## Contract: renderer smoke for staged npcdens surface data. It consumes a
  ## plot-data object and must not re-enter public estimators, bandwidth
  ## constructors, bootstrap helpers, or target builders.
  if (!is.list(plot.data) || is.null(plot.data$cd1) ||
      !(inherits(plot.data$cd1, "condensity") || inherits(plot.data$cd1, "condistribution")))
    stop("renderer prototype requires plot-data with a conditional density/distribution 'cd1' element", call. = FALSE)
  cd <- plot.data$cd1
  value.name <- if (inherits(cd, "condistribution")) "condist" else "condens"
  value.label <- if (identical(value.name, "condist")) "Conditional distribution" else "Conditional density"
  x <- unique(as.numeric(cd$xeval))
  y <- unique(as.numeric(cd$yeval))
  z <- matrix(as.numeric(cd[[value.name]]), nrow = length(x), ncol = length(y), byrow = FALSE)
  if (length(x) * length(y) != length(cd[[value.name]]) ||
      any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(z))) {
    stop("renderer prototype requires a finite rectangular conditional surface", call. = FALSE)
  }

  xlab <- if (is.null(xlab)) cd$xnames[1L] else xlab
  ylab <- if (is.null(ylab)) cd$ynames[1L] else ylab
  zlab <- if (is.null(zlab)) value.label else zlab
  main <- if (is.null(main)) value.label else main

  if (isTRUE(perspective)) {
    graphics::persp(
      x = x,
      y = y,
      z = z,
      theta = theta,
      phi = phi,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      main = main,
      col = col,
      ...
    )
  } else {
    graphics::image(
      x = x,
      y = y,
      z = z,
      xlab = xlab,
      ylab = ylab,
      main = main,
      col = grDevices::hcl.colors(64L, "YlOrRd", rev = TRUE),
      ...
    )
  }
  invisible(plot.data)
}

.np_plot_proto_check_npreg_fixed_surface <- function(bws,
                                                     xdat,
                                                     ydat,
                                                     neval,
                                                     xtrim) {
  if (!inherits(bws, "rbandwidth"))
    stop("prototype route requires a regression bandwidth object", call. = FALSE)
  regtype.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  if (!is.element(regtype.engine, c("lc", "lp")))
    stop("prototype route currently supports regtype='lc', 'll', or 'lp' only", call. = FALSE)
  if (!identical(as.character(bws$type), "fixed"))
    stop("prototype route currently supports fixed bandwidths only", call. = FALSE)
  if (bws$ndim != 2L)
    stop("prototype route currently supports two-dimensional regression surfaces only", call. = FALSE)
  if (bws$nuno != 0L)
    stop("prototype route currently supports continuous/ordered surface variables only", call. = FALSE)
  if (bws$ncon + bws$nord != 2L)
    stop("prototype route requires a two-dimensional regression surface", call. = FALSE)
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("prototype route requires scalar neval >= 2", call. = FALSE)
  invisible(TRUE)
}

.np_plot_proto_clean_regression_data <- function(xdat, ydat) {
  xdat <- toFrame(xdat)
  ydat <- as.vector(ydat)
  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE
  if (!any(keep.rows))
    stop("Data has no rows without NAs")
  list(
    xdat = xdat[keep.rows, , drop = FALSE],
    ydat = ydat[keep.rows]
  )
}

.np_plot_proto_regression_surface_grid <- function(bws,
                                                   xdat,
                                                   neval,
                                                   xtrim) {
  xtrim <- double(bws$ndim) + xtrim

  if (is.ordered(xdat[, 1L])) {
    x1.eval <- bws$xdati$all.ulev[[1L]]
    x1.neval <- length(x1.eval)
  } else {
    x1.neval <- as.integer(neval)
    qi <- trim.quantiles(xdat[, 1L], xtrim[1L])
    x1.eval <- seq(qi[1L], qi[2L], length.out = x1.neval)
  }

  if (is.ordered(xdat[, 2L])) {
    x2.eval <- bws$xdati$all.ulev[[2L]]
    x2.neval <- length(x2.eval)
  } else {
    x2.neval <- as.integer(neval)
    qi <- trim.quantiles(xdat[, 2L], xtrim[2L])
    x2.eval <- seq(qi[1L], qi[2L], length.out = x2.neval)
  }

  x.eval <- expand.grid(x1.eval, x2.eval)
  if (is.ordered(xdat[, 1L]))
    x1.eval <- bws$xdati$all.dlev[[1L]][as.integer(x1.eval)]
  if (is.ordered(xdat[, 2L]))
    x2.eval <- bws$xdati$all.dlev[[2L]][as.integer(x2.eval)]

  list(
    x.eval = x.eval,
    x1.eval = x1.eval,
    x2.eval = x2.eval,
    x1.neval = x1.neval,
    x2.neval = x2.neval
  )
}

.np_plot_proto_npreg_fixed_data <- function(bws,
                                            xdat,
                                            ydat,
                                            neval = 50,
                                            xtrim = 0.0,
                                            plot.errors.method = c("none", "asymptotic", "bootstrap"),
                                            plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                            plot.errors.boot.nonfixed = c("exact", "frozen"),
                                            plot.errors.boot.wild = c("rademacher", "mammen"),
                                            plot.errors.boot.blocklen = NULL,
                                            plot.errors.boot.num = 399L,
                                            plot.errors.center = c("estimate", "bias-corrected"),
                                            plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                 "simultaneous", "all"),
                                            plot.errors.alpha = 0.05,
                                            return.stages = FALSE) {
  if (missing(xdat) || missing(ydat))
    stop("prototype route requires explicit xdat and ydat", call. = FALSE)
  plot.errors.method <- match.arg(plot.errors.method)
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.boot.wild <- match.arg(plot.errors.boot.wild)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)

  dat <- .np_plot_proto_clean_regression_data(xdat = xdat, ydat = ydat)
  xdat <- dat$xdat
  ydat <- dat$ydat
  .np_plot_proto_check_npreg_fixed_surface(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim
  )
  grid <- .np_plot_proto_regression_surface_grid(
    bws = bws,
    xdat = xdat,
    neval = neval,
    xtrim = xtrim
  )
  fit <- .np_plot_regression_eval(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = grid$x.eval,
    gradients = FALSE,
    need.asymptotic = identical(plot.errors.method, "asymptotic")
  )

  terr <- matrix(data = fit$merr, nrow = nrow(grid$x.eval), ncol = 3L)
  terr[, 3L] <- NA
  interval <- NULL
  bootstrap <- NULL
  treg <- matrix(data = fit$mean,
                 nrow = grid$x1.neval,
                 ncol = grid$x2.neval,
                 byrow = FALSE)
  if (identical(plot.errors.method, "bootstrap")) {
    bootstrap <- compute.bootstrap.errors(
      xdat = xdat,
      ydat = ydat,
      exdat = grid$x.eval,
      gradients = FALSE,
      gradient.order = 1L,
      slice.index = 0L,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
      plot.errors.boot.wild = plot.errors.boot.wild,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      progress.target = NULL,
      bws = bws
    )
    terr <- bootstrap$boot.err
    interval <- list(
      err = bootstrap$boot.err[, 1:2, drop = FALSE],
      all.err = bootstrap$boot.all.err
    )
  } else if (identical(plot.errors.method, "asymptotic")) {
    interval <- .np_plot_asymptotic_error_from_se(
      se = fit$merr,
      alpha = plot.errors.alpha,
      band.type = plot.errors.type,
      m = nrow(grid$x.eval)
    )
    terr[, 1:2] <- interval$err
  }

  r1 <- npregression(
    bws = bws,
    eval = grid$x.eval,
    mean = as.double(treg),
    merr = terr[, 1:2, drop = FALSE],
    ntrain = nrow(xdat)
  )
  r1$bias <- NA
  if (identical(plot.errors.center, "bias-corrected"))
    r1$bias <- terr[, 3L] - treg

  plot.data <- list(r1 = r1)
  if (!isTRUE(return.stages))
    return(plot.data)

  list(
    state = list(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      ntrain = nrow(xdat),
      family = "npreg",
      gradients = FALSE
    ),
    target_grid = grid,
    evaluator = fit,
    intervals = if (is.null(interval)) NULL else list(
      method = plot.errors.method,
      type = plot.errors.type,
      alpha = plot.errors.alpha,
      err = interval$err,
      all.err = interval$all.err
    ),
    bootstrap = if (is.null(bootstrap)) NULL else list(
      method = plot.errors.boot.method,
      nonfixed = plot.errors.boot.nonfixed,
      wild = plot.errors.boot.wild,
      blocklen = plot.errors.boot.blocklen,
      B = plot.errors.boot.num,
      center = plot.errors.center,
      boot.err = bootstrap$boot.err,
      boot.all.err = bootstrap$boot.all.err,
      bxp = bootstrap$bxp
    ),
    plot_data = plot.data
  )
}

.np_plot_proto_npreg_fixed_none_data <- function(bws,
                                                 xdat,
                                                 ydat,
                                                 neval = 50,
                                                 xtrim = 0.0,
                                                 return.stages = FALSE) {
  .np_plot_proto_npreg_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    plot.errors.method = "none",
    return.stages = return.stages
  )
}

.np_plot_proto_npreg_fixed_asymptotic_data <- function(bws,
                                                       xdat,
                                                       ydat,
                                                       neval = 50,
                                                       xtrim = 0.0,
                                                       plot.errors.type = c("pmzsd", "pointwise",
                                                                            "bonferroni", "simultaneous",
                                                                            "all"),
                                                       plot.errors.alpha = 0.05,
                                                       return.stages = FALSE) {
  .np_plot_proto_npreg_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    plot.errors.method = "asymptotic",
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    return.stages = return.stages
  )
}

.np_plot_proto_npreg_fixed_bootstrap_data <- function(bws,
                                                     xdat,
                                                     ydat,
                                                     neval = 50,
                                                     xtrim = 0.0,
                                                     plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                                     plot.errors.boot.nonfixed = c("exact", "frozen"),
                                                     plot.errors.boot.wild = c("rademacher", "mammen"),
                                                     plot.errors.boot.blocklen = NULL,
                                                     plot.errors.boot.num = 399L,
                                                     plot.errors.center = c("estimate", "bias-corrected"),
                                                     plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                          "simultaneous", "all"),
                                                     plot.errors.alpha = 0.05,
                                                     return.stages = FALSE) {
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.boot.wild <- match.arg(plot.errors.boot.wild)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)
  .np_plot_proto_npreg_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
    plot.errors.boot.wild = plot.errors.boot.wild,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    return.stages = return.stages
  )
}

.np_plot_proto_check_npindex_fixed_curve <- function(bws,
                                                    xdat,
                                                    ydat,
                                                    neval) {
  if (!inherits(bws, "sibandwidth"))
    stop("prototype route requires a single-index bandwidth object", call. = FALSE)
  regtype.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  if (!is.element(regtype.engine, c("lc", "lp")))
    stop("prototype route currently supports regtype='lc', 'll', or 'lp' only", call. = FALSE)
  if (!identical(as.character(bws$type), "fixed"))
    stop("prototype route currently supports fixed bandwidths only", call. = FALSE)
  if (bws$nuno != 0L)
    stop("prototype route currently supports continuous/ordered index variables only", call. = FALSE)
  if (bws$ncon + bws$nord < 1L)
    stop("prototype route requires at least one continuous/ordered index variable", call. = FALSE)
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("prototype route requires scalar neval >= 2", call. = FALSE)
  invisible(TRUE)
}

.np_plot_proto_npindex_payload <- function(bws,
                                           index,
                                           mean,
                                           ntrain,
                                           trainiseval = FALSE,
                                           merr = NA,
                                           bias = NULL) {
  payload <- singleindex(
    bws = bws,
    index = index,
    mean = mean,
    merr = merr,
    ntrain = ntrain,
    trainiseval = trainiseval,
    gradients = FALSE
  )
  if (!is.null(bias))
    payload$bias <- bias
  payload
}

.np_plot_proto_npindex_fixed_data <- function(bws,
                                             xdat,
                                             ydat,
                                             neval = 50,
                                             plot.errors.method = c("none", "asymptotic", "bootstrap"),
                                             plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                             plot.errors.boot.nonfixed = c("exact", "frozen"),
                                             plot.errors.boot.wild = c("rademacher", "mammen"),
                                             plot.errors.boot.blocklen = NULL,
                                             plot.errors.boot.num = 399L,
                                             plot.errors.center = c("estimate", "bias-corrected"),
                                             plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                  "simultaneous", "all"),
                                             plot.errors.alpha = 0.05,
                                             return.stages = FALSE) {
  if (missing(xdat) || missing(ydat))
    stop("prototype route requires explicit xdat and ydat", call. = FALSE)
  plot.errors.method <- match.arg(plot.errors.method)
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.boot.wild <- match.arg(plot.errors.boot.wild)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)

  dat <- .np_plot_proto_clean_regression_data(xdat = xdat, ydat = ydat)
  xdat <- dat$xdat
  ydat <- dat$ydat
  .np_plot_proto_check_npindex_fixed_curve(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval
  )
  eval.info <- .np_plot_singleindex_eval_grid(
    bws = bws,
    xdat = xdat,
    neval = neval,
    trim = 0.0,
    where = "plot.sibandwidth"
  )
  maxneval <- nrow(eval.info$idx.eval)
  fit <- if (identical(plot.errors.method, "asymptotic")) {
    .np_plot_singleindex_asymptotic_eval(
      bws = bws,
      txdat = xdat,
      tydat = ydat,
      gradients = FALSE,
      index.eval = eval.info$index.eval
    )
  } else {
    .np_plot_singleindex_local_eval(
      bws = bws,
      idx.train = eval.info$idx.train,
      idx.eval = eval.info$idx.eval,
      ydat = ydat,
      gradients = FALSE
    )
  }

  temp.err <- matrix(data = NA_real_, nrow = maxneval, ncol = 3L)
  temp.mean <- as.vector(fit$mean)
  temp.all.err <- NULL
  interval <- NULL
  bootstrap <- NULL
  if (identical(plot.errors.method, "bootstrap")) {
    bootstrap <- compute.bootstrap.errors(
      xdat = xdat,
      ydat = ydat,
      gradients = FALSE,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
      plot.errors.boot.wild = plot.errors.boot.wild,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      idx.eval = eval.info$idx.eval,
      bws = bws
    )
    temp.err[,] <- bootstrap$boot.err
    temp.all.err <- bootstrap$boot.all.err
    interval <- list(
      err = bootstrap$boot.err[, 1:2, drop = FALSE],
      all.err = bootstrap$boot.all.err
    )
  } else if (identical(plot.errors.method, "asymptotic")) {
    interval <- .np_plot_asymptotic_error_from_se(
      se = fit$merr,
      alpha = plot.errors.alpha,
      band.type = plot.errors.type,
      m = maxneval
    )
    temp.err[, 1:2] <- interval$err
    temp.all.err <- interval$all.err
  }

  merr <- NA
  bias <- NULL
  if (!identical(plot.errors.method, "none")) {
    merr <- cbind(-temp.err[, 1L], temp.err[, 2L])
    if (identical(plot.errors.center, "bias-corrected"))
      bias <- temp.err[, 3L] - temp.mean
  }
  plot.data <- list(si1 = .np_plot_proto_npindex_payload(
    bws = bws,
    index = fit$index,
    mean = fit$mean,
    ntrain = nrow(xdat),
    trainiseval = eval.info$trainiseval,
    merr = merr,
    bias = bias
  ))
  if (!isTRUE(return.stages))
    return(plot.data)

  list(
    state = list(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      ntrain = nrow(xdat),
      family = "npindex",
      gradients = FALSE
    ),
    target_grid = eval.info,
    evaluator = fit,
    intervals = if (is.null(interval)) NULL else list(
      method = plot.errors.method,
      type = plot.errors.type,
      alpha = plot.errors.alpha,
      err = interval$err,
      all.err = temp.all.err
    ),
    bootstrap = if (is.null(bootstrap)) NULL else list(
      method = plot.errors.boot.method,
      nonfixed = plot.errors.boot.nonfixed,
      wild = plot.errors.boot.wild,
      blocklen = plot.errors.boot.blocklen,
      B = plot.errors.boot.num,
      center = plot.errors.center,
      boot.err = bootstrap$boot.err,
      boot.all.err = bootstrap$boot.all.err,
      bxp = bootstrap$bxp
    ),
    plot_data = plot.data
  )
}

.np_plot_proto_npindex_fixed_none_data <- function(bws,
                                                  xdat,
                                                  ydat,
                                                  neval = 50,
                                                  return.stages = FALSE) {
  .np_plot_proto_npindex_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    plot.errors.method = "none",
    return.stages = return.stages
  )
}

.np_plot_proto_npindex_fixed_asymptotic_data <- function(bws,
                                                        xdat,
                                                        ydat,
                                                        neval = 50,
                                                        plot.errors.type = c("pmzsd", "pointwise",
                                                                             "bonferroni", "simultaneous",
                                                                             "all"),
                                                        plot.errors.alpha = 0.05,
                                                        return.stages = FALSE) {
  .np_plot_proto_npindex_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    plot.errors.method = "asymptotic",
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    return.stages = return.stages
  )
}

.np_plot_proto_npindex_fixed_bootstrap_data <- function(bws,
                                                       xdat,
                                                       ydat,
                                                       neval = 50,
                                                       plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                                       plot.errors.boot.nonfixed = c("exact", "frozen"),
                                                       plot.errors.boot.wild = c("rademacher", "mammen"),
                                                       plot.errors.boot.blocklen = NULL,
                                                       plot.errors.boot.num = 399L,
                                                       plot.errors.center = c("estimate", "bias-corrected"),
                                                       plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                            "simultaneous", "all"),
                                                       plot.errors.alpha = 0.05,
                                                       return.stages = FALSE) {
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.boot.wild <- match.arg(plot.errors.boot.wild)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)
  .np_plot_proto_npindex_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
    plot.errors.boot.wild = plot.errors.boot.wild,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    return.stages = return.stages
  )
}

.np_plot_proto_check_npscoef_fixed_surface <- function(bws,
                                                       xdat,
                                                       ydat,
                                                       zdat,
                                                       neval,
                                                       xtrim,
                                                       ztrim) {
  if (!inherits(bws, "scbandwidth"))
    stop("prototype route requires a smooth-coefficient bandwidth object", call. = FALSE)
  regtype <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  if (!is.element(regtype, c("lc", "ll", "lp")))
    stop("prototype route currently supports regtype='lc', 'll', or 'lp' only", call. = FALSE)
  if (!identical(as.character(bws$type), "fixed"))
    stop("prototype route currently supports fixed bandwidths only", call. = FALSE)
  if (ncol(xdat) != 1L || ncol(zdat) != 1L)
    stop("prototype route currently supports one x variable and one z variable", call. = FALSE)
  if (bws$xdati$iuno[1L] || bws$zdati$iuno[1L])
    stop("prototype route currently supports continuous/ordered surface variables only", call. = FALSE)
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("prototype route requires scalar neval >= 2", call. = FALSE)
  invisible(TRUE)
}

.np_plot_proto_clean_scoef_data <- function(xdat, ydat, zdat) {
  xdat <- toFrame(xdat)
  zdat <- toFrame(zdat)
  ydat <- as.vector(ydat)
  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat, zdat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE
  if (!any(keep.rows))
    stop("Data has no rows without NAs")
  list(
    xdat = xdat[keep.rows, , drop = FALSE],
    ydat = ydat[keep.rows],
    zdat = zdat[keep.rows, , drop = FALSE]
  )
}

.np_plot_proto_scoef_surface_grid <- function(bws,
                                              xdat,
                                              zdat,
                                              neval,
                                              xtrim,
                                              ztrim) {
  xtrim <- double(ncol(xdat)) + xtrim
  ztrim <- double(ncol(zdat)) + ztrim

  if (is.ordered(xdat[, 1L])) {
    x1.eval <- bws$xdati$all.ulev[[1L]]
    x1.neval <- length(x1.eval)
  } else {
    x1.neval <- as.integer(neval)
    qi <- trim.quantiles(xdat[, 1L], xtrim[1L])
    x1.eval <- seq(qi[1L], qi[2L], length.out = x1.neval)
  }

  if (is.ordered(zdat[, 1L])) {
    x2.eval <- bws$zdati$all.ulev[[1L]]
    x2.neval <- length(x2.eval)
  } else {
    x2.neval <- as.integer(neval)
    qi <- trim.quantiles(zdat[, 1L], ztrim[1L])
    x2.eval <- seq(qi[1L], qi[2L], length.out = x2.neval)
  }

  x.eval <- expand.grid(x1.eval, x2.eval)
  colnames(x.eval) <- c(colnames(xdat)[1L], colnames(zdat)[1L])
  if (is.ordered(xdat[, 1L]))
    x1.eval <- bws$xdati$all.dlev[[1L]][as.integer(x1.eval)]
  if (is.ordered(zdat[, 1L]))
    x2.eval <- bws$zdati$all.dlev[[1L]][as.integer(x2.eval)]

  list(
    x.eval = x.eval,
    exdat = x.eval[, 1L, drop = FALSE],
    ezdat = x.eval[, 2L, drop = FALSE],
    x1.eval = x1.eval,
    x2.eval = x2.eval,
    x1.neval = x1.neval,
    x2.neval = x2.neval
  )
}

.np_plot_proto_npscoef_fixed_data <- function(bws,
                                             xdat,
                                             ydat,
                                             zdat,
                                             neval = 50,
                                             xtrim = 0.0,
                                             ztrim = 0.0,
                                             plot.errors.method = c("none", "asymptotic", "bootstrap"),
                                             plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                             plot.errors.boot.nonfixed = c("exact", "frozen"),
                                             plot.errors.boot.wild = c("rademacher", "mammen"),
                                             plot.errors.boot.blocklen = NULL,
                                             plot.errors.boot.num = 399L,
                                             plot.errors.center = c("estimate", "bias-corrected"),
                                             plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                  "simultaneous", "all"),
                                             plot.errors.alpha = 0.05,
                                             return.stages = FALSE) {
  if (missing(xdat) || missing(ydat) || missing(zdat))
    stop("prototype route requires explicit xdat, ydat, and zdat", call. = FALSE)
  plot.errors.method <- match.arg(plot.errors.method)
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.boot.wild <- match.arg(plot.errors.boot.wild)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)

  dat <- .np_plot_proto_clean_scoef_data(xdat = xdat, ydat = ydat, zdat = zdat)
  xdat <- dat$xdat
  ydat <- dat$ydat
  zdat <- dat$zdat
  .np_plot_proto_check_npscoef_fixed_surface(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim
  )
  grid <- .np_plot_proto_scoef_surface_grid(
    bws = bws,
    xdat = xdat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim
  )
  fit <- .np_scoef_fit_internal(
    bws = bws,
    txdat = xdat,
    tydat = ydat,
    tzdat = zdat,
    exdat = grid$exdat,
    ezdat = grid$ezdat,
    iterate = FALSE,
    errors = identical(plot.errors.method, "asymptotic"),
    betas = FALSE
  )

  treg <- matrix(data = fit$mean,
                 nrow = grid$x1.neval,
                 ncol = grid$x2.neval,
                 byrow = FALSE)
  terr <- matrix(data = fit$merr, nrow = nrow(grid$x.eval), ncol = 3L)
  terr[, 3L] <- NA_real_
  interval <- NULL
  bootstrap <- NULL
  if (identical(plot.errors.method, "bootstrap")) {
    bootstrap <- compute.bootstrap.errors(
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      exdat = grid$exdat,
      ezdat = grid$ezdat,
      gradients = FALSE,
      slice.index = 0L,
      progress.target = NULL,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
      plot.errors.boot.wild = plot.errors.boot.wild,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      bws = bws
    )
    terr <- bootstrap$boot.err
    interval <- list(
      err = bootstrap$boot.err[, 1:2, drop = FALSE],
      all.err = bootstrap$boot.all.err
    )
  } else if (identical(plot.errors.method, "asymptotic")) {
    interval <- .np_plot_asymptotic_error_from_se(
      se = fit$merr,
      alpha = plot.errors.alpha,
      band.type = plot.errors.type,
      m = nrow(grid$x.eval)
    )
    terr[, 1:2] <- interval$err
  }

  r1 <- smoothcoefficient(
    bws = bws,
    eval = list(exdat = grid$exdat, ezdat = grid$ezdat),
    mean = as.double(treg),
    merr = if (identical(plot.errors.method, "none")) {
      matrix(NA, nrow = nrow(grid$x.eval), ncol = 2L)
    } else {
      terr[, 1:2, drop = FALSE]
    },
    ntrain = nrow(xdat)
  )
  r1$bias <- NA
  if (identical(plot.errors.center, "bias-corrected"))
    r1$bias <- terr[, 3L] - treg

  plot.data <- list(r1 = r1)
  if (!isTRUE(return.stages))
    return(plot.data)

  list(
    state = list(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      ntrain = nrow(xdat),
      family = "npscoef",
      gradients = FALSE,
      coef = FALSE
    ),
    target_grid = grid,
    evaluator = fit,
    intervals = if (is.null(interval)) NULL else list(
      method = plot.errors.method,
      type = plot.errors.type,
      alpha = plot.errors.alpha,
      err = interval$err,
      all.err = interval$all.err
    ),
    bootstrap = if (is.null(bootstrap)) NULL else list(
      method = plot.errors.boot.method,
      nonfixed = plot.errors.boot.nonfixed,
      wild = plot.errors.boot.wild,
      blocklen = plot.errors.boot.blocklen,
      B = plot.errors.boot.num,
      center = plot.errors.center,
      boot.err = bootstrap$boot.err,
      boot.all.err = bootstrap$boot.all.err,
      bxp = bootstrap$bxp
    ),
    plot_data = plot.data
  )
}

.np_plot_proto_npscoef_fixed_none_data <- function(bws,
                                                  xdat,
                                                  ydat,
                                                  zdat,
                                                  neval = 50,
                                                  xtrim = 0.0,
                                                  ztrim = 0.0,
                                                  return.stages = FALSE) {
  .np_plot_proto_npscoef_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim,
    plot.errors.method = "none",
    return.stages = return.stages
  )
}

.np_plot_proto_npscoef_fixed_asymptotic_data <- function(bws,
                                                        xdat,
                                                        ydat,
                                                        zdat,
                                                        neval = 50,
                                                        xtrim = 0.0,
                                                        ztrim = 0.0,
                                                        plot.errors.type = c("pmzsd", "pointwise",
                                                                             "bonferroni", "simultaneous",
                                                                             "all"),
                                                        plot.errors.alpha = 0.05,
                                                        return.stages = FALSE) {
  .np_plot_proto_npscoef_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim,
    plot.errors.method = "asymptotic",
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    return.stages = return.stages
  )
}

.np_plot_proto_npscoef_fixed_bootstrap_data <- function(bws,
                                                       xdat,
                                                       ydat,
                                                       zdat,
                                                       neval = 50,
                                                       xtrim = 0.0,
                                                       ztrim = 0.0,
                                                       plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                                       plot.errors.boot.nonfixed = c("exact", "frozen"),
                                                       plot.errors.boot.wild = c("rademacher", "mammen"),
                                                       plot.errors.boot.blocklen = NULL,
                                                       plot.errors.boot.num = 399L,
                                                       plot.errors.center = c("estimate", "bias-corrected"),
                                                       plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                            "simultaneous", "all"),
                                                       plot.errors.alpha = 0.05,
                                                       return.stages = FALSE) {
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.boot.wild <- match.arg(plot.errors.boot.wild)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)
  .np_plot_proto_npscoef_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
    plot.errors.boot.wild = plot.errors.boot.wild,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    return.stages = return.stages
  )
}

.np_plot_proto_check_npplreg_fixed_surface <- function(bws,
                                                       xdat,
                                                       ydat,
                                                       zdat,
                                                       neval,
                                                       xtrim,
                                                       ztrim) {
  if (!inherits(bws, "plbandwidth"))
    stop("prototype route requires a partially linear bandwidth object", call. = FALSE)
  regtype <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  if (!is.element(regtype, c("lc", "ll", "lp")))
    stop("prototype route currently supports regtype='lc', 'll', or 'lp' only", call. = FALSE)
  if (!identical(as.character(bws$type), "fixed"))
    stop("prototype route currently supports fixed bandwidths only", call. = FALSE)
  if (ncol(xdat) != 1L || ncol(zdat) != 1L)
    stop("prototype route currently supports one x variable and one z variable", call. = FALSE)
  if (bws$xdati$iuno[1L] || bws$zdati$iuno[1L])
    stop("prototype route currently supports continuous/ordered surface variables only", call. = FALSE)
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("prototype route requires scalar neval >= 2", call. = FALSE)
  invisible(TRUE)
}

.np_plot_proto_npplreg_fixed_data <- function(bws,
                                             xdat,
                                             ydat,
                                             zdat,
                                             neval = 50,
                                             xtrim = 0.0,
                                             ztrim = 0.0,
                                             plot.errors.method = c("none", "asymptotic", "bootstrap"),
                                             plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                             plot.errors.boot.nonfixed = c("exact", "frozen"),
                                             plot.errors.boot.wild = c("rademacher", "mammen"),
                                             plot.errors.boot.blocklen = NULL,
                                             plot.errors.boot.num = 399L,
                                             plot.errors.center = c("estimate", "bias-corrected"),
                                             plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                  "simultaneous", "all"),
                                             plot.errors.alpha = 0.05,
                                             return.stages = FALSE) {
  if (missing(xdat) || missing(ydat) || missing(zdat))
    stop("prototype route requires explicit xdat, ydat, and zdat", call. = FALSE)
  plot.errors.method <- match.arg(plot.errors.method)
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.boot.wild <- match.arg(plot.errors.boot.wild)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)

  dat <- .np_plot_proto_clean_scoef_data(xdat = xdat, ydat = ydat, zdat = zdat)
  xdat <- dat$xdat
  ydat <- dat$ydat
  zdat <- dat$zdat
  .np_plot_proto_check_npplreg_fixed_surface(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim
  )
  grid <- .np_plot_proto_scoef_surface_grid(
    bws = bws,
    xdat = xdat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim
  )
  fit <- if (identical(plot.errors.method, "asymptotic")) {
    .np_plot_plreg_asymptotic_fit(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      exdat = grid$exdat,
      ezdat = grid$ezdat
    )
  } else {
    .np_plot_plreg_local_fit(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      exdat = grid$exdat,
      ezdat = grid$ezdat
    )
  }

  treg <- matrix(data = fit$mean,
                 nrow = grid$x1.neval,
                 ncol = grid$x2.neval,
                 byrow = FALSE)
  terr <- matrix(data = NA, nrow = nrow(grid$x.eval), ncol = 3L)
  interval <- NULL
  bootstrap <- NULL
  if (identical(plot.errors.method, "bootstrap")) {
    bootstrap <- compute.bootstrap.errors(
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      exdat = grid$exdat,
      ezdat = grid$ezdat,
      gradients = FALSE,
      slice.index = 0L,
      progress.target = NULL,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
      plot.errors.boot.wild = plot.errors.boot.wild,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      bws = bws
    )
    terr <- bootstrap$boot.err
    interval <- list(
      err = bootstrap$boot.err[, 1:2, drop = FALSE],
      all.err = bootstrap$boot.all.err
    )
  } else if (identical(plot.errors.method, "asymptotic")) {
    interval <- .np_plot_asymptotic_error_from_se(
      se = fit$merr,
      alpha = plot.errors.alpha,
      band.type = plot.errors.type,
      m = nrow(grid$x.eval)
    )
    terr[, 1:2] <- interval$err
  }

  r1 <- plregression(
    bws = bws,
    xcoef = fit$xcoef,
    xcoefvcov = vcov(fit),
    xcoeferr = fit$xcoeferr,
    evalx = grid$exdat[, 1L],
    evalz = grid$ezdat[, 1L],
    mean = fit$mean,
    ntrain = nrow(xdat),
    trainiseval = FALSE,
    xtra = c(fit$RSQ, fit$MSE, 0, 0, 0, 0)
  )
  r1$merr <- NA
  r1$bias <- NA
  if (!identical(plot.errors.method, "none"))
    r1$merr <- terr[, 1:2, drop = FALSE]
  if (identical(plot.errors.center, "bias-corrected"))
    r1$bias <- terr[, 3L] - treg

  plot.data <- list(r1 = r1)
  if (!isTRUE(return.stages))
    return(plot.data)

  list(
    state = list(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      ntrain = nrow(xdat),
      family = "npplreg",
      gradients = FALSE,
      coef = FALSE
    ),
    target_grid = grid,
    evaluator = fit,
    intervals = if (is.null(interval)) NULL else list(
      method = plot.errors.method,
      type = plot.errors.type,
      alpha = plot.errors.alpha,
      err = interval$err,
      all.err = interval$all.err
    ),
    bootstrap = if (is.null(bootstrap)) NULL else list(
      method = plot.errors.boot.method,
      nonfixed = plot.errors.boot.nonfixed,
      wild = plot.errors.boot.wild,
      blocklen = plot.errors.boot.blocklen,
      B = plot.errors.boot.num,
      center = plot.errors.center,
      boot.err = bootstrap$boot.err,
      boot.all.err = bootstrap$boot.all.err,
      bxp = bootstrap$bxp
    ),
    plot_data = plot.data
  )
}

.np_plot_proto_npplreg_fixed_none_data <- function(bws,
                                                  xdat,
                                                  ydat,
                                                  zdat,
                                                  neval = 50,
                                                  xtrim = 0.0,
                                                  ztrim = 0.0,
                                                  return.stages = FALSE) {
  .np_plot_proto_npplreg_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim,
    plot.errors.method = "none",
    return.stages = return.stages
  )
}

.np_plot_proto_npplreg_fixed_asymptotic_data <- function(bws,
                                                        xdat,
                                                        ydat,
                                                        zdat,
                                                        neval = 50,
                                                        xtrim = 0.0,
                                                        ztrim = 0.0,
                                                        plot.errors.type = c("pmzsd", "pointwise",
                                                                             "bonferroni", "simultaneous",
                                                                             "all"),
                                                        plot.errors.alpha = 0.05,
                                                        return.stages = FALSE) {
  .np_plot_proto_npplreg_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim,
    plot.errors.method = "asymptotic",
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    return.stages = return.stages
  )
}

.np_plot_proto_npplreg_fixed_bootstrap_data <- function(bws,
                                                       xdat,
                                                       ydat,
                                                       zdat,
                                                       neval = 50,
                                                       xtrim = 0.0,
                                                       ztrim = 0.0,
                                                       plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                                       plot.errors.boot.nonfixed = c("exact", "frozen"),
                                                       plot.errors.boot.wild = c("rademacher", "mammen"),
                                                       plot.errors.boot.blocklen = NULL,
                                                       plot.errors.boot.num = 399L,
                                                       plot.errors.center = c("estimate", "bias-corrected"),
                                                       plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                            "simultaneous", "all"),
                                                       plot.errors.alpha = 0.05,
                                                       return.stages = FALSE) {
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.boot.wild <- match.arg(plot.errors.boot.wild)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)
  .np_plot_proto_npplreg_fixed_data(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    neval = neval,
    xtrim = xtrim,
    ztrim = ztrim,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
    plot.errors.boot.wild = plot.errors.boot.wild,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    return.stages = return.stages
  )
}

.np_plot_proto_check_npqreg_fixed_slices <- function(object,
                                                     xdat,
                                                     ydat,
                                                     neval,
                                                     xtrim) {
  if (!inherits(object, "qregression"))
    stop("prototype route requires a quantile regression object", call. = FALSE)
  bws <- object$bws
  if (!(inherits(bws, "condbandwidth") || inherits(bws, "conbandwidth")))
    stop("prototype route requires a conditional distribution bandwidth object", call. = FALSE)
  if (!identical(as.character(bws$type), "fixed"))
    stop("prototype route currently supports fixed bandwidths only", call. = FALSE)
  if (bws$xnuno != 0L)
    stop("prototype route currently supports continuous/ordered slice variables only", call. = FALSE)
  if (bws$xndim < 1L || bws$yndim != 1L)
    stop("prototype route requires at least one x variable and one y variable", call. = FALSE)
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("prototype route requires scalar neval >= 2", call. = FALSE)
  invisible(TRUE)
}

.np_plot_proto_npqreg_slice_targets <- function(bws,
                                                xdat,
                                                ydat,
                                                neval,
                                                xq,
                                                yq,
                                                xtrim) {
  xq <- double(bws$xndim) + xq
  xtrim <- double(bws$xndim) + xtrim
  x.ev <- xdat[1L, , drop = FALSE]
  y.ev <- ydat[1L, , drop = FALSE]
  for (i in seq_len(bws$xndim))
    x.ev[1L, i] <- uocquantile(xdat[, i], prob = xq[i])
  for (i in seq_len(bws$yndim))
    y.ev[1L, i] <- uocquantile(ydat[, i], prob = yq)

  maxneval <- max(c(sapply(xdat, nlevels), neval))
  exdat.base <- xdat[rep(1L, maxneval), , drop = FALSE]
  eydat.base <- as.data.frame(matrix(data = 0, nrow = maxneval, ncol = bws$yndim))
  colnames(eydat.base) <- colnames(ydat)
  for (i in seq_len(bws$xndim))
    exdat.base[, i] <- x.ev[1L, i]
  for (i in seq_len(bws$yndim))
    eydat.base[, i] <- y.ev[1L, i]

  slices <- vector("list", bws$xndim)
  for (i in seq_len(bws$xndim)) {
    if (is.factor(xdat[, i])) {
      ei <- levels(xdat[, i])
      ei <- factor(ei, levels = ei)
      xi.neval <- length(ei)
    } else {
      xi.neval <- as.integer(neval)
      qi <- trim.quantiles(xdat[, i], xtrim[i])
      ei <- seq(qi[1L], qi[2L], length.out = neval)
    }
    exdat.i <- subcol(exdat.base, ei, i)[seq_len(xi.neval), , drop = FALSE]
    eydat.i <- eydat.base[seq_len(xi.neval), , drop = FALSE]
    slices[[i]] <- list(
      index = i,
      eval = ei,
      neval = xi.neval,
      exdat = exdat.i,
      eydat = eydat.i
    )
  }

  list(
    base_x = x.ev,
    base_y = y.ev,
    maxneval = maxneval,
    slices = slices
  )
}

.np_plot_proto_npqreg_fixed_slices_data <- function(object,
                                                   xdat,
                                                   ydat,
                                                   neval = 50,
                                                   xq = 0.5,
                                                   yq = 0.5,
                                                   xtrim = 0.0,
                                                   plot.errors.method = c("none", "asymptotic", "bootstrap"),
                                                   plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                                   plot.errors.boot.nonfixed = c("exact", "frozen"),
                                                   plot.errors.boot.blocklen = NULL,
                                                   plot.errors.boot.num = 399L,
                                                   plot.errors.center = c("estimate", "bias-corrected"),
                                                   plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                        "simultaneous", "all"),
                                                   plot.errors.alpha = 0.05,
                                                   gradients = FALSE,
                                                   return.stages = FALSE) {
  if (missing(xdat) || missing(ydat))
    stop("prototype route requires explicit xdat and ydat", call. = FALSE)
  plot.errors.method <- match.arg(plot.errors.method)
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  bws <- object$bws
  dat <- .np_plot_proto_clean_conditional_data(xdat = xdat, ydat = ydat)
  xdat <- dat$xdat
  ydat <- dat$ydat
  .np_plot_proto_check_npqreg_fixed_slices(
    object = object,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim
  )
  targets <- .np_plot_proto_npqreg_slice_targets(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xq = xq,
    yq = yq,
    xtrim = xtrim
  )

  plot.data <- vector("list", length(targets$slices))
  evaluators <- vector("list", length(targets$slices))
  intervals <- vector("list", length(targets$slices))
  bootstraps <- vector("list", length(targets$slices))
  for (ii in seq_along(targets$slices)) {
    target <- targets$slices[[ii]]
    fit <- .np_plot_quantile_eval(
      txdat = xdat,
      tydat = ydat,
      exdat = target$exdat,
      tau = object$tau,
      gradients = gradients,
      bws = bws
    )
    temp.err <- matrix(data = NA, nrow = targets$maxneval, ncol = 3L)
    temp.quant <- rep(NA_real_, targets$maxneval)
    temp.quant[seq_len(target$neval)] <- fit$quantile
    interval <- NULL
    bootstrap <- NULL

    if (identical(plot.errors.method, "bootstrap") && !isTRUE(gradients)) {
      bootstrap <- compute.bootstrap.errors(
        xdat = xdat,
        ydat = ydat,
        exdat = target$exdat,
        eydat = target$eydat,
        cdf = TRUE,
        quantreg = TRUE,
        tau = object$tau,
        gradients = FALSE,
        gradient.index = 0L,
        slice.index = target$index,
        plot.errors.boot.method = plot.errors.boot.method,
        plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
        plot.errors.boot.blocklen = plot.errors.boot.blocklen,
        plot.errors.boot.num = plot.errors.boot.num,
        plot.errors.center = plot.errors.center,
        plot.errors.type = plot.errors.type,
        plot.errors.alpha = plot.errors.alpha,
        progress.target = .np_plot_conditional_bootstrap_target_label(
          bws = bws,
          slice.index = target$index,
          gradients = FALSE,
          gradient.index = 0L
        ),
        bws = bws
      )
      temp.err[seq_len(target$neval), ] <- bootstrap$boot.err
      interval <- list(
        err = bootstrap$boot.err[, 1:2, drop = FALSE],
        all.err = bootstrap$boot.all.err
      )
    } else if (identical(plot.errors.method, "bootstrap") && isTRUE(gradients)) {
      bootstrap <- vector("list", bws$xndim)
      interval <- vector("list", bws$xndim)
      for (jj in seq_len(bws$xndim)) {
        temp.err.j <- matrix(data = NA, nrow = targets$maxneval, ncol = 3L)
        bootstrap.j <- compute.bootstrap.errors(
          xdat = xdat,
          ydat = ydat,
          exdat = target$exdat,
          eydat = target$eydat,
          cdf = TRUE,
          quantreg = TRUE,
          tau = object$tau,
          gradients = TRUE,
          gradient.index = jj,
          slice.index = target$index,
          plot.errors.boot.method = plot.errors.boot.method,
          plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
          plot.errors.boot.blocklen = plot.errors.boot.blocklen,
          plot.errors.boot.num = plot.errors.boot.num,
          plot.errors.center = plot.errors.center,
          plot.errors.type = plot.errors.type,
          plot.errors.alpha = plot.errors.alpha,
          progress.target = .np_plot_conditional_bootstrap_target_label(
            bws = bws,
            slice.index = target$index,
            gradients = TRUE,
            gradient.index = jj
          ),
          bws = bws
        )
        temp.err.j[seq_len(target$neval), ] <- bootstrap.j$boot.err
        bootstrap[[jj]] <- bootstrap.j
        interval[[jj]] <- list(
          err = bootstrap.j$boot.err[, 1:2, drop = FALSE],
          all.err = bootstrap.j$boot.all.err
        )
      }
    } else if (identical(plot.errors.method, "asymptotic")) {
      asym <- .np_plot_asymptotic_error_from_se(
        se = fit$quanterr,
        alpha = plot.errors.alpha,
        band.type = plot.errors.type,
        m = target$neval
      )
      temp.err[seq_len(target$neval), 1:2] <- asym$err
      interval <- list(err = asym$err, all.err = asym$all.err)
    }

    payload <- fit
    if (!identical(plot.errors.method, "none") && !isTRUE(gradients)) {
      payload$quanterr <- na.omit(cbind(
        -temp.err[seq_len(target$neval), 1L],
        temp.err[seq_len(target$neval), 2L]
      ))
      payload$bias <- na.omit(
        temp.quant[seq_len(target$neval)] - temp.err[seq_len(target$neval), 3L]
      )
      payload$bxp <- if (is.null(bootstrap)) list() else bootstrap$bxp
    } else if (identical(plot.errors.method, "bootstrap") && isTRUE(gradients)) {
      for (jj in seq_len(bws$xndim)) {
        temp.err.j <- matrix(data = NA, nrow = targets$maxneval, ncol = 3L)
        temp.err.j[seq_len(target$neval), ] <- bootstrap[[jj]]$boot.err
        temp.grad.j <- rep(NA_real_, targets$maxneval)
        temp.grad.j[seq_len(target$neval)] <- fit$quantgrad[, jj]
        payload[[paste0("gc", jj, "err")]] <- na.omit(cbind(
          -temp.err.j[, 1L],
          temp.err.j[, 2L]
        ))
        payload[[paste0("gc", jj, "bias")]] <- na.omit(
          temp.grad.j - temp.err.j[, 3L]
        )
        payload$bxp <- bootstrap[[jj]]$bxp
      }
    }

    plot.data[[ii]] <- payload
    evaluators[[ii]] <- fit
    intervals[[ii]] <- interval
    bootstraps[[ii]] <- bootstrap
  }
  names(plot.data) <- paste0("cd", seq_along(plot.data))

  if (!isTRUE(return.stages))
    return(plot.data)

  list(
    state = list(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      ntrain = nrow(xdat),
      family = "npqreg",
      tau = object$tau,
      gradients = gradients
    ),
    target_grid = targets,
    evaluator = evaluators,
    intervals = intervals,
    bootstrap = bootstraps,
    plot_data = plot.data
  )
}

.np_plot_proto_npqreg_fixed_none_data <- function(object,
                                                 xdat,
                                                 ydat,
                                                 neval = 50,
                                                 xq = 0.5,
                                                 yq = 0.5,
                                                 xtrim = 0.0,
                                                 gradients = FALSE,
                                                 return.stages = FALSE) {
  .np_plot_proto_npqreg_fixed_slices_data(
    object = object,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xq = xq,
    yq = yq,
    xtrim = xtrim,
    plot.errors.method = "none",
    gradients = gradients,
    return.stages = return.stages
  )
}

.np_plot_proto_npqreg_fixed_asymptotic_data <- function(object,
                                                       xdat,
                                                       ydat,
                                                       neval = 50,
                                                       xq = 0.5,
                                                       yq = 0.5,
                                                       xtrim = 0.0,
                                                       plot.errors.type = c("pmzsd", "pointwise",
                                                                            "bonferroni", "simultaneous",
                                                                            "all"),
                                                       plot.errors.alpha = 0.05,
                                                       gradients = FALSE,
                                                       return.stages = FALSE) {
  .np_plot_proto_npqreg_fixed_slices_data(
    object = object,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xq = xq,
    yq = yq,
    xtrim = xtrim,
    plot.errors.method = "asymptotic",
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    gradients = gradients,
    return.stages = return.stages
  )
}

.np_plot_proto_npqreg_fixed_bootstrap_data <- function(object,
                                                      xdat,
                                                      ydat,
                                                      neval = 50,
                                                      xq = 0.5,
                                                      yq = 0.5,
                                                      xtrim = 0.0,
                                                      plot.errors.boot.method = c("wild", "inid", "fixed", "geom"),
                                                      plot.errors.boot.nonfixed = c("exact", "frozen"),
                                                      plot.errors.boot.blocklen = NULL,
                                                      plot.errors.boot.num = 399L,
                                                      plot.errors.center = c("estimate", "bias-corrected"),
                                                      plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                           "simultaneous", "all"),
                                                      plot.errors.alpha = 0.05,
                                                      gradients = FALSE,
                                                      return.stages = FALSE) {
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)
  .np_plot_proto_npqreg_fixed_slices_data(
    object = object,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xq = xq,
    yq = yq,
    xtrim = xtrim,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = plot.errors.boot.method,
    plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
    plot.errors.boot.blocklen = plot.errors.boot.blocklen,
    plot.errors.boot.num = plot.errors.boot.num,
    plot.errors.center = plot.errors.center,
    plot.errors.type = plot.errors.type,
    plot.errors.alpha = plot.errors.alpha,
    gradients = gradients,
    return.stages = return.stages
  )
}

.np_plot_proto_check_unconditional_fixed_surface <- function(bws,
                                                            xdat,
                                                            neval,
                                                            xtrim,
                                                            cdf = FALSE) {
  expected.class <- if (isTRUE(cdf)) "dbandwidth" else "bandwidth"
  expected.label <- if (isTRUE(cdf)) "unconditional distribution" else "unconditional density"
  if (!inherits(bws, expected.class))
    stop(sprintf("prototype route requires an %s bandwidth object", expected.label), call. = FALSE)
  if (!identical(as.character(bws$type), "fixed"))
    stop("prototype route currently supports fixed bandwidths only", call. = FALSE)
  if (bws$ndim != 2L)
    stop(sprintf("prototype route currently supports two-dimensional %s surfaces", expected.label), call. = FALSE)
  if (bws$nuno != 0L)
    stop("prototype route currently supports continuous/ordered surface variables only", call. = FALSE)
  if (bws$ncon + bws$nord != 2L)
    stop(sprintf("prototype route requires a two-dimensional %s surface", expected.label), call. = FALSE)
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("prototype route requires scalar neval >= 2", call. = FALSE)
  invisible(TRUE)
}

.np_plot_proto_clean_unconditional_data <- function(xdat) {
  xdat <- toFrame(xdat)
  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(xdat), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE
  if (!any(keep.rows))
    stop("Data has no rows without NAs")
  list(xdat = xdat[keep.rows, , drop = FALSE])
}

.np_plot_proto_unconditional_fixed_data <- function(bws,
                                                   xdat,
                                                   neval = 50,
                                                   xtrim = 0.0,
                                                   plot.errors.method = c("none", "asymptotic", "bootstrap"),
                                                   plot.errors.boot.method = c("inid", "fixed", "geom"),
                                                   plot.errors.boot.nonfixed = c("exact", "frozen"),
                                                   plot.errors.boot.blocklen = NULL,
                                                   plot.errors.boot.num = 1999,
                                                   plot.errors.center = c("estimate", "bias-corrected"),
                                                   plot.errors.type = c("pmzsd", "pointwise", "bonferroni",
                                                                        "simultaneous", "all"),
                                                   plot.errors.alpha = 0.05,
                                                   return.stages = FALSE,
                                                   cdf = FALSE) {
  if (missing(xdat))
    stop("prototype route requires explicit xdat", call. = FALSE)
  plot.errors.method <- match.arg(plot.errors.method)
  plot.errors.boot.method <- match.arg(plot.errors.boot.method)
  plot.errors.boot.nonfixed <- match.arg(plot.errors.boot.nonfixed)
  plot.errors.center <- match.arg(plot.errors.center)
  plot.errors.type <- match.arg(plot.errors.type)

  dat <- .np_plot_proto_clean_unconditional_data(xdat = xdat)
  xdat <- dat$xdat
  .np_plot_proto_check_unconditional_fixed_surface(
    bws = bws,
    xdat = xdat,
    neval = neval,
    xtrim = xtrim,
    cdf = cdf
  )
  grid <- .np_plot_proto_regression_surface_grid(
    bws = bws,
    xdat = xdat,
    neval = neval,
    xtrim = xtrim
  )
  fit <- .np_plot_unconditional_eval(
    xdat = xdat,
    exdat = grid$x.eval,
    bws = bws,
    cdf = isTRUE(cdf),
    need.asymptotic = identical(plot.errors.method, "asymptotic")
  )

  estimate <- if (isTRUE(cdf)) fit$dist else fit$dens
  se <- if (isTRUE(cdf)) fit$derr else fit$derr
  terr <- matrix(data = se, nrow = length(estimate), ncol = 3L)
  terr[, 3L] <- NA_real_
  interval <- NULL
  bootstrap <- NULL
  if (identical(plot.errors.method, "bootstrap")) {
    bootstrap.args <- list(
      xdat = xdat,
      exdat = grid$x.eval,
      slice.index = 0L,
      plot.errors.boot.method = plot.errors.boot.method,
      plot.errors.boot.nonfixed = plot.errors.boot.nonfixed,
      plot.errors.boot.blocklen = plot.errors.boot.blocklen,
      plot.errors.boot.num = plot.errors.boot.num,
      plot.errors.center = plot.errors.center,
      plot.errors.type = plot.errors.type,
      plot.errors.alpha = plot.errors.alpha,
      bws = bws
    )
    if (!isTRUE(cdf))
      bootstrap.args$cdf <- FALSE
    bootstrap <- do.call(compute.bootstrap.errors, bootstrap.args)
    terr <- bootstrap$boot.err
    interval <- list(
      err = bootstrap$boot.err[, 1:2, drop = FALSE],
      all.err = bootstrap$boot.all.err
    )
  } else if (identical(plot.errors.method, "asymptotic")) {
    interval <- .np_plot_asymptotic_error_from_se(
      se = se,
      alpha = plot.errors.alpha,
      band.type = plot.errors.type,
      m = nrow(grid$x.eval)
    )
    terr[, 1:2] <- interval$err
  }

  d1 <- if (isTRUE(cdf)) {
    npdistribution(
      bws = bws,
      eval = grid$x.eval,
      dist = estimate,
      derr = terr[, 1:2, drop = FALSE],
      ntrain = nrow(xdat)
    )
  } else {
    npdensity(
      bws = bws,
      eval = grid$x.eval,
      dens = estimate,
      derr = terr[, 1:2, drop = FALSE],
      ntrain = nrow(xdat)
    )
  }
  d1$bias <- NA
  if (identical(plot.errors.center, "bias-corrected"))
    d1$bias <- terr[, 3L] - estimate

  plot.data <- list(d1 = d1)
  if (!isTRUE(return.stages))
    return(plot.data)

  list(
    state = list(
      bws = bws,
      xdat = xdat,
      ntrain = nrow(xdat),
      family = if (isTRUE(cdf)) "npudist" else "npudens",
      cdf = isTRUE(cdf)
    ),
    target_grid = grid,
    evaluator = fit,
    intervals = if (is.null(interval)) NULL else list(
      method = plot.errors.method,
      type = plot.errors.type,
      alpha = plot.errors.alpha,
      err = interval$err,
      all.err = interval$all.err
    ),
    bootstrap = if (is.null(bootstrap)) NULL else list(
      method = plot.errors.boot.method,
      nonfixed = plot.errors.boot.nonfixed,
      blocklen = plot.errors.boot.blocklen,
      B = plot.errors.boot.num,
      center = plot.errors.center,
      boot.err = bootstrap$boot.err,
      boot.all.err = bootstrap$boot.all.err,
      bxp = bootstrap$bxp
    ),
    plot_data = plot.data
  )
}

.np_plot_proto_npudens_fixed_none_data <- function(bws,
                                                   xdat,
                                                   neval = 50,
                                                   xtrim = 0.0,
                                                   return.stages = FALSE) {
  .np_plot_proto_unconditional_fixed_data(
    bws = bws,
    xdat = xdat,
    neval = neval,
    xtrim = xtrim,
    plot.errors.method = "none",
    return.stages = return.stages
  )
}

.np_plot_proto_npudist_fixed_none_data <- function(bws,
                                                   xdat,
                                                   neval = 50,
                                                   xtrim = 0.0,
                                                   return.stages = FALSE) {
  .np_plot_proto_unconditional_fixed_data(
    bws = bws,
    xdat = xdat,
    neval = neval,
    xtrim = xtrim,
    plot.errors.method = "none",
    return.stages = return.stages,
    cdf = TRUE
  )
}

.np_plot_proto_rectangular_surface_base_render <- function(plot.data,
                                                          perspective = TRUE,
                                                          main = NULL,
                                                          xlab = NULL,
                                                          ylab = NULL,
                                                          zlab = NULL,
                                                          col = "lightblue",
                                                          theta = 0.0,
                                                          phi = 20.0,
                                                          ...) {
  if (!is.list(plot.data))
    stop("renderer prototype requires a plot-data list", call. = FALSE)

  obj <- NULL
  value.name <- NULL
  value.label <- NULL
  if (!is.null(plot.data$r1) && inherits(plot.data$r1, "npregression")) {
    obj <- plot.data$r1
    value.name <- "mean"
    value.label <- "Conditional mean"
  } else if (!is.null(plot.data$d1) && inherits(plot.data$d1, "npdensity")) {
    obj <- plot.data$d1
    value.name <- "dens"
    value.label <- "Density"
  } else if (!is.null(plot.data$d1) && inherits(plot.data$d1, "npdistribution")) {
    obj <- plot.data$d1
    value.name <- "dist"
    value.label <- "Distribution"
  }
  if (is.null(obj))
    stop("renderer prototype requires npregression, npdensity, or npdistribution plot data", call. = FALSE)

  xeval <- obj$eval
  if (!is.data.frame(xeval) || ncol(xeval) != 2L)
    stop("renderer prototype requires a two-column rectangular evaluation grid", call. = FALSE)
  x <- unique(as.numeric(xeval[, 1L]))
  y <- unique(as.numeric(xeval[, 2L]))
  z <- matrix(as.numeric(obj[[value.name]]), nrow = length(x), ncol = length(y), byrow = FALSE)
  if (length(x) * length(y) != length(obj[[value.name]]) ||
      any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(z))) {
    stop("renderer prototype requires a finite rectangular surface", call. = FALSE)
  }

  xlab <- if (is.null(xlab)) names(xeval)[1L] else xlab
  ylab <- if (is.null(ylab)) names(xeval)[2L] else ylab
  zlab <- if (is.null(zlab)) value.label else zlab
  main <- if (is.null(main)) value.label else main

  if (isTRUE(perspective)) {
    graphics::persp(
      x = x,
      y = y,
      z = z,
      theta = theta,
      phi = phi,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      main = main,
      col = col,
      ...
    )
  } else {
    graphics::image(
      x = x,
      y = y,
      z = z,
      xlab = xlab,
      ylab = ylab,
      main = main,
      col = grDevices::hcl.colors(64L, "YlOrRd", rev = TRUE),
      ...
    )
  }
  invisible(plot.data)
}
