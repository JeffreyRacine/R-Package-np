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
  if (!is.list(plot.data) || is.null(plot.data$cd1) || !inherits(plot.data$cd1, "condensity"))
    stop("renderer prototype requires plot-data with a condensity 'cd1' element", call. = FALSE)
  cd <- plot.data$cd1
  x <- unique(as.numeric(cd$xeval))
  y <- unique(as.numeric(cd$yeval))
  z <- matrix(as.numeric(cd$condens), nrow = length(x), ncol = length(y), byrow = FALSE)
  if (length(x) * length(y) != length(cd$condens) ||
      any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(z))) {
    stop("renderer prototype requires a finite rectangular conditional-density surface", call. = FALSE)
  }

  xlab <- if (is.null(xlab)) cd$xnames[1L] else xlab
  ylab <- if (is.null(ylab)) cd$ynames[1L] else ylab
  zlab <- if (is.null(zlab)) "Conditional density" else zlab
  main <- if (is.null(main)) "Conditional Density" else main

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
