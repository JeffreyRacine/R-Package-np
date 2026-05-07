.np_plot_proto_check_npcdens_lc_fixed_none <- function(bws,
                                                       xdat,
                                                       ydat,
                                                       neval,
                                                       xtrim,
                                                       ytrim) {
  if (!inherits(bws, "conbandwidth"))
    stop("prototype route requires a conditional density bandwidth object", call. = FALSE)
  if (!identical(if (is.null(bws$regtype)) "lc" else as.character(bws$regtype), "lc"))
    stop("prototype route currently supports regtype='lc' only", call. = FALSE)
  if (!identical(as.character(bws$type), "fixed"))
    stop("prototype route currently supports bwtype='fixed' only", call. = FALSE)
  if (bws$xndim != 1L || bws$yndim != 1L)
    stop("prototype route currently supports one x variable and one y variable", call. = FALSE)
  if (bws$xnuno + bws$ynuno != 0L)
    stop("prototype route currently supports continuous/ordered surface variables only", call. = FALSE)
  if (bws$xncon + bws$xnord + bws$yncon + bws$ynord != 2L)
    stop("prototype route requires a two-dimensional conditional density surface", call. = FALSE)
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) || neval < 2L)
    stop("prototype route requires scalar neval >= 2", call. = FALSE)
  invisible(TRUE)
}

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

.np_plot_proto_npcdens_lc_fixed_none_data <- function(bws,
                                                      xdat,
                                                      ydat,
                                                      neval = 50,
                                                      xtrim = 0.0,
                                                      ytrim = 0.0,
                                                      return.stages = FALSE) {
  ## Contract: first plot-runtime vertical slice. It reproduces the current
  ## npcdens LC/fixed/no-error plot.behavior='data' surface object while keeping
  ## state cleanup, target construction, evaluator, and plot-data assembly
  ## inspectable as separate stages.
  if (missing(xdat) || missing(ydat))
    stop("prototype route requires explicit xdat and ydat", call. = FALSE)
  dat <- .np_plot_proto_clean_conditional_data(xdat = xdat, ydat = ydat)
  xdat <- dat$xdat
  ydat <- dat$ydat
  .np_plot_proto_check_npcdens_lc_fixed_none(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    neval = neval,
    xtrim = xtrim,
    ytrim = ytrim
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
    cdf = FALSE,
    gradients = FALSE,
    proper = FALSE
  )

  cd1 <- condensity(
    bws = bws,
    xeval = grid$x.eval[, 1L],
    yeval = grid$x.eval[, 2L],
    ntrain = nrow(xdat),
    condens = fit$condens,
    conderr = matrix(fit$conderr, nrow = length(fit$condens), ncol = 3L)[, 1:2, drop = FALSE],
    proper.requested = fit$proper.requested,
    proper.applied = fit$proper.applied,
    proper.method = fit$proper.method,
    condens.raw = fit$condens.raw,
    proper.info = fit$proper.info
  )
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
      family = "npcdens",
      cdf = FALSE,
      gradients = FALSE
    ),
    target_grid = grid,
    evaluator = fit,
    plot_data = plot.data
  )
}
