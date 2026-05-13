## Kernel copula estimation.  The worker preserves the historical data-frame
## return columns while the public wrapper adds formula routing and S3 metadata.

npcopula <- function(bws, ...) {
  if (!missing(bws)) {
    UseMethod("npcopula", bws)
  } else {
    UseMethod("npcopula", NULL)
  }
}

.npcopula_validate_target <- function(target) {
  match.arg(target[1L], c("distribution", "density"))
}

.npcopula_validate_evaluation <- function(evaluation) {
  match.arg(evaluation[1L], c("grid", "sample"))
}

.npcopula_validate_neval <- function(neval) {
  if (!is.numeric(neval) || length(neval) != 1L || is.na(neval) ||
      !is.finite(neval) || neval != floor(neval) || neval < 2L)
    stop("'neval' must be a positive integer greater than one")
  as.integer(neval)
}

.npcopula_make_auto_u <- function(xnames, neval) {
  out <- as.data.frame(
    replicate(length(xnames), seq(0, 1, length.out = neval), simplify = FALSE)
  )
  names(out) <- xnames
  out
}

.npcopula_grid_dim <- function(u, num.var) {
  if (is.null(u))
    return(NULL)
  rep.int(nrow(as.data.frame(u)), num.var)
}

.npcopula_progress_begin <- function(target, evaluation, total) {
  state <- .np_progress_begin(
    label = sprintf("Copula %s %s", target, evaluation),
    total = total,
    domain = "general",
    surface = "copula"
  )
  state$start_note_grace_sec <- 0
  .np_progress_show_now(state, done = if (!is.null(total)) 0L else NULL)
}

.npcopula_progress_step <- function(state, done, detail) {
  if (is.null(state))
    return(state)
  .np_progress_step_at(
    state = state,
    now = .np_progress_now(),
    done = done,
    detail = detail,
    force = TRUE
  )
}

.npcopula_progress_total <- function(density, u.provided, num.var) {
  if (!u.provided)
    return(1L + num.var)
  if (density)
    return(2L * num.var + 2L)
  num.var + 2L
}

.npcopula_object <- function(result,
                             bws,
                             target,
                             evaluation,
                             u.provided,
                             u.auto,
                             grid.dim,
                             neval,
                             n.quasi.inv,
                             er.quasi.inv,
                             timing = NULL) {
  class(result) <- c("npcopula", "data.frame")
  attr(result, "bws") <- bws
  attr(result, "target") <- target
  attr(result, "density") <- identical(target, "density")
  attr(result, "xnames") <- bws$xnames
  attr(result, "ntrain") <- bws$nobs
  attr(result, "u.provided") <- isTRUE(u.provided)
  attr(result, "u.auto") <- isTRUE(u.auto)
  attr(result, "evaluation") <- evaluation
  attr(result, "grid.dim") <- grid.dim
  attr(result, "neval") <- neval
  attr(result, "n.quasi.inv") <- n.quasi.inv
  attr(result, "er.quasi.inv") <- er.quasi.inv
  attr(result, "timing") <- timing
  result
}

print.npcopula <- function(x, ...) {
  target <- attr(x, "target")
  evaluation <- attr(x, "evaluation")
  ntrain <- attr(x, "ntrain")
  xnames <- attr(x, "xnames")
  grid.dim <- attr(x, "grid.dim")

  cat("\nKernel Copula ", if (identical(target, "density")) "Density" else "Distribution",
      ": ", ntrain, " training points, in ", length(xnames),
      " variable(s)\n", sep = "")
  cat("Evaluation: ", evaluation, sep = "")
  if (!is.null(grid.dim))
    cat(" (grid ", paste(grid.dim, collapse = " x "), ")", sep = "")
  cat("\n\n")
  print(utils::head(as.data.frame(x), ...))
  invisible(x)
}

summary.npcopula <- function(object, ...) {
  print(object)
  cat("\nBandwidth summary:\n")
  summary(attr(object, "bws"))
  cat("\nCopula value summary:\n")
  print(summary(object$copula))
  invisible(object)
}

fitted.npcopula <- function(object, ...) {
  object$copula
}

plot.npcopula <- function(x,
                          perspective = TRUE,
                          view = c("surface", "contour", "image"),
                          renderer = c("base", "rgl"),
                          theta = 45,
                          phi = 30,
                          xlab = "u1",
                          ylab = "u2",
                          zlab = NULL,
                          main = NULL,
                          col = NULL,
                          border = "black",
                          zlim = NULL,
                          ...) {
  view <- match.arg(view)
  renderer <- .np_plot_match_renderer(renderer)
  dots <- list(...)
  target <- attr(x, "target")
  target.label <- if (identical(target, "density")) "Copula Density" else "Copula"
  if (is.null(zlab))
    zlab <- target.label
  if (is.null(main))
    main <- target.label

  xnames <- attr(x, "xnames")
  if (length(xnames) != 2L)
    stop("plot.npcopula currently supports two-dimensional copula displays")

  evaluation <- attr(x, "evaluation")
  if (identical(evaluation, "sample")) {
    if (identical(renderer, "rgl"))
      stop("renderer='rgl' is supported only for grid surface displays in plot.npcopula",
           call. = FALSE)
    do.call(graphics::plot, c(list(x = x$u1, y = x$u2,
                                   xlab = xlab, ylab = ylab, main = main),
                              dots))
    return(invisible(x))
  }

  grid.dim <- attr(x, "grid.dim")
  if (is.null(grid.dim) || length(grid.dim) != 2L ||
      prod(grid.dim) != nrow(x))
    stop("npcopula grid output cannot be reshaped into a two-dimensional surface")

  u1 <- sort(unique(x$u1))
  u2 <- sort(unique(x$u2))
  if (length(u1) != grid.dim[1L] || length(u2) != grid.dim[2L])
    stop("npcopula grid output is not rectangular")

  z <- matrix(x$copula, nrow = length(u1), ncol = length(u2))
  if (isTRUE(perspective) && identical(view, "surface")) {
    renderer <- .np_plot_validate_renderer_request(
      renderer = renderer,
      route = "plot.npcopula",
      perspective = perspective,
      supported.route = TRUE,
      view = "fixed"
    )
    if (identical(renderer, "rgl")) {
      rgl.view <- .np_plot_rgl_view_angles(theta = theta, phi = phi)
      rgl.out <- .np_plot_render_surface_rgl(
        x = u1,
        y = u2,
        z = z,
        xlab = xlab,
        ylab = ylab,
        zlab = zlab,
        main = main,
        theta = rgl.view$theta,
        phi = rgl.view$phi,
        col = col,
        border = border,
        zlim = zlim,
        par3d.args = .np_plot_collect_rgl_args(dots, "rgl.par3d", "rgl.par3d."),
        view3d.args = .np_plot_collect_rgl_args(dots, "rgl.view3d", "rgl.view3d."),
        persp3d.args = .np_plot_collect_rgl_args(dots, "rgl.persp3d", "rgl.persp3d."),
        grid3d.args = .np_plot_collect_rgl_args(dots, "rgl.grid3d", "rgl.grid3d."),
        widget.args = .np_plot_collect_rgl_args(dots, "rgl.widget", "rgl.widget.")
      )
      return(.np_plot_rgl_finalize(rgl.out = rgl.out))
    }
    persp.args <- .np_plot_user_args(dots, "persp")
    persp.col <- grDevices::adjustcolor(
      .np_plot_persp_surface_colors(z = z, col = col),
      alpha.f = 0.5
    )
    persp.call <- list(x = u1, y = u2, z = z,
                       theta = theta, phi = phi,
                       xlab = xlab, ylab = ylab, zlab = zlab,
                       main = main, col = persp.col, border = border)
    if (!is.null(zlim))
      persp.call$zlim <- zlim
    do.call(graphics::persp, .np_plot_merge_user_args(persp.call, persp.args))
  } else if (identical(view, "contour")) {
    if (identical(renderer, "rgl"))
      stop("renderer='rgl' is supported only for view='surface' in plot.npcopula",
           call. = FALSE)
    do.call(graphics::contour, c(list(x = u1, y = u2, z = z,
                                      xlab = xlab, ylab = ylab, main = main),
                                 dots))
  } else {
    if (identical(renderer, "rgl"))
      stop("renderer='rgl' is supported only for view='surface' in plot.npcopula",
           call. = FALSE)
    image.col <- if (is.null(col)) {
      grDevices::hcl.colors(100L, palette = "viridis")
    } else {
      col
    }
    do.call(graphics::image, c(list(x = u1, y = u2, z = z,
                                    xlab = xlab, ylab = ylab, main = main,
                                    col = image.col),
                               dots))
  }
  invisible(x)
}

npcopula.formula <- function(bws,
                             data = NULL,
                             u = NULL,
                             target = c("distribution", "density"),
                             evaluation = c("grid", "sample"),
                             neval = 30,
                             n.quasi.inv = 1000,
                             er.quasi.inv = 1,
                             ...) {
  target <- .npcopula_validate_target(target)
  evaluation <- .npcopula_validate_evaluation(evaluation)
  neval <- .npcopula_validate_neval(neval)

  mf <- stats::model.frame(bws, data = data)
  dat <- mf[, attr(attr(mf, "terms"), "term.labels"), drop = FALSE]
  num.var <- ncol(dat)

  u.auto <- FALSE
  if (is.null(u) && identical(evaluation, "grid")) {
    if (num.var != 2L)
      stop("automatic copula probability grids are supported only for two variables; supply 'u' or use evaluation='sample'")
    u <- .npcopula_make_auto_u(names(dat), neval)
    u.auto <- TRUE
  }

  bw.fun <- if (identical(target, "density")) npudensbw else npudistbw
  bw.args <- c(list(formula = bws), if (is.null(data)) list() else list(data = data), list(...))
  bw <- .np_progress_select_bandwidth_enhanced(
    if (identical(target, "density")) "Selecting copula density bandwidth" else "Selecting copula distribution bandwidth",
    do.call(bw.fun, bw.args)
  )

  npcopula.default(
    bws = bw,
    data = dat,
    u = u,
    target = target,
    evaluation = if (is.null(u)) "sample" else "grid",
    neval = neval,
    n.quasi.inv = n.quasi.inv,
    er.quasi.inv = er.quasi.inv,
    u.auto = u.auto
  )
}

npcopula.default <- function(bws,
                             data,
                             u = NULL,
                             target = NULL,
                             evaluation = c("sample", "grid"),
                             neval = 30,
                             n.quasi.inv = 1000,
                             er.quasi.inv = 1,
                             u.auto = FALSE,
                             ...) {
  evaluation <- .npcopula_validate_evaluation(evaluation)
  neval <- .npcopula_validate_neval(neval)
  start.time <- proc.time()[3]

  if(missing(data)) stop("You must provide a data frame")
  if(!is.data.frame(data)) stop("Object `data' must be a data frame")
  if(missing(bws)) stop("You must provide a bandwidth object")
  if(!isa(bws,"dbandwidth") && !isa(bws,"bandwidth"))
    stop("you must provide a density (npudensbw) or distribution (npudistbw) object")

  density <- isa(bws, "bandwidth")
  inferred.target <- if (density) "density" else "distribution"
  if (!is.null(target)) {
    target <- .npcopula_validate_target(target)
    if (!identical(target, inferred.target))
      stop("'target' conflicts with the supplied bandwidth object")
  } else {
    target <- inferred.target
  }

  result <- .npcopula_eval(
    bws = bws,
    data = data,
    u = u,
    target = target,
    evaluation = evaluation,
    neval = neval,
    n.quasi.inv = n.quasi.inv,
    er.quasi.inv = er.quasi.inv,
    u.auto = u.auto
  )
  attr(result, "timing") <- proc.time()[3] - start.time
  result
}

.npcopula_eval <- function(bws,
                           data,
                           u = NULL,
                           target,
                           evaluation,
                           neval,
                           n.quasi.inv,
                           er.quasi.inv,
                           u.auto = FALSE) {
  density <- identical(target, "density")
  if(!is.null(u)) {
    if(is.vector(u) && !is.list(u)) u <- matrix(u,1,length(u))
    if(anyNA(u)) stop("u must not contain missing values")
    if(any(u>1 | u<0, na.rm = TRUE)) stop("u must lie in [0,1]")
  }
  num.var <- length(bws$xnames)
  if(!is.null(u) && (ncol(u)!=num.var)) stop("u and bws are incompatible")
  if(n.quasi.inv < 1) stop("n.quasi.inv must be a positive integer")
  if(any(is.na(data))) stop("NA values present in data - recompute bw object with na.omit on data")
  if(bws$nuno>0) stop("unordered factors not suitable for copula estimation")

  u.provided <- !is.null(u)
  grid.dim <- .npcopula_grid_dim(u, num.var)

  if(!is.null(u)) {
    if(ncol(u) != num.var) stop(paste("matrix u must have ", num.var," columns",sep=""))
  }

  progress.total <- .npcopula_progress_total(density = density,
                                             u.provided = u.provided,
                                             num.var = num.var)
  progress <- .npcopula_progress_begin(target = target,
                                       evaluation = if (u.provided) "grid" else "sample",
                                       total = progress.total)
  on.exit(.np_progress_end(progress), add = TRUE)
  stage <- 0L

  if(is.null(u)) {
    stage <- stage + 1L
    progress <- .npcopula_progress_step(
      progress, stage,
      if(!density) "joint distribution at sample realizations" else "joint density at sample realizations"
    )
    if(!density) {
      copula <- fitted(npudist(bws=bws,data=data))
    } else {
      copula <- fitted(npudens(bws=bws,data=data))
    }

    u <- matrix(NA,bws$nobs,num.var)
    for (j in seq_len(num.var)) {
      stage <- stage + 1L
      progress <- .npcopula_progress_step(
        progress, stage,
        sprintf("marginal %s at sample realizations", bws$xnames[j])
      )
      bws.F.marginal <- npudistbw(formula(paste("~",bws$xnames[j])),
                         bws=bws$bw[j],
                         bandwidth.compute=FALSE,
                         bwtype=bws$type,
                         ckerorder=bws$ckerorder,
                         ckertype=bws$ckertype,
                         okertype=bws$okertype,
                         data=data)

      u[,j] <- fitted(npudist(bws=bws.F.marginal,data=data))
      if(density) {
        bws.f.marginal <- npudensbw(formula(paste("~",bws$xnames[j])),
                           bws=bws$bw[j],
                           bandwidth.compute=FALSE,
                           bwtype=bws$type,
                           ckerorder=bws$ckerorder,
                           ckertype=bws$ckertype,
                           okertype=bws$okertype,
                           ukertype=bws$ukertype,
                           data=data)
        copula <- copula/NZD(fitted(npudens(bws=bws.f.marginal,data=data)))

      }
    }
  } else {
    n.u <- nrow(u)
    x.u <- data.frame(matrix(NA,n.u,num.var))
    names(x.u) <- bws$xnames
    for (j in seq_len(num.var)) {
      stage <- stage + 1L
      progress <- .npcopula_progress_step(
        progress, stage,
        sprintf("quasi-inverse marginal %s", bws$xnames[j])
      )
      x.marginal <- data[[bws$xnames[j]]]
      quantile.seq <- seq(0,1,length=round(n.quasi.inv/2))
      if(is.numeric(x.marginal)) {
        x.er <- extendrange(x.marginal,f=er.quasi.inv)
        x.q <- quantile(x.marginal,quantile.seq)
        x.eval <- sort(c(seq(x.er[1],x.er[2],length=round(n.quasi.inv/2)),x.q))
      } else {
        x.u[,j] <- ordered(x.u[,j],levels=levels(x.marginal))
        x.q <- sapply(seq_len(round(n.quasi.inv/2)), function(i) { uocquantile(x.marginal, quantile.seq[i]) })
        x.eval <- sort(ordered(c(as.character(x.q),as.character(x.q)),levels=levels(x.marginal)))
      }
      F <- fitted(npudist(tdat=x.marginal,
                          edat=x.eval,
                          bws=bws$bw[j],
                          bwtype=bws$type,
                          ckerorder=bws$ckerorder,
                          ckertype=bws$ckertype,
                          okertype=bws$okertype,
                          ukertype=bws$ukertype,data=data))
      for (i in seq_len(n.u)) {
        u[u[,j]<min(F),j] <- min(F)
        u[u[,j]>max(F),j] <- max(F)
        x.u[i,j] <-  min(x.eval[F>=u[i,j]])
      }
    }

    stage <- stage + 1L
    progress <- .npcopula_progress_step(progress, stage, "expand u grid")
    x.u <- expand.grid(data.frame(x.u))
    for (k in seq_len(ncol(x.u))) {
      if(is.ordered(data[,k])) x.u[,k] <- ordered(x.u[,k],levels=levels(data[,k]))
    }
    unit.weights <- rep_len(1, nrow(data))
    if(!density) {
      stage <- stage + 1L
      progress <- .npcopula_progress_step(progress, stage, "joint distribution on expanded grid")
      copula <- npudisthat(
        bws = bws,
        tdat = data,
        edat = x.u,
        y = unit.weights,
        output = "apply"
      )
    } else {
      stage <- stage + 1L
      progress <- .npcopula_progress_step(progress, stage, "joint density on expanded grid")
      copula <- npudenshat(
        bws = bws,
        tdat = data,
        edat = x.u,
        y = unit.weights,
        output = "apply"
      )
      for (j in seq_len(num.var)) {
        stage <- stage + 1L
        progress <- .npcopula_progress_step(
          progress, stage,
          sprintf("density marginal %s on expanded grid", bws$xnames[j])
        )
        bws.f.marginal <- npudensbw(formula(paste("~",bws$xnames[j])),
                           bws=bws$bw[j],
                           bandwidth.compute=FALSE,
                           bwtype=bws$type,
                           ckerorder=bws$ckerorder,
                           ckertype=bws$ckertype,
                           okertype=bws$okertype,
                           ukertype=bws$ukertype,
                           data=data)
        xeval <- data.frame(x.u[,j])
        names(xeval) <- bws$xnames[j]
        copula <- copula/NZD(npudenshat(
          bws = bws.f.marginal,
          tdat = data[, bws$xnames[j], drop = FALSE],
          edat = xeval,
          y = unit.weights,
          output = "apply"
        ))
      }
    }
  }

  if(!u.provided) {
    u <- data.frame(u)
    names(u) <- paste("u", seq_len(num.var), sep = "")
    out <- data.frame(copula,u)
  } else {
    u <- expand.grid(data.frame(u))
    names(u) <- paste("u", seq_len(num.var), sep = "")
    out <- data.frame(copula,u,x.u)
  }

  .npcopula_object(
    result = out,
    bws = bws,
    target = target,
    evaluation = if (u.provided) "grid" else "sample",
    u.provided = u.provided,
    u.auto = u.auto,
    grid.dim = grid.dim,
    neval = neval,
    n.quasi.inv = n.quasi.inv,
    er.quasi.inv = er.quasi.inv
  )
}
