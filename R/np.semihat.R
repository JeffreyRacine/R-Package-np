.np_semihat_require_class <- function(bws, class.name, fn.name) {
  if (!inherits(bws, class.name))
    stop(sprintf("argument 'bws' must inherit from class '%s' in %s()", class.name, fn.name))
}

.npRmpi_with_local_hat_helper <- function(expr) {
  old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  options(npRmpi.autodispatch.context = TRUE)
  options(npRmpi.autodispatch.disable = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  force(expr)
}

.np_semihat_make_regbw_args <- function(source, xdat, ydat, bw) {
  xdat <- toFrame(xdat)
  ncon.target <- sum(untangle(xdat)$icon)
  spec <- npCanonicalConditionalRegSpec(
    regtype = if (is.null(source$regtype)) "lc" else as.character(source$regtype),
    basis = if (is.null(source$basis)) "glp" else as.character(source$basis),
    degree = source$degree,
    bernstein.basis = isTRUE(source$bernstein.basis),
    ncon = ncon.target,
    where = "npscoef"
  )

  collapse_bound <- function(v, nm) {
    if (is.null(v))
      return(NULL)
    vv <- as.double(v)
    if (ncon.target <= 0L)
      return(NULL)
    if (length(vv) <= 1L || length(vv) == ncon.target)
      return(vv)
    if (ncon.target == 1L) {
      uu <- unique(vv)
      if (length(uu) == 1L)
        return(uu)
      stop(sprintf("cannot collapse %s with %d distinct values to scalar helper bound", nm, length(uu)),
           call. = FALSE)
    }
    stop(sprintf("%s length (%d) is incompatible with helper continuous dimension (%d)",
                 nm, length(vv), ncon.target),
         call. = FALSE)
  }

  args <- list(
    xdat = xdat,
    ydat = ydat,
    bws = as.double(bw),
    regtype = spec$regtype.engine,
    bandwidth.compute = FALSE
  )

  # Forward the full compatible estimator/bandwidth option contract.
  args$basis <- spec$basis.engine
  args$degree <- spec$degree.engine
  args$bernstein.basis <- spec$bernstein.basis.engine
  if (!is.null(source$method) && source$method %in% c("cv.ls", "cv.aic"))
    args$bwmethod <- source$method
  if (!is.null(source$scaling))
    args$bwscaling <- isTRUE(source$scaling)
  if (!is.null(source$type))
    args$bwtype <- source$type
  if (!is.null(source$ckertype))
    args$ckertype <- source$ckertype
  if (!is.null(source$ckerorder))
    args$ckerorder <- source$ckerorder
  if (!is.null(source$ckerbound))
    args$ckerbound <- source$ckerbound
  if (!is.null(source$ckerlb))
    args$ckerlb <- collapse_bound(source$ckerlb, "ckerlb")
  if (!is.null(source$ckerub))
    args$ckerub <- collapse_bound(source$ckerub, "ckerub")
  if (!is.null(source$ukertype))
    args$ukertype <- source$ukertype
  if (!is.null(source$okertype))
    args$okertype <- source$okertype

  args
}

.np_semihat_make_regbw_state <- function(source, xdat, ydat, bw) {
  resolve_choice <- function(value, choices, default = choices[1L]) {
    if (is.null(value) || !length(value))
      return(default)
    cand <- as.character(value)[1L]
    if (!cand %in% choices)
      return(default)
    cand
  }

  xdat <- toFrame(xdat)
  if (!(is.vector(ydat) || is.factor(ydat)))
    stop("'ydat' must be a vector")

  args <- .np_semihat_make_regbw_args(source = source, xdat = xdat, ydat = ydat, bw = bw)
  xdati <- untangle(xdat)
  ydati <- untangle(data.frame(ydat))
  xmat <- toMatrix(xdat)
  rcon <- xmat[, xdati$icon, drop = FALSE]
  nobs <- nrow(xdat)
  ncon <- sum(xdati$icon)
  method <- resolve_choice(if (!is.null(args$bwmethod)) args$bwmethod else source$method,
                           c("cv.ls", "cv.aic"))
  bwtype <- resolve_choice(if (!is.null(args$bwtype)) args$bwtype else source$type,
                           c("fixed", "generalized_nn", "adaptive_nn"))
  ckertype <- resolve_choice(if (!is.null(args$ckertype)) args$ckertype else source$ckertype,
                             c("gaussian", "truncated gaussian", "epanechnikov", "uniform"))
  ckerorder <- if (!is.null(args$ckerorder)) args$ckerorder else source$ckerorder
  ckerbound <- resolve_choice(if (!is.null(args$ckerbound)) args$ckerbound else source$ckerbound,
                              c("none", "range", "fixed"))
  ukertype <- resolve_choice(if (!is.null(args$ukertype)) args$ukertype else source$ukertype,
                             c("aitchisonaitken", "liracine"))
  okertype <- resolve_choice(if (!is.null(args$okertype)) args$okertype else source$okertype,
                             c("liracine", "wangvanryzin", "racineliyan", "nliracine"))
  scaling <- if (!is.null(args$bwscaling)) isTRUE(args$bwscaling) else FALSE
  bw.vec <- stats::setNames(as.double(args$bws), names(xdat))
  sfactor <- bandwidth <- bw.vec
  nconfac <- nobs^(-1.0 / (2.0 * ckerorder + ncon))
  ncatfac <- nobs^(-2.0 / (2.0 * ckerorder + ncon))
  sdev <- if (ncon > 0L) EssDee(rcon) else numeric(0)

  if (sum(xdati$iuno) > 0L) {
    if (scaling) {
      bandwidth[xdati$iuno] <- bandwidth[xdati$iuno] * ncatfac
    } else {
      sfactor[xdati$iuno] <- sfactor[xdati$iuno] / ncatfac
    }
  }

  if (sum(xdati$iord) > 0L) {
    if (scaling) {
      bandwidth[xdati$iord] <- bandwidth[xdati$iord] * ncatfac
    } else {
      sfactor[xdati$iord] <- sfactor[xdati$iord] / ncatfac
    }
  }

  if (ncon > 0L) {
    dfactor <- sdev * nconfac
    if (scaling) {
      bandwidth[xdati$icon] <- bandwidth[xdati$icon] * dfactor
    } else {
      sfactor[xdati$icon] <- sfactor[xdati$icon] / dfactor
    }
  }

  out <- rbandwidth(
    bw = bw.vec,
    regtype = args$regtype,
    basis = args$basis,
    degree = args$degree,
    bernstein.basis = args$bernstein.basis,
    bwmethod = method,
    bwscaling = scaling,
    bwtype = bwtype,
    ckertype = ckertype,
    ckerorder = ckerorder,
    ckerbound = ckerbound,
    ckerlb = args$ckerlb,
    ckerub = args$ckerub,
    ukertype = ukertype,
    okertype = okertype,
    fval = NA,
    ifval = NA,
    num.feval = NA,
    num.feval.fast = NA,
    fval.history = NA,
    eval.history = NA,
    invalid.history = NA,
    nobs = nobs,
    xdati = xdati,
    ydati = ydati,
    xnames = names(xdat),
    ynames = if (!is.null(source$ynames) && length(source$ynames)) source$ynames else "y",
    sfactor = sfactor,
    bandwidth = bandwidth,
    rows.omit = NA,
    nconfac = nconfac,
    ncatfac = ncatfac,
    sdev = sdev,
    bandwidth.compute = FALSE,
    timing = NA,
    total.time = c(elapsed = 0.0)
  )
  out$call <- call(".np_semihat_make_regbw_state")
  out
}

.npscoef_canonical_spec <- function(source, zdat, where = "npscoef") {
  zdat <- toFrame(zdat)
  npCanonicalConditionalRegSpec(
    regtype = if (is.null(source$regtype)) "lc" else as.character(source$regtype),
    basis = if (is.null(source$basis)) "glp" else as.character(source$basis),
    degree = source$degree,
    bernstein.basis = isTRUE(source$bernstein.basis),
    ncon = sum(untangle(zdat)$icon),
    where = where
  )
}

.np_indexhat_rbw <- function(bws, idx.train) {
  .np_semihat_make_regbw_state(
    source = bws,
    xdat = idx.train,
    ydat = rep.int(0.0, nrow(idx.train)),
    bw = bws$bw
  )
}

.np_indexhat_ll_owner_rbw <- function(bws, idx.train) {
  base <- .np_indexhat_rbw(bws = bws, idx.train = idx.train)
  npregbw(
    xdat = idx.train,
    ydat = rep.int(0.0, nrow(idx.train)),
    bws = base$bw,
    regtype = "ll",
    bwtype = base$type,
    bandwidth.compute = FALSE,
    ckertype = base$ckertype,
    ckerorder = base$ckerorder,
    ckerbound = base$ckerbound,
    ckerlb = base$ckerlb,
    ckerub = base$ckerub,
    ukertype = base$ukertype,
    okertype = base$okertype
  )
}

.np_indexhat_numeric_y <- function(y) {
  if (is.factor(y) || is.vector(y))
    return(matrix(as.double(y), ncol = 1L))

  as.matrix(y)
}

.np_indexhat_kbw <- function(bws, idx.train) {
  resolve_choice <- function(value, choices, default = choices[1L]) {
    if (is.null(value) || !length(value))
      return(default)
    cand <- as.character(value)[1L]
    if (!cand %in% choices)
      return(default)
    cand
  }

  collapse_bound <- function(v, nm) {
    if (is.null(v))
      return(NULL)
    vv <- as.double(v)
    if (length(vv) <= 1L)
      return(vv)
    uu <- unique(vv)
    if (length(uu) == 1L)
      return(uu)
    stop(sprintf("cannot collapse %s with %d distinct values to scalar helper bound",
                 nm, length(uu)),
         call. = FALSE)
  }

  bwtype <- resolve_choice(bws$type, c("fixed", "generalized_nn", "adaptive_nn"))
  ckertype <- resolve_choice(bws$ckertype, c("gaussian", "truncated gaussian", "epanechnikov", "uniform"))
  ckerbound <- resolve_choice(bws$ckerbound, c("none", "range", "fixed"))
  ukertype <- resolve_choice(bws$ukertype, c("aitchisonaitken", "liracine"))
  okertype <- resolve_choice(bws$okertype, c("liracine", "wangvanryzin", "racineliyan", "nliracine"))

  kbandwidth.numeric(
    bw = c(bws$bw),
    bwtype = bwtype,
    bwscaling = FALSE,
    ckertype = ckertype,
    ckerorder = bws$ckerorder,
    ckerbound = ckerbound,
    ckerlb = collapse_bound(bws$ckerlb, "ckerlb"),
    ckerub = collapse_bound(bws$ckerub, "ckerub"),
    ukertype = ukertype,
    okertype = okertype,
    nobs = nrow(idx.train),
    xdati = untangle(idx.train),
    ydati = bws$ydati,
    xnames = "index",
    ynames = bws$ynames,
    xdat = idx.train
  )
}

.np_indexhat_core <- function(bws,
                              idx.train,
                              idx.eval,
                              y = NULL,
                              output = c("matrix", "apply"),
                              ridge = 0.0) {
  output <- match.arg(output)
  kbw <- .np_indexhat_kbw(bws = bws, idx.train = idx.train)
  spec <- .npindex_resolve_spec(bws, where = "npindexhat")
  regtype <- spec$regtype.engine

  kw <- .np_kernel_weights_direct(
    bws = kbw,
    txdat = idx.train,
    exdat = idx.eval,
    bandwidth.divide = TRUE,
    kernel.pow = 1.0
  )

  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(idx.train))

  ntrain <- nrow(idx.train)
  neval <- nrow(idx.eval)
  if (nrow(kw) != ntrain || ncol(kw) != neval)
    stop("single-index hat kernel-weight matrix shape mismatch", call. = FALSE)

  if (!is.null(y)) {
    y <- .np_indexhat_numeric_y(y)
    if (nrow(y) != ntrain)
      stop("number of rows in 'y' must equal number of training rows")
  }

  if (identical(regtype, "lc")) {
    den <- pmax(colSums(kw), .Machine$double.eps)
    if (identical(output, "matrix"))
      return(sweep(t(kw), 1L, den, "/", check.margin = FALSE))

    out <- sweep(t(kw), 1L, den, "/", check.margin = FALSE) %*% y
    return(if (ncol(out) == 1L) as.vector(out) else out)
  }

  degree <- if (identical(regtype, "ll")) {
    1L
  } else {
    spec$degree.engine
  }

  W <- W.lp(
    xdat = idx.train,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  W.eval <- W.lp(
    xdat = idx.train,
    exdat = idx.eval,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )

  if (!is.matrix(W))
    W <- matrix(W, nrow = ntrain)
  if (!is.matrix(W.eval))
    W.eval <- matrix(W.eval, nrow = neval)
  if (nrow(W.eval) != neval)
    W.eval <- matrix(W.eval, nrow = neval, byrow = FALSE)

  if (identical(output, "matrix")) {
    H <- matrix(NA_real_, nrow = neval, ncol = ntrain)
  } else {
    out <- matrix(0.0, nrow = neval, ncol = ncol(y))
  }

  for (i in seq_len(neval)) {
    solve.out <- .npreghat_solve_eval(
      W = W,
      w.eval = W.eval[i, ],
      k = kw[, i],
      ridge.base = ridge
    )

    if (is.null(solve.out))
      stop(sprintf("failed to solve single-index hat system at evaluation row %d", i),
           call. = FALSE)

    h.row <- kw[, i] * drop(W %*% solve.out$v)

    if (identical(output, "matrix")) {
      H[i, ] <- h.row
    } else {
      out[i, ] <- drop(crossprod(h.row, y))
    }
  }

  if (identical(output, "matrix"))
    H
  else if (ncol(out) == 1L)
    as.vector(out)
  else
    out
}

.np_indexhat_lc_kernel_weights <- function(bws, idx.train, idx.eval) {
  kw <- .npRmpi_with_local_regression(npksum(
    txdat = idx.train,
    exdat = idx.eval,
    bws = bws$bw,
    bwtype = bws$type,
    ckertype = bws$ckertype,
    ckerorder = bws$ckerorder,
    return.kernel.weights = TRUE
  ))$kw

  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(idx.train))

  if (nrow(kw) != nrow(idx.train) || ncol(kw) != nrow(idx.eval))
    stop("single-index lc kernel-weight matrix shape mismatch", call. = FALSE)

  kw
}

.np_indexhat_lc_mean <- function(bws,
                                 idx.train,
                                 idx.eval,
                                 y = NULL,
                                 output = c("matrix", "apply")) {
  output <- match.arg(output)
  kw <- .np_indexhat_lc_kernel_weights(
    bws = bws,
    idx.train = idx.train,
    idx.eval = idx.eval
  )
  H <- sweep(
    t(kw),
    1L,
    pmax(colSums(kw), .Machine$double.eps),
    "/",
    check.margin = FALSE
  )

  if (identical(output, "matrix"))
    return(H)

  if (is.null(y))
    stop("argument 'y' is required when output='apply'")

  y <- .np_indexhat_numeric_y(y)
  if (nrow(y) != nrow(idx.train))
    stop("number of rows in 'y' must equal number of training rows")

  out <- H %*% y
  if (ncol(out) == 1L) as.vector(out) else out
}

.np_indexhat_lc_derivative <- function(bws,
                                       idx.train,
                                       idx.eval,
                                       y = NULL,
                                       output = c("matrix", "apply")) {
  output <- match.arg(output)
  out <- .npRmpi_with_local_hat_helper(npreghat(
    bws = .np_indexhat_rbw(bws = bws, idx.train = idx.train),
    txdat = idx.train,
    exdat = idx.eval,
    y = y,
    output = output,
    s = 1L
  ))

  if (!identical(output, "matrix"))
    return(out)

  matrix(
    as.double(out),
    nrow = nrow(out),
    ncol = ncol(out),
    dimnames = dimnames(out)
  )
}

.np_indexhat_lp_mean_matrix <- function(bws, idx.train, idx.eval) {
  out <- .npRmpi_with_local_hat_helper(npreghat(
    bws = .np_indexhat_rbw(bws = bws, idx.train = idx.train),
    txdat = idx.train,
    exdat = idx.eval,
    output = "matrix",
    s = 0L
  ))

  matrix(
    as.double(out),
    nrow = nrow(out),
    ncol = ncol(out),
    dimnames = dimnames(out)
  )
}

.np_indexhat_gradient_matrix <- function(bws, idx.train, idx.eval) {
  regtype <- if (!is.null(bws$regtype)) as.character(bws$regtype)[1L] else "lc"
  owner.bw <- if (identical(regtype, "ll")) {
    .np_indexhat_ll_owner_rbw(bws = bws, idx.train = idx.train)
  } else {
    .np_indexhat_rbw(bws = bws, idx.train = idx.train)
  }
  out <- .npRmpi_with_local_hat_helper(npreghat(
    bws = owner.bw,
    txdat = idx.train,
    exdat = idx.eval,
    output = "matrix",
    s = 1L
  ))

  matrix(
    as.double(out),
    nrow = nrow(out),
    ncol = ncol(out),
    dimnames = dimnames(out)
  )
}

.np_indexhat_exact <- function(bws,
                               idx.train,
                               idx.eval,
                               y = NULL,
                               output = c("matrix", "apply"),
                               s = 0L) {
  output <- match.arg(output)
  spec <- .npindex_resolve_spec(bws, where = "npindexhat")
  regtype <- if (!is.null(bws$regtype)) as.character(bws$regtype)[1L] else spec$regtype.engine
  regtype.engine <- spec$regtype.engine
  if (identical(regtype, "ll")) {
    out <- .npRmpi_with_local_hat_helper(npreghat(
      bws = .np_indexhat_ll_owner_rbw(bws = bws, idx.train = idx.train),
      txdat = idx.train,
      exdat = idx.eval,
      y = y,
      output = output,
      s = s
    ))

    if (!identical(output, "matrix"))
      return(out)

    return(matrix(
      as.double(out),
      nrow = nrow(out),
      ncol = ncol(out),
      dimnames = dimnames(out)
    ))
  }

  if (identical(regtype.engine, "lc")) {
    if (s == 1L) {
      return(.np_indexhat_lc_derivative(
        bws = bws,
        idx.train = idx.train,
        idx.eval = idx.eval,
        y = y,
        output = output
      ))
    }

    return(.np_indexhat_lc_mean(
      bws = bws,
      idx.train = idx.train,
      idx.eval = idx.eval,
      y = y,
      output = output
    ))
  }

  rbw <- .np_indexhat_rbw(bws = bws, idx.train = idx.train)

  lp.mean.owner.safe <- identical(regtype.engine, "lp") &&
    identical(output, "matrix") &&
    s == 0L &&
    !identical(bws$type, "fixed") &&
    !(
      identical(bws$type, "generalized_nn") &&
        all(as.integer(spec$degree.engine) == 1L) &&
        (
          !identical(spec$basis.engine, "glp") ||
            isTRUE(spec$bernstein.basis.engine)
        )
    )

  if (lp.mean.owner.safe) {
    return(.np_indexhat_lp_mean_matrix(
      bws = bws,
      idx.train = idx.train,
      idx.eval = idx.eval
    ))
  }

  lp.grad.owner.safe <- identical(output, "matrix") &&
    s == 1L &&
    (
      identical(regtype, "ll") ||
        (
          identical(regtype.engine, "lp") &&
            !(
              identical(bws$type, "generalized_nn") &&
                all(as.integer(spec$degree.engine) == 1L) &&
                (
                  !identical(spec$basis.engine, "glp") ||
                    isTRUE(spec$bernstein.basis.engine)
                )
            )
        )
    )

  if (lp.grad.owner.safe) {
    return(.np_indexhat_gradient_matrix(
      bws = bws,
      idx.train = idx.train,
      idx.eval = idx.eval
    ))
  }

  if (s == 1L) {
    fit_one <- function(ycol) {
      fit <- .npRmpi_with_local_regression(.np_regression_direct(
        bws = rbw,
        txdat = idx.train,
        tydat = ycol,
        exdat = idx.eval,
        gradients = TRUE,
        gradient.order = 1L
      ))
      fit$grad[, 1L]
    }

    if (identical(output, "apply")) {
      if (is.null(y))
        stop("argument 'y' is required when output='apply'")

      y <- .np_indexhat_numeric_y(y)
      if (nrow(y) != nrow(idx.train))
        stop("number of rows in 'y' must equal number of training rows")

      out <- vapply(
        seq_len(ncol(y)),
        function(j) fit_one(y[, j]),
        numeric(nrow(idx.eval))
      )

      if (is.null(dim(out)))
        return(as.vector(out))
      if (ncol(out) == 1L)
        return(as.vector(out[, 1L]))
      return(out)
    }

    neval <- nrow(idx.eval)
    ntrain <- nrow(idx.train)
    H <- matrix(NA_real_, nrow = neval, ncol = ntrain)
    for (j in seq_len(ntrain)) {
      yj <- numeric(ntrain)
      yj[j] <- 1.0
      H[, j] <- fit_one(yj)
    }
    return(H)
  }

  fit_one <- function(ycol) {
    fit <- .npRmpi_with_local_regression(.np_regression_direct(
      bws = rbw,
      txdat = idx.train,
      tydat = ycol,
      exdat = idx.eval,
      gradients = FALSE
    ))
    fit$mean
  }

  if (identical(output, "matrix")) {
    neval <- nrow(idx.eval)
    ntrain <- nrow(idx.train)
    H <- matrix(NA_real_, nrow = neval, ncol = ntrain)
    for (j in seq_len(ntrain)) {
      yj <- numeric(ntrain)
      yj[j] <- 1.0
      H[, j] <- fit_one(yj)
    }
    return(H)
  }

  if (is.null(y))
    stop("argument 'y' is required when output='apply'")

  y <- .np_indexhat_numeric_y(y)
  if (nrow(y) != nrow(idx.train))
    stop("number of rows in 'y' must equal number of training rows")

  out <- matrix(0.0, nrow = nrow(idx.eval), ncol = ncol(y))
  for (j in seq_len(ncol(y)))
    out[, j] <- fit_one(y[, j])

  if (ncol(out) == 1L) as.vector(out) else out
}

.npscoef_make_regbw <- function(bws, zdat, bw = bws$bw) {
  .np_semihat_make_regbw_state(
    source = bws,
    xdat = zdat,
    ydat = rep.int(0.0, nrow(zdat)),
    bw = bw
  )
}

.npscoef_lp_state <- function(bws, tzdat, ezdat, leave.one.out = FALSE, where = "npscoef") {
  tzdat <- toFrame(tzdat)
  ezdat <- toFrame(ezdat)
  leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")
  spec <- .npscoef_canonical_spec(source = bws, zdat = tzdat, where = where)
  same.eval <- isTRUE(all.equal(tzdat, ezdat, check.attributes = FALSE))

  if (leave.one.out && !same.eval) {
    stop("leave.one.out=TRUE requires evaluation 'z' data to match training 'z' data")
  }

  state <- list(
    spec = spec,
    leave.one.out = leave.one.out
  )

  if (identical(spec$regtype.engine, "lc")) {
    state$z.train <- adjustLevels(tzdat, bws$zdati)
    state$z.eval <- if (leave.one.out) state$z.train else adjustLevels(ezdat, bws$zdati, allowNewCells = TRUE)
    return(state)
  }

  rbw <- .npscoef_make_regbw(bws = bws, zdat = tzdat)
  z.train <- adjustLevels(tzdat, rbw$xdati)
  z.eval <- if (leave.one.out) z.train else adjustLevels(ezdat, rbw$xdati, allowNewCells = TRUE)
  W.train <- W.lp(
    xdat = z.train,
    degree = spec$degree.engine,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  W.eval <- W.lp(
    xdat = z.train,
    exdat = if (leave.one.out) NULL else z.eval,
    degree = spec$degree.engine,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  neval <- nrow(z.eval)
  if (!is.matrix(W.eval))
    W.eval <- matrix(W.eval, nrow = neval)
  if (nrow(W.eval) != neval)
    W.eval <- matrix(W.eval, nrow = neval, byrow = FALSE)

  state$rbw <- rbw
  state$z.train <- z.train
  state$z.eval <- z.eval
  state$W.train <- W.train
  state$W.eval <- W.eval
  state
}

.npscoef_row_tensor_design <- function(X, Z) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  if (nrow(X) != nrow(Z))
    stop("tensor design inputs must have the same number of rows", call. = FALSE)
  do.call(cbind, lapply(seq_len(ncol(X)), function(j) Z * X[, j]))
}

.npscoef_lp1_largeh_global_fit <- function(bws,
                                           tzdat,
                                           ezdat,
                                           W.train,
                                           tydat,
                                           u2 = NULL,
                                           leave.one.out = FALSE,
                                           where = "npscoef",
                                           solver,
                                           ksum_fun = npksum) {
  if (missing(solver) || !is.function(solver))
    stop("LP1 large-h global fit requires a solver function", call. = FALSE)
  if (!is.function(ksum_fun))
    stop("LP1 large-h global fit requires a ksum function", call. = FALSE)

  lp_state <- .npscoef_lp_state(
    bws = bws,
    tzdat = tzdat,
    ezdat = ezdat,
    leave.one.out = leave.one.out,
    where = where
  )
  spec <- lp_state$spec
  if (!identical(spec$regtype.engine, "lp") ||
      !identical(spec$basis.engine, "glp") ||
      isTRUE(spec$bernstein.basis.engine) ||
      any(as.integer(spec$degree.engine) != 1L)) {
    stop("LP1 large-h global fit called outside canonical degree-1 raw glp regime",
         call. = FALSE)
  }

  tensor.train <- .npscoef_row_tensor_design(W.train, lp_state$W.train)
  z.eval.one <- lp_state$z.eval[1L, , drop = FALSE]

  ytensor <- cbind(tydat, tensor.train)
  main.ks <- ksum_fun(
    txdat = lp_state$z.train,
    tydat = ytensor,
    weights = ytensor,
    exdat = z.eval.one,
    bws = lp_state$rbw,
    bandwidth.divide = TRUE
  )$ksum
  tyw <- as.double(main.ks[-1L, 1L, 1L])
  tww <- main.ks[-1L, -1L, 1L, drop = TRUE]

  solve.out <- solver(tyw, tww)
  theta <- as.double(solve.out$coef)
  theta.mat <- matrix(theta,
                      nrow = ncol(lp_state$W.eval),
                      ncol = ncol(W.train))
  coef <- t(lp_state$W.eval %*% theta.mat)

  out <- list(
    lp_state = lp_state,
    coef = coef,
    theta = theta,
    ridge = solve.out$ridge,
    tww = tww
  )

  if (!is.null(u2)) {
    out$s <- ksum_fun(
      txdat = lp_state$z.train,
      tydat = tensor.train,
      weights = tensor.train * as.double(u2),
      exdat = z.eval.one,
      bws = lp_state$rbw,
      bandwidth.divide = TRUE,
      kernel.pow = 2
    )$ksum[, , 1L, drop = TRUE]
  } else {
    out$s <- NULL
  }

  out
}

.npscoef_effective_weight_state <- function(bws, tzdat, ezdat, leave.one.out = FALSE) {
  state <- .npscoef_lp_state(
    bws = bws,
    tzdat = tzdat,
    ezdat = ezdat,
    leave.one.out = leave.one.out,
    where = "npscoefhat"
  )
  state$bws <- bws
  state$regtype <- state$spec$regtype.engine
  state
}

.npscoef_effective_weight_one <- function(state, eval.index) {
  .npscoef_effective_weight_chunk(state = state, eval.indices = eval.index)[, 1L]
}

.npscoef_effective_weight_chunk_size <- function(ntrain, neval) {
  chunk.opt <- getOption("np.scoef.weight.chunk.size")
  if (!is.null(chunk.opt)) {
    chunk.opt <- as.integer(chunk.opt)[1L]
    if (is.na(chunk.opt) || chunk.opt < 1L)
      stop("option 'np.scoef.weight.chunk.size' must be a positive integer")
    return(min(as.integer(neval), chunk.opt))
  }

  target.bytes <- 32L * 1024L * 1024L
  chunk <- as.integer(floor(target.bytes / (8 * max(1L, as.integer(ntrain)))))
  if (!is.finite(chunk) || is.na(chunk) || chunk < 1L)
    chunk <- 1L
  min(as.integer(neval), max(1L, chunk))
}

.npscoef_effective_weight_chunk <- function(state, eval.indices) {
  eval.indices <- as.integer(eval.indices)
  if (!length(eval.indices))
    return(matrix(0.0, nrow = nrow(state$z.train), ncol = 0L))
  if (anyNA(eval.indices) || any(!is.finite(eval.indices)) ||
      any(eval.indices < 1L) || any(eval.indices > nrow(state$z.eval))) {
    stop("invalid evaluation index for smooth-coefficient effective weights", call. = FALSE)
  }

  eval.rows <- state$z.eval[eval.indices, , drop = FALSE]
  bw.use <- if (identical(state$regtype, "lc")) state$bws else state$rbw
  kw <- .np_kernel_weights_direct(
    txdat = state$z.train,
    exdat = eval.rows,
    bws = bw.use,
    bandwidth.divide = TRUE,
    kernel.pow = 1.0
  )

  if (isTRUE(state$leave.one.out)) {
    for (jj in seq_along(eval.indices))
      kw[eval.indices[[jj]], jj] <- 0.0
  }

  if (identical(state$regtype, "lc"))
    return(kw)

  out <- matrix(0.0, nrow = nrow(kw), ncol = ncol(kw))
  for (jj in seq_along(eval.indices)) {
    ii <- eval.indices[[jj]]
    solve.out <- .npreghat_solve_eval(
      W = state$W.train,
      w.eval = state$W.eval[ii, ],
      k = kw[, jj],
      ridge.base = 0.0
    )
    if (is.null(solve.out)) {
      stop(sprintf("failed to solve smooth-coefficient effective-weight system at evaluation row %d", ii))
    }
    out[, jj] <- kw[, jj] * drop(state$W.train %*% solve.out$v)
  }

  out
}

.npscoef_weight_matrix <- function(bws, tzdat, ezdat, leave.one.out = FALSE) {
  state <- .npscoef_effective_weight_state(
    bws = bws,
    tzdat = tzdat,
    ezdat = ezdat,
    leave.one.out = leave.one.out
  )
  neval <- nrow(state$z.eval)
  ntrain <- nrow(state$z.train)
  H <- matrix(0.0, nrow = ntrain, ncol = neval)
  for (ii in seq_len(neval))
    H[, ii] <- .npscoef_effective_weight_one(state = state, eval.index = ii)
  H
}

npindexhat <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           exdat = txdat,
           y = NULL,
           output = c("matrix", "apply"),
           s = 0L,
           fd.step = NULL,
           ...){

    output <- match.arg(output)
    .np_semihat_require_class(bws, "sibandwidth", "npindexhat")
    s <- as.integer(s)
    if (length(s) != 1L || is.na(s) || !(s %in% c(0L, 1L)))
      stop("argument 's' must be 0 (fit) or 1 (index derivative)")
    if (!is.null(fd.step)) {
      fd.step <- as.double(fd.step)
      if (length(fd.step) != 1L || is.na(fd.step) || !is.finite(fd.step) || fd.step <= 0)
        stop("argument 'fd.step' must be a positive finite scalar")
    }

    txdat <- toFrame(txdat)
    exdat <- toFrame(exdat)

    npKernelBoundsCheckEval(exdat, bws$xdati$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")

    txdat <- adjustLevels(txdat, bws$xdati)
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)

    txm <- toMatrix(txdat)
    exm <- toMatrix(exdat)

    index.train <- as.vector(txm %*% bws$beta)
    index.eval <- as.vector(exm %*% bws$beta)

    idx.train <- data.frame(index = index.train)
    idx.eval <- data.frame(index = index.eval)
    spec <- .npindex_resolve_spec(bws, where = "npindexhat")

    if (identical(spec$regtype.engine, "lc")) {
      return(.np_indexhat_exact(
        bws = bws,
        idx.train = idx.train,
        idx.eval = idx.eval,
        y = y,
        output = output,
        s = s
      ))
    }
    if (s == 1L) {
      return(.np_indexhat_exact(
        bws = bws,
        idx.train = idx.train,
        idx.eval = idx.eval,
        y = y,
        output = output,
        s = s
      ))
    }
    if (!identical(bws$type, "fixed")) {
      return(.np_indexhat_exact(
        bws = bws,
        idx.train = idx.train,
        idx.eval = idx.eval,
        y = y,
        output = output,
        s = s
      ))
    }

    if (s == 0L) {
      return(.np_indexhat_core(
        bws = bws,
        idx.train = idx.train,
        idx.eval = idx.eval,
        y = y,
        output = output,
        ridge = 0.0
      ))
    }

    step <- if (is.null(fd.step)) {
      span <- diff(range(index.train))
      span.scale <- if (is.finite(span) && span > 0) span else 1.0
      max(1.0e-6, 1.0e-5 * span.scale)
    } else {
      fd.step <- as.double(fd.step)
      if (length(fd.step) != 1L || is.na(fd.step) || !is.finite(fd.step) || fd.step <= 0)
        stop("argument 'fd.step' must be a positive finite scalar")
      fd.step
    }
    if (is.na(step) || !is.finite(step) || step <= 0)
      stop("argument 'fd.step' must be a positive finite scalar")

    if (identical(output, "matrix")) {
      H.plus <- .np_indexhat_core(
        bws = bws,
        idx.train = idx.train,
        idx.eval = data.frame(index = index.eval + step),
        output = "matrix",
        ridge = 0.0
      )
      H.minus <- .np_indexhat_core(
        bws = bws,
        idx.train = idx.train,
        idx.eval = data.frame(index = index.eval - step),
        output = "matrix",
        ridge = 0.0
      )
      return((H.plus - H.minus) / (2.0 * step))
    }

    if (is.null(y))
      stop("argument 'y' is required when output='apply'")

    out.plus <- .np_indexhat_core(
      bws = bws,
      idx.train = idx.train,
      idx.eval = data.frame(index = index.eval + step),
      y = y,
      output = "apply",
      ridge = 0.0
    )
    out.minus <- .np_indexhat_core(
      bws = bws,
      idx.train = idx.train,
      idx.eval = data.frame(index = index.eval - step),
      y = y,
      output = "apply",
      ridge = 0.0
    )
    (out.plus - out.minus) / (2.0 * step)
  }

npplreghat <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tzdat = stop("training data 'tzdat' missing"),
           exdat = txdat,
           ezdat = tzdat,
           y = NULL,
           output = c("apply", "matrix"),
           ...){

    output <- match.arg(output)
    .np_semihat_require_class(bws, "plbandwidth", "npplreghat")
    txdat <- toFrame(txdat)
    tzdat <- toFrame(tzdat)
    exdat <- toFrame(exdat)
    ezdat <- toFrame(ezdat)

    if (nrow(txdat) != nrow(tzdat))
      stop("training data 'txdat' and 'tzdat' must have same number of rows")
    if (nrow(exdat) != nrow(ezdat))
      stop("evaluation data 'exdat' and 'ezdat' must have same number of rows")
    if (ncol(txdat) != ncol(exdat))
      stop("'txdat' and 'exdat' must have same number of columns")

    npKernelBoundsCheckEval(ezdat, bws$zdati$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")

    txdat <- adjustLevels(txdat, bws$xdati)
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    tzdat <- adjustLevels(tzdat, bws$zdati)
    ezdat <- adjustLevels(ezdat, bws$zdati, allowNewCells = TRUE)

    if (is.null(y) && identical(output, "apply"))
      stop("argument 'y' is required when output='apply'")

    x.train.num <- matrix(0.0, nrow = nrow(txdat), ncol = ncol(txdat))
    x.eval.num <- matrix(0.0, nrow = nrow(exdat), ncol = ncol(txdat))

    for (j in seq_len(ncol(txdat))) {
      bw.xj <- bws$bw[[j + 1L]]
      if (is.factor(txdat[[j]])) {
        trj <- adjustLevels(txdat[, j, drop = FALSE], bw.xj$ydati)
        evj <- adjustLevels(exdat[, j, drop = FALSE], bw.xj$ydati, allowNewCells = TRUE)
        lev <- bws$bw[[j + 1L]]$ydati$all.dlev[[1L]]
        x.train.num[, j] <- lev[as.integer(trj[, 1L])]
        x.eval.num[, j] <- lev[as.integer(evj[, 1L])]
      } else {
        x.train.num[, j] <- as.double(txdat[[j]])
        x.eval.num[, j] <- as.double(exdat[[j]])
      }
    }

    n <- nrow(txdat)
    p <- ncol(txdat)
    m <- nrow(exdat)

    resx.train <- matrix(0.0, nrow = n, ncol = p)
    resx.eval <- matrix(0.0, nrow = m, ncol = p)

    if (identical(output, "matrix")) {
      H.y.eval <- npreghat(
        bws = bws$bw$yzbw,
        txdat = tzdat,
        exdat = ezdat,
        output = "matrix"
      )

      for (j in seq_len(p)) {
        xhat.train <- as.vector(npreghat(
          bws = bws$bw[[j + 1L]],
          txdat = tzdat,
          y = x.train.num[, j],
          output = "apply"
        ))
        xhat.eval <- as.vector(npreghat(
          bws = bws$bw[[j + 1L]],
          txdat = tzdat,
          exdat = ezdat,
          y = x.train.num[, j],
          output = "apply"
        ))

        resx.train[, j] <- x.train.num[, j] - xhat.train
        resx.eval[, j] <- x.eval.num[, j] - xhat.eval
      }

      qrR <- qr(resx.train, tol = .Machine$double.eps)
      H.y.train <- npreghat(
        bws = bws$bw$yzbw,
        txdat = tzdat,
        output = "matrix"
      )
      H.y.train[] <- -H.y.train
      diag(H.y.train) <- diag(H.y.train) + 1.0
      A <- qr.coef(qrR, H.y.train)
      A[is.na(A)] <- 0.0
      H <- H.y.eval + resx.eval %*% A
      return(H)
    }

    yy <- as.matrix(y)
    if (nrow(yy) != n)
      stop("number of rows in 'y' must equal number of training rows")

    out <- .npRmpi_with_local_hat_helper({
      for (j in seq_len(p)) {
        xhat.train <- npreghat(
          bws = bws$bw[[j + 1L]],
          txdat = tzdat,
          y = x.train.num[, j],
          output = "apply"
        )
        xhat.eval <- npreghat(
          bws = bws$bw[[j + 1L]],
          txdat = tzdat,
          exdat = ezdat,
          y = x.train.num[, j],
          output = "apply"
        )
        resx.train[, j] <- x.train.num[, j] - as.vector(xhat.train)
        resx.eval[, j] <- x.eval.num[, j] - as.vector(xhat.eval)
      }

      qrR <- qr(resx.train, tol = .Machine$double.eps)

      Hy.eval <- npreghat(
        bws = bws$bw$yzbw,
        txdat = tzdat,
        exdat = ezdat,
        y = yy,
        output = "apply"
      )
      if (!is.matrix(Hy.eval))
        Hy.eval <- matrix(Hy.eval, ncol = ncol(yy))

      Hy.train <- npreghat(
        bws = bws$bw$yzbw,
        txdat = tzdat,
        y = yy,
        output = "apply"
      )
      if (!is.matrix(Hy.train))
        Hy.train <- matrix(Hy.train, ncol = ncol(yy))

      B <- qr.coef(qrR, yy - Hy.train)
      B[is.na(B)] <- 0.0

      Hy.eval + resx.eval %*% B
    })

    if (ncol(out) == 1L) as.vector(out) else out
  }

npscoefhat <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tzdat = NULL,
           exdat = txdat,
           ezdat = tzdat,
           y = NULL,
           output = c("matrix", "apply"),
           ridge = 0.0,
           iterate = FALSE,
           leave.one.out = FALSE,
           ...){

    output <- match.arg(output)
    .np_semihat_require_class(bws, "scbandwidth", "npscoefhat")
    iterate <- npValidateScalarLogical(iterate, "iterate")
    leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")
    ridge <- as.double(ridge)
    if (length(ridge) != 1L || is.na(ridge) || !is.finite(ridge) || ridge < 0)
      stop("argument 'ridge' must be a non-negative finite scalar")
    if (iterate)
      stop("iterate=TRUE is not supported in npscoefhat")

    miss.z <- missing(tzdat) || is.null(tzdat)

    txdat <- toFrame(txdat)
    if (miss.z)
      tzdat <- txdat
    else
      tzdat <- toFrame(tzdat)

    exdat <- toFrame(exdat)
    if (missing(ezdat) || is.null(ezdat))
      ezdat <- if (miss.z) exdat else tzdat
    else
      ezdat <- toFrame(ezdat)

    if (nrow(txdat) != nrow(tzdat))
      stop("'txdat' and 'tzdat' must have same number of training rows")
    if (nrow(exdat) != nrow(ezdat))
      stop("'exdat' and 'ezdat' must have same number of evaluation rows")
    if (ncol(txdat) != ncol(exdat))
      stop("'txdat' and 'exdat' must have same number of columns")

    npKernelBoundsCheckEval(ezdat, bws$zdati$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")

    txdat <- adjustLevels(txdat, bws$xdati)
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    if (!miss.z) {
      tzdat <- adjustLevels(tzdat, bws$zdati)
      ezdat <- adjustLevels(ezdat, bws$zdati, allowNewCells = TRUE)
    }

    X.train <- toMatrix(txdat)
    X.eval <- toMatrix(exdat)

    W.train <- cbind(1.0, X.train)
    W.eval <- cbind(1.0, X.eval)

    n <- nrow(W.train)
    m <- nrow(W.eval)
    spec <- .npscoef_canonical_spec(source = bws, zdat = tzdat, where = "npscoefhat")

    if (identical(output, "apply")) {
      if (is.null(y))
        stop("argument 'y' is required when output='apply'")
      yy <- as.matrix(y)
      if (nrow(yy) != n)
        stop("number of rows in 'y' must equal number of training rows")
      out <- matrix(0.0, nrow = m, ncol = ncol(yy))

      if (identical(spec$regtype.engine, "lc")) {
        state <- .npscoef_effective_weight_state(
          bws = bws,
          tzdat = tzdat,
          ezdat = ezdat,
          leave.one.out = leave.one.out
        )
        chunk.size <- .npscoef_effective_weight_chunk_size(ntrain = n, neval = m)

        for (start in seq.int(1L, m, by = chunk.size)) {
          stopi <- min(m, start + chunk.size - 1L)
          idx <- seq.int(start, stopi)
          kw.chunk <- .npscoef_effective_weight_chunk(state = state, eval.indices = idx)

          for (jj in seq_along(idx)) {
            ii <- idx[[jj]]
            solve.out <- .npreghat_solve_eval(
              W = W.train,
              w.eval = W.eval[ii, ],
              k = kw.chunk[, jj],
              ridge.base = ridge
            )
            if (is.null(solve.out))
              stop(sprintf("failed to solve local hat system at evaluation row %d", ii))
            h.row <- kw.chunk[, jj] * drop(W.train %*% solve.out$v)
            out[ii, ] <- drop(crossprod(h.row, yy))
          }
        }
      } else {
        lp_state <- .npscoef_lp_state(
          bws = bws,
          tzdat = tzdat,
          ezdat = ezdat,
          leave.one.out = leave.one.out,
          where = "npscoefhat"
        )
        tensor.train <- .npscoef_row_tensor_design(W.train, lp_state$W.train)
        tensor.eval <- .npscoef_row_tensor_design(W.eval, lp_state$W.eval)
        chunk.size <- .npscoef_effective_weight_chunk_size(ntrain = n, neval = m)

        for (start in seq.int(1L, m, by = chunk.size)) {
          stopi <- min(m, start + chunk.size - 1L)
          idx <- seq.int(start, stopi)
          kw.chunk <- .np_kernel_weights_direct(
            txdat = lp_state$z.train,
            exdat = lp_state$z.eval[idx, , drop = FALSE],
            bws = lp_state$rbw,
            bandwidth.divide = TRUE,
            kernel.pow = 1.0
          )

          if (leave.one.out) {
            for (jj in seq_along(idx))
              kw.chunk[idx[[jj]], jj] <- 0.0
          }

          for (jj in seq_along(idx)) {
            ii <- idx[[jj]]
            solve.out <- .npreghat_solve_eval(
              W = tensor.train,
              w.eval = tensor.eval[ii, ],
              k = kw.chunk[, jj],
              ridge.base = ridge
            )
            if (is.null(solve.out))
              stop(sprintf("failed to solve smooth-coefficient local hat system at evaluation row %d", ii))
            h.row <- kw.chunk[, jj] * drop(tensor.train %*% solve.out$v)
            out[ii, ] <- drop(crossprod(h.row, yy))
          }
        }
      }

      if (ncol(out) == 1L) return(as.vector(out))
      return(out)
    }

    H <- matrix(0.0, nrow = m, ncol = n)

    if (identical(spec$regtype.engine, "lc")) {
      state <- .npscoef_effective_weight_state(
        bws = bws,
        tzdat = tzdat,
        ezdat = ezdat,
        leave.one.out = leave.one.out
      )
      chunk.size <- .npscoef_effective_weight_chunk_size(ntrain = n, neval = m)

      for (start in seq.int(1L, m, by = chunk.size)) {
        stopi <- min(m, start + chunk.size - 1L)
        idx <- seq.int(start, stopi)
        kw.chunk <- .npscoef_effective_weight_chunk(state = state, eval.indices = idx)

        for (jj in seq_along(idx)) {
          ii <- idx[[jj]]
          solve.out <- .npreghat_solve_eval(
            W = W.train,
            w.eval = W.eval[ii, ],
            k = kw.chunk[, jj],
            ridge.base = ridge
          )
          if (is.null(solve.out))
            stop(sprintf("failed to solve local hat system at evaluation row %d", ii))
          H[ii, ] <- kw.chunk[, jj] * drop(W.train %*% solve.out$v)
        }
      }
    } else {
      lp_state <- .npscoef_lp_state(
        bws = bws,
        tzdat = tzdat,
        ezdat = ezdat,
        leave.one.out = leave.one.out,
        where = "npscoefhat"
      )
      tensor.train <- .npscoef_row_tensor_design(W.train, lp_state$W.train)
      tensor.eval <- .npscoef_row_tensor_design(W.eval, lp_state$W.eval)
      chunk.size <- .npscoef_effective_weight_chunk_size(ntrain = n, neval = m)

      for (start in seq.int(1L, m, by = chunk.size)) {
        stopi <- min(m, start + chunk.size - 1L)
        idx <- seq.int(start, stopi)
        kw.chunk <- .np_kernel_weights_direct(
          txdat = lp_state$z.train,
          exdat = lp_state$z.eval[idx, , drop = FALSE],
          bws = lp_state$rbw,
          bandwidth.divide = TRUE,
          kernel.pow = 1.0
        )

        if (leave.one.out) {
          for (jj in seq_along(idx))
            kw.chunk[idx[[jj]], jj] <- 0.0
        }

        for (jj in seq_along(idx)) {
          ii <- idx[[jj]]
          solve.out <- .npreghat_solve_eval(
            W = tensor.train,
            w.eval = tensor.eval[ii, ],
            k = kw.chunk[, jj],
            ridge.base = ridge
          )
          if (is.null(solve.out))
            stop(sprintf("failed to solve smooth-coefficient local hat system at evaluation row %d", ii))
          H[ii, ] <- kw.chunk[, jj] * drop(tensor.train %*% solve.out$v)
        }
      }
    }

    H
  }
