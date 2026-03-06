.np_semihat_require_class <- function(bws, class.name, fn.name) {
  if (!inherits(bws, class.name))
    stop(sprintf("argument 'bws' must inherit from class '%s' in %s()", class.name, fn.name))
}

.np_semihat_make_regbw_args <- function(source, xdat, ydat, bw) {
  regtype <- if (is.null(source$regtype)) "lc" else as.character(source$regtype)
  xdat <- toFrame(xdat)
  ncon.target <- sum(untangle(xdat)$icon)

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
    regtype = regtype,
    bandwidth.compute = FALSE
  )

  # Forward the full compatible estimator/bandwidth option contract.
  if (!is.null(source$basis))
    args$basis <- source$basis
  if (!is.null(source$degree))
    args$degree <- source$degree
  if (!is.null(source$bernstein.basis))
    args$bernstein.basis <- source$bernstein.basis
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

.np_indexhat_rbw <- function(bws, idx.train) {
  do.call(npregbw, .np_semihat_make_regbw_args(
    source = bws,
    xdat = idx.train,
    ydat = rep.int(0.0, nrow(idx.train)),
    bw = bws$bw
  ))
}

.npscoef_make_regbw <- function(bws, zdat, bw = bws$bw) {
  do.call(npregbw, .np_semihat_make_regbw_args(
    source = bws,
    xdat = zdat,
    ydat = rep.int(0.0, nrow(zdat)),
    bw = bw
  ))
}

.npscoef_effective_weight_state <- function(bws, tzdat, ezdat, leave.one.out = FALSE) {
  tzdat <- toFrame(tzdat)
  ezdat <- toFrame(ezdat)
  leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")
  regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype
  same.eval <- isTRUE(all.equal(tzdat, ezdat, check.attributes = FALSE))

  if (leave.one.out && !same.eval) {
    stop("leave.one.out=TRUE requires evaluation 'z' data to match training 'z' data")
  }

  state <- list(
    bws = bws,
    regtype = regtype,
    leave.one.out = leave.one.out
  )

  if (identical(regtype, "lc")) {
    state$z.train <- adjustLevels(tzdat, bws$zdati)
    state$z.eval <- if (leave.one.out) state$z.train else adjustLevels(ezdat, bws$zdati, allowNewCells = TRUE)
    return(state)
  }

  rbw <- .npscoef_make_regbw(bws = bws, zdat = tzdat)
  regtype.rbw <- if (is.null(rbw$regtype)) "lc" else as.character(rbw$regtype)
  degree.rbw <- if (identical(regtype.rbw, "lc")) {
    rep.int(0L, rbw$ncon)
  } else if (identical(regtype.rbw, "ll")) {
    rep.int(1L, rbw$ncon)
  } else {
    npValidateGlpDegree(regtype = "lp", degree = rbw$degree, ncon = rbw$ncon)
  }
  basis.rbw <- npValidateLpBasis(
    regtype = "lp",
    basis = if (is.null(rbw$basis)) "glp" else rbw$basis
  )
  bernstein.rbw <- npValidateGlpBernstein(
    regtype = "lp",
    bernstein.basis = isTRUE(rbw$bernstein.basis)
  )

  z.train <- adjustLevels(tzdat, rbw$xdati)
  z.eval <- if (leave.one.out) z.train else adjustLevels(ezdat, rbw$xdati, allowNewCells = TRUE)
  W.train <- W.lp(
    xdat = z.train,
    degree = degree.rbw,
    basis = basis.rbw,
    bernstein.basis = bernstein.rbw
  )
  W.eval <- W.lp(
    xdat = z.train,
    exdat = if (leave.one.out) NULL else z.eval,
    degree = degree.rbw,
    basis = basis.rbw,
    bernstein.basis = bernstein.rbw
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

    txdat <- toFrame(txdat)
    exdat <- toFrame(exdat)

    txdat <- adjustLevels(txdat, bws$xdati)
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)

    txm <- toMatrix(txdat)
    exm <- toMatrix(exdat)

    index.train <- as.vector(txm %*% bws$beta)
    index.eval <- as.vector(exm %*% bws$beta)

    idx.train <- data.frame(index = index.train)
    idx.eval <- data.frame(index = index.eval)
    rbw <- .np_indexhat_rbw(bws = bws, idx.train = idx.train)

    if (s == 0L) {
      return(npreghat(
        bws = rbw,
        txdat = idx.train,
        exdat = idx.eval,
        y = y,
        output = output,
        ...
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

    H.plus <- npreghat(
      bws = rbw,
      txdat = idx.train,
      exdat = data.frame(index = index.eval + step),
      output = "matrix",
      ...
    )
    H.minus <- npreghat(
      bws = rbw,
      txdat = idx.train,
      exdat = data.frame(index = index.eval - step),
      output = "matrix",
      ...
    )
    H.deriv <- (H.plus - H.minus) / (2.0 * step)

    if (identical(output, "matrix"))
      return(H.deriv)

    if (is.null(y))
      stop("argument 'y' is required when output='apply'")
    yy <- as.matrix(y)
    if (nrow(yy) != nrow(txdat))
      stop("number of rows in 'y' must equal number of training rows")
    out <- H.deriv %*% yy
    if (ncol(out) == 1L) as.vector(out) else out
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

    txdat <- adjustLevels(txdat, bws$xdati)
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    tzdat <- adjustLevels(tzdat, bws$zdati)
    ezdat <- adjustLevels(ezdat, bws$zdati, allowNewCells = TRUE)

    if (is.null(y) && identical(output, "apply"))
      stop("argument 'y' is required when output='apply'")

    x.train.num <- matrix(0.0, nrow = nrow(txdat), ncol = ncol(txdat))
    x.eval.num <- matrix(0.0, nrow = nrow(exdat), ncol = ncol(txdat))

    for (j in seq_len(ncol(txdat))) {
      if (is.factor(txdat[[j]])) {
        trj <- adjustLevels(txdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati)
        evj <- adjustLevels(exdat[, j, drop = FALSE], bws$bw[[j + 1L]]$ydati, allowNewCells = TRUE)
        lev <- bws$bw[[j + 1L]]$ydati$all.dlev[[1L]]
        x.train.num[, j] <- lev[as.integer(trj[, 1L])]
        x.eval.num[, j] <- lev[as.integer(evj[, 1L])]
      } else {
        x.train.num[, j] <- as.double(txdat[[j]])
        x.eval.num[, j] <- as.double(exdat[[j]])
      }
    }

    H.y.eval <- npreghat(
      bws = bws$bw$yzbw,
      txdat = tzdat,
      exdat = ezdat,
      output = "matrix"
    )

    n <- nrow(txdat)
    p <- ncol(txdat)
    m <- nrow(exdat)

    resx.train <- matrix(0.0, nrow = n, ncol = p)
    resx.eval <- matrix(0.0, nrow = m, ncol = p)

    for (j in seq_len(p)) {
      H.x.eval <- npreghat(
        bws = bws$bw[[j + 1L]],
        txdat = tzdat,
        exdat = ezdat,
        output = "matrix"
      )
      xhat.train <- as.vector(npreghat(
        bws = bws$bw[[j + 1L]],
        txdat = tzdat,
        y = x.train.num[, j],
        output = "apply"
      ))
      xhat.eval <- as.vector(H.x.eval %*% x.train.num[, j])

      resx.train[, j] <- x.train.num[, j] - xhat.train
      resx.eval[, j] <- x.eval.num[, j] - xhat.eval
    }

    qrR <- qr(resx.train, tol = .Machine$double.eps)

    if (identical(output, "matrix")) {
      H.y.train <- npreghat(
        bws = bws$bw$yzbw,
        txdat = tzdat,
        output = "matrix"
      )
      A <- qr.coef(qrR, diag(n) - H.y.train)
      A[is.na(A)] <- 0.0
      H <- H.y.eval + resx.eval %*% A
      return(H)
    }

    yy <- as.matrix(y)
    if (nrow(yy) != n)
      stop("number of rows in 'y' must equal number of training rows")

    Hy.eval <- H.y.eval %*% yy
    if (is.null(dim(Hy.eval)))
      Hy.eval <- matrix(Hy.eval, ncol = 1L)

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

    out <- Hy.eval + resx.eval %*% B
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

    kw <- .npscoef_weight_matrix(
      bws = bws,
      tzdat = tzdat,
      ezdat = ezdat,
      leave.one.out = leave.one.out
    )

    n <- nrow(W.train)
    m <- ncol(kw)

    if (identical(output, "matrix")) {
      H <- matrix(0.0, nrow = m, ncol = n)
    } else {
      if (is.null(y))
        stop("argument 'y' is required when output='apply'")
      yy <- as.matrix(y)
      if (nrow(yy) != n)
        stop("number of rows in 'y' must equal number of training rows")
      out <- matrix(0.0, nrow = m, ncol = ncol(yy))
    }

    for (i in seq_len(m)) {
      solve.out <- .npreghat_solve_eval(
        W = W.train,
        w.eval = W.eval[i, ],
        k = kw[, i],
        ridge.base = ridge
      )
      if (is.null(solve.out))
        stop(sprintf("failed to solve local hat system at evaluation row %d", i))
      h.row <- kw[, i] * drop(W.train %*% solve.out$v)

      if (identical(output, "matrix")) {
        H[i, ] <- h.row
      } else {
        out[i, ] <- drop(crossprod(h.row, yy))
      }
    }

    if (identical(output, "matrix"))
      return(H)
    if (ncol(out) == 1L) as.vector(out) else out
  }
