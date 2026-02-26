.np_semihat_require_class <- function(bws, class.name, fn.name) {
  if (!inherits(bws, class.name))
    stop(sprintf("argument 'bws' must inherit from class '%s' in %s()", class.name, fn.name))
}

.np_indexhat_rbw <- function(bws, idx.train) {
  regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype
  args <- list(
    xdat = idx.train,
    ydat = rep.int(0.0, nrow(idx.train)),
    bws = as.double(bws$bw),
    regtype = regtype,
    bandwidth.compute = FALSE,
    ckertype = bws$ckertype,
    ckerorder = bws$ckerorder
  )
  if (!is.null(bws$basis))
    args$basis <- bws$basis
  if (!is.null(bws$degree))
    args$degree <- bws$degree
  if (!is.null(bws$bernstein.basis))
    args$bernstein.basis <- bws$bernstein.basis
  do.call(npregbw, args)
}

.npscoef_make_regbw <- function(bws, zdat, bw = bws$bw) {
  regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype
  args <- list(
    xdat = zdat,
    ydat = rep.int(0.0, nrow(zdat)),
    bws = as.double(bw),
    regtype = regtype,
    bwscaling = isTRUE(bws$scaling),
    bandwidth.compute = FALSE,
    bwtype = if (is.null(bws$type)) "fixed" else bws$type,
    ckertype = if (is.null(bws$ckertype)) "gaussian" else bws$ckertype,
    ckerorder = if (is.null(bws$ckerorder)) 2L else bws$ckerorder,
    ckerbound = if (is.null(bws$ckerbound)) "none" else bws$ckerbound,
    ckerlb = bws$ckerlb,
    ckerub = bws$ckerub,
    ukertype = if (is.null(bws$ukertype)) "aitchisonaitken" else bws$ukertype,
    okertype = if (is.null(bws$okertype)) "liracine" else bws$okertype
  )
  if (!is.null(bws$basis))
    args$basis <- bws$basis
  if (!is.null(bws$degree))
    args$degree <- bws$degree
  if (!is.null(bws$bernstein.basis))
    args$bernstein.basis <- bws$bernstein.basis
  do.call(npregbw, args)
}

.npscoef_weight_matrix <- function(bws, tzdat, ezdat, leave.one.out = FALSE) {
  tzdat <- toFrame(tzdat)
  ezdat <- toFrame(ezdat)
  leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")
  regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype

  if (identical(regtype, "lc")) {
    return(.np_kernel_weights_direct(
      bws = bws,
      txdat = tzdat,
      exdat = ezdat,
      leave.one.out = leave.one.out,
      bandwidth.divide = TRUE,
      kernel.pow = 1.0
    ))
  }

  rbw <- .npscoef_make_regbw(bws = bws, zdat = tzdat)
  H <- npreghat(
    bws = rbw,
    txdat = tzdat,
    exdat = ezdat,
    output = "matrix",
    leave.one.out = leave.one.out
  )
  if (!is.matrix(H))
    H <- matrix(H, nrow = nrow(ezdat))
  t(H)
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
           ridge = 1.0e-12,
           iterate = FALSE,
           leave.one.out = FALSE,
           ...){

    output <- match.arg(output)
    .np_semihat_require_class(bws, "scbandwidth", "npscoefhat")
    iterate <- npValidateScalarLogical(iterate, "iterate")
    leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")
    ridge <- as.double(ridge)
    if (length(ridge) != 1L || is.na(ridge) || !is.finite(ridge) || ridge <= 0)
      stop("argument 'ridge' must be a positive finite scalar")
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
