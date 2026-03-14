.npcdhat_make_xbw <- function(bws, txdat) {
  regtype <- if (is.null(bws$regtype.engine)) bws$regtype else bws$regtype.engine
  basis <- if (is.null(bws$basis.engine)) bws$basis else bws$basis.engine
  degree <- if (is.null(bws$degree.engine)) bws$degree else bws$degree.engine
  bernstein <- if (is.null(bws$bernstein.basis.engine)) bws$bernstein.basis else bws$bernstein.basis.engine

  xbw.args <- list(
    xdat = txdat,
    ydat = rep.int(0.0, nrow(txdat)),
    bws = bws$xbw,
    regtype = regtype,
    bwtype = bws$type,
    bandwidth.compute = FALSE,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = bws$cxkerlb,
    ckerub = bws$cxkerub,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype
  )

  if (!is.null(basis))
    xbw.args$basis <- basis
  if (!is.null(degree))
    xbw.args$degree <- degree
  if (!is.null(bernstein))
    xbw.args$bernstein.basis <- bernstein

  do.call(npregbw, xbw.args)
}

.npcdhat_make_xkbw <- function(bws, txdat) {
  kbandwidth.numeric(
    bw = bws$xbw,
    bwscaling = FALSE,
    bwtype = bws$type,
    ckertype = bws$cxkertype,
    ckerorder = bws$cxkerorder,
    ckerbound = bws$cxkerbound,
    ckerlb = bws$cxkerlb,
    ckerub = bws$cxkerub,
    ukertype = bws$uxkertype,
    okertype = bws$oxkertype,
    nobs = nrow(txdat),
    xdati = untangle(txdat),
    ydati = NULL,
    xnames = names(txdat),
    ynames = NULL
  )
}

.npcdhat_make_ybw <- function(bws, tydat) {
  kbandwidth.numeric(
    bw = bws$ybw,
    bwscaling = FALSE,
    bwtype = bws$type,
    ckertype = bws$cykertype,
    ckerorder = bws$cykerorder,
    ckerbound = bws$cykerbound,
    ckerlb = bws$cykerlb,
    ckerub = bws$cykerub,
    ukertype = bws$uykertype,
    okertype = bws$oykertype,
    nobs = nrow(tydat),
    xdati = untangle(tydat),
    ydati = NULL,
    xnames = names(tydat),
    ynames = NULL
  )
}

.np_direct_operator_matrix <- function(kbw, txdat, exdat, operator, where) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  op.info <- .np_operator_kernel_weight_scale(
    bws = kbw,
    operator = operator,
    nvars = ncol(txdat),
    where = where
  )
  kw <- .np_kernel_weights_direct(
    bws = op.info$bws,
    txdat = txdat,
    exdat = exdat,
    bandwidth.divide = TRUE,
    operator = op.info$operator
  )

  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(txdat))
  if (nrow(kw) != nrow(txdat) || ncol(kw) != nrow(exdat))
    stop(sprintf("%s returned unexpected operator shape", where))

  t(kw) / op.info$scale
}

.np_direct_operator_apply <- function(kbw, txdat, exdat, operator, rhs, where) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  rhs <- as.matrix(rhs)
  storage.mode(rhs) <- "double"

  if (nrow(rhs) != nrow(txdat))
    stop(sprintf("%s received RHS with unexpected number of rows", where))

  op.info <- .np_operator_kernel_weight_scale(
    bws = kbw,
    operator = operator,
    nvars = ncol(txdat),
    where = where
  )
  kw <- .np_kernel_weights_direct(
    bws = op.info$bws,
    txdat = txdat,
    exdat = exdat,
    bandwidth.divide = TRUE,
    operator = op.info$operator
  )

  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(txdat))
  if (nrow(kw) != nrow(txdat) || ncol(kw) != nrow(exdat))
    stop(sprintf("%s returned unexpected operator shape", where))

  crossprod(kw, rhs) / op.info$scale
}

.np_exact_operator_apply <- function(kbw, txdat, exdat, operator, rhs, where) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  rhs <- as.matrix(rhs)
  storage.mode(rhs) <- "double"

  if (nrow(rhs) != nrow(txdat))
    stop(sprintf("%s received RHS with unexpected number of rows", where))

  out <- npksum(
    bws = kbw,
    txdat = txdat,
    exdat = exdat,
    tydat = rhs,
    operator = operator,
    bandwidth.divide = TRUE
  )$ksum

  if (!is.matrix(out))
    out <- matrix(out, nrow = nrow(exdat), ncol = ncol(rhs))

  if (ncol(out) == 1L)
    out[, 1L, drop = FALSE]
  else
    out
}

.np_exact_operator_matrix <- function(kbw, txdat, exdat, operator, where) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  .np_exact_operator_apply(
    kbw = kbw,
    txdat = txdat,
    exdat = exdat,
    operator = operator,
    rhs = diag(nrow(txdat)),
    where = where
  )
}

.npcdhat_make_kernel_matrix <- function(kbw, txdat, exdat, operator) {
  if (!identical(kbw$type, "fixed")) {
    return(.np_exact_operator_matrix(
      kbw = kbw,
      txdat = txdat,
      exdat = exdat,
      operator = operator,
      where = "conditional hat exact kernel matrix"
    ))
  }

  .np_direct_operator_matrix(
    kbw = kbw,
    txdat = txdat,
    exdat = exdat,
    operator = operator,
    where = "conditional hat direct kernel matrix"
  )
}

.npcdhat_ratio_matrix <- function(bws, txdat, tydat, exdat, eydat, operator) {
  xkbw <- .npcdhat_make_xkbw(bws = bws, txdat = txdat)
  ybw <- .npcdhat_make_ybw(bws = bws, tydat = tydat)
  Kx <- .npcdhat_make_kernel_matrix(
    kbw = xkbw,
    txdat = txdat,
    exdat = exdat,
    operator = rep.int("normal", ncol(txdat))
  )
  Ky <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = tydat,
    exdat = eydat,
    operator = rep.int(operator, ncol(tydat))
  )

  denom <- rowSums(Kx) / nrow(txdat)
  sweep((Kx * Ky) / nrow(txdat), 1L, pmax(denom, .Machine$double.eps), "/")
}

.npcdhat_exact_matrix <- function(bws, txdat, tydat, exdat, eydat, operator) {
  if (identical(bws$type, "adaptive_nn")) {
    return(.npcdhat_ratio_matrix(
      bws = bws,
      txdat = txdat,
      tydat = tydat,
      exdat = exdat,
      eydat = eydat,
      operator = operator
    ))
  }

  xbw <- .npcdhat_make_xbw(bws = bws, txdat = txdat)
  ybw <- .npcdhat_make_ybw(bws = bws, tydat = tydat)
  Gy <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = tydat,
    exdat = eydat,
    operator = rep.int(operator, ncol(tydat))
  )

  H <- matrix(NA_real_, nrow = nrow(exdat), ncol = nrow(txdat))
  for (i in seq_len(nrow(exdat))) {
    x.row <- .npreghat_exact_matrix_from_core(
      bws = xbw,
      txdat = txdat,
      exdat = exdat[i, , drop = FALSE]
    )
    H[i, ] <- as.vector(x.row[1L, ]) * as.vector(Gy[i, ])
  }

  H
}

.npcdhat_exact_apply <- function(bws, txdat, tydat, exdat, eydat, rhs, operator) {
  if (identical(bws$type, "adaptive_nn")) {
    H <- .npcdhat_ratio_matrix(
      bws = bws,
      txdat = txdat,
      tydat = tydat,
      exdat = exdat,
      eydat = eydat,
      operator = operator
    )
    return(H %*% rhs)
  }

  xbw <- .npcdhat_make_xbw(bws = bws, txdat = txdat)
  ybw <- .npcdhat_make_ybw(bws = bws, tydat = tydat)
  Gy <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = tydat,
    exdat = eydat,
    operator = rep.int(operator, ncol(tydat))
  )

  out <- matrix(0.0, nrow = nrow(exdat), ncol = ncol(rhs))
  for (i in seq_len(nrow(exdat))) {
    x.row <- .npreghat_exact_matrix_from_core(
      bws = xbw,
      txdat = txdat,
      exdat = exdat[i, , drop = FALSE]
    )
    h.row <- as.vector(x.row[1L, ]) * as.vector(Gy[i, ])
    out[i, ] <- drop(crossprod(h.row, rhs))
  }

  out
}

.npcdhat_core <- function(bws,
                          txdat,
                          tydat,
                          exdat,
                          eydat,
                          y,
                          output,
                          operator,
                          class_name,
                          where) {
  output <- match.arg(output, c("matrix", "apply"))

  if (xor(is.null(exdat), is.null(eydat)))
    stop("evaluation data must be supplied for both 'exdat' and 'eydat'")

  no.exy <- is.null(exdat)

  txdat <- toFrame(txdat)
  tydat <- toFrame(tydat)

  if (!no.exy) {
    exdat <- toFrame(exdat)
    eydat <- toFrame(eydat)

    if (!(txdat %~% exdat))
      stop("'txdat' and 'exdat' are not similar data frames!")
    if (!(tydat %~% eydat))
      stop("'tydat' and 'eydat' are not similar data frames!")
  }

  if (!is.null(y)) {
    if (is.factor(y) || is.vector(y)) {
      y <- matrix(as.double(y), ncol = 1L)
    } else {
      y <- as.matrix(y)
      storage.mode(y) <- "double"
    }

    if (nrow(y) != nrow(txdat))
      stop("number of rows in 'y' must match the number of training rows in 'txdat'")
  }

  keep.rows <- rep_len(TRUE, nrow(txdat))
  rows.omit.train <- attr(stats::na.omit(data.frame(txdat, tydat)), "na.action")
  if (!is.null(y))
    rows.omit.train <- union(rows.omit.train, attr(stats::na.omit(as.data.frame(y)), "na.action"))

  if (length(rows.omit.train) > 0L)
    keep.rows[as.integer(rows.omit.train)] <- FALSE

  if (!any(keep.rows))
    stop("Training data has no rows without NAs")

  txdat <- txdat[keep.rows, , drop = FALSE]
  tydat <- tydat[keep.rows, , drop = FALSE]
  if (!is.null(y))
    y <- y[keep.rows, , drop = FALSE]

  if (!no.exy) {
    keep.eval <- rep_len(TRUE, nrow(exdat))
    rows.omit.eval <- attr(stats::na.omit(data.frame(exdat, eydat)), "na.action")
    if (length(rows.omit.eval) > 0L)
      keep.eval[as.integer(rows.omit.eval)] <- FALSE

    if (!any(keep.eval))
      stop("Evaluation data has no rows without NAs")

    exdat <- exdat[keep.eval, , drop = FALSE]
    eydat <- eydat[keep.eval, , drop = FALSE]
  } else {
    exdat <- txdat
    eydat <- tydat
  }

  if (identical(output, "apply")) {
    if (is.null(y))
      stop("argument 'y' is required when output='apply'")

    out <- .npcdhat_exact_apply(
      bws = bws,
      txdat = txdat,
      tydat = tydat,
      exdat = exdat,
      eydat = eydat,
      rhs = y,
      operator = operator
    )
    if (ncol(out) == 1L)
      return(as.vector(out))
    return(out)
  }

  H <- .npcdhat_exact_matrix(
    bws = bws,
    txdat = txdat,
    tydat = tydat,
    exdat = exdat,
    eydat = eydat,
    operator = operator
  )

  class(H) <- c(class_name, "matrix")
  attr(H, "bws") <- bws
  attr(H, "txdat") <- txdat
  attr(H, "tydat") <- tydat
  attr(H, "exdat") <- exdat
  attr(H, "eydat") <- eydat
  attr(H, "trainiseval") <- no.exy
  attr(H, "rows.omit.train") <- rows.omit.train
  attr(H, "call") <- match.call(expand.dots = FALSE)

  if (!is.null(y)) {
    Hy <- H %*% y
    if (ncol(Hy) == 1L)
      Hy <- as.vector(Hy)
    attr(H, "Hy") <- Hy
  }

  H
}
