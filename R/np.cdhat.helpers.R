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

.npcdhat_make_xhat_matrix <- function(bws, txdat, exdat) {
  xbw <- .npcdhat_make_xbw(bws = bws, txdat = txdat)
  spec <- npCanonicalConditionalRegSpec(
    regtype = if (is.null(xbw$regtype)) "lc" else as.character(xbw$regtype),
    basis = if (is.null(xbw$basis)) "glp" else xbw$basis,
    degree = xbw$degree,
    bernstein.basis = isTRUE(xbw$bernstein.basis),
    ncon = xbw$ncon,
    where = "npcdhat"
  )
  regtype <- if (is.null(xbw$regtype.engine)) {
    as.character(spec$regtype.engine)
  } else {
    as.character(xbw$regtype.engine)
  }
  basis <- if (is.null(xbw$basis.engine)) spec$basis.engine else xbw$basis.engine
  degree <- if (is.null(xbw$degree.engine)) as.integer(spec$degree.engine) else as.integer(xbw$degree.engine)
  bernstein <- if (is.null(xbw$bernstein.basis.engine)) {
    isTRUE(spec$bernstein.basis.engine)
  } else {
    isTRUE(xbw$bernstein.basis.engine)
  }

  if (identical(regtype, "lc")) {
    return(.npreghat_exact_lc_matrix_from_kernel_weights(
      bws = xbw,
      txdat = txdat,
      exdat = exdat
    ))
  }

  if (identical(regtype, "ll")) {
    return(.npreghat_exact_ll_matrix_from_kernel_weights(
      bws = xbw,
      txdat = txdat,
      exdat = exdat,
      s = integer(0)
    ))
  }

  if (identical(regtype, "lp")) {
    if (identical(xbw$type, "generalized_nn") && any(degree > 1L)) {
      H <- matrix(NA_real_, nrow = nrow(exdat), ncol = nrow(txdat))
      for (i in seq_len(nrow(exdat))) {
        H[i, ] <- .npreghat_exact_matrix_from_core(
          bws = xbw,
          txdat = txdat,
          exdat = exdat[i, , drop = FALSE]
        )[1L, ]
      }
      return(H)
    }

    return(.npreghat_exact_lp_matrix_from_kernel_weights(
      bws = xbw,
      txdat = txdat,
      exdat = exdat,
      s = integer(0),
      basis = basis,
      degree = degree,
      bernstein.basis = bernstein
    ))
  }

  .npreghat_exact_matrix_from_core(
    bws = xbw,
    txdat = txdat,
    exdat = exdat
  )
}

.np_operator_kernel_weight_scale <- function(bws, operator, nvars, where) {
  if (!isa(bws, "kbandwidth"))
    bws <- kbandwidth(bws)

  operator <- as.character(operator)
  if (length(operator) == 1L)
    operator <- rep.int(operator, nvars)
  if (length(operator) != nvars)
    stop(sprintf("%s requires one operator per column", where))

  bw.scale <- 1.0
  if (bws$ncon > 0L) {
    con.ops <- operator[bws$icon]
    if (any(con.ops == "normal"))
      bw.scale <- prod(bws$bw[bws$icon][con.ops == "normal"])
  }

  list(bws = bws, scale = bw.scale, operator = operator)
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

.np_local_operator_ksum <- function(kbw,
                                    txdat,
                                    exdat,
                                    operator,
                                    rhs,
                                    return.kernel.weights = FALSE,
                                    where) {
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
  bws <- op.info$bws

  txdat <- adjustLevels(txdat, bws$xdati, allowNewCells = TRUE)
  exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
  npKernelBoundsCheckEval(exdat, bws$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")

  txm <- toMatrix(txdat)
  exm <- toMatrix(exdat)
  tuno <- txm[, bws$iuno, drop = FALSE]
  tcon <- txm[, bws$icon, drop = FALSE]
  tord <- txm[, bws$iord, drop = FALSE]
  euno <- exm[, bws$iuno, drop = FALSE]
  econ <- exm[, bws$icon, drop = FALSE]
  eord <- exm[, bws$iord, drop = FALSE]

  tnrow <- nrow(txdat)
  enrow <- nrow(exdat)
  nkw <- if (isTRUE(return.kernel.weights)) tnrow * enrow else 0L

  operator.num <- ALL_OPERATORS[op.info$operator]
  myopti <- list(
    num_obs_train = tnrow,
    num_obs_eval = enrow,
    num_uno = bws$nuno,
    num_ord = bws$nord,
    num_con = bws$ncon,
    int_LARGE_SF = SF_ARB,
    BANDWIDTH_reg_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN
    ),
    int_MINIMIZE_IO = if (isTRUE(getOption("np.messages"))) IO_MIN_FALSE else IO_MIN_TRUE,
    kerneval = switch(bws$ckertype,
      gaussian = CKER_GAUSS + bws$ckerorder / 2 - 1,
      epanechnikov = CKER_EPAN + bws$ckerorder / 2 - 1,
      uniform = CKER_UNI,
      "truncated gaussian" = CKER_TGAUSS
    ),
    ukerneval = switch(bws$ukertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR
    ),
    okerneval = switch(bws$okertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR,
      nliracine = OKER_NLR,
      racineliyan = OKER_RLY
    ),
    miss.ex = FALSE,
    leave.one.out = FALSE,
    bandwidth.divide = TRUE,
    mcv.numRow = attr(bws$xmcv, "num.row"),
    wncol = 0L,
    yncol = ncol(rhs),
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    return.kernel.weights = isTRUE(return.kernel.weights),
    permutation.operator = PERMUTATION_OPERATORS[["none"]],
    compute.score = FALSE,
    compute.ocg = FALSE,
    suppress.parallel = 1L
  )

  cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

  asDouble <- function(data) {
    if (is.null(data)) as.double(0.0) else as.double(data)
  }

  myout <- .Call(
    "C_np_kernelsum",
    asDouble(tuno), asDouble(tord), asDouble(tcon),
    asDouble(rhs), as.double(0.0),
    asDouble(euno), asDouble(eord), asDouble(econ),
    as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
    as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
    as.integer(c(operator.num[bws$icon], operator.num[bws$iuno], operator.num[bws$iord])),
    as.integer(myopti), as.double(1.0),
    as.integer(enrow),
    as.integer(0L),
    as.integer(nkw),
    as.double(cker.bounds.c$lb),
    as.double(cker.bounds.c$ub),
    PACKAGE = "npRmpi"
  )

  out <- myout[["ksum"]]
  if (!is.matrix(out))
    out <- matrix(out, nrow = enrow, ncol = ncol(rhs))

  kw <- NULL
  if (isTRUE(return.kernel.weights))
    kw <- matrix(as.double(myout[["kernel.weights"]]), nrow = tnrow, ncol = enrow)

  list(ksum = out, kernel.weights = kw)
}

.np_exact_operator_apply <- function(kbw, txdat, exdat, operator, rhs, where) {
  rhs <- as.matrix(rhs)
  storage.mode(rhs) <- "double"

  out <- .np_local_operator_ksum(
    kbw = kbw,
    txdat = txdat,
    exdat = exdat,
    operator = operator,
    rhs = rhs,
    return.kernel.weights = FALSE,
    where = where
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

  probe <- .np_local_operator_ksum(
    kbw = kbw,
    txdat = txdat,
    exdat = exdat,
    operator = operator,
    rhs = matrix(1.0, nrow = nrow(txdat), ncol = 1L),
    return.kernel.weights = TRUE,
    where = where
  )

  kw <- probe$kernel.weights
  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(txdat))
  if (nrow(kw) != nrow(txdat) || ncol(kw) != nrow(exdat))
    stop(sprintf("%s returned unexpected kernel-weight shape", where))

  H.raw <- t(kw)
  row.sum <- drop(H.raw %*% rep.int(1.0, nrow(txdat)))
  row.sum[abs(row.sum) < .Machine$double.xmin] <- .Machine$double.xmin
  row.scale <- as.vector(probe$ksum[, 1L]) / row.sum

  sweep(H.raw, 1L, row.scale, "*")
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

  ybw <- .npcdhat_make_ybw(bws = bws, tydat = tydat)
  Hx <- .npcdhat_make_xhat_matrix(
    bws = bws,
    txdat = txdat,
    exdat = exdat
  )
  Gy <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = tydat,
    exdat = eydat,
    operator = rep.int(operator, ncol(tydat))
  )

  Hx * Gy
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

  ybw <- .npcdhat_make_ybw(bws = bws, tydat = tydat)
  Hx <- .npcdhat_make_xhat_matrix(
    bws = bws,
    txdat = txdat,
    exdat = exdat
  )
  Gy <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = tydat,
    exdat = eydat,
    operator = rep.int(operator, ncol(tydat))
  )

  (Hx * Gy) %*% rhs
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
