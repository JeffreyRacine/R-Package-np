.npcdhat_make_xbw <- function(bws, txdat) {
  spec <- npConditionalRegEngineSpec(
    bws,
    where = "conditional hat"
  )

  xbw.args <- list(
    xdat = txdat,
    ydat = rep.int(0.0, nrow(txdat)),
    bws = bws$xbw,
    regtype = spec$reg.engine,
    basis = spec$basis.engine,
    degree = spec$degree.engine,
    bernstein.basis = spec$bernstein.engine,
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

.npcdhat_resolve_x_s <- function(bws, txdat, s = NULL, deriv = NULL, where = "npcdhat") {
  if (is.null(s) && is.null(deriv))
    return(NULL)

  txdat <- toFrame(txdat)
  ncon <- sum(bws$ixcon)
  con.names <- names(txdat)[bws$ixcon]
  if (!is.null(s) && !is.null(names(s))) {
    sout <- integer(ncon)
    names(sout) <- con.names
    keep <- intersect(names(s), con.names)
    if (length(keep))
      sout[keep] <- as.integer(s[keep])
    s <- sout
  }
  .npreghat_resolve_s(s = s, deriv = deriv, ncon = ncon, con.names = con.names)
}

.npcdhat_has_x_derivative <- function(s) {
  length(s) > 0L && any(s > 0L)
}

.npcdhat_validate_lp_x_s <- function(degree, s, ncon, where = "npcdhat") {
  if (!.npcdhat_has_x_derivative(s))
    return(invisible(NULL))

  degree <- as.integer(degree)
  s <- as.integer(s)
  ncon <- as.integer(ncon)
  if (length(degree) != ncon || length(s) != ncon) {
    stop(sprintf("%s received incoherent local-polynomial derivative metadata", where),
         call. = FALSE)
  }

  first.derivative.request <- (sum(s) == 1L) && all(s %in% c(0L, 1L))
  lp.degree0.lc.derivative.route <- first.derivative.request &&
    npGlpDegree0FirstDerivativeLcOk(
      regtype.engine = "lp",
      degree.engine = degree,
      gradient.order = rep.int(1L, ncon),
      ncon = ncon
    )

  if (any(s > degree) && !lp.degree0.lc.derivative.route) {
    stop(sprintf("%s requested derivative order in 's' exceeds local polynomial degree", where),
         call. = FALSE)
  }

  invisible(NULL)
}

.npcdhat_use_adaptive_ratio <- function(bws, x.s = NULL) {
  spec <- npConditionalRegEngineSpec(
    bws,
    where = "conditional hat adaptive ratio"
  )
  identical(bws$type, "adaptive_nn") &&
    !.npcdhat_has_x_derivative(x.s) &&
    (identical(spec$reg.engine, "lc") ||
       all(spec$degree.engine == 0L))
}

.npcdhat_make_xhat_matrix <- function(bws,
                                      txdat,
                                      exdat,
                                      s = NULL,
                                      train.is.eval = FALSE) {
  train.is.eval <- npValidateScalarLogical(train.is.eval, "train.is.eval")
  if (train.is.eval && nrow(exdat) != nrow(txdat))
    stop("conditional X hat received inconsistent training-identity rows")
  eval.arg <- if (train.is.eval) NULL else exdat
  xbw <- .npcdhat_make_xbw(bws = bws, txdat = txdat)
  spec <- npConditionalRegEngineSpec(
    xbw,
    where = "conditional hat regression bandwidth",
    ncon.field = "ncon"
  )
  regtype <- spec$reg.engine
  basis <- spec$basis.engine
  degree <- spec$degree.engine
  bernstein <- spec$bernstein.engine

  if (npIsCanonicalLp0(
        regtype.engine = regtype,
        degree.engine = degree,
        ncon = xbw$ncon
      )) {
    if (.npcdhat_has_x_derivative(s)) {
      if (!identical(xbw[["ckertype", exact = TRUE]], "beta")) {
        if (identical(xbw[["ckertype", exact = TRUE]], "uniform"))
          .np_warning("ignoring kernel order specified with uniform kernel type")
        return(.npreghat_exact_lp_matrix_from_regression_core(
          bws = xbw,
          txdat = txdat,
          exdat = eval.arg,
          s = s,
          basis = basis,
          degree = degree,
          bernstein.basis = bernstein
        ))
      }
      return(.npreghat_exact_lc_derivative_matrix_from_npksum_chunked(
        bws = xbw,
        txdat = txdat,
        exdat = eval.arg,
        s = s
      ))
    }

    return(.npreghat_exact_lc_matrix_from_kernel_weights(
      bws = xbw,
      txdat = txdat,
      exdat = eval.arg
    ))
  }

  if (identical(regtype, "lp")) {
    .npcdhat_validate_lp_x_s(
      degree = degree,
      s = s,
      ncon = xbw$ncon,
      where = "npcdhat"
    )

    if (identical(xbw$type, "generalized_nn") &&
        any(degree > 1L) && train.is.eval) {
      return(.npreghat_exact_matrix_from_core(
        bws = xbw,
        txdat = txdat,
        exdat = NULL,
        s = s
      ))
    }

    if (identical(xbw$type, "generalized_nn") && any(degree > 1L)) {
      H <- matrix(NA_real_, nrow = nrow(exdat), ncol = nrow(txdat))
      for (i in seq_len(nrow(exdat))) {
        H[i, ] <- .npreghat_exact_matrix_from_core(
          bws = xbw,
          txdat = txdat,
          exdat = exdat[i, , drop = FALSE],
          s = s
        )[1L, ]
      }
      return(H)
    }

    return(.npreghat_exact_lp_matrix_from_kernel_weights(
      bws = xbw,
      txdat = txdat,
      exdat = eval.arg,
      s = s,
      basis = basis,
      degree = degree,
      bernstein.basis = bernstein
    ))
  }

  .npreghat_exact_matrix_from_core(
    bws = xbw,
    txdat = txdat,
    exdat = eval.arg,
    s = s
  )
}

.np_direct_operator_matrix <- function(kbw,
                                       txdat,
                                       exdat,
                                       operator,
                                       where,
                                       train.is.eval = FALSE) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  train.is.eval <- npValidateScalarLogical(train.is.eval, "train.is.eval")
  if (train.is.eval && nrow(exdat) != nrow(txdat))
    stop(sprintf("%s received inconsistent training-identity rows", where))
  op.info <- .np_operator_kernel_weight_scale(
    bws = kbw,
    operator = operator,
    nvars = ncol(txdat),
    where = where
  )
  kw <- .np_kernel_weights_direct(
    bws = op.info$bws,
    txdat = txdat,
    exdat = if (train.is.eval) NULL else exdat,
    bandwidth.divide = op.info$bandwidth.divide,
    operator = op.info$operator
  )

  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(txdat))
  if (nrow(kw) != nrow(txdat) || ncol(kw) != nrow(exdat))
    stop(sprintf("%s returned unexpected operator shape", where))

  t(kw) / op.info$scale
}

.np_direct_operator_apply <- function(kbw,
                                      txdat,
                                      exdat,
                                      operator,
                                      rhs,
                                      where,
                                      train.is.eval = FALSE) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  train.is.eval <- npValidateScalarLogical(train.is.eval, "train.is.eval")
  if (train.is.eval && nrow(exdat) != nrow(txdat))
    stop(sprintf("%s received inconsistent training-identity rows", where))
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
    exdat = if (train.is.eval) NULL else exdat,
    bandwidth.divide = op.info$bandwidth.divide,
    operator = op.info$operator
  )

  if (!is.matrix(kw))
    kw <- matrix(kw, nrow = nrow(txdat))
  if (nrow(kw) != nrow(txdat) || ncol(kw) != nrow(exdat))
    stop(sprintf("%s returned unexpected operator shape", where))

  crossprod(kw, rhs) / op.info$scale
}

.np_native_output_extent <- function(..., where) {
  factors <- as.double(c(...))
  extent <- prod(factors)

  if (!length(factors) || any(!is.finite(factors)) ||
      any(factors < 0.0) || any(factors != floor(factors)) ||
      !is.finite(extent) || extent > .Machine$integer.max) {
    stop(sprintf("%s exceeds the native output-size capacity", where),
         call. = FALSE)
  }

  as.integer(extent)
}

.np_local_operator_ksum <- function(kbw,
                                    txdat,
                                    exdat,
                                    operator,
                                    rhs,
                                    return.kernel.weights = FALSE,
                                    where,
                                    train.is.eval = FALSE) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)
  train.is.eval <- npValidateScalarLogical(train.is.eval, "train.is.eval")
  if (train.is.eval && nrow(exdat) != nrow(txdat))
    stop(sprintf("%s received inconsistent training-identity rows", where))
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
  tuno <- txm[, bws$iuno, drop = FALSE]
  tcon <- txm[, bws$icon, drop = FALSE]
  tord <- txm[, bws$iord, drop = FALSE]
  if (train.is.eval) {
    euno <- data.frame()
    econ <- data.frame()
    eord <- data.frame()
  } else {
    exm <- toMatrix(exdat)
    euno <- exm[, bws$iuno, drop = FALSE]
    econ <- exm[, bws$icon, drop = FALSE]
    eord <- exm[, bws$iord, drop = FALSE]
  }

  tnrow <- nrow(txdat)
  enrow <- if (train.is.eval) tnrow else nrow(exdat)
  nksum <- .np_native_output_extent(
    enrow,
    max(1L, ncol(rhs)),
    where = sprintf("%s kernel-sum result", where)
  )
  nkw <- if (isTRUE(return.kernel.weights)) {
    .np_native_output_extent(
      tnrow,
      enrow,
      where = sprintf("%s kernel-weight result", where)
    )
  } else {
    0L
  }

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
      beta = CKER_COORDINATE
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
    miss.ex = train.is.eval,
    leave.one.out = FALSE,
    bandwidth.divide = op.info$bandwidth.divide,
    mcv.numRow = attr(bws$xmcv, "num.row"),
    wncol = 0L,
    yncol = ncol(rhs),
    int_do_tree = if (identical(bws[["ckertype", exact = TRUE]], "beta")) {
      DO_TREE_NO
    } else {
      npDoTreeOrCategoricalCompress(
        ncon = bws$ncon,
        ncat = bws$nuno + bws$nord,
        bws = bws
      )
    },
    return.kernel.weights = isTRUE(return.kernel.weights),
    permutation.operator = PERMUTATION_OPERATORS[["none"]],
    compute.score = FALSE,
    compute.ocg = FALSE
  )
  myopti <- c(myopti, npContinuousKernelDescriptorOptions(bws))
  myopti <- c(myopti, list(divide.returned.kernel.weights = FALSE))

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
    nksum,
    as.integer(0L),
    as.integer(nkw),
    as.double(cker.bounds.c$lb),
    as.double(cker.bounds.c$ub),
    PACKAGE = "np"
  )

  out <- myout[["ksum"]]
  if (ncol(rhs) == 1L) {
    out <- matrix(out, nrow = enrow, ncol = 1L)
  } else {
    out <- t(matrix(out, nrow = ncol(rhs), ncol = enrow))
  }

  kw <- NULL
  if (isTRUE(return.kernel.weights))
    kw <- matrix(as.double(myout[["kernel.weights"]]), nrow = tnrow, ncol = enrow)

  list(ksum = out, kernel.weights = kw)
}

.np_exact_operator_apply <- function(kbw,
                                     txdat,
                                     exdat,
                                     operator,
                                     rhs,
                                     where,
                                     train.is.eval = FALSE) {
  rhs <- as.matrix(rhs)
  storage.mode(rhs) <- "double"

  out <- .np_local_operator_ksum(
    kbw = kbw,
    txdat = txdat,
    exdat = exdat,
    operator = operator,
    rhs = rhs,
    return.kernel.weights = FALSE,
    where = where,
    train.is.eval = train.is.eval
  )$ksum

  if (!is.matrix(out))
    out <- matrix(out, nrow = nrow(exdat), ncol = ncol(rhs))

  if (ncol(out) == 1L)
    out[, 1L, drop = FALSE]
  else
    out
}

.np_exact_operator_matrix <- function(kbw,
                                      txdat,
                                      exdat,
                                      operator,
                                      where,
                                      train.is.eval = FALSE) {
  txdat <- toFrame(txdat)
  exdat <- toFrame(exdat)

  probe <- .np_local_operator_ksum(
    kbw = kbw,
    txdat = txdat,
    exdat = exdat,
    operator = operator,
    rhs = matrix(1.0, nrow = nrow(txdat), ncol = 1L),
    return.kernel.weights = TRUE,
    where = where,
    train.is.eval = train.is.eval
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

.npcdhat_make_kernel_matrix <- function(kbw,
                                        txdat,
                                        exdat,
                                        operator,
                                        train.is.eval = FALSE) {
  if (!identical(kbw$type, "fixed")) {
    return(.np_exact_operator_matrix(
      kbw = kbw,
      txdat = txdat,
      exdat = exdat,
      operator = operator,
      where = "conditional hat exact kernel matrix",
      train.is.eval = train.is.eval
    ))
  }

  .np_direct_operator_matrix(
    kbw = kbw,
    txdat = txdat,
    exdat = exdat,
    operator = operator,
    where = "conditional hat direct kernel matrix",
    train.is.eval = train.is.eval
  )
}

.npcdhat_ratio_matrix <- function(bws,
                                  txdat,
                                  tydat,
                                  exdat,
                                  eydat,
                                  operator,
                                  train.is.eval = FALSE) {
  xkbw <- .npcdhat_make_xkbw(bws = bws, txdat = txdat)
  ybw <- .npcdhat_make_ybw(bws = bws, tydat = tydat)
  Kx <- .npcdhat_make_kernel_matrix(
    kbw = xkbw,
    txdat = txdat,
    exdat = exdat,
    operator = rep.int("normal", ncol(txdat)),
    train.is.eval = train.is.eval
  )
  Ky <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = tydat,
    exdat = eydat,
    operator = rep.int(operator, ncol(tydat)),
    train.is.eval = train.is.eval
  )

  denom <- rowSums(Kx) / nrow(txdat)
  if (length(denom) &&
      !anyNA(denom) &&
      min(denom) > -Inf &&
      max(denom) < Inf &&
      !any(denom == 0.0))
    return(sweep((Kx * Ky) / nrow(txdat), 1L, denom, "/"))

  finite.nonzero <- is.finite(denom) & denom != 0.0
  divisor <- pmax(denom, .Machine$double.eps)
  divisor[finite.nonzero] <- denom[finite.nonzero]
  sweep((Kx * Ky) / nrow(txdat), 1L, divisor, "/")
}

.npcdhat_exact_matrix <- function(bws,
                                  txdat,
                                  tydat,
                                  exdat,
                                  eydat,
                                  operator,
                                  x.s = NULL,
                                  train.is.eval = FALSE) {
  if (.npcdhat_use_adaptive_ratio(bws = bws, x.s = x.s)) {
    return(.npcdhat_ratio_matrix(
      bws = bws,
      txdat = txdat,
      tydat = tydat,
      exdat = exdat,
      eydat = eydat,
      operator = operator,
      train.is.eval = train.is.eval
    ))
  }

  ybw <- .npcdhat_make_ybw(bws = bws, tydat = tydat)
  Hx <- .npcdhat_make_xhat_matrix(
    bws = bws,
    txdat = txdat,
    exdat = exdat,
    s = x.s,
    train.is.eval = train.is.eval
  )
  Gy <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = tydat,
    exdat = eydat,
    operator = rep.int(operator, ncol(tydat)),
    train.is.eval = train.is.eval
  )

  Hx * Gy
}

.npcdhat_exact_apply <- function(bws,
                                 txdat,
                                 tydat,
                                 exdat,
                                 eydat,
                                 rhs,
                                 operator,
                                 x.s = NULL,
                                 train.is.eval = FALSE) {
  if (.npcdhat_use_adaptive_ratio(bws = bws, x.s = x.s)) {
    H <- .npcdhat_ratio_matrix(
      bws = bws,
      txdat = txdat,
      tydat = tydat,
      exdat = exdat,
      eydat = eydat,
      operator = operator,
      train.is.eval = train.is.eval
    )
    return(H %*% rhs)
  }

  ybw <- .npcdhat_make_ybw(bws = bws, tydat = tydat)
  Hx <- .npcdhat_make_xhat_matrix(
    bws = bws,
    txdat = txdat,
    exdat = exdat,
    s = x.s,
    train.is.eval = train.is.eval
  )
  Gy <- .npcdhat_make_kernel_matrix(
    kbw = ybw,
    txdat = tydat,
    exdat = eydat,
    operator = rep.int(operator, ncol(tydat)),
    train.is.eval = train.is.eval
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
                          x.deriv = NULL,
                          x.s = NULL,
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

  if (!is.null(x.s) || !is.null(x.deriv)) {
    x.s <- .npcdhat_resolve_x_s(
      bws = bws,
      txdat = txdat,
      s = x.s,
      deriv = x.deriv,
      where = where
    )
  } else {
    x.s <- NULL
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
      operator = operator,
      x.s = x.s,
      train.is.eval = no.exy
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
    operator = operator,
    x.s = x.s,
    train.is.eval = no.exy
  )

  class(H) <- c(class_name, "matrix")
  attr(H, "bws") <- bws
  attr(H, "txdat") <- txdat
  attr(H, "tydat") <- tydat
  attr(H, "exdat") <- exdat
  attr(H, "eydat") <- eydat
  attr(H, "trainiseval") <- no.exy
  attr(H, "rows.omit.train") <- rows.omit.train
  attr(H, "s") <- x.s
  attr(H, "call") <- match.call(expand.dots = FALSE)

  if (!is.null(y)) {
    Hy <- H %*% y
    if (ncol(Hy) == 1L)
      Hy <- as.vector(Hy)
    attr(H, "Hy") <- Hy
  }

  H
}
npConditionalCategoricalFirstDifferences <- function(hat.fun,
                                                      bws,
                                                      txdat,
                                                      tydat,
                                                      exdat = NULL,
                                                      eydat = NULL,
                                                      where) {
  if (!is.function(hat.fun))
    stop(sprintf("%s received an invalid conditional hat evaluator", where),
         call. = FALSE)

  eval.x <- if (is.null(exdat)) txdat else exdat
  eval.y <- if (is.null(eydat)) tydat else eydat
  out <- matrix(NA_real_, nrow = nrow(eval.x), ncol = bws$xndim)
  cat.idx <- which(bws$ixuno | bws$ixord)
  if (!length(cat.idx))
    return(out)

  rhs <- rep.int(1.0, nrow(txdat))
  eval.hat <- function(z) {
    as.vector(hat.fun(
      bws = bws,
      txdat = txdat,
      tydat = tydat,
      exdat = z,
      eydat = eval.y,
      y = rhs,
      output = "apply"
    ))
  }

  for (jj in cat.idx) {
    frames <- npCategoricalFirstDifferenceFrames(
      exdat = eval.x,
      index = jj,
      where = where
    )
    out[, jj] <- eval.hat(frames$upper) - eval.hat(frames$lower)
  }
  out
}
