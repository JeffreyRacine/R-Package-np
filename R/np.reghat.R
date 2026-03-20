npreghat <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npreghat", bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npreghat", bws$call)
        else if (!is.call(bws))
          UseMethod("npreghat", bws)
        else
          UseMethod("npreghat", NULL)
      } else {
        UseMethod("npreghat", NULL)
      }
    } else {
      UseMethod("npreghat", NULL)
    }
  }

.npreghat_resolve_s <- function(s, deriv, ncon, con.names) {
  if (ncon == 0L)
    return(integer(0))

  if (is.null(s)) {
    if (is.null(deriv)) {
      s <- integer(ncon)
    } else {
      s <- deriv
    }
  }

  if (length(s) == 1L)
    s <- c(as.integer(s), rep.int(0L, ncon - 1L))

  if (!is.null(names(s))) {
    sout <- integer(ncon)
    names(sout) <- con.names
    keep <- intersect(names(s), con.names)
    if (length(keep))
      sout[keep] <- as.integer(s[keep])
    s <- sout
  }

  if (length(s) != ncon)
    stop("argument 's' must have length equal to the number of continuous predictors")

  if (anyNA(s) || any(s < 0) || any(s != as.integer(s)))
    stop("argument 's' must be a non-negative integer vector")

  as.integer(s)
}

.npreghat_solve_eval <- function(W, w.eval, k, ridge.base) {
  XtWX <- crossprod(W, W * k)
  p <- nrow(XtWX)
  diag.loc <- cbind(seq_len(p), seq_len(p))
  XtWX.diag <- XtWX[diag.loc]
  ridge.grid <- npRidgeSequenceFromBase(
    n.train = nrow(W),
    ridge.base = max(0.0, as.double(ridge.base)),
    cap = 1.0
  )
  ridge <- ridge.grid[1L]
  solved <- FALSE

  for (ridge.try in ridge.grid) {
    ridge <- ridge.try
    A <- XtWX
    if (ridge > 0)
      A[diag.loc] <- XtWX.diag + ridge

    v <- tryCatch(
      drop(solve(t(A), matrix(w.eval, ncol = 1L))),
      error = function(e) NULL
    )

    if (!is.null(v) && all(is.finite(v))) {
      solved <- TRUE
      break
    }
  }

  if (!solved)
    return(NULL)

  list(v = v, ridge = ridge)
}

.npreghat_exact_matrix_from_core <- function(bws, txdat, exdat = NULL, s = NULL) {
  miss.ex <- is.null(exdat)
  neval <- if (miss.ex) nrow(txdat) else nrow(exdat)
  ntrain <- nrow(txdat)
  H <- matrix(NA_real_, nrow = neval, ncol = ntrain)
  if (is.null(s))
    s <- integer(bws$ncon)

  want.grad <- length(s) > 0L && any(s > 0L)
  target.col <- NULL
  if (want.grad) {
    if (!(sum(s) == 1L && all(s %in% c(0L, 1L))))
      stop("exact direct hat matrix supports only mean and first derivatives")
    target.cont <- which(s == 1L)
    target.col <- which(bws$icon)[target.cont]
  }

  fit_one <- function(ycol) {
    direct.args <- list(
      bws = bws,
      txdat = txdat,
      tydat = ycol,
      exdat = exdat,
      gradients = want.grad,
      gradient.order = 1L
    )
    fit <- .npRmpi_with_local_regression(do.call(.np_regression_direct, direct.args))
    if (is.null(target.col)) {
      fit$mean
    } else {
      fit$grad[, target.col]
    }
  }

  for (j in seq_len(ntrain)) {
    yj <- numeric(ntrain)
    yj[j] <- 1.0
    H[, j] <- fit_one(yj)
  }

  H
}

.npreghat_exact_lc_matrix_from_kernel_weights <- function(bws, txdat, exdat = NULL) {
  kw <- .np_kernel_weights_direct(
    bws = bws,
    txdat = txdat,
    exdat = exdat,
    leave.one.out = FALSE,
    bandwidth.divide = TRUE,
    kernel.pow = 1.0
  )

  denom <- colSums(kw)
  denom[denom == 0.0] <- .Machine$double.xmin
  sweep(t(kw), 1L, denom, "/")
}

.npreghat_exact_lc_derivative_matrix_from_npksum_chunked <- function(bws,
                                                                     txdat,
                                                                     exdat = NULL,
                                                                     s) {
  no.ex <- is.null(exdat)
  txdat <- toFrame(txdat)
  if (!no.ex) {
    exdat <- toFrame(exdat)
    if (!(txdat %~% exdat))
      stop("'txdat' and 'exdat' are not similar data frames!")
  }

  target.cont <- which(s == 1L)
  if (length(target.cont) != 1L)
    stop("lc derivative hat matrix requires exactly one continuous first derivative")

  eval.data <- if (no.ex) txdat else exdat
  ntrain <- nrow(txdat)
  neval <- nrow(eval.data)

  if (identical(bws$type, "adaptive_nn")) {
    out <- npksum.default(
      bws = bws,
      txdat = txdat,
      exdat = if (no.ex) txdat else eval.data,
      bandwidth.divide = TRUE,
      return.kernel.weights = TRUE,
      return.derivative.kernel.weights = TRUE,
      permutation.operator = "derivative"
    )

    kw <- out$kw
    pkw <- out$p.kw
    if (is.null(kw) || is.null(pkw))
      stop("adaptive lc derivative hat matrix requires kernel weights and derivative kernel weights")

    if (!is.matrix(kw))
      kw <- matrix(kw, nrow = ntrain, ncol = neval)

    if (length(dim(pkw)) == 3L) {
      pkw <- pkw[, , target.cont, drop = TRUE]
    } else {
      pkw <- as.matrix(pkw)
    }

    sk <- as.vector(out$ksum)
    dsk <- out$p.ksum
    if (length(dim(dsk)) == 3L) {
      dsk <- dsk[, 1L, target.cont, drop = TRUE]
    } else if (is.matrix(dsk)) {
      dsk <- dsk[, target.cont, drop = TRUE]
    }
    dsk <- as.vector(dsk)

    return(t(
      sweep(pkw, 2L, sk, "/") -
      sweep(kw, 2L, dsk / (sk^2), "*")
    ))
  }

  block.size <- min(512L, ntrain)
  ones <- rep.int(1.0, ntrain)
  H <- matrix(0.0, nrow = neval, ncol = ntrain)

  for (start in seq.int(1L, ntrain, by = block.size)) {
    stop.col <- min(ntrain, start + block.size - 1L)
    cols <- seq.int(start, stop.col)
    ib <- length(cols)
    W <- matrix(0.0, nrow = ntrain, ncol = ib + 1L)
    W[cbind(cols, seq_len(ib))] <- 1.0
    W[, ib + 1L] <- 1.0

    # The derivative call already returns both ksum and p.ksum.
    out <- npksum.default(
      bws = bws,
      txdat = txdat,
      exdat = if (no.ex) txdat else eval.data,
      tydat = ones,
      weights = W,
      bandwidth.divide = TRUE,
      permutation.operator = "derivative"
    )

    ks <- out$ksum
    ps <- out$p.ksum

    ks <- if (is.null(dim(ks))) {
      matrix(ks, nrow = ib + 1L, ncol = neval)
    } else {
      as.matrix(ks)
    }
    ps <- if (is.null(dim(ps))) {
      matrix(ps, nrow = ib + 1L, ncol = neval)
    } else if (length(dim(ps)) == 3L && dim(ps)[3L] == 1L) {
      ps[, , 1L, drop = TRUE]
    } else {
      as.matrix(ps)
    }

    sk <- ks[nrow(ks), ]
    dsk <- ps[nrow(ps), ]

    H[, cols] <- t(
      (ps[seq_len(ib), , drop = FALSE] / rep(sk, each = ib)) -
      (ks[seq_len(ib), , drop = FALSE] * rep(dsk / (sk^2), each = ib))
    )
  }

  H
}

.npreghat_exact_ll_matrix_from_kernel_weights <- function(bws, txdat, exdat = NULL, s = NULL) {
  miss.ex <- is.null(exdat)
  eval.data <- if (miss.ex) txdat else exdat
  ntrain <- nrow(txdat)
  neval <- nrow(eval.data)
  kw <- .np_kernel_weights_direct(
    bws = bws,
    txdat = txdat,
    exdat = if (miss.ex) NULL else eval.data,
    leave.one.out = FALSE,
    bandwidth.divide = identical(bws$type, "adaptive_nn"),
    kernel.pow = 1.0
  )

  xcon.train <- as.matrix(txdat[, bws$icon, drop = FALSE])
  xcon.eval <- as.matrix(eval.data[, bws$icon, drop = FALSE])
  design <- cbind(1.0, xcon.train)
  H <- matrix(NA_real_, nrow = neval, ncol = ntrain)

  want.grad <- length(s) > 0L && any(s > 0L)
  target.cont <- if (want.grad) which(s == 1L) else integer(0)

  for (j in seq_len(neval)) {
    w <- kw[, j]
    A.base <- crossprod(design, design * w)
    rhs.base <- t(design * w)
    solved <- tryCatch(solve(A.base, rhs.base), error = function(e) NULL)

    if (is.null(solved) || !all(is.finite(solved))) {
      eps <- 1.0 / ntrain
      nepsilon <- 0.0
      A.try <- A.base
      rhs.try <- rhs.base
      diag(A.try) <- diag(A.try) + eps
      solved <- tryCatch(solve(A.try, rhs.try), error = function(e) NULL)

      while (is.null(solved) || !all(is.finite(solved))) {
        diag(A.try) <- diag(A.try) + eps
        nepsilon <- nepsilon + eps
        solved <- tryCatch(solve(A.try, rhs.try), error = function(e) NULL)
      }

      if (nepsilon > 0.0) {
        sumw <- sum(w)
        if (sumw == 0.0)
          sumw <- .Machine$double.xmin
        rhs.try[1L, ] <- rhs.try[1L, ] + nepsilon * (w / sumw)
        solved <- solve(A.try, rhs.try)
      }
    }

    if (!want.grad) {
      fit.weights <- c(1.0, xcon.eval[j, , drop = TRUE])
      H[j, ] <- drop(fit.weights %*% solved)
    } else {
      H[j, ] <- solved[1L + target.cont, ]
    }
  }

  H
}

.npreghat_exact_lp_matrix_from_kernel_weights <- function(bws,
                                                          txdat,
                                                          exdat = NULL,
                                                          s = NULL,
                                                          basis = "glp",
                                                          degree = integer(0),
                                                          bernstein.basis = FALSE) {
  miss.ex <- is.null(exdat)
  eval.data <- if (miss.ex) txdat else exdat
  ntrain <- nrow(txdat)
  neval <- nrow(eval.data)
  want.grad <- length(s) > 0L && any(s > 0L)

  if (identical(bws$type, "generalized_nn") &&
      any(degree > 1L) &&
      isTRUE(getOption("np.tree"))) {
    return(.npreghat_exact_matrix_from_core(
      bws = bws,
      txdat = txdat,
      exdat = if (miss.ex) NULL else exdat,
      s = s
    ))
  }

  if (identical(bws$type, "generalized_nn") &&
      any(degree > 1L) &&
      !isTRUE(getOption("np.tree"))) {
    return(.npreghat_exact_lp_matrix_from_regression_core_chunked(
      bws = bws,
      txdat = txdat,
      exdat = if (miss.ex) NULL else exdat,
      basis = basis,
      degree = degree,
      bernstein.basis = bernstein.basis,
      s = s
    ))
  }

  kw <- .np_kernel_weights_direct(
    bws = bws,
    txdat = txdat,
    exdat = if (miss.ex) NULL else eval.data,
    leave.one.out = FALSE,
    bandwidth.divide = TRUE,
    kernel.pow = 1.0
  )

  W.train <- W.lp(
    xdat = txdat[, bws$icon, drop = FALSE],
    degree = degree,
    basis = basis,
    bernstein.basis = bernstein.basis
  )
  W.eval <- W.lp(
    xdat = txdat[, bws$icon, drop = FALSE],
    exdat = if (miss.ex) NULL else eval.data[, bws$icon, drop = FALSE],
    degree = degree,
    gradient.vec = if (want.grad) s else NULL,
    basis = basis,
    bernstein.basis = bernstein.basis
  )

  H <- matrix(NA_real_, nrow = neval, ncol = ntrain)
  eps <- 1.0 / max(1L, ntrain)

  for (j in seq_len(neval)) {
    w <- kw[, j]
    A.base <- crossprod(W.train, W.train * w)
    rhs <- W.eval[j, ]
    solved <- tryCatch(solve(A.base, rhs), error = function(e) NULL)

    if (is.null(solved) || !all(is.finite(solved))) {
      A.try <- A.base
      nepsilon <- 0.0

      repeat {
        diag(A.try) <- diag(A.try) + eps
        nepsilon <- nepsilon + eps
        solved <- tryCatch(solve(A.try, rhs), error = function(e) NULL)
        if (!is.null(solved) && all(is.finite(solved)))
          break
      }

      denom <- A.try[1L, 1L]
      if (!is.finite(denom) || abs(denom) < .Machine$double.xmin)
        denom <- .Machine$double.xmin
      solved[1L] <- solved[1L] * (1.0 + nepsilon / denom)
    }

    H[j, ] <- w * drop(W.train %*% solved)
  }

  H
}

.npreghat_exact_lp_matrix_from_regression_core_chunked <- function(bws,
                                                                   txdat,
                                                                   exdat = NULL,
                                                                   basis = "glp",
                                                                   degree = integer(0),
                                                                   bernstein.basis = FALSE,
                                                                   s = NULL,
                                                                   chunk.size = 128L) {
  no.ex <- is.null(exdat)

  txdat <- toFrame(txdat)
  if (!no.ex) {
    exdat <- toFrame(exdat)
    if (!(txdat %~% exdat))
      stop("'txdat' and 'exdat' are not similar data frames!")
  }

  if (length(bws$bw) != length(txdat))
    stop("length of bandwidth vector does not match number of columns of 'txdat'")

  if ((any(bws$icon) &&
       !all(vapply(txdat[, bws$icon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
      (any(bws$iord) &&
       !all(vapply(txdat[, bws$iord, drop = FALSE], inherits, logical(1), "ordered"))) ||
      (any(bws$iuno) &&
       !all(vapply(txdat[, bws$iuno, drop = FALSE], inherits, logical(1), "factor")))) {
    stop("supplied bandwidths do not match 'txdat' in type")
  }

  txdat <- adjustLevels(txdat, bws$xdati)
  if (!no.ex) {
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    npKernelBoundsCheckEval(exdat, bws$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")
  }

  txmat <- toMatrix(txdat)
  tuno <- txmat[, bws$iuno, drop = FALSE]
  tcon <- txmat[, bws$icon, drop = FALSE]
  tord <- txmat[, bws$iord, drop = FALSE]
  reg.c <- npRegtypeToC(regtype = "lp",
                        degree = degree,
                        ncon = bws$ncon,
                        context = ".npreghat_exact_lp_matrix_from_regression_core_chunked")

  npCheckRegressionDesignCondition(reg.code = reg.c$code,
                                   xcon = tcon,
                                   basis = basis,
                                   degree = degree,
                                   bernstein.basis = bernstein.basis,
                                   where = ".npreghat_exact_lp_matrix_from_regression_core_chunked")

  if (!no.ex) {
    exmat <- toMatrix(exdat)
    euno <- exmat[, bws$iuno, drop = FALSE]
    econ <- exmat[, bws$icon, drop = FALSE]
    eord <- exmat[, bws$iord, drop = FALSE]
  } else {
    euno <- tuno
    econ <- tcon
    eord <- tord
  }

  tuno <- if (bws$nuno > 0L) as.matrix(tuno) else matrix(double(), nrow = nrow(txdat), ncol = 0L)
  tcon <- if (bws$ncon > 0L) as.matrix(tcon) else matrix(double(), nrow = nrow(txdat), ncol = 0L)
  tord <- if (bws$nord > 0L) as.matrix(tord) else matrix(double(), nrow = nrow(txdat), ncol = 0L)
  euno <- if (bws$nuno > 0L) as.matrix(euno) else matrix(double(), nrow = if (no.ex) nrow(txdat) else nrow(exdat), ncol = 0L)
  econ <- if (bws$ncon > 0L) as.matrix(econ) else matrix(double(), nrow = if (no.ex) nrow(txdat) else nrow(exdat), ncol = 0L)
  eord <- if (bws$nord > 0L) as.matrix(eord) else matrix(double(), nrow = if (no.ex) nrow(txdat) else nrow(exdat), ncol = 0L)

  bwtype.c <- switch(bws$type,
    fixed = BW_FIXED,
    generalized_nn = BW_GEN_NN,
    adaptive_nn = BW_ADAP_NN
  )
  kernel.x.c <- switch(bws$ckertype,
    gaussian = CKER_GAUSS + bws$ckerorder / 2 - 1,
    epanechnikov = CKER_EPAN + bws$ckerorder / 2 - 1,
    uniform = CKER_UNI,
    "truncated gaussian" = CKER_TGAUSS
  )
  kernel.xu.c <- switch(bws$ukertype,
    aitchisonaitken = UKER_AIT,
    liracine = UKER_LR
  )
  kernel.xo.c <- switch(bws$okertype,
    wangvanryzin = OKER_WANG,
    liracine = OKER_LR,
    racineliyan = OKER_RLY
  )

  ntrain <- nrow(txdat)
  neval <- if (no.ex) ntrain else nrow(exdat)
  H <- matrix(0.0, nrow = neval, ncol = ntrain)
  chunk.size <- max(1L, min(as.integer(chunk.size), ntrain))
  bw.vec <- as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord]))
  tree.flag <- isTRUE(getOption("np.tree"))
  grad.vec <- if (length(s) && any(s > 0L)) as.integer(s) else integer(0L)
  on.exit(
    tryCatch(.Call("C_np_shadow_reset_state", PACKAGE = "npRmpi"),
             error = function(e) NULL),
    add = TRUE
  )

  for (start in seq.int(1L, ntrain, by = chunk.size)) {
    stop.col <- min(ntrain, start + chunk.size - 1L)
    cols <- seq.int(start, stop.col)
    rhs <- matrix(0.0, nrow = ntrain, ncol = length(cols))
    rhs[cbind(cols, seq_along(cols))] <- 1.0
    H[, cols] <- .Call(
      "C_np_regression_lp_apply_conditional",
      tuno,
      tord,
      tcon,
      euno,
      eord,
      econ,
      rhs,
      bw.vec,
      as.integer(bwtype.c),
      as.integer(kernel.x.c),
        as.integer(kernel.xu.c),
        as.integer(kernel.xo.c),
        as.logical(tree.flag),
        as.integer(degree),
        grad.vec,
        as.integer(isTRUE(bernstein.basis)),
        as.integer(npLpBasisCode(basis)),
        PACKAGE = "npRmpi"
      )
  }

  H
}

.np_kernel_weights_direct <- function(bws,
                                      txdat,
                                      exdat = NULL,
                                      leave.one.out = FALSE,
                                      bandwidth.divide = TRUE,
                                      kernel.pow = 1.0,
                                      operator = NULL) {
  miss.ex <- is.null(exdat)
  txdat <- toFrame(txdat)
  if (!miss.ex) {
    exdat <- toFrame(exdat)
    if (!(txdat %~% exdat))
      stop("'txdat' and 'exdat' are not similar data frames!")
  }

  if (!isa(bws, "kbandwidth"))
    bws <- kbandwidth(bws)

  if (length(bws$bw) != length(txdat))
    stop("length of bandwidth vector does not match number of columns of 'txdat'")

  leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")
  bandwidth.divide <- npValidateScalarLogical(bandwidth.divide, "bandwidth.divide")
  if (!miss.ex && leave.one.out)
    stop("you may not specify 'leave.one.out = TRUE' and provide evaluation data")

  if (is.null(operator)) {
    operator <- rep.int("normal", length(txdat))
  } else {
    if (length(operator) == 1L)
      operator <- rep.int(operator, length(txdat))
    if (length(operator) != length(txdat))
      stop("operator not specified for all variables")
  }

  uo.operators <- c("normal", "convolution", "integral")
  if (!all(operator %in% names(ALL_OPERATORS)))
    stop("invalid operator specification")
  if (!all(operator[bws$iuno | bws$iord] %in% uo.operators))
    stop("unordered and ordered variables may only use normal, convolution, or integral operators")
  operator.num <- ALL_OPERATORS[operator]

  txdat <- adjustLevels(txdat, bws$xdati, allowNewCells = TRUE)
  if (!miss.ex) {
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    npKernelBoundsCheckEval(exdat, bws$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")
  }

  txm <- toMatrix(txdat)
  tuno <- txm[, bws$iuno, drop = FALSE]
  tcon <- txm[, bws$icon, drop = FALSE]
  tord <- txm[, bws$iord, drop = FALSE]

  if (!miss.ex) {
    exm <- toMatrix(exdat)
    euno <- exm[, bws$iuno, drop = FALSE]
    econ <- exm[, bws$icon, drop = FALSE]
    eord <- exm[, bws$iord, drop = FALSE]
  } else {
    euno <- data.frame()
    eord <- data.frame()
    econ <- data.frame()
  }

  tnrow <- nrow(txdat)
  enrow <- if (miss.ex) tnrow else nrow(exdat)
  nkw <- tnrow * enrow

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
    miss.ex = miss.ex,
    leave.one.out = leave.one.out,
    bandwidth.divide = bandwidth.divide,
    mcv.numRow = attr(bws$xmcv, "num.row"),
    wncol = 0L,
    yncol = 0L,
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    return.kernel.weights = TRUE,
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
    as.double(0.0), as.double(0.0),
    asDouble(euno), asDouble(eord), asDouble(econ),
    as.double(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
    as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
    as.integer(c(operator.num[bws$icon], operator.num[bws$iuno], operator.num[bws$iord])),
    as.integer(myopti), as.double(kernel.pow),
    as.integer(enrow),
    as.integer(0L),
    as.integer(nkw),
    as.double(cker.bounds.c$lb),
    as.double(cker.bounds.c$ub),
    PACKAGE = "npRmpi"
  )

  kw <- matrix(as.double(myout[["kernel.weights"]]), nrow = tnrow, ncol = enrow)
  if (leave.one.out && miss.ex && nrow(kw) == ncol(kw))
    diag(kw) <- 0.0
  kw
}

.np_regression_lc_mean_from_kernel_weights <- function(bws,
                                                       txdat,
                                                       tydat,
                                                       exdat = NULL) {
  kw <- .np_kernel_weights_direct(
    bws = bws,
    txdat = txdat,
    exdat = exdat,
    leave.one.out = FALSE,
    bandwidth.divide = identical(bws$type, "adaptive_nn"),
    kernel.pow = 1.0
  )
  kw.sum <- colSums(kw)
  kw.sum[kw.sum == 0.0] <- .Machine$double.xmin
  drop(crossprod(kw, as.double(tydat)) / kw.sum)
}

.np_regression_lc_gradient_from_kernel_weights <- function(bws,
                                                           txdat,
                                                           tydat,
                                                           exdat = NULL) {
  no.ex <- is.null(exdat)
  txdat <- toFrame(txdat)
  if (!no.ex) {
    exdat <- toFrame(exdat)
    if (!(txdat %~% exdat))
      stop("'txdat' and 'exdat' are not similar data frames!")
  }

  if (!(is.vector(tydat) || is.factor(tydat)))
    stop("'tydat' must be a vector or a factor")
  if (nrow(txdat) != length(tydat))
    stop("number of explanatory data 'txdat' and dependent data 'tydat' do not match")

  if (is.factor(tydat)) {
    tydat <- adjustLevels(data.frame(tydat), bws$ydati)[, 1L]
    tydat <- (bws$ydati$all.dlev[[1L]])[as.integer(tydat)]
  } else {
    tydat <- as.double(tydat)
  }

  eval.data <- if (no.ex) txdat else exdat
  kw.obj <- npksum.default(
    bws = bws,
    txdat = txdat,
    exdat = eval.data,
    tydat = rep.int(1.0, nrow(txdat)),
    weights = cbind(tydat, 1.0),
    bandwidth.divide = TRUE,
    permutation.operator = "derivative"
  )

  ks <- kw.obj$ksum
  ps <- kw.obj$p.ksum

  if (is.null(dim(ks)))
    ks <- matrix(ks, nrow = 2L)
  else
    ks <- as.matrix(ks)

  if (is.null(dim(ps))) {
    ps <- array(ps, dim = c(2L, nrow(eval.data), 1L))
  } else if (length(dim(ps)) == 2L) {
    ps <- array(ps, dim = c(dim(ps), 1L))
  }

  denom <- ks[2L, ]
  denom[denom == 0.0] <- .Machine$double.xmin
  numer <- ks[1L, ]

  grad <- matrix(0.0, nrow = nrow(eval.data), ncol = bws$ncon)
  for (j in seq_len(bws$ncon)) {
    ps.j <- ps[, , j, drop = TRUE]
    grad[, j] <- (ps.j[1L, ] / denom) - (numer * ps.j[2L, ] / (denom^2))
  }

  grad
}

.np_regression_direct <- function(bws,
                                  txdat,
                                  tydat,
                                  exdat = NULL,
                                  gradients = FALSE,
                                  gradient.order = 1L,
                                  local.mode = FALSE) {
  no.ex <- is.null(exdat)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  local.mode <- npValidateScalarLogical(local.mode, "local.mode")

  txdat <- toFrame(txdat)
  if (!no.ex) {
    exdat <- toFrame(exdat)
    if (!(txdat %~% exdat))
      stop("'txdat' and 'exdat' are not similar data frames!")
  }

  if (length(bws$bw) != length(txdat))
    stop("length of bandwidth vector does not match number of columns of 'txdat'")

  regtype <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  basis <- npValidateLpBasis(regtype = regtype, basis = bws$basis)
  degree <- npValidateGlpDegree(regtype = regtype,
                                degree = bws$degree,
                                ncon = bws$ncon)
  bernstein.basis <- npValidateGlpBernstein(regtype = regtype,
                                            bernstein.basis = bws$bernstein.basis)
  reg.spec <- npCanonicalConditionalRegSpec(
    regtype = regtype,
    basis = basis,
    degree = degree,
    bernstein.basis = bernstein.basis,
    ncon = bws$ncon,
    where = ".np_regression_direct"
  )
  glp.gradient.order <- if (identical(reg.spec$regtype.engine, "lp")) {
    if (identical(regtype, "lp")) {
      npValidateGlpGradientOrder(regtype = regtype,
                                 gradient.order = gradient.order,
                                 ncon = bws$ncon)
    } else if (bws$ncon > 0L) {
      rep.int(1L, bws$ncon)
    } else {
      integer(0)
    }
  } else {
    NULL
  }
  if (isTRUE(gradients) &&
      identical(reg.spec$regtype.engine, "lp") &&
      (bws$ncon > 0L) &&
      all(reg.spec$degree.engine == 0L)) {
    stop("regtype='lp' with degree=0 does not support derivatives; use gradients=FALSE for fitted/predicted values")
  }

  reg.c <- npRegtypeToC(regtype = reg.spec$regtype.engine,
                        degree = reg.spec$degree.engine,
                        ncon = bws$ncon,
                        context = ".np_regression_direct")
  degree.c <- if (bws$ncon > 0L) {
    as.integer(if (is.null(reg.c$degree)) rep.int(0L, bws$ncon) else reg.c$degree)
  } else {
    integer(1L)
  }

  if ((any(bws$icon) &&
       !all(vapply(txdat[, bws$icon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
      (any(bws$iord) &&
       !all(vapply(txdat[, bws$iord, drop = FALSE], inherits, logical(1), "ordered"))) ||
      (any(bws$iuno) &&
       !all(vapply(txdat[, bws$iuno, drop = FALSE], inherits, logical(1), "factor")))) {
    stop("supplied bandwidths do not match 'txdat' in type")
  }

  if (!(is.vector(tydat) || is.factor(tydat)))
    stop("'tydat' must be a vector or a factor")
  if (nrow(txdat) != length(tydat))
    stop("number of explanatory data 'txdat' and dependent data 'tydat' do not match")

  if (is.factor(tydat)) {
    tydat <- adjustLevels(data.frame(tydat), bws$ydati)[, 1L]
    tydat <- (bws$ydati$all.dlev[[1L]])[as.integer(tydat)]
  } else {
    tydat <- as.double(tydat)
  }

  mean.override <- !isTRUE(gradients) &&
    identical(regtype, "lc") &&
    identical(bws$type, "adaptive_nn")
  grad.override <- isTRUE(gradients) &&
    identical(regtype, "lc") &&
    identical(bws$type, "adaptive_nn") &&
    (bws$ncon > 0L)
  txdat.frame <- txdat
  exdat.frame <- if (no.ex) NULL else exdat

  txdat <- adjustLevels(txdat, bws$xdati)
  if (!no.ex) {
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    npKernelBoundsCheckEval(exdat, bws$icon, bws$ckerlb, bws$ckerub, argprefix = "cker")
  }

  txmat <- toMatrix(txdat)
  tuno <- txmat[, bws$iuno, drop = FALSE]
  tcon <- txmat[, bws$icon, drop = FALSE]
  tord <- txmat[, bws$iord, drop = FALSE]

  npCheckRegressionDesignCondition(reg.code = reg.c$code,
                                   xcon = tcon,
                                   basis = reg.spec$basis.engine,
                                   degree = reg.spec$degree.engine,
                                   bernstein.basis = reg.spec$bernstein.basis.engine,
                                   where = ".np_regression_direct")

  if (!no.ex) {
    exmat <- toMatrix(exdat)
    euno <- exmat[, bws$iuno, drop = FALSE]
    econ <- exmat[, bws$icon, drop = FALSE]
    eord <- exmat[, bws$iord, drop = FALSE]
  } else {
    euno <- data.frame()
    econ <- data.frame()
    eord <- data.frame()
  }

  tnrow <- nrow(txdat)
  enrow <- if (no.ex) tnrow else nrow(exdat)
  ncol.x <- ncol(txdat)

  myopti <- list(
    num_obs_train = tnrow,
    num_obs_eval = enrow,
    num_uno = bws$nuno,
    num_ord = bws$nord,
    num_con = bws$ncon,
    int_LARGE_SF = if (bws$scaling) SF_NORMAL else SF_ARB,
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
      "racineliyan" = OKER_RLY
    ),
    ey_is_ty = TRUE,
    do_grad = gradients,
    regtype = reg.c$code,
    no.ex = no.ex,
    mcv.numRow = attr(bws$xmcv, "num.row"),
    int_do_tree = if (isTRUE(getOption("np.tree"))) DO_TREE_YES else DO_TREE_NO,
    old.reg = FALSE
  )

  cker.bounds.c <- npKernelBoundsMarshal(bws$ckerlb[bws$icon], bws$ckerub[bws$icon])

  asDouble <- function(data) {
    if (is.null(data)) as.double(0.0) else as.double(data)
  }

  glp.gradient.order.c <- if (bws$ncon > 0L) {
    as.integer(if (is.null(glp.gradient.order)) rep.int(1L, bws$ncon) else glp.gradient.order)
  } else {
    integer(1L)
  }

  if (isTRUE(local.mode) && .npRmpi_has_active_slave_pool(comm = 1L)) {
    old.mode <- .Call("C_np_set_local_regression_mode", TRUE, PACKAGE = "npRmpi")
    on.exit(.Call("C_np_set_local_regression_mode", old.mode, PACKAGE = "npRmpi"), add = TRUE)
  }

  myout <- .Call(
    "C_np_regression",
    asDouble(tuno), asDouble(tord), asDouble(tcon), as.double(tydat),
    asDouble(euno), asDouble(eord), asDouble(econ), as.double(double()),
    asDouble(c(bws$bw[bws$icon], bws$bw[bws$iuno], bws$bw[bws$iord])),
    asDouble(bws$xmcv), asDouble(attr(bws$xmcv, "pad.num")),
    asDouble(bws$nconfac), asDouble(bws$ncatfac), asDouble(bws$sdev),
    as.integer(myopti),
    as.integer(degree.c),
    as.integer(glp.gradient.order.c),
    as.integer(isTRUE(reg.spec$bernstein.basis.engine)),
    as.integer(npLpBasisCode(reg.spec$basis.engine)),
    as.integer(enrow),
    as.integer(ncol.x),
    as.logical(gradients),
    as.double(cker.bounds.c$lb),
    as.double(cker.bounds.c$ub),
    PACKAGE = "npRmpi"
  )

  mean.out <- as.double(myout$mean)
  if (mean.override) {
    mean.out <- .np_regression_lc_mean_from_kernel_weights(
      bws = bws,
      txdat = txdat.frame,
      tydat = tydat,
      exdat = exdat.frame
    )
  }

  out <- list(mean = mean.out)

  if (gradients) {
    grad <- matrix(data = myout$g, nrow = enrow, ncol = ncol.x, byrow = FALSE)
    rorder <- numeric(ncol.x)
    ord.idx <- seq_len(ncol.x)
    rorder[c(ord.idx[bws$icon], ord.idx[bws$iuno], ord.idx[bws$iord])] <- ord.idx
    grad <- as.matrix(grad[, rorder, drop = FALSE])

    if (identical(regtype, "lp")) {
      cont.idx <- which(bws$icon)
      if (length(cont.idx)) {
        invalid.order <- glp.gradient.order > degree
        if (any(invalid.order))
          grad[, cont.idx[invalid.order]] <- NA_real_
      }
    }

    if (grad.override) {
      grad[, bws$icon] <- .np_regression_lc_gradient_from_kernel_weights(
        bws = bws,
        txdat = txdat.frame,
        tydat = tydat,
        exdat = exdat.frame
      )
    }

    out$grad <- grad
  }

  out
}

npreghat.formula <-
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1, m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    mf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    y <- model.response(mf)
    txdat <- mf[, attr(attr(mf, "terms"), "term.labels"), drop = FALSE]

    has.eval <- !is.null(newdata)
    if (has.eval) {
      et <- delete.response(tt)
      emf <- model.frame(et, data = newdata)
      exdat <- emf[, attr(attr(emf, "terms"), "term.labels"), drop = FALSE]
    }

    hat.args <- list(bws = bws, txdat = txdat)
    if (!is.null(y))
      hat.args$y <- y
    if (has.eval)
      hat.args$exdat <- exdat

    ev <- do.call(npreghat, c(hat.args, list(...)))
    attr(ev, "call") <- match.call(expand.dots = FALSE)
    ev
  }

npreghat.call <-
  function(bws, ...) {
    ev <- npreghat(
      txdat = .np_eval_bws_call_arg(bws, "xdat"),
      y = .np_eval_bws_call_arg(bws, "ydat"),
      bws = bws,
      ...
    )
    attr(ev, "call") <- match.call(expand.dots = FALSE)
    ev
  }

npreghat.npregression <-
  function(bws, txdat, y, ...){
    model <- bws
    rbw <- model$bws

    if (missing(txdat))
      txdat <- .np_eval_bws_call_arg(rbw, "xdat")

    if (missing(y))
      y <- .np_eval_bws_call_arg(rbw, "ydat")

    npreghat(bws = rbw, txdat = txdat, y = y, ...)
  }

npreghat.rbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           exdat,
           y = NULL,
           output = c("matrix", "apply"),
           basis = NULL,
           bernstein.basis = NULL,
           degree = NULL,
           deriv = NULL,
           leave.one.out = FALSE,
           ridge = 0.0,
           s = NULL,
           ...){

    no.ex <- missing(exdat)
    regtype.raw <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
    basis.raw <- npValidateLpBasis(regtype = regtype.raw, basis = bws$basis)
    degree.raw <- npValidateGlpDegree(regtype = regtype.raw,
                                      degree = bws$degree,
                                      ncon = bws$ncon)
    bernstein.raw <- npValidateGlpBernstein(regtype = regtype.raw,
                                            bernstein.basis = bws$bernstein.basis)
    reg.spec.raw <- npCanonicalConditionalRegSpec(
      regtype = regtype.raw,
      basis = basis.raw,
      degree = degree.raw,
      bernstein.basis = bernstein.raw,
      ncon = bws$ncon,
      where = "npreghat"
    )
    world.size <- .npRmpi_safe_int(mpi.comm.size(0))
    world.size <- if (is.na(world.size)) 1L else as.integer(world.size)
    use.spawn.local.exec <- .npRmpi_autodispatch_active() &&
      !.npRmpi_autodispatch_in_context() &&
      !.npRmpi_autodispatch_called_from_bcast() &&
      (world.size <= 1L) &&
      identical(bws$type, "generalized_nn") &&
      identical(reg.spec.raw$regtype.engine, "lp") &&
      identical(reg.spec.raw$basis.engine, "glp") &&
      !isTRUE(reg.spec.raw$bernstein.basis.engine) &&
      (bws$ncon > 0L) &&
      all(reg.spec.raw$degree.engine == 1L)
    if (use.spawn.local.exec) {
      old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
      old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
      options(npRmpi.autodispatch.context = TRUE)
      options(npRmpi.autodispatch.disable = TRUE)
      on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
      on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
      if (no.ex) {
        return(get("npreghat.rbandwidth", envir = asNamespace("npRmpi"), inherits = FALSE)(
          bws = .npRmpi_autodispatch_untag(bws),
          txdat = txdat,
          y = y,
          output = output,
          s = s,
          deriv = deriv,
          degree = degree,
          basis = basis,
          bernstein.basis = bernstein.basis,
          ridge = ridge,
          leave.one.out = leave.one.out
        ))
      }
      return(get("npreghat.rbandwidth", envir = asNamespace("npRmpi"), inherits = FALSE)(
        bws = .npRmpi_autodispatch_untag(bws),
        txdat = txdat,
        exdat = exdat,
        y = y,
        output = output,
        s = s,
        deriv = deriv,
        degree = degree,
        basis = basis,
        bernstein.basis = bernstein.basis,
        ridge = ridge,
        leave.one.out = leave.one.out
      ))
    }
    if (.npRmpi_has_active_slave_pool(comm = 1L) &&
        !.npRmpi_autodispatch_in_context() &&
        !.npRmpi_autodispatch_called_from_bcast()) {
      expr <- substitute({
        old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
        old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
        options(npRmpi.autodispatch.context = TRUE)
        options(npRmpi.autodispatch.disable = TRUE)
        on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
        on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
        if (NOEX) {
          get("npreghat.rbandwidth", envir = asNamespace("npRmpi"), inherits = FALSE)(
            bws = BWS,
            txdat = TXDAT,
            y = YDAT,
            output = OUTPUT,
            s = SVAL,
            deriv = DERIV,
            degree = DEGREE,
            basis = BASIS,
            bernstein.basis = BERN,
            ridge = RIDGE,
            leave.one.out = LOO
          )
        } else {
          get("npreghat.rbandwidth", envir = asNamespace("npRmpi"), inherits = FALSE)(
            bws = BWS,
            txdat = TXDAT,
            exdat = EXDAT,
            y = YDAT,
            output = OUTPUT,
            s = SVAL,
            deriv = DERIV,
            degree = DEGREE,
            basis = BASIS,
            bernstein.basis = BERN,
            ridge = RIDGE,
            leave.one.out = LOO
          )
        }
      }, list(
        NOEX = no.ex,
        BWS = bws,
        TXDAT = txdat,
        EXDAT = if (no.ex) NULL else exdat,
        YDAT = y,
        OUTPUT = output,
        SVAL = s,
        DERIV = deriv,
        DEGREE = degree,
        BASIS = basis,
        BERN = bernstein.basis,
        RIDGE = ridge,
        LOO = leave.one.out
      ))
      return(.npRmpi_bcast_cmd_expr(expr, comm = 1L, caller.execute = TRUE))
    }

    output <- match.arg(output)
    dots <- list(...)
    npRejectLegacyLpArgs(names(dots), where = "npreghat")

    txdat <- toFrame(txdat)
    leave.one.out <- npValidateScalarLogical(leave.one.out, "leave.one.out")

    if (!no.ex) {
      exdat <- toFrame(exdat)
      if (!(txdat %~% exdat))
        stop("'txdat' and 'exdat' are not similar data frames!")
    }
    if (!no.ex && leave.one.out)
      stop("you may not specify 'leave.one.out = TRUE' and provide evaluation data")

    if (length(bws$bw) != ncol(txdat))
      stop("length of bandwidth vector does not match number of columns of 'txdat'")

    npValidateRegressionNnLowerBound(bws, where = "npreghat")

    if ((any(bws$icon) &&
         !all(vapply(txdat[, bws$icon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$iord) &&
         !all(vapply(txdat[, bws$iord, drop = FALSE], inherits, logical(1), "ordered"))) ||
        (any(bws$iuno) &&
         !all(vapply(txdat[, bws$iuno, drop = FALSE], inherits, logical(1), "factor"))))
      stop("supplied bandwidths do not match 'txdat' in type")

    if (!is.null(y)) {
      if (is.factor(y) || is.vector(y)) {
        if (length(y) != nrow(txdat))
          stop("length of 'y' must match the number of training rows in 'txdat'")
        y <- matrix(as.double(y), ncol = 1L)
      } else {
        y <- as.matrix(y)
        if (nrow(y) != nrow(txdat))
          stop("number of rows in 'y' must match the number of training rows in 'txdat'")
      }
    }

    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(txdat), "na.action")
    if (!is.null(y))
      rows.omit <- union(rows.omit, attr(na.omit(as.data.frame(y)), "na.action"))

    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows, , drop = FALSE]
    if (!is.null(y))
      y <- y[keep.rows, , drop = FALSE]

    if (!no.ex) {
      keep.eval <- rep_len(TRUE, nrow(exdat))
      rows.omit.eval <- attr(na.omit(exdat), "na.action")
      if (length(rows.omit.eval) > 0L)
        keep.eval[as.integer(rows.omit.eval)] <- FALSE
      if (!any(keep.eval))
        stop("Evaluation data has no rows without NAs")
      exdat <- exdat[keep.eval, , drop = FALSE]
    }

    regtype <- if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
    ncon <- bws$ncon
    con.names <- names(txdat)[which(bws$icon)]

    degree <- if (identical(regtype, "lc")) {
      rep.int(0L, ncon)
    } else if (identical(regtype, "ll")) {
      rep.int(1L, ncon)
    } else {
      npValidateGlpDegree(regtype = "lp",
                          degree = if (is.null(degree)) bws$degree else degree,
                          ncon = ncon)
    }

    basis <- npValidateLpBasis(
      regtype = "lp",
      basis = if (is.null(basis)) {
        if (is.null(bws$basis)) "glp" else bws$basis
      } else {
        basis
      }
    )

    bernstein.basis <- npValidateGlpBernstein(
      regtype = "lp",
      bernstein.basis = if (is.null(bernstein.basis)) {
        isTRUE(bws$bernstein.basis)
      } else {
        bernstein.basis
      }
    )

    s <- .npreghat_resolve_s(s = s, deriv = deriv, ncon = ncon, con.names = con.names)
    reg.spec <- npCanonicalConditionalRegSpec(
      regtype = regtype,
      basis = basis,
      degree = degree,
      bernstein.basis = bernstein.basis,
      ncon = ncon,
      where = "npreghat"
    )

    first.derivative.request <- (sum(s) == 1L) && all(s %in% c(0L, 1L))
    simple.operator.request <- (sum(s) == 0L) || first.derivative.request

    direct.apply <- identical(output, "apply") &&
      !is.null(y) &&
      !isTRUE(leave.one.out) &&
      (ncol(y) == 1L) &&
      simple.operator.request

    lc.derivative.exact.route <- identical(regtype, "lc") &&
      first.derivative.request

    exact.lc.kernel.route <- !isTRUE(leave.one.out) &&
      !any(s > 0L) &&
      identical(regtype, "lc") &&
      (bws$ncon > 0L)

    exact.ll.kernel.route <- !isTRUE(leave.one.out) &&
      simple.operator.request &&
      identical(regtype, "ll") &&
      (bws$ncon > 0L)

    exact.lp.kernel.route <- !isTRUE(leave.one.out) &&
      simple.operator.request &&
      identical(regtype, "lp") &&
      identical(reg.spec$regtype.engine, "lp") &&
      (bws$ncon > 0L)

    exact.core.route <- !isTRUE(leave.one.out) &&
      simple.operator.request &&
      (
        exact.lc.kernel.route ||
        exact.ll.kernel.route ||
        exact.lp.kernel.route ||
        lc.derivative.exact.route ||
        FALSE
      )

    if (direct.apply) {
      direct.args <- list(
        bws = bws,
        txdat = txdat,
        tydat = as.vector(y[, 1L]),
        exdat = if (no.ex) NULL else exdat,
        gradients = any(s > 0L),
        gradient.order = 1L
      )
      direct.out <- .npRmpi_with_local_regression(do.call(.np_regression_direct, direct.args))

      if (!any(s > 0L))
        return(as.vector(direct.out$mean))

      target.cont <- which(s == 1L)
      target.col <- which(bws$icon)[target.cont]
      return(as.vector(direct.out$grad[, target.col]))
    }

    if (exact.core.route) {
      H <- if (lc.derivative.exact.route) {
        .npRmpi_with_local_regression(.npreghat_exact_lc_derivative_matrix_from_npksum_chunked(
          bws = bws,
          txdat = txdat,
          exdat = if (no.ex) NULL else exdat,
          s = s
        ))
      } else if (exact.lc.kernel.route) {
        .npRmpi_with_local_regression(.npreghat_exact_lc_matrix_from_kernel_weights(
          bws = bws,
          txdat = txdat,
          exdat = if (no.ex) NULL else exdat
        ))
      } else if (exact.ll.kernel.route) {
        .npRmpi_with_local_regression(.npreghat_exact_ll_matrix_from_kernel_weights(
          bws = bws,
          txdat = txdat,
          exdat = if (no.ex) NULL else exdat,
          s = s
        ))
      } else if (exact.lp.kernel.route) {
        .npRmpi_with_local_regression(.npreghat_exact_lp_matrix_from_kernel_weights(
          bws = bws,
          txdat = txdat,
          exdat = if (no.ex) NULL else exdat,
          s = s,
          basis = reg.spec$basis.engine,
          degree = reg.spec$degree.engine,
          bernstein.basis = reg.spec$bernstein.basis.engine
        ))
      } else {
        .npreghat_exact_matrix_from_core(
          bws = bws,
          txdat = txdat,
          exdat = if (no.ex) NULL else exdat,
          s = s
        )
      }

      if (identical(output, "apply")) {
        if (is.null(y))
          stop("argument 'y' is required when output='apply'")
        out <- H %*% y
        if (ncol(out) == 1L)
          return(as.vector(out))
        return(out)
      }

      class(H) <- c("npreghat", "matrix")
      attr(H, "bws") <- bws
      attr(H, "txdat") <- txdat
      attr(H, "exdat") <- if (no.ex) txdat else exdat
      attr(H, "trainiseval") <- no.ex
      attr(H, "regtype") <- regtype
      attr(H, "degree") <- degree
      attr(H, "basis") <- basis
      attr(H, "bernstein.basis") <- bernstein.basis
      attr(H, "s") <- s
      attr(H, "leave.one.out") <- leave.one.out
      attr(H, "ridge.used") <- rep.int(0.0, nrow(H))
      attr(H, "rows.omit") <- rows.omit
      attr(H, "call") <- match.call(expand.dots = FALSE)

      if (!is.null(y)) {
        Hy <- H %*% y
        if (ncol(Hy) == 1L)
          Hy <- as.vector(Hy)
        attr(H, "Hy") <- Hy
      }

      return(H)
    }

    if (any(s > degree))
      stop("requested derivative order in 's' exceeds local polynomial degree")

    kw.args <- list(
      txdat = txdat,
      bws = bws,
      return.kernel.weights = TRUE,
      bandwidth.divide = TRUE,
      leave.one.out = leave.one.out
    )
    if (!no.ex)
      kw.args$exdat <- exdat

    use.local.kw <- identical(bws$type, "generalized_nn") &&
      identical(reg.spec$regtype.engine, "lp") &&
      .npRmpi_has_active_slave_pool(comm = 1L)

    kw <- if (use.local.kw) {
      direct.kw.args <- kw.args[intersect(
        names(kw.args),
        c("bws", "txdat", "exdat", "leave.one.out", "bandwidth.divide")
      )]
      do.call(.np_kernel_weights_direct, direct.kw.args)
    } else {
      kw.obj <- do.call(npksum.default, kw.args)
      kw.obj$kw
    }

    if (is.null(kw))
      stop("kernel-weight extraction failed")

    if (!is.matrix(kw))
      kw <- matrix(kw, nrow = nrow(txdat))
    if (leave.one.out && no.ex && nrow(kw) == ncol(kw))
      diag(kw) <- 0.0

    ntrain <- nrow(txdat)
    neval <- ncol(kw)

    W <- W.lp(
      xdat = txdat,
      degree = degree,
      basis = basis,
      bernstein.basis = bernstein.basis
    )

    eval.data <- if (no.ex) txdat else exdat
    W.eval <- W.lp(
      xdat = txdat,
      exdat = if (no.ex) NULL else eval.data,
      degree = degree,
      gradient.vec = if (any(s > 0L)) s else NULL,
      basis = basis,
      bernstein.basis = bernstein.basis
    )

    if (!is.matrix(W.eval))
      W.eval <- matrix(W.eval, nrow = neval)

    if (nrow(W.eval) != neval)
      W.eval <- matrix(W.eval, nrow = neval, byrow = FALSE)

    ridge.used <- rep.int(0.0, neval)

    if (identical(output, "matrix")) {
      H <- matrix(NA_real_, nrow = neval, ncol = ntrain)
    } else {
      if (is.null(y))
        stop("argument 'y' is required when output='apply'")
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
        stop(sprintf("failed to solve local hat system at evaluation row %d", i))

      ridge.used[i] <- solve.out$ridge
      h.row <- kw[, i] * drop(W %*% solve.out$v)

      if (identical(output, "matrix")) {
        H[i, ] <- h.row
      } else {
        out[i, ] <- drop(crossprod(h.row, y))
      }
    }

    if (identical(output, "apply")) {
      if (ncol(out) == 1L)
        return(as.vector(out))
      return(out)
    }

    class(H) <- c("npreghat", "matrix")
    attr(H, "bws") <- bws
    attr(H, "txdat") <- txdat
    attr(H, "exdat") <- if (no.ex) txdat else exdat
    attr(H, "trainiseval") <- no.ex
    attr(H, "regtype") <- regtype
    attr(H, "degree") <- degree
    attr(H, "basis") <- basis
    attr(H, "bernstein.basis") <- bernstein.basis
    attr(H, "s") <- s
    attr(H, "leave.one.out") <- leave.one.out
    attr(H, "ridge.used") <- ridge.used
    attr(H, "rows.omit") <- rows.omit
    attr(H, "call") <- match.call(expand.dots = FALSE)

    if (!is.null(y)) {
      Hy <- H %*% y
      if (ncol(Hy) == 1L)
        Hy <- as.vector(Hy)
      attr(H, "Hy") <- Hy
    }

    H
  }

npreghat.default <-
  function(bws, ...){
    stop("unsupported 'bws' type; provide an 'rbandwidth' or 'npregression' object")
  }

npreghat.NULL <- npreghat.default

predict.npreghat <-
  function(object,
           newdata = NULL,
           y = NULL,
           output = c("matrix", "apply"),
           s = attr(object, "s"),
           leave.one.out = attr(object, "leave.one.out"),
           deriv = NULL,
           ...){
    output <- match.arg(output)
    bws <- attr(object, "bws")
    txdat <- attr(object, "txdat")
    dots <- list(...)

    if (is.null(bws) || is.null(txdat))
      stop("object does not carry 'bws' and training data attributes")

    leave.one.out <- if (is.null(leave.one.out)) FALSE else leave.one.out
    call.args <- list(
      bws = bws,
      txdat = txdat,
      y = y,
      output = output,
      s = s,
      deriv = deriv,
      degree = attr(object, "degree"),
      basis = attr(object, "basis"),
      bernstein.basis = attr(object, "bernstein.basis"),
      leave.one.out = leave.one.out
    )
    if (!is.null(newdata) || !isTRUE(leave.one.out))
      call.args$exdat <- if (is.null(newdata)) attr(object, "exdat") else newdata
    do.call(npreghat, c(call.args, dots))
  }

print.npreghat <- function(x, ...) {
  s <- attr(x, "s")
  cat("\nNonparametric hat operator:", nrow(x), "evaluation x", ncol(x), "training\n")
  cat("regtype:", attr(x, "regtype"),
      "| basis:", attr(x, "basis"),
      "| s:", if (length(s)) paste(s, collapse = ",") else "none", "\n")
  invisible(x)
}
