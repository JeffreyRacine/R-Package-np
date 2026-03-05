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
  solved <- FALSE
  ridge <- ridge.grid[1L]

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

.np_kernel_weights_direct <- function(bws,
                                      txdat,
                                      exdat = NULL,
                                      leave.one.out = FALSE,
                                      bandwidth.divide = TRUE,
                                      kernel.pow = 1.0) {
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
    compute.ocg = FALSE
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
    as.integer(rep.int(ALL_OPERATORS[["normal"]], length(txdat))),
    as.integer(myopti), as.double(kernel.pow),
    as.integer(enrow),
    as.integer(0L),
    as.integer(nkw),
    as.double(cker.bounds.c$lb),
    as.double(cker.bounds.c$ub),
    PACKAGE = "np"
  )

  kw <- matrix(as.double(myout[["kernel.weights"]]), nrow = tnrow, ncol = enrow)
  if (leave.one.out && miss.ex && nrow(kw) == ncol(kw))
    diag(kw) <- 0.0
  kw
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

    output <- match.arg(output)
    dots <- list(...)
    npRejectLegacyLpArgs(names(dots), where = "npreghat")

    txdat <- toFrame(txdat)
    no.ex <- missing(exdat)
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

    if (any(s > degree))
      stop("requested derivative order in 's' exceeds local polynomial degree")

    if (all(degree == 0L) && any(s > 0L))
      stop("local-constant derivative hat operators are not available yet; use local polynomial degree >= 1")

    kw <- .np_kernel_weights_direct(
      bws = bws,
      txdat = txdat,
      exdat = if (no.ex) NULL else exdat,
      leave.one.out = leave.one.out,
      bandwidth.divide = TRUE,
      kernel.pow = 1.0
    )

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
