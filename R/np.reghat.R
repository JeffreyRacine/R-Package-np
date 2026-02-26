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
  ridge <- max(0, as.double(ridge.base))
  solved <- FALSE

  for (attempt in 0:8) {
    A <- XtWX
    if (ridge > 0)
      diag(A) <- diag(A) + ridge

    v <- tryCatch(
      qr.solve(t(A), w.eval, tol = .Machine$double.eps),
      error = function(e) NULL
    )

    if (!is.null(v) && all(is.finite(v))) {
      solved <- TRUE
      break
    }

    ridge <- if (ridge > 0) ridge * 10 else 1e-12
  }

  if (!solved)
    return(NULL)

  list(v = v, ridge = ridge)
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
           s = NULL,
           deriv = NULL,
           degree = NULL,
           basis = NULL,
           bernstein.basis = NULL,
           ridge = 1.0e-12,
           ...){

    no.ex <- missing(exdat)
    if (.npRmpi_autodispatch_active() &&
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
          npRmpi:::npreghat.rbandwidth(
            bws = BWS,
            txdat = TXDAT,
            y = YDAT,
            output = OUTPUT,
            s = SVAL,
            deriv = DERIV,
            degree = DEGREE,
            basis = BASIS,
            bernstein.basis = BERN,
            ridge = RIDGE
          )
        } else {
          npRmpi:::npreghat.rbandwidth(
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
            ridge = RIDGE
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
        RIDGE = ridge
      ))
      return(.npRmpi_bcast_cmd_expr(expr, comm = 1L, caller.execute = TRUE))
    }

    output <- match.arg(output)
    dots <- list(...)
    npRejectLegacyLpArgs(names(dots), where = "npreghat")

    txdat <- toFrame(txdat)

    if (!no.ex) {
      exdat <- toFrame(exdat)
      if (!(txdat %~% exdat))
        stop("'txdat' and 'exdat' are not similar data frames!")
    }

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

    kw.args <- list(
      txdat = txdat,
      bws = bws,
      return.kernel.weights = TRUE,
      bandwidth.divide = TRUE,
      leave.one.out = FALSE
    )
    if (!no.ex)
      kw.args$exdat <- exdat

    kw.obj <- do.call(npksum.default, kw.args)
    kw <- kw.obj$kw

    if (is.null(kw))
      stop("kernel-weight extraction failed")

    if (!is.matrix(kw))
      kw <- matrix(kw, nrow = nrow(txdat))

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
           deriv = NULL,
           ...){
    output <- match.arg(output)
    bws <- attr(object, "bws")
    txdat <- attr(object, "txdat")

    if (is.null(bws) || is.null(txdat))
      stop("object does not carry 'bws' and training data attributes")

    exdat <- if (is.null(newdata)) attr(object, "exdat") else newdata

    npreghat(
      bws = bws,
      txdat = txdat,
      exdat = exdat,
      y = y,
      output = output,
      s = s,
      deriv = deriv,
      degree = attr(object, "degree"),
      basis = attr(object, "basis"),
      bernstein.basis = attr(object, "bernstein.basis"),
      ...
    )
  }

print.npreghat <- function(x, ...) {
  s <- attr(x, "s")
  cat("\nNonparametric hat operator:", nrow(x), "evaluation x", ncol(x), "training\n")
  cat("regtype:", attr(x, "regtype"),
      "| basis:", attr(x, "basis"),
      "| s:", if (length(s)) paste(s, collapse = ",") else "none", "\n")
  invisible(x)
}
