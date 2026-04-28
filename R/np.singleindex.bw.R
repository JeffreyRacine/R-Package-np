
npindexbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npindexbw")
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat"),
                                     eval_env = parent.frame())
    UseMethod("npindexbw", target)
  }

npindexbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    formula.terms <- terms(formula)
    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = formula.terms,
                        data = environment(formula),
                        eval_env = environment(formula))
    else .np_terms_ts_mask(terms_obj = formula.terms,
                           data = data,
                           eval_env = environment(formula))

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    mf[[1]] <- as.name("model.frame")

    if(all(orig.ts)){
      args <- (as.list(attr(formula.terms, "variables"))[-1])
      formula <- formula.terms
      attr(formula, "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
      mf[["formula"]] <- formula
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(formula.terms, "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      formula <- formula.terms
      attr(formula, "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
      mf[["formula"]] <- formula
    }

    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())

    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]

    tbw <- npindexbw(xdat = xdat, ydat = ydat, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)

    tbw
  }

npindexbw.NULL <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npindexbw()")
    mc <- match.call(expand.dots = FALSE)
    dots <- list(...)
    dot.names <- names(dots)
    nomad.requested <- "nomad" %in% dot.names && isTRUE(dots$nomad)
    degree.select.value <- if ("degree.select" %in% dot.names) {
      match.arg(as.character(dots$degree.select[[1L]]),
                c("manual", "coordinate", "exhaustive"))
    } else {
      "manual"
    }
    automatic.degree.search <- isTRUE(nomad.requested) ||
      !identical(degree.select.value, "manual")
    regtype.value <- if ("regtype" %in% dot.names) {
      match.arg(as.character(dots$regtype[[1L]]), c("lc", "ll", "lp"))
    } else if (isTRUE(nomad.requested)) {
      "lp"
    } else {
      "lc"
    }
    search.engine.value <- if ("search.engine" %in% dot.names) {
      match.arg(as.character(dots$search.engine[[1L]]), c("nomad+powell", "cell", "nomad"))
    } else if (isTRUE(nomad.requested)) {
      "nomad+powell"
    } else {
      "nomad+powell"
    }
    .np_nomad_validate_inner_multistart(
      call_names = names(mc),
      dot.args = dots,
      regtype = regtype.value,
      automatic.degree.search = automatic.degree.search,
      search.engine = search.engine.value
    )
    if (.npRmpi_autodispatch_active() && !isTRUE(automatic.degree.search))
      return(.npRmpi_autodispatch_call(mc, parent.frame()))

    xdat <- toFrame(xdat)

    bws <- double(ncol(xdat)+1)

    tbw <- npindexbw.default(xdat = xdat,
                             ydat = ydat,
                             bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <- updateBwNameMetadata(nameList =
                                list(ynames = deparse(substitute(ydat))),
                                bws = tbw)

    tbw
  }

.npindex_resolve_spec <- function(source, where = "npindex") {
  if (!is.null(source$regtype.engine)) {
    return(list(
      regtype = if (is.null(source$regtype)) "lc" else as.character(source$regtype),
      basis = if (is.null(source$basis)) "glp" else as.character(source$basis),
      degree = if (is.null(source$degree)) integer(0) else as.integer(source$degree),
      bernstein.basis = isTRUE(source$bernstein.basis),
      regtype.engine = as.character(source$regtype.engine),
      basis.engine = if (is.null(source$basis.engine)) "glp" else as.character(source$basis.engine),
      degree.engine = if (is.null(source$degree.engine)) integer(0) else as.integer(source$degree.engine),
      bernstein.basis.engine = isTRUE(source$bernstein.basis.engine)
    ))
  }

  npCanonicalConditionalRegSpec(
    regtype = if (is.null(source$regtype)) "lc" else as.character(source$regtype),
    basis = if (is.null(source$basis)) "glp" else as.character(source$basis),
    degree = source$degree,
    bernstein.basis = isTRUE(source$bernstein.basis),
    ncon = 1L,
    where = where
  )
}

.npindex_nn_candidate_bandwidth <- function(h, bwtype, nobs) {
  if (identical(bwtype, "fixed")) {
    return(list(ok = is.finite(h) && (h > 0), value = as.double(h)))
  }

  if (!is.finite(h)) {
    return(list(ok = FALSE, value = NA_real_))
  }

  lower <- 2L
  upper <- max(1L, as.integer(nobs) - 1L)
  k <- .np_round_half_to_even(h)
  list(ok = (k >= lower) && (k <= upper), value = as.double(k))
}

.npindexbw_h_start_controls <- function(scale.factor.init.lower = 0.1,
                                        scale.factor.init.upper = 2.0,
                                        scale.factor.init = 0.5,
                                        scale.factor.search.lower = 0,
                                        where = "npindexbw") {
  cont.start <- npContinuousSearchStartControls(
    scale.factor.init.lower,
    scale.factor.init.upper,
    scale.factor.init,
    scale.factor.search.lower,
    where = where
  )
  cont.start$scale.factor.search.lower <- as.double(scale.factor.search.lower)
  cont.start
}

.npindex_start_bandwidth_scale <- function(fit, nobs) {
  EssDee(fit) * nobs^(-1 / 5)
}

.npindex_default_start_bandwidth <- function(fit,
                                             bwtype,
                                             nobs,
                                             start.controls = .npindexbw_h_start_controls()) {
  if (identical(bwtype, "fixed")) {
    return(start.controls$scale.factor.init * .npindex_start_bandwidth_scale(fit = fit, nobs = nobs))
  }

  lower <- 2L
  max(lower, min(max(1L, as.integer(nobs) - 1L), .np_round_half_to_even(sqrt(nobs))))
}

.npindex_random_start_bandwidth <- function(fit,
                                            bwtype,
                                            nobs,
                                            start.controls = .npindexbw_h_start_controls()) {
  if (identical(bwtype, "fixed")) {
    return(runif(1, min = start.controls$scale.factor.init.lower, max = start.controls$scale.factor.init.upper) *
             .npindex_start_bandwidth_scale(fit = fit, nobs = nobs))
  }

  upper <- max(1L, as.integer(nobs) - 1L)
  runif(1, min = 2, max = max(2L, upper))
}

.npindex_finalize_bandwidth <- function(h,
                                        bwtype,
                                        nobs,
                                        lower = NULL,
                                        where = "npindexbw") {
  candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bwtype, nobs = nobs)
  if (!candidate$ok) {
    if (identical(bwtype, "fixed")) {
      stop(sprintf("%s: bandwidth must be positive and finite", where), call. = FALSE)
    }
    stop(
      sprintf(
        "%s: nearest-neighbor bandwidth candidate must map to an integer in [2, %d]",
        where,
        max(2L, as.integer(nobs) - 1L)
      ),
      call. = FALSE
    )
  }

  if (identical(bwtype, "fixed") && !is.null(lower) && candidate$value < lower) {
    stop(sprintf("%s: bandwidth is below the continuous scale-factor lower bound", where),
         call. = FALSE)
  }

  candidate$value
}

.npindexbw_nomad_fixed_h_scale <- function(fit,
                                           h.start.raw,
                                           nobs,
                                           start.controls = .npindexbw_h_start_controls()) {
  scale <- .npindex_default_start_bandwidth(
    fit = fit,
    bwtype = "fixed",
    nobs = nobs,
    start.controls = start.controls
  )
  if (!is.finite(scale) || scale <= 0) {
    scale <- max(
      if (is.finite(h.start.raw)) abs(as.double(h.start.raw)) else 0,
      1e-3
    )
  }
  as.double(scale)
}

.npindexbw_nomad_fixed_start_setup <- function(xmat,
                                               ydat,
                                               baseline.bws,
                                               degree.search,
                                               nmulti,
                                               random.seed,
                                               h.start.controls = .npindexbw_h_start_controls()) {
  p <- ncol(xmat)
  nobs <- nrow(xmat)
  beta.free <- if (p > 1L) seq_len(p - 1L) else integer(0)
  h.col <- length(beta.free) + 1L
  degree.col <- h.col + 1L

  ols.fit <- lm(ydat ~ xmat, x = TRUE)
  fit.proxy <- as.double(fitted(ols.fit))
  fit.proxy[!is.finite(fit.proxy)] <- mean(ydat[is.finite(ydat)])
  fit.proxy[!is.finite(fit.proxy)] <- 0

  ols.beta <- if (length(beta.free)) {
    as.double(coef(ols.fit)[3:ncol(ols.fit$x)])
  } else {
    numeric(0)
  }
  if (length(ols.beta))
    ols.beta[!is.finite(ols.beta)] <- 0

  beta.start.raw <- if (length(beta.free)) {
    beta.user <- as.double(baseline.bws$beta[beta.free + 1L])
    if (setequal(beta.user, c(0))) ols.beta else beta.user
  } else {
    numeric(0)
  }
  if (length(beta.start.raw))
    beta.start.raw[!is.finite(beta.start.raw)] <- 0

  if (isTRUE(all.equal(as.double(baseline.bws$bw[1L]), 0))) {
    h.start.raw <- .npindex_default_start_bandwidth(
      fit = fit.proxy,
      bwtype = "fixed",
      nobs = nobs,
      start.controls = h.start.controls
    )
  } else {
    h.lower.raw <- h.start.controls$scale.factor.search.lower *
      .npindex_start_bandwidth_scale(fit = fit.proxy, nobs = nobs)
    h.start.raw <- tryCatch(
      .npindex_finalize_bandwidth(
        h = baseline.bws$bw[1L],
        bwtype = "fixed",
        nobs = nobs,
        lower = h.lower.raw,
        where = "npindexbw"
      ),
      error = function(e) .npindex_default_start_bandwidth(
        fit = fit.proxy,
        bwtype = "fixed",
        nobs = nobs,
        start.controls = h.start.controls
      )
    )
  }
  if (!is.finite(h.start.raw) || h.start.raw <= 0)
    h.start.raw <- 1e-3

  h.scale <- .npindexbw_nomad_fixed_h_scale(
    fit = fit.proxy,
    h.start.raw = h.start.raw,
    nobs = nobs,
    start.controls = h.start.controls
  )
  degree.starts <- .np_lp_nomad_build_degree_starts(
    initial = degree.search$start.degree,
    lower = degree.search$lower,
    upper = degree.search$upper,
    basis = degree.search$basis,
    nobs = degree.search$nobs,
    nmulti = nmulti,
    random.seed = random.seed,
    user_supplied = isTRUE(degree.search$start.user)
  )

  start_matrix.raw <- matrix(0, nrow = nmulti, ncol = degree.col)
  if (length(beta.free))
    start_matrix.raw[1L, beta.free] <- beta.start.raw
  start_matrix.raw[1L, h.col] <- h.start.raw
  start_matrix.raw[, degree.col] <- as.integer(degree.starts[, 1L])

  if (nmulti > 1L) {
    seed.state <- .np_seed_enter(random.seed)
    on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)

    for (j in 2:nmulti) {
      if (length(beta.free)) {
        beta.rand <- runif(length(ols.beta), min = 0.5, max = 1.5) * ols.beta
        beta.rand[!is.finite(beta.rand)] <- 0
        start_matrix.raw[j, beta.free] <- beta.rand
      }

      h.rand <- .npindex_random_start_bandwidth(
        fit = fit.proxy,
        bwtype = "fixed",
        nobs = nobs,
        start.controls = h.start.controls
      )
      if (!is.finite(h.rand) || h.rand <= 0)
        h.rand <- h.start.raw
      start_matrix.raw[j, h.col] <- h.rand
    }
  }

  start_matrix.point <- start_matrix.raw
  start_matrix.point[, h.col] <- start_matrix.point[, h.col] / h.scale

  list(
    beta.free = beta.free,
    ols.beta = ols.beta,
    h.col = h.col,
    degree.col = degree.col,
    h.scale = h.scale,
    start_matrix.raw = start_matrix.raw,
    start_matrix.point = start_matrix.point
  )
}

.npindexbw_fast_eligible <- function(h, bws, eval.index) {
  if (!identical(bws$type, "fixed"))
    return(FALSE)

  fast_largeh_tol <- getOption("np.largeh.rel.tol", 1e-3)
  if (!is.numeric(fast_largeh_tol) || length(fast_largeh_tol) != 1L ||
      is.na(fast_largeh_tol) || !is.finite(fast_largeh_tol) ||
      fast_largeh_tol <= 0 || fast_largeh_tol >= 0.1)
    fast_largeh_tol <- 1e-3

  cont_utol <- switch(
    bws$ckertype,
    gaussian = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
    "truncated gaussian" = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
    epanechnikov = sqrt(fast_largeh_tol),
    uniform = 1.0 - 32.0 * .Machine$double.eps,
    0.0
  )

  if (!is.finite(cont_utol) || cont_utol <= 0 || !is.finite(h) || h <= 0)
    return(FALSE)

  vals <- as.double(eval.index)
  vals <- vals[is.finite(vals)]
  if (!length(vals))
    return(FALSE)

  h >= (diff(range(vals)) / cont_utol)
}

.npindex_selector_with_local_kernelsum <- function(expr) {
  if (!isTRUE(.npRmpi_autodispatch_active()))
    return(force(expr))
  if (!isTRUE(.npRmpi_autodispatch_called_from_bcast()))
    return(force(expr))
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.local <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old.local), add = TRUE)
  old.mode <- .Call("C_np_set_local_regression_mode", TRUE, PACKAGE = "npRmpi")
  on.exit(.Call("C_np_set_local_regression_mode", old.mode, PACKAGE = "npRmpi"), add = TRUE)
  force(expr)
}

.npindexbw_build_lp_regression_leaf <- function(index,
                                                ydat,
                                                h,
                                                bws,
                                                spec) {
  index.df <- data.frame(index = as.double(index))
  reg.args <- list(
    regtype = "lp",
    basis = as.character(spec$basis.engine),
    degree = as.integer(spec$degree.engine),
    bernstein.basis = isTRUE(spec$bernstein.basis.engine),
    bwmethod = "cv.ls",
    bwtype = bws$type,
    ckertype = bws$ckertype,
    ckerorder = bws$ckerorder,
    ckerbound = bws$ckerbound,
    ckerlb = bws$ckerlb,
    ckerub = bws$ckerub,
    ukertype = if (is.null(bws$ukertype)) "aitchisonaitken" else bws$ukertype,
    okertype = if (is.null(bws$okertype)) "liracine" else bws$okertype
  )

  list(
    xdat = index.df,
    bws = .npregbw_build_rbandwidth(
      xdat = index.df,
      ydat = ydat,
      bws = c(h),
      bandwidth.compute = FALSE,
      reg.args = reg.args,
      yname = if (is.null(bws$ynames)) "y" else as.character(bws$ynames[1L])
    )
  )
}

.npindexbw_eval_ichimura_lp_via_npreg <- function(index,
                                                  ydat,
                                                  h,
                                                  bws,
                                                  spec,
                                                  invalid.penalty) {
  leaf <- .npindexbw_build_lp_regression_leaf(
    index = index,
    ydat = ydat,
    h = h,
    bws = bws,
    spec = spec
  )

  out <- tryCatch(
    .npregbw_eval_only(
      xdat = leaf$xdat,
      ydat = ydat,
      bws = leaf$bws,
      invalid.penalty = "baseline",
      penalty.multiplier = 10
    ),
    error = function(e) NULL
  )

  if (is.null(out) || !is.finite(out$objective[1L]))
    return(list(objective = as.numeric(invalid.penalty), num.feval.fast = 0L))

  list(
    objective = as.numeric(out$objective[1L]),
    num.feval.fast = as.numeric(out$num.feval.fast[1L])
  )
}

.npindex_lp_loo_fit_block <- function(idx,
                                      index,
                                      ydat,
                                      h,
                                      bws,
                                      spec) {
  idx <- as.integer(idx)
  if (!length(idx))
    return(list(index = idx, fitted = numeric(0L)))

  index <- as.double(index)
  ydat <- as.double(ydat)
  n <- length(index)
  index.df <- data.frame(index = index)
  kbw <- kbandwidth(
    bw = c(h),
    bwtype = bws$type,
    ckertype = bws$ckertype,
    ckerorder = bws$ckerorder,
    ckerbound = bws$ckerbound,
    ckerlb = bws$ckerlb,
    ckerub = bws$ckerub,
    nobs = n,
    xdati = untangle(index.df),
    ydati = untangle(data.frame(ydat)),
    xnames = "index",
    ynames = bws$ynames
  )
  degree <- as.integer(spec$degree.engine)
  W <- W.lp(
    xdat = index.df,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  W.eval.block <- W.lp(
    xdat = index.df,
    exdat = index.df[idx, , drop = FALSE],
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  if (!is.matrix(W))
    W <- matrix(W, nrow = n)
  if (!is.matrix(W.eval.block))
    W.eval.block <- matrix(W.eval.block, nrow = length(idx))

  kw <- .np_kernel_weights_direct(
    bws = kbw,
    txdat = index.df,
    exdat = index.df[idx, , drop = FALSE],
    bandwidth.divide = TRUE,
    kernel.pow = 1.0
  )
  kw[cbind(idx, seq_along(idx))] <- 0.0

  ridge.grid <- npRidgeSequenceFromBase(n.train = n, ridge.base = 0.0, cap = 1.0)
  fit.block <- double(length(idx))

  for (jj in seq_along(idx)) {
    w <- kw[, jj]
    XtWy0 <- sum(w * ydat)
    lc.mean <- XtWy0 / NZD(sum(w))
    A <- crossprod(W, W * w)
    z <- crossprod(W, ydat * w)
    solved <- FALSE

    for (ridge in ridge.grid) {
      Ar <- A
      zr <- z
      if (ridge > 0) {
        diag(Ar) <- diag(Ar) + ridge
        zr[1L] <- XtWy0 + ridge * XtWy0 / NZD(A[1L, 1L])
      }

      coef <- tryCatch(chol2inv(chol(Ar)) %*% zr, error = function(e) NULL)
      if (!is.null(coef) && all(is.finite(coef))) {
        alpha <- min(1.0, ridge)
        fit.block[jj] <- (1.0 - alpha) * drop(W.eval.block[jj, , drop = FALSE] %*% coef) +
          alpha * lc.mean
        solved <- TRUE
        break
      }
    }

    if (!solved)
      fit.block[jj] <- lc.mean
  }

  list(index = idx, fitted = fit.block)
}

.npindex_lp_loo_fit <- function(index,
                                ydat,
                                h,
                                bws,
                                spec,
                                chunk.size = 128L) {
  index <- as.double(index)
  ydat <- as.double(ydat)
  n <- length(index)
  if (length(ydat) != n)
    stop("index and ydat must have the same length")

  if (.npRmpi_autodispatch_active() &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    slave.num <- mpi.comm.size(1L) - 1L
    if (slave.num > 0L) {
      worker.spec <- list(
        degree.engine = spec$degree.engine,
        basis.engine = spec$basis.engine,
        bernstein.basis.engine = spec$bernstein.basis.engine
      )
      blocks <- Filter(length, .splitIndices(n, min(as.integer(n), as.integer(slave.num))))
      pieces <- mpi.apply(
        blocks,
        .npindex_lp_loo_fit_block,
        index = index,
        ydat = ydat,
        h = h,
        bws = bws,
        spec = worker.spec,
        comm = 1L
      )
      fit.loo <- double(n)
      for (piece in pieces) {
        if (inherits(piece, "try-error"))
          stop(conditionMessage(attr(piece, "condition", exact = TRUE)), call. = FALSE)
        fit.loo[as.integer(piece$index)] <- as.double(piece$fitted)
      }
      return(fit.loo)
    }
  }

  index.df <- data.frame(index = index)
  kbw <- kbandwidth(
    bw = c(h),
    bwtype = bws$type,
    ckertype = bws$ckertype,
    ckerorder = bws$ckerorder,
    ckerbound = bws$ckerbound,
    ckerlb = bws$ckerlb,
    ckerub = bws$ckerub,
    nobs = n,
    xdati = untangle(index.df),
    ydati = untangle(data.frame(ydat)),
    xnames = "index",
    ynames = bws$ynames
  )
  degree <- as.integer(spec$degree.engine)
  W <- W.lp(
    xdat = index.df,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  W.eval <- W.lp(
    xdat = index.df,
    exdat = index.df,
    degree = degree,
    basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.basis.engine
  )
  if (!is.matrix(W))
    W <- matrix(W, nrow = n)
  if (!is.matrix(W.eval))
    W.eval <- matrix(W.eval, nrow = n)

  ridge.grid <- npRidgeSequenceFromBase(n.train = n, ridge.base = 0.0, cap = 1.0)
  fit.loo <- double(n)
  chunk.size <- max(1L, npValidatePositiveInteger(chunk.size, "chunk.size"))

  for (start in seq.int(1L, n, by = chunk.size)) {
    idx <- seq.int(start, min(n, start + chunk.size - 1L))
    kw <- .np_kernel_weights_direct(
      bws = kbw,
      txdat = index.df,
      exdat = index.df[idx, , drop = FALSE],
      bandwidth.divide = TRUE,
      kernel.pow = 1.0
    )

    kw[cbind(idx, seq_along(idx))] <- 0.0
    W.eval.block <- W.eval[idx, , drop = FALSE]

    for (jj in seq_along(idx)) {
      w <- kw[, jj]
      XtWy0 <- sum(w * ydat)
      lc.mean <- XtWy0 / NZD(sum(w))
      A <- crossprod(W, W * w)
      z <- crossprod(W, ydat * w)
      solved <- FALSE

      for (ridge in ridge.grid) {
        Ar <- A
        zr <- z
        if (ridge > 0) {
          diag(Ar) <- diag(Ar) + ridge
          zr[1L] <- XtWy0 + ridge * XtWy0 / NZD(A[1L, 1L])
        }

        coef <- tryCatch(chol2inv(chol(Ar)) %*% zr, error = function(e) NULL)
        if (!is.null(coef) && all(is.finite(coef))) {
          alpha <- min(1.0, ridge)
          fit.loo[idx[jj]] <- (1.0 - alpha) * drop(W.eval.block[jj, , drop = FALSE] %*% coef) +
            alpha * lc.mean
          solved <- TRUE
          break
        }
      }

      if (!solved)
        fit.loo[idx[jj]] <- lc.mean
    }
  }

  fit.loo
}

.npindexbw_build_sibandwidth <- function(xdat,
                                         ydat,
                                         bws,
                                         template,
                                         bandwidth.compute,
                                         reg.args) {
  p <- ncol(xdat)
  out <- do.call(
    sibandwidth,
    c(
      list(
        beta = bws[seq_len(p)],
        h = bws[p + 1L],
        method = template$method,
        ckertype = template$ckertype,
        ckerorder = template$ckerorder,
        ckerbound = template$ckerbound,
        ckerlb = template$ckerlb,
        ckerub = template$ckerub,
        bwtype = template$type,
        nobs = nrow(xdat),
        xdati = untangle(xdat),
        ydati = untangle(data.frame(ydat)),
        xnames = names(xdat),
        ynames = template$ynames,
        bandwidth = bws[p + 1L],
        bandwidth.compute = bandwidth.compute
      ),
      reg.args
    )
  )
  if (!is.null(reg.args$scale.factor.search.lower))
    out$scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      reg.args$scale.factor.search.lower
    )
  out
}

.npindexbw_run_fixed_degree <- function(xdat, ydat, bws, template, reg.args, opt.args) {
  tbw <- .npindexbw_build_sibandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    template = template,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args
  )

  do.call(npindexbw.sibandwidth, c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args))
}

.npindexbw_eval_objective <- function(param,
                                      xmat,
                                      ydat,
                                      bws,
                                      spec) {
  p <- ncol(xmat)
  beta.idx <- if (p > 1L) seq_len(p - 1L) else integer(0)
  beta <- if (length(beta.idx)) as.double(param[beta.idx]) else numeric(0)
  h <- as.double(param[p])
  nobs <- nrow(xmat)

  if (identical(bws$method, "ichimura")) {
    invalid.penalty <- 10 * mean(ydat^2)
  } else {
    invalid.penalty <- sqrt(.Machine$double.xmax)
  }

  h.candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bws$type, nobs = nobs)
  if (!h.candidate$ok)
    return(list(objective = invalid.penalty, num.feval.fast = 0L))
  h <- h.candidate$value

  index <- xmat %*% c(1, beta)
  wmat <- cbind(ydat, 1.0)

  if (!identical(spec$regtype.engine, "lc") &&
      identical(bws$method, "ichimura")) {
    return(.npindexbw_eval_ichimura_lp_via_npreg(
      index = index,
      ydat = ydat,
      h = h,
      bws = bws,
      spec = spec,
      invalid.penalty = invalid.penalty
    ))
  }

  fit.loo <- if (identical(spec$regtype.engine, "lc")) {
    tww <- tryCatch(
      .npRmpi_with_local_regression(
        npksum(
          txdat = index,
          tydat = wmat,
          weights = wmat,
          leave.one.out = TRUE,
          bandwidth.divide = TRUE,
          bws = c(h),
          bwtype = bws$type,
          ckertype = bws$ckertype,
          ckerorder = bws$ckerorder,
          ckerbound = bws$ckerbound,
          ckerlb = bws$ckerlb,
          ckerub = bws$ckerub
        )
      )$ksum,
      error = function(e) NULL
    )
    if (is.null(tww))
      return(list(objective = invalid.penalty, num.feval.fast = 0L))
    tww[1, 2, ] / NZD(tww[2, 2, ])
  } else {
    ok.design <- tryCatch({
      npCheckRegressionDesignCondition(
        reg.code = REGTYPE_LP,
        xcon = data.frame(index = index),
        basis = spec$basis.engine,
        degree = spec$degree.engine,
        bernstein.basis = spec$bernstein.basis.engine,
        where = "npindexbw"
      )
      TRUE
    }, error = function(e) FALSE)

    if (!ok.design)
      return(list(objective = invalid.penalty, num.feval.fast = 0L))

    if (identical(bws$method, "kleinspady") &&
        .npRmpi_has_active_slave_pool(comm = 1L) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
      .npRmpi_with_local_regression(
        .npindex_lp_loo_fit(
          index = index,
          ydat = ydat,
          h = h,
          bws = bws,
          spec = spec
        )
      )
    } else {
      .npindex_lp_loo_fit(
        index = index,
        ydat = ydat,
        h = h,
        bws = bws,
        spec = spec
      )
    }
  }

  if (any(!is.finite(fit.loo)))
    return(list(objective = invalid.penalty, num.feval.fast = 0L))

  if (identical(bws$method, "ichimura")) {
    return(list(
      objective = mean((ydat - fit.loo)^2),
      num.feval.fast = if (.npindexbw_fast_eligible(h = h, bws = bws, eval.index = index)) 1L else 0L
    ))
  }

  floor <- sqrt(.Machine$double.eps)
  fit.loo[fit.loo < floor] <- floor
  fit.loo[fit.loo > 1 - floor] <- 1 - floor
  list(
    objective = -mean(ydat * log(fit.loo) + (1.0 - ydat) * log1p(-fit.loo)),
    num.feval.fast = if (.npindexbw_fast_eligible(h = h, bws = bws, eval.index = index)) 1L else 0L
  )
}

.npindexbw_nomad_search <- function(xdat,
                                    ydat,
                                    bws,
                                    template,
                                    reg.args,
                                    opt.args,
                                    degree.search,
                                    nomad.inner.nmulti = 0L) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis
  template.reg.args$regtype.engine <- "lp"
  template.reg.args$degree.engine <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

  baseline.bws <- .npindexbw_build_sibandwidth(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    template = template,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args
  )

  keep.rows <- rep_len(TRUE, nrow(xdat))
  rows.omit <- attr(na.omit(data.frame(xdat, ydat)), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE
  x.clean <- toMatrix(xdat[keep.rows, , drop = FALSE])
  y.clean <- ydat[keep.rows]
  if (is.factor(y.clean))
    y.clean <- dlev(y.clean)[as.integer(y.clean)]
  else
    y.clean <- as.double(y.clean)

  p <- ncol(x.clean)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(ncol(xdat)) else npValidateNmulti(opt.args$nmulti[1L])
  scale.factor.search.lower <- npResolveScaleFactorLowerBound(opt.args$scale.factor.search.lower)
  h.start.controls <- .npindexbw_h_start_controls(
    scale.factor.init.lower = if (is.null(opt.args$scale.factor.init.lower)) 0.1 else opt.args$scale.factor.init.lower,
    scale.factor.init.upper = if (is.null(opt.args$scale.factor.init.upper)) 2.0 else opt.args$scale.factor.init.upper,
    scale.factor.init = if (is.null(opt.args$scale.factor.init)) 0.5 else opt.args$scale.factor.init,
    scale.factor.search.lower = scale.factor.search.lower,
    where = "npindexbw"
  )
  fixed.nomad <- identical(baseline.bws$type, "fixed")
  fixed.setup <- if (fixed.nomad) {
    .npindexbw_nomad_fixed_start_setup(
      xmat = x.clean,
      ydat = y.clean,
      baseline.bws = baseline.bws,
      degree.search = degree.search,
      nmulti = nomad.nmulti,
      random.seed = if (!is.null(opt.args$random.seed)) opt.args$random.seed else 42L,
      h.start.controls = h.start.controls
    )
  } else {
    NULL
  }
  beta.free <- if (fixed.nomad) fixed.setup$beta.free else if (p > 1L) seq_len(p - 1L) else integer(0)
  beta.start <- if (fixed.nomad) {
    if (length(beta.free)) as.double(fixed.setup$start_matrix.raw[1L, beta.free]) else numeric(0)
  } else {
    if (length(beta.free)) as.double(baseline.bws$beta[beta.free + 1L]) else numeric(0)
  }
  h.start.raw <- if (fixed.nomad) {
    as.double(fixed.setup$start_matrix.raw[1L, fixed.setup$h.col])
  } else {
    as.double(baseline.bws$bw[1L])
  }
  beta.lower <- if (length(beta.start)) {
    if (fixed.nomad) {
      pmin(0.5 * fixed.setup$ols.beta, 1.5 * fixed.setup$ols.beta)
    } else {
      -pmax(10, 10 * abs(beta.start))
    }
  } else {
    numeric(0)
  }
  beta.upper <- if (length(beta.start)) {
    if (fixed.nomad) {
      pmax(0.5 * fixed.setup$ols.beta, 1.5 * fixed.setup$ols.beta)
    } else {
      pmax(10, 10 * abs(beta.start))
    }
  } else {
    numeric(0)
  }
  h.lower.raw <- if (fixed.nomad) 1e-3 else 2
  h.upper.raw <- if (fixed.nomad) {
    max(1e6, abs(h.start.raw) * 1e3, 1)
  } else {
    max(2L, as.integer(nrow(x.clean)) - 1L)
  }
  h.lower <- if (fixed.nomad) {
    h.start.controls$scale.factor.search.lower
  } else {
    h.lower.raw
  }
  h.upper <- if (fixed.nomad) {
    h.start.controls$scale.factor.init.upper
  } else {
    h.upper.raw
  }
  h.integer <- !identical(baseline.bws$type, "fixed")

  x0 <- if (fixed.nomad) {
    as.numeric(fixed.setup$start_matrix.point[1L, ])
  } else {
    c(beta.start, h.start.raw, as.integer(degree.search$start.degree))
  }
  lb <- c(beta.lower, h.lower, degree.search$lower)
  ub <- c(beta.upper, h.upper, degree.search$upper)
  bbin <- c(rep.int(0L, length(beta.start)), if (isTRUE(h.integer)) 1L else 0L, 1L)
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0

  .np_nomad_baseline_note(degree.search$start.degree)

  point_h_to_raw <- function(h.point) {
    h.raw <- if (fixed.nomad) as.double(h.point) * fixed.setup$h.scale else as.double(h.point)
    as.double(h.raw)
  }

  point_to_public <- function(point) {
    point <- as.numeric(point)
    if (length(point) >= (length(beta.free) + 1L))
      point[length(beta.free) + 1L] <- point_h_to_raw(point[length(beta.free) + 1L])
    point
  }

  bandwidth_point_to_public <- function(point) {
    point <- as.numeric(point)
    if (length(point))
      point[length(point)] <- point_h_to_raw(point[length(point)])
    point
  }

  point_to_param <- function(point) {
    beta.tail <- if (length(beta.free)) as.double(point[seq_along(beta.free)]) else numeric(0)
    h <- point_h_to_raw(point[length(beta.free) + 1L])
    h <- .npindex_finalize_bandwidth(
      h = h,
      bwtype = baseline.bws$type,
      nobs = nrow(x.clean),
      lower = if (fixed.nomad) h.start.controls$scale.factor.search.lower * fixed.setup$h.scale else NULL,
      where = "npindexbw"
    )
    c(beta.tail, h)
  }

  point_to_bws <- function(point) {
    param <- point_to_param(point)
    beta.tail <- if (length(beta.free)) param[seq_along(beta.free)] else numeric(0)
    h <- param[length(beta.free) + 1L]
    c(1.0, beta.tail, h)
  }

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- .np_degree_clip_to_grid(
      as.integer(round(point[length(point)])),
      degree.search$candidates
    )
    eval.spec <- reg.args
    eval.spec$regtype.engine <- "lp"
    eval.spec$degree.engine <- degree
    eval.spec$bernstein.basis.engine <- degree.search$bernstein.basis
    eval.spec$basis.engine <- reg.args$basis.engine
    objective <- .npindexbw_eval_objective(
      param = point_to_param(point),
      xmat = x.clean,
      ydat = y.clean,
      bws = baseline.bws,
      spec = eval.spec
    )
    nomad.num.feval.total <<- nomad.num.feval.total + 1L
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(objective$num.feval.fast[1L])
    list(
      objective = as.numeric(objective$objective[1L]),
      degree = as.integer(degree),
      num.feval = 1,
      num.feval.fast = as.numeric(objective$num.feval.fast[1L])
    )
  }

  build_payload <- function(point, best_record, solution, interrupted) {
    point <- as.numeric(point)
    degree <- as.integer(best_record$degree)
    bw.vec <- point_to_bws(point)
    powell.elapsed <- NA_real_

    build_direct_payload <- function() {
      final.reg.args <- reg.args
      final.reg.args$regtype <- "lp"
      final.reg.args$degree <- degree
      final.reg.args$bernstein.basis <- degree.search$bernstein.basis
      final.reg.args$regtype.engine <- "lp"
      final.reg.args$degree.engine <- degree
      final.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis

      tbw <- .npindexbw_build_sibandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = bw.vec,
        template = baseline.bws,
        bandwidth.compute = FALSE,
        reg.args = final.reg.args
      )
      tbw$fval <- as.numeric(best_record$objective)
      tbw$num.feval <- as.numeric(nomad.num.feval.total)
      tbw$num.feval.fast <- as.numeric(nomad.num.feval.fast.total)
      tbw$total.time <- NA_real_
      payload <- npindexbw.sibandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = tbw,
        bandwidth.compute = FALSE
      )
      if (!is.null(payload$method) && length(payload$method))
        payload$pmethod <- bwmToPrint(as.character(payload$method[1L]))
      payload
    }

    direct.payload <- build_direct_payload()
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.reg.args$regtype.engine <- "lp"
      hot.reg.args$degree.engine <- degree
      hot.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis
      hot.opt.args <- opt.args
      hot.opt.args$nmulti <- .np_nomad_powell_hotstart_nmulti("single_iteration")
      powell.start <- proc.time()[3L]
      hot.payload <- .np_nomad_with_powell_progress(
        degree,
        .npindexbw_run_fixed_degree(
          xdat = xdat,
          ydat = ydat,
          bws = bw.vec,
          template = template,
          reg.args = hot.reg.args,
          opt.args = hot.opt.args
        )
      )
      powell.elapsed <- proc.time()[3L] - powell.start
      direct.payload$num.feval <- as.numeric(nomad.num.feval.total) + as.numeric(hot.payload$num.feval[1L])
      direct.payload$num.feval.fast <- as.numeric(nomad.num.feval.fast.total) + as.numeric(hot.payload$num.feval.fast[1L])
      hot.payload$num.feval <- direct.payload$num.feval
      hot.payload$num.feval.fast <- direct.payload$num.feval.fast
      if (!is.null(hot.payload$method) && length(hot.payload$method))
        hot.payload$pmethod <- bwmToPrint(as.character(hot.payload$method[1L]))
      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = "min")) {
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  search.result <- .np_nomad_search(
    engine = degree.search$engine,
    baseline_record = baseline.record,
    start_degree = degree.search$start.degree,
    x0 = x0,
    bbin = bbin,
    lb = lb,
    ub = ub,
    eval_fun = eval_fun,
    build_payload = build_payload,
    direction = "min",
    objective_name = "fval",
    nmulti = nomad.nmulti,
    nomad.inner.nmulti = nomad.inner.nmulti,
    random.seed = if (!is.null(opt.args$random.seed)) opt.args$random.seed else 42L,
    degree_spec = list(
      initial = degree.search$start.degree,
      lower = degree.search$lower,
      upper = degree.search$upper,
      basis = degree.search$basis,
      nobs = degree.search$nobs,
      user_supplied = degree.search$start.user
    )
  )

  if (fixed.nomad) {
    if (!is.null(search.result$restart.starts)) {
      search.result$restart.starts <- lapply(search.result$restart.starts, point_to_public)
    }
    if (!is.null(search.result$restart.bandwidth.starts)) {
      search.result$restart.bandwidth.starts <- lapply(search.result$restart.bandwidth.starts, bandwidth_point_to_public)
    }
    if (!is.null(search.result$best_point) && length(search.result$best_point)) {
      search.result$best_point <- point_to_public(search.result$best_point)
    }
    if (!is.null(search.result$restart.results) && length(search.result$restart.results)) {
      search.result$restart.results <- lapply(search.result$restart.results, function(rec) {
        if (!is.null(rec$start) && length(rec$start))
          rec$start <- point_to_public(rec$start)
        if (!is.null(rec$solution) && length(rec$solution))
          rec$solution <- point_to_public(rec$solution)
        rec
      })
    }
  }

  search.result
}

.npindexbw_degree_search_controls <- function(regtype,
                                              regtype.named,
                                              bandwidth.compute,
                                              nobs,
                                              basis,
                                              degree.select,
                                              search.engine,
                                              degree.min,
                                              degree.max,
                                              degree.start,
                                              degree.restarts,
                                              degree.max.cycles,
                                              degree.verify,
                                              bernstein.basis,
                                              bernstein.named) {
  degree.select <- match.arg(degree.select, c("manual", "coordinate", "exhaustive"))
  if (identical(degree.select, "manual"))
    return(NULL)
  search.engine <- .np_degree_search_engine_controls(search.engine)

  regtype.requested <- if (isTRUE(regtype.named)) match.arg(regtype, c("lc", "ll", "lp")) else "lc"
  if (!identical(regtype.requested, "lp"))
    stop("automatic degree search currently requires regtype='lp'")
  if (!isTRUE(bandwidth.compute))
    stop("automatic degree search requires bandwidth.compute=TRUE")

  bern.auto <- if (isTRUE(bernstein.named)) bernstein.basis else TRUE
  bern.auto <- npValidateGlpBernstein(regtype = "lp", bernstein.basis = bern.auto)

  bounds <- .np_degree_normalize_bounds(
    ncon = 1L,
    degree.min = degree.min,
    degree.max = degree.max,
    default.max = 3L
  )

  baseline.degree <- 0L
  default.start.degree <- if (identical(search.engine, "cell")) baseline.degree else 1L
  start.degree <- if (is.null(degree.start)) {
    pmax(bounds$lower, pmin(bounds$upper, default.start.degree))
  } else {
    start.raw <- npValidateGlpDegree(
      regtype = "lp",
      degree = degree.start,
      ncon = 1L,
      argname = "degree.start"
    )
    if (!(start.raw[1L] %in% bounds$candidates[[1L]]))
      stop("degree.start must lie within the searched degree candidates")
    start.raw
  }

  list(
    method = if (identical(search.engine, "cell")) degree.select else search.engine,
    engine = search.engine,
    candidates = bounds$candidates,
    lower = bounds$lower,
    upper = bounds$upper,
    baseline.degree = as.integer(baseline.degree),
    start.degree = as.integer(start.degree),
    start.user = !is.null(degree.start),
    basis = if (missing(basis) || is.null(basis)) "glp" else as.character(basis[1L]),
    nobs = as.integer(nobs[1L]),
    restarts = npValidateNonNegativeInteger(degree.restarts, "degree.restarts"),
    max.cycles = npValidatePositiveInteger(degree.max.cycles, "degree.max.cycles"),
    verify = npValidateScalarLogical(degree.verify, "degree.verify"),
    bernstein.basis = bern.auto
  )
}

.npindexbw_attach_degree_search <- function(bws, search_result) {
  metadata <- list(
    mode = search_result$method,
    verify = isTRUE(search_result$verify),
    completed = isTRUE(search_result$completed),
    certified = isTRUE(search_result$certified),
    interrupted = isTRUE(search_result$interrupted),
    baseline.degree = search_result$baseline$degree,
    baseline.fval = search_result$baseline$objective,
    best.degree = search_result$best$degree,
    best.fval = search_result$best$objective,
    nomad.time = search_result$nomad.time,
    powell.time = search_result$powell.time,
    optim.time = search_result$optim.time,
    n.unique = search_result$n.unique,
    n.visits = search_result$n.visits,
    n.cached = search_result$n.cached,
    grid.size = search_result$grid.size,
    best.restart = search_result$best.restart,
    restart.starts = search_result$restart.starts,
    restart.degree.starts = search_result$restart.degree.starts,
    restart.bandwidth.starts = search_result$restart.bandwidth.starts,
    restart.start.info = search_result$restart.start.info,
    restart.results = search_result$restart.results,
    trace = search_result$trace
  )

  if (!is.null(search_result$nomad.time))
    bws$nomad.time <- as.numeric(search_result$nomad.time[1L])
  if (!is.null(search_result$powell.time))
    bws$powell.time <- as.numeric(search_result$powell.time[1L])
  if (!is.null(search_result$optim.time) && is.finite(search_result$optim.time))
    bws$total.time <- as.numeric(search_result$optim.time[1L])
  if (is.null(bws$timing.profile) && is.list(search_result$best$timing.profile))
    bws$timing.profile <- search_result$best$timing.profile
  bws <- .np_attach_nomad_restart_summary(bws, search_result)
  bws$degree.search <- metadata
  bws
}

npindexbw.default <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws,
           bandwidth.compute = TRUE,
           basis = c("glp", "additive", "tensor"),
           bernstein.basis = FALSE,
           degree = NULL,
           degree.select = c("manual", "coordinate", "exhaustive"),
           search.engine = c("nomad+powell", "cell", "nomad"),
           nomad = FALSE,
           nomad.nmulti = 0L,
           degree.min = NULL,
           degree.max = NULL,
           degree.start = NULL,
           degree.restarts = 0L,
           degree.max.cycles = 20L,
           degree.verify = FALSE,
           nmulti,
           only.optimize.beta,
           optim.abstol,
           optim.maxattempts,
           optim.maxit,
           optim.method,
           optim.reltol,
           random.seed,
           regtype = c("lc", "ll", "lp"),
           scale.factor.init.lower = 0.1,
           scale.factor.init.upper = 2.0,
           scale.factor.init = 0.5,
           scale.factor.search.lower = NULL,
           ...){
    .npRmpi_require_active_slave_pool(where = "npindexbw()")
    search.mc.names <- names(match.call(expand.dots = FALSE))
    degree.select.value <- if (isTRUE(npValidateScalarLogical(nomad, "nomad"))) {
      "coordinate"
    } else if ("degree.select" %in% search.mc.names) {
      degree.select
    } else {
      "manual"
    }
    automatic.degree.search <- !identical(match.arg(degree.select.value, c("manual", "coordinate", "exhaustive")), "manual")
    if (.npRmpi_autodispatch_active() && !isTRUE(automatic.degree.search))
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    xdat <- toFrame(xdat)

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (ncol(xdat) < 2) {
      if (coarseclass(xdat[,1]) != "numeric")
        stop("xdat must contain at least one continuous variable")

      .np_warning(paste("xdat has one dimension. Using a single index model to reduce",
                    "dimensionality is unnecessary."))
    }

    if (coarseclass(bws) != "numeric" || length(bws) != ncol(xdat)+1)
      stop(paste("manually specified 'bws' must be a numeric vector of length ncol(xdat)+1.",
                 "See documentation for details."))

    p <- ncol(xdat)
    mc <- match.call(expand.dots = FALSE)
    mc.names <- names(mc)
    dots <- list(...)
    dot.names <- names(dots)
    nomad.shortcut <- .np_prepare_nomad_shortcut(
      nomad = nomad,
      call_names = unique(c(mc.names, dot.names)),
      preset = list(
        regtype = "lp",
        search.engine = "nomad+powell",
        degree.select = "coordinate",
        bernstein.basis = TRUE,
        degree.min = 0L,
        degree.max = 10L,
        degree.verify = FALSE,
        bwtype = "fixed"
      ),
      values = list(
        regtype = if ("regtype" %in% mc.names) regtype else NULL,
        search.engine = if ("search.engine" %in% mc.names) search.engine else NULL,
        degree.select = if ("degree.select" %in% mc.names) degree.select else NULL,
        bernstein.basis = if ("bernstein.basis" %in% mc.names) bernstein.basis else NULL,
        degree.min = if ("degree.min" %in% mc.names) degree.min else NULL,
        degree.max = if ("degree.max" %in% mc.names) degree.max else NULL,
        degree.verify = if ("degree.verify" %in% mc.names) degree.verify else NULL,
        bwtype = if ("bwtype" %in% dot.names) dots$bwtype else NULL,
        degree = if ("degree" %in% mc.names) degree else NULL
      ),
      where = "npindexbw"
    )

    if (isTRUE(nomad.shortcut$enabled)) {
      if ("degree" %in% mc.names)
        stop("nomad=TRUE does not support an explicit degree; remove degree or set nomad=FALSE")
      if ("regtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp")))[1L], "lp"))
        stop("nomad=TRUE requires regtype='lp'")
      if ("bwtype" %in% dot.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$bwtype, c("fixed", "generalized_nn", "adaptive_nn")))[1L], "fixed"))
        stop("nomad=TRUE currently requires bwtype='fixed'")
      if ("degree.select" %in% mc.names &&
          identical(as.character(match.arg(nomad.shortcut$values$degree.select, c("manual", "coordinate", "exhaustive")))[1L], "manual"))
        stop("nomad=TRUE requires automatic degree search; use degree.select='coordinate' or 'exhaustive'")
      if ("search.engine" %in% mc.names &&
          !(as.character(match.arg(nomad.shortcut$values$search.engine, c("nomad+powell", "cell", "nomad")))[1L] %in%
              c("nomad", "nomad+powell")))
        stop("nomad=TRUE requires search.engine='nomad' or 'nomad+powell'")
      if ("bernstein.basis" %in% mc.names &&
          !isTRUE(npValidateGlpBernstein(regtype = "lp",
                                        bernstein.basis = nomad.shortcut$values$bernstein.basis)))
        stop("nomad=TRUE currently requires bernstein.basis=TRUE")
      if ("degree.verify" %in% mc.names &&
          isTRUE(npValidateScalarLogical(nomad.shortcut$values$degree.verify, "degree.verify")))
        stop("nomad=TRUE currently requires degree.verify=FALSE")
    }

    random.seed.value <- if ("random.seed" %in% mc.names) {
      npValidateNonNegativeInteger(random.seed, "random.seed")
    } else {
      42L
    }
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else degree.select.value
    degree.search <- .npindexbw_degree_search_controls(
      regtype = if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else regtype,
      regtype.named = isTRUE(nomad.shortcut$enabled) || ("regtype" %in% mc.names),
      bandwidth.compute = bandwidth.compute,
      nobs = NROW(xdat),
      basis = basis,
      degree.select = degree.select.value,
      search.engine = if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell",
      degree.min = nomad.shortcut$values$degree.min,
      degree.max = nomad.shortcut$values$degree.max,
      degree.start = if ("degree.start" %in% mc.names) degree.start else NULL,
      degree.restarts = if ("degree.restarts" %in% mc.names) degree.restarts else 0L,
      degree.max.cycles = if ("degree.max.cycles" %in% mc.names) degree.max.cycles else 20L,
      degree.verify = if (!is.null(nomad.shortcut$values$degree.verify)) nomad.shortcut$values$degree.verify else FALSE,
      bernstein.basis = if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else bernstein.basis,
      bernstein.named = isTRUE(nomad.shortcut$enabled) || ("bernstein.basis" %in% mc.names)
    )
    nomad.inner <- .np_nomad_validate_inner_multistart(
      call_names = mc.names,
      dot.args = dots,
      nomad.nmulti = nomad.nmulti,
      regtype = if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else regtype,
      automatic.degree.search = !is.null(degree.search),
      search.engine = if (is.null(degree.search)) "" else degree.search$engine
    )
    nomad.inner.named <- nomad.inner$named
    nomad.inner.nmulti <- nomad.inner$nmulti
    degree.setup <- npSetupGlpDegree(
      regtype = if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else regtype,
      degree = degree,
      ncon = 1L,
      degree.select = degree.select.value
    )
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(scale.factor.search.lower)
    spec <- npResolveCanonicalConditionalRegSpec(
      mc.names = mc.names,
      regtype = if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else regtype,
      basis = basis,
      degree = degree.setup,
      bernstein.basis = if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else bernstein.basis,
      ncon = 1L,
      where = "npindexbw"
    )
    if (!is.null(degree.search)) {
      spec$bernstein.basis <- degree.search$bernstein.basis
      spec$bernstein.basis.engine <- degree.search$bernstein.basis
    }
    tbw <- sibandwidth(beta = bws[seq_len(p)],
                       h = bws[p+1L], ...,
                       regtype = spec$regtype,
                       basis = spec$basis,
                       degree = spec$degree,
                       bernstein.basis = spec$bernstein.basis,
                       nobs = dim(xdat)[1],
                       xdati = untangle(xdat),
                       ydati = untangle(data.frame(ydat)),
                       xnames = names(xdat),
                       ynames = deparse(substitute(ydat)),
                       bandwidth = bws[p+1L],
                       bandwidth.compute = bandwidth.compute)
    tbw <- npSetScaleFactorSearchLower(tbw, scale.factor.search.lower)

    if (tbw$method == "kleinspady" && !setequal(ydat,c(0,1)))
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    margs <- c("nmulti","random.seed", "optim.method", "optim.maxattempts",
               "optim.reltol", "optim.abstol", "optim.maxit", "only.optimize.beta",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "scale.factor.search.lower")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)
    bwsel.args <- list(xdat = xdat, ydat = ydat, bws = tbw)
    if (any.m) {
      nms <- mc.names[m]
      bwsel.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }

    if (bandwidth.compute && !is.null(degree.search)) {
      reg.args <- list(
        regtype = spec$regtype,
        basis = spec$basis,
        degree = spec$degree,
        bernstein.basis = spec$bernstein.basis,
        regtype.engine = spec$regtype.engine,
        basis.engine = spec$basis.engine,
        degree.engine = spec$degree.engine,
        bernstein.basis.engine = spec$bernstein.basis.engine,
        scale.factor.search.lower = scale.factor.search.lower
      )
      opt.args <- c(
        list(bandwidth.compute = bandwidth.compute),
        bwsel.args[setdiff(names(bwsel.args), c("xdat", "ydat", "bws"))]
      )
      opt.args$scale.factor.search.lower <- scale.factor.search.lower

      if (identical(degree.search$engine, "cell")) {
        eval_fun <- function(degree.vec) {
          cell.reg.args <- reg.args
          cell.reg.args$regtype <- "lp"
          cell.reg.args$degree <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
          cell.reg.args$regtype.engine <- "lp"
          cell.reg.args$degree.engine <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis.engine <- degree.search$bernstein.basis
          cell.bws <- .npindexbw_run_fixed_degree(
            xdat = xdat,
            ydat = ydat,
            bws = bws,
            template = tbw,
            reg.args = cell.reg.args,
            opt.args = opt.args
          )
          list(
            objective = as.numeric(cell.bws$fval[1L]),
            payload = cell.bws,
            num.feval = if (!is.null(cell.bws$num.feval)) as.numeric(cell.bws$num.feval[1L]) else NA_real_
          )
        }

        search.result <- .np_degree_search(
          method = degree.search$method,
          candidates = degree.search$candidates,
          baseline_degree = degree.search$baseline.degree,
          start_degree = degree.search$start.degree,
          restarts = degree.search$restarts,
          max_cycles = degree.search$max.cycles,
          verify = degree.search$verify,
          eval_fun = eval_fun,
          direction = "min",
          trace_level = "full",
          objective_name = "fval"
        )
      } else {
        search.result <- .npindexbw_nomad_search(
          xdat = xdat,
          ydat = ydat,
          bws = bws,
          template = tbw,
          reg.args = reg.args,
          opt.args = utils::modifyList(opt.args, list(random.seed = random.seed.value)),
          degree.search = degree.search,
          nomad.inner.nmulti = nomad.inner.nmulti
        )
      }
      tbw <- .npindexbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
    } else if (bandwidth.compute) {
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting single-index bandwidth",
        do.call(npindexbw.sibandwidth, bwsel.args)
      )
    }

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    return(tbw)
  }

npindexbw.sibandwidth <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws,
           bandwidth.compute = TRUE,
           nmulti,
           only.optimize.beta = FALSE,
           optim.abstol = .Machine$double.eps,
           optim.maxattempts = 10,
           optim.maxit = 500,
           optim.method = c("Nelder-Mead", "BFGS", "CG"),
           optim.reltol = sqrt(.Machine$double.eps),
           random.seed = 42,
           scale.factor.init.lower = 0.1,
           scale.factor.init.upper = 2.0,
           scale.factor.init = 0.5,
           scale.factor.search.lower = NULL,
           ...){
    ## Save seed prior to setting

    seed.state <- .np_seed_enter(random.seed)


    xdat = toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(ncol(xdat))
    }
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    only.optimize.beta <- npValidateScalarLogical(only.optimize.beta, "only.optimize.beta")
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)
    optim.maxattempts <- npValidatePositiveInteger(optim.maxattempts, "optim.maxattempts")
    optim.maxit <- npValidatePositiveInteger(optim.maxit, "optim.maxit")
    optim.reltol <- npValidatePositiveFiniteNumeric(optim.reltol, "optim.reltol")
    optim.abstol <- npValidatePositiveFiniteNumeric(optim.abstol, "optim.abstol")
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
    )
    h.start.controls <- .npindexbw_h_start_controls(
      scale.factor.init.lower = scale.factor.init.lower,
      scale.factor.init.upper = scale.factor.init.upper,
      scale.factor.init = scale.factor.init,
      scale.factor.search.lower = scale.factor.search.lower,
      where = "npindexbw"
    )
    .npRmpi_require_active_slave_pool(where = "npindexbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    if (bws$method == "kleinspady" && !setequal(ydat,c(0,1)))
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    if (ncol(xdat) < 2) {
      if (coarseclass(xdat[,1]) != "numeric")
        stop("xdat must contain at least one continuous variable")

      .np_warning(paste("xdat has one dimension. Using a single index model to reduce",
                    "dimensionality is unnecessary."))
    }

    optim.method <- match.arg(optim.method)

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    rows.omit <- attr(na.omit(data.frame(xdat,ydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]

    ## convert to numeric
    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat = toMatrix(xdat)
    p <- ncol(xdat)
    beta.idx <- if (p > 1L) seq_len(p - 1L) else integer(0)
    nobs <- nrow(xdat)
    spec <- .npindex_resolve_spec(bws, where = "npindexbw")

    total.time <-
      system.time({

        if(bandwidth.compute){

          ## Invariant objects used by objective evaluations.
          xmat <- xdat
          wmat <- cbind(ydat, 1.0)
          bandwidth_eval_count <- 0L

          bandwidth_progress_step <- function() {
            bandwidth_eval_count <<- bandwidth_eval_count + 1L
            .np_progress_bandwidth_activity_step(done = bandwidth_eval_count)
            invisible(NULL)
          }
          fixed.h.lower <- NULL

          ## Note - there are two methods currently implemented, Ichimura's
          ## least squares approach and Klein and Spady's likelihood approach.

          ##We define ichimura's objective function. Since we normalize beta_1
          ##to equal 1, we only worry about beta_2...beta_k in the index
          ##function for these internals, and when computing the leave-one-out
          ##objective function use c(1,beta). However, we do indeed return
          ##c(1,beta) which can be used in the index.model function above.

          ichimuraMaxPenalty <- 10*mean(ydat^2)

          ichimura <- function(param) {
            bandwidth_progress_step()

            ##Define the leave-one-out objective function, sum (y - \hat
            ## G(X\hat\beta))^2. We let beta denote beta_2...beta_k (first k-1
            ## parameters in `param') and then let h denote the kth column.
            beta <- param[beta.idx]
            h <- param[p]
            h.candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bws$type, nobs = nobs)
            if (!h.candidate$ok)
              return(ichimuraMaxPenalty)
            h <- h.candidate$value
            if (!is.null(fixed.h.lower) && h < fixed.h.lower)
              return(ichimuraMaxPenalty)

            ## Next we define the sum of squared leave-one-out residuals

            sum.squares.leave.one.out <- function(beta,h) {

              ## Normalize beta_1 = 1 hence multiply X by c(1,beta)

              index <- xmat %*% c(1, beta)

              fit.loo <- if (identical(spec$regtype.engine, "lc")) {
                ## One call to npksum to avoid repeated computation of the
                ## product kernel (the expensive part)
                tww <- .npindex_selector_with_local_kernelsum(
                  npksum(txdat=index,
                         tydat=wmat,
                         weights=wmat,
                         leave.one.out=TRUE,
                         bandwidth.divide=TRUE,
                         bws=c(h),
                         bwtype = bws$type,
                         ckertype = bws$ckertype,
                         ckerorder = bws$ckerorder,
                         ckerbound = bws$ckerbound,
                         ckerlb = bws$ckerlb,
                         ckerub = bws$ckerub)
                )$ksum
                tww[1,2,]/NZD(tww[2,2,])
              } else {
                ok.design <- tryCatch({
                  npCheckRegressionDesignCondition(
                    reg.code = REGTYPE_LP,
                    xcon = data.frame(index = index),
                    basis = spec$basis.engine,
                    degree = spec$degree.engine,
                    bernstein.basis = spec$bernstein.basis.engine,
                    where = "npindexbw"
                  )
                  TRUE
                }, error = function(e) FALSE)
                if (!ok.design)
                  return(ichimuraMaxPenalty)

                .npindex_lp_loo_fit(
                  index = index,
                  ydat = ydat,
                  h = h,
                  bws = bws,
                  spec = spec
                )
              }

              t.ret <- mean((ydat-fit.loo)^2)
              if (is.finite(t.ret) && .npindexbw_fast_eligible(h = h, bws = bws, eval.index = index))
                num.feval.fast.overall <<- num.feval.fast.overall + 1L
              return(t.ret)

            }

            ## For the objective function, we require a positive bandwidth, so
            ## return an infinite penalty for negative h

            if(h > 0) {
              return(sum.squares.leave.one.out(beta,h))
            } else {
              return(ichimuraMaxPenalty)
            }

          }

          ichimura.nobw <- function(param,h){ return(ichimura(c(param,h))) }

          ## We define ichimura's objective function. Since we normalize beta_1
          ## to equal 1, we only worry about beta_2...beta_k in the index
          ## function for these internals, and when computing the leave-one-out
          ## objective function use c(1,beta). However, we do indeed return
          ## c(1,beta) which can be used in the index.model function above.

          kleinspadyFloor <- sqrt(.Machine$double.eps)

          kleinspady <- function(param) {
            bandwidth_progress_step()

            ## Define the leave-one-out objective function, sum (y - \hat
            ## G(X\hat\beta))^2. We let beta denote beta_2...beta_k (first k-1
            ## parameters in `param') and then let h denote the kth column.
            beta <- param[beta.idx]
            h <- param[p]
            h.candidate <- .npindex_nn_candidate_bandwidth(h = h, bwtype = bws$type, nobs = nobs)
            if (!h.candidate$ok)
              return(sqrt(.Machine$double.xmax))
            h <- h.candidate$value
            if (!is.null(fixed.h.lower) && h < fixed.h.lower)
              return(sqrt(.Machine$double.xmax))

            ## Next we define the sum of logs

            sum.log.leave.one.out <- function(beta,h) {

              ## Normalize beta_1 = 1 hence multiply X by c(1,beta)

              index <- xmat %*% c(1, beta)

              ks.loo <- if (identical(spec$regtype.engine, "lc")) {
                ## One call to npksum to avoid repeated computation of the
                ## product kernel (the expensive part)
                tww <- .npindex_selector_with_local_kernelsum(
                  npksum(txdat=index,
                         tydat=wmat,
                         weights=wmat,
                         leave.one.out=TRUE,
                         bandwidth.divide=TRUE,
                         bws=c(h),
                         bwtype = bws$type,
                         ckertype = bws$ckertype,
                         ckerorder = bws$ckerorder,
                         ckerbound = bws$ckerbound,
                         ckerlb = bws$ckerlb,
                         ckerub = bws$ckerub)
                )$ksum
                tww[1,2,]/NZD(tww[2,2,])
              } else {
                ok.design <- tryCatch({
                  npCheckRegressionDesignCondition(
                    reg.code = REGTYPE_LP,
                    xcon = data.frame(index = index),
                    basis = spec$basis.engine,
                    degree = spec$degree.engine,
                    bernstein.basis = spec$bernstein.basis.engine,
                    where = "npindexbw"
                  )
                  TRUE
                }, error = function(e) FALSE)
                if (!ok.design)
                  return(sqrt(.Machine$double.xmax))

                .npindex_lp_loo_fit(
                  index = index,
                  ydat = ydat,
                  h = h,
                  bws = bws,
                  spec = spec
                )
              }

              ## Avoid taking log of zero (ks.loo = 0 or 1 since we take
              ## the log of ks.loo and the log of 1-ks.loo)

              ks.loo[ks.loo < kleinspadyFloor] <- kleinspadyFloor
              ks.loo[ks.loo > 1 - kleinspadyFloor] <- 1 - kleinspadyFloor

              ## Maximize the log likelihood, therefore minimize minus.
              ## Here ydat is binary (0/1), so this is equivalent to
              ## ifelse(ydat==1, log(ks.loo), log(1-ks.loo)) but avoids
              ## branchy ifelse and uses stable log1p for the second term.
              t.ret <- -mean(ydat * log(ks.loo) + (1.0 - ydat) * log1p(-ks.loo))
              if (is.finite(t.ret) && .npindexbw_fast_eligible(h = h, bws = bws, eval.index = index))
                num.feval.fast.overall <<- num.feval.fast.overall + 1L
              return(t.ret)

            }

            ## For the objective function, we require a positive bandwidth, so
            ## return an infinite penalty for negative h

            if(h > 0) {
              return(sum.log.leave.one.out(beta,h))
            } else {
              ## No natural counterpart to var of y here, unlike Ichimura above...
              return(sqrt(.Machine$double.xmax))
            }

          }

          kleinspady.nobw <- function(param,h){ return(kleinspady(c(param,h))) }
          ## Now we implement multistarting

          fval.min <- .Machine$double.xmax
          numimp <- 0
          fval.value <- numeric(nmulti)
          num.feval.overall <- 0
          num.feval.fast.overall <- 0L

          if(bws$method == "ichimura"){
            optim.fn <- if(only.optimize.beta) ichimura.nobw else ichimura
            optim.control <- list(abstol=optim.abstol,
                                  reltol=optim.reltol,
                                  maxit=optim.maxit)
          } else if(bws$method == "kleinspady"){
            optim.fn <- if(only.optimize.beta) kleinspady.nobw else  kleinspady
            optim.control <- list(reltol=optim.reltol,maxit=optim.maxit)
          }

          for (i in seq_len(nmulti)) {
            ## We use the nlm command to minimize the objective function using
            ## starting values. Note that since we normalize beta_1=1 here beta
            ## is the k-1 vector containing beta_2...beta_k

            if(i == 1) {

              ## Initial values taken from OLS fit with a constant used for
              ## multistart 1
              ols.fit <- lm(ydat~xdat,x=TRUE)
              fit <- fitted(ols.fit)
              fixed.h.lower <- if (identical(bws$type, "fixed")) {
                h.start.controls$scale.factor.search.lower * .npindex_start_bandwidth_scale(fit = fit, nobs = nobs)
              } else {
                NULL
              }

              if (p != 1L){
                if (setequal(bws$beta[2:p], c(0)))
                  beta <- coef(ols.fit)[3:ncol(ols.fit$x)]
                else
                  beta = bws$beta[2:p]
              } else { beta = numeric(0) }

              if (bws$bw == 0)
                h <- .npindex_default_start_bandwidth(
                  fit = fit,
                  bwtype = bws$type,
                  nobs = nobs,
                  start.controls = h.start.controls
                )
              else
                h <- tryCatch(
                  .npindex_finalize_bandwidth(
                    h = bws$bw,
                    bwtype = bws$type,
                    nobs = nobs,
                    lower = fixed.h.lower,
                    where = "npindexbw"
                  ),
                  error = function(e) .npindex_default_start_bandwidth(
                    fit = fit,
                    bwtype = bws$type,
                    nobs = nobs,
                    start.controls = h.start.controls
                  )
                )
            } else {
              ## Random initialization used for remaining multistarts

              beta.length <- length(coef(ols.fit)[3:ncol(ols.fit$x)])
              beta <- runif(beta.length,min=0.5,max=1.5)*coef(ols.fit)[3:ncol(ols.fit$x)]
              if (!only.optimize.beta)
                h <- .npindex_random_start_bandwidth(
                  fit = fit,
                  bwtype = bws$type,
                  nobs = nobs,
                  start.controls = h.start.controls
                )
            }

            optim.parm <- if(only.optimize.beta) beta else c(beta,h)

            optim.base.args <- list(
              par = optim.parm,
              fn = optim.fn,
              gr = NULL,
              method = optim.method,
              control = optim.control
            )
            if (only.optimize.beta) {
              optim.base.args$h <- h
            }
            suppressWarnings(optim.return <- do.call(optim, optim.base.args))
            if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
              num.feval.overall <- num.feval.overall + optim.return$counts[1]
            attempts <- 0
            while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
              attempts <- attempts + 1
              beta.length <- length(coef(ols.fit)[3:ncol(ols.fit$x)])
              beta <- runif(beta.length,min=0.5,max=1.5)*coef(ols.fit)[3:ncol(ols.fit$x)]
              if(!only.optimize.beta)
                h <- .npindex_random_start_bandwidth(
                  fit = fit,
                  bwtype = bws$type,
                  nobs = nobs,
                  start.controls = h.start.controls
                )

              if(optim.return$convergence == 1){
                if(optim.control$maxit < (2^32/10))
                  optim.control$maxit <- 10*optim.control$maxit
                else
                  stop(paste("optim failed to converge after optim.maxattempts = ", optim.maxattempts, " iterations."))
              }

              if(optim.return$convergence == 10){
                optim.control$reltol <-  10.0*optim.control$reltol
                if(!is.null(optim.control$abstol))
                  optim.control$abstol <-  10.0*optim.control$abstol
              }

              optim.parm <- if(only.optimize.beta) beta else c(beta,h)
              optim.base.args$par <- optim.parm
              optim.base.args$control <- optim.control
              if (!only.optimize.beta && ("h" %in% names(optim.base.args))) {
                optim.base.args$h <- NULL
              }
              if (only.optimize.beta) {
                optim.base.args$h <- h
              }
              suppressWarnings(optim.return <- do.call(optim, optim.base.args))
              if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                num.feval.overall <- num.feval.overall + optim.return$counts[1]
            }

            if(optim.return$convergence != 0)
              stop(paste("optim failed to converge after optim.maxattempts = ", optim.maxattempts, " iterations."))

            fval.value[i] <- optim.return$value
            if(optim.return$value < fval.min) {
              param <- if(only.optimize.beta) c(optim.return$par,h) else optim.return$par
              fval.min <- optim.return$value
              numimp <- numimp + 1
              best <- i
            }

            .np_progress_bandwidth_multistart_step(done = i, total = nmulti)

          }

          bws$beta <- c(1.0, param[beta.idx])
          bws$bw <- .npindex_finalize_bandwidth(
            h = param[p],
            bwtype = bws$type,
            nobs = nobs,
            lower = fixed.h.lower,
            where = "npindexbw"
          )
          bws$fval <- fval.min
          bws$ifval <- best
          bws$num.feval <- num.feval.overall
          bws$num.feval.fast <- num.feval.fast.overall
          bws$numimp <- numimp
          bws$fval.vector <- fval.value
        }
      })[["elapsed"]]
    ## Return a list with beta (we append the restricted value of
    ## beta_1=1), the bandwidth h, the value of the objective function at
    ## its minimum, the number of restarts that resulted in an improved
    ## value of the objective function, the restart that resulted in the
    ## smallest value, and the vector of objective function values.

    ## Restore seed

    .np_seed_exit(seed.state)

    bws <- sibandwidth(beta = bws$beta,
                       h = bws$bw,
                       method = bws$method,
                       regtype = bws$regtype,
                       basis = bws$basis,
                       degree = bws$degree,
                       bernstein.basis = bws$bernstein.basis,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       ckerbound = bws$ckerbound,
                       ckerlb = bws$ckerlb,
                       ckerub = bws$ckerub,
                       bwtype = bws$type,
                       fval = bws$fval,
                       ifval = bws$ifval,
                       num.feval = bws$num.feval,
                       num.feval.fast = bws$num.feval.fast,
                       numimp = bws$numimp,
                       fval.vector = bws$fval.vector,
                       nobs = bws$nobs,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       bandwidth = bws$bw,
                       rows.omit = rows.omit,
                       bandwidth.compute = bandwidth.compute,
                       optim.method = optim.method,
                       only.optimize.beta = only.optimize.beta,
                       total.time = total.time)
    bws <- npSetScaleFactorSearchLower(bws, scale.factor.search.lower)

    bws

  }
