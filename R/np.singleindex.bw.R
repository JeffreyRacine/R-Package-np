
npindexbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
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

.npindex_default_start_bandwidth <- function(fit, bwtype, nobs) {
  if (identical(bwtype, "fixed")) {
    return(1.059224 * EssDee(fit) * nobs^(-1 / 5))
  }

  lower <- 2L
  max(lower, min(max(1L, as.integer(nobs) - 1L), .np_round_half_to_even(sqrt(nobs))))
}

.npindex_random_start_bandwidth <- function(fit, bwtype, nobs) {
  if (identical(bwtype, "fixed")) {
    return(runif(1, min = 0.5, max = 1.5) * EssDee(fit) * nobs^(-1 / 5))
  }

  upper <- max(1L, as.integer(nobs) - 1L)
  runif(1, min = 2, max = max(2L, upper))
}

.npindex_finalize_bandwidth <- function(h, bwtype, nobs, where = "npindexbw") {
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

  candidate$value
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
    return(as.numeric(invalid.penalty))

  as.numeric(out$objective[1L])
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
      blocks <- Filter(length, .splitIndices(n, min(as.integer(n), as.integer(slave.num))))
      pieces <- mpi.apply(
        blocks,
        .npindex_lp_loo_fit_block,
        index = index,
        ydat = ydat,
        h = h,
        bws = bws,
        spec = spec,
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
  do.call(
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
    return(invalid.penalty)
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
      )$ksum,
      error = function(e) NULL
    )
    if (is.null(tww))
      return(invalid.penalty)
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
      return(invalid.penalty)

    .npindex_lp_loo_fit(
      index = index,
      ydat = ydat,
      h = h,
      bws = bws,
      spec = spec
    )
  }

  if (any(!is.finite(fit.loo)))
    return(invalid.penalty)

  if (identical(bws$method, "ichimura")) {
    return(mean((ydat - fit.loo)^2))
  }

  floor <- sqrt(.Machine$double.eps)
  fit.loo[fit.loo < floor] <- floor
  fit.loo[fit.loo > 1 - floor] <- 1 - floor
  -mean(ydat * log(fit.loo) + (1.0 - ydat) * log1p(-fit.loo))
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
  beta.free <- if (p > 1L) seq_len(p - 1L) else integer(0)
  beta.start <- if (length(beta.free)) as.double(baseline.bws$beta[beta.free + 1L]) else numeric(0)
  h.start <- as.double(baseline.bws$bw[1L])
  beta.bound <- if (length(beta.start)) pmax(10, 10 * abs(beta.start)) else numeric(0)
  h.lower <- if (identical(baseline.bws$type, "fixed")) 1e-3 else 2
  h.upper <- if (identical(baseline.bws$type, "fixed")) {
    max(1e6, abs(h.start) * 1e3, 1)
  } else {
    max(2L, as.integer(nrow(x.clean)) - 1L)
  }
  h.integer <- !identical(baseline.bws$type, "fixed")

  x0 <- c(beta.start, h.start, as.integer(degree.search$start.degree))
  lb <- c(-beta.bound, h.lower, degree.search$lower)
  ub <- c(beta.bound, h.upper, degree.search$upper)
  bbin <- c(rep.int(0L, length(beta.start)), if (isTRUE(h.integer)) 1L else 0L, 1L)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(ncol(xdat)) else max(1L, as.integer(opt.args$nmulti[1L]))
  baseline.record <- NULL

  .np_nomad_baseline_note(degree.search$start.degree)

  point_to_bws <- function(point) {
    beta.tail <- if (length(beta.free)) as.double(point[seq_along(beta.free)]) else numeric(0)
    h <- as.double(point[length(beta.free) + 1L])
    h <- .npindex_finalize_bandwidth(h = h, bwtype = baseline.bws$type, nobs = nrow(x.clean), where = "npindexbw")
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
      param = point_to_bws(point),
      xmat = x.clean,
      ydat = y.clean,
      bws = baseline.bws,
      spec = eval.spec
    )
    list(
      objective = as.numeric(objective),
      degree = as.integer(degree),
      num.feval = 1
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
      tbw$num.feval <- if (!is.null(solution$bbe)) as.numeric(solution$bbe) else as.numeric(best_record$num.feval)
      tbw$total.time <- NA_real_
      npindexbw.sibandwidth(
        xdat = xdat,
        ydat = ydat,
        bws = tbw,
        bandwidth.compute = FALSE
      )
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
      hot.opt.args$nmulti <- 1L
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
      hot.objective <- as.numeric(hot.payload$fval[1L])
      if (is.finite(hot.objective) &&
          .np_degree_better(hot.objective, direct.objective, direction = "min")) {
        return(list(payload = hot.payload, objective = hot.objective, powell.time = powell.elapsed))
      }
    }

    list(payload = direct.payload, objective = direct.objective, powell.time = powell.elapsed)
  }

  .np_nomad_search(
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

    if (tbw$method == "kleinspady" && !setequal(ydat,c(0,1)))
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    bwsel.args <- list(xdat = xdat, ydat = ydat, bws = tbw)
    if (!missing(nmulti)) bwsel.args$nmulti <- nmulti
    if (!missing(random.seed)) bwsel.args$random.seed <- random.seed
    if (!missing(optim.method)) bwsel.args$optim.method <- optim.method
    if (!missing(optim.maxattempts)) bwsel.args$optim.maxattempts <- optim.maxattempts
    if (!missing(optim.reltol)) bwsel.args$optim.reltol <- optim.reltol
    if (!missing(optim.abstol)) bwsel.args$optim.abstol <- optim.abstol
    if (!missing(optim.maxit)) bwsel.args$optim.maxit <- optim.maxit
    if (!missing(only.optimize.beta)) bwsel.args$only.optimize.beta <- only.optimize.beta

    if (bandwidth.compute && !is.null(degree.search)) {
      reg.args <- list(
        regtype = spec$regtype,
        basis = spec$basis,
        degree = spec$degree,
        bernstein.basis = spec$bernstein.basis,
        regtype.engine = spec$regtype.engine,
        basis.engine = spec$basis.engine,
        degree.engine = spec$degree.engine,
        bernstein.basis.engine = spec$bernstein.basis.engine
      )
      opt.args <- c(
        list(bandwidth.compute = bandwidth.compute),
        bwsel.args[setdiff(names(bwsel.args), c("xdat", "ydat", "bws"))]
      )

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
      opt.args <- bwsel.args
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting single-index bandwidth",
        do.call(npindexbw.sibandwidth, opt.args)
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
           ...){
    ## Save seed prior to setting

    seed.state <- .np_seed_enter(random.seed)


    xdat = toFrame(xdat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(ncol(xdat))
    }
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    only.optimize.beta <- npValidateScalarLogical(only.optimize.beta, "only.optimize.beta")
    nmulti <- npValidateNonNegativeInteger(nmulti, "nmulti")
    .np_progress_bandwidth_set_total(nmulti)
    optim.maxattempts <- npValidatePositiveInteger(optim.maxattempts, "optim.maxattempts")
    optim.maxit <- npValidatePositiveInteger(optim.maxit, "optim.maxit")
    optim.reltol <- npValidatePositiveFiniteNumeric(optim.reltol, "optim.reltol")
    optim.abstol <- npValidatePositiveFiniteNumeric(optim.abstol, "optim.abstol")
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

              if (p != 1L){
                if (setequal(bws$beta[2:p], c(0)))
                  beta <- coef(ols.fit)[3:ncol(ols.fit$x)]
                else
                  beta = bws$beta[2:p]
              } else { beta = numeric(0) }

              if (bws$bw == 0)
                h <- .npindex_default_start_bandwidth(fit = fit, bwtype = bws$type, nobs = nobs)
              else
                h <- .npindex_finalize_bandwidth(h = bws$bw, bwtype = bws$type, nobs = nobs)
            } else {
              ## Random initialization used for remaining multistarts

              beta.length <- length(coef(ols.fit)[3:ncol(ols.fit$x)])
              beta <- runif(beta.length,min=0.5,max=1.5)*coef(ols.fit)[3:ncol(ols.fit$x)]
              if (!only.optimize.beta)
                h <- .npindex_random_start_bandwidth(fit = fit, bwtype = bws$type, nobs = nobs)
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
                h <- .npindex_random_start_bandwidth(fit = fit, bwtype = bws$type, nobs = nobs)

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
          bws$bw <- .npindex_finalize_bandwidth(h = param[p], bwtype = bws$type, nobs = nobs)
          bws$fval <- fval.min
          bws$ifval <- best
          bws$num.feval <- num.feval.overall
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

    bws

  }
