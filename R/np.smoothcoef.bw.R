npscoefbw <-
  function(...){
    mc <- match.call(expand.dots = FALSE)
    npRejectRenamedScaleFactorSearchArgs(names(mc$...), where = "npscoefbw")
    target <- .np_bw_dispatch_target(dots = mc$...,
                                     data_arg_names = c("xdat", "ydat", "zdat"),
                                     eval_env = parent.frame())
    UseMethod("npscoefbw", target)
  }

npscoefbw.formula <-
  function(formula, data, subset, na.action, call, ...){
    orig.ts <- if (missing(data))
      .np_terms_ts_mask(terms_obj = terms(formula),
                        data = environment(formula),
                        eval_env = environment(formula))
    else .np_terms_ts_mask(terms_obj = terms(formula, data = data),
                           data = data,
                           eval_env = environment(formula))

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]

    formula.call <- .np_bw_formula_from_call(call_obj = call, eval_env = parent.frame())
    if (!is.null(formula.call))
      mf[[2]] <- formula.call

    mf[[1]] <- as.name("model.frame")

    formula.obj <- .np_bw_resolve_formula(formula_obj = formula,
                                        formula_call = formula.call,
                                        eval_env = parent.frame())
    chromoly <- explodePipe(formula.obj, env = environment(formula))

    bronze <- sapply(chromoly, paste, collapse = " + ")
    mf[["formula"]] <-
      as.formula(paste(bronze[1]," ~ ",
                       paste(bronze[2:length(bronze)],
                             collapse =" + ")),
                 env = environment(formula))

    mf[["formula"]] <- terms(mf[["formula"]])
    if(all(orig.ts)){
      args <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      attr(mf[["formula"]], "predvars") <- as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), args))))
    }else if(any(orig.ts)){
      arguments <- (as.list(attr(mf[["formula"]], "variables"))[-1])
      arguments.normal <- arguments[which(!orig.ts)]
      arguments.timeseries <- arguments[which(orig.ts)]

      ix <- sort(c(which(orig.ts),which(!orig.ts)),index.return = TRUE)$ix
      attr(mf[["formula"]], "predvars") <- bquote(.(as.call(c(quote(cbind),as.call(c(quote(as.data.frame),as.call(c(quote(ts.intersect), arguments.timeseries)))),arguments.normal,check.rows = TRUE)))[,.(ix)])
    }

    mf.args <- as.list(mf[-1L])
    mf <- do.call(stats::model.frame, mf.args, envir = parent.frame())

    ydat <- model.response(mf)
    xdat <- mf[, chromoly[[2]], drop = FALSE]
    miss.z <- !(length(chromoly) == 3)
    if (!miss.z)
      zdat <- mf[, chromoly[[3]], drop = FALSE]

    bw.args <- list(xdat = xdat, ydat = ydat)
    if (!miss.z)
      bw.args$zdat <- zdat
    tbw <- do.call(npscoefbw, c(bw.args, list(...)))

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")
    tbw$chromoly <- chromoly

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)

    tbw
  }

npscoefbw.NULL <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws, ...){
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
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

    miss.z <- missing(zdat)

    xdat <- toFrame(xdat)

    if(!miss.z)
      zdat <- toFrame(zdat)

    n.bw <- if (miss.z) ncol(xdat) else ncol(zdat)
    bws <- double(n.bw)

    bw.args <- list(xdat = xdat, ydat = ydat, bws = bws)
    if (!miss.z)
      bw.args$zdat <- zdat
    tbw <- do.call(npscoefbw.default, c(bw.args, list(...)))

    ## clean up (possible) inconsistencies due to recursion ...
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <-
      updateBwNameMetadata(nameList = list(ynames = deparse(substitute(ydat))),
                           bws = tbw)

    tbw
  }

.npscoefbw_build_scbandwidth <- function(xdat,
                                         ydat,
                                         zdat,
                                         bws,
                                         bandwidth.compute,
                                         reg.args) {
  miss.z <- is.null(zdat)
  zdati <- if (miss.z) NULL else untangle(zdat)
  znames <- if (miss.z) NULL else names(zdat)

  sbw.args <- c(
    list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      zdati = zdati,
      xnames = names(xdat),
      ynames = deparse(substitute(ydat)),
      znames = znames,
      bandwidth.compute = bandwidth.compute
    ),
    reg.args
  )

  out <- do.call(scbandwidth, sbw.args)
  if (!is.null(reg.args$scale.factor.search.lower))
    out$scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      reg.args$scale.factor.search.lower
    )
  out
}

.npscoefbw_run_fixed_degree <- function(xdat, ydat, zdat, bws, reg.args, opt.args) {
  tbw <- .npscoefbw_build_scbandwidth(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = bws,
    bandwidth.compute = opt.args$bandwidth.compute,
    reg.args = reg.args
  )

  scbw.args <- c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args)
  if (!is.null(zdat))
    scbw.args$zdat <- zdat
  do.call(npscoefbw.scbandwidth, scbw.args)
}

.npscoefbw_worker_opt_args <- function(opt.args) {
  opt.args[setdiff(names(opt.args), "zdat")]
}

.npscoefbw_run_fixed_degree_bcast_payload <- function(xdat, ydat, zdat, bws, reg.args, opt.args) {
  old.messages <- getOption("np.messages")
  rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)

  if (!isTRUE(rank == 0L))
    options(np.messages = FALSE)

  on.exit(options(np.messages = old.messages), add = TRUE)

  .npRmpi_with_local_regression(
    .npscoefbw_run_fixed_degree(
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      bws = bws,
      reg.args = reg.args,
      opt.args = opt.args
    )
  )
}

.npscoefbw_run_fixed_degree_collective <- function(xdat,
                                                   ydat,
                                                   zdat,
                                                   bws,
                                                   reg.args,
                                                   opt.args,
                                                   comm = 1L) {
  worker.opt.args <- .npscoefbw_worker_opt_args(opt.args)
  if (.npRmpi_has_active_slave_pool(comm = comm) &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    mc <- substitute(
      get(".npscoefbw_run_fixed_degree_bcast_payload", envir = asNamespace("npRmpi"), inherits = FALSE)(
        XDAT,
        YDAT,
        ZDAT,
        BWS,
        REGARGS,
        OPTARGS
      ),
      list(
        XDAT = xdat,
        YDAT = ydat,
        ZDAT = zdat,
        BWS = bws,
        REGARGS = reg.args,
        OPTARGS = worker.opt.args
      )
    )
    return(.npRmpi_bcast_cmd_expr(mc, comm = comm, caller.execute = TRUE))
  }

  .npscoefbw_run_fixed_degree(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = bws,
    reg.args = reg.args,
    opt.args = worker.opt.args
  )
}

.npscoefbw_nomad_controls <- function(search.engine) {
  .np_degree_search_engine_controls(search.engine)
}

.npscoefbw_fast_eligible <- function(sbw, eval.zdat) {
  if (!identical(sbw$type, "fixed"))
    return(FALSE)

  tdati <- if (is.null(sbw$zdati)) sbw$xdati else sbw$zdati
  eval.zdat <- toFrame(eval.zdat)

  fast_largeh_tol <- getOption("np.largeh.rel.tol", 1e-3)
  if (!is.numeric(fast_largeh_tol) || length(fast_largeh_tol) != 1L ||
      is.na(fast_largeh_tol) || !is.finite(fast_largeh_tol) ||
      fast_largeh_tol <= 0 || fast_largeh_tol >= 0.1)
    fast_largeh_tol <- 1e-3

  fast_disc_tol <- getOption("np.disc.upper.rel.tol", 1e-2)
  if (!is.numeric(fast_disc_tol) || length(fast_disc_tol) != 1L ||
      is.na(fast_disc_tol) || !is.finite(fast_disc_tol) ||
      fast_disc_tol <= 0 || fast_disc_tol >= 0.5)
    fast_disc_tol <- 1e-2

  cont_utol <- switch(
    sbw$ckertype,
    gaussian = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
    "truncated gaussian" = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
    epanechnikov = sqrt(fast_largeh_tol),
    uniform = 1.0 - 32.0 * .Machine$double.eps,
    0.0
  )

  cont_hmin <- numeric(0)
  if (any(sbw$icon) && is.finite(cont_utol) && cont_utol > 0) {
    zcon <- eval.zdat[, sbw$icon, drop = FALSE]
    cont_hmin <- vapply(zcon, function(col) {
      vals <- as.double(col)
      vals <- vals[is.finite(vals)]
      if (!length(vals))
        return(Inf)
      diff(range(vals)) / cont_utol
    }, numeric(1))
  }

  disc_upper_tol <- function(upper) {
    max(fast_disc_tol * abs(upper),
        16.0 * .Machine$double.eps * max(1.0, abs(upper)))
  }

  uno_upper <- numeric(0)
  if (any(sbw$iuno)) {
    uno_idx <- which(sbw$iuno)
    uno_upper <- vapply(uno_idx, function(i) {
      uMaxL(tdati$all.nlev[[i]], kertype = sbw$ukertype)
    }, numeric(1))
  }

  ord_upper <- numeric(0)
  if (any(sbw$iord)) {
    ord_idx <- which(sbw$iord)
    ord_upper <- vapply(ord_idx, function(i) {
      oMaxL(tdati$all.nlev[[i]], kertype = sbw$okertype)
    }, numeric(1))
  }

  bwv <- sbw$bandwidth[[1L]]
  if (!length(bwv) || length(bwv) != length(sbw$icon))
    return(FALSE)

  if (any(sbw$icon)) {
    bw_cont <- bwv[sbw$icon]
    if (any(!is.finite(bw_cont)) || any(bw_cont <= 0) ||
        any(bw_cont < cont_hmin))
      return(FALSE)
  }

  if (any(sbw$iuno)) {
    bw_uno <- bwv[sbw$iuno]
    ok_uno <- mapply(function(bw, upper) {
      is.finite(bw) && abs(bw - upper) <= disc_upper_tol(upper)
    }, bw = bw_uno, upper = uno_upper, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    if (!all(ok_uno))
      return(FALSE)
  }

  if (any(sbw$iord)) {
    bw_ord <- bwv[sbw$iord]
    ok_ord <- mapply(function(bw, upper) {
      is.finite(bw) && abs(bw - upper) <= disc_upper_tol(upper)
    }, bw = bw_ord, upper = ord_upper, SIMPLIFY = TRUE, USE.NAMES = FALSE)
    if (!all(ok_ord))
      return(FALSE)
  }

  TRUE
}

.npscoefbw_nomad_lp_context_local <- function(xdat, ydat, zdat = NULL) {
  miss.z <- is.null(zdat)
  xdat <- toFrame(xdat)

  if (!(is.vector(ydat) || is.factor(ydat)))
    stop("'ydat' must be a vector or a factor")

  if (!miss.z)
    zdat <- toFrame(zdat)

  keep.rows <- rep_len(TRUE, nrow(xdat))
  train.df <- data.frame(xdat, ydat)
  if (!miss.z)
    train.df <- data.frame(train.df, zdat)
  rows.omit <- attr(na.omit(train.df), "na.action")
  if (length(rows.omit) > 0L)
    keep.rows[as.integer(rows.omit)] <- FALSE

  if (!any(keep.rows))
    stop("Data has no rows without NAs")

  xdat <- xdat[keep.rows, , drop = FALSE]
  ydat <- ydat[keep.rows]
  if (!miss.z)
    zdat <- zdat[keep.rows, , drop = FALSE]

  yvec <- if (is.factor(ydat)) {
    dlev(ydat)[as.integer(ydat)]
  } else {
    as.double(ydat)
  }

  xmat <- toMatrix(xdat)
  if (qr(xmat)$rank < ncol(xmat))
    stop("columns of the independent variable (xdat) are linearly dependent")

  zdf <- if (miss.z) xdat else if (is.data.frame(zdat)) zdat else as.data.frame(zdat)

  list(
    xdat = xdat,
    xmat = xmat,
    ydat = yvec,
    zdat.df = zdf,
    W = cbind(1.0, xmat),
    n = nrow(xmat)
  )
}

.npscoefbw_normalize_nomad_scbw <- function(scbw, eval.zdat, bw = scbw$bw) {
  rbw <- .npscoef_make_regbw(
    bws = scbw,
    zdat = eval.zdat,
    bw = bw
  )
  scbw$bw <- rbw$bw
  scbw$bandwidth[[1L]] <- rbw$bandwidth[[1L]]
  scbw$sfactor[[1L]] <- rbw$sfactor[[1L]]
  scbw$nconfac <- rbw$nconfac
  scbw$ncatfac <- rbw$ncatfac
  scbw$sdev <- rbw$sdev
  scbw
}

.npscoefbw_nomad_lp_npksum <- function(args, localize = TRUE) {
  if (isTRUE(localize))
    return(.npRmpi_with_local_regression(do.call(npksum, args)))
  do.call(npksum, args)
}

.npscoefbw_nomad_solve_cv_moment_system <- function(tyw,
                                                    tww,
                                                    W.eval.design,
                                                    maxPenalty,
                                                    n.train,
                                                    Wz.eval = NULL) {
  neval <- ncol(tyw)
  ncoef <- nrow(tyw)
  pcoef <- ncol(W.eval.design)
  coef.out <- matrix(maxPenalty, nrow = pcoef, ncol = neval)
  ridge.grid <- npRidgeSequenceAdditive(n.train = n.train, cap = 1.0)
  ridge <- rep.int(ridge.grid[1L], neval)
  ridge.idx <- rep.int(1L, neval)
  doridge <- rep.int(TRUE, neval)

  while (any(doridge)) {
    iloo <- seq_len(neval)[doridge]
    for (ii in iloo) {
      doridge[ii] <- FALSE
      ridge.val <- ridge[ii] * tyw[, ii][1L] / NZD(tww[, , ii][1L, 1L])
      theta.ii <- tryCatch(
        solve(
          tww[, , ii] + diag(rep(ridge[ii], ncoef)),
          tyw[, ii] + c(ridge.val, rep(0, ncoef - 1L))
        ),
        error = function(e) e
      )
      if (inherits(theta.ii, "error")) {
        ridge.idx[ii] <- ridge.idx[ii] + 1L
        if (ridge.idx[ii] <= length(ridge.grid)) {
          ridge[ii] <- ridge.grid[ridge.idx[ii]]
          doridge[ii] <- TRUE
        }
        theta.ii <- rep(maxPenalty, ncoef)
      }

      if (is.null(Wz.eval)) {
        coef.out[, ii] <- theta.ii
      } else {
        coef.out[, ii] <- as.vector(crossprod(
          Wz.eval[ii, ],
          matrix(theta.ii, nrow = ncol(Wz.eval), ncol = pcoef)
        ))
      }
    }
  }

  coef.out
}

.npscoefbw_nomad_lp_eval_subset <- function(ctx,
                                            bws,
                                            idx,
                                            maxPenalty,
                                            localize = TRUE) {
  idx <- as.integer(idx)
  if (!length(idx))
    return(list(sse = 0.0, invalid = 0L))

  lp_state <- .npscoef_lp_state(
    bws = bws,
    tzdat = ctx$zdat.df,
    ezdat = ctx$zdat.df[idx, , drop = FALSE],
    leave.one.out = FALSE,
    where = "npscoefbw"
  )
  tensor.train <- .npscoef_row_tensor_design(ctx$W, lp_state$W.train)
  ytensor <- cbind(ctx$ydat, tensor.train)
  main.ks <- .npscoefbw_nomad_lp_npksum(
    args = list(
      txdat = lp_state$z.train,
      exdat = lp_state$z.eval,
      tydat = ytensor,
      weights = ytensor,
      bws = lp_state$rbw,
      leave.one.out = FALSE,
      bandwidth.divide = TRUE
    ),
    localize = localize
  )$ksum
  kw.self <- .np_kernel_weights_direct(
    bws = lp_state$rbw,
    txdat = lp_state$z.train,
    exdat = lp_state$z.eval,
    bandwidth.divide = TRUE,
    kernel.pow = 1.0
  )
  bwv <- lp_state$rbw$bandwidth[[1L]]
  bw.divisor <- if (any(lp_state$rbw$icon))
    prod(as.double(bwv[lp_state$rbw$icon]))
  else
    1.0

  for (jj in seq_along(idx)) {
    main.ks[, , jj] <- main.ks[, , jj] -
      (kw.self[idx[jj], jj] / bw.divisor) * tcrossprod(ytensor[idx[jj], ], ytensor[idx[jj], ])
  }

  tyw <- main.ks[-1L, 1L, , drop = FALSE]
  if (length(dim(tyw)) == 3L)
    dim(tyw) <- c(dim(tyw)[1L], dim(tyw)[3L])
  tww <- main.ks[-1L, -1L, , drop = FALSE]
  coef.block <- .npscoefbw_nomad_solve_cv_moment_system(
    tyw = tyw,
    tww = tww,
    W.eval.design = ctx$W[idx, , drop = FALSE],
    Wz.eval = lp_state$W.eval,
    maxPenalty = maxPenalty,
    n.train = ctx$n
  )
  mean.block <- rowSums(ctx$W[idx, , drop = FALSE] * t(coef.block))

  if (any(!is.finite(mean.block)) || any(mean.block == maxPenalty)) {
    return(list(sse = 0.0, invalid = 1L))
  }

  list(
    sse = as.numeric(sum((ctx$ydat[idx] - mean.block)^2)),
    invalid = 0L
  )
}

.npscoefbw_nomad_lp_eval_collective_local <- function(ctx,
                                                      bws,
                                                      invalid.penalty = c("baseline", "large"),
                                                      penalty.multiplier = 10,
                                                      comm = 1L) {
  invalid.penalty <- match.arg(invalid.penalty)
  base.penalty <- switch(
    invalid.penalty,
    baseline = if (is.finite(bws$fval[1L])) as.numeric(bws$fval[1L]) else 1,
    large = 1
  )
  base.penalty <- max(abs(base.penalty), 1)
  penalty <- penalty.multiplier * base.penalty

  if (!is.list(ctx) || is.null(ctx$W) || is.null(ctx$ydat) || is.null(ctx$zdat.df))
    stop("invalid NOMAD smooth-coefficient context")

  if (!validateBandwidthTF(bws)) {
    return(list(
      objective = penalty,
      num.feval = 1L,
      num.feval.fast = 0L
    ))
  }

  rank <- tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) 0L)
  size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) 1L)
  if (!is.finite(rank) || is.na(rank) || rank < 0L)
    rank <- 0L
  if (!is.finite(size) || is.na(size) || size < 1L)
    size <- 1L

  assignments <- .splitIndices(ctx$n, size)
  local.indices <- if (length(assignments) >= (rank + 1L))
    as.integer(assignments[[rank + 1L]])
  else
    integer(0L)

  local.out <- tryCatch(
    .npscoefbw_nomad_lp_eval_subset(
      ctx = ctx,
      bws = bws,
      idx = local.indices,
      maxPenalty = sqrt(.Machine$double.xmax),
      localize = TRUE
    ),
    error = function(e) list(sse = 0.0, invalid = 1L)
  )

  totals <- mpi.allreduce(
    as.double(c(local.out$sse, local.out$invalid)),
    type = 2,
    op = "sum",
    comm = comm
  )

  if (totals[2L] > 0) {
    return(list(
      objective = penalty,
      num.feval = 1L,
      num.feval.fast = 0L
    ))
  }

  list(
    objective = as.numeric(totals[1L] / ctx$n),
    num.feval = 1L,
    num.feval.fast = if (.npscoefbw_fast_eligible(bws, eval.zdat = ctx$zdat.df)) 1L else 0L
  )
}

.npscoefbw_nomad_lp_eval_direct <- function(ctx,
                                            bws,
                                            invalid.penalty = c("baseline", "large"),
                                            penalty.multiplier = 10,
                                            localize = TRUE) {
  invalid.penalty <- match.arg(invalid.penalty)
  base.penalty <- switch(
    invalid.penalty,
    baseline = if (is.finite(bws$fval[1L])) as.numeric(bws$fval[1L]) else 1,
    large = 1
  )
  base.penalty <- max(abs(base.penalty), 1)
  penalty <- penalty.multiplier * base.penalty

  if (!is.list(ctx) || is.null(ctx$W) || is.null(ctx$ydat) || is.null(ctx$zdat.df))
    stop("invalid NOMAD smooth-coefficient context")

  if (!validateBandwidthTF(bws))
    return(list(
      objective = penalty,
      num.feval = 1L,
      num.feval.fast = 0L
    ))

  maxPenalty <- sqrt(.Machine$double.xmax)

  result <- tryCatch({
    lp_state <- .npscoef_lp_state(
      bws = bws,
      tzdat = ctx$zdat.df,
      ezdat = ctx$zdat.df,
      leave.one.out = TRUE,
      where = "npscoefbw"
    )
    tensor.train <- .npscoef_row_tensor_design(ctx$W, lp_state$W.train)
    ytensor <- cbind(ctx$ydat, tensor.train)
    main.ks <- .npscoefbw_nomad_lp_npksum(
      args = list(
        txdat = lp_state$z.train,
        tydat = ytensor,
        weights = ytensor,
        bws = lp_state$rbw,
        leave.one.out = TRUE,
        bandwidth.divide = TRUE
      ),
      localize = localize
    )$ksum
    tyw <- main.ks[-1L, 1L, , drop = FALSE]
    if (length(dim(tyw)) == 3L)
      dim(tyw) <- c(dim(tyw)[1L], dim(tyw)[3L])
    tww <- main.ks[-1L, -1L, , drop = FALSE]
    coef.loo <- .npscoefbw_nomad_solve_cv_moment_system(
      tyw = tyw,
      tww = tww,
      W.eval.design = ctx$W,
      Wz.eval = lp_state$W.eval,
      maxPenalty = maxPenalty,
      n.train = ctx$n
    )
    mean.loo <- rowSums(ctx$W * t(coef.loo))

    if (any(!is.finite(mean.loo)) || any(mean.loo == maxPenalty)) {
      list(objective = penalty, num.feval = 1L, num.feval.fast = 0L)
    } else {
      list(
        objective = as.numeric(mean((ctx$ydat - mean.loo)^2)),
        num.feval = 1L,
        num.feval.fast = if (.npscoefbw_fast_eligible(bws, eval.zdat = ctx$zdat.df)) 1L else 0L
      )
    }
  }, error = function(e) {
    list(objective = penalty, num.feval = 1L, num.feval.fast = 0L)
  })

  result
}

.npscoefbw_nomad_context_prepare <- function(xdat, ydat, zdat = NULL) {
  ctx <- .npscoefbw_nomad_lp_context_local(xdat = xdat, ydat = ydat, zdat = zdat)
  if (.npRmpi_autodispatch_active() &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    mc <- match.call()
    mc[[1L]] <- get(".npscoefbw_nomad_context_prepare", envir = asNamespace("npRmpi"), inherits = FALSE)
    return(.npRmpi_autodispatch_call(mc, parent.frame()))
  }
  ctx
}

.npscoefbw_nomad_context_cleanup <- function(ctx, comm = 1L) {
  ref <- .npRmpi_autodispatch_remote_ref(ctx)
  if (!is.null(ref))
    .npRmpi_autodispatch_cleanup(ref, comm = comm)
  invisible(NULL)
}

.npscoefbw_nomad_pool_start <- function(ctx, comm = 1L) {
  if (!.npRmpi_has_active_slave_pool(comm = comm) ||
      isTRUE(.npRmpi_autodispatch_called_from_bcast()) ||
      isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    return(NULL)
  }

  ref <- .npRmpi_autodispatch_remote_ref(ctx)
  if (is.null(ref))
    stop("invalid NOMAD smooth-coefficient context: missing remote reference")

  req.base <- 61100L
  res.base <- 61200L

  mc <- substitute({
    get(".npscoefbw_nomad_slave_loop", envir = asNamespace("npRmpi"), inherits = FALSE)(
      REF = REF,
      COMM = COMM,
      REQ_BASE = REQ_BASE,
      RES_BASE = RES_BASE
    )
  }, list(
    REF = ref,
    COMM = comm,
    REQ_BASE = req.base,
    RES_BASE = res.base
  ))
  .npRmpi_bcast_cmd_expr(mc, comm = comm, caller.execute = FALSE)

  list(
    comm = comm,
    req.base = req.base,
    res.base = res.base,
    nslaves = max(0L, .npRmpi_safe_int(mpi.comm.size(comm)) - 1L)
  )
}

.npscoefbw_nomad_pool_stop <- function(pool) {
  if (is.null(pool) || !is.list(pool) || is.null(pool$nslaves) || pool$nslaves < 1L)
    return(invisible(NULL))

  for (rk in seq_len(pool$nslaves)) {
    mpi.send.Robj(
      obj = list(type = "stop"),
      dest = rk,
      tag = pool$req.base + rk,
      comm = pool$comm
    )
  }

  invisible(NULL)
}

.npscoefbw_nomad_slave_loop <- function(REF,
                                        COMM = 1L,
                                        REQ_BASE = 61100L,
                                        RES_BASE = 61200L) {
  rank <- tryCatch(as.integer(mpi.comm.rank(COMM)), error = function(e) 0L)
  if (isTRUE(rank == 0L))
    return(invisible(NULL))

  old.messages <- getOption("np.messages")
  options(np.messages = FALSE)
  on.exit(options(np.messages = old.messages), add = TRUE)

  ctx <- get(REF, envir = .GlobalEnv, inherits = FALSE)

  repeat {
    msg <- mpi.recv.Robj(source = 0L, tag = REQ_BASE + rank, comm = COMM)
    if (is.list(msg) && identical(msg$type, "stop"))
      break

    out <- tryCatch(
      .npscoefbw_nomad_lp_eval_subset(
        ctx = ctx,
        bws = msg$bws,
        idx = as.integer(msg$idx),
        maxPenalty = sqrt(.Machine$double.xmax),
        localize = TRUE
      ),
      error = function(e) list(sse = 0.0, invalid = 1L)
    )

    mpi.send.Robj(
      obj = out,
      dest = 0L,
      tag = RES_BASE + rank,
      comm = COMM
    )
  }

  invisible(NULL)
}

.npscoefbw_eval_pool <- function(ctx,
                                 bws,
                                 pool,
                                 invalid.penalty = c("baseline", "large"),
                                 penalty.multiplier = 10) {
  invalid.penalty <- match.arg(invalid.penalty)
  base.penalty <- switch(
    invalid.penalty,
    baseline = if (is.finite(bws$fval[1L])) as.numeric(bws$fval[1L]) else 1,
    large = 1
  )
  base.penalty <- max(abs(base.penalty), 1)
  penalty <- penalty.multiplier * base.penalty

  if (!validateBandwidthTF(bws)) {
    return(list(
      objective = penalty,
      num.feval = 1L,
      num.feval.fast = 0L
    ))
  }

  if (is.null(pool) || !is.list(pool) || is.null(pool$nslaves) || pool$nslaves < 1L) {
    return(.npscoefbw_nomad_lp_eval_direct(
      ctx = ctx,
      bws = bws,
      invalid.penalty = invalid.penalty,
      penalty.multiplier = penalty.multiplier,
      localize = TRUE
    ))
  }

  size <- pool$nslaves + 1L
  assignments <- .splitIndices(ctx$n, size)
  local.indices <- if (length(assignments)) as.integer(assignments[[1L]]) else integer(0L)

  for (rk in seq_len(pool$nslaves)) {
    mpi.send.Robj(
      obj = list(
        type = "eval",
        bws = bws,
        idx = as.integer(assignments[[rk + 1L]])
      ),
      dest = rk,
      tag = pool$req.base + rk,
      comm = pool$comm
    )
  }

  local.out <- tryCatch(
    .npscoefbw_nomad_lp_eval_subset(
      ctx = ctx,
      bws = bws,
      idx = local.indices,
      maxPenalty = sqrt(.Machine$double.xmax),
      localize = TRUE
    ),
    error = function(e) list(sse = 0.0, invalid = 1L)
  )

  worker.out <- lapply(seq_len(pool$nslaves), function(rk) {
    mpi.recv.Robj(
      source = rk,
      tag = pool$res.base + rk,
      comm = pool$comm
    )
  })

  invalid.total <- local.out$invalid + sum(vapply(worker.out, `[[`, numeric(1), "invalid"))
  if (invalid.total > 0) {
    return(list(
      objective = penalty,
      num.feval = 1L,
      num.feval.fast = 0L
    ))
  }

  total.sse <- local.out$sse + sum(vapply(worker.out, `[[`, numeric(1), "sse"))

  list(
    objective = as.numeric(total.sse / ctx$n),
    num.feval = 1L,
    num.feval.fast = if (.npscoefbw_fast_eligible(bws, eval.zdat = ctx$zdat.df)) 1L else 0L
  )
}

.npscoefbw_eval_only <- function(xdat,
                                 ydat,
                                 zdat,
                                 bws,
                                 invalid.penalty = c("baseline", "large"),
                                 penalty.multiplier = 10) {
  invalid.penalty <- match.arg(invalid.penalty)
  base.penalty <- switch(
    invalid.penalty,
    baseline = if (is.finite(bws$fval[1L])) as.numeric(bws$fval[1L]) else 1,
    large = 1
  )
  base.penalty <- max(abs(base.penalty), 1)
  penalty <- penalty.multiplier * base.penalty

  fit <- tryCatch(
    .npRmpi_with_local_regression(
      npscoef(
        bws = bws,
        txdat = xdat,
        tydat = ydat,
        tzdat = zdat,
        leave.one.out = TRUE,
        iterate = FALSE,
        betas = FALSE,
        errors = FALSE,
        .np_fit_progress_allow = FALSE
      )
    ),
    error = function(e) e
  )

  if (inherits(fit, "error") || is.null(fit$mean) || any(!is.finite(fit$mean))) {
    return(list(objective = penalty, num.feval = 1L, num.feval.fast = 0L))
  }

  list(
    objective = as.numeric(mean((as.double(ydat) - as.double(fit$mean))^2)),
    num.feval = 1L,
    num.feval.fast = if (.npscoefbw_fast_eligible(bws, eval.zdat = if (is.null(zdat)) xdat else zdat)) 1L else 0L
  )
}

.npscoefbw_eval_collective <- function(data,
                                       bws,
                                       invalid.penalty = c("baseline", "large"),
                                       penalty.multiplier = 10,
                                       comm = 1L) {
  if (!is.list(data))
    stop("invalid NOMAD smooth-coefficient context")

  if (.npRmpi_has_active_slave_pool(comm = comm) &&
      !isTRUE(.npRmpi_autodispatch_called_from_bcast()) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE))) {
    ref <- .npRmpi_autodispatch_remote_ref(data)
    if (is.null(ref))
      stop("invalid NOMAD smooth-coefficient context: missing remote reference")
    mc <- substitute({
      old.messages <- getOption("np.messages")
      rank <- tryCatch(as.integer(mpi.comm.rank(COMM)), error = function(e) 0L)
      if (!isTRUE(rank == 0L))
        options(np.messages = FALSE)
      on.exit(options(np.messages = old.messages), add = TRUE)
      ctx.local <- get(REF, envir = .GlobalEnv, inherits = FALSE)
      get(".npscoefbw_nomad_lp_eval_collective_local", envir = asNamespace("npRmpi"), inherits = FALSE)(
        ctx = ctx.local,
        bws = BWS,
        invalid.penalty = INVALID,
        penalty.multiplier = PMULT,
        comm = COMM
      )
    }, list(
      REF = ref,
      BWS = bws,
      INVALID = invalid.penalty,
      PMULT = penalty.multiplier,
      COMM = comm
    ))
    return(.npRmpi_bcast_cmd_expr(mc, comm = comm, caller.execute = TRUE))
  }

  if (.npRmpi_safe_int(mpi.comm.size(comm)) > 1L) {
    return(.npscoefbw_nomad_lp_eval_collective_local(
      ctx = data,
      bws = bws,
      invalid.penalty = invalid.penalty,
      penalty.multiplier = penalty.multiplier,
      comm = comm
    ))
  }

  .npscoefbw_nomad_lp_eval_direct(
    ctx = data,
    bws = bws,
    invalid.penalty = invalid.penalty,
    penalty.multiplier = penalty.multiplier,
    localize = TRUE
  )
}

.npscoefbw_nomad_search <- function(xdat,
                                    ydat,
                                    zdat,
                                    bws,
                                    reg.args,
                                    opt.args,
                                    degree.search,
                                    nomad.inner.nmulti = 0L) {
  if (isTRUE(degree.search$verify))
    stop("automatic degree search with search.engine='nomad' does not support degree.verify")

  if (!identical(opt.args$bandwidth.compute, TRUE))
    stop("automatic degree search with search.engine='nomad' requires bandwidth.compute=TRUE")

  eval.zdat <- if (is.null(zdat)) xdat else zdat

  template.reg.args <- reg.args
  template.reg.args$regtype <- "lp"
  template.reg.args$degree <- as.integer(degree.search$start.degree)
  template.reg.args$bernstein.basis <- degree.search$bernstein.basis

  template <- .npscoefbw_build_scbandwidth(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = bws,
    bandwidth.compute = FALSE,
    reg.args = template.reg.args
  )

  if (!identical(template$type, "fixed"))
    stop("automatic degree search with search.engine='nomad' currently requires bwtype='fixed'")
  setup <- .npregbw_nomad_bw_setup(xdat = eval.zdat, template = template)
  ncon <- length(setup$cont_idx)
  ncat <- length(setup$cat_idx)
  ndeg <- length(degree.search$start.degree)
  nomad.nmulti <- if (is.null(opt.args$nmulti)) npDefaultNmulti(NCOL(eval.zdat)) else npValidateNmulti(opt.args$nmulti[1L])

  cont_lower <- npGetScaleFactorSearchLower(
    template,
    argname = "template$scale.factor.search.lower"
  )
  bw_lower <- c(rep.int(cont_lower, ncon), rep.int(0, ncat))
  bw_upper <- c(rep.int(1e6, ncon), setup$cat_upper * setup$bandwidth.scale.categorical)

  x0 <- c(
    .np_nomad_complete_start_point(
      point = if (all(template$bw == 0)) NULL else .npregbw_nomad_bw_to_point(template$bw, template = template, setup = setup),
      lower = bw_lower,
      upper = bw_upper,
      ncont = ncon
    ),
    as.integer(degree.search$start.degree)
  )
  lb <- c(bw_lower, degree.search$lower)
  ub <- c(bw_upper, degree.search$upper)
  bbin <- c(rep.int(0L, ncon + ncat), rep.int(1L, ndeg))
  baseline.record <- NULL
  nomad.num.feval.total <- 0
  nomad.num.feval.fast.total <- 0
  ctx <- .npscoefbw_nomad_context_prepare(xdat = xdat, ydat = ydat, zdat = zdat)
  on.exit(.npscoefbw_nomad_context_cleanup(ctx, comm = 1L), add = TRUE)
  pool <- .npscoefbw_nomad_pool_start(ctx, comm = 1L)
  on.exit(.npscoefbw_nomad_pool_stop(pool), add = TRUE)

  .np_nomad_baseline_note(degree.search$start.degree)

  eval_fun <- function(point) {
    point <- as.numeric(point)
    degree <- as.integer(round(point[ncon + ncat + seq_len(ndeg)]))
    degree <- .np_degree_clip_to_grid(degree, degree.search$candidates)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)

    eval.reg.args <- reg.args
    eval.reg.args$regtype <- "lp"
    eval.reg.args$degree <- degree
    eval.reg.args$bernstein.basis <- degree.search$bernstein.basis

    tbw <- .npscoefbw_build_scbandwidth(
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      bws = bw_vec,
      bandwidth.compute = FALSE,
      reg.args = eval.reg.args
    )
    tbw <- .npscoefbw_normalize_nomad_scbw(
      scbw = tbw,
      eval.zdat = eval.zdat,
      bw = bw_vec
    )

    out <- .npscoefbw_eval_pool(
      ctx = ctx,
      bws = tbw,
      pool = pool,
      invalid.penalty = "baseline",
      penalty.multiplier = if (is.null(opt.args$penalty.multiplier)) 10 else opt.args$penalty.multiplier
    )
    nomad.num.feval.total <<- nomad.num.feval.total + as.numeric(out$num.feval[1L])
    nomad.num.feval.fast.total <<- nomad.num.feval.fast.total + as.numeric(out$num.feval.fast[1L])

    list(
      objective = out$objective,
      degree = degree,
      num.feval = out$num.feval,
      num.feval.fast = out$num.feval.fast
    )
  }

  build_payload <- function(point, best_record, solution, interrupted) {
    point <- as.numeric(point)
    degree <- as.integer(best_record$degree)
    bw_vec <- .npregbw_nomad_point_to_bw(point[seq_len(ncon + ncat)], template = template, setup = setup)
    powell.elapsed <- NA_real_

    build_direct_payload <- function() {
      final.reg.args <- reg.args
      final.reg.args$regtype <- "lp"
      final.reg.args$degree <- degree
      final.reg.args$bernstein.basis <- degree.search$bernstein.basis

      tbw <- .npscoefbw_build_scbandwidth(
        xdat = xdat,
        ydat = ydat,
        zdat = zdat,
        bws = bw_vec,
        bandwidth.compute = FALSE,
        reg.args = final.reg.args
      )
      tbw <- .npscoefbw_normalize_nomad_scbw(
        scbw = tbw,
        eval.zdat = eval.zdat,
        bw = bw_vec
      )
      tbw$fval <- as.numeric(best_record$objective)
      tbw$ifval <- as.numeric(best_record$objective)
      tbw$num.feval <- as.numeric(nomad.num.feval.total)
      tbw$num.feval.fast <- as.numeric(nomad.num.feval.fast.total)
      tbw$numimp <- 0
      tbw$fval.vector <- as.numeric(best_record$objective)
      tbw$total.time <- NA_real_
      if (!is.null(tbw$method) && length(tbw$method))
        tbw$pmethod <- bwmToPrint(as.character(tbw$method[1L]))
      tbw
    }

    direct.payload <- build_direct_payload()
    if (is.null(direct.payload$timing.profile) && is.list(best_record$timing.profile))
      direct.payload$timing.profile <- best_record$timing.profile
    direct.objective <- as.numeric(best_record$objective)

    if (identical(degree.search$engine, "nomad+powell")) {
      hot.reg.args <- reg.args
      hot.reg.args$regtype <- "lp"
      hot.reg.args$degree <- degree
      hot.reg.args$bernstein.basis <- degree.search$bernstein.basis
      hot.opt.args <- opt.args
      hot.opt.args$nmulti <- .np_nomad_powell_hotstart_nmulti("single_iteration")
      powell.start <- proc.time()[3L]
      hot.payload <- .np_nomad_with_powell_progress(
        degree,
        .npscoefbw_run_fixed_degree_collective(
          xdat = xdat,
          ydat = ydat,
          zdat = zdat,
          bws = bw_vec,
          reg.args = hot.reg.args,
          opt.args = hot.opt.args
        )
      )
      powell.elapsed <- proc.time()[3L] - powell.start
      direct.payload$num.feval <- as.numeric(direct.payload$num.feval[1L]) + as.numeric(hot.payload$num.feval[1L])
      direct.payload$num.feval.fast <- as.numeric(direct.payload$num.feval.fast[1L]) + as.numeric(hot.payload$num.feval.fast[1L])
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

.npscoefbw_degree_search_controls <- function(regtype,
                                              regtype.named,
                                              cv.iterate,
                                              cv.iterate.named,
                                              bandwidth.compute,
                                              ncon,
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
  search.engine <- .npscoefbw_nomad_controls(search.engine)

  regtype.requested <- if (isTRUE(regtype.named)) match.arg(regtype, c("lc", "ll", "lp")) else "lc"
  if (!identical(regtype.requested, "lp"))
    stop("automatic degree search currently requires regtype='lp'")
  if (!isTRUE(bandwidth.compute))
    stop("automatic degree search requires bandwidth.compute=TRUE")
  if (isTRUE(cv.iterate.named) && isTRUE(cv.iterate))
    stop("automatic degree search currently requires cv.iterate=FALSE")
  if (ncon < 1L)
    stop("automatic degree search requires at least one continuous smoothing predictor")

  bern.auto <- if (isTRUE(bernstein.named)) bernstein.basis else TRUE
  bern.auto <- npValidateGlpBernstein(regtype = "lp", bernstein.basis = bern.auto)

  bounds <- .np_degree_normalize_bounds(
    ncon = ncon,
    degree.min = degree.min,
    degree.max = degree.max,
    default.max = 3L
  )

  baseline.degree <- rep.int(0L, ncon)
  default.start.degree <- if (identical(search.engine, "cell")) {
    baseline.degree
  } else {
    rep.int(1L, ncon)
  }
  start.degree <- if (is.null(degree.start)) {
    pmax(bounds$lower, pmin(bounds$upper, default.start.degree))
  } else {
    start.raw <- npValidateGlpDegree(regtype = "lp", degree = degree.start, ncon = ncon, argname = "degree.start")
    out.of.range <- vapply(seq_len(ncon), function(j) !(start.raw[j] %in% bounds$candidates[[j]]), logical(1))
    if (any(out.of.range))
      stop("degree.start must lie within the searched degree candidates for every continuous smoothing predictor")
    start.raw
  }

  list(
    method = if (identical(search.engine, "cell")) degree.select else search.engine,
    engine = search.engine,
    candidates = bounds$candidates,
    lower = bounds$lower,
    upper = bounds$upper,
    baseline.degree = baseline.degree,
    start.degree = start.degree,
    start.user = !is.null(degree.start),
    basis = if (missing(basis) || is.null(basis)) "glp" else as.character(basis[1L]),
    nobs = as.integer(nobs[1L]),
    restarts = npValidateNonNegativeInteger(degree.restarts, "degree.restarts"),
    max.cycles = npValidatePositiveInteger(degree.max.cycles, "degree.max.cycles"),
    verify = npValidateScalarLogical(degree.verify, "degree.verify"),
    bernstein.basis = bern.auto
  )
}

.npscoefbw_attach_degree_search <- function(bws, search_result) {
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

.npscoef_nn_candidate_bandwidth <- function(param, bwtype, nobs) {
  if (identical(bwtype, "fixed"))
    return(as.double(param))

  lower <- 2L
  upper <- max(1L, as.integer(nobs) - 1L)
  vapply(param, function(h) {
    if (!is.finite(h))
      return(NA_real_)
    as.double(max(lower, min(upper, .np_round_half_to_even(h))))
  }, numeric(1))
}

.npscoefbw_start_controls <- function(scale.factor.init.lower = 0.1,
                                     scale.factor.init.upper = 2.0,
                                     scale.factor.init = 0.5,
                                     lbd.init = 0.5,
                                     hbd.init = 1.5,
                                     dfac.init = 1.0,
                                     scale.factor.search.lower = 0,
                                     where = "npscoefbw") {
  cont.start <- npContinuousSearchStartControls(
    scale.factor.init.lower,
    scale.factor.init.upper,
    scale.factor.init,
    scale.factor.search.lower,
    where = where
  )
  lbd.init <- npValidatePositiveFiniteNumeric(lbd.init, "lbd.init")
  hbd.init <- npValidatePositiveFiniteNumeric(hbd.init, "hbd.init")
  dfac.init <- npValidatePositiveFiniteNumeric(dfac.init, "dfac.init")

  if (hbd.init < lbd.init) {
    stop(sprintf("%s: 'hbd.init' must be greater than or equal to 'lbd.init'", where),
         call. = FALSE)
  }
  if (lbd.init > 2 || hbd.init > 2 || dfac.init > 2) {
    stop(sprintf("%s: categorical start factors must be less than or equal to 2", where),
         call. = FALSE)
  }

  list(
    scale.factor.init.lower = cont.start$scale.factor.init.lower,
    scale.factor.init.upper = cont.start$scale.factor.init.upper,
    scale.factor.init = cont.start$scale.factor.init,
    scale.factor.search.lower = as.double(scale.factor.search.lower),
    lbd.init = as.double(lbd.init),
    hbd.init = as.double(hbd.init),
    dfac.init = as.double(dfac.init)
  )
}

.npscoef_start_factor_vector <- function(param,
                                         icon = NULL,
                                         iord = NULL,
                                         iuno = NULL,
                                         continuous.factor,
                                         categorical.factor,
                                         where = "npscoefbw") {
  ndim <- length(param)
  if (is.null(icon) || is.null(iord) || is.null(iuno))
    return(rep.int(as.double(continuous.factor), ndim))

  icon <- as.logical(icon)
  iord <- as.logical(iord)
  iuno <- as.logical(iuno)

  if (length(icon) != ndim || length(iord) != ndim || length(iuno) != ndim) {
    stop(sprintf("%s: invalid fixed-start coordinate map", where), call. = FALSE)
  }

  cat.mask <- iord | iuno
  factor <- rep_len(NA_real_, ndim)
  factor[icon] <- as.double(continuous.factor)
  factor[cat.mask] <- as.double(categorical.factor)

  if (anyNA(factor)) {
    stop(sprintf("%s: unsupported fixed-start coordinate type", where), call. = FALSE)
  }

  factor
}

.npscoef_default_start_bandwidth <- function(param,
                                             bwtype,
                                             nobs,
                                             start.controls = .npscoefbw_start_controls(),
                                             icon = NULL,
                                             iord = NULL,
                                             iuno = NULL) {
  if (identical(bwtype, "fixed"))
    return(as.double(
      param * .npscoef_start_factor_vector(
        param = param,
        icon = icon,
        iord = iord,
        iuno = iuno,
        continuous.factor = start.controls$scale.factor.init,
        categorical.factor = start.controls$dfac.init
      )
    ))

  lower <- 2L
  upper <- max(1L, as.integer(nobs) - 1L)
  start <- max(lower, min(upper, .np_round_half_to_even(sqrt(nobs))))
  rep.int(as.double(start), length(param))
}

.npscoef_random_start_bandwidth <- function(param,
                                            bwtype,
                                            nobs,
                                            start.controls = .npscoefbw_start_controls(),
                                            icon = NULL,
                                            iord = NULL,
                                            iuno = NULL) {
  if (identical(bwtype, "fixed")) {
    draws <- .npscoef_start_factor_vector(
      param = param,
      icon = icon,
      iord = iord,
      iuno = iuno,
      continuous.factor = start.controls$scale.factor.init.lower,
      categorical.factor = start.controls$lbd.init
    )

    for (ii in seq_along(draws)) {
      if (!is.null(icon) && !is.null(iord) && !is.null(iuno) &&
          !isTRUE(as.logical(icon)[ii])) {
        draws[ii] <- runif(1L, min = start.controls$lbd.init, max = start.controls$hbd.init)
      } else {
        draws[ii] <- runif(1L, min = start.controls$scale.factor.init.lower, max = start.controls$scale.factor.init.upper)
      }
    }

    return(as.double(draws * param))
  }

  upper <- max(1L, as.integer(nobs) - 1L)
  .npscoef_nn_candidate_bandwidth(
    param = runif(length(param), min = 2, max = max(2L, upper)),
    bwtype = bwtype,
    nobs = nobs
  )
}

.npscoef_candidate_is_admissible <- function(param,
                                            bwtype,
                                            nobs,
                                            lower = NULL) {
  candidate <- .npscoef_nn_candidate_bandwidth(
    param = param,
    bwtype = bwtype,
    nobs = nobs
  )
  if (any(!is.finite(candidate)))
    return(FALSE)
  if (identical(bwtype, "fixed")) {
    if (!is.null(lower))
      return(all(candidate >= lower))
    return(all(candidate > 0))
  }
  TRUE
}

.npscoef_finalize_bandwidth <- function(param,
                                        bwtype,
                                        nobs,
                                        lower = NULL,
                                        where = "npscoefbw") {
  candidate <- .npscoef_nn_candidate_bandwidth(param = param, bwtype = bwtype, nobs = nobs)
  if (any(!is.finite(candidate))) {
    if (identical(bwtype, "fixed")) {
      stop(sprintf("%s: bandwidth must be finite", where), call. = FALSE)
    }
    stop(
      sprintf(
        "%s: nearest-neighbor bandwidth must be an integer vector in [2, %d]",
        where,
        max(2L, as.integer(nobs) - 1L)
      ),
      call. = FALSE
    )
  }
  if (identical(bwtype, "fixed") && any(candidate <= 0)) {
    stop(sprintf("%s: bandwidth must be strictly positive", where), call. = FALSE)
  }
  if (identical(bwtype, "fixed") && !is.null(lower) && any(candidate < lower)) {
    stop(sprintf("%s: bandwidth is below the continuous scale-factor lower bound", where),
         call. = FALSE)
  }
  as.double(candidate)
}

npscoefbw.scbandwidth <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws,
           backfit.iterate = FALSE,
           backfit.maxiter = 100,
           backfit.tol = .Machine$double.eps,
           bandwidth.compute = TRUE,
           cv.iterate = FALSE,
           cv.num.iterations = 1,
           nmulti,
           optim.abstol = .Machine$double.eps,
           optim.maxattempts = 10,
           optim.maxit = 500,
           optim.method = c("Nelder-Mead", "BFGS", "CG"),
           optim.reltol = sqrt(.Machine$double.eps),
           random.seed = 42,
           scale.factor.init.lower = 0.1,
           scale.factor.init.upper = 2.0,
           scale.factor.init = 0.5,
           lbd.init = 0.5,
           hbd.init = 1.5,
           dfac.init = 1.0,
           scale.factor.search.lower = NULL,
           ...){
    ## Save seed prior to setting

    seed.state <- .np_seed_enter(random.seed)


    miss.z <- missing(zdat)

    xdat <- toFrame(xdat)

    if (!miss.z)
      zdat <- toFrame(zdat)

    if (missing(nmulti)){
      nmulti <- npDefaultNmulti(if (miss.z) NCOL(xdat) else NCOL(zdat))
    }
    regtype <- if (is.null(bws$regtype)) "lc" else bws$regtype
    cv.iterate <- npValidateScalarLogical(cv.iterate, "cv.iterate")
    backfit.iterate <- npValidateScalarLogical(backfit.iterate, "backfit.iterate")
    bandwidth.compute <- npValidateScalarLogical(bandwidth.compute, "bandwidth.compute")
    nmulti <- npValidateNmulti(nmulti)
    .np_progress_bandwidth_set_total(nmulti)
    backfit.maxiter <- npValidatePositiveInteger(backfit.maxiter, "backfit.maxiter")
    backfit.tol <- npValidatePositiveFiniteNumeric(backfit.tol, "backfit.tol")
    optim.maxattempts <- npValidatePositiveInteger(optim.maxattempts, "optim.maxattempts")
    optim.maxit <- npValidatePositiveInteger(optim.maxit, "optim.maxit")
    optim.reltol <- npValidatePositiveFiniteNumeric(optim.reltol, "optim.reltol")
    optim.abstol <- npValidatePositiveFiniteNumeric(optim.abstol, "optim.abstol")
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(
      if (is.null(scale.factor.search.lower)) npGetScaleFactorSearchLower(bws) else scale.factor.search.lower
    )
    start.controls <- .npscoefbw_start_controls(
      scale.factor.init.lower = scale.factor.init.lower,
      scale.factor.init.upper = scale.factor.init.upper,
      scale.factor.init = scale.factor.init,
      lbd.init = lbd.init,
      hbd.init = hbd.init,
      dfac.init = dfac.init,
      scale.factor.search.lower = scale.factor.search.lower,
      where = "npscoefbw"
    )
    if (cv.iterate)
      cv.num.iterations <- npValidatePositiveInteger(cv.num.iterations, "cv.num.iterations")
    spec <- .npscoef_canonical_spec(
      source = bws,
      zdat = if (miss.z) xdat else zdat,
      where = "npscoefbw"
    )
    reg.engine <- spec$regtype.engine
    if (!identical(reg.engine, "lc") && cv.iterate)
      stop("cv.iterate currently supports regtype='lc' for npscoefbw")
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (.npRmpi_autodispatch_active())
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector or a factor")

    if (miss.z) {
      bwMatch(xdat, bws$xdati)
    } else {
      bwMatch(zdat, bws$zdati)
    }

    if (dim(xdat)[1] != length(ydat))
      stop("number of regression data and response data do not match")

    if (ncol(xdat) == 1 && missing(cv.iterate))
      cv.iterate = FALSE

    if (!all(bws$xdati$icon))
      stop("Only continuous 'x' regressors are supported in this version.")

    optim.method <- match.arg(optim.method)

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(xdat))
    train.df <- data.frame(xdat, ydat)
    if (!miss.z)
      train.df <- data.frame(train.df, zdat)
    rows.omit <- attr(na.omit(train.df), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Data has no rows without NAs")

    xdat <- xdat[keep.rows,,drop = FALSE]
    ydat <- ydat[keep.rows]

    if(!miss.z)
      zdat <- zdat[keep.rows,, drop = FALSE]

    nrow = dim(xdat)[1]
    ncol = dim(xdat)[2]

    ## at this stage, data to be sent to the c routines must be converted to
    ## numeric type.

    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)

    xdat <- toMatrix(xdat)

    ## if (!miss.z)
    ##  zdat <- toMatrix(zdat)

    ## bad data
    if (qr(xdat)$rank < ncol(xdat)){
      stop("columns of the independent variable (xdat) are linearly dependent")
    }

    n <- nrow(xdat)

    ## ... do bandwidth selection

    ## construct 'W' matrix
    ## in the future one will be able to use a switch to npksum
    ## to emulate W

    W <- cbind(1.0, xdat)
    yW <- cbind(ydat, W)

    if (miss.z){
      zdat <- xdat
      dati <- bws$xdati
    }
    else
      dati <- bws$zdati
    zdat.df <- if (is.data.frame(zdat)) zdat else as.data.frame(zdat)

    mysd <- EssDee(zdat[, dati$icon, drop = FALSE])
    nconfac <- n^(-1.0/(2.0*bws$ckerorder+bws$ncon))
    ncatfac <- n^(-2.0/(2.0*bws$ckerorder+bws$ncon))

    bws$sdev <- mysd
    bws$nconfac <- nconfac
    bws$ncatfac <- ncatfac
    bw.scale.multiplier <- NULL
    if (bws$scaling) {
      bw.scale.multiplier <- rep(ncatfac, bws$ndim)
      if (any(bws$icon)) {
        icon.cumsum <- cumsum(dati$icon)
        bw.scale.multiplier[bws$icon] <- nconfac * bws$sdev[icon.cumsum[bws$icon]]
      }
    }
    apply_bw_to_scbw <- function(scbw, param) {
      param <- .npscoef_nn_candidate_bandwidth(param = param, bwtype = scbw$type, nobs = n)
      scbw$bw <- param
      if (scbw$scaling)
        scbw$bandwidth[[1]] <- scbw$bw * bw.scale.multiplier
      else
        scbw$bandwidth[[1]] <- scbw$bw
      scbw
    }

    fast_largeh_tol <- getOption("np.largeh.rel.tol", 1e-3)
    if (!is.numeric(fast_largeh_tol) || length(fast_largeh_tol) != 1L ||
        is.na(fast_largeh_tol) || !is.finite(fast_largeh_tol) ||
        fast_largeh_tol <= 0 || fast_largeh_tol >= 0.1)
      fast_largeh_tol <- 1e-3

    fast_disc_tol <- getOption("np.disc.upper.rel.tol", 1e-2)
    if (!is.numeric(fast_disc_tol) || length(fast_disc_tol) != 1L ||
        is.na(fast_disc_tol) || !is.finite(fast_disc_tol) ||
        fast_disc_tol <= 0 || fast_disc_tol >= 0.5)
      fast_disc_tol <- 1e-2

    cont_utol <- switch(
      bws$ckertype,
      gaussian = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
      "truncated gaussian" = sqrt(-2.0 * log(1.0 - fast_largeh_tol)),
      epanechnikov = sqrt(fast_largeh_tol),
      uniform = 1.0 - 32.0 * .Machine$double.eps,
      0.0
    )

    cont_hmin <- numeric(0)
    if (any(dati$icon) && is.finite(cont_utol) && cont_utol > 0) {
      zcon <- zdat.df[, dati$icon, drop = FALSE]
      cont_hmin <- vapply(zcon, function(col) {
        vals <- as.double(col)
        vals <- vals[is.finite(vals)]
        if (!length(vals))
          return(Inf)
        diff(range(vals)) / cont_utol
      }, numeric(1))
    }

    disc_upper_tol <- function(upper) {
      max(fast_disc_tol * abs(upper),
          16.0 * .Machine$double.eps * max(1.0, abs(upper)))
    }

    uno_upper <- numeric(0)
    if (any(dati$iuno)) {
      uno_idx <- which(dati$iuno)
      uno_upper <- vapply(uno_idx, function(i) {
        uMaxL(dati$all.nlev[[i]], kertype = bws$ukertype)
      }, numeric(1))
    }

    ord_upper <- numeric(0)
    if (any(dati$iord)) {
      ord_idx <- which(dati$iord)
      ord_upper <- vapply(ord_idx, function(i) {
        oMaxL(dati$all.nlev[[i]], kertype = bws$okertype)
      }, numeric(1))
    }

    npscoef_fast_eligible <- function(sbw) {
      .npscoefbw_fast_eligible(sbw = sbw, eval.zdat = zdat.df)
    }

    solve_cv_moment_system <- function(tyw, tww, W.eval.design, maxPenalty, Wz.eval = NULL) {
      neval <- ncol(tyw)
      ncoef <- nrow(tyw)
      pcoef <- ncol(W.eval.design)
      coef.out <- matrix(maxPenalty, nrow = pcoef, ncol = neval)
      ridge.grid <- npRidgeSequenceAdditive(n.train = n, cap = 1.0)
      ridge <- rep.int(ridge.grid[1L], neval)
      ridge.idx <- rep.int(1L, neval)
      doridge <- rep.int(TRUE, neval)

      while(any(doridge)){
        iloo <- seq_len(neval)[doridge]
        for (ii in iloo) {
          doridge[ii] <- FALSE
          ridge.val <- ridge[ii]*tyw[,ii][1]/NZD(tww[,,ii][1,1])
          theta.ii <- tryCatch(
            solve(tww[,,ii] + diag(rep(ridge[ii], ncoef)),
                  tyw[,ii] + c(ridge.val, rep(0, ncoef - 1))),
            error = function(e) e
          )
          if (inherits(theta.ii, "error")) {
            ridge.idx[ii] <- ridge.idx[ii] + 1L
            if (ridge.idx[ii] <= length(ridge.grid)) {
              ridge[ii] <- ridge.grid[ridge.idx[ii]]
              doridge[ii] <- TRUE
            }
            theta.ii <- rep(maxPenalty, ncoef)
          }

          if (is.null(Wz.eval)) {
            coef.out[,ii] <- theta.ii
          } else {
            coef.out[,ii] <- as.vector(crossprod(
              Wz.eval[ii,],
              matrix(theta.ii, nrow = ncol(Wz.eval), ncol = pcoef)
            ))
          }
        }
      }

      coef.out
    }

    lp_full_coef <- function(sbw, leave.one.out.eval) {
      lp_state <- .npscoef_lp_state(
        bws = sbw,
        tzdat = zdat.df,
        ezdat = zdat.df,
        leave.one.out = leave.one.out.eval,
        where = "npscoefbw"
      )
      tensor.train <- .npscoef_row_tensor_design(W, lp_state$W.train)
      ytensor <- cbind(ydat, tensor.train)
      ksum.args <- list(
        txdat = lp_state$z.train,
        tydat = ytensor,
        weights = ytensor,
        bws = lp_state$rbw,
        leave.one.out = leave.one.out.eval,
        bandwidth.divide = TRUE
      )
      main.ks <- do.call(npksum, ksum.args)$ksum
      tyw <- main.ks[-1L, 1L, , drop = FALSE]
      if (length(dim(tyw)) == 3L)
        dim(tyw) <- c(dim(tyw)[1L], dim(tyw)[3L])
      tww <- main.ks[-1L, -1L, , drop = FALSE]
      solve_cv_moment_system(
        tyw = tyw,
        tww = tww,
        W.eval.design = W,
        Wz.eval = lp_state$W.eval,
        maxPenalty = maxPenalty
      )
    }

    lp_partial_coef <- function(sbw, wj, partial.y, leave.one.out.eval) {
      lp_state <- .npscoef_lp_state(
        bws = sbw,
        tzdat = zdat.df,
        ezdat = zdat.df,
        leave.one.out = leave.one.out.eval,
        where = "npscoefbw"
      )
      U <- lp_state$W.train * wj
      yU <- cbind(partial.y, U)
      ksum.args <- list(
        txdat = lp_state$z.train,
        tydat = yU,
        weights = yU,
        bws = lp_state$rbw,
        leave.one.out = leave.one.out.eval,
        bandwidth.divide = TRUE
      )
      main.ks <- do.call(npksum, ksum.args)$ksum
      tyw <- main.ks[-1L, 1L, , drop = FALSE]
      if (length(dim(tyw)) == 3L)
        dim(tyw) <- c(dim(tyw)[1L], dim(tyw)[3L])
      tww <- main.ks[-1L, -1L, , drop = FALSE]
      as.vector(solve_cv_moment_system(
        tyw = tyw,
        tww = tww,
        W.eval.design = matrix(1.0, nrow = n, ncol = 1L),
        Wz.eval = lp_state$W.eval,
        maxPenalty = maxPenalty
      ))
    }

    total.time <-
      system.time({
        if (bandwidth.compute){
          maxPenalty <- sqrt(.Machine$double.xmax)
          cv_state <- new.env(parent = emptyenv())
          cv_state$fast_total <- 0L
          cv_state$objective_fast <- FALSE
          cv_state$optim_progress <- NULL
          cv_state$optim_eval <- 0L
          cv_state$multistart_index <- NA_integer_
          cv_state$partial_progress <- NULL
          cv_state$partial_eval <- 0L
          cv_state$backfit_iteration <- NA_integer_
          cv_state$partial_index <- NA_integer_

          cv_progress_detail <- function(ridging = FALSE) {
            detail <- sprintf("multistart %d", cv_state$multistart_index)
            if (isTRUE(ridging)) {
              paste(detail, "near-singular system encountered, ridging", sep = ", ")
            } else {
              detail
            }
          }

          cv_progress_begin <- function() {
            cv_state$optim_eval <- 0L
            cv_state$optim_progress <- .np_progress_begin("Optimizing smooth coefficient bandwidth")
            invisible(NULL)
          }

          cv_progress_step <- function(ridging = FALSE) {
            cv_state$optim_eval <- cv_state$optim_eval + 1L
            .np_progress_bandwidth_activity_step(done = cv_state$optim_eval)
            cv_state$optim_progress <- .np_progress_step(
              state = cv_state$optim_progress,
              done = cv_state$optim_eval,
              detail = cv_progress_detail(ridging = ridging)
            )
            invisible(NULL)
          }

          cv_progress_end <- function(state) {
            if (is.null(state))
              return(invisible(NULL))

            if (isTRUE(state$known_total) && identical(state$last_done, state$total))
              return(invisible(NULL))

            state$last_emit <- -Inf
            .np_progress_end(state)
            invisible(NULL)
          }

          cv_progress_finish <- function(ridging = FALSE) {
            if (is.null(cv_state$optim_progress))
              return(invisible(NULL))

            cv_state$optim_progress$last_emit <- -Inf
            cv_state$optim_progress <- .np_progress_end(
              cv_state$optim_progress,
              detail = cv_progress_detail(ridging = ridging)
            )
            cv_state$optim_progress <- NULL
            invisible(NULL)
          }

          partial_progress_detail <- function(fv = NULL) {
            detail <- sprintf(
              "backfitting iteration %d of %d, partial residual %d of %d",
              cv_state$backfit_iteration,
              cv.num.iterations,
              cv_state$partial_index,
              ncol(W)
            )
            if (!is.null(fv)) {
              detail <- paste(
                detail,
                sprintf(
                  "fval %s",
                  format(signif(fv, digits = getOption("digits", 7L)), trim = TRUE)
                ),
                sep = ", "
              )
            }
            detail
          }

          partial_progress_begin <- function(iteration, partial.index) {
            cv_state$backfit_iteration <- iteration
            cv_state$partial_index <- partial.index
            cv_state$partial_eval <- 0L
            cv_state$partial_progress <- .np_progress_begin("Optimizing partial residual bandwidth")
            invisible(NULL)
          }

          partial_progress_step <- function(fv) {
            cv_state$partial_eval <- cv_state$partial_eval + 1L
            cv_state$partial_progress <- .np_progress_step(
              state = cv_state$partial_progress,
              done = cv_state$partial_eval,
              detail = partial_progress_detail(fv = fv)
            )
            invisible(NULL)
          }

          partial_progress_finish <- function(fv = NULL) {
            if (is.null(cv_state$partial_progress))
              return(invisible(NULL))

            cv_state$partial_progress$last_emit <- -Inf
            cv_state$partial_progress <- .np_progress_end(
              cv_state$partial_progress,
              detail = partial_progress_detail(fv = fv)
            )
            cv_state$partial_progress <- NULL
            invisible(NULL)
          }

          overall.cv.ls <- function(param) {
            cv_state$objective_fast <- FALSE
            sbw <- apply_bw_to_scbw(bws, param)
            if (!validateBandwidthTF(sbw) ||
                (!is.null(fixed.lower) && any(param < fixed.lower)) ||
                ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            cv_state$objective_fast <- npscoef_fast_eligible(sbw)

            if (identical(reg.engine, "lc")) {
              tww <- npksum(txdat = zdat, tydat = yW, weights = yW, bws = sbw,
                            leave.one.out = TRUE)$ksum

              mean.loo <- rep(maxPenalty,n)
              ridge.grid <- npRidgeSequenceAdditive(n.train = n, cap = 1.0)
              ridge <- rep.int(ridge.grid[1L], n)
              ridge.idx <- rep.int(1L, n)
              doridge <- rep.int(TRUE, n)

              nc <- ncol(tww[-1,-1,1])

              while(any(doridge)){
                iloo <- which(doridge)
                for (ii in iloo) {
                  doridge[ii] <- FALSE
                  ridge.val <- ridge[ii]*tww[-1,1,ii][1]/NZD(tww[-1,-1,ii][1,1])
                  beta.ii <- tryCatch(
                    solve(tww[-1,-1,ii] + diag(rep(ridge[ii], nc)),
                          tww[-1,1,ii] + c(ridge.val, rep(0, nc - 1))),
                    error = function(e) e
                  )
                  if (inherits(beta.ii, "error")) {
                    ridge.idx[ii] <- ridge.idx[ii] + 1L
                    if (ridge.idx[ii] <= length(ridge.grid)) {
                      ridge[ii] <- ridge.grid[ridge.idx[ii]]
                      doridge[ii] <- TRUE
                    }
                    beta.ii <- rep(maxPenalty, nc)
                  }
                  mean.loo[ii] <- W[ii,, drop = FALSE] %*% beta.ii
                }
              }
            } else {
              coef.loo <- lp_full_coef(sbw = sbw, leave.one.out.eval = TRUE)
              mean.loo <- rowSums(W * t(coef.loo))
            }

            stopifnot(all(is.finite(mean.loo)))

            if(!any(mean.loo == maxPenalty)){
              fv <- sum((ydat-mean.loo)^2)/n
              cv_progress_step()
            } else {
              cv_progress_step(ridging = TRUE)
              fv <- maxPenalty
            }

            if (isTRUE(cv_state$objective_fast))
              cv_state$fast_total <- cv_state$fast_total + 1L

            return((if (is.finite(fv)) fv else maxPenalty))

          }

          scoef.loo.args <- list(
            bws = bws, txdat = xdat, tydat = ydat,
            leave.one.out = TRUE, iterate = TRUE,
            maxiter = backfit.maxiter, tol = backfit.tol,
            betas = TRUE,
            .np_fit_progress_allow = FALSE
          )
          if (!miss.z)
            scoef.loo.args$tzdat <- zdat

          partial.cv.ls <- function(param, partial.index) {
            cv_state$objective_fast <- FALSE
            sbw <- apply_bw_to_scbw(bws, param)

            if (!validateBandwidthTF(sbw) ||
                (!is.null(fixed.lower) && any(param < fixed.lower)) ||
                ((bws$nord+bws$nuno > 0) && any(param[!bws$icon] > 2.0*x.scale[!bws$icon])))
              return(maxPenalty)
            cv_state$objective_fast <- npscoef_fast_eligible(sbw)

            if (backfit.iterate){
              bws$bw.fitted[,partial.index] <- param
              scoef.loo <- do.call(npscoef, scoef.loo.args)
              partial.loo <- W[,partial.index]*scoef.loo$beta[,partial.index]
            } else {
              wj <- W[,partial.index]
              if (identical(reg.engine, "lc")) {
                tww <- npksum(txdat=zdat,
                              tydat=cbind(partial.orig * wj, wj * wj),
                              weights=cbind(partial.orig * wj, 1),
                              bws=sbw,
                              leave.one.out=TRUE)$ksum

                partial.loo <- wj * tww[1,2,]/NZD(tww[2,2,])
              } else {
                partial.loo <- wj * lp_partial_coef(
                  sbw = sbw,
                  wj = wj,
                  partial.y = partial.orig,
                  leave.one.out.eval = TRUE
                )
              }
            }


            fv <- sum((partial.orig - partial.loo)^2)/n

            if (isTRUE(cv_state$objective_fast))
              cv_state$fast_total <- cv_state$fast_total + 1L

            partial_progress_step(fv = fv)
            return((if (is.finite(fv)) fv else maxPenalty))
          }

          ## Now we implement multistarting

          fval.min <- .Machine$double.xmax
          have_best <- FALSE
          numimp <- 0
          value.overall <- numeric(nmulti)
          num.feval.overall <- 0

          x.scale <- sapply(seq_len(bws$ndim), function(i){
            if (dati$icon[i]){
              return(1.059224*((if (bws$scaling) 1.0 else mysd[sum(dati$icon[seq_len(i)])]*nconfac)))
            }

            if (dati$iord[i])
              return(0.5*oMaxL(dati$all.nlev[[i]], kertype = bws$okertype)*
                     (if (bws$scaling) ncatfac else 1.0))

            if (dati$iuno[i])
              return(0.5*uMaxL(dati$all.nlev[[i]], kertype = bws$ukertype)*
                     (if (bws$scaling) ncatfac else 1.0))
          })
          fixed.lower <- if (identical(bws$type, "fixed")) {
            out <- rep.int(0, length(x.scale))
            out[dati$icon] <- x.scale[dati$icon] * start.controls$scale.factor.search.lower
            out
          } else {
            NULL
          }

          optim.control <- list(abstol = optim.abstol,
                                reltol = optim.reltol,
                                maxit = optim.maxit)

          for (i in seq_len(nmulti)) {

            cv_state$multistart_index <- i
            cv_progress_begin()

            if (i == 1) {
              tbw <- .npscoef_default_start_bandwidth(
                param = x.scale,
                bwtype = bws$type,
                nobs = n,
                start.controls = start.controls,
                icon = dati$icon,
                iord = dati$iord,
                iuno = dati$iuno
              )
              if (all(bws$bw != 0) &&
                  .npscoef_candidate_is_admissible(param = bws$bw, bwtype = bws$type, nobs = n,
                                                   lower = fixed.lower)) {
                tbw <- .npscoef_finalize_bandwidth(
                  param = bws$bw,
                  bwtype = bws$type,
                  nobs = n,
                  lower = fixed.lower,
                  where = "npscoefbw"
                )
              }
            } else {
              tbw <- .npscoef_random_start_bandwidth(
                param = x.scale,
                bwtype = bws$type,
                nobs = n,
                start.controls = start.controls,
                icon = dati$icon,
                iord = dati$iord,
                iuno = dati$iuno
              )
            }

            suppressWarnings(optim.return <- optim(tbw,
                                                   fn = overall.cv.ls,
                                                   control = optim.control))
            if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
              num.feval.overall <- num.feval.overall + optim.return$counts[1]
            attempts <- 0
            while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
              attempts <- attempts + 1
              tbw <- .npscoef_random_start_bandwidth(
                param = x.scale,
                bwtype = bws$type,
                nobs = n,
                start.controls = start.controls,
                icon = dati$icon,
                iord = dati$iord,
                iuno = dati$iuno
              )
              optim.control <- lapply(optim.control, '*', 10.0)
              suppressWarnings(optim.return <- optim(tbw,
                                                     fn = overall.cv.ls,
                                                     control = optim.control))
              if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                num.feval.overall <- num.feval.overall + optim.return$counts[1]

            }

            cv_progress_finish()

            value.overall[i] <- optim.return$value

            if (.npscoef_candidate_is_admissible(
              param = optim.return$par,
              bwtype = bws$type,
              nobs = n,
              lower = fixed.lower
            ) && (!have_best || optim.return$value < fval.min)) {
              param <- .npscoef_finalize_bandwidth(
                param = optim.return$par,
                bwtype = bws$type,
                nobs = n,
                lower = fixed.lower,
                where = "npscoefbw"
              )
              min.overall <- optim.return$value
              fval.min <- min.overall ## Added by jracine Jul 22 2010
              numimp.overall <- numimp + 1
              best.overall <- i
              have_best <- TRUE
            }

            .np_progress_bandwidth_multistart_step(done = i, total = nmulti)
          }

          if (!have_best) {
            if (identical(bws$type, "fixed")) {
              stop("npscoefbw: no feasible fixed bandwidths found", call. = FALSE)
            }
            stop("npscoefbw: no feasible bandwidths found", call. = FALSE)
          }

          param.overall <- bws$bw <- .npscoef_finalize_bandwidth(
            param = param,
            bwtype = bws$type,
            nobs = n,
            where = "npscoefbw"
          )
          bws <- apply_bw_to_scbw(bws, bws$bw)

          if(cv.iterate){
            n.part <- (ncol(xdat)+1)
            backfit.progress <- .np_progress_begin(
              "Backfitting smooth coefficient bandwidth",
              total = cv.num.iterations
            )
            on.exit(cv_progress_end(backfit.progress), add = TRUE)
            bws$bw.fitted <- matrix(data = bws$bw, nrow = length(bws$bw), ncol = n.part)
            ## obtain matrix of alpha.hat | h0 and beta.hat | h0

            scoef.args <- list(
              bws = bws,
              txdat = xdat,
              tydat = ydat,
              iterate = FALSE,
              betas = TRUE,
              .np_fit_progress_allow = FALSE
            )
            if (!miss.z)
              scoef.args$tzdat <- zdat
            scoef <- do.call(npscoef, scoef.args)

            resid.full <- ydat - scoef$mean


            for (i in seq_len(cv.num.iterations)) {
              backfit.progress <- .np_progress_step(
                state = backfit.progress,
                done = i,
                detail = sprintf("iteration %d of %d", i, cv.num.iterations)
              )

              for (j in seq_len(n.part)) {
                ## estimate partial residuals
                partial.orig <- W[,j] * scoef$beta[,j] + resid.full
                partial_progress_begin(iteration = i, partial.index = j)

                ## minimise
                suppressWarnings(optim.return <-
                                 optim(tbw, fn = partial.cv.ls,
                                       control = optim.control,
                                       partial.index = j))
                if(!is.null(optim.return$counts) && length(optim.return$counts) > 0)
                  num.feval.overall <- num.feval.overall + optim.return$counts[1]
                partial_progress_finish(fv = optim.return$value)

                ## grab parameter
                bws$bw.fitted[,j] <- optim.return$par

                if (backfit.iterate){
                  ## re-estimate all betas
                  scoef.args <- list(
                    bws = bws, txdat = xdat, tydat = ydat,
                    iterate = TRUE, maxiter = backfit.maxiter,
                    tol = backfit.tol, betas = TRUE,
                    .np_fit_progress_allow = FALSE
                  )
                  if (!miss.z)
                    scoef.args$tzdat <- zdat
                  scoef <- do.call(npscoef, scoef.args)
                  resid.full <- ydat - scoef$mean
                } else {
                  bws$bw <- bws$bw.fitted[,j]
                  ## estimate new beta.hats

                  bws <- apply_bw_to_scbw(bws, bws$bw)

                  if (identical(reg.engine, "lc")) {
                    wj <- W[,j]
                    tww <- npksum(txdat=zdat,
                                  tydat=cbind(partial.orig * wj, wj * wj),
                                  weights=cbind(partial.orig * wj, 1),
                                  bws=bws)$ksum
                    scoef$beta[,j] <- tww[1,2,]/NZD(tww[2,2,])
                  } else {
                    wj <- W[,j]
                    scoef$beta[,j] <- lp_partial_coef(
                      sbw = bws,
                      wj = wj,
                      partial.y = partial.orig,
                      leave.one.out.eval = FALSE
                    )
                  }

                  bws$bw <- param.overall
                  bws <- apply_bw_to_scbw(bws, bws$bw)
                  ## estimate new full residuals
                  resid.full <- partial.orig - W[,j] * scoef$beta[,j]
                }
              }
            }
            scoef.loo.args <- list(
              bws = bws, txdat = xdat, tydat = ydat,
              iterate = TRUE, maxiter = backfit.maxiter,
              tol = backfit.tol, leave.one.out = TRUE,
              .np_fit_progress_allow = FALSE
            )
            if (!miss.z)
              scoef.loo.args$tzdat <- zdat
            scoef.loo <- do.call(npscoef, scoef.loo.args)$mean
            bws$fval.fitted <- sum((ydat - scoef.loo)^2)/n
          }

          bws$fval = min.overall
          bws$ifval = best.overall
          bws$num.feval = num.feval.overall
          bws$num.feval.fast = cv_state$fast_total
          bws$numimp = numimp.overall
          bws$fval.vector = value.overall
        }
      })[["elapsed"]]

    bws$sfactor <- bws$bandwidth <- bws$bw
    nfactor <- nrow^(-2.0/(2.0*bws$ckerorder+bws$ncon))
    dfactor <- EssDee(zdat[, dati$icon, drop = FALSE])*nrow^(-1.0/(2.0*bws$ckerorder+sum(dati$icon)))

    if (bws$scaling) {
      bws$bandwidth[dati$icon] <- bws$bandwidth[dati$icon]*dfactor

      if(bws$nuno > 0)
        bws$bandwidth[dati$iuno] <- bws$bandwidth[dati$iuno]*nfactor

      if(bws$nord > 0)
        bws$bandwidth[dati$iord] <- bws$bandwidth[dati$iord]*nfactor

    } else {
      bws$sfactor[dati$icon] <- bws$sfactor[dati$icon]/dfactor

      if(bws$nuno > 0)
        bws$sfactor[dati$iuno] <- bws$sfactor[dati$iuno]/nfactor

      if(bws$nord > 0)
        bws$sfactor[dati$iord] <- bws$sfactor[dati$iord]/nfactor
    }

    ## Restore seed

    .np_seed_exit(seed.state)

    bws <- scbandwidth(bw = bws$bw,
                       regtype = regtype,
                       basis = if (is.null(bws$basis)) "glp" else bws$basis,
                       degree = bws$degree,
                       bernstein.basis = bws$bernstein.basis,
                       bwmethod = bws$method,
                       bwscaling = bws$scaling,
                       bwtype = bws$type,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       ckerbound = bws$ckerbound,
                       ckerlb = bws$ckerlb,
                       ckerub = bws$ckerub,
                       ukertype = bws$ukertype,
                       okertype = bws$okertype,
                       fval = bws$fval,
                       ifval = bws$ifval,
                       num.feval = bws$num.feval,
                       num.feval.fast = bws$num.feval.fast,
                       numimp = bws$numimp,
                       fval.vector = bws$fval.vector,
                       nobs = bws$nobs,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       zdati = bws$zdati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       znames = bws$znames,
                       sfactor = bws$sfactor,
                       bandwidth = bws$bandwidth,
                       rows.omit = rows.omit,
                       bandwidth.compute = bandwidth.compute,
                       optim.method = optim.method,
                       total.time = total.time)
    bws <- npSetScaleFactorSearchLower(bws, scale.factor.search.lower)

    bws
  }

npscoefbw.default <-
  function(xdat = stop("invoked without data 'xdat'"),
           ydat = stop("invoked without data 'ydat'"),
           zdat = NULL,
           bws,
           backfit.iterate,
           backfit.maxiter,
           backfit.tol,
           bandwidth.compute = TRUE,
           basis,
           bernstein.basis,
           bwmethod,
           bwscaling,
           bwtype,
           ckerbound,
           ckerlb,
           ckerorder,
           ckertype,
           ckerub,
           cv.iterate,
           cv.num.iterations,
           degree,
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
           okertype,
           optim.abstol,
           optim.maxattempts,
           optim.maxit,
           optim.method,
           optim.reltol,
           random.seed,
           regtype,
           ukertype,
           scale.factor.init.lower = 0.1,
           scale.factor.init.upper = 2.0,
           scale.factor.init = 0.5,
           lbd.init = 0.5,
           hbd.init = 1.5,
           dfac.init = 1.0,
           scale.factor.search.lower = NULL,
           ...){
    .npRmpi_require_active_slave_pool(where = "npscoefbw()")
    if (!missing(bwmethod) && identical(match.arg(bwmethod, c("cv.ls", "manual")), "manual") &&
        missing(bws))
      stop("bwmethod='manual' requires argument 'bws'")

    miss.z <- missing(zdat)
    xdat <- toFrame(xdat)

    if (!(is.vector(ydat) || is.factor(ydat)))
      stop("'ydat' must be a vector or a factor")

    if(!miss.z)
      zdat <- toFrame(zdat)

    mc <- match.call(expand.dots = FALSE)
    mc.names <- names(mc)
    dots <- list(...)
    nomad.shortcut <- .np_prepare_nomad_shortcut(
      nomad = nomad,
      call_names = mc.names,
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
        bwtype = if ("bwtype" %in% mc.names) bwtype else NULL,
        degree = if ("degree" %in% mc.names) degree else NULL
      ),
      where = "npscoefbw"
    )

    if (isTRUE(nomad.shortcut$enabled)) {
      if ("degree" %in% mc.names)
        stop("nomad=TRUE does not support an explicit degree; remove degree or set nomad=FALSE")
      if ("regtype" %in% mc.names &&
          !identical(as.character(match.arg(nomad.shortcut$values$regtype, c("lc", "ll", "lp")))[1L], "lp"))
        stop("nomad=TRUE requires regtype='lp'")
      if ("bwtype" %in% mc.names &&
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

    regtype.named <- isTRUE(nomad.shortcut$enabled) || any(mc.names == "regtype")
    bernstein.named <- isTRUE(nomad.shortcut$enabled) || any(mc.names == "bernstein.basis")
    cv.iterate.named <- any(mc.names == "cv.iterate")
    regtype.value <- if (!is.null(nomad.shortcut$values$regtype)) nomad.shortcut$values$regtype else "lc"
    degree.select.value <- if (!is.null(nomad.shortcut$values$degree.select)) nomad.shortcut$values$degree.select else "manual"
    degree.search <- .npscoefbw_degree_search_controls(
      regtype = regtype.value,
      regtype.named = regtype.named,
      cv.iterate = cv.iterate,
      cv.iterate.named = cv.iterate.named,
      bandwidth.compute = bandwidth.compute,
      ncon = sum(if (miss.z) untangle(xdat)$icon else untangle(zdat)$icon),
      nobs = if (miss.z) NROW(xdat) else NROW(zdat),
      basis = if ("basis" %in% mc.names) basis else "glp",
      degree.select = degree.select.value,
      search.engine = if (!is.null(nomad.shortcut$values$search.engine)) nomad.shortcut$values$search.engine else "nomad+powell",
      degree.min = nomad.shortcut$values$degree.min,
      degree.max = nomad.shortcut$values$degree.max,
      degree.start = if ("degree.start" %in% mc.names) degree.start else NULL,
      degree.restarts = if ("degree.restarts" %in% mc.names) degree.restarts else 0L,
      degree.max.cycles = if ("degree.max.cycles" %in% mc.names) degree.max.cycles else 20L,
      degree.verify = if (!is.null(nomad.shortcut$values$degree.verify)) nomad.shortcut$values$degree.verify else FALSE,
      bernstein.basis = if (!is.null(nomad.shortcut$values$bernstein.basis)) nomad.shortcut$values$bernstein.basis else bernstein.basis,
      bernstein.named = bernstein.named
    )
    nomad.inner <- .np_nomad_validate_inner_multistart(
      call_names = mc.names,
      dot.args = dots,
      nomad.nmulti = nomad.nmulti,
      regtype = regtype.value,
      automatic.degree.search = !is.null(degree.search),
      search.engine = if (is.null(degree.search)) "" else degree.search$engine
    )
    nomad.inner.named <- nomad.inner$named
    nomad.inner.nmulti <- nomad.inner$nmulti
    degree.setup <- npSetupGlpDegree(
      regtype = regtype.value,
      degree = if ("degree" %in% mc.names) degree else NULL,
      ncon = sum(if (miss.z) untangle(xdat)$icon else untangle(zdat)$icon),
      degree.select = degree.select.value
    )
    scale.factor.search.lower <- npResolveScaleFactorLowerBound(scale.factor.search.lower)
    if (.npRmpi_autodispatch_active() && is.null(degree.search))
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    ## first grab dummy args for scbandwidth() and perform 'bootstrap'
    ## bandwidth call

    margs <- c("regtype", "basis", "degree", "bernstein.basis",
               "bwmethod", "bwscaling", "bwtype", "ckertype", "ckerorder",
               "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype",
               "scale.factor.search.lower")

    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    sbw.args <- list(
      bw = bws,
      nobs = dim(xdat)[1],
      xdati = untangle(xdat),
      ydati = untangle(data.frame(ydat)),
      zdati = untangle(zdat),
      xnames = names(xdat),
      ynames = deparse(substitute(ydat)),
      znames = names(zdat),
      bandwidth.compute = bandwidth.compute
    )
    if (any.m) {
      nms <- mc.names[m]
      sbw.args[nms] <- mget(nms, envir = environment(), inherits = FALSE)
    }
    if (!("degree" %in% names(sbw.args)) && !is.null(degree.setup))
      sbw.args$degree <- degree.setup
    reg.args <- sbw.args[setdiff(names(sbw.args), c("bw", "nobs", "xdati", "ydati", "zdati", "xnames", "ynames", "znames", "bandwidth.compute"))]
    if (!is.null(degree.search))
      reg.args$bernstein.basis <- degree.search$bernstein.basis
    tbw <- do.call(scbandwidth, sbw.args)
    tbw <- npSetScaleFactorSearchLower(tbw, scale.factor.search.lower)

    ## next grab dummies for actual bandwidth selection and perform call
    margs <- c("zdat",
               "nmulti",
               "random.seed",
               "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
               "lbd.init", "hbd.init", "dfac.init",
               "scale.factor.search.lower",
               "cv.iterate",
               "cv.num.iterations",
               "backfit.iterate",
               "backfit.maxiter",
               "backfit.tol",
               "optim.method", "optim.maxattempts",
               "optim.reltol", "optim.abstol", "optim.maxit")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    if (any.m) {
      nms <- mc.names[m]
      opt.args <- mget(nms, envir = environment(), inherits = FALSE)
    } else {
      opt.args <- list()
    }
    opt.args <- c(list(bandwidth.compute = bandwidth.compute), opt.args)
    reg.args$scale.factor.search.lower <- scale.factor.search.lower
    opt.args$scale.factor.search.lower <- scale.factor.search.lower

    if (!is.null(degree.search)) {
      if (identical(degree.search$engine, "cell")) {
        eval_fun <- function(degree.vec) {
          cell.reg.args <- reg.args
          cell.reg.args$regtype <- "lp"
          cell.reg.args$degree <- as.integer(degree.vec)
          cell.reg.args$bernstein.basis <- degree.search$bernstein.basis
          cell.bws <- .npscoefbw_run_fixed_degree(
            xdat = xdat,
            ydat = ydat,
            zdat = if (miss.z) NULL else zdat,
            bws = bws,
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
        search.result <- .npscoefbw_nomad_search(
          xdat = xdat,
          ydat = ydat,
          zdat = if (miss.z) NULL else zdat,
          bws = bws,
          reg.args = reg.args,
          opt.args = opt.args,
          degree.search = degree.search,
          nomad.inner.nmulti = nomad.inner.nmulti
        )
      }
      tbw <- .npscoefbw_attach_degree_search(
        bws = search.result$best_payload,
        search_result = search.result
      )
    } else {
      scbw.args <- c(list(xdat = xdat, ydat = ydat, bws = tbw), opt.args)
      if (!miss.z)
        scbw.args$zdat <- zdat
      tbw <- .np_progress_select_bandwidth_enhanced(
        "Selecting smooth coefficient bandwidth",
        do.call(npscoefbw.scbandwidth, scbw.args)
      )
    }

    environment(mc) <- parent.frame()
    tbw$call <- mc
    tbw <- .np_attach_nomad_shortcut(tbw, nomad.shortcut$metadata)

    return(tbw)

  }
