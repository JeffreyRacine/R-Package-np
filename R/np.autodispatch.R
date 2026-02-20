.npRmpi_autodispatch_option_keys <- function() {
  c("np.messages", "np.tree", "np.largeh.rel.tol",
    "np.disc.upper.rel.tol", "np.groupcv.fast")
}

.npRmpi_warn_pkg_conflict_once <- function() {
  if (!isTRUE(getOption("npRmpi.conflicts.warn", TRUE)))
    return(invisible(FALSE))
  if (!("package:np" %in% search()))
    return(invisible(FALSE))
  if (isTRUE(getOption("npRmpi.conflicts.warned", FALSE)))
    return(invisible(TRUE))
  warning("both packages 'npRmpi' and 'np' are attached: use explicit npRmpi:: calls to avoid masked-function ambiguity")
  options(npRmpi.conflicts.warned = TRUE)
  invisible(TRUE)
}

.npRmpi_warn_rmpi_conflict_once <- function() {
  if (!isTRUE(getOption("npRmpi.conflicts.warn", TRUE)))
    return(invisible(FALSE))
  if (!("package:Rmpi" %in% search()))
    return(invisible(FALSE))
  if (isTRUE(getOption("npRmpi.conflicts.warned.rmpi", FALSE)))
    return(invisible(TRUE))
  warning("both packages 'npRmpi' and 'Rmpi' are attached: prefer explicit npRmpi:: calls in user code to avoid dispatch ambiguity")
  options(npRmpi.conflicts.warned.rmpi = TRUE)
  invisible(TRUE)
}

.npRmpi_bcast_cmd_expr <- function(expr, comm = 1L, caller.execute = TRUE) {
  eval(substitute(mpi.bcast.cmd(cmd = EXPR, comm = COMM, caller.execute = CE),
                  list(EXPR = expr, COMM = comm, CE = caller.execute)),
       envir = parent.frame())
}

.npRmpi_bcast_robj_by_name <- function(name, caller_env = parent.frame()) {
  expr <- parse(text = sprintf("mpi.bcast.Robj2slave(%s)", name))[[1L]]
  eval(expr, envir = caller_env)
}

.npRmpi_autodispatch_active <- function() {
  isTRUE(getOption("npRmpi.autodispatch", FALSE)) &&
    !isTRUE(getOption("npRmpi.autodispatch.disable", FALSE))
}

.npRmpi_autodispatch_in_context <- function() {
  isTRUE(getOption("npRmpi.autodispatch.context", FALSE))
}

.npRmpi_eval_without_dispatch <- function(mc, caller_env) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  mc.eval <- mc
  if (is.call(mc.eval) && length(mc.eval) >= 1L && is.symbol(mc.eval[[1L]])) {
    fname <- as.character(mc.eval[[1L]])
    in.caller <- exists(fname, envir = caller_env, mode = "function", inherits = TRUE)
    in.ns <- exists(fname, envir = asNamespace("npRmpi"), mode = "function", inherits = FALSE)
    if (!in.caller && in.ns) {
      mc.eval[[1L]] <- get(fname, envir = asNamespace("npRmpi"), mode = "function", inherits = FALSE)
    }
  }
  eval(mc.eval, envir = caller_env)
}

.npRmpi_autodispatch_called_from_bcast <- function() {
  calls <- sys.calls()
  if (!length(calls)) return(FALSE)
  any(vapply(calls, function(cl) {
    if (!is.call(cl) || length(cl) < 1) return(FALSE)
    fn <- cl[[1]]
    if (is.symbol(fn) && identical(as.character(fn), "mpi.bcast.cmd"))
      return(TRUE)
    if (is.call(fn) && length(fn) >= 3 && is.symbol(fn[[1]]) &&
        as.character(fn[[1]]) %in% c("::", ":::")) {
      tgt <- fn[[3]]
      return(is.symbol(tgt) && identical(as.character(tgt), "mpi.bcast.cmd"))
    }
    FALSE
  }, logical(1)))
}

.npRmpi_autodispatch_warn_nested <- function() {
  if (!isTRUE(getOption("npRmpi.autodispatch.warned.nested", FALSE))) {
    warning("detected active mpi.bcast.cmd context; skipping nested auto-dispatch for this call")
    options(npRmpi.autodispatch.warned.nested = TRUE)
  }
}

.npRmpi_has_active_slave_pool <- function(comm = 1L) {
  size <- try(mpi.comm.size(comm), silent = TRUE)
  rank <- try(mpi.comm.rank(comm), silent = TRUE)
  if (inherits(size, "try-error") || inherits(rank, "try-error") ||
      is.na(size) || is.na(rank))
    return(FALSE)
  size >= 2L
}

.npRmpi_require_active_slave_pool <- function(comm = 1L, where = "this call") {
  if (.npRmpi_has_active_slave_pool(comm = comm))
    return(invisible(TRUE))
  stop(sprintf("%s requires an active MPI slave pool; call npRmpi.start(nslaves=...) first", where))
}

.npRmpi_autodispatch_preflight <- function(comm = 1L) {
  strict <- isTRUE(getOption("npRmpi.autodispatch.strict", TRUE))
  if (!.npRmpi_has_active_slave_pool(comm = comm)) {
    msg <- "npRmpi auto-dispatch requires an active slave pool; call npRmpi.start(nslaves=...) first"
    if (strict) stop(msg)
    warning(msg)
    return(FALSE)
  }

  TRUE
}

.npRmpi_autodispatch_sync_options <- function(comm = 1L) {
  keys <- .npRmpi_autodispatch_option_keys()
  vals <- lapply(keys, getOption)
  cmd.sync <- substitute({
    for (i in seq_along(KEYS)) {
      options(structure(list(VALS[[i]]), names = KEYS[[i]]))
    }
  }, list(KEYS = keys, VALS = vals))
  .npRmpi_bcast_cmd_expr(cmd.sync, comm = comm, caller.execute = TRUE)

  if (isTRUE(getOption("npRmpi.autodispatch.strict", TRUE))) {
    for (i in seq_along(keys)) {
      k <- keys[[i]]
      mval <- vals[[i]]
      sval <- mpi.remote.exec(getOption, k, simplify = TRUE, comm = comm, ret = TRUE)
      svals <- unname(unlist(sval, recursive = TRUE, use.names = FALSE))
      if (!all(vapply(as.list(svals), identical, logical(1), mval)))
        stop(sprintf("failed to synchronize option '%s' across MPI ranks", k))
    }
  }

  invisible(TRUE)
}

.npRmpi_autodispatch_target_args <- function() {
  c("formula", "data", "bws",
    "dat", "tdat", "edat",
    "xdat", "ydat", "zdat", "txdat", "tydat", "tzdat",
    "exdat", "eydat", "ezdat", "newdata",
    "gydat", "wdat", "gdata",
    "data.x", "data.y", "model",
    "y", "z", "w", "x", "zeval", "weval", "xeval", "bw",
    "gradients", "residuals", "errors", "gradient.order")
}

.npRmpi_autodispatch_failfast_formula_data <- function(mc, caller_env) {
  invisible(FALSE)
}

.npRmpi_autodispatch_eval_arg <- function(expr, caller_env) {
  val <- try(eval(expr, envir = caller_env), silent = TRUE)
  if (!inherits(val, "try-error"))
    return(val)

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    env_i <- frames[[i]]
    if (identical(env_i, caller_env))
      next
    val_i <- try(eval(expr, envir = env_i), silent = TRUE)
    if (!inherits(val_i, "try-error"))
      return(val_i)
  }

  stop(attr(val, "condition")$message)
}

.npRmpi_autodispatch_materialize_call <- function(mc, caller_env, comm = 1L) {
  arg.list <- as.list(mc)
  nms <- names(arg.list)
  targets <- .npRmpi_autodispatch_target_args()

  out <- mc
  if (is.call(out) && length(out) >= 1L && is.symbol(out[[1L]])) {
    fname <- as.character(out[[1L]])
    if (exists(fname, envir = asNamespace("npRmpi"), mode = "function", inherits = FALSE)) {
      out[[1L]] <- get(fname, envir = asNamespace("npRmpi"), mode = "function", inherits = FALSE)
    }
  }
  tmpnames <- character(0)
  tmpvals <- list()
  idx <- 0L

  has_data_inputs <- !is.null(nms) && any(nms %in% c("data", "xdat", "ydat", "txdat", "tydat", "zdat"))

  formula.expr <- NULL
  if (!is.null(nms) && any(nms == "formula")) {
    formula.expr <- arg.list[[which(nms == "formula")[1L]]]
  } else if (!is.null(nms) && any(nms == "bws")) {
    bexpr <- arg.list[[which(nms == "bws")[1L]]]
    bval <- try(.npRmpi_autodispatch_eval_arg(bexpr, caller_env = caller_env), silent = TRUE)
    if (!inherits(bval, "try-error") && inherits(bval, "formula"))
      formula.expr <- bexpr
  } else if (length(arg.list) >= 2L) {
    nm2 <- if (!is.null(nms) && length(nms) >= 2L) nms[[2L]] else ""
    if (is.null(nm2) || identical(nm2, "")) {
      fval <- try(.npRmpi_autodispatch_eval_arg(arg.list[[2L]], caller_env = caller_env), silent = TRUE)
      if (!inherits(fval, "try-error") && inherits(fval, "formula"))
        formula.expr <- arg.list[[2L]]
    }
  }

  if (!has_data_inputs && !is.null(formula.expr)) {
    fval <- .npRmpi_autodispatch_eval_arg(formula.expr, caller_env = caller_env)
    vars <- all.vars(fval)
    if (length(vars)) {
      dlist <- setNames(lapply(vars, function(v) {
        .npRmpi_autodispatch_eval_arg(as.name(v), caller_env = caller_env)
      }), vars)
      idx <- idx + 1L
      tmp <- sprintf(".__npRmpi_autod_data_%d", idx)
      out[["data"]] <- as.name(tmp)
      tmpnames <- c(tmpnames, tmp)
      tmpvals[[tmp]] <- as.data.frame(dlist, stringsAsFactors = FALSE)
    }
  }

  for (i in seq_along(arg.list)) {
    if (i == 1L) next
    nm <- nms[i]
    if (is.null(nm) || identical(nm, "")) next
    if (!nm %in% targets) next

    val <- .npRmpi_autodispatch_eval_arg(arg.list[[i]], caller_env = caller_env)
    idx <- idx + 1L
    tmp <- sprintf(".__npRmpi_autod_%s_%d", nm, idx)
    .GlobalEnv[[tmp]] <- val
    .npRmpi_bcast_robj_by_name(tmp, caller_env = .GlobalEnv)

    out[[i]] <- as.name(tmp)
    tmpnames <- c(tmpnames, tmp)
    tmpvals[[tmp]] <- val
  }

  list(call = out, tmpnames = unique(tmpnames), tmpvals = tmpvals)
}

.npRmpi_autodispatch_cleanup <- function(tmpnames, comm = 1L) {
  if (!length(tmpnames))
    return(invisible(TRUE))
  rm(list = tmpnames, envir = .GlobalEnv)
  cmd.rm <- substitute(rm(list = TMPS, envir = .GlobalEnv),
                       list(TMPS = tmpnames))
  .npRmpi_bcast_cmd_expr(cmd.rm, comm = comm, caller.execute = FALSE)
  invisible(TRUE)
}

.npRmpi_autodispatch_replace_tmps <- function(x, tmpvals) {
  if (!length(tmpvals))
    return(x)

  if (is.symbol(x)) {
    nm <- as.character(x)
    if (!is.null(tmpvals[[nm]]))
      return(tmpvals[[nm]])
    return(x)
  }

  if (is.call(x)) {
    for (i in seq_len(length(x))) {
      xi <- try(x[[i]], silent = TRUE)
      if (!inherits(xi, "try-error"))
        x[[i]] <- .npRmpi_autodispatch_replace_tmps(xi, tmpvals = tmpvals)
    }
    return(x)
  }

  if (is.pairlist(x)) {
    for (i in seq_len(length(x))) {
      xi <- try(x[[i]], silent = TRUE)
      if (!inherits(xi, "try-error"))
        x[[i]] <- .npRmpi_autodispatch_replace_tmps(xi, tmpvals = tmpvals)
    }
    return(x)
  }

  if (inherits(x, "formula"))
    return(x)

  x
}

.npRmpi_autodispatch_replace_tmp_calls <- function(x, tmpvals) {
  if (!is.list(x) || !length(tmpvals))
    return(x)

  nms <- names(x)
  for (i in seq_along(x)) {
    xi <- x[[i]]
    nm <- if (!is.null(nms)) nms[[i]] else ""

    if (is.list(xi)) {
      x[i] <- list(.npRmpi_autodispatch_replace_tmp_calls(xi, tmpvals = tmpvals))
      next
    }

    if ((is.call(xi) || is.pairlist(xi)) && identical(nm, "call")) {
      x[i] <- list(.npRmpi_autodispatch_replace_tmps(xi, tmpvals = tmpvals))
    }
  }

  x
}

.npRmpi_autodispatch_tag_result <- function(x, mode = "auto") {
  if (is.list(x) || is.environment(x)) {
    attr(x, "npRmpi.dispatch.mode") <- mode
  }
  x
}

.npRmpi_autodispatch_untag <- function(x) {
  if (is.list(x) && !is.pairlist(x) && !is.call(x)) {
    for (i in seq_len(length(x))) {
      xi <- try(x[[i]], silent = TRUE)
      if (!inherits(xi, "try-error"))
        x[i] <- list(.npRmpi_autodispatch_untag(xi))
    }
  }
  attr(x, "npRmpi.dispatch.mode") <- NULL
  x
}

.npRmpi_guard_no_auto_object_in_manual_bcast <- function(obj, where = "this call") {
  if (.npRmpi_autodispatch_in_context())
    return(invisible(FALSE))
  if (!.npRmpi_autodispatch_called_from_bcast())
    return(invisible(FALSE))
  mode <- try(attr(obj, "npRmpi.dispatch.mode", exact = TRUE), silent = TRUE)
  if (inherits(mode, "try-error") || is.null(mode))
    return(invisible(FALSE))
  if (identical(mode, "auto")) {
    stop(sprintf("%s received an object created under npRmpi.autodispatch and cannot be executed inside mpi.bcast.cmd(...): avoid mixing dispatch modes; either rerun the full workflow in manual broadcast mode or keep this call outside mpi.bcast.cmd", where))
  }
  invisible(FALSE)
}

.npRmpi_autodispatch_call <- function(mc, caller_env = parent.frame(), comm = 1L) {
  .npRmpi_warn_pkg_conflict_once()
  .npRmpi_warn_rmpi_conflict_once()
  if (!.npRmpi_autodispatch_active())
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  if (.npRmpi_autodispatch_in_context())
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  rank <- try(mpi.comm.rank(comm), silent = TRUE)
  if (!inherits(rank, "try-error") && !is.na(rank) && rank != 0L)
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  if (.npRmpi_autodispatch_called_from_bcast()) {
    .npRmpi_autodispatch_warn_nested()
    return(.npRmpi_eval_without_dispatch(mc, caller_env))
  }

  .npRmpi_autodispatch_failfast_formula_data(mc, caller_env = caller_env)

  if (!.npRmpi_autodispatch_preflight(comm = comm))
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  .npRmpi_autodispatch_sync_options(comm = comm)
  prepared <- .npRmpi_autodispatch_materialize_call(mc = mc, caller_env = caller_env, comm = comm)
  on.exit(.npRmpi_autodispatch_cleanup(prepared$tmpnames, comm = comm), add = TRUE)

  if (length(prepared$tmpnames)) {
    for (nm in prepared$tmpnames) {
      val <- .npRmpi_autodispatch_untag(prepared$tmpvals[[nm]])
      prepared$tmpvals[[nm]] <- val
      .GlobalEnv[[nm]] <- val
      .npRmpi_bcast_robj_by_name(nm, caller_env = .GlobalEnv)
    }
  }

  cmd <- substitute({
    old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
    old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
    options(npRmpi.autodispatch.context = TRUE)
    options(npRmpi.autodispatch.disable = TRUE)
    on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
    on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
    eval(CALL, envir = .GlobalEnv)
  }, list(CALL = prepared$call))

  result <- .npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE)

  if (is.list(result))
    return(.npRmpi_autodispatch_tag_result(.npRmpi_autodispatch_replace_tmp_calls(result, tmpvals = prepared$tmpvals), mode = "auto"))

  if (is.call(result))
    return(.npRmpi_autodispatch_tag_result(.npRmpi_autodispatch_replace_tmps(result, tmpvals = prepared$tmpvals), mode = "auto"))

  .npRmpi_autodispatch_tag_result(result, mode = "auto")
}

.npRmpi_manual_distributed_call <- function(mc, caller_env = parent.frame(), comm = 1L) {
  .npRmpi_warn_pkg_conflict_once()
  .npRmpi_warn_rmpi_conflict_once()

  if (.npRmpi_autodispatch_in_context())
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  rank <- try(mpi.comm.rank(comm), silent = TRUE)
  if (!inherits(rank, "try-error") && !is.na(rank) && rank != 0L)
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  if (.npRmpi_autodispatch_called_from_bcast())
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  .npRmpi_autodispatch_failfast_formula_data(mc, caller_env = caller_env)

  if (!.npRmpi_autodispatch_preflight(comm = comm))
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  .npRmpi_autodispatch_sync_options(comm = comm)
  prepared <- .npRmpi_autodispatch_materialize_call(mc = mc, caller_env = caller_env, comm = comm)
  on.exit(.npRmpi_autodispatch_cleanup(prepared$tmpnames, comm = comm), add = TRUE)

  if (length(prepared$tmpnames)) {
    for (nm in prepared$tmpnames) {
      val <- .npRmpi_autodispatch_untag(prepared$tmpvals[[nm]])
      prepared$tmpvals[[nm]] <- val
      .GlobalEnv[[nm]] <- val
      .npRmpi_bcast_robj_by_name(nm, caller_env = .GlobalEnv)
    }
  }

  cmd <- substitute({
    old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
    old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
    options(npRmpi.autodispatch.context = TRUE)
    options(npRmpi.autodispatch.disable = TRUE)
    on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
    on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
    CALL
  }, list(CALL = prepared$call))

  result <- .npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE)

  if (is.list(result))
    return(.npRmpi_autodispatch_tag_result(.npRmpi_autodispatch_replace_tmp_calls(result, tmpvals = prepared$tmpvals), mode = "auto"))

  if (is.call(result))
    return(.npRmpi_autodispatch_tag_result(.npRmpi_autodispatch_replace_tmps(result, tmpvals = prepared$tmpvals), mode = "auto"))

  .npRmpi_autodispatch_tag_result(result, mode = "auto")
}

.npRmpi_bootstrap_make_indices <- function(n,
                                           boot.num,
                                           boot.method = c("inid", "fixed", "geom"),
                                           boot.blocklen = NULL) {
  boot.method <- match.arg(boot.method)

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 2L)
    stop("invalid bootstrap sample size")
  if (!is.numeric(boot.num) || length(boot.num) != 1L || is.na(boot.num) || boot.num < 1L)
    stop("invalid bootstrap replication count")

  n <- as.integer(n)
  boot.num <- as.integer(boot.num)

  if (boot.method == "inid") {
    idx <- matrix(sample.int(n, size = n * boot.num, replace = TRUE),
                  nrow = boot.num, ncol = n)
    return(idx)
  }

  if (is.null(boot.blocklen) || !is.numeric(boot.blocklen) ||
      length(boot.blocklen) != 1L || is.na(boot.blocklen) || boot.blocklen < 1)
    stop("invalid block length for bootstrap method")

  boot.blocklen <- as.integer(boot.blocklen)
  idx.ts <- tsboot(tseries = seq_len(n),
                   statistic = function(tsb) tsb,
                   R = boot.num,
                   l = boot.blocklen,
                   sim = boot.method)
  idx <- idx.ts$t
  if (!is.matrix(idx))
    idx <- matrix(idx, nrow = boot.num, byrow = TRUE)
  storage.mode(idx) <- "integer"
  idx
}

.npRmpi_bootstrap_compute_payload <- function(payload, comm = 1L) {
  .npRmpi_require_active_slave_pool(comm = comm, where = "bootstrap payload computation")

  if (!is.list(payload))
    stop("bootstrap payload must be a list")

  tmp <- ".__npRmpi_boot_payload"
  .GlobalEnv[[tmp]] <- payload
  .npRmpi_bcast_robj_by_name(tmp, caller_env = .GlobalEnv)
  on.exit({
    rm(list = tmp, envir = .GlobalEnv)
    cmd.rm <- substitute(rm(list = TMP, envir = .GlobalEnv), list(TMP = tmp))
    .npRmpi_bcast_cmd_expr(cmd.rm, comm = comm, caller.execute = FALSE)
  }, add = TRUE)

  cmd <- substitute({
    old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
    old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
    options(npRmpi.autodispatch.context = TRUE)
    options(npRmpi.autodispatch.disable = TRUE)
    on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
    on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)

    p <- .__npRmpi_boot_payload
    B <- nrow(p$indices)
    if (is.null(p$family) || is.null(p$out.field))
      stop("invalid bootstrap payload: missing family/out.field")

    fit_one <- function(ix) {
      switch(
        p$family,
        npreg = suppressWarnings(npreg(
          txdat = p$xdat[ix, , drop = FALSE],
          tydat = p$ydat[ix],
          exdat = p$exdat,
          bws = p$bws,
          gradients = p$gradients,
          gradient.order = p$gradient.order,
          warn.glp.gradient = FALSE)),
        npcdens = suppressWarnings(npcdens(
          txdat = p$xdat[ix, , drop = FALSE],
          tydat = p$ydat[ix, , drop = FALSE],
          exdat = p$exdat,
          eydat = p$eydat,
          bws = p$bws,
          gradients = p$gradients)),
        npcdist = suppressWarnings(npcdist(
          txdat = p$xdat[ix, , drop = FALSE],
          tydat = p$ydat[ix, , drop = FALSE],
          exdat = p$exdat,
          eydat = p$eydat,
          bws = p$bws,
          gradients = p$gradients)),
        npqreg = suppressWarnings(npqreg(
          txdat = p$xdat[ix, , drop = FALSE],
          tydat = p$ydat[ix, , drop = FALSE],
          exdat = p$exdat,
          tau = p$tau,
          bws = p$bws,
          gradients = p$gradients)),
        npscoef = suppressWarnings(npscoef(
          txdat = p$xdat[ix, , drop = FALSE],
          tydat = p$ydat[ix],
          tzdat = if (is.null(p$zdat)) NULL else p$zdat[ix, , drop = FALSE],
          exdat = p$exdat,
          ezdat = p$ezdat,
          bws = p$bws)),
        npplreg = suppressWarnings(npplreg(
          txdat = p$xdat[ix, , drop = FALSE],
          tydat = p$ydat[ix],
          tzdat = p$zdat[ix, , drop = FALSE],
          exdat = p$exdat,
          ezdat = p$ezdat,
          bws = p$bws)),
        npindex = suppressWarnings(npindex(
          txdat = p$xdat[ix, , drop = FALSE],
          tydat = p$ydat[ix],
          exdat = p$exdat,
          bws = p$bws,
          gradients = p$gradients)),
        npudens = suppressWarnings(npudens(
          tdat = p$xdat[ix, , drop = FALSE],
          edat = p$exdat,
          bws = p$bws)),
        npudist = suppressWarnings(npudist(
          tdat = p$xdat[ix, , drop = FALSE],
          edat = p$exdat,
          bws = p$bws)),
        stop("unsupported bootstrap payload family"))
    }

    extract_stat <- function(obj) {
      out <- obj[[p$out.field]]
      if (!is.null(p$out.index))
        out <- out[, p$out.index]
      as.numeric(out)
    }

    base.fit <- fit_one(seq_len(nrow(p$xdat)))
    t0 <- extract_stat(base.fit)

    t.mat <- matrix(NA_real_, nrow = B, ncol = length(t0))
    for (b in seq_len(B)) {
      idx <- as.integer(p$indices[b, ])
      fit.b <- fit_one(idx)
      t.mat[b, ] <- extract_stat(fit.b)
    }

    list(t0 = as.numeric(t0), t = t.mat)
  })

  .npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE)
}
