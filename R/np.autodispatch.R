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

.npRmpi_bcast_cmd_expr <- function(expr, comm = 1L, caller.execute = TRUE) {
  eval(substitute(mpi.bcast.cmd(cmd = EXPR, comm = COMM, caller.execute = CE),
                  list(EXPR = expr, COMM = comm, CE = caller.execute)),
       envir = parent.frame())
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
  eval(mc, envir = caller_env)
}

.npRmpi_autodispatch_called_from_bcast <- function() {
  calls <- sys.calls()
  if (!length(calls)) return(FALSE)
  any(vapply(calls, function(cl) {
    if (!is.call(cl) || length(cl) < 1) return(FALSE)
    fn <- cl[[1]]
    is.symbol(fn) && identical(as.character(fn), "mpi.bcast.cmd")
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
      out[[1L]] <- as.call(list(as.name(":::"), as.name("npRmpi"), as.name(fname)))
    }
  }
  tmpnames <- character(0)
  idx <- 0L

  for (i in seq_along(arg.list)) {
    if (i == 1L) next
    nm <- nms[i]
    if (is.null(nm) || identical(nm, "")) next
    if (!nm %in% targets) next

    val <- .npRmpi_autodispatch_eval_arg(arg.list[[i]], caller_env = caller_env)
    idx <- idx + 1L
    tmp <- sprintf(".__npRmpi_autod_%s_%d", nm, idx)
    assign(tmp, val, envir = .GlobalEnv)
    cmd.assign <- substitute(assign(TMP, VAL, envir = .GlobalEnv),
                             list(TMP = tmp, VAL = val))
    .npRmpi_bcast_cmd_expr(cmd.assign, comm = comm, caller.execute = FALSE)

    out[[i]] <- as.name(tmp)
    tmpnames <- c(tmpnames, tmp)
  }

  list(call = out, tmpnames = unique(tmpnames))
}

.npRmpi_autodispatch_cleanup <- function(tmpnames, comm = 1L) {
  if (!length(tmpnames)) return(invisible(TRUE))

  rm(list = tmpnames, envir = .GlobalEnv)
  cmd.rm <- substitute(rm(list = TMPS, envir = .GlobalEnv),
                       list(TMPS = tmpnames))
  .npRmpi_bcast_cmd_expr(cmd.rm, comm = comm, caller.execute = FALSE)
  invisible(TRUE)
}

.npRmpi_autodispatch_call <- function(mc, caller_env = parent.frame(), comm = 1L) {
  .npRmpi_warn_pkg_conflict_once()
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

  if (!.npRmpi_autodispatch_preflight(comm = comm))
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  .npRmpi_autodispatch_sync_options(comm = comm)
  prepared <- .npRmpi_autodispatch_materialize_call(mc = mc, caller_env = caller_env, comm = comm)
  on.exit(.npRmpi_autodispatch_cleanup(prepared$tmpnames, comm = comm), add = TRUE)

  cmd <- substitute({
    old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
    old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
    options(npRmpi.autodispatch.context = TRUE)
    options(npRmpi.autodispatch.disable = TRUE)
    on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
    on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
    eval(CALL, envir = .GlobalEnv)
  }, list(CALL = prepared$call))

  .npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE)
}
