.npRmpi_autodispatch_option_keys <- function() {
  c("np.messages", "np.tree", "np.largeh.rel.tol",
    "np.disc.upper.rel.tol", "np.groupcv.fast")
}

.npRmpi_autodispatch_option_snapshot <- function(keys = .npRmpi_autodispatch_option_keys()) {
  setNames(lapply(keys, getOption), keys)
}

.npRmpi_autodispatch_should_sync_options <- function(snapshot) {
  mode <- getOption("npRmpi.autodispatch.option.sync", "onchange")
  if (!is.character(mode) || length(mode) != 1L || is.na(mode))
    mode <- "onchange"
  mode <- match.arg(mode, choices = c("onchange", "always", "never"))

  if (identical(mode, "never"))
    return(FALSE)

  if (identical(mode, "always")) {
    options(npRmpi.autodispatch.option.snapshot = snapshot)
    return(TRUE)
  }

  cached <- getOption("npRmpi.autodispatch.option.snapshot", NULL)
  if (!identical(cached, snapshot)) {
    options(npRmpi.autodispatch.option.snapshot = snapshot)
    return(TRUE)
  }

  FALSE
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

.npRmpi_abort_if_rmpi_attached <- function(where = "this call") {
  if (!("package:Rmpi" %in% search()))
    return(invisible(FALSE))
  stop(
    sprintf("%s cannot run because package 'Rmpi' is attached. Detach 'Rmpi' and use npRmpi APIs directly.", where),
    call. = FALSE
  )
}

.npRmpi_bcast_cmd_expr <- function(expr, comm = 1L, caller.execute = TRUE) {
  target <- get("mpi.bcast.cmd",
                envir = parent.frame(),
                mode = "function",
                inherits = TRUE)
  do.call(target,
          list(cmd = expr, comm = comm, caller.execute = caller.execute),
          envir = parent.frame())
}

.npRmpi_bcast_robj_by_name <- function(name, caller_env = parent.frame()) {
  fn <- get("mpi.bcast.Robj2slave",
            envir = asNamespace("npRmpi"),
            mode = "function",
            inherits = FALSE)
  call <- as.call(list(fn, as.name(name)))
  .npRmpi_eval_scmd(call, envir = caller_env)
}

.npRmpi_autodispatch_next_remote_name <- function() {
  id <- getOption("npRmpi.autodispatch.remote.counter", 0L)
  id <- as.integer(id) + 1L
  options(npRmpi.autodispatch.remote.counter = id)
  sprintf(".__npRmpi_autod_ret_%d", id)
}

.npRmpi_autodispatch_remote_ref <- function(x) {
  ref <- tryCatch(attr(x, "npRmpi.autodispatch.remote", exact = TRUE), error = function(e) NULL)
  if (!is.character(ref) || length(ref) != 1L || !nzchar(ref))
    return(NULL)
  ref
}

.npRmpi_autodispatch_large_arg_threshold <- function() {
  thr <- getOption("npRmpi.autodispatch.arg.broadcast.threshold", 4096L)
  if (!is.numeric(thr) || length(thr) != 1L || is.na(thr) || thr < 0)
    return(4096L)
  as.integer(thr)
}

.npRmpi_autodispatch_as_generic_call <- function(generic, mc) {
  args <- as.list(mc)[-1L]
  as.call(c(list(as.name(generic)), args))
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
    not_found <- new.env(parent = emptyenv())
    caller_fun <- get0(fname, envir = caller_env, mode = "function", inherits = TRUE, ifnotfound = not_found)
    ns_fun <- get0(fname, envir = asNamespace("npRmpi"), mode = "function", inherits = FALSE, ifnotfound = not_found)
    if (identical(caller_fun, not_found) && !identical(ns_fun, not_found))
      mc.eval[[1L]] <- ns_fun
  }
  .npRmpi_eval_scmd(mc.eval, envir = caller_env)
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
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
    return(FALSE)
  size <- tryCatch(mpi.comm.size(comm), error = function(e) NA_integer_)
  rank <- tryCatch(mpi.comm.rank(comm), error = function(e) NA_integer_)
  if (is.null(size) || is.null(rank) ||
      is.na(size) || is.na(rank))
    return(FALSE)
  size >= 2L
}

.npRmpi_master_only_mode <- function(comm = 1L) {
  if (!isTRUE(getOption("npRmpi.master.only", FALSE)))
    return(FALSE)
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
    return(FALSE)
  size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) NA_integer_)
  rank <- tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) NA_integer_)
  if (is.na(size) || is.na(rank) || size < 1L || rank < 0L) {
    size <- tryCatch(as.integer(mpi.comm.size(0)), error = function(e) NA_integer_)
    rank <- tryCatch(as.integer(mpi.comm.rank(0)), error = function(e) NA_integer_)
  }
  if (is.null(size) || is.null(rank) || is.na(size) || is.na(rank))
    return(FALSE)
  as.integer(size) == 1L && as.integer(rank) == 0L
}

.npRmpi_require_active_slave_pool <- function(comm = 1L,
                                              where = "this call",
                                              allow.master.only = TRUE) {
  .npRmpi_abort_if_rmpi_attached(where = where)
  if (.npRmpi_has_active_slave_pool(comm = comm))
    return(invisible(TRUE))
  if (isTRUE(allow.master.only) && .npRmpi_master_only_mode(comm = comm))
    return(invisible(TRUE))
  stop(sprintf("%s requires an active MPI slave pool; call npRmpi.init(...) first", where))
}

.npRmpi_autodispatch_preflight <- function(comm = 1L) {
  .npRmpi_abort_if_rmpi_attached(where = "npRmpi auto-dispatch")
  strict <- isTRUE(getOption("npRmpi.autodispatch.strict", TRUE))
  if (!.npRmpi_has_active_slave_pool(comm = comm)) {
    if (.npRmpi_master_only_mode(comm = comm)) {
      warning("npRmpi auto-dispatch is disabled in master-only mode (nslaves=0); executing call locally")
      return(FALSE)
    }
    msg <- "npRmpi auto-dispatch requires an active slave pool; call npRmpi.init(...) first"
    if (strict) stop(msg)
    warning(msg)
    return(FALSE)
  }

  TRUE
}

.npRmpi_autodispatch_target_args <- function() {
  c("formula", "data", "bws",
     "dat", "tdat", "edat",
     "xdat", "ydat", "zdat", "txdat", "tydat", "tzdat",
     "exdat", "eydat", "ezdat", "newdata",
     "gydat", "gdat", "wdat", "gdata",
     "data.x", "data.y", "model",
     "y", "z", "w", "x", "zeval", "weval", "xeval", "bw",
     "regtype", "basis", "degree", "bernstein.basis",
     "nmulti", "bandwidth.compute", "remin", "itmax", "ftol", "tol", "small",
     "gradients", "residuals", "errors", "gradient.order")
}

.npRmpi_autodispatch_failfast_formula_data <- function(mc, caller_env) {
  invisible(FALSE)
}

.npRmpi_autodispatch_eval_arg <- function(expr, caller_env) {
  val <- tryCatch(.npRmpi_eval_scmd(expr, envir = caller_env), error = function(e) e)
  if (!inherits(val, "error"))
    return(val)

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    env_i <- frames[[i]]
    if (identical(env_i, caller_env))
      next
    val_i <- tryCatch(.npRmpi_eval_scmd(expr, envir = env_i), error = function(e) e)
    if (!inherits(val_i, "error"))
      return(val_i)
  }

  stop(conditionMessage(val), call. = FALSE)
}

.npRmpi_autodispatch_materialize_call <- function(mc, caller_env, comm = 1L) {
  arg.list <- as.list(mc)
  nms <- names(arg.list)
  targets <- .npRmpi_autodispatch_target_args()

  out <- mc
  if (is.call(out) && length(out) >= 1L && is.symbol(out[[1L]])) {
    fname <- as.character(out[[1L]])
    if (grepl(".", fname, fixed = TRUE)) {
      # For non-exported S3 methods (e.g., npcdens.conbandwidth), broadcast
      # the generic symbol to avoid serializing closure heads in MPI commands.
      generic <- strsplit(fname, ".", fixed = TRUE)[[1L]][1L]
      generic_fun <- get0(generic, envir = asNamespace("npRmpi"), mode = "function", inherits = FALSE)
      if (nzchar(generic) && is.function(generic_fun)) {
        out[[1L]] <- as.name(generic)
      }
    }
  }
  tmpnames <- character(0)
  tmpvals <- list()
  prepublish <- list()
  idx <- 0L

  has_data_inputs <- !is.null(nms) && any(nms %in% c("data", "xdat", "ydat", "txdat", "tydat", "zdat"))

  formula.expr <- NULL
  if (!is.null(nms) && any(nms == "formula")) {
    formula.expr <- arg.list[[which(nms == "formula")[1L]]]
  } else if (!is.null(nms) && any(nms == "bws")) {
    bexpr <- arg.list[[which(nms == "bws")[1L]]]
    bval <- tryCatch(.npRmpi_autodispatch_eval_arg(bexpr, caller_env = caller_env),
                     error = function(e) NULL)
    if (!is.null(bval) && inherits(bval, "formula"))
      formula.expr <- bexpr
  } else if (length(arg.list) >= 2L) {
    nm2 <- if (!is.null(nms) && length(nms) >= 2L) nms[[2L]] else ""
    if (is.null(nm2) || identical(nm2, "")) {
      fval <- tryCatch(.npRmpi_autodispatch_eval_arg(arg.list[[2L]], caller_env = caller_env),
                       error = function(e) NULL)
      if (!is.null(fval) && inherits(fval, "formula"))
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

    # Prefer forcing already-bound formal arguments in the caller frame.
    # This avoids unresolved placeholders such as `..1` from match.call()
    # when methods forward `...` through nested generic/method dispatch.
    not_found <- new.env(parent = emptyenv())
    val <- tryCatch(get0(nm, envir = caller_env, inherits = FALSE, ifnotfound = not_found),
                    error = function(e) not_found)
    if (identical(val, not_found))
      val <- .npRmpi_autodispatch_eval_arg(arg.list[[i]], caller_env = caller_env)
    ref <- .npRmpi_autodispatch_remote_ref(val)
    # `bws` objects can be post-processed locally (e.g. formula methods
    # rewriting call/formula metadata). Reusing a stale remote reference can
    # reintroduce unresolved temporary symbols during downstream plotting.
    if (!is.null(ref) && !identical(nm, "bws")) {
      out[[i]] <- as.name(ref)
      next
    }
    idx <- idx + 1L
    tmp <- sprintf(".__npRmpi_autod_%s_%d", nm, idx)

    out[[i]] <- as.name(tmp)
    # Preserve formula-method dispatch semantics for generics that use
    # formal name `bws` in the generic but `formula` in the method.
    if (identical(nm, "bws") && inherits(val, "formula")) {
      out.nms <- names(out)
      if (!is.null(out.nms) && length(out.nms) >= i)
        out.nms[i] <- ""
      names(out) <- out.nms
    }
    tmpnames <- c(tmpnames, tmp)
    if (as.numeric(object.size(val)) >= .npRmpi_autodispatch_large_arg_threshold()) {
      prepublish[[tmp]] <- val
    } else {
      tmpvals[[tmp]] <- val
    }
  }

  list(call = out, tmpnames = unique(tmpnames), tmpvals = tmpvals, prepublish = prepublish)
}

.npRmpi_autodispatch_cleanup <- function(tmpnames, comm = 1L) {
  if (!length(tmpnames))
    return(invisible(TRUE))
  .npRmpi_rm_existing(tmpnames, envir = .GlobalEnv)
  cmd.rm <- substitute(get(".npRmpi_rm_existing", envir = asNamespace("npRmpi"), inherits = FALSE)(TMPS, envir = .GlobalEnv),
                       list(TMPS = tmpnames))
  .npRmpi_bcast_cmd_expr(cmd.rm, comm = comm, caller.execute = FALSE)
  invisible(TRUE)
}

.npRmpi_autodispatch_call_name <- function(mc) {
  clean <- function(x) sub("\\..*$", "", x)
  if (!is.call(mc) || length(mc) < 1L)
    return("autodispatch-call")
  hd <- mc[[1L]]
  if (is.symbol(hd))
    return(clean(as.character(hd)))
  if (is.character(hd) && length(hd))
    return(clean(as.character(hd)[1L]))
  if (is.call(hd) && length(hd) >= 3L &&
      is.symbol(hd[[1L]]) && as.character(hd[[1L]]) %in% c("::", ":::")) {
    lhs <- hd[[2L]]
    rhs <- hd[[3L]]
    if (is.symbol(lhs) && is.symbol(rhs))
      return(sprintf("%s%s%s",
                     as.character(lhs),
                     as.character(hd[[1L]]),
                     as.character(rhs)))
  }
  "autodispatch-call"
}

.npRmpi_autodispatch_attach_timing <- function(obj, rec) {
  if (!is.list(obj) || !is.list(rec))
    return(obj)

  if (is.null(obj$timing.profile))
    obj$timing.profile <- rec
  if (is.list(obj$bws) && is.null(obj$bws$timing.profile))
    obj$bws$timing.profile <- rec
  obj
}

.npRmpi_is_missing_call_arg <- function(arg) {
  if (missing(arg))
    return(TRUE)
  is.symbol(arg) && identical(as.character(arg), "")
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
    xlist <- as.list(x)
    for (i in seq_along(xlist)) {
      xi <- xlist[[i]]
      if (.npRmpi_is_missing_call_arg(xi))
        next
      xlist[[i]] <- .npRmpi_autodispatch_replace_tmps(xi, tmpvals = tmpvals)
    }
    return(as.call(xlist))
  }

  if (is.pairlist(x)) {
    xlist <- as.list(x)
    for (i in seq_along(xlist)) {
      xi <- xlist[[i]]
      if (.npRmpi_is_missing_call_arg(xi))
        next
      xlist[[i]] <- .npRmpi_autodispatch_replace_tmps(xi, tmpvals = tmpvals)
    }
    return(as.pairlist(xlist))
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

.npRmpi_rm_existing <- function(nms, envir = .GlobalEnv) {
  if (!length(nms))
    return(invisible(character(0)))
  nms <- unique(as.character(nms))
  present <- nms[vapply(nms, exists, logical(1), envir = envir, inherits = FALSE)]
  if (length(present))
    rm(list = present, envir = envir)
  invisible(present)
}

.npRmpi_autodispatch_tag_result <- function(x, mode = "auto", remote = NULL) {
  if (!is.null(x)) {
    attr(x, "npRmpi.dispatch.mode") <- mode
    if (is.character(remote) && length(remote) == 1L && nzchar(remote))
      attr(x, "npRmpi.autodispatch.remote") <- remote
  }
  x
}

.npRmpi_autodispatch_untag <- function(x) {
  if (is.list(x) && !is.pairlist(x) && !is.call(x)) {
    for (i in seq_len(length(x))) {
      xi <- tryCatch(x[[i]], error = function(e) NULL)
      if (!is.null(xi))
        x[i] <- list(.npRmpi_autodispatch_untag(xi))
    }
  }
  attr(x, "npRmpi.dispatch.mode") <- NULL
  attr(x, "npRmpi.autodispatch.remote") <- NULL
  x
}

.npRmpi_guard_no_auto_object_in_manual_bcast <- function(obj, where = "this call") {
  if (.npRmpi_autodispatch_in_context())
    return(invisible(FALSE))
  if (!.npRmpi_autodispatch_called_from_bcast())
    return(invisible(FALSE))
  mode <- tryCatch(attr(obj, "npRmpi.dispatch.mode", exact = TRUE), error = function(e) NULL)
  if (is.null(mode))
    return(invisible(FALSE))
  if (identical(mode, "auto")) {
    stop(sprintf("%s received an object created under npRmpi.autodispatch and cannot be executed inside mpi.bcast.cmd(...): avoid mixing dispatch modes; either rerun the full workflow in manual broadcast mode or keep this call outside mpi.bcast.cmd", where))
  }
  invisible(FALSE)
}

.npRmpi_distributed_call_impl <- function(mc,
                                          caller_env = parent.frame(),
                                          comm = 1L,
                                          warn_nested = FALSE) {
  t.start <- proc.time()
  start.wall <- Sys.time()
  comm.elapsed <- 0.0
  comm.calls <- 0L
  note.vec <- character(0)
  rec.comm <- function(expr, note = NA_character_) {
    t0 <- proc.time()[["elapsed"]]
    out <- force(expr)
    dt <- proc.time()[["elapsed"]] - t0
    if (!is.finite(dt) || dt < 0)
      dt <- 0.0
    # Intentional <<- : closure-level accumulator for communication timing.
    comm.elapsed <<- comm.elapsed + as.double(dt)
    # Intentional <<- : closure-level counter for communication events.
    comm.calls <<- as.integer(comm.calls) + 1L
    if (!is.na(note) && nzchar(as.character(note)[1L]) &&
        identical(getOption("npRmpi.profile.level", "basic"), "detailed")) {
      # Intentional <<- : append detailed event notes in parent accumulator.
      note.vec <<- c(note.vec, as.character(note)[1L])
    }
    out
  }
  rec.step <- function(expr, note = NA_character_) {
    t0 <- proc.time()[["elapsed"]]
    out <- force(expr)
    dt <- proc.time()[["elapsed"]] - t0
    if (!is.finite(dt) || dt < 0)
      dt <- 0.0
    if (!is.na(note) && nzchar(as.character(note)[1L]) &&
        identical(getOption("npRmpi.profile.level", "basic"), "detailed")) {
      # Intentional <<- : append detailed step notes in parent accumulator.
      note.vec <<- c(note.vec, as.character(note)[1L])
    }
    out
  }
  make.rec <- function() {
    dt <- proc.time() - t.start
    wall <- unname(as.double(dt[["elapsed"]]))
    if (!is.finite(wall) || wall < 0)
      wall <- 0.0
    comm <- as.double(comm.elapsed)
    if (!is.finite(comm) || comm < 0)
      comm <- 0.0
    compute <- max(0.0, wall - comm)
    denom <- comm + compute
    ratio <- if (denom > 0) min(1.0, max(0.0, comm / denom)) else NA_real_
    rec <- list(
      profile_kind = "call",
      where = .npRmpi_autodispatch_call_name(mc),
      method = NA_character_,
      B = NA_integer_,
      ntrain = NA_integer_,
      neval = NA_integer_,
      rank = tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) NA_integer_),
      size = tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) NA_integer_),
      via_bcast = isTRUE(.npRmpi_autodispatch_called_from_bcast()),
      wall_elapsed_sec = wall,
      comm_elapsed_sec = comm,
      compute_elapsed_sec = compute,
      comm_ratio = ratio,
      comm_calls = as.integer(comm.calls),
      user_sec = unname(as.double(dt[["user.self"]])),
      system_sec = unname(as.double(dt[["sys.self"]])),
      timestamp_start = start.wall,
      timestamp_end = Sys.time()
    )
    if (identical(getOption("npRmpi.profile.level", "basic"), "detailed"))
      rec$comm_notes <- note.vec
    rec
  }

  if (.npRmpi_autodispatch_called_from_bcast()) {
    if (warn_nested)
      .npRmpi_autodispatch_warn_nested()
    return(.npRmpi_eval_without_dispatch(mc, caller_env))
  }

  if (.npRmpi_autodispatch_in_context())
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  rank <- tryCatch(mpi.comm.rank(comm), error = function(e) NA_integer_)
  if (!is.na(rank) && rank != 0L)
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  .npRmpi_autodispatch_failfast_formula_data(mc, caller_env = caller_env)

  if (!.npRmpi_autodispatch_preflight(comm = comm))
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  prepared <- .npRmpi_autodispatch_materialize_call(mc = mc, caller_env = caller_env, comm = comm)
  cleaned <- FALSE
  if (length(prepared$tmpnames))
    on.exit({
      if (!isTRUE(cleaned))
        rec.comm(.npRmpi_autodispatch_cleanup(prepared$tmpnames, comm = comm),
                 note = "autodispatch.cleanup")
    }, add = TRUE)
  if (length(prepared$tmpnames))
    for (nm in prepared$tmpnames)
      if (!is.null(prepared$tmpvals[[nm]]))
        prepared$tmpvals[[nm]] <- .npRmpi_autodispatch_untag(prepared$tmpvals[[nm]])
      else if (!is.null(prepared$prepublish[[nm]]))
        prepared$prepublish[[nm]] <- .npRmpi_autodispatch_untag(prepared$prepublish[[nm]])

  if (length(prepared$prepublish)) {
    for (nm in names(prepared$prepublish)) {
      .GlobalEnv[[nm]] <- prepared$prepublish[[nm]]
      rec.comm(.npRmpi_bcast_robj_by_name(nm, caller_env = .GlobalEnv),
               note = "mpi.bcast.Robj2slave")
    }
  }

  opt.keys <- .npRmpi_autodispatch_option_keys()
  opt.snapshot <- .npRmpi_autodispatch_option_snapshot(opt.keys)
  opt.sync <- .npRmpi_autodispatch_should_sync_options(opt.snapshot)
  if (opt.sync) {
    opt.keys <- names(opt.snapshot)
    opt.vals <- unname(opt.snapshot)
  } else {
    opt.keys <- character(0)
    opt.vals <- list()
  }
  opt.verify <- opt.sync && isTRUE(getOption("npRmpi.autodispatch.verify.options", FALSE))
  remote.name <- .npRmpi_autodispatch_next_remote_name()

  cmd <- substitute({
    for (i in seq_along(OPT_KEYS))
      options(structure(list(OPT_VALS[[i]]), names = OPT_KEYS[[i]]))

    if (OPT_VERIFY) {
      for (i in seq_along(OPT_KEYS)) {
        lval <- getOption(OPT_KEYS[[i]])
        if (!identical(lval, OPT_VALS[[i]]))
          stop(sprintf("failed to synchronize option '%s' across MPI ranks", OPT_KEYS[[i]]))
      }
    }

    tmpvals <- TMPVALS
    if (length(tmpvals)) {
      for (nm in names(tmpvals))
        .GlobalEnv[[nm]] <- tmpvals[[nm]]
    }

    old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
    old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
    options(npRmpi.autodispatch.context = TRUE)
    options(npRmpi.autodispatch.disable = TRUE)
    on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
    on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
    if (length(TMP_NAMES))
      on.exit(get(".npRmpi_rm_existing", envir = asNamespace("npRmpi"), inherits = FALSE)(TMP_NAMES, envir = .GlobalEnv), add = TRUE)

    res <- CALL
    .GlobalEnv[[REMOTE_NAME]] <- res
    res
  }, list(CALL = prepared$call,
          TMPVALS = prepared$tmpvals,
          TMP_NAMES = prepared$tmpnames,
          OPT_KEYS = opt.keys,
          OPT_VALS = opt.vals,
          OPT_VERIFY = opt.verify,
          REMOTE_NAME = remote.name))

  result <- rec.step(.npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE),
                     note = "mpi.bcast.cmd.execute")
  if (length(prepared$tmpnames)) {
    rec.comm(.npRmpi_autodispatch_cleanup(prepared$tmpnames, comm = comm),
             note = "autodispatch.cleanup")
    cleaned <- TRUE
  }
  tmpreplace <- c(prepared$tmpvals, prepared$prepublish)
  rec <- make.rec()

  if (is.list(result))
    return(.npRmpi_autodispatch_attach_timing(
      .npRmpi_autodispatch_tag_result(
        .npRmpi_autodispatch_replace_tmp_calls(result, tmpvals = tmpreplace),
        mode = "auto", remote = remote.name
      ),
      rec
    ))

  if (is.call(result))
    return(.npRmpi_autodispatch_tag_result(
      .npRmpi_autodispatch_replace_tmps(result, tmpvals = tmpreplace),
      mode = "auto", remote = remote.name
    ))

  .npRmpi_autodispatch_tag_result(result, mode = "auto", remote = remote.name)
}

.npRmpi_autodispatch_call <- function(mc, caller_env = parent.frame(), comm = 1L) {
  .npRmpi_warn_pkg_conflict_once()
  .npRmpi_warn_rmpi_conflict_once()
  caller_env <- parent.frame()
  if (!.npRmpi_autodispatch_active())
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  .npRmpi_distributed_call_impl(mc = mc, caller_env = caller_env, comm = comm, warn_nested = TRUE)
}

.npRmpi_manual_distributed_call <- function(mc, caller_env = parent.frame(), comm = 1L) {
  .npRmpi_warn_pkg_conflict_once()
  .npRmpi_warn_rmpi_conflict_once()
  caller_env <- parent.frame()
  .npRmpi_distributed_call_impl(mc = mc, caller_env = caller_env, comm = comm, warn_nested = FALSE)
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
  .npRmpi_require_active_slave_pool(
    comm = comm,
    where = "bootstrap payload computation",
    allow.master.only = FALSE
  )

  if (!is.list(payload))
    stop("bootstrap payload must be a list")

  tmp <- ".__npRmpi_boot_payload"
  .GlobalEnv[[tmp]] <- payload
  .npRmpi_bcast_robj_by_name(tmp, caller_env = .GlobalEnv)
  on.exit({
    .npRmpi_rm_existing(tmp, envir = .GlobalEnv)
    cmd.rm <- substitute(get(".npRmpi_rm_existing", envir = asNamespace("npRmpi"), inherits = FALSE)(TMP, envir = .GlobalEnv), list(TMP = tmp))
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
