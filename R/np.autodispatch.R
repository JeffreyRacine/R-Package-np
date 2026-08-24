.npRmpi_autodispatch_option_keys <- function() {
  c("np.messages", "np.tree", "np.categorical.compress", "np.largeh",
    "np.objective.cache", "np.largelambda", "np.extendednn", "np.largeh.rel.tol",
    "np.disc.upper.rel.tol", "np.nomad.degree.start.policy")
}

.npRmpi_autodispatch_option_snapshot <- function(keys = .npRmpi_autodispatch_option_keys()) {
  setNames(lapply(keys, getOption), keys)
}

.npRmpi_autodispatch_validate_option_snapshot <- function(snapshot) {
  if (is.list(snapshot) && "np.objective.cache" %in% names(snapshot))
    npObjectiveCacheEnabled(snapshot[["np.objective.cache"]])
  invisible(TRUE)
}

.npRmpi_autodispatch_should_sync_options <- function(snapshot) {
  .npRmpi_autodispatch_validate_option_snapshot(snapshot)
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
  .np_warning("both packages 'npRmpi' and 'np' are attached: use explicit npRmpi:: calls to avoid masked-function ambiguity")
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
  .np_warning("both packages 'npRmpi' and 'Rmpi' are attached: prefer explicit npRmpi:: calls in user code to avoid dispatch ambiguity")
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

.npRmpi_lease_process_state <- new.env(parent = emptyenv())
.npRmpi_lease_process_state$counter <- 0

.npRmpi_lease_state <- new.env(parent = emptyenv())
.npRmpi_lease_state$session <- NULL
.npRmpi_lease_state$status <- "open"
.npRmpi_lease_state$store <- new.env(parent = emptyenv())
.npRmpi_lease_state$registry <- new.env(parent = emptyenv())
.npRmpi_lease_state$pending <- new.env(parent = emptyenv())

.npRmpi_lease_next_id <- function(prefix) {
  .npRmpi_lease_process_state$counter <- .npRmpi_lease_process_state$counter + 1
  sprintf("%s_p%d_c%016.0f", prefix, Sys.getpid(), .npRmpi_lease_process_state$counter)
}

.npRmpi_lease_capability_table <- function() {
  data.frame(
    producer = c("npregbw", "nplsqregbw", "npplregbw", "npindexbw", "npscoefbw",
                 ".npscoefbw_nomad_context_prepare"),
    kind = c(rep("leased_public", 5L), "scoped_internal"),
    capability = c("bw.npreg", "bw.nplsqreg", "bw.npplreg", "bw.npindex", "bw.npscoef",
                   "internal.npscoef.nomad_context"),
    consumer = c("npreg", "nplsqreg", "npplreg", "npindex", "npscoef",
                 ".npscoefbw_nomad_pool_start"),
    argument = c(rep("bws", 5L), "ctx"),
    expected_class = c("rbandwidth", "lsqregressionbandwidth", "plbandwidth",
                       "sibandwidth", "scbandwidth", "list"),
    fingerprint_policy = c(rep("bandwidth.v1", 5L), "context.v1"),
    stringsAsFactors = FALSE
  )
}

.npRmpi_lease_context_symbol <- function() {
  as.name("npRmpi_internal_npscoefbw_nomad_context_prepare")
}

.npRmpi_lease_is_context_call <- function(mc) {
  if (!is.call(mc) || length(mc) < 1L)
    return(FALSE)
  head <- mc[[1L]]
  if (is.symbol(head))
    return(identical(head, .npRmpi_lease_context_symbol()))
  if (!is.function(head))
    return(FALSE)
  target <- get0(".npscoefbw_nomad_context_prepare",
                 envir = asNamespace("npRmpi"), mode = "function", inherits = FALSE)
  is.function(target) && identical(head, target)
}

.npRmpi_lease_producer <- function(mc) {
  if (.npRmpi_lease_is_context_call(mc))
    return(".npscoefbw_nomad_context_prepare")
  sub("\\..*$", "", .npRmpi_autodispatch_call_name(mc))
}

.npRmpi_lease_require_open <- function() {
  if (!identical(.npRmpi_lease_state$status, "open"))
    stop(sprintf("npRmpi autodispatch lease session is not reusable: state=%s",
                 .npRmpi_lease_state$status), call. = FALSE)
  invisible(TRUE)
}

.npRmpi_lease_ensure_session <- function() {
  .npRmpi_lease_require_open()
  if (is.null(.npRmpi_lease_state$session))
    .npRmpi_lease_state$session <- .npRmpi_lease_next_id("session")
  .npRmpi_lease_state$session
}

.npRmpi_lease_remove_local <- function(ids) {
  ids <- unique(as.character(ids))
  ids <- ids[nzchar(ids)]
  for (id in ids) {
    meta <- if (exists(id, envir = .npRmpi_lease_state$registry, inherits = FALSE))
      get(id, envir = .npRmpi_lease_state$registry, inherits = FALSE) else NULL
    if (is.list(meta) && identical(meta$kind, "scoped_internal") &&
        is.character(meta$remote_name) && length(meta$remote_name) == 1L &&
        exists(meta$remote_name, envir = .GlobalEnv, inherits = FALSE))
      rm(list = meta$remote_name, envir = .GlobalEnv)
    if (exists(id, envir = .npRmpi_lease_state$store, inherits = FALSE))
      rm(list = id, envir = .npRmpi_lease_state$store)
    if (exists(id, envir = .npRmpi_lease_state$registry, inherits = FALSE))
      rm(list = id, envir = .npRmpi_lease_state$registry)
  }
  invisible(ids)
}

.npRmpi_lease_reset_local <- function() {
  .npRmpi_lease_remove_local(ls(.npRmpi_lease_state$registry, all.names = TRUE))
  pending <- ls(.npRmpi_lease_state$pending, all.names = TRUE)
  if (length(pending))
    rm(list = pending, envir = .npRmpi_lease_state$pending)
  .npRmpi_lease_state$session <- NULL
  .npRmpi_lease_state$status <- "open"
  invisible(TRUE)
}

.npRmpi_lease_poison <- function() {
  .npRmpi_lease_state$status <- "poisoned"
  invisible(FALSE)
}

.npRmpi_lease_publication_plan <- function(mc) {
  producer <- .npRmpi_lease_producer(mc)
  policy <- .npRmpi_lease_capability_table()
  hit <- which(policy$producer == producer)
  if (length(hit) != 1L)
    return(list(kind = "none", producer = producer))
  row <- policy[hit, , drop = FALSE]
  session <- .npRmpi_lease_ensure_session()
  id <- paste(session, .npRmpi_lease_next_id("lease"), sep = "__")
  list(
    kind = row$kind[[1L]],
    producer = row$producer[[1L]],
    capability = row$capability[[1L]],
    expected_class = row$expected_class[[1L]],
    fingerprint_policy = row$fingerprint_policy[[1L]],
    session = session,
    id = id,
    remote_name = if (identical(row$kind[[1L]], "scoped_internal"))
      paste0(".__npRmpi_scoped_", id) else id
  )
}

.npRmpi_lease_validate_plan <- function(plan) {
  if (!is.list(plan) || !is.character(plan$kind) || length(plan$kind) != 1L)
    stop("invalid autodispatch publication plan", call. = FALSE)
  if (identical(plan$kind, "none"))
    return(invisible(TRUE))
  tab <- .npRmpi_lease_capability_table()
  hit <- which(tab$producer == plan$producer & tab$kind == plan$kind &
               tab$capability == plan$capability &
               tab$expected_class == plan$expected_class &
               tab$fingerprint_policy == plan$fingerprint_policy)
  if (length(hit) != 1L)
    stop("autodispatch publication plan is not capability-authorized", call. = FALSE)
  req <- c("session", "id", "remote_name")
  if (!all(vapply(req, function(nm) {
    val <- plan[[nm]]
    is.character(val) && length(val) == 1L && !is.na(val) && nzchar(val)
  }, logical(1))))
    stop("autodispatch publication plan has invalid identity", call. = FALSE)
  invisible(TRUE)
}

.npRmpi_lease_adopt_session <- function(session) {
  .npRmpi_lease_require_open()
  if (is.null(.npRmpi_lease_state$session)) {
    .npRmpi_lease_state$session <- session
  } else if (!identical(.npRmpi_lease_state$session, session)) {
    stop("autodispatch publication session-generation mismatch", call. = FALSE)
  }
  invisible(TRUE)
}

.npRmpi_lease_prepare_local <- function(value, plan) {
  .npRmpi_lease_validate_plan(plan)
  if (identical(plan$kind, "none"))
    return(list(kind = "none", id = "", capability = "", fingerprint = ""))
  .npRmpi_lease_adopt_session(plan$session)
  class.ok <- if (identical(plan$expected_class, "list")) is.list(value) else inherits(value, plan$expected_class)
  if (!isTRUE(class.ok))
    stop(sprintf("autodispatch capability '%s' returned unexpected class", plan$capability), call. = FALSE)
  if (exists(plan$id, envir = .npRmpi_lease_state$registry, inherits = FALSE) ||
      exists(plan$id, envir = .npRmpi_lease_state$store, inherits = FALSE) ||
      exists(plan$remote_name, envir = .GlobalEnv, inherits = FALSE))
    stop("autodispatch publication identity collision", call. = FALSE)
  fingerprint <- .npRmpi_autodispatch_fingerprint(value, policy = plan$fingerprint_policy)
  if (!is.character(fingerprint) || length(fingerprint) != 1L || is.na(fingerprint))
    stop("autodispatch canonical fingerprint failed", call. = FALSE)
  assign(plan$id, value, envir = .npRmpi_lease_state$store)
  assign(plan$id, list(
    id = plan$id, session = plan$session, kind = plan$kind,
    capability = plan$capability, fingerprint = fingerprint,
    remote_name = plan$remote_name, state = "provisional",
    estimated_bytes = as.numeric(object.size(value))
  ), envir = .npRmpi_lease_state$registry)
  list(kind = plan$kind, id = plan$id, capability = plan$capability,
       fingerprint = fingerprint)
}

.npRmpi_lease_commit_local <- function(plan) {
  .npRmpi_lease_validate_plan(plan)
  .npRmpi_lease_adopt_session(plan$session)
  if (!exists(plan$id, envir = .npRmpi_lease_state$registry, inherits = FALSE) ||
      !exists(plan$id, envir = .npRmpi_lease_state$store, inherits = FALSE))
    stop("autodispatch commit is missing provisional state", call. = FALSE)
  meta <- get(plan$id, envir = .npRmpi_lease_state$registry, inherits = FALSE)
  if (!identical(meta$state, "provisional") || !identical(meta$capability, plan$capability))
    stop("autodispatch commit found contradictory provisional state", call. = FALSE)
  if (identical(plan$kind, "scoped_internal")) {
    if (exists(plan$remote_name, envir = .GlobalEnv, inherits = FALSE))
      stop("autodispatch scoped publication would overwrite an existing binding", call. = FALSE)
    assign(plan$remote_name,
           get(plan$id, envir = .npRmpi_lease_state$store, inherits = FALSE),
           envir = .GlobalEnv)
    rm(list = plan$id, envir = .npRmpi_lease_state$store)
  }
  meta$state <- "live"
  assign(plan$id, meta, envir = .npRmpi_lease_state$registry)
  invisible(TRUE)
}

.npRmpi_lease_rollback_local <- function(ids) {
  .npRmpi_lease_remove_local(ids)
  invisible(TRUE)
}

.npRmpi_lease_session_drain_local <- function() {
  .npRmpi_lease_remove_local(ls(.npRmpi_lease_state$registry, all.names = TRUE))
  pending <- ls(.npRmpi_lease_state$pending, all.names = TRUE)
  if (length(pending))
    rm(list = pending, envir = .npRmpi_lease_state$pending)
  .npRmpi_lease_state$session <- NULL
  invisible(TRUE)
}

.npRmpi_lease_finalizer_factory <- function(id, pending) {
  force(id)
  force(pending)
  function(e) {
    if (!exists(id, envir = pending, inherits = FALSE))
      assign(id, TRUE, envir = pending)
    invisible(NULL)
  }
}

.npRmpi_lease_attach_handle <- function(value, plan, fingerprint) {
  handle <- new.env(parent = emptyenv())
  handle$id <- plan$id
  handle$session <- plan$session
  handle$capability <- plan$capability
  handle$fingerprint <- fingerprint
  reg.finalizer(
    handle,
    .npRmpi_lease_finalizer_factory(plan$id, .npRmpi_lease_state$pending),
    onexit = FALSE
  )
  attr(value, "npRmpi.autodispatch.lease") <- handle
  value
}

.npRmpi_lease_handle_fields <- function(value) {
  handle <- tryCatch(attr(value, "npRmpi.autodispatch.lease", exact = TRUE),
                     error = function(e) NULL)
  if (!is.environment(handle))
    return(NULL)
  needed <- c("id", "session", "capability", "fingerprint")
  if (!all(vapply(needed, exists, logical(1), envir = handle, inherits = FALSE)))
    return(NULL)
  setNames(lapply(needed, get, envir = handle, inherits = FALSE), needed)
}

.npRmpi_lease_resolve <- function(value, consumer, argument) {
  if (!identical(.npRmpi_lease_state$status, "open"))
    return(NULL)
  fields <- .npRmpi_lease_handle_fields(value)
  if (is.null(fields) || is.null(.npRmpi_lease_state$session) ||
      !identical(fields$session, .npRmpi_lease_state$session))
    return(NULL)
  tab <- .npRmpi_lease_capability_table()
  hit <- which(tab$kind == "leased_public" & tab$capability == fields$capability &
               tab$consumer == consumer & tab$argument == argument)
  if (length(hit) != 1L || !inherits(value, tab$expected_class[[hit]]))
    return(NULL)
  if ((!is.null(value$call) && .npRmpi_autodispatch_has_tmp_symbols(value$call)) ||
      (!is.null(value$formula) && .npRmpi_autodispatch_has_tmp_symbols(value$formula)))
    return(NULL)
  fingerprint <- .npRmpi_autodispatch_fingerprint(value, policy = tab$fingerprint_policy[[hit]])
  if (!identical(fingerprint, fields$fingerprint) ||
      !exists(fields$id, envir = .npRmpi_lease_state$registry, inherits = FALSE) ||
      !exists(fields$id, envir = .npRmpi_lease_state$store, inherits = FALSE))
    return(NULL)
  meta <- get(fields$id, envir = .npRmpi_lease_state$registry, inherits = FALSE)
  if (!identical(meta$state, "live") || !identical(meta$session, fields$session) ||
      !identical(meta$capability, fields$capability) ||
      !identical(meta$fingerprint, fields$fingerprint))
    return(NULL)
  list(id = fields$id, capability = fields$capability,
       fingerprint = fields$fingerprint)
}

.npRmpi_lease_diagnostics <- function() {
  ids <- ls(.npRmpi_lease_state$registry, all.names = TRUE)
  rows <- lapply(ids, function(id) {
    meta <- get(id, envir = .npRmpi_lease_state$registry, inherits = FALSE)
    data.frame(id = id, capability = meta$capability, kind = meta$kind,
               state = meta$state, estimated_bytes = meta$estimated_bytes,
               stringsAsFactors = FALSE)
  })
  list(
    session = if (is.null(.npRmpi_lease_state$session)) NA_character_ else .npRmpi_lease_state$session,
    state = .npRmpi_lease_state$status,
    pending_count = length(ls(.npRmpi_lease_state$pending, all.names = TRUE)),
    entries = if (length(rows)) do.call(rbind, rows) else
      data.frame(id = character(), capability = character(), kind = character(),
                 state = character(), estimated_bytes = numeric(), stringsAsFactors = FALSE)
  )
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

.npRmpi_autodispatch_raw_md5 <- function(raw) {
  path <- tempfile("npRmpi-autodispatch-fingerprint-", fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  writeBin(raw, path, useBytes = TRUE)
  unname(as.character(tools::md5sum(path)))
}

.npRmpi_autodispatch_is_degree_trace <- function(x) {
  if (!is.data.frame(x) || is.null(names(x)) || anyDuplicated(names(x)))
    return(FALSE)
  objective <- intersect(c("objective", "fval"), names(x))
  required <- c(
    "trace_id", "eval_id", "degree", "status", "cached", "message",
    "elapsed", "num.feval", "num.feval.fast"
  )
  if (length(objective) != 1L || !all(required %in% names(x)))
    return(FALSE)
  is.numeric(x$trace_id) &&
    is.numeric(x$eval_id) &&
    is.character(x$degree) &&
    is.numeric(x[[objective]]) &&
    is.character(x$status) &&
    is.logical(x$cached) &&
    is.character(x$message) &&
    is.numeric(x$elapsed) &&
    is.numeric(x$num.feval) &&
    is.numeric(x$num.feval.fast)
}

.npRmpi_autodispatch_is_restart_record <- function(x) {
  required <- c(
    "restart", "start", "elapsed", "status", "message", "objective",
    "bbe", "iterations", "solution"
  )
  if (!is.list(x) || is.data.frame(x) || is.null(names(x)) ||
      anyDuplicated(names(x)) || !all(required %in% names(x)))
    return(FALSE)
  is.numeric(x$restart) && length(x$restart) == 1L &&
    is.numeric(x$start) && length(x$start) > 0L &&
    is.numeric(x$elapsed) && length(x$elapsed) == 1L &&
    is.character(x$status) && length(x$status) == 1L &&
    (is.null(x$message) ||
       (is.character(x$message) && length(x$message) == 1L)) &&
    is.numeric(x$objective) && length(x$objective) == 1L &&
    is.numeric(x$bbe) && length(x$bbe) == 1L &&
    is.numeric(x$iterations) && length(x$iterations) == 1L &&
    (is.null(x$solution) || is.numeric(x$solution))
}

.npRmpi_autodispatch_fingerprint <- function(x, policy = "bandwidth.v1") {
  if (!is.character(policy) || length(policy) != 1L ||
      !policy %in% c("bandwidth.v1", "context.v1"))
    stop("unsupported autodispatch fingerprint policy", call. = FALSE)
  timing.fields <- c(
    "timing", "timing.profile", "total.time", "optim.time", "fit.time",
    "nomad.time", "powell.time", "verify.time", "child.nomad.time",
    "child.powell.time"
  )
  restart.collections <- c(
    "nomad.restarts", "nomad.restart.results", "restart.results"
  )
  normalize <- function(z) {
    if (is.environment(z))
      return(sprintf("<environment:%s>", environmentName(z)))
    if (inherits(z, "formula")) {
      environment(z) <- emptyenv()
      return(z)
    }
    if (is.call(z)) {
      out <- lapply(as.list(z), normalize)
      names(out) <- names(as.list(z))
      return(as.call(out))
    }
    if (is.pairlist(z)) {
      out <- lapply(as.list(z), normalize)
      names(out) <- names(as.list(z))
      return(as.pairlist(out))
    }
    if (is.list(z)) {
      nms <- names(z)
      if (!inherits(z, "data.frame") && !is.null(nms)) {
        # Ignore elapsed time only inside recognized package telemetry records.
        for (degree.index in which(nms == "degree.search")) {
          degree.search <- z[[degree.index]]
          degree.names <- if (is.list(degree.search)) names(degree.search) else NULL
          if (!is.null(degree.names)) {
            for (trace.index in which(degree.names == "trace")) {
              trace <- degree.search[[trace.index]]
              if (.npRmpi_autodispatch_is_degree_trace(trace)) {
                trace[["elapsed"]] <- NULL
                degree.search[trace.index] <- list(trace)
              }
            }
            z[degree.index] <- list(degree.search)
          }
        }
        for (collection.index in which(nms %in% restart.collections)) {
          records <- z[[collection.index]]
          if (is.list(records)) {
            for (record.index in seq_along(records)) {
              record <- records[[record.index]]
              if (.npRmpi_autodispatch_is_restart_record(record)) {
                record[["elapsed"]] <- NULL
                records[record.index] <- list(record)
              }
            }
            z[collection.index] <- list(records)
          }
        }
      }
      if (!inherits(z, "data.frame") && !is.null(nms))
        z[intersect(nms, timing.fields)] <- NULL
      for (i in seq_along(z))
        z[i] <- list(normalize(z[[i]]))
    }
    attrs <- attributes(z)
    if (length(attrs)) {
      attrs[intersect(names(attrs), c(
        "npRmpi.dispatch.mode", "npRmpi.autodispatch.remote",
        "npRmpi.autodispatch.fingerprint", "npRmpi.autodispatch.lease",
        "timing.profile"
      ))] <- NULL
      attrs <- lapply(attrs, normalize)
      attributes(z) <- attrs
    }
    z
  }
  y <- x
  attr(y, "npRmpi.dispatch.mode") <- NULL
  attr(y, "npRmpi.autodispatch.remote") <- NULL
  attr(y, "npRmpi.autodispatch.fingerprint") <- NULL
  attr(y, "npRmpi.autodispatch.lease") <- NULL
  y <- normalize(y)
  raw <- tryCatch(serialize(y, connection = NULL, xdr = FALSE), error = function(e) raw())
  if (!length(raw))
    return(NA_character_)
  digest <- tryCatch(.npRmpi_autodispatch_raw_md5(raw), error = function(e) NA_character_)
  if (!is.character(digest) || length(digest) != 1L || is.na(digest) || !nzchar(digest))
    return(NA_character_)
  sprintf("md5:%d:%s", length(raw), digest)
}

.npRmpi_autodispatch_has_tmp_symbols <- function(x) {
  if (is.symbol(x)) {
    nm <- as.character(x)
    return(grepl("^\\.__npRmpi_autod_", nm) || grepl("^\\.\\.[0-9]+$", nm))
  }

  if (is.call(x) || is.pairlist(x))
    return(any(vapply(as.list(x), .npRmpi_autodispatch_has_tmp_symbols, logical(1))))

  if (inherits(x, "formula"))
    return(.npRmpi_autodispatch_has_tmp_symbols(as.list(x)))

  FALSE
}

.npRmpi_autodispatch_ref_is_current <- function(val) {
  if (!identical(.npRmpi_lease_state$status, "open"))
    return(FALSE)
  fields <- .npRmpi_lease_handle_fields(val)
  if (is.null(fields) || is.null(.npRmpi_lease_state$session) ||
      !identical(fields$session, .npRmpi_lease_state$session))
    return(FALSE)
  tab <- .npRmpi_lease_capability_table()
  hit <- which(tab$kind == "leased_public" &
               tab$capability == fields$capability)
  if (length(hit) != 1L || !inherits(val, tab$expected_class[[hit]]))
    return(FALSE)
  fingerprint <- .npRmpi_autodispatch_fingerprint(
    val,
    policy = tab$fingerprint_policy[[hit]]
  )
  if (!identical(fingerprint, fields$fingerprint) ||
      !exists(fields$id, envir = .npRmpi_lease_state$registry, inherits = FALSE) ||
      !exists(fields$id, envir = .npRmpi_lease_state$store, inherits = FALSE))
    return(FALSE)
  meta <- get(fields$id, envir = .npRmpi_lease_state$registry, inherits = FALSE)
  identical(meta$state, "live") && identical(meta$session, fields$session) &&
    identical(meta$capability, fields$capability) &&
    identical(meta$fingerprint, fields$fingerprint)
}

.npRmpi_autodispatch_can_reuse_bws_ref <- function(val, call.base) {
  !is.null(.npRmpi_lease_resolve(val, consumer = call.base, argument = "bws"))
}

.npRmpi_autodispatch_large_arg_threshold <- function() {
  thr <- getOption("npRmpi.autodispatch.arg.broadcast.threshold", 4096L)
  if (!is.numeric(thr) || length(thr) != 1L || is.na(thr) || thr < 0)
    return(4096L)
  as.integer(thr)
}

.npRmpi_autodispatch_large_arg_threshold_for_call <- function(mc) {
  thr <- .npRmpi_autodispatch_large_arg_threshold()
  call.name <- .npRmpi_autodispatch_call_name(mc)
  call.name <- sub("\\..*$", "", call.name)
  ud.calls <- c("npudens", "npudensbw", "npudist", "npudistbw")
  if (call.name %in% ud.calls) {
    uthr <- getOption("npRmpi.autodispatch.arg.broadcast.threshold.unconditional", 32768L)
    if (!is.numeric(uthr) || length(uthr) != 1L || is.na(uthr) || uthr < 0)
      return(thr)
    return(as.integer(uthr))
  }

  regression.calls <- c(
    "npreg", "npregbw", "npplreg", "npplregbw",
    "nplsqreg", "nplsqregbw",
    "npqreg", "npqregbw", "npscoef", "npscoefbw",
    "npindex", "npindexbw"
  )
  if (!call.name %in% regression.calls)
    return(thr)

  rthr <- getOption("npRmpi.autodispatch.arg.broadcast.threshold.regression", 32768L)
  if (!is.numeric(rthr) || length(rthr) != 1L || is.na(rthr) || rthr < 0)
    return(thr)
  as.integer(rthr)
}

.npRmpi_autodispatch_store_tmp <- function(name,
                                           value,
                                           tmpvals,
                                           prepublish,
                                           threshold,
                                           force.inline = FALSE) {
  if (isTRUE(force.inline) || as.numeric(object.size(value)) < threshold) {
    tmpvals[[name]] <- value
  } else {
    prepublish[[name]] <- value
  }
  list(tmpvals = tmpvals, prepublish = prepublish)
}

.npRmpi_autodispatch_as_generic_call <- function(generic, mc) {
  mc <- .npRmpi_autodispatch_expand_dots_call(mc)
  args <- as.list(mc)[-1L]
  as.call(c(list(as.name(generic)), args))
}

.npRmpi_autodispatch_expand_dots_call <- function(mc) {
  if (!is.call(mc) || length(mc) < 2L)
    return(mc)

  head <- mc[[1L]]
  args <- as.list(mc)[-1L]
  arg.names <- names(args)

  if (is.null(arg.names) || !any(arg.names == "..."))
    return(mc)

  out.args <- list()
  out.names <- character(0)

  append_arg <- function(val, nm) {
    out.args[length(out.args) + 1L] <<- list(val)
    out.names[[length(out.names) + 1L]] <<- if (is.null(nm)) "" else nm
  }

  for (i in seq_along(args)) {
    nm <- arg.names[[i]]
    val <- args[[i]]

    dots.like <- identical(nm, "...") && (
      is.pairlist(val) ||
      is.list(val) ||
      (is.call(val) &&
         length(val) >= 1L &&
         is.symbol(val[[1L]]) &&
         as.character(val[[1L]]) %in% c("pairlist", "list"))
    )

    if (dots.like) {
      dot.args <- if (is.call(val)) as.list(val)[-1L] else as.list(val)
      dot.names <- names(dot.args)
      if (!length(dot.args))
        next
      for (j in seq_along(dot.args))
        append_arg(dot.args[[j]], if (is.null(dot.names)) "" else dot.names[[j]])
      next
    }

    append_arg(val, nm)
  }

  names(out.args) <- out.names
  as.call(c(list(head), out.args))
}

.npRmpi_autodispatch_active <- function() {
  isTRUE(getOption("npRmpi.autodispatch", FALSE)) &&
    !isTRUE(getOption("npRmpi.autodispatch.disable", FALSE))
}

.npRmpi_autodispatch_in_context <- function() {
  isTRUE(getOption("npRmpi.autodispatch.context", FALSE))
}

.npRmpi_manual_bcast_in_context <- function() {
  isTRUE(getOption("npRmpi.manual.bcast.context", FALSE))
}

.npRmpi_with_manual_bcast_context <- function(expr) {
  old <- getOption("npRmpi.manual.bcast.context", FALSE)
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  options(npRmpi.manual.bcast.context = TRUE)
  options(npRmpi.autodispatch.disable = TRUE)
  on.exit({
    options(npRmpi.autodispatch.disable = old.disable)
    options(npRmpi.manual.bcast.context = old)
  }, add = TRUE)
  force(expr)
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
  if (isTRUE(.npRmpi_manual_bcast_in_context()))
    return(TRUE)
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
  if (isTRUE(getOption("npRmpi.profile.active", FALSE)) &&
      isTRUE(.npRmpi_manual_bcast_in_context())) {
    return(invisible(FALSE))
  }
  if (!isTRUE(getOption("npRmpi.autodispatch.warned.nested", FALSE))) {
    .np_warning("detected active mpi.bcast.cmd context; skipping nested auto-dispatch for this call")
    options(npRmpi.autodispatch.warned.nested = TRUE)
  }
  invisible(TRUE)
}

.npRmpi_has_active_slave_pool <- function(comm = 1L) {
  if (!isTRUE(getOption("npRmpi.mpi.initialized", FALSE)))
    return(FALSE)
  if (!isTRUE(getOption("npRmpi.pool.active", FALSE)))
    return(FALSE)
  size <- tryCatch(mpi.comm.size(comm), error = function(e) NA_integer_)
  rank <- tryCatch(mpi.comm.rank(comm), error = function(e) NA_integer_)
  if (is.null(size) || is.null(rank) ||
      is.na(size) || is.na(rank))
    return(FALSE)
  size >= 2L
}

.npRmpi_rank_local_regression_context <- function(comm = 1L) {
  isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) ||
    isTRUE(.npRmpi_autodispatch_in_context()) ||
    isTRUE(.npRmpi_manual_bcast_in_context()) ||
    isTRUE(.npRmpi_autodispatch_called_from_bcast()) ||
    isTRUE(.npRmpi_has_active_slave_pool(comm = comm))
}

.npRmpi_regression_native_local_mode <- function(local.mode = FALSE,
                                                 comm = 1L) {
  local.mode <- npValidateScalarLogical(local.mode, "local.mode")
  if (local.mode)
    return(TRUE)

  ## A native regression call may use MPI collectives only when the current
  ## R expression is executing symmetrically on every rank.  An active pool
  ## alone is not sufficient: internal fit/predict helpers are also called
  ## directly on the master.  Classify those calls positively as local so a
  ## rank-local helper can never enter a collective without its workers.
  if (isTRUE(.npRmpi_autodispatch_in_context()) ||
      isTRUE(.npRmpi_autodispatch_called_from_bcast()))
    return(FALSE)

  .npRmpi_has_active_slave_pool(comm = comm)
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
                                              where = "this call") {
  .npRmpi_abort_if_rmpi_attached(where = where)
  if (isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
    return(invisible(TRUE))
  if (isTRUE(.npRmpi_autodispatch_in_context()) ||
      isTRUE(.npRmpi_manual_bcast_in_context()))
    return(invisible(TRUE))
  if (.npRmpi_has_active_slave_pool(comm = comm))
    return(invisible(TRUE))
  stop(sprintf("%s requires an active MPI slave pool; call npRmpi.init(...) first", where))
}

.npRmpi_spmd_seq_get <- function() {
  seq.id <- getOption("npRmpi.spmd.seq_id", 0L)
  if (!is.numeric(seq.id) || length(seq.id) != 1L || is.na(seq.id) || seq.id < 0)
    return(0L)
  as.integer(seq.id)
}

.npRmpi_spmd_seq_set <- function(seq.id) {
  if (!is.numeric(seq.id) || length(seq.id) != 1L || is.na(seq.id) ||
      !is.finite(seq.id) || seq.id < 0 || seq.id != floor(seq.id))
    stop("'seq.id' must be a non-negative integer", call. = FALSE)
  options(npRmpi.spmd.seq_id = as.integer(seq.id))
  invisible(as.integer(seq.id))
}

.npRmpi_spmd_make_envelope <- function(opcode,
                                       args_ref = NULL,
                                       timeout_class = "default",
                                       seq_id = NULL) {
  if (!is.character(opcode) || length(opcode) != 1L || is.na(opcode) || !nzchar(opcode))
    stop("'opcode' must be a non-empty character scalar", call. = FALSE)
  if (!is.character(timeout_class) || length(timeout_class) != 1L ||
      is.na(timeout_class) || !nzchar(timeout_class))
    stop("'timeout_class' must be a non-empty character scalar", call. = FALSE)

  if (is.null(seq_id)) {
    seq_id <- .npRmpi_spmd_seq_get() + 1L
  } else {
    if (!is.numeric(seq_id) || length(seq_id) != 1L || is.na(seq_id) ||
        !is.finite(seq_id) || seq_id < 1 || seq_id != floor(seq_id))
      stop("'seq_id' must be a positive integer", call. = FALSE)
    seq_id <- as.integer(seq_id)
  }

  list(
    seq_id = as.integer(seq_id),
    opcode = opcode,
    args_ref = args_ref,
    timeout_class = timeout_class
  )
}

.npRmpi_spmd_validate_envelope <- function(envelope) {
  if (!is.list(envelope))
    stop("SPMD envelope must be a list", call. = FALSE)
  req <- c("seq_id", "opcode", "args_ref", "timeout_class")
  miss <- setdiff(req, names(envelope))
  if (length(miss))
    stop(sprintf("SPMD envelope missing required fields: %s", paste(miss, collapse = ", ")), call. = FALSE)

  seq.id <- envelope$seq_id
  if (!is.numeric(seq.id) || length(seq.id) != 1L || is.na(seq.id) ||
      !is.finite(seq.id) || seq.id < 1 || seq.id != floor(seq.id))
    stop("SPMD envelope has invalid 'seq_id'", call. = FALSE)
  if (!is.character(envelope$opcode) || length(envelope$opcode) != 1L ||
      is.na(envelope$opcode) || !nzchar(envelope$opcode))
    stop("SPMD envelope has invalid 'opcode'", call. = FALSE)
  if (!is.character(envelope$timeout_class) || length(envelope$timeout_class) != 1L ||
      is.na(envelope$timeout_class) || !nzchar(envelope$timeout_class))
    stop("SPMD envelope has invalid 'timeout_class'", call. = FALSE)

  invisible(TRUE)
}

.npRmpi_spmd_assert_sequence <- function(envelope,
                                         last_seq = .npRmpi_spmd_seq_get(),
                                         where = "SPMD step") {
  .npRmpi_spmd_validate_envelope(envelope)
  expected <- as.integer(last_seq) + 1L
  got <- as.integer(envelope$seq_id)
  if (!identical(got, expected)) {
    stop(
      sprintf("%s sequence mismatch: expected seq_id=%d but received seq_id=%d [opcode=%s]",
              where, expected, got, envelope$opcode),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.npRmpi_spmd_registry <- new.env(parent = emptyenv())

.npRmpi_spmd_register_opcode <- function(opcode, handler) {
  if (!is.character(opcode) || length(opcode) != 1L || is.na(opcode) || !nzchar(opcode))
    stop("'opcode' must be a non-empty character scalar", call. = FALSE)
  if (!is.function(handler))
    stop("'handler' must be a function", call. = FALSE)
  assign(opcode, handler, envir = .npRmpi_spmd_registry)
  invisible(TRUE)
}

.npRmpi_spmd_eval_payload <- function(payload, envelope) {
  if (!is.list(payload))
    stop("SPMD payload must be a list", call. = FALSE)

  call.obj <- payload$call
  if (!(is.call(call.obj) || is.expression(call.obj)))
    stop("SPMD payload missing executable call object", call. = FALSE)

  opt.keys <- payload$opt.keys
  opt.vals <- payload$opt.vals
  opt.verify <- isTRUE(payload$opt.verify)
  tmpvals <- payload$tmpvals
  tmpnames <- payload$tmpnames
  prepublish.names <- payload$prepublish.names
  lease.bindings <- payload$lease.bindings
  publication <- payload$publication
  if (is.null(publication))
    publication <- list(kind = "none", producer = "unknown")

  if (!is.null(opt.keys) && length(opt.keys)) {
    for (i in seq_along(opt.keys))
      options(structure(list(opt.vals[[i]]), names = opt.keys[[i]]))
    if (opt.verify) {
      for (i in seq_along(opt.keys)) {
        lval <- getOption(opt.keys[[i]])
        if (!identical(lval, opt.vals[[i]]))
          stop(sprintf("failed to synchronize option '%s' across MPI ranks", opt.keys[[i]]), call. = FALSE)
      }
    }
  }

  if (!is.null(tmpvals) && length(tmpvals)) {
    for (nm in names(tmpvals))
      .GlobalEnv[[nm]] <- tmpvals[[nm]]
  }
  if (!is.null(lease.bindings) && length(lease.bindings)) {
    for (nm in names(lease.bindings)) {
      id <- as.character(lease.bindings[[nm]])[1L]
      if (exists(nm, envir = .GlobalEnv, inherits = FALSE))
        stop(sprintf("autodispatch lease binding collision for '%s'", nm), call. = FALSE)
      if (!exists(id, envir = .npRmpi_lease_state$store, inherits = FALSE))
        stop(sprintf("autodispatch lease '%s' is unavailable", id), call. = FALSE)
      meta <- if (exists(id, envir = .npRmpi_lease_state$registry, inherits = FALSE))
        get(id, envir = .npRmpi_lease_state$registry, inherits = FALSE) else NULL
      if (!is.list(meta) || !identical(meta$state, "live"))
        stop(sprintf("autodispatch lease '%s' is not live", id), call. = FALSE)
      assign(nm, get(id, envir = .npRmpi_lease_state$store, inherits = FALSE),
             envir = .GlobalEnv)
    }
  }

  old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  options(npRmpi.autodispatch.context = TRUE)
  options(npRmpi.autodispatch.disable = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  if (!is.null(tmpnames) && length(tmpnames))
    on.exit(get(".npRmpi_rm_existing", envir = asNamespace("npRmpi"), inherits = FALSE)(tmpnames, envir = .GlobalEnv), add = TRUE)
  if (!is.null(lease.bindings) && length(lease.bindings))
    on.exit(get(".npRmpi_rm_existing", envir = asNamespace("npRmpi"), inherits = FALSE)(names(lease.bindings), envir = .GlobalEnv), add = TRUE)

  res <- .npRmpi_eval_scmd(call.obj, envir = .GlobalEnv)
  tmpreplace <- tmpvals
  if (!is.null(prepublish.names) && length(prepublish.names)) {
    prepublish.names <- unique(as.character(prepublish.names))
    for (nm in prepublish.names) {
      if (nzchar(nm) && exists(nm, envir = .GlobalEnv, inherits = FALSE))
        tmpreplace[[nm]] <- get(nm, envir = .GlobalEnv, inherits = FALSE)
    }
  }
  if (!is.null(tmpreplace) && length(tmpreplace))
    res <- .npRmpi_autodispatch_sanitize_object(res, tmpvals = tmpreplace)
  publication.ack <- .npRmpi_lease_prepare_local(res, publication)
  if (identical(publication$kind, "leased_public"))
    .npRmpi_lease_commit_local(publication)
  structure(
    list(result = res, publication = publication.ack),
    class = "npRmpi_spmd_internal_result"
  )
}

.npRmpi_lease_lifecycle_handler <- function(payload, envelope) {
  if (!is.list(payload))
    stop("autodispatch lifecycle payload must be a list", call. = FALSE)
  switch(
    envelope$opcode,
    "autodispatch.lifecycle.commit" = .npRmpi_lease_commit_local(payload$plan),
    "autodispatch.lifecycle.rollback" = .npRmpi_lease_rollback_local(payload$ids),
    "autodispatch.lifecycle.retire" = .npRmpi_lease_rollback_local(payload$ids),
    "autodispatch.lifecycle.session_drain" = .npRmpi_lease_session_drain_local(),
    stop("unsupported autodispatch lifecycle opcode", call. = FALSE)
  )
}

.npRmpi_spmd_locked_opcodes <- function() {
  c(
    "autodispatch.npregbw.cv_lllp",
    "autodispatch.npscoefbw.cv_lllp",
    "autodispatch.npplregbw.cv_lllp",
    "autodispatch.npindexbw.core",
    "autodispatch.npreg.core",
    "autodispatch.npscoef.core",
    "autodispatch.npplreg.core",
    "autodispatch.npindex.core",
    "autodispatch.npudens.core",
    "autodispatch.npudist.core",
    "autodispatch.npcdens.core",
    "autodispatch.npcdist.core",
    "autodispatch.npudensbw.cv",
    "autodispatch.npudistbw.cv",
    "autodispatch.npcdensbw.cv",
    "autodispatch.npcdistbw.cv",
    "autodispatch.nplsqreg.core",
    "autodispatch.npqreg.core",
    "autodispatch.npconmode.core",
    "autodispatch.npksum.core",
    "autodispatch.npregiv.core",
    "autodispatch.npregivderiv.core",
    "autodispatch.npcmstest.core",
    "autodispatch.npqcmstest.core",
    "autodispatch.npdeneqtest.core",
    "autodispatch.npdeptest.core",
    "autodispatch.npsdeptest.core",
    "autodispatch.npsigtest.core",
    "autodispatch.npsymtest.core",
    "autodispatch.npunitest.core"
    ,"autodispatch.npscoefbw.nomad_context"
    ,"autodispatch.lifecycle.commit"
    ,"autodispatch.lifecycle.rollback"
    ,"autodispatch.lifecycle.retire"
    ,"autodispatch.lifecycle.session_drain"
  )
}

.npRmpi_spmd_eval_payload_call_guard <- function(payload,
                                                 envelope,
                                                 allowed_calls,
                                                 where = "SPMD call guard") {
  call.obj <- payload$call
  if (is.expression(call.obj)) {
    if (length(call.obj) != 1L)
      stop(sprintf("%s expected single expression call", where), call. = FALSE)
    call.obj <- call.obj[[1L]]
  }
  if (!is.call(call.obj))
    stop(sprintf("%s missing call payload", where), call. = FALSE)

  call.name <- .npRmpi_autodispatch_call_name(call.obj)
  if (!length(call.name) || !nzchar(call.name) || !(call.name %in% allowed_calls)) {
    stop(
      sprintf("%s opcode '%s' restricted to {%s}; received '%s'",
              where,
              as.character(envelope$opcode)[1L],
              paste(allowed_calls, collapse = ", "),
              if (length(call.name) && nzchar(call.name)) call.name else "<unknown>"),
      call. = FALSE
    )
  }

  payload$call <- call.obj
  .npRmpi_spmd_eval_payload(payload = payload, envelope = envelope)
}

.npRmpi_spmd_get_opcode <- function(opcode) {
  if (!exists(opcode, envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    locked <- .npRmpi_spmd_locked_opcodes()
    if (is.character(opcode) && length(opcode) == 1L && opcode %in% locked) {
      stop(sprintf("SPMD locked opcode '%s' is not registered", opcode), call. = FALSE)
    }
    if (is.character(opcode) && length(opcode) == 1L &&
        startsWith(opcode, "autodispatch.")) {
      .npRmpi_spmd_register_opcode(
        opcode,
        function(payload, envelope) .npRmpi_spmd_eval_payload(payload = payload, envelope = envelope)
      )
    } else {
      stop(sprintf("SPMD opcode '%s' is not registered", opcode), call. = FALSE)
    }
  }
  get(opcode, envir = .npRmpi_spmd_registry, inherits = FALSE)
}

.npRmpi_spmd_register_defaults <- function() {
  if (!exists("spmd.ping", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode("spmd.ping", function(payload, envelope) {
      rank <- tryCatch(as.integer(mpi.comm.rank(1L)), error = function(e) 0L)
      list(
        opcode = envelope$opcode,
        seq_id = as.integer(envelope$seq_id),
        rank = as.integer(rank),
        payload = payload
      )
    })
  }
  if (!exists("autodispatch.npscoefbw.nomad_context", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npscoefbw.nomad_context",
      function(payload, envelope) {
        if (!is.list(payload) || !.npRmpi_lease_is_context_call(payload$call))
          stop("SPMD npscoef NOMAD context opcode received an unauthorized call", call. = FALSE)
        payload$call[[1L]] <- get(
          ".npscoefbw_nomad_context_prepare",
          envir = asNamespace("npRmpi"),
          mode = "function",
          inherits = FALSE
        )
        .npRmpi_spmd_eval_payload(payload = payload, envelope = envelope)
      }
    )
  }
  for (opcode in c("autodispatch.lifecycle.commit",
                   "autodispatch.lifecycle.rollback",
                   "autodispatch.lifecycle.retire",
                   "autodispatch.lifecycle.session_drain")) {
    if (!exists(opcode, envir = .npRmpi_spmd_registry, inherits = FALSE))
      .npRmpi_spmd_register_opcode(opcode, .npRmpi_lease_lifecycle_handler)
  }
  if (!exists("autodispatch.npregbw.cv_lllp", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npregbw.cv_lllp",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npregbw",
        where = "SPMD regression opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npscoefbw.cv_lllp", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npscoefbw.cv_lllp",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npscoefbw",
        where = "SPMD regression opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npplregbw.cv_lllp", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npplregbw.cv_lllp",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npplregbw",
        where = "SPMD regression opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npindexbw.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npindexbw.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npindexbw", "npindexbw.default", "npindexbw.formula", "npindexbw.sibandwidth"),
        where = "SPMD regression opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npreg.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npreg.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npreg", "npreg.default", "npreg.formula", "npreg.call", "npreg.rbandwidth"),
        where = "SPMD estimator opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npscoef.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npscoef.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npscoef", "npscoef.default", "npscoef.formula", "npscoef.call", "npscoef.scbandwidth"),
        where = "SPMD estimator opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npplreg.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npplreg.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npplreg", "npplreg.default", "npplreg.formula", "npplreg.call", "npplreg.plbandwidth"),
        where = "SPMD estimator opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npindex.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npindex.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npindex", "npindex.default", "npindex.formula", "npindex.call", "npindex.sibandwidth"),
        where = "SPMD estimator opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npudens.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npudens.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npudens", "npudens.default", "npudens.formula", "npudens.call", "npudens.bandwidth"),
        where = "SPMD estimator opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npudist.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npudist.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npudist", "npudist.default", "npudist.formula", "npudist.call", "npudist.dbandwidth"),
        where = "SPMD estimator opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npcdens.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npcdens.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npcdens", "npcdens.default", "npcdens.formula", "npcdens.call", "npcdens.conbandwidth"),
        where = "SPMD estimator opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npcdist.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npcdist.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npcdist", "npcdist.default", "npcdist.formula", "npcdist.call", "npcdist.condbandwidth"),
        where = "SPMD estimator opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npudensbw.cv", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npudensbw.cv",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npudensbw",
        where = "SPMD density opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npudistbw.cv", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npudistbw.cv",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npudistbw",
        where = "SPMD density opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npcdensbw.cv", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npcdensbw.cv",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npcdensbw",
        where = "SPMD density opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npcdistbw.cv", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npcdistbw.cv",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npcdistbw",
        where = "SPMD density opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npqreg.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npqreg.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npqreg", "npqreg.default", "npqreg.formula", "npqreg.call", "npqreg.condbandwidth", "npqreg.conbandwidth"),
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.nplsqreg.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.nplsqreg.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("nplsqreg", "nplsqreg.default", "nplsqreg.formula", "nplsqreg.lsqregressionbandwidth",
                          "nplsqregbw", "nplsqregbw.default", "nplsqregbw.formula", "nplsqregbw.lsqregressionbandwidth"),
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npconmode.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npconmode.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npconmode", "npconmode.default", "npconmode.formula", "npconmode.call", "npconmode.conbandwidth"),
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npksum.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npksum.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npksum", "npksum.formula", "npksum.default", "npksum.numeric", "npksum.integer"),
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npregiv.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npregiv.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npregiv",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npregivderiv.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npregivderiv.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npregivderiv",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npcmstest.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npcmstest.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npcmstest",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npqcmstest.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npqcmstest.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npqcmstest",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npdeneqtest.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npdeneqtest.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npdeneqtest",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npdeptest.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npdeptest.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npdeptest",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npsdeptest.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npsdeptest.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npsdeptest",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npsigtest.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npsigtest.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = c("npsigtest", "npsigtest.formula", "npsigtest.call", "npsigtest.npregression", "npsigtest.rbandwidth"),
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npsymtest.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npsymtest.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npsymtest",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  if (!exists("autodispatch.npunitest.core", envir = .npRmpi_spmd_registry, inherits = FALSE)) {
    .npRmpi_spmd_register_opcode(
      "autodispatch.npunitest.core",
      function(payload, envelope) .npRmpi_spmd_eval_payload_call_guard(
        payload = payload,
        envelope = envelope,
        allowed_calls = "npunitest",
        where = "SPMD non-core opcode guard"
      )
    )
  }
  invisible(TRUE)
}

.npRmpi_spmd_execute_local <- function(envelope,
                                       payload = NULL,
                                       where = "SPMD step local execute") {
  .npRmpi_spmd_register_defaults()
  last.seq <- .npRmpi_spmd_seq_get()
  .npRmpi_spmd_assert_sequence(envelope, last_seq = last.seq, where = where)
  .npRmpi_spmd_seq_set(envelope$seq_id)

  handler <- .npRmpi_spmd_get_opcode(envelope$opcode)
  out <- tryCatch(
    handler(payload, envelope),
    error = function(e) {
      stop(
        sprintf("%s opcode failure [opcode=%s seq_id=%d]: %s",
                where, envelope$opcode, envelope$seq_id, conditionMessage(e)),
        call. = FALSE
      )
    }
  )
  publication <- list(kind = "", id = "", capability = "", fingerprint = "")
  if (inherits(out, "npRmpi_spmd_internal_result")) {
    publication <- out$publication
    out <- out$result
  }

  list(
    ok = TRUE,
    ack = list(
      seq_id = as.integer(envelope$seq_id),
      opcode = envelope$opcode,
      status = "ACK",
      publication_kind = as.character(publication$kind)[1L],
      publication_id = as.character(publication$id)[1L],
      publication_capability = as.character(publication$capability)[1L],
      publication_fingerprint = as.character(publication$fingerprint)[1L]
    ),
    result = out
  )
}

.npRmpi_spmd_try_execute_local <- function(envelope,
                                           payload = NULL,
                                           where = "SPMD step local execute") {
  tryCatch(
    .npRmpi_spmd_execute_local(envelope = envelope, payload = payload, where = where),
    error = function(e) list(
      ok = FALSE,
      ack = list(
        seq_id = if (is.list(envelope) && !is.null(envelope$seq_id)) as.integer(envelope$seq_id) else NA_integer_,
        opcode = if (is.list(envelope) && !is.null(envelope$opcode)) as.character(envelope$opcode)[1L] else NA_character_,
        status = "ERR",
        error = conditionMessage(e),
        publication_kind = "",
        publication_id = "",
        publication_capability = "",
        publication_fingerprint = ""
      ),
      error = conditionMessage(e)
    )
  )
}

.npRmpi_spmd_timeout_seconds <- function(timeout_class = "default") {
  cls <- as.character(timeout_class)[1L]
  if (!is.character(cls) || length(cls) != 1L || is.na(cls) || !nzchar(cls))
    cls <- "default"
  opt.key <- paste0("npRmpi.spmd.timeout.", gsub("[^A-Za-z0-9_.-]+", "_", cls))
  timeout <- getOption(opt.key, NULL)
  if (is.null(timeout))
    timeout <- getOption("npRmpi.spmd.timeout.default", 0)
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) ||
      !is.finite(timeout) || timeout <= 0)
    return(0)
  as.numeric(timeout)
}

.npRmpi_spmd_with_timeout <- function(timeout_sec, thunk) {
  if (!is.function(thunk))
    stop("'thunk' must be a function", call. = FALSE)
  if (!is.numeric(timeout_sec) || length(timeout_sec) != 1L || is.na(timeout_sec) ||
      !is.finite(timeout_sec) || timeout_sec <= 0)
    return(thunk())

  base::setTimeLimit(elapsed = as.numeric(timeout_sec), transient = TRUE)
  on.exit(base::setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  thunk()
}

.npRmpi_spmd_collective_ack <- function(local,
                                        envelope,
                                        comm = 1L,
                                        where = "SPMD step") {
  .npRmpi_spmd_validate_envelope(envelope)
  if (!is.list(local) || is.null(local$ack))
    stop(sprintf("%s local ACK payload is malformed", where), call. = FALSE)

  size <- tryCatch(as.integer(mpi.comm.size(comm)), error = function(e) 1L)
  if (is.na(size) || size < 2L) {
    if (!isTRUE(local$ok)) {
      stop(
        sprintf("%s local failure [opcode=%s seq_id=%d]: %s",
                where, envelope$opcode, as.integer(envelope$seq_id),
                as.character(local$error)[1L]),
        call. = FALSE
      )
    }
    return(local)
  }

  rank <- tryCatch(as.integer(mpi.comm.rank(comm)), error = function(e) NA_integer_)
  local.seq <- if (is.list(local$ack) && !is.null(local$ack$seq_id)) as.integer(local$ack$seq_id) else NA_integer_
  local.op <- if (is.list(local$ack) && !is.null(local$ack$opcode)) as.character(local$ack$opcode)[1L] else as.character(envelope$opcode)[1L]
  local.status <- if (isTRUE(local$ok)) "ACK" else "ERR"
  local.error <- if (isTRUE(local$ok)) "" else as.character(local$error)[1L]
  publication.field <- function(name) {
    val <- local$ack[[name]]
    if (is.null(val) || !length(val) || is.na(val[[1L]])) "" else as.character(val[[1L]])
  }
  ack.local <- c(
    rank = as.character(rank),
    seq_id = as.character(local.seq),
    opcode = local.op,
    status = local.status,
    error = local.error,
    publication_kind = publication.field("publication_kind"),
    publication_id = publication.field("publication_id"),
    publication_capability = publication.field("publication_capability"),
    publication_fingerprint = publication.field("publication_fingerprint")
  )

  timeout.sec <- .npRmpi_spmd_timeout_seconds(timeout_class = envelope$timeout_class)
  ack.gather <- tryCatch(
    .npRmpi_spmd_with_timeout(
      timeout_sec = timeout.sec,
      thunk = function() mpi.allgather.Robj(ack.local, comm = comm)
    ),
    error = function(e) {
      stop(
        sprintf("%s ACK gather failed [opcode=%s seq_id=%d]: %s",
                where, envelope$opcode, as.integer(envelope$seq_id), conditionMessage(e)),
        call. = FALSE
      )
    }
  )

  flat <- as.character(ack.gather)
  expected <- 9L * size
  if (length(flat) != expected) {
    stop(
      sprintf("%s ACK gather shape mismatch [opcode=%s seq_id=%d]: expected %d fields, received %d",
              where, envelope$opcode, as.integer(envelope$seq_id), expected, length(flat)),
      call. = FALSE
    )
  }

  ack.mat <- matrix(flat, nrow = 9L, byrow = FALSE,
                    dimnames = list(c("rank", "seq_id", "opcode", "status", "error",
                                      "publication_kind", "publication_id",
                                      "publication_capability", "publication_fingerprint"), NULL))
  seq.vec <- suppressWarnings(as.integer(ack.mat["seq_id", ]))
  rank.vec <- suppressWarnings(as.integer(ack.mat["rank", ]))
  status.vec <- ack.mat["status", ]
  opcode.vec <- ack.mat["opcode", ]
  bad <- which(is.na(seq.vec) |
                 (seq.vec != as.integer(envelope$seq_id)) |
                 (status.vec != "ACK") |
                 (opcode.vec != as.character(envelope$opcode)[1L]))
  if (!length(bad)) {
    publication.rows <- c("publication_kind", "publication_id",
                          "publication_capability", "publication_fingerprint")
    for (row in publication.rows)
      bad <- union(bad, which(ack.mat[row, ] != ack.mat[row, 1L]))
  }

  if (length(bad)) {
    detail.index <- unique(c(1L, bad))
    details <- vapply(
      detail.index,
      function(i) {
        sprintf(paste0(
          "rank=%s status=%s seq_id=%s opcode=%s err=%s ",
          "publication_kind=%s publication_id=%s publication_capability=%s ",
          "publication_fingerprint=%s"
        ),
                if (is.na(rank.vec[[i]])) "NA" else as.character(rank.vec[[i]]),
                as.character(status.vec[[i]]),
                as.character(ack.mat["seq_id", i]),
                as.character(opcode.vec[[i]]),
                as.character(ack.mat["error", i]),
                as.character(ack.mat["publication_kind", i]),
                as.character(ack.mat["publication_id", i]),
                as.character(ack.mat["publication_capability", i]),
                as.character(ack.mat["publication_fingerprint", i]))
      },
      character(1)
    )
    if (!is.na(rank) && rank != 0L)
      return(local)
    stop(
      sprintf("%s ACK mismatch [opcode=%s seq_id=%d]: %s",
              where, envelope$opcode, as.integer(envelope$seq_id),
              paste(details, collapse = "; ")),
      call. = FALSE
    )
  }

  local
}

.npRmpi_spmd_execute_step <- function(envelope,
                                      payload = NULL,
                                      comm = 1L,
                                      where = "SPMD step") {
  local <- .npRmpi_spmd_try_execute_local(
    envelope = envelope,
    payload = payload,
    where = where
  )
  .npRmpi_spmd_collective_ack(
    local = local,
    envelope = envelope,
    comm = comm,
    where = where
  )
}

.npRmpi_spmd_tiny_smoke <- function(label = "spmd-smoke", comm = 1L) {
  .npRmpi_require_active_slave_pool(comm = comm, where = "SPMD tiny smoke")
  envelope <- .npRmpi_spmd_make_envelope(
    opcode = "spmd.ping",
    timeout_class = "smoke"
  )
  payload <- list(label = as.character(label)[1L], comm = as.integer(comm))
  cmd <- substitute({
    exec.step <- get(".npRmpi_spmd_execute_step", envir = asNamespace("npRmpi"), inherits = FALSE)
    exec.step(
      envelope = ENVELOPE,
      payload = PAYLOAD,
      comm = COMM,
      where = "SPMD tiny smoke"
    )
  }, list(ENVELOPE = envelope, PAYLOAD = payload, COMM = comm))

  .npRmpi_bcast_cmd_expr(expr = cmd, comm = comm, caller.execute = TRUE)
}

.npRmpi_set_validated_backend_seed <- function(seed, comm = 1L) {
  # Calls already executing on every rank must remain local. In particular,
  # the command broadcast below re-enters npseed() under this context.
  if (isTRUE(.npRmpi_autodispatch_in_context()) ||
      isTRUE(.npRmpi_manual_bcast_in_context())) {
    .np_backend_seed_local(seed)
    return(invisible())
  }

  if (!isTRUE(getOption("npRmpi.pool.active", FALSE))) {
    .np_backend_seed_local(seed)
    return(invisible())
  }
  if (!isTRUE(.npRmpi_has_active_slave_pool(comm = comm))) {
    stop(
      "npseed() cannot synchronize the C backend because the active MPI pool is inconsistent",
      call. = FALSE
    )
  }

  # The validated seed write is total and non-allocating. Reuse the canonical
  # command path so every rank receives the write before any following command;
  # the manual-broadcast context above prevents nested dispatch.
  cmd <- substitute(npseed(SEED), list(SEED = seed))

  .npRmpi_bcast_cmd_expr(expr = cmd, comm = comm, caller.execute = TRUE)
  invisible()
}

.npRmpi_autodispatch_eval_char_arg <- function(mc, caller_env, argname) {
  arg.list <- as.list(mc)
  nms <- names(arg.list)
  if (is.null(nms) || !any(nms == argname))
    return(NULL)
  idx <- which(nms == argname)[1L]
  val <- tryCatch(.npRmpi_autodispatch_eval_arg(arg.list[[idx]], caller_env = caller_env),
                  error = function(e) NULL)
  if (is.null(val))
    return(NULL)
  val <- as.character(val)
  if (!length(val) || is.na(val[[1L]]) || !nzchar(val[[1L]]))
    return(NULL)
  tolower(val[[1L]])
}

.npRmpi_autodispatch_eval_logical_arg <- function(mc, caller_env, argname, default = NULL) {
  arg.list <- as.list(mc)
  nms <- names(arg.list)
  if (is.null(nms) || !any(nms == argname))
    return(default)
  idx <- which(nms == argname)[1L]
  val <- tryCatch(.npRmpi_autodispatch_eval_arg(arg.list[[idx]], caller_env = caller_env),
                  error = function(e) NULL)
  if (is.null(val) || length(val) != 1L)
    return(default)
  val <- as.logical(val[[1L]])
  if (is.na(val))
    return(default)
  isTRUE(val)
}

.npRmpi_autodispatch_is_bw_lllp_cv <- function(mc, caller_env, call_name) {
  if (!identical(.npRmpi_autodispatch_call_name(mc), call_name))
    return(FALSE)
  regtype <- .npRmpi_autodispatch_eval_char_arg(mc, caller_env, "regtype")
  bwmethod <- .npRmpi_autodispatch_eval_char_arg(mc, caller_env, "bwmethod")
  if (is.null(regtype) || is.null(bwmethod))
    return(FALSE)

  regtype %in% c("ll", "lp") && bwmethod %in% c("cv.ls", "cv.aic")
}

.npRmpi_autodispatch_is_density_bw_cv <- function(mc, caller_env, call_name) {
  if (!identical(.npRmpi_autodispatch_call_name(mc), call_name))
    return(FALSE)
  bw.compute <- .npRmpi_autodispatch_eval_logical_arg(
    mc = mc,
    caller_env = caller_env,
    argname = "bandwidth.compute",
    default = TRUE
  )
  isTRUE(bw.compute)
}

.npRmpi_autodispatch_is_npregbw_lllp_cv <- function(mc, caller_env) {
  .npRmpi_autodispatch_is_bw_lllp_cv(mc = mc, caller_env = caller_env, call_name = "npregbw")
}

.npRmpi_autodispatch_is_npscoefbw_lllp_cv <- function(mc, caller_env) {
  .npRmpi_autodispatch_is_bw_lllp_cv(mc = mc, caller_env = caller_env, call_name = "npscoefbw")
}

.npRmpi_autodispatch_is_npplregbw_lllp_cv <- function(mc, caller_env) {
  .npRmpi_autodispatch_is_bw_lllp_cv(mc = mc, caller_env = caller_env, call_name = "npplregbw")
}

.npRmpi_autodispatch_is_npindexbw_core <- function(mc, caller_env) {
  call.name <- .npRmpi_autodispatch_call_name(mc)
  if (!(call.name %in% c("npindexbw", "npindexbw.default", "npindexbw.formula", "npindexbw.sibandwidth")))
    return(FALSE)
  bw.compute <- .npRmpi_autodispatch_eval_logical_arg(
    mc = mc,
    caller_env = caller_env,
    argname = "bandwidth.compute",
    default = TRUE
  )
  isTRUE(bw.compute)
}

.npRmpi_autodispatch_call_name_in <- function(mc, valid_names) {
  call.name <- .npRmpi_autodispatch_call_name(mc)
  call.name %in% valid_names
}

.npRmpi_autodispatch_is_npreg_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npreg", "npreg.default", "npreg.formula", "npreg.call", "npreg.rbandwidth"))
}

.npRmpi_autodispatch_is_npscoef_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npscoef", "npscoef.default", "npscoef.formula", "npscoef.call", "npscoef.scbandwidth"))
}

.npRmpi_autodispatch_is_npplreg_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npplreg", "npplreg.default", "npplreg.formula", "npplreg.call", "npplreg.plbandwidth"))
}

.npRmpi_autodispatch_is_npindex_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npindex", "npindex.default", "npindex.formula", "npindex.call", "npindex.sibandwidth"))
}

.npRmpi_autodispatch_is_npudens_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npudens", "npudens.default", "npudens.formula", "npudens.call", "npudens.bandwidth"))
}

.npRmpi_autodispatch_is_npudist_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npudist", "npudist.default", "npudist.formula", "npudist.call", "npudist.dbandwidth"))
}

.npRmpi_autodispatch_is_npcdens_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npcdens", "npcdens.default", "npcdens.formula", "npcdens.call", "npcdens.conbandwidth"))
}

.npRmpi_autodispatch_is_npcdist_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npcdist", "npcdist.default", "npcdist.formula", "npcdist.call", "npcdist.condbandwidth"))
}

.npRmpi_autodispatch_is_npqreg_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npqreg", "npqreg.default", "npqreg.formula", "npqreg.call", "npqreg.condbandwidth", "npqreg.conbandwidth"))
}

.npRmpi_autodispatch_is_nplsqreg_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("nplsqreg", "nplsqreg.default", "nplsqreg.formula", "nplsqreg.lsqregressionbandwidth",
                                          "nplsqregbw", "nplsqregbw.default", "nplsqregbw.formula", "nplsqregbw.lsqregressionbandwidth"))
}

.npRmpi_autodispatch_is_npconmode_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npconmode", "npconmode.default", "npconmode.formula", "npconmode.call", "npconmode.conbandwidth"))
}

.npRmpi_autodispatch_is_npksum_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npksum", "npksum.formula", "npksum.default", "npksum.numeric", "npksum.integer"))
}

.npRmpi_autodispatch_is_npregiv_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npregiv")
}

.npRmpi_autodispatch_is_npregivderiv_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npregivderiv")
}

.npRmpi_autodispatch_is_npcmstest_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npcmstest")
}

.npRmpi_autodispatch_is_npqcmstest_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npqcmstest")
}

.npRmpi_autodispatch_is_npdeneqtest_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npdeneqtest")
}

.npRmpi_autodispatch_is_npdeptest_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npdeptest")
}

.npRmpi_autodispatch_is_npsdeptest_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npsdeptest")
}

.npRmpi_autodispatch_is_npsigtest_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, c("npsigtest", "npsigtest.formula", "npsigtest.call", "npsigtest.npregression", "npsigtest.rbandwidth"))
}

.npRmpi_autodispatch_is_npsymtest_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npsymtest")
}

.npRmpi_autodispatch_is_npunitest_core <- function(mc) {
  .npRmpi_autodispatch_call_name_in(mc, "npunitest")
}

.npRmpi_spmd_opcode_from_call <- function(mc, caller_env) {
  if (.npRmpi_lease_is_context_call(mc))
    return("autodispatch.npscoefbw.nomad_context")
  if (.npRmpi_autodispatch_is_npregbw_lllp_cv(mc = mc, caller_env = caller_env))
    return("autodispatch.npregbw.cv_lllp")
  if (.npRmpi_autodispatch_is_npscoefbw_lllp_cv(mc = mc, caller_env = caller_env))
    return("autodispatch.npscoefbw.cv_lllp")
  if (.npRmpi_autodispatch_is_npplregbw_lllp_cv(mc = mc, caller_env = caller_env))
    return("autodispatch.npplregbw.cv_lllp")
  if (.npRmpi_autodispatch_is_npindexbw_core(mc = mc, caller_env = caller_env))
    return("autodispatch.npindexbw.core")
  if (.npRmpi_autodispatch_is_npreg_core(mc = mc))
    return("autodispatch.npreg.core")
  if (.npRmpi_autodispatch_is_npscoef_core(mc = mc))
    return("autodispatch.npscoef.core")
  if (.npRmpi_autodispatch_is_npplreg_core(mc = mc))
    return("autodispatch.npplreg.core")
  if (.npRmpi_autodispatch_is_npindex_core(mc = mc))
    return("autodispatch.npindex.core")
  if (.npRmpi_autodispatch_is_npudens_core(mc = mc))
    return("autodispatch.npudens.core")
  if (.npRmpi_autodispatch_is_npudist_core(mc = mc))
    return("autodispatch.npudist.core")
  if (.npRmpi_autodispatch_is_npcdens_core(mc = mc))
    return("autodispatch.npcdens.core")
  if (.npRmpi_autodispatch_is_npcdist_core(mc = mc))
    return("autodispatch.npcdist.core")
  if (.npRmpi_autodispatch_is_nplsqreg_core(mc = mc))
    return("autodispatch.nplsqreg.core")
  if (.npRmpi_autodispatch_is_npqreg_core(mc = mc))
    return("autodispatch.npqreg.core")
  if (.npRmpi_autodispatch_is_npconmode_core(mc = mc))
    return("autodispatch.npconmode.core")
  if (.npRmpi_autodispatch_is_npksum_core(mc = mc))
    return("autodispatch.npksum.core")
  if (.npRmpi_autodispatch_is_npregiv_core(mc = mc))
    return("autodispatch.npregiv.core")
  if (.npRmpi_autodispatch_is_npregivderiv_core(mc = mc))
    return("autodispatch.npregivderiv.core")
  if (.npRmpi_autodispatch_is_npcmstest_core(mc = mc))
    return("autodispatch.npcmstest.core")
  if (.npRmpi_autodispatch_is_npqcmstest_core(mc = mc))
    return("autodispatch.npqcmstest.core")
  if (.npRmpi_autodispatch_is_npdeneqtest_core(mc = mc))
    return("autodispatch.npdeneqtest.core")
  if (.npRmpi_autodispatch_is_npdeptest_core(mc = mc))
    return("autodispatch.npdeptest.core")
  if (.npRmpi_autodispatch_is_npsdeptest_core(mc = mc))
    return("autodispatch.npsdeptest.core")
  if (.npRmpi_autodispatch_is_npsigtest_core(mc = mc))
    return("autodispatch.npsigtest.core")
  if (.npRmpi_autodispatch_is_npsymtest_core(mc = mc))
    return("autodispatch.npsymtest.core")
  if (.npRmpi_autodispatch_is_npunitest_core(mc = mc))
    return("autodispatch.npunitest.core")
  if (.npRmpi_autodispatch_is_density_bw_cv(mc = mc, caller_env = caller_env, call_name = "npudensbw"))
    return("autodispatch.npudensbw.cv")
  if (.npRmpi_autodispatch_is_density_bw_cv(mc = mc, caller_env = caller_env, call_name = "npudistbw"))
    return("autodispatch.npudistbw.cv")
  if (.npRmpi_autodispatch_is_density_bw_cv(mc = mc, caller_env = caller_env, call_name = "npcdensbw"))
    return("autodispatch.npcdensbw.cv")
  if (.npRmpi_autodispatch_is_density_bw_cv(mc = mc, caller_env = caller_env, call_name = "npcdistbw"))
    return("autodispatch.npcdistbw.cv")

  call.name <- .npRmpi_autodispatch_call_name(mc)
  paste0("autodispatch.", gsub("[^A-Za-z0-9_]+", "_", call.name))
}

.npRmpi_spmd_timeout_class_from_opcode <- function(opcode) {
  if (identical(opcode, "autodispatch.npregbw.cv_lllp") ||
      identical(opcode, "autodispatch.npscoefbw.cv_lllp") ||
      identical(opcode, "autodispatch.npplregbw.cv_lllp") ||
      identical(opcode, "autodispatch.npindexbw.core")) {
    return("cv-regression")
  }
  if (identical(opcode, "autodispatch.npudensbw.cv") ||
      identical(opcode, "autodispatch.npudistbw.cv") ||
      identical(opcode, "autodispatch.npcdensbw.cv") ||
      identical(opcode, "autodispatch.npcdistbw.cv")) {
    return("cv-density")
  }
  "default"
}

.npRmpi_autodispatch_preflight <- function(comm = 1L) {
  .npRmpi_abort_if_rmpi_attached(where = "npRmpi auto-dispatch")
  strict <- isTRUE(getOption("npRmpi.autodispatch.strict", TRUE))
  if (!.npRmpi_has_active_slave_pool(comm = comm)) {
    msg <- "npRmpi auto-dispatch requires an active slave pool; call npRmpi.init(...) first"
    if (strict) stop(msg)
    .np_warning(msg)
    return(FALSE)
  }

  TRUE
}

.npRmpi_autodispatch_target_arg_names <-
  c("formula", "data", "bws",
     "dat", "tdat", "edat",
     "xdat", "ydat", "zdat", "txdat", "tydat", "tzdat",
     "exdat", "eydat", "ezdat", "newdata",
     "gydat", "gdat", "wdat", "gdata",
     "data.x", "data.y", "model",
     "y", "z", "w", "x", "zeval", "weval", "xeval", "bw",
     "method", "distribution", "boot.method", "boot.type", "bootstrap", "boot.num", "random.seed",
     "tau", "tau.search", "delta", "delta.bounds", "scale",
     "regtype.pilot", "nomad", "nomad.pilot", "pilot.args",
     "optim.control",
     "bwydat", "joint", "lag.num", "index", "pivot", "density.weighted",
     "weights", "bandwidth.divide", "compute.ocg", "compute.score", "kernel.pow",
     "leave.one.out", "operator", "permutation.operator", "return.kernel.weights",
     "regtype", "basis", "degree", "bernstein.basis",
     "degree.select", "search.engine", "degree.min", "degree.max",
     "degree.start", "degree.restarts", "degree.max.cycles", "degree.verify",
     "nomad.nmulti", "nomad.opts",
     "bwmethod", "bwtype", "bwscaling", "do.full.integral", "ngrid",
     "cvls.quadrature.grid", "cvls.quadrature.extend.factor",
     "cvls.quadrature.points", "cvls.quadrature.ratios",
     "ckertype", "ckerorder", "ckerbound", "ckerlb", "ckerub",
     "cxkertype", "cxkerorder", "cxkerbound", "cxkerlb", "cxkerub",
     "cykertype", "cykerorder", "cykerbound", "cykerlb", "cykerub",
     "ukertype", "okertype", "uxkertype", "oxkertype", "uykertype", "oykertype",
     "nmulti", "bandwidth.compute", "nomad.remin", "powell.remin", "itmax", "ftol", "tol", "small",
     "gradients", "residuals", "se", "errors", "gradient.order",
     "probabilities", "level", "y.eval", "na.action",
     "proper", "proper.method", "proper.control",
     "alpha", "alpha.iter", "alpha.max", "alpha.min", "alpha.tol",
     "constant", "iterate.diff.tol", "iterate.max", "iterate.Tikhonov", "iterate.Tikhonov.num",
     "optim.abstol", "optim.maxattempts", "optim.maxit", "optim.method", "optim.reltol",
     "p", "penalize.iteration", "return.weights.phi", "return.weights.phi.deriv.1",
     "return.weights.phi.deriv.2", "smooth.residuals", "start.from", "starting.values",
     "stop.on.increase", "iterate.break", "n.quasi.inv", "er.quasi.inv",
     "bw.x", "bw.y", "only.optimize.beta",
     "scale.factor.init.lower", "scale.factor.init.upper", "scale.factor.init",
     "scale.factor.search.lower", "betas", "iterate", "maxiter")

.npRmpi_autodispatch_target_args <- function() {
  .npRmpi_autodispatch_target_arg_names
}

.npRmpi_autodispatch_failfast_formula_data <- function(mc, caller_env) {
  invisible(FALSE)
}

.npRmpi_autodispatch_eval_arg <- function(expr, caller_env) {
  .npRmpi_eval_scmd(expr, envir = caller_env)
}

.npRmpi_autodispatch_owner_dot_names <- function(caller_env) {
  dots.call <- tryCatch(
    evalq(substitute(list(...)), envir = caller_env),
    error = function(e) NULL
  )
  if (!is.call(dots.call) || length(dots.call) < 1L ||
      !is.symbol(dots.call[[1L]]) ||
      !identical(as.character(dots.call[[1L]]), "list"))
    return(character(0))

  dots <- as.list(dots.call)[-1L]
  dot.names <- names(dots)
  if (is.null(dot.names))
    dot.names <- rep.int("", length(dots))
  dot.names[is.na(dot.names)] <- ""
  dot.names
}

.npRmpi_autodispatch_owner_dot_value <- function(argname,
                                                  caller_env,
                                                  dot.names) {
  if (!is.character(argname) || length(argname) != 1L ||
      is.na(argname) || !nzchar(argname) || !length(dot.names))
    return(list(found = FALSE, value = NULL))

  idx <- which(dot.names == argname)
  if (!length(idx))
    return(list(found = FALSE, value = NULL))
  if (length(idx) > 1L)
    stop(sprintf("autodispatch cannot materialize duplicate '%s' arguments", argname),
         call. = FALSE)

  dot.expr <- as.name(sprintf("..%d", idx[[1L]]))
  list(
    found = TRUE,
    value = .npRmpi_autodispatch_eval_arg(dot.expr, caller_env = caller_env)
  )
}

.npRmpi_autodispatch_owner_binding <- function(argname, caller_env) {
  if (!is.character(argname) || length(argname) != 1L ||
      is.na(argname) || !nzchar(argname) ||
      !exists(argname, envir = caller_env, inherits = FALSE))
    return(list(found = FALSE, value = NULL))

  list(
    found = TRUE,
    value = get(argname, envir = caller_env, inherits = FALSE)
  )
}

.npRmpi_autodispatch_resolve_owned_arg <- function(expr,
                                                    argname,
                                                    caller_env,
                                                    dot.names) {
  dot.res <- .npRmpi_autodispatch_owner_dot_value(
    argname = argname,
    caller_env = caller_env,
    dot.names = dot.names
  )
  if (isTRUE(dot.res$found))
    return(dot.res$value)

  forwarded <- is.symbol(expr) &&
    grepl("^\\.\\.[0-9]+$", as.character(expr))
  if (forwarded) {
    binding.res <- .npRmpi_autodispatch_owner_binding(
      argname,
      caller_env = caller_env
    )
    if (!isTRUE(binding.res$found))
      stop(sprintf("autodispatch cannot resolve forwarded argument '%s' in its owner frame",
                   argname), call. = FALSE)
    return(binding.res$value)
  }

  eval.res <- tryCatch(
    list(ok = TRUE,
         value = .npRmpi_autodispatch_eval_arg(expr, caller_env = caller_env)),
    error = function(e) list(ok = FALSE, error = e)
  )
  if (isTRUE(eval.res$ok))
    return(eval.res$value)

  binding.res <- .npRmpi_autodispatch_owner_binding(
    argname,
    caller_env = caller_env
  )
  if (isTRUE(binding.res$found))
    return(binding.res$value)
  stop(conditionMessage(eval.res$error), call. = FALSE)
}

.npRmpi_autodispatch_materialize_call <- function(mc, caller_env, comm = 1L) {
  mc <- .npRmpi_autodispatch_expand_dots_call(mc)
  arg.list <- as.list(mc)
  nms <- names(arg.list)
  targets <- .npRmpi_autodispatch_target_args()
  call.base <- sub("\\..*$", "", .npRmpi_autodispatch_call_name(mc))

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
  lease.bindings <- list()
  lease.replacements <- list()
  idx <- 0L
  large.arg.threshold <- .npRmpi_autodispatch_large_arg_threshold_for_call(mc)
  forwarded <- vapply(arg.list, function(expr) {
    is.symbol(expr) && grepl("^\\.\\.[0-9]+$", as.character(expr))
  }, logical(1L))
  owner.dot.names <- if (any(forwarded)) {
    .npRmpi_autodispatch_owner_dot_names(caller_env)
  } else {
    character(0)
  }

  has_data_inputs <- !is.null(nms) && any(nms %in% c("data", "xdat", "ydat", "txdat", "tydat", "zdat"))

  formula.expr <- NULL
  if (!is.null(nms) && any(nms == "formula")) {
    formula.expr <- arg.list[[which(nms == "formula")[1L]]]
  } else if (!is.null(nms) && any(nms == "bws")) {
    bexpr <- arg.list[[which(nms == "bws")[1L]]]
    bval <- tryCatch(
      .npRmpi_autodispatch_resolve_owned_arg(
        expr = bexpr,
        argname = "bws",
        caller_env = caller_env,
        dot.names = owner.dot.names
      ),
      error = function(e) NULL
    )
    if (!is.null(bval) && inherits(bval, "formula"))
      formula.expr <- bexpr
  } else if (length(arg.list) >= 2L) {
    nm2 <- if (!is.null(nms) && length(nms) >= 2L) nms[[2L]] else ""
    if (is.null(nm2) || identical(nm2, "")) {
      fval <- tryCatch(
        .npRmpi_autodispatch_resolve_owned_arg(
          expr = arg.list[[2L]],
          argname = "bws",
          caller_env = caller_env,
          dot.names = owner.dot.names
        ),
        error = function(e) NULL
      )
      if (!is.null(fval) && inherits(fval, "formula"))
        formula.expr <- arg.list[[2L]]
    }
  }

  if (!has_data_inputs && !is.null(formula.expr)) {
    formula.argname <- if (!is.null(nms) && any(nms == "formula")) {
      "formula"
    } else {
      "bws"
    }
    fval <- .npRmpi_autodispatch_resolve_owned_arg(
      expr = formula.expr,
      argname = formula.argname,
      caller_env = caller_env,
      dot.names = owner.dot.names
    )
    vars <- all.vars(fval)
    if (length(vars)) {
      formula.env <- environment(fval)
      if (!is.environment(formula.env))
        formula.env <- caller_env
      dlist <- setNames(lapply(vars, function(v) {
        .npRmpi_autodispatch_eval_arg(as.name(v), caller_env = formula.env)
      }), vars)
      idx <- idx + 1L
      tmp <- sprintf(".__npRmpi_autod_data_%d", idx)
      out[["data"]] <- as.name(tmp)
      tmpnames <- c(tmpnames, tmp)
      stored <- .npRmpi_autodispatch_store_tmp(
        name = tmp,
        value = as.data.frame(dlist, stringsAsFactors = FALSE),
        tmpvals = tmpvals,
        prepublish = prepublish,
        threshold = large.arg.threshold
      )
      tmpvals <- stored$tmpvals
      prepublish <- stored$prepublish
    }
  }

  for (i in seq_along(arg.list)) {
    if (i == 1L) next
    nm <- nms[i]
    if (is.null(nm) || identical(nm, "")) next
    if (!nm %in% targets) next

    expr_i <- arg.list[[i]]
    # S3 forwarding can encode a method formal or named `...` promise as a
    # `..N` placeholder whose numeric position belongs to an outer generic.
    # Materialize from the method frame that owns the promise; never search
    # unrelated dynamic frames for a same-named binding.
    val <- .npRmpi_autodispatch_resolve_owned_arg(
      expr = expr_i,
      argname = nm,
      caller_env = caller_env,
      dot.names = owner.dot.names
    )
    lease <- if (identical(nm, "bws"))
      .npRmpi_lease_resolve(val, consumer = call.base, argument = nm) else NULL
    if (!is.null(lease)) {
      binding <- paste0(".__npRmpi_lease_binding_", .npRmpi_lease_next_id("call"))
      out[[i]] <- as.name(binding)
      lease.bindings[[binding]] <- lease$id
      # Keep the public call reconstructible without adding the bandwidth
      # object back to the worker payload.  The leased binding is private
      # transport state; returned estimator calls must retain the same
      # untagged bandwidth value as the ordinary serialized route.
      lease.replacements[[binding]] <- .npRmpi_autodispatch_untag(val)
      next
    }
    if (is.null(val)) {
      out[i] <- list(NULL)
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
    # npsigtest formulas store master-only symbols in bws$call; forcing the
    # bws object through inline payload avoids worker-side prepublish divergence.
    force.inline <- identical(call.base, "npsigtest") && identical(nm, "bws")
    stored <- .npRmpi_autodispatch_store_tmp(
      name = tmp,
      value = val,
      tmpvals = tmpvals,
      prepublish = prepublish,
      threshold = large.arg.threshold,
      force.inline = force.inline
    )
    tmpvals <- stored$tmpvals
    prepublish <- stored$prepublish
  }

  list(call = out, tmpnames = unique(tmpnames), tmpvals = tmpvals,
       prepublish = prepublish, lease.bindings = lease.bindings,
       lease.replacements = lease.replacements)
}

.npRmpi_lease_run_lifecycle <- function(opcode, payload, comm = 1L, where = opcode) {
  envelope <- .npRmpi_spmd_make_envelope(
    opcode = opcode,
    args_ref = if (is.list(payload)) payload$ids else NULL,
    timeout_class = "default"
  )
  cmd <- substitute({
    run <- get(".npRmpi_spmd_execute_step", envir = asNamespace("npRmpi"), inherits = FALSE)
    run(envelope = ENVELOPE, payload = PAYLOAD, comm = COMM, where = WHERE)
  }, list(ENVELOPE = envelope, PAYLOAD = payload, COMM = comm, WHERE = where))
  .npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE)
}

.npRmpi_lease_drain_pending <- function(comm = 1L) {
  .npRmpi_lease_require_open()
  ids <- ls(.npRmpi_lease_state$pending, all.names = TRUE)
  if (!length(ids))
    return(invisible(character()))
  out <- try(
    .npRmpi_lease_run_lifecycle(
      "autodispatch.lifecycle.retire", list(ids = ids), comm = comm,
      where = "autodispatch pending lease retirement"
    ),
    silent = TRUE
  )
  if (inherits(out, "try-error")) {
    .npRmpi_lease_poison()
    stop(sprintf("autodispatch pending lease retirement failed: %s",
                 paste(as.character(out), collapse = " ")), call. = FALSE)
  }
  present <- ids[vapply(ids, exists, logical(1), envir = .npRmpi_lease_state$pending,
                        inherits = FALSE)]
  if (length(present))
    rm(list = present, envir = .npRmpi_lease_state$pending)
  invisible(ids)
}

.npRmpi_lease_session_drain <- function(comm = 1L) {
  if (!identical(.npRmpi_lease_state$status, "open"))
    return(invisible(FALSE))
  out <- try(
    .npRmpi_lease_run_lifecycle(
      "autodispatch.lifecycle.session_drain", list(ids = NULL), comm = comm,
      where = "autodispatch session lease drain"
    ),
    silent = TRUE
  )
  if (inherits(out, "try-error")) {
    .npRmpi_lease_poison()
    stop(sprintf("autodispatch session lease drain failed: %s",
                 paste(as.character(out), collapse = " ")), call. = FALSE)
  }
  invisible(TRUE)
}

.npRmpi_autodispatch_cleanup <- function(tmpnames, comm = 1L) {
  if (!length(tmpnames))
    return(invisible(TRUE))
  tmpnames <- unique(as.character(tmpnames))
  registry.ids <- ls(.npRmpi_lease_state$registry, all.names = TRUE)
  scoped.ids <- registry.ids[vapply(registry.ids, function(id) {
    meta <- get(id, envir = .npRmpi_lease_state$registry, inherits = FALSE)
    is.list(meta) && identical(meta$kind, "scoped_internal") &&
      is.character(meta$remote_name) && meta$remote_name %in% tmpnames
  }, logical(1))]
  scoped.names <- if (length(scoped.ids)) vapply(scoped.ids, function(id) {
    get(id, envir = .npRmpi_lease_state$registry, inherits = FALSE)$remote_name
  }, character(1)) else character()
  if (length(scoped.ids))
    .npRmpi_lease_run_lifecycle(
      "autodispatch.lifecycle.retire", list(ids = scoped.ids), comm = comm,
      where = "autodispatch scoped context retirement"
    )
  ordinary <- setdiff(tmpnames, scoped.names)
  if (length(ordinary)) {
    .npRmpi_rm_existing(ordinary, envir = .GlobalEnv)
    cmd.rm <- substitute(get(".npRmpi_rm_existing", envir = asNamespace("npRmpi"), inherits = FALSE)(TMPS, envir = .GlobalEnv),
                         list(TMPS = ordinary))
    .npRmpi_bcast_cmd_expr(cmd.rm, comm = comm, caller.execute = FALSE)
  }
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

  # Single-bracket list assignment preserves explicit NULL entries and arity.
  if (is.call(x)) {
    xlist <- as.list(x)
    for (i in seq_along(xlist)) {
      xi <- xlist[[i]]
      if (.npRmpi_is_missing_call_arg(xi))
        next
      if (is.null(xi)) {
        xlist[i] <- list(NULL)
        next
      }
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
      if (is.null(xi)) {
        xlist[i] <- list(NULL)
        next
      }
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

.npRmpi_autodispatch_sanitize_object <- function(x, tmpvals) {
  if (is.list(x))
    return(.npRmpi_autodispatch_untag(
      .npRmpi_autodispatch_replace_tmp_calls(x, tmpvals = tmpvals)
    ))

  if (is.call(x))
    return(.npRmpi_autodispatch_replace_tmps(x, tmpvals = tmpvals))

  .npRmpi_autodispatch_untag(x)
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

.npRmpi_autodispatch_tag_result <- function(x, mode = "auto", remote = NULL,
                                            publication = NULL) {
  if (!is.null(x)) {
    attr(x, "npRmpi.dispatch.mode") <- mode
    if (is.list(publication) && !identical(publication$kind, "none")) {
      if (!exists(publication$id, envir = .npRmpi_lease_state$registry, inherits = FALSE))
        stop("autodispatch publication metadata is unavailable after commit", call. = FALSE)
      meta <- get(publication$id, envir = .npRmpi_lease_state$registry, inherits = FALSE)
      if (!identical(meta$state, "live"))
        stop("autodispatch publication is not live after commit", call. = FALSE)
      current.fingerprint <- .npRmpi_autodispatch_fingerprint(
        x, policy = publication$fingerprint_policy
      )
      if (!identical(current.fingerprint, meta$fingerprint)) {
        .npRmpi_lease_poison()
        stop("coordinator public-result fingerprint contradicts committed ranks", call. = FALSE)
      }
      attr(x, "npRmpi.autodispatch.remote") <- publication$remote_name
      attr(x, "npRmpi.autodispatch.fingerprint") <- meta$fingerprint
      if (identical(publication$kind, "leased_public"))
        x <- .npRmpi_lease_attach_handle(x, publication, meta$fingerprint)
      return(x)
    }
    if (is.character(remote) && length(remote) == 1L && nzchar(remote)) {
      attr(x, "npRmpi.autodispatch.remote") <- remote
      attr(x, "npRmpi.autodispatch.fingerprint") <- .npRmpi_autodispatch_fingerprint(x)
    }
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
  attr(x, "npRmpi.autodispatch.fingerprint") <- NULL
  attr(x, "npRmpi.autodispatch.lease") <- NULL
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
  mc <- .npRmpi_autodispatch_expand_dots_call(mc)
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
    comm.time <- as.double(comm.elapsed)
    if (!is.finite(comm.time) || comm.time < 0)
      comm.time <- 0.0
    compute <- max(0.0, wall - comm.time)
    denom <- comm.time + compute
    ratio <- if (denom > 0) min(1.0, max(0.0, comm.time / denom)) else NA_real_
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
      comm_elapsed_sec = comm.time,
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

  .npRmpi_autodispatch_validate_option_snapshot(.npRmpi_autodispatch_option_snapshot())
  .npRmpi_autodispatch_failfast_formula_data(mc, caller_env = caller_env)

  if (!.npRmpi_autodispatch_preflight(comm = comm))
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  rec.comm(.npRmpi_lease_drain_pending(comm = comm),
           note = "autodispatch.lease.drain")
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
  opcode <- .npRmpi_spmd_opcode_from_call(mc = mc, caller_env = caller_env)
  publication <- .npRmpi_lease_publication_plan(mc)
  publication$opcode <- opcode
  timeout.class <- .npRmpi_spmd_timeout_class_from_opcode(opcode)
  envelope <- .npRmpi_spmd_make_envelope(
    opcode = opcode,
    args_ref = prepared$tmpnames,
    timeout_class = timeout.class
  )
  payload <- list(
    call = prepared$call,
    tmpvals = prepared$tmpvals,
    tmpnames = prepared$tmpnames,
    prepublish.names = names(prepared$prepublish),
    opt.keys = opt.keys,
    opt.vals = opt.vals,
    opt.verify = opt.verify,
    lease.bindings = prepared$lease.bindings,
    publication = publication
  )

  cmd <- substitute({
    exec.step <- get(".npRmpi_spmd_execute_step", envir = asNamespace("npRmpi"), inherits = FALSE)
    ans <- exec.step(
      envelope = ENVELOPE,
      payload = PAYLOAD,
      comm = COMM,
      where = WHERE
    )
    ans$result
  }, list(ENVELOPE = envelope,
          PAYLOAD = payload,
          COMM = comm,
          WHERE = .npRmpi_autodispatch_call_name(mc)))

  eval.out <- try(
    rec.step(.npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE),
             note = "mpi.bcast.cmd.execute"),
    silent = TRUE
  )
  if (inherits(eval.out, "try-error")) {
    if (!identical(publication$kind, "none")) {
      rolled <- try(
        rec.comm(.npRmpi_lease_run_lifecycle(
          "autodispatch.lifecycle.rollback", list(ids = publication$id),
          comm = comm, where = "autodispatch publication rollback"
        ), note = "autodispatch.lease.rollback"),
        silent = TRUE
      )
      if (inherits(rolled, "try-error"))
        .npRmpi_lease_poison()
    }
    stop(paste(as.character(eval.out), collapse = " "), call. = FALSE)
  }
  result <- eval.out
  if (identical(publication$kind, "scoped_internal")) {
    committed <- try(
      rec.comm(.npRmpi_lease_run_lifecycle(
        "autodispatch.lifecycle.commit", list(plan = publication),
        comm = comm, where = "autodispatch publication commit"
      ), note = "autodispatch.lease.commit"),
      silent = TRUE
    )
    if (inherits(committed, "try-error")) {
      rolled <- try(
        rec.comm(.npRmpi_lease_run_lifecycle(
          "autodispatch.lifecycle.rollback", list(ids = publication$id),
          comm = comm, where = "autodispatch failed-commit rollback"
        ), note = "autodispatch.lease.rollback"),
        silent = TRUE
      )
      if (inherits(rolled, "try-error"))
        .npRmpi_lease_poison()
      stop(paste(as.character(committed), collapse = " "), call. = FALSE)
    }
  }
  if (length(prepared$tmpnames))
    cleaned <- TRUE
  tmpreplace <- c(prepared$tmpvals, prepared$prepublish,
                  prepared$lease.replacements)
  rec <- make.rec()

  if (is.list(result))
    return(.npRmpi_autodispatch_attach_timing(
      .npRmpi_autodispatch_tag_result(
        .npRmpi_autodispatch_sanitize_object(result, tmpvals = tmpreplace),
        mode = "auto", publication = publication
      ),
      rec
    ))

  if (is.call(result))
    return(.npRmpi_autodispatch_tag_result(
      .npRmpi_autodispatch_replace_tmps(result, tmpvals = tmpreplace),
      mode = "auto", publication = publication
    ))

  .npRmpi_autodispatch_tag_result(result, mode = "auto", publication = publication)
}

.npRmpi_autodispatch_call <- function(mc, caller_env = parent.frame(), comm = 1L) {
  .npRmpi_warn_pkg_conflict_once()
  .npRmpi_warn_rmpi_conflict_once()
  if (missing(caller_env) || !is.environment(caller_env))
    caller_env <- parent.frame()
  if (!.npRmpi_autodispatch_active())
    return(.npRmpi_eval_without_dispatch(mc, caller_env))

  .npRmpi_distributed_call_impl(mc = mc, caller_env = caller_env, comm = comm, warn_nested = TRUE)
}

.npRmpi_manual_distributed_call <- function(mc, caller_env = parent.frame(), comm = 1L) {
  .npRmpi_warn_pkg_conflict_once()
  .npRmpi_warn_rmpi_conflict_once()
  if (missing(caller_env) || !is.environment(caller_env))
    caller_env <- parent.frame()
  .npRmpi_distributed_call_impl(mc = mc, caller_env = caller_env, comm = comm, warn_nested = FALSE)
}

.npRmpi_bootstrap_compute_payload <- function(payload, comm = 1L) {
  .npRmpi_require_active_slave_pool(
    comm = comm,
    where = "bootstrap payload computation"
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
