# Package-owned loaded identity. This does not certify subsequent tracing or
# namespace mutation, nor hand-patched binaries outside ordinary installation.
.npRmpi_worker_identity_state <- new.env(parent = emptyenv())

.npRmpi_worker_identity_load <- function(ns) {
  # Force the installed lazy-load bindings while their source database still
  # corresponds to the loaded namespace. No user callback is evaluated here.
  for (name in ls(ns, all.names = TRUE)) {
    if (bindingIsActive(name, ns))
      stop("npRmpi cannot bind source identity for an active namespace binding", call. = FALSE)
    get(name, envir = ns, inherits = FALSE)
  }
  source <- get(".npRmpi_worker_source", envir = ns, inherits = FALSE)
  native <- .Call(C_np_worker_native_identity)
  .npRmpi_worker_identity_state$record <- list(
    ok = identical(source[["native"]], native), source = source,
    native = native, library = getNamespaceInfo(ns, "path"),
    platform = R.version$platform, r.version = as.character(getRversion()))
  invisible(NULL)
}

.npRmpi_worker_identity_record <- function() {
  .npRmpi_worker_identity_state$record
}

.npRmpi_worker_identity_problems <- function(records) {
  valid <- function(x) {
    scalar <- function(v) is.character(v) && length(v) == 1L && !is.na(v) && nzchar(v)
    hash <- function(v) scalar(v) && grepl("^[0-9a-f]{32}$", v)
    is.list(x) && isTRUE(x[["ok"]]) && is.list(x[["source"]]) &&
      identical(x[["source"]][["schema"]], 1L) &&
      scalar(x[["source"]][["version"]]) && hash(x[["source"]][["r"]]) &&
      hash(x[["native"]]) && identical(x[["source"]][["native"]], x[["native"]])
  }
  good <- vapply(records, valid, logical(1))
  problems <- rep("", length(records))
  problems[!good] <- "missing, malformed, or inconsistent loaded source identity"
  if (length(records) && good[[1L]]) {
    same <- vapply(seq_along(records), function(i)
      good[[i]] && identical(records[[i]][["source"]], records[[1L]][["source"]]), logical(1))
    problems[good & !same] <- "loaded package source differs from coordinator"
  } else problems[good] <- "coordinator has no valid loaded source identity"
  problems
}

.npRmpi_worker_identity_collect <- function(comm) {
  # Literal inspection uses only primitives present on older compatible workers.
  # Never throw a record error before every rank reaches the existing collective.
  command <- substitute(base::local({
    result <- tryCatch({
      ns <- asNamespace("npRmpi")
      getter <- get0(".npRmpi_worker_identity_record", envir = ns,
                     mode = "function", inherits = FALSE)
      if (is.null(getter)) list(ok = FALSE) else getter()
    }, error = function(e) list(ok = FALSE))
    bytes <- tryCatch(serialize(result, NULL, version = 2L),
                      error = function(e) serialize(list(ok = FALSE), NULL, version = 2L))
    if (length(bytes) > 4096L)
      bytes <- serialize(list(ok = FALSE), NULL, version = 2L)
    buffer <- raw(4096L)
    buffer[seq_along(bytes)] <- bytes
    get("mpi.gather", asNamespace("npRmpi"))(
      buffer, type = 4L,
      rdata = if (get("mpi.comm.rank", asNamespace("npRmpi"))(COMM) == 0L)
        raw(4096L * get("mpi.comm.size", asNamespace("npRmpi"))(COMM)) else raw(1L),
      root = 0L, comm = COMM)
  }), list(COMM = comm))
  bytes <- .npRmpi_bcast_cmd_expr(command, comm = comm, caller.execute = TRUE)
  count <- mpi.comm.size(comm)
  if (!is.raw(bytes) || length(bytes) != count * 4096L)
    stop("npRmpi received an invalid worker identity record buffer", call. = FALSE)
  lapply(seq_len(count), function(i) {
    segment <- bytes[seq.int((i - 1L) * 4096L + 1L, i * 4096L)]
    tryCatch(unserialize(segment), error = function(e) list(ok = FALSE))
  })
}

.npRmpi_worker_identity_admit <- function(comm) {
  records <- .npRmpi_worker_identity_collect(comm)
  problems <- .npRmpi_worker_identity_problems(records)
  bad <- which(nzchar(problems))
  if (length(bad)) {
    # Reuse the existing closed-session disposition; do not change transport,
    # worker loops, or cleanup. Workers have already finished the collective.
    .npRmpi_lease_poison()
    stop(paste0("npRmpi worker code identity mismatch: ",
                paste(sprintf("rank %d: %s", bad - 1L, problems[bad]), collapse = "; "),
                ". Close this pool and restart all ranks with matching package sources."),
         call. = FALSE)
  }
  invisible(TRUE)
}
