# One fixed-coordinate density-ratio uncertainty owner, shared by fitting,
# prediction and asymptotic plot payloads. No work is done for se = FALSE.
.npcopula_density_se_local <- function(bws, data, xgrid, progress.context = NULL) {
  if (!identical(bws[["type", exact = TRUE]], "fixed"))
    stop("internal copula influence requires fixed bandwidths", call. = FALSE)
  n <- nrow(data)
  m <- nrow(xgrid)
  p <- length(bws[["xnames", exact = TRUE]])
  if (n < 1L || p < 1L)
    stop("internal copula influence requires nonempty training data", call. = FALSE)
  if (m == 0L)
    return(numeric(0L))

  # Account for all p resident marginal tiles, not just a single matrix.
  tile <- max(1L, min(64L, floor(1048576 / (as.double(n) * p))))
  marginal.bws <- lapply(seq_len(p), function(j)
    .npcopula_marginal_bw(bws, data, j, target = "density", kbandwidth = TRUE))
  training <- lapply(seq_len(p), function(j)
    .npcopula_marginal_data(bws, data, j))
  evaluation <- lapply(seq_len(p), function(j)
    .npcopula_marginal_eval_data(bws, xgrid, j))
  out <- numeric(m)
  for (start in seq.int(1L, m, by = tile)) {
    rows <- seq.int(start, min(m, start + tile - 1L))
    weights <- lapply(seq_len(p), function(j) {
      npksum(bws = marginal.bws[[j]], txdat = training[[j]],
             exdat = evaluation[[j]][rows, , drop = FALSE],
             operator = "normal", bandwidth.divide = FALSE,
             return.kernel.weights = TRUE)[["kw", exact = TRUE]]
    })
    out[rows] <- .Call("C_np_copula_density_se", weights, PACKAGE = "npRmpi")
    weights <- NULL
    if (!is.null(progress.context))
      progress.context$state <- .np_progress_step_at(
        progress.context$state, now = .np_progress_now(),
        done = if (isTRUE(progress.context$use.bootstrap.done)) max(rows) else
          progress.context$done,
        detail = sprintf("asymptotic standard errors (%d/%d rows)", max(rows), m))
  }
  out
}

.npcopula_density_se_fixed <- function(bws, data, xgrid, progress.context = NULL) {
  if (is.null(progress.context)) {
    progress.context <- new.env(parent = emptyenv())
    progress.context$state <- .np_progress_begin(
      "Copula density standard errors", total = nrow(xgrid), surface = "copula")
    progress.context$use.bootstrap.done <- TRUE
    on.exit(.np_progress_end(progress.context$state), add = TRUE)
  }
  out <- .npRmpi_copula_density_se(bws, data, xgrid, progress.context)
  if (anyNA(out))
    .np_warning(paste0(
      "Copula density standard errors are unavailable for ",
      sum(is.na(out)), " evaluation row(s): a marginal density is zero ",
      "or the ratio influence is nonfinite. These standard errors are NA."
    ), call. = FALSE)
  out
}

# Each worker owns complete query rows; marginal weight tiles never travel.
.npRmpi_copula_density_se <- function(bws, data, xgrid, progress.context = NULL) {
  .npRmpi_require_active_slave_pool(where = "npcopula standard errors")
  if (isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
    return(.npcopula_density_se_local(bws, data, xgrid, progress.context))
  m <- nrow(xgrid)
  if (.npRmpi_autodispatch_called_from_bcast()) {
    rank <- as.integer(mpi.comm.rank(comm = 1L))
    size <- as.integer(mpi.comm.size(comm = 1L))
    base <- m %/% size
    extra <- m %% size
    count <- base + as.integer(rank < extra)
    start <- rank * base + min(rank, extra) + 1L
    rows <- if (count > 0L) seq.int(start, length.out = count) else integer(0L)
    local <- .npRmpi_capture_local_work(.npRmpi_with_local_regression(
      .npcopula_density_se_local(bws, data, xgrid[rows, , drop = FALSE], progress.context)))
    parts <- mpi.allgather.Robj(list(list(rows = rows, value = local)), comm = 1L)
    out <- numeric(m)
    seen <- integer(m)
    for (wrapped in parts) {
      part <- if (is.list(wrapped) && length(wrapped) == 1L &&
                  is.list(wrapped[[1L]])) wrapped[[1L]] else wrapped
      if (!is.list(part) || is.null(part[["rows"]]) || is.null(part[["value"]]))
        stop("copula SE gathered a malformed result", call. = FALSE)
      if (inherits(part[["value"]], "npRmpi_local_failure"))
        .npRmpi_raise_completed_failure(part[["value"]])
      ids <- part[["rows"]]
      if (!is.integer(ids) || anyNA(ids) || any(ids < 1L | ids > m) ||
          anyDuplicated(ids) || !is.double(part[["value"]]) ||
          length(part[["value"]]) != length(ids))
        stop("copula SE gathered invalid row ownership", call. = FALSE)
      out[ids] <- part[["value"]]
      seen[ids] <- seen[ids] + 1L
    }
    if (any(seen != 1L))
      stop("copula SE row ownership was incomplete or duplicated", call. = FALSE)
    return(out)
  }

  # The established task runner also serves hat-operator evaluation rows.
  # No random seeds or bootstrap data are generated for these deterministic tasks.
  tasks <- .npRmpi_bootstrap_chunk_tasks(
    B = m, chunk.size = 64L, with.seeds = FALSE)
  snapshot <- .npRmpi_autodispatch_option_snapshot()
  worker <- function(task, bws, data, xgrid, snapshot) {
    old <- options(snapshot)
    on.exit(options(old), add = TRUE)
    rows <- seq.int(task$start, length.out = task$bsz)
    values <- .npRmpi_with_local_regression(
      .npcopula_density_se_local(bws, data, xgrid[rows, , drop = FALSE]))
    matrix(values, ncol = 1L)
  }
  # The one-task mpi.apply branch draws a transport tag. This deterministic
  # uncertainty calculation must not advance (or create) the user's RNG state.
  seed.state <- list(
    exists.seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    save.seed = get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)
  result <- .npRmpi_bootstrap_run_fanout(
    tasks = tasks, worker = worker, ncol.out = 1L,
    what = "npcopula-density-se", progress.label = "Copula density standard errors",
    profile.where = "npcopula-density-se", comm = 1L,
    prefer.local.single_worker = FALSE, master_local_chunk = TRUE,
    bws = bws, data = data, xgrid = xgrid, snapshot = snapshot,
    progress.context = progress.context)
  as.numeric(result[, 1L])
}
