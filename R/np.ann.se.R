# Same-sample random-radius influence for unconditional ANN estimates.
# This owner is entered only when uncertainty was requested. Point estimates,
# bandwidths, kernel evaluation and computational options are never changed.
.np_ann_se_domain <- function(bws, n, density) {
  if (n < 4L)
    return("the sample is too small for a two-sided NN spacing pilot")
  if (!identical(bws[["ckerbound", exact = TRUE]], "none"))
    return("ANN uncertainty with finite kernel bounds is not yet implemented")
  kernel <- bws[["ckertype", exact = TRUE]]
  if (!kernel %in% c("gaussian", "epanechnikov", "uniform") ||
      (density && identical(kernel, "uniform")))
    return("the kernel's ANN support-boundary influence is not yet implemented")
  k <- as.double(unlist(bws[["bandwidth", exact = TRUE]], use.names = FALSE))[bws[["icon", exact = TRUE]]]
  if (!length(k) || any(!is.finite(k)) || any(k != floor(k)) ||
      any(k < 2 | k > n - 3L))
    return("the NN rank has no regular two-sided interior spacing pilot")
  NULL
}

.np_ann_se_geometry <- function(bws, data) {
  x <- as.matrix(data[, bws[["icon", exact = TRUE]], drop = FALSE])
  storage.mode(x) <- "double"
  k <- as.double(unlist(bws[["bandwidth", exact = TRUE]], use.names = FALSE))[bws[["icon", exact = TRUE]]]
  orders <- lapply(seq_len(ncol(x)), function(j) order(x[, j]))
  geometry <- lapply(seq_len(ncol(x)), function(j)
    .Call("C_np_ann_geometry", x[orders[[j]], j], as.integer(k[j]),
          PACKAGE = "npRmpi"))
  list(x = x, order = orders, geometry = geometry,
       valid = all(vapply(geometry, function(g) all(g[["status"]] == 0L), logical(1L))))
}

.np_ann_unconditional_se_local <- function(bws, data, evaluation, density,
                                            progress.context = NULL) {
  n <- nrow(data)
  m <- nrow(evaluation)
  if (!m) return(numeric(0L))
  state <- .np_ann_se_geometry(bws, data)
  if (!state$valid) return(rep.int(NA_real_, m))
  p <- ncol(state$x)
  e <- as.matrix(evaluation[, bws[["icon", exact = TRUE]], drop = FALSE])
  storage.mode(e) <- "double"
  kbw <- kbandwidth(bws)
  # Bound the entire resident normal + permutation tensor, not each separately.
  tile <- max(1L, min(64L, floor(1048576 / (as.double(n) * (1L + p)))))
  out <- numeric(m)
  for (start in seq.int(1L, m, by = tile)) {
    rows <- seq.int(start, min(m, start + tile - 1L))
    weights <- npksum(bws = kbw, txdat = data,
      exdat = evaluation[rows, , drop = FALSE],
      operator = if (density) "normal" else "integral",
      permutation.operator = if (density) "derivative" else "normal",
      bandwidth.divide = TRUE, return.kernel.weights = TRUE,
      return.derivative.kernel.weights = TRUE,
      .np.internal.bandwidth.divide.weights = TRUE)
    out[rows] <- .Call("C_np_ann_variance", state$geometry, state$order,
      state$x, e[rows, , drop = FALSE], weights[["kw", exact = TRUE]],
      weights[["p.kw", exact = TRUE]], density, PACKAGE = "npRmpi")
    weights <- NULL
    if (!is.null(progress.context))
      progress.context$state <- .np_progress_step_at(
        progress.context$state, now = .np_progress_now(), done = max(rows),
        detail = sprintf("asymptotic standard errors (%d/%d rows)", max(rows), m))
  }
  out
}

.np_ann_unconditional_se <- function(bws, data, evaluation, density) {
  reason <- .np_ann_se_domain(bws, nrow(data), density)
  if (!is.null(reason)) {
    .np_warning(paste0("Adaptive-NN standard errors are unavailable: ", reason,
                      ". Point estimates are retained; standard errors are NA."),
                call. = FALSE)
    return(rep.int(NA_real_, nrow(evaluation)))
  }
  progress.context <- new.env(parent = emptyenv())
  progress.context$state <- .np_progress_begin(
    "Adaptive-NN standard errors", total = nrow(evaluation), surface = "fit")
  on.exit(.np_progress_end(progress.context$state), add = TRUE)
  out <- .npRmpi_ann_unconditional_se(bws, data, evaluation, density, progress.context)
  if (anyNA(out))
    .np_warning(paste0("Adaptive-NN standard errors are unavailable for ",
      sum(is.na(out)), " evaluation row(s): a NN spacing pilot is degenerate ",
      "or the estimated influence is nonfinite. Point estimates are retained; ",
      "these standard errors are NA."), call. = FALSE)
  out
}

# Existing query-row task/completion protocol; no donor-weight tensor is sent.
.npRmpi_ann_unconditional_se <- function(bws, data, evaluation, density,
                                         progress.context = NULL) {
  .npRmpi_require_active_slave_pool(where = "adaptive-NN standard errors")
  if (isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
    return(.np_ann_unconditional_se_local(bws, data, evaluation, density, progress.context))
  m <- nrow(evaluation)
  if (.npRmpi_autodispatch_called_from_bcast()) {
    rank <- as.integer(mpi.comm.rank(comm = 1L))
    size <- as.integer(mpi.comm.size(comm = 1L))
    base <- m %/% size
    extra <- m %% size
    count <- base + as.integer(rank < extra)
    start <- rank * base + min(rank, extra) + 1L
    rows <- if (count > 0L) seq.int(start, length.out = count) else integer(0L)
    local <- .npRmpi_capture_local_work(.npRmpi_with_local_regression(
      .np_ann_unconditional_se_local(bws, data,
        evaluation[rows, , drop = FALSE], density, progress.context)))
    parts <- mpi.allgather.Robj(list(list(rows = rows, value = local)), comm = 1L)
    out <- numeric(m)
    seen <- integer(m)
    for (wrapped in parts) {
      part <- if (is.list(wrapped) && length(wrapped) == 1L &&
                  is.list(wrapped[[1L]])) wrapped[[1L]] else wrapped
      if (!is.list(part) || is.null(part[["rows"]]) || is.null(part[["value"]]))
        stop("ANN SE gathered a malformed result", call. = FALSE)
      if (inherits(part[["value"]], "npRmpi_local_failure"))
        .npRmpi_raise_completed_failure(part[["value"]])
      ids <- part[["rows"]]
      if (!is.integer(ids) || anyNA(ids) || any(ids < 1L | ids > m) ||
          anyDuplicated(ids) || !is.double(part[["value"]]) ||
          length(part[["value"]]) != length(ids))
        stop("ANN SE gathered invalid row ownership", call. = FALSE)
      out[ids] <- part[["value"]]
      seen[ids] <- seen[ids] + 1L
    }
    if (any(seen != 1L))
      stop("ANN SE row ownership was incomplete or duplicated", call. = FALSE)
    return(out)
  }
  tasks <- .npRmpi_bootstrap_chunk_tasks(B = m, chunk.size = 64L, with.seeds = FALSE)
  snapshot <- .npRmpi_autodispatch_option_snapshot()
  worker <- function(task, bws, data, evaluation, density, snapshot) {
    old <- options(snapshot)
    on.exit(options(old), add = TRUE)
    rows <- seq.int(task$start, length.out = task$bsz)
    values <- .npRmpi_with_local_regression(
      .np_ann_unconditional_se_local(bws, data,
        evaluation[rows, , drop = FALSE], density))
    matrix(values, ncol = 1L)
  }
  seed.state <- list(
    exists.seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    save.seed = get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)
  result <- .npRmpi_bootstrap_run_fanout(
    tasks = tasks, worker = worker, ncol.out = 1L,
    what = "ann-unconditional-se", progress.label = "Adaptive-NN standard errors",
    profile.where = "ann-unconditional-se", comm = 1L,
    prefer.local.single_worker = FALSE, master_local_chunk = TRUE,
    bws = bws, data = data, evaluation = evaluation, density = density,
    snapshot = snapshot, progress.context = progress.context)
  as.numeric(result[, 1L])
}
