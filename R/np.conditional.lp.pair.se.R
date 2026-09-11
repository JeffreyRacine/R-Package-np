# Complete endpoint pairs are rank-owned. Kernels, bases and influence vectors
# remain local; only the six-column result crosses the completion boundary.
.np_conditional_lp_pair_se <- function(bws, txdat, tydat, upper, lower,
                                       upper.y, lower.y = upper.y, cdf = FALSE,
                                       quantile = FALSE) {
  .npRmpi_require_active_slave_pool(where = "conditional LP contrast standard errors")
  if (isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
    return(.np_conditional_lp_pair_se_local(bws, txdat, tydat, upper, lower,
                                           upper.y, lower.y, cdf, quantile))
  m <- nrow(upper)
  if (.npRmpi_autodispatch_called_from_bcast()) {
    rank <- as.integer(mpi.comm.rank(comm = 1L))
    size <- as.integer(mpi.comm.size(comm = 1L))
    base <- m %/% size
    extra <- m %% size
    count <- base + as.integer(rank < extra)
    first <- rank * base + min(rank, extra) + 1L
    rows <- if (count) seq.int(first, length.out = count) else integer(0L)
    value <- .npRmpi_capture_local_work(.npRmpi_with_local_regression(
      .np_conditional_lp_pair_se_local(bws, txdat, tydat,
        upper[rows, , drop = FALSE], lower[rows, , drop = FALSE],
        upper.y[rows, , drop = FALSE], lower.y[rows, , drop = FALSE], cdf, quantile)))
    parts <- mpi.allgather.Robj(list(list(rows = rows, value = value)), comm = 1L)
    out <- matrix(NA_real_, m, 6L)
    seen <- integer(m)
    for (wrapped in parts) {
      part <- if (is.list(wrapped) && length(wrapped) == 1L &&
                  is.list(wrapped[[1L]])) wrapped[[1L]] else wrapped
      if (!is.list(part) || is.null(part[["rows", exact = TRUE]]) ||
          is.null(part[["value", exact = TRUE]]))
        stop("conditional LP covariance gathered a malformed result", call. = FALSE)
      ids <- part[["rows", exact = TRUE]]
      values <- part[["value", exact = TRUE]]
      if (inherits(values, "npRmpi_local_failure"))
        .npRmpi_raise_completed_failure(values)
      if (!is.integer(ids) || anyNA(ids) || any(ids < 1L | ids > m) ||
          anyDuplicated(ids) || !is.double(values) ||
          !identical(dim(values), c(length(ids), 6L)))
        stop("conditional LP covariance gathered invalid row ownership", call. = FALSE)
      out[ids, ] <- values
      seen[ids] <- seen[ids] + 1L
    }
    if (any(seen != 1L))
      stop("conditional LP covariance row ownership was incomplete or duplicated", call. = FALSE)
    return(out)
  }
  if (!m) return(matrix(numeric(0L), 0L, 6L))
  tasks <- .npRmpi_bootstrap_chunk_tasks(B = m, chunk.size = 8L, with.seeds = FALSE)
  snapshot <- .npRmpi_autodispatch_option_snapshot()
  worker <- function(task, bws, txdat, tydat, upper, lower, upper.y, lower.y,
                     cdf, quantile, snapshot) {
    old <- options(snapshot)
    on.exit(options(old), add = TRUE)
    rows <- seq.int(task$start, length.out = task$bsz)
    .npRmpi_with_local_regression(.np_conditional_lp_pair_se_local(
      bws, txdat, tydat, upper[rows, , drop = FALSE], lower[rows, , drop = FALSE],
      upper.y[rows, , drop = FALSE], lower.y[rows, , drop = FALSE], cdf, quantile))
  }
  seed.state <- list(
    exists.seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    save.seed = get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)
  .npRmpi_bootstrap_run_fanout(tasks = tasks, worker = worker, ncol.out = 6L,
    what = "conditional-lp-pair-se", progress.label = "Categorical contrast standard errors",
    profile.where = "conditional-lp-pair-se", comm = 1L,
    prefer.local.single_worker = FALSE, master_local_chunk = TRUE,
    bws = bws, txdat = txdat, tydat = tydat, upper = upper, lower = lower,
    upper.y = upper.y, lower.y = lower.y, cdf = cdf, quantile = quantile,
    snapshot = snapshot)
}

# Requested-only paired empirical-mass covariance for the realized conditional
# LP fit. Endpoint construction and published point results remain with callers.
.np_conditional_lp_pair_se_local <- function(bws, txdat, tydat, upper, lower,
                                       upper.y, lower.y = upper.y, cdf = FALSE,
                                       quantile = FALSE) {
  m <- nrow(upper)
  if (nrow(lower) != m || nrow(upper.y) != m || nrow(lower.y) != m)
    stop("internal conditional LP covariance endpoint mismatch", call. = FALSE)
  out <- matrix(NA_real_, m, 6L)
  if (!m) return(out)
  xbw <- .npcdhat_make_xbw(bws, txdat)
  ybw <- .npcdhat_make_ybw(bws, tydat)
  spec <- npConditionalRegEngineSpec(xbw, where = "conditional LP covariance",
                                    ncon.field = "ncon")
  z <- as.matrix(W.lp(txdat[, xbw$icon, drop = FALSE],
    degree = spec$degree.engine, basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.engine))
  storage.mode(z) <- "double"
  # At most eight endpoint pairs per tile; additionally cap the kernel rows
  # near one MiB each. No observation-by-all-evaluation result is retained.
  tile <- max(1L, min(8L, floor(131072 / nrow(txdat))))
  for (first in seq.int(1L, m, by = tile)) {
    rows <- seq.int(first, min(m, first + tile - 1L))
    ex <- rbind(upper[rows, , drop = FALSE], lower[rows, , drop = FALSE])
    ey <- rbind(upper.y[rows, , drop = FALSE], lower.y[rows, , drop = FALSE])
    d <- as.matrix(W.lp(txdat[, xbw$icon, drop = FALSE],
      exdat = ex[, xbw$icon, drop = FALSE], degree = spec$degree.engine,
      basis = spec$basis.engine, bernstein.basis = spec$bernstein.engine))
    storage.mode(d) <- "double"
    w <- .np_kernel_weights_direct(xbw, txdat, exdat = ex,
      bandwidth.divide = TRUE,
      int.do.tree = .npreg_fit_tree_code(xbw, xbw$ncon, xbw$nuno + xbw$nord))
    r <- t(.npcdhat_make_kernel_matrix(ybw, tydat, ey,
      operator = if (cdf) "integral" else "normal"))
    slope <- if (quantile) t(.npcdhat_make_kernel_matrix(
      ybw, tydat, ey, operator = "normal")) else NULL
    out[rows, ] <- .Call("C_np_conditional_lp_pair_se", z, d,
      as.matrix(w), r, slope, PACKAGE = "npRmpi")
  }
  out
}
