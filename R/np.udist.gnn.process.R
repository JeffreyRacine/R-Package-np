# Exact resampling of the GNN statistic; no bootstrap studentizing pilot.
.np_udist_gnn_chunk <- function(prepared, base, first, last, B, law, blocklen) {
  prepared$zero.limit <- law %in% c("inid", "fixed")
  counts <- .np_bias_process_innovations(prepared$n, first, last, B, law, blocklen) + 1
  out <- matrix(0, length(base), last - first + 1L)
  for (j in seq_len(ncol(counts))) {
    value <- tryCatch(.np_udist_gnn_center(prepared, counts[, j])$center,
      error = function(e) stop(sprintf("exact GNN bootstrap draw %d unavailable: %s",
        first + j - 1L, conditionMessage(e)), call. = FALSE))
    out[, j] <- value - base
  }
  out
}

.np_udist_gnn_draws <- function(prepared, base, B, law, blocklen) {
  n <- prepared$n
  .npRmpi_bootstrap_fanout_enabled(n = n, B = B, chunk.size = 32L,
    what = "exact GNN CDF process")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B, 32L, with.seeds = FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks, function(task)
    invisible(.np_bias_process_innovations(n, task$start,
      task$start + task$bsz - 1L, B, law, blocklen)))
  worker <- function(task, prepared, base, B, law, blocklen) {
    .npRmpi_bootstrap_use_task_rng(task)
    t(.np_udist_gnn_chunk(prepared, base, task$start,
      task$start + task$bsz - 1L, B, law, blocklen))
  }
  t(.npRmpi_bootstrap_run_fanout(tasks, worker, ncol.out = length(base),
    what = "exact GNN CDF process", profile.where = "mpi.applyLB:udist-gnn-process",
    progress.label = "Calibrating exact GNN CDF process", prepared = prepared,
    base = base, B = B, law = law, blocklen = blocklen))
}

.np_udist_gnn_stream <- function(state, grid, law, alpha, B, corrected, blocklen, start) {
  m <- nrow(grid)
  budget <- .np_udist_gnn_budget(state, m, B, 1L + .npRmpi_bootstrap_worker_count())
  center <- original <- se <- numeric(m)
  for (block in .np_bias_target_blocks(m, 1L, anchors = 25L)) {
    ii <- block$grid
    main <- .np_udist_gnn_main(state, grid[ii, , drop = FALSE], corrected)
    center[ii] <- main$center; original[ii] <- main$original.center
  }
  limits <- .np_bias_limits_prepare(rep("cdf", m), alpha, FALSE, B)
  end <- NULL; done <- 0L
  for (block in .np_udist_gnn_blocks(state, grid, budget$marginal.tile)) {
    prepared <- .np_udist_gnn_block(state, block, corrected)
    # Coupled evaluation of identical sums preserves exactly zero displacement
    # for an unchanged sample, while the displayed observed center stays intact.
    base <- .np_udist_gnn_center(prepared, rep.int(1, state$n))$center
    stream <- .np_bias_stream(start,
      .np_udist_gnn_draws(prepared, base, B, law, blocklen))
    if (!is.null(end) && !identical(end, stream$end))
      stop("exact GNN CDF tile RNG mismatch", call. = FALSE)
    end <- stream$end
    se[block$rows] <- sqrt(apply(stream$value, 1L, stats::var))
    limits <- .np_bias_limits_update(limits, abs(stream$value), block$rows)
    done <- done + length(block$rows)
    .np_plot_progress_notify(stage = "Calibrating exact GNN CDF process",
      done = done, total = m, unit = "targets")
  }
  list(center = center, original = original, se = se, limits = limits,
    end = end, budget = budget)
}
