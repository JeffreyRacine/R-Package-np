# Exact resampling of the GNN statistic; no bootstrap studentizing pilot.
.np_udist_gnn_chunk <- function(prepared, base, first, last, B, law, blocklen) {
  prepared$zero.limit <- identical(law, "inid")
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
  out <- matrix(0, length(base), B)
  for (first in seq.int(1L, B, by = 32L)) {
    last <- min(B, first + 31L)
    out[, first:last] <- .np_udist_gnn_chunk(prepared, base, first, last, B, law, blocklen)
  }
  out
}

.np_udist_gnn_stream <- function(state, grid, law, alpha, B, corrected, blocklen, start) {
  m <- nrow(grid)
  budget <- .np_udist_gnn_budget(state, m, B, 1L)
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
