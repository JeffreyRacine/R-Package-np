# Exact expanded-observation NN refits, with immutable counts and lambdas.
.np_udist_bc_nn_indices <- function(n, first, last, B, law, blocklen) {
  if (law == "inid" || blocklen == 1L) {
    N <- .np_bias_process_innovations(n, first, last, B, law, blocklen) + 1
    return(lapply(seq_len(ncol(N)), function(j) rep.int(seq_len(n), N[, j])))
  }
  draw <- .np_block_replicate_drawer(n, n, blocklen, law, TRUE)
  lapply(seq_len(last-first+1L), function(j) draw())
}

.np_udist_bc_nn_chunk <- function(state, grid, base, first, last, B, law, blocklen, corrected) {
  index <- .np_udist_bc_nn_indices(state$n, first, last, B, law, blocklen)
  out <- matrix(0, nrow(grid), length(index)); degenerate <- integer(nrow(grid))
  for (b in seq_along(index)) {
    sample <- state; sample$x <- state$x[index[[b]], , drop = FALSE]
    for (j in seq_along(sample$category))
      sample$category[[j]]$training <- state$category[[j]]$training[index[[b]]]
    z <- .np_udist_bc_gnn_rows(sample, grid, corrected)
    variance <- if (law == "inid") rowSums(z$factor^2) else
      .np_bias_score_variance(z$factor, blocklen, law)
    numerator <- z$center-base$center
    if (any(!is.finite(variance)) || any(variance < 0)) stop("invalid exact NN CDF variance", call. = FALSE)
    zero <- which(variance == 0)
    if (length(zero)) for (j in zero) {
      if (!all(z$corrected[j, ] == z$corrected[j, 1L]) || numerator[j] == 0)
        stop("undefined or uncertified exact NN CDF pivot", call. = FALSE)
      degenerate[j] <- degenerate[j] + 1L
    }
    out[, b] <- abs(numerator/sqrt(variance))
    if (anyNA(out[, b])) stop("undefined exact NN CDF pivot", call. = FALSE)
  }
  list(absolute = out, degeneracy = degenerate, cancellation = integer(nrow(grid)))
}

.np_udist_bc_nn_dispatch <- function(state, grid, base, B, law, blocklen, corrected) {
  .npRmpi_bootstrap_fanout_enabled(n = state$n, B = B, chunk.size = 32L, what = "exact NN CDF process")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B, 32L, with.seeds = FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks, function(task)
    invisible(.np_udist_bc_nn_indices(state$n, task$start, task$start+task$bsz-1L, B, law, blocklen)))
  worker <- function(task, state, grid, base, B, law, blocklen, corrected) {
    .npRmpi_bootstrap_use_task_rng(task)
    z <- .np_udist_bc_nn_chunk(state, grid, base, task$start, task$start+task$bsz-1L, B, law, blocklen, corrected)
    out <- t(z$absolute); attr(out, "cdf.pivot") <- z[c("degeneracy", "cancellation")]
    out
  }
  reducer <- function(out, parts, tasks) {
    d <- lapply(parts, attr, which = "cdf.pivot", exact = TRUE)
    if (any(vapply(d, is.null, logical(1L)))) stop("CDF NN worker omitted diagnostics", call. = FALSE)
    attr(out, "cdf.pivot") <- list(degeneracy = Reduce(`+`, lapply(d, `[[`, "degeneracy")),
      cancellation = Reduce(`+`, lapply(d, `[[`, "cancellation")))
    out
  }
  deps <- list(.np_udist_bc_nn_chunk = .np_udist_bc_nn_chunk,
    .np_udist_bc_nn_indices = .np_udist_bc_nn_indices, .np_udist_bc_gnn_rows = .np_udist_bc_gnn_rows,
    .np_udist_bc_gnn_geometry = .np_udist_bc_gnn_geometry, .np_udist_bc_grid = .np_udist_bc_grid,
    .np_udist_bc_ordered = .np_udist_bc_ordered)
  out <- .npRmpi_bootstrap_run_fanout(tasks, worker, ncol.out = nrow(grid),
    what = "exact NN CDF process", profile.where = "mpi.applyLB:cdf-nn",
    progress.label = "Calibrating exact NN CDF process", progress.unit = "draws",
    required.bindings = deps, metadata.reducer = reducer, state = state, grid = grid,
    base = base, B = B, law = law, blocklen = blocklen, corrected = corrected)
  d <- attr(out, "cdf.pivot", exact = TRUE)
  list(absolute = t(out), degeneracy = d$degeneracy, cancellation = d$cancellation)
}
