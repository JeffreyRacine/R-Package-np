# Fit-local reuse of exact existing block innovations. The cap is independent
# of target count; larger jobs retain the incumbent bounded replay owner.
.np_udist_cache_plan <- function(n, B, m, tile, blocklen, budget, max.bytes,
                                seed.length) {
  payload <- 8 * as.double(n) * B
  # Tasks are disjoint across ranks: this bounds all retained/serialized
  # innovation copies across the pool, and also any single process.
  extra <- 64 * as.double(n) * B +
    128 * ceiling(B / 32) * (seed.length + 16)
  if (m <= tile || blocklen <= 1L || !is.finite(payload) ||
      payload > 16 * 1024^2 || !is.finite(extra) ||
      budget$bytes + extra > max.bytes)
    return(NULL)
  list(payload.bytes = payload, extra.bytes = extra,
       bytes = budget$bytes + extra,
       cap.bytes = 16 * 1024^2, draw.chunk = 32L)
}

.np_udist_cache_prepare <- function(n, B, law, blocklen) {
  starts <- seq.int(1L, B, by = 32L)
  tasks <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    first <- starts[i]; size <- min(32L, B - first + 1L)
    entry <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
    weights <- .np_bias_process_innovations(n, first, first + size - 1L,
                                           B, law, blocklen)
    tasks[[i]] <- list(start = first, bsz = size, seed = NA_integer_,
      rng_state = entry, weights = weights,
      rng_end = get(".Random.seed", .GlobalEnv, inherits = FALSE))
  }
  attr(tasks, "rng_final_state") <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  tasks
}

.np_udist_cache_chunk <- function(F, task) {
  value <- F %*% task$weights
  assign(".Random.seed", task$rng_end, .GlobalEnv)
  value
}

.np_udist_cache_draws <- function(F, tasks, B, law, blocklen) {
  .npRmpi_bootstrap_fanout_enabled(n = ncol(F), B = B, chunk.size = 32L,
                                 what = "corrected process")
  worker <- function(task, factors) {
    .npRmpi_bootstrap_use_task_rng(task)
    t(.np_udist_cache_chunk(factors, task))
  }
  t(.npRmpi_bootstrap_run_fanout(tasks, worker, ncol.out = nrow(F),
    what = "corrected process", profile.where = "mpi.applyLB:bias-process",
    progress.label = "Calibrating corrected process",
    required.bindings = list(.np_udist_cache_chunk = .np_udist_cache_chunk), factors = F))
}
