# Q6 transport owns only primitive draw replay and scalar replicate maxima.
# Logical blocks are independent of rank count. All donors precede all noise.
.np_udens_q6_plan <- function(start, n, B, chunk = 32L) {
  .np_udens_q6_budget(n, 1L, B, chunk)
  if (B < 1L || as.double(n) * chunk > .Machine$integer.max)
    stop("invalid Q6 replicate count or block", call. = FALSE)
  n <- as.integer(n); B <- as.integer(B); chunk <- as.integer(chunk)
  .np_bias_replay_check("inid", NULL, TRUE, seed = start)
  .np_bias_replay_check("gaussian", NULL, TRUE, seed = start)
  starts <- seq.int(1L, B, by = chunk)
  stream <- .np_bias_stream(start, {
    tasks <- vector("list", length(starts))
    for (i in seq_along(starts)) {
      size <- min(chunk, B - starts[i] + 1L)
      state <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
      invisible(sample.int(n, n * size, replace = TRUE))
      tasks[[i]] <- list(start = starts[i], bsz = size, seed = NA_integer_,
        rng_state = state, donor = state,
        donor.end = get(".Random.seed", .GlobalEnv, inherits = FALSE))
    }
    for (i in seq_along(tasks)) {
      tasks[[i]]$noise <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
      invisible(rnorm(n * tasks[[i]]$bsz))
      tasks[[i]]$noise.end <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
      .np_plot_progress_notify(stage = "Preparing Q6 conditional draws",
        done = i, total = length(tasks), unit = "blocks")
    }
    tasks
  })
  attr(stream$value, "rng_final_state") <- stream$end
  stream$value
}

.np_udens_q6_samples <- function(task, x, h) {
  n <- length(x)
  donor <- .np_bias_stream(task$donor,
    sample.int(n, n * task$bsz, replace = TRUE))
  noise <- .np_bias_stream(task$noise, rnorm(n * task$bsz))
  if (!identical(donor$end, task$donor.end) ||
      !identical(noise$end, task$noise.end))
    stop("Q6 draw replay consumed inconsistent streams", call. = FALSE)
  probability <- pnorm(noise$value)
  if (any(!is.finite(probability)) || any(probability <= 0 | probability >= 1))
    stop("Q6 normal transform reached an unrepresentable probability", call. = FALSE)
  value <- x[donor$value] + h * sqrt(2/3) * qt(probability, df = 6)
  if (any(!is.finite(value))) stop("nonfinite Q6 conditional sample", call. = FALSE)
  matrix(value, n, task$bsz)
}

.np_udens_q6_chunk <- function(task, reference) {
  Y <- .np_udens_q6_samples(task, reference$x, reference$h)
  z <- .np_udens_q6_maximum(Y, reference)
  cbind(maximum = z$maximum, refinement = rep(z$refinement, task$bsz))
}

.np_udens_q6_dispatch <- function(tasks, reference) {
  .npRmpi_bootstrap_fanout_enabled(n = length(reference$x),
    B = sum(vapply(tasks, function(z) z$bsz, 0L)), chunk.size = 32L,
    what = "Q6 density simultaneous process")
  worker <- function(task, reference) .np_udens_q6_chunk(task, reference)
  .npRmpi_bootstrap_run_fanout(tasks, worker, ncol.out = 2L,
    what = "Q6 density simultaneous process",
    profile.where = "mpi.applyLB:q6-density",
    progress.label = "Constructing Q6 simultaneous bands",
    required.bindings = list(
      .np_udens_q6_chunk = .np_udens_q6_chunk,
      .np_udens_q6_samples = .np_udens_q6_samples,
      .np_udens_q6_maximum = .np_udens_q6_maximum,
      .np_udens_q6_curve = .np_udens_q6_curve,
      .np_udens_q6_convolution = .np_udens_q6_convolution,
      .np_udens_q6_rule = .np_udens_q6_rule,
      .np_udens_q6_kernel = .np_udens_q6_kernel,
      .np_udens_q6_check = .np_udens_q6_check),
    reference = reference)
}
