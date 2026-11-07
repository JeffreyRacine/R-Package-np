# Fixed-pilot complete corrected curves; direct ranks, without SD transport.
.np_udens_smoothed_streams <- function(start) {
  if (start[1L] %% 100L == 7L)
    return(list(count = start, jitter = parallel::nextRNGStream(start)))
  seed <- 104729
  for (v in as.double(start)) seed <- (65599 * seed + v %% 2147483646) %% 2147483646
  jitter <- .np_bias_stream(start, {
    set.seed(as.integer(seed + 1))
    get(".Random.seed", .GlobalEnv, inherits = FALSE)
  })$value
  list(count = start, jitter = jitter)
}

# Whole observations share one donor. Numeric order is occurrence, coordinate,
# replicate, so chunk boundaries cannot select different innovations.
.np_udens_smoothed_samples <- function(x, h, first, size, B, law, blocklen, streams, layout = NULL) {
  n <- nrow(x); q <- ncol(x)
  counts <- .np_bias_stream(streams$count,
    .np_bias_process_innovations(n, first, first + size - 1L, B, law, blocklen) + 1)
  if (!is.null(layout)) {
    jitter <- .np_bias_stream(streams$jitter, {
      value <- array(0, c(n, q, size))
      for (b in seq_len(size)) {
        donor <- rep(seq_len(n), counts$value[, b])
        for (j in seq_len(q)) {
          base <- x[donor, j]
          if (layout$numeric[j]) value[, j, b] <- base + h[j] * rnorm(n) else {
            category <- layout$category[[match(j, vapply(layout$category, function(z) z$index, integer(1L)))]]
            lambda <- category$lambda
            if (category$ordered) {
              value[, j, b] <- base + rgeom(n, 1-lambda) - rgeom(n, 1-lambda)
            } else {
              u <- runif(n); selected <- match(base, category$support)
              change <- u >= 1-lambda
              if (lambda > 0 && any(change)) {
                other <- floor((u[change]-(1-lambda))/lambda*(category$k-1)) + 1L
                selected[change] <- other + (other >= selected[change])
              }
              value[, j, b] <- category$support[selected]
            }
          }
        }
      }
      value
    })
    return(list(y = jitter$value, end = list(count = counts$end, jitter = jitter$end)))
  }
  jitter <- .np_bias_stream(streams$jitter, array(rnorm(n * q * size), c(n, q, size)))
  y <- jitter$value
  for (b in seq_len(size)) {
    donor <- rep(seq_len(n), counts$value[, b])
    for (j in seq_len(q)) y[, j, b] <- x[donor, j] + h[j] * y[, j, b]
  }
  list(y = y, end = list(count = counts$end, jitter = jitter$end))
}

.np_udens_smoothed_plan <- function(x, h, B, law, blocklen, streams, chunk, cache, layout = NULL) {
  .np_bias_replay_check(law, blocklen, !cache, seed = streams$count)
  normal.jitter <- if (is.null(layout)) ncol(x) > 0L else
    any(layout$numeric) || any(vapply(layout$category,
      function(z) isTRUE(z$ordered) && isTRUE(z$lambda > 0), logical(1L)))
  .np_bias_replay_check("gaussian", NULL, !cache && normal.jitter,
    seed = streams$jitter)
  uniform.category <- !is.null(layout) && any(vapply(layout$category,
    function(z) !isTRUE(z$ordered) || isTRUE(z$lambda > 0), logical(1L)))
  .np_bias_replay_check("inid", NULL, !cache && uniform.category,
    seed = streams$jitter)
  starts <- seq.int(1L, B, by = chunk)
  tasks <- vector("list", length(starts)); current <- streams
  for (i in seq_along(starts)) {
    first <- starts[i]; size <- min(chunk, B - first + 1L)
    z <- .np_udens_smoothed_samples(x, h, first, size, B, law, blocklen, current, layout)
    tasks[[i]] <- list(start = first, bsz = size, seed = NA_integer_,
      rng_state = current$count, streams = current, end = z$end,
      y = if (cache) z$y else NULL)
    current <- z$end
  }
  attr(tasks, "rng_final_state") <- current$count
  attr(tasks, "streams_end") <- current
  tasks
}

# At most tile * n * chunk contributions live at once. The complete numeric
# correction is K - K * (sum(u^2) - q)/2, not a product of corrected kernels.
.np_udens_smoothed_chunk <- function(task, x, h, grid, B, law, blocklen, layout = NULL) {
  y <- task$y
  if (is.null(y)) {
    draw <- .np_udens_smoothed_samples(x, h, task$start, task$bsz, B, law, blocklen, task$streams, layout)
    if (!identical(draw$end, task$end))
      stop("smoothed density replay consumed inconsistent streams", call. = FALSE)
    y <- draw$y
  }
  n <- nrow(x); numeric <- if (is.null(layout)) rep(TRUE, ncol(x)) else layout$numeric
  q <- sum(numeric); m <- nrow(grid)
  K <- matrix(1, m, n * task$bsz); radius <- matrix(0, m, n * task$bsz)
  for (j in which(numeric)) {
    u <- outer(grid[[j]], as.vector(y[, j, ]), "-") / h[j]
    K <- K * dnorm(u) / h[j]
    radius <- radius + u * u
  }
  L <- K - K * (radius - q) / 2
  if (!is.null(layout)) for (category in layout$category) {
    j <- category$index
    category$training <- as.vector(y[, j, ])
    value <- .np_udens_bc_categorical_rows(category, grid[[j]])
    L <- L * value$corrected
  }
  out <- matrix(0, m, task$bsz)
  for (j in seq_len(m)) out[j, ] <- colMeans(matrix(L[j, ], n, task$bsz))
  if (any(!is.finite(out))) stop("nonfinite smoothed density contributions", call. = FALSE)
  out
}

.np_udens_smoothed_draws <- function(tasks, x, h, grid, B, law, blocklen, layout = NULL) {
  out <- matrix(0, nrow(grid), B)
  for (task in tasks) {
    jj <- seq.int(task$start, task$start + task$bsz - 1L)
    out[, jj] <- .np_udens_smoothed_chunk(task, x, h, grid, B, law, blocklen, layout)
  }
  out
}
