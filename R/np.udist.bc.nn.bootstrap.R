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
  out <- matrix(0, nrow(grid), B); degenerate <- cancellation <- integer(nrow(grid))
  for (first in seq.int(1L, B, by = 32L)) {
    last <- min(B, first+31L)
    z <- .np_udist_bc_nn_chunk(state, grid, base, first, last, B, law, blocklen, corrected)
    out[, first:last] <- z$absolute
    degenerate <- degenerate+z$degeneracy; cancellation <- cancellation+z$cancellation
    .np_plot_progress_notify(stage = "Calibrating exact NN CDF process", done = last, total = B, unit = "draws")
  }
  list(absolute = out, degeneracy = degenerate, cancellation = cancellation)
}
