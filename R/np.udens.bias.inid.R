# Literal multinomial bootstrap-t. Extended pivots are private to density.
.np_udens_bc_inid_pivot <- function(F, raw, counts, accurate = FALSE) {
  if (!is.matrix(F) || !is.numeric(F) || !is.matrix(raw) ||
      !identical(dim(F), dim(raw)) || !is.matrix(counts) ||
      ncol(F) < 2L || nrow(counts) != ncol(F) ||
      any(!is.finite(F)) || any(!is.finite(raw)) || any(!is.finite(counts)) ||
      any(counts < 0) || any(counts != floor(counts)) ||
      any(colSums(counts) != ncol(F)))
    stop("invalid density inid contributions or counts", call. = FALSE)
  n <- ncol(F)
  numerator <- F %*% counts
  moment <- (F^2) %*% counts
  variance <- moment - numerator^2/n
  if (any(!is.finite(variance)) || any(!is.finite(numerator)))
    stop("nonfinite density replicate variance or numerator", call. = FALSE)
  certified <- matrix(FALSE, nrow(F), ncol(counts))
  cancellation <- integer(nrow(F))
  # This is a cancellation diagnostic, not a variance threshold/floor. Every
  # suspect cell is recomputed from centered raw observations and its counts.
  # Mixed rows can concentrate on few profiles. Recompute ill-conditioned
  # subtractions to an accuracy target; never change or floor their variance.
  guard <- 128*.Machine$double.eps / if (accurate) 1e-12 else 1
  suspect <- which(variance <= guard*(moment + numerator^2/n), arr.ind = TRUE)
  if (nrow(suspect)) for (i in seq_len(nrow(suspect))) {
    r <- suspect[i, 1L]; b <- suspect[i, 2L]
    active <- which(counts[, b] > 0)
    values <- raw[r, active]; w <- counts[active, b]
    if (all(values == values[1L])) {
      variance[r, b] <- 0
      if (numerator[r, b] == 0)
        stop("undefined zero-over-zero density bootstrap pivot", call. = FALSE)
      certified[r, b] <- TRUE
    } else {
      dev <- (values - values[1L])/n
      mean.dev <- sum(w*dev)/n
      direct <- sum(w*(dev-mean.dev)^2)
      if (!is.finite(direct) || direct <= 0)
        stop("uncertified nonpositive density replicate variance", call. = FALSE)
      variance[r, b] <- direct
      cancellation[r] <- cancellation[r] + 1L
    }
  }
  if (any(variance[!certified] <= 0))
    stop("negative density replicate variance", call. = FALSE)
  absolute <- matrix(Inf, nrow(F), ncol(counts))
  absolute[!certified] <- abs(numerator[!certified]/sqrt(variance[!certified]))
  if (any(!is.finite(absolute[!certified])))
    stop("uncertified nonfinite density bootstrap pivot", call. = FALSE)
  list(absolute = absolute, degeneracy = rowSums(certified),
       cancellation = cancellation)
}

.np_udens_bc_inid_chunk <- function(F, raw, first, last, B, accurate = FALSE) {
  counts <- .np_bias_process_innovations(ncol(F), first, last, B, "inid", 1L) + 1
  .np_udens_bc_inid_pivot(F, raw, counts, accurate)
}

.np_udens_bc_inid_dispatch <- function(F, raw, B, accurate = FALSE) {
  absolute <- matrix(0, nrow(F), B)
  degenerate <- cancellation <- integer(nrow(F))
  for (first in seq.int(1L, B, by = 32L)) {
    last <- min(B, first + 31L)
    z <- .np_udens_bc_inid_chunk(F, raw, first, last, B, accurate)
    absolute[, first:last] <- z$absolute
    degenerate <- degenerate + z$degeneracy
    cancellation <- cancellation + z$cancellation
    .np_plot_progress_notify(stage = "Calibrating density inid process",
                             done = last, total = B, unit = "draws")
  }
  list(absolute = absolute, degeneracy = degenerate, cancellation = cancellation)
}

.np_udens_bc_limits_update <- function(state, absolute, rows, certified) {
  if (nrow(absolute) != length(rows) || ncol(absolute) != ncol(state$maximum) ||
      anyNA(absolute) || any(absolute < 0) || length(certified) != nrow(absolute) ||
      any(rowSums(is.infinite(absolute)) != certified))
    stop("invalid or uncertified density extended-real draws", call. = FALSE)
  for (g in seq_along(state$index)) {
    local <- which(rows %in% state$index[[g]])
    if (!length(local)) next
    maximum <- state$maximum[g, ]
    for (j in local) maximum <- pmax(maximum, absolute[j, ])
    state$maximum[g, ] <- maximum
    for (j in local) {
      q <- unname(quantile(absolute[j, ],
        c(1-state$alpha, 1-state$alpha/length(state$index[[g]])), type = 1L))
      state$point[rows[j], ] <- q[seq_along(state$alpha)]
      state$bonf[rows[j], ] <- q[length(state$alpha) + seq_along(state$alpha)]
    }
  }
  state
}

.np_udens_bc_inid <- function(state, grid, alpha, B, tile, corrected, start) {
  m <- nrow(grid)
  limits <- .np_bias_limits_prepare(rep("density", m), alpha, FALSE, B)
  center <- original <- se <- numeric(m)
  degenerate <- cancellation <- integer(m); visited <- logical(m); end <- NULL
  for (block in .np_bias_target_blocks(m, 1L, anchors = tile)) {
    ii <- block$grid
    if (any(visited[ii])) stop("density target updated twice", call. = FALSE)
    value <- .np_udens_bc_contribution(state, grid[ii, , drop = FALSE], corrected)
    if (any(!is.finite(value$se)) || any(value$se <= 0))
      stop("density inference requires positive finite observed scales", call. = FALSE)
    stream <- .np_bias_stream(start, .np_udens_bc_inid_dispatch(value$factor, value$corrected, B,
      accurate = state$r + state$s > 0L))
    if (!is.null(end) && !identical(end, stream$end))
      stop("density inid tiles consumed inconsistent streams", call. = FALSE)
    end <- stream$end; z <- stream$value
    limits <- .np_udens_bc_limits_update(limits, z$absolute, ii, z$degeneracy)
    center[ii] <- value$center; original[ii] <- value$original.center; se[ii] <- value$se
    degenerate[ii] <- z$degeneracy; cancellation[ii] <- z$cancellation
    visited[ii] <- TRUE
  }
  if (!all(visited)) stop("density target was not updated", call. = FALSE)
  raw.bonf <- limits$bonf
  finished <- .np_bias_limits_finish(limits)
  critical <- lapply(finished$limits, function(z) matrix(z[, 2L, ], m, length(alpha)))
  list(center = center, original.center = original, se = se, covariance = NULL,
       critical = critical, raw.max = drop(finished$raw_max_critical),
       raw.bonf = raw.bonf, degeneracy = degenerate, cancellation = cancellation, end = end)
}
