# Reuse the qualified covariance calibrator; no incumbent-SE rescaling.
.np_udens_bc_gaussian <- function(state, grid, alpha, B, tile, corrected, start) {
  moments <- .np_udens_bc_covariance(state, grid, tile, corrected)
  m <- nrow(grid)
  output <- .np_plot_gaussian_calibrate_alphas(moments$covariance,
    alpha = alpha, draws = B, max.bytes = state$max.bytes)
  end <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  point <- matrix(rep(qnorm(alpha/2, lower.tail = FALSE), each = m), m)
  bonf <- matrix(rep(qnorm(alpha/(2*m), lower.tail = FALSE), each = m), m)
  joint <- matrix(rep(vapply(output, function(z) z$critical, numeric(1L)), each = m), m)
  list(center = moments$center, original.center = moments$original.center,
       se = moments$se, covariance = moments$covariance,
       critical = list(pointwise = point, simultaneous = joint, bonferroni = bonf),
       raw.max = vapply(output, function(z) z$critical.raw, numeric(1L)),
       raw.bonf = bonf, degeneracy = integer(m), cancellation = integer(m),
       gaussian = output, end = end)
}

.np_udens_bc_factor <- function(state, grid, alpha, B, tile, corrected, start) {
  m <- nrow(grid)
  .np_bias_normal_replay_check("gaussian", NULL, m > tile, seed = start)
  limits <- .np_bias_limits_prepare(rep("density", m), alpha, TRUE, B)
  center <- original <- se <- numeric(m); end <- NULL
  for (block in .np_bias_target_blocks(m, 1L, anchors = tile)) {
    ii <- block$grid
    value <- .np_udens_bc_contribution(state, grid[ii, , drop = FALSE], corrected)
    if (any(!is.finite(value$se)) || any(value$se <= 0))
      stop("density inference requires positive finite observed scales", call. = FALSE)
    stream <- .np_bias_stream(start, .np_bias_process_draws(value$factor, B, "gaussian", NULL))
    if (!is.null(end) && !identical(end, stream$end))
      stop("Gaussian factor tiles consumed inconsistent streams", call. = FALSE)
    end <- stream$end
    limits <- .np_bias_limits_update(limits, abs(stream$value/value$se), ii)
    center[ii] <- value$center; original[ii] <- value$original.center; se[ii] <- value$se
  }
  finished <- .np_bias_limits_finish(limits)
  critical <- lapply(finished$limits, function(z) matrix(z[, 2L, ], m, length(alpha)))
  list(center = center, original.center = original, se = se, covariance = NULL,
       critical = critical, raw.max = drop(finished$raw_max_critical),
       raw.bonf = critical$bonferroni, degeneracy = integer(m),
       cancellation = integer(m), end = end)
}
