# Private recovery provider. F3 center, F4 marginal scales, complete empirical
# corrected correlation. This intentionally differs from F4 pilot correlation.
.np_udens_reference_scale <- function(state, grid) {
  vapply(seq_len(nrow(grid)), function(j) {
    target <- grid[j, , drop = FALSE]
    mean1 <- mean(.np_udens_b2_pilot_moment(state, target))
    second1 <- mean(.np_udens_b2_pilot_moment(state, target, target))
    mean2 <- mean(.np_udens_b4_pilot_moment(state, target, tau = 2))
    second2 <- mean(.np_udens_b4_pilot_moment(state, target, target, tau = 2))
    variance <- (second1 - mean1^2) / state$n
    wide <- (second2 - mean2^2) / state$n
    if (!is.finite(variance) || !is.finite(wide) || variance <= 0 || wide <= 0)
      stop("fixed-reference density requires positive finite pilot variances", call. = FALSE)
    .np_plot_progress_notify(stage = "Correcting density scales", done = j,
                             total = nrow(grid), unit = "targets")
    variance / sqrt(wide)
  }, numeric(1L))
}

.np_udens_reference_inference <- function(state, grid, law = c("gaussian", "inid"),
    alpha = .05, B = 9999L, tile = 25L,
    bands = c("pointwise", "simultaneous", "bonferroni")) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (existed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind(); success <- FALSE
  on.exit({
    if (!success) {
      do.call(RNGkind, as.list(kind))
      if (existed) assign(".Random.seed", previous, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  law <- match.arg(law)
  bands <- match.arg(bands, c("pointwise", "simultaneous", "bonferroni"), several.ok = TRUE)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= .5) || anyDuplicated(alpha))
    stop("density alpha values must be unique and in (0, .5)", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "density draws", .Machine$integer.max)
  if (B < 2L) stop("density calibration requires at least two draws", call. = FALSE)
  .np_udens_bc_grid(state, grid)
  m <- nrow(grid)
  bonferroni.available <- law == "gaussian" || all(alpha * B / m >= 1)
  if (!bonferroni.available && "bonferroni" %in% bands)
    stop("fixed-reference density bootstrap Bonferroni is unqualified at this grid/draw resolution; no calibration law or draw count was changed", call. = FALSE)
  covariance.mode <- law == "gaussian" && m <= 512L
  budget <- .np_udens_b2_budget(state, m, B, tile, covariance.mode)
  scale <- .np_udens_reference_scale(state, grid)
  if (any(!is.finite(scale)) || any(scale <= 0))
    stop("fixed-reference density requires positive finite observed scales", call. = FALSE)
  if (!existed) runif(1L)
  start <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  covariance <- NULL; end <- NULL
  if (covariance.mode) {
    model <- .np_udens_bc_covariance(state, grid, tile)
    center <- model$center; original <- model$original.center; empirical <- model$se
    covariance <- model$covariance * tcrossprod(scale / empirical)
    calibrated <- lapply(alpha, function(a) .np_bias_stream(start,
      .np_plot_gaussian_calibrate(covariance, alpha = a, draws = B,
                                 max.bytes = state$max.bytes)))
    end <- calibrated[[1L]]$end
    if (!all(vapply(calibrated, function(z) identical(z$end, end), logical(1L))))
      stop("density Gaussian alpha replays consumed inconsistent streams", call. = FALSE)
    critical <- list(
      pointwise = matrix(rep(qnorm(alpha/2, lower.tail = FALSE), each = m), m),
      simultaneous = matrix(rep(vapply(calibrated, function(z) z$value$critical,
                                      numeric(1L)), each = m), m),
      bonferroni = matrix(rep(qnorm(alpha/(2*m), lower.tail = FALSE), each = m), m))
    raw.critical <- critical
    raw.critical$simultaneous <- matrix(rep(vapply(calibrated,
      function(z) z$value$critical.raw, numeric(1L)), each = m), m)
  } else {
    limits <- .np_bias_limits_prepare(rep("density", m), alpha, law == "gaussian", B)
    center <- original <- empirical <- numeric(m)
    for (block in .np_bias_target_blocks(m, 1L, anchors = tile)) {
      ii <- block$grid
      value <- .np_udens_bc_contribution(state, grid[ii, , drop = FALSE])
      if (any(!is.finite(value$se)) || any(value$se <= 0))
        stop("fixed-reference density requires positive finite empirical scales", call. = FALSE)
      stream <- .np_bias_stream(start,
        .np_bias_process_draws(value$factor, B, law, NULL))
      if (!is.null(end) && !identical(end, stream$end))
        stop("density process tiles consumed inconsistent streams", call. = FALSE)
      end <- stream$end
      limits <- .np_bias_limits_update(limits, abs(stream$value/value$se), ii)
      center[ii] <- value$center; original[ii] <- value$original.center
      empirical[ii] <- value$se
    }
    finished <- .np_bias_limits_finish(limits)
    critical <- lapply(finished$limits, function(z)
      matrix(z[, 2L, ], m, length(alpha)))
    raw.critical <- list(
      pointwise = if (law == "gaussian") critical$pointwise else limits$point,
      simultaneous = matrix(rep(finished$raw_max_critical[, 1L], each = m), m),
      bonferroni = if (law == "gaussian") critical$bonferroni else limits$bonf)
  }
  available <- if (bonferroni.available) names(critical) else
    setdiff(names(critical), "bonferroni")
  intervals <- lapply(critical[available], function(k) {
    out <- array(0, c(m, 2L, length(alpha)),
                 dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    out[, 1L, ] <- center - scale * k
    out[, 2L, ] <- center + scale * k
    if (any(!is.finite(out))) stop("nonfinite fixed-reference density endpoints", call. = FALSE)
    out
  })
  result <- list(center = center, original.center = original, se = scale,
    empirical.se = empirical, covariance = covariance, bands = intervals,
    intervals = intervals, absolute.intervals = intervals,
    critical = critical[available], raw.critical = raw.critical,
    settings = list(construction = "density-fixed-reference-v1", law = law,
      corrected = TRUE, h = state$h, n = state$n, q = state$q, r = state$r,
      s = state$s, alpha = alpha, B = B, grid.size = m, budget = budget,
      observed.scale = "F4 positive marginal correction",
      covariance = "F4 marginal congruence of complete empirical corrected correlation",
      bootstrap.statistic = "observed variance transport of literal centered multinomial process",
      representation = if (covariance.mode) "bounded covariance" else "bounded observation factors",
      available.bands = available, bonferroni.tail.expected = alpha * B / m,
      bonferroni.raw.rank = ceiling(B * (1-alpha/m)),
      calibration.status = "private recovery candidate; coverage unqualified"),
    rng = list(kind = kind, entry = previous, start = start, end = end))
  assign(".Random.seed", end, .GlobalEnv)
  success <- TRUE
  result
}
