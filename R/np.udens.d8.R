# Private plot-inference provider. Public default activation is separately gated.
# The original density estimator, bandwidth selector and estimator SEs are not used
# as a surrogate for the complete covariance of this corrected plot curve.
.np_udens_d8_prepare <- function(bws, x, max.bytes = 256 * 1024^2) {
  state <- .np_udens_bc_prepare(bws, x, max.bytes)
  if (state$r != 0L || state$s != 0L || !state$q %in% 1:3)
    stop("D8 plot inference requires one to three continuous coordinates", call. = FALSE)
  state$construction <- "density-d8-gaussian-v1"
  state
}

# Total-degree radial Gaussian correction; not a product of Gaussian8 kernels.
.np_udens_d8_rows <- function(state, grid) {
  t <- matrix(0, nrow(grid), state$n)
  for (j in seq_len(state$q))
    t <- t + (outer(grid[[j]], state$x[[j]], "-") / state$h[j])^2 / 2
  a <- 1 + state$q / 2
  log.normal <- -state$q * log(2*pi) / 2 - sum(log(state$h))
  gaussian <- exp(log.normal - t)
  original <- rowMeans(gaussian)
  # Keep ordinary evaluation close to the frozen arithmetic oracle. Beyond 100,
  # factor out t^3 and evaluate in logs, avoiding 0*Inf and premature underflow.
  ordinary <- is.finite(t) & t <= 100
  value <- matrix(0, nrow(grid), state$n)
  u <- t[ordinary]
  value[ordinary] <- gaussian[ordinary] *
    (a*(a+1)*(a+2)/6 - (a+1)*(a+2)*u/2 + (a+2)*u^2/2 - u^3/6)
  far <- is.finite(t) & !ordinary
  u <- t[far]
  reciprocal <- 1/u
  polynomial <- -1/6 + (a+2)*reciprocal/2 -
    (a+1)*(a+2)*reciprocal^2/2 + a*(a+1)*(a+2)*reciprocal^3/6
  value[far] <- sign(polynomial) *
    exp(log.normal - u + 3*log(u) + log(abs(polynomial)))
  if (any(!is.finite(value)) || any(!is.finite(original)))
    stop("nonfinite complete D8 density contributions", call. = FALSE)
  center <- rowMeans(value)
  factor <- (value - center) / state$n
  scale <- sqrt(rowSums(factor^2))
  if (any(!is.finite(factor)) || any(!is.finite(scale)))
    stop("nonfinite complete D8 density covariance", call. = FALSE)
  list(center = center, original.center = original, factor = factor, se = scale)
}

# The shared forecast covers factor/draw work. Charge the density-local
# multi-alpha endpoint/width payload and construction temporaries separately.
.np_udens_d8_budget <- function(state, m, B, tile, levels) {
  budget <- .np_udens_b2_budget(state, m, B, tile, FALSE)
  extra <- 8 * as.double(m) * (16 + 12 * as.double(levels))
  budget$per.process <- budget$per.process + extra
  budget$bytes <- budget$per.process * budget$processes
  budget$d8.payload.per.process <- extra
  if (!is.finite(budget$bytes) || budget$bytes > state$max.bytes)
    stop("D8 density payload exceeds the requested byte budget", call. = FALSE)
  budget
}

.np_udens_d8_inference <- function(state, grid, alpha = .05, B = 9999L,
                                    tile = 25L) {
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
  if (!identical(state$construction, "density-d8-gaussian-v1") ||
      state$r != 0L || state$s != 0L || !state$q %in% 1:3)
    stop("invalid D8 plot-inference state", call. = FALSE)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= .5) || anyDuplicated(alpha))
    stop("density alpha values must be unique and in (0, 0.5)", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "density draws", .Machine$integer.max)
  .np_udens_bc_grid(state, grid)
  m <- nrow(grid)
  budget <- .np_udens_d8_budget(state, m, B, tile, length(alpha))
  tile <- budget$tile
  if (m > 1L) for (level in alpha) .np_plot_gaussian_draws(B, level, m)
  if (!existed) sample.int(1L, 0L)
  start <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  end <- NULL; maximum <- numeric(B)
  center <- original <- scale <- numeric(m)
  for (first in seq.int(1L, m, by = tile)) {
    ii <- seq.int(first, min(m, first + tile - 1L))
    z <- .np_udens_d8_rows(state, grid[ii, , drop = FALSE])
    stream <- .np_bias_stream(start,
      .np_bias_process_draws(z$factor, B, "gaussian", NULL))
    if (!identical(dim(stream$value), c(length(ii), B)) ||
        any(!is.finite(stream$value)))
      stop("invalid D8 Gaussian process draws", call. = FALSE)
    if (!is.null(end) && !identical(end, stream$end))
      stop("D8 target tiles consumed inconsistent Gaussian streams", call. = FALSE)
    end <- stream$end
    for (j in seq_along(ii)) maximum <- pmax(maximum, abs(stream$value[j, ]))
    center[ii] <- z$center; original[ii] <- z$original.center; scale[ii] <- z$se
    .np_plot_progress_notify(stage = "Constructing corrected density covariance",
      done = max(ii), total = m, unit = "targets")
  }
  point <- outer(scale, qnorm(1-alpha/2))
  raw <- unname(quantile(maximum, 1-alpha, type = 7L))
  critical <- pmax(raw, qnorm(1-alpha/2)*max(scale))
  if (m == 1L) critical <- qnorm(1-alpha/2)*scale
  simultaneous <- matrix(critical, m, length(alpha), byrow = TRUE)
  bonferroni <- pmax(outer(scale, qnorm(1-alpha/(2*m))), simultaneous)
  widths <- list(pointwise = point, simultaneous = simultaneous, bonferroni = bonferroni)
  intervals <- lapply(widths, function(width) {
    out <- array(0, c(m, 2L, length(alpha)),
      dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    out[, 1L, ] <- center - width
    out[, 2L, ] <- center + width
    if (any(!is.finite(out))) stop("nonfinite D8 density endpoints", call. = FALSE)
    out
  })
  assign(".Random.seed", end, .GlobalEnv)
  result <- list(center = center, original.center = original, se = scale,
    empirical.se = scale, covariance = NULL, bands = intervals,
    intervals = intervals, absolute.intervals = intervals, half.width = widths,
    maximum = maximum, raw.maximum.quantile = raw, maximum.half.width = critical,
    settings = list(construction = "density-d8-gaussian-v1", law = "gaussian",
      corrected = TRUE, h = state$h, n = state$n, q = state$q, r = 0L, s = 0L,
      alpha = alpha, B = B, grid.size = m, budget = budget,
      observed.scale = "complete centered D8 observation contributions / n",
      covariance = "F F-transpose, denominator n-squared",
      representation = "bounded observation factors", quantile.type = 7L,
      maximum.statistic = "unstudentized absolute Gaussian maximum",
      family = "one supplied finite grid", available.bands = names(intervals),
      calibration.status = "private provider; automatic activation held by pointwise coverage admission"),
    rng = list(kind = kind, entry = previous, start = start, end = end))
  success <- TRUE
  result
}
