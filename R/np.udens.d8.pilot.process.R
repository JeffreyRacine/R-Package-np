# Explicit seeded private plot entry. Public dispatch must call this directly,
# never from inside the incumbent .np_with_seed wrapper under Box-Muller.
.np_udens_d8_pilot_seeded <- function(state, grid, alpha = .05, B = 9999L,
                                      tile = NULL, panel = 32L, batch = 256L,
                                      pointwise.only = FALSE, random.seed = 42L) {
  npValidateNonNegativeInteger(random.seed, "random.seed")
  kind <- RNGkind()
  if (kind[2L] != "Box-Muller")
    return(.np_with_seed(random.seed, .np_udens_d8_pilot_inference(state, grid,
      alpha, B, tile, panel, batch, pointwise.only)))
  request <- .np_plot_normal_request(random.seed, kind, 0)
  .np_udens_d8_pilot_inference(state, grid, alpha, B, tile, panel, batch,
    pointwise.only, .normal.request = request)
}

# Actual innovations are generated once and transported as values. Seed-only
# replay cannot reproduce a cached Box-Muller normal and is not used here.
.np_udens_d8_pilot_dispatch <- function(factor, innovations) {
  factor %*% innovations
}

.np_udens_d8_pilot_inference <- function(state, grid, alpha = .05, B = 9999L,
                                         tile = NULL, panel = 32L, batch = 256L,
                                         pointwise.only = FALSE, .normal.request = NULL) {
  if (!identical(state$construction, "density-d8-positive-pilot-v2") ||
      state$r != 0L || state$s != 0L || !state$q %in% 1:3)
    stop("invalid D8 pilot-inference state", call. = FALSE)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 1e-4 | alpha >= .5) || anyDuplicated(alpha))
    stop("D8 pilot alpha values must be unique and in (0.0001, 0.5)", call. = FALSE)
  if (!is.logical(pointwise.only) || length(pointwise.only) != 1L || is.na(pointwise.only))
    stop("invalid D8 pointwise request", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "density draws", .Machine$integer.max)
  .np_udens_bc_grid(state, grid)
  m <- nrow(grid)
  budget <- .np_udens_d8_pilot_budget(state, m, B, length(alpha), tile, panel, batch)
  if (!pointwise.only && m > 1L)
    for (level in alpha) .np_plot_gaussian_draws(B, level, m)
  moments <- .np_udens_d8_pilot_moments(state, grid)
  kind <- RNGkind()
  entry <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  factor <- NULL; maxima <- numeric(); remainder <- 0; delta <- 1e-4
  normal <- qnorm(1-alpha/2); critical <- normal
  if (!pointwise.only && m > 1L) {
    store <- .np_udens_d8_store(m)
    on.exit(.np_udens_d8_store_close(store), add = TRUE)
    factor <- .np_udens_d8_pilot_factor(state, grid, moments, store, budget$panel)
    budget <- .np_udens_d8_pilot_draw_budget(state, m, factor$rank, budget)
    .np_udens_d8_store_for_draws(store, budget$tile)
    normal.stream <- NULL
    if (!is.null(.normal.request)) {
      normal.stream <- .np_plot_normal_stream(.normal.request$seed,
        .normal.request$kind, as.double(factor$rank) * B)
      on.exit(.np_plot_normal_close(normal.stream), add = TRUE)
    }
    maxima <- numeric(B)
    for (first in seq.int(1L, B, by = budget$batch)) {
      jj <- seq.int(first, min(B, first+budget$batch-1L))
      innovations <- matrix(if (is.null(normal.stream)) rnorm(factor$rank * length(jj))
        else .np_plot_normal_read(normal.stream, factor$rank * length(jj)), factor$rank)
      maximum <- numeric(length(jj))
      for (target in seq.int(1L, m, by = budget$tile)) {
        last <- min(m, target+budget$tile-1L)
        block <- .np_udens_d8_store_rows(store, target, last)
        values <- .np_udens_d8_pilot_dispatch(block, innovations)
        if (!identical(dim(values), c(nrow(block), length(jj))) || any(!is.finite(values)))
          stop("invalid D8 pilot Gaussian innovations", call. = FALSE)
        for (row in seq_len(nrow(values))) maximum <- pmax(maximum, abs(values[row, ]))
      }
      maxima[jj] <- maximum
      .np_plot_progress_notify(stage = "Calibrating full-grid D8 Gaussian maximum",
        done = max(jj), total = B, unit = "draws")
    }
    if (!is.null(normal.stream) && normal.stream$remaining != 0)
      stop("incomplete D8 seeded normal consumption", call. = FALSE)
    remainder <- sqrt(2 * max(factor$residual) * log(2*m/delta))
    critical <- pmax(as.vector(quantile(maxima, 1-alpha+delta, type = 7L)) + remainder, normal)
  }
  widths <- list(pointwise = outer(moments$se, normal))
  if (!pointwise.only) {
    widths$simultaneous <- outer(moments$se, critical)
    widths$bonferroni <- pmax(outer(moments$se, qnorm(1-alpha/(2*m))), widths$simultaneous)
  }
  intervals <- lapply(widths, function(width) {
    result <- array(0, c(m, 2L, length(alpha)),
      dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    result[, 1L, ] <- moments$center - width
    result[, 2L, ] <- moments$center + width
    if (any(!is.finite(result))) stop("nonfinite D8 pilot endpoints", call. = FALSE)
    result
  })
  list(center = moments$center, original.center = moments$original.center, se = moments$se,
    covariance = NULL, bands = intervals, intervals = intervals,
    absolute.intervals = intervals, half.width = widths, maximum = maxima,
    critical = critical, remainder.bound = remainder, factor = factor,
    settings = list(construction = state$construction, law = "gaussian", h = state$h,
      n = state$n, q = state$q, alpha = alpha, B = B, grid.size = m,
      budget = budget, covariance = "complete D8 covariance under fixed positive Gaussian pilot",
      pilot.bandwidth = "unchanged fitted bandwidth", delta = delta,
      maximum.statistic = "standardized absolute full-grid Gaussian maximum",
      quantile.type = 7L, family = "one supplied finite grid", available.bands = names(intervals),
      calibration.status = "private provider; public activation separately gated"),
    rng = list(kind = kind, entry = entry,
      end = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL))
}
