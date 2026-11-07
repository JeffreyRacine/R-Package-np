# CDF absolute-process calibration. Critical values are in CDF units.
# The complete centered observation factor is shared across every target tile.
.np_udist_joint_inference <- function(state, grid, law = "gaussian", alpha = .05,
                                     B = 9999L, tile = 25L,
                                     bands = c("pointwise", "simultaneous", "bonferroni"),
                                     corrected = TRUE, blocklen = NULL, selection = NULL,
                                     calibration = NULL) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  entry <- if (existed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind(); success <- FALSE
  on.exit(if (!success) {
    do.call(RNGkind, as.list(kind))
    if (existed) assign(".Random.seed", entry, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  law <- match.arg(law, c("gaussian", "inid", "fixed", "geom"))
  if (identical(calibration, "default")) calibration <- NULL
  if (!is.null(calibration) && (!is.character(calibration) || length(calibration) != 1L ||
      is.na(calibration) || !calibration %in% c("absolute", "separate-tails", "tail-envelope", "score-envelope")))
    stop("invalid CDF calibration", call. = FALSE)
  signed.calibration <- !is.null(calibration) && calibration %in% c("separate-tails", "tail-envelope")
  score.calibration <- identical(calibration, "score-envelope")
  if (!is.null(calibration) && ((!score.calibration && law == "gaussian") ||
      (score.calibration && (law != "gaussian" || !isTRUE(corrected) ||
       !identical(state$type, "fixed") || !state$q %in% 1:2 || state$s != 0L))))
    stop("CDF calibration is unavailable for the selected law or contribution geometry", call. = FALSE)
  if (signed.calibration && (!isTRUE(corrected) || !state$q %in% 0:2 || !state$s %in% 0:2))
    stop("signed CDF calibration requires the admitted corrected contribution geometry", call. = FALSE)
  if (signed.calibration && !identical(state$type, "fixed") && state$q > 0L)
    stop("asymmetric numeric nearest-neighbor CDF calibration remains unqualified", call. = FALSE)
  dependent <- law %in% c("fixed", "geom")
  if (dependent) {
    if (is.null(selection)) selection <- .np_udist_bc_block_request(state$x, law, blocklen)
    if (!is.list(selection) || !identical(selection$law, law))
      stop("invalid CDF block selection", call. = FALSE)
    .np_udist_bc_block_request(state$x, law, selection$length)
    blocklen <- selection$length
  } else if (!is.null(selection) || !is.null(blocklen))
    stop("CDF block controls require a dependent bootstrap law", call. = FALSE)
  if (!state$type %in% c("fixed", "generalized_nn", "adaptive_nn") ||
      (identical(state$type, "adaptive_nn") && state$q > 0L))
    stop("absolute CDF process has no qualified numeric ANN construction", call. = FALSE)
  nearest <- !identical(state$type, "fixed")
  if (nearest && !law %in% c("gaussian", "inid", "fixed", "geom"))
    stop("selected nearest-neighbor CDF bootstrap law remains unqualified", call. = FALSE)
  bands <- match.arg(bands, c("pointwise", "simultaneous", "bonferroni"), several.ok = TRUE)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= .5) || anyDuplicated(alpha))
    stop("invalid CDF alpha", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "CDF draws", .Machine$integer.max)
  if (B < 2L && nearest && state$q > 0L && law != "gaussian")
    stop("numeric GNN bootstrap empirical standard deviation requires at least two CDF draws", call. = FALSE)
  .np_udist_bc_grid(state, grid)
  m <- nrow(grid); budget <- .np_udist_bc_budget(state, m, B, tile)
  if (signed.calibration || score.calibration) {
    extra <- 192*as.double(m)*length(alpha)
    if (signed.calibration)
      extra <- extra + 96*as.double(B) + 48*as.double(m)*length(alpha) +
        8*min(tile, m)*as.double(B)
    budget$per.process <- budget$per.process + extra
    budget$bytes <- budget$per.process * budget$processes
    if (!is.finite(budget$bytes) || budget$bytes > state$max.bytes)
      stop("CDF asymmetric workspace exceeds the requested budget", call. = FALSE)
    if (signed.calibration) {
      map <- data.frame(row = seq_len(m), target = "cdf")
      tails <- .np_band_tail_prepare(B, map, alpha, byte.limit = extra)
    }
  }
  stochastic <- law != "gaussian" || any(bands != "pointwise")
  if (stochastic && !existed) runif(1L)
  start <- if (stochastic) get(".Random.seed", .GlobalEnv, inherits = FALSE) else entry
  end <- NULL
  cached <- NULL; cache.plan <- NULL
  if (signed.calibration && dependent) {
    cache.plan <- .np_udist_cache_plan(state$n, B, m, tile, blocklen,
      budget, state$max.bytes, length(start))
    if (!is.null(cache.plan)) {
      cached <- .np_bias_stream(start,
        .np_udist_cache_prepare(state$n, B, law, blocklen))$value
      budget$per.process <- budget$per.process + cache.plan$extra.bytes
      budget$bytes <- cache.plan$bytes
      budget$innovation.cache <- cache.plan
    } else if (m > tile && .np_plot_replay_uniform_state_absent(kind)) {
      # Multiple target tiles must share the same observation-count draws.
      # A header-only custom seed cannot replay uncached dependent draws.
      stop("unsafe RNG replay: resolved signed CDF calibration requires custom uniform state absent from .Random.seed when its innovation cache is unavailable; explicitly select a qualified built-in generator before plotting", call. = FALSE)
    }
  }
  center <- original <- se <- numeric(m)
  limits <- .np_bias_limits_prepare(rep("cdf", m), alpha, law == "gaussian", B)
  gnn.bootstrap <- nearest && state$q > 0L && law != "gaussian"
  if (!gnn.bootstrap && !signed.calibration && identical(law, "fixed") &&
      m > tile && .np_plot_replay_uniform_state_absent(kind))
    stop("unsafe RNG replay: resolved unsigned fixed-block CDF calibration requires custom uniform state absent from .Random.seed; explicitly select a qualified built-in generator before plotting", call. = FALSE)
  if (!gnn.bootstrap)
    .np_bias_replay_check(law, blocklen, stochastic && m > tile,
      retained = !is.null(cached), seed = start)
  if (gnn.bootstrap) {
    exact <- .np_udist_gnn_stream(state, grid, law, alpha, B, corrected, blocklen, start)
    center <- exact$center; original <- exact$original; se <- exact$se
    limits <- exact$limits; end <- exact$end; budget <- exact$budget
  } else for (block in .np_bias_target_blocks(m, 1L, anchors = tile)) {
    ii <- block$grid
    rows <- .np_udist_bc_rows(state, grid[ii, , drop = FALSE], corrected)
    center[ii] <- rows$center; original[ii] <- rows$original.center
    se[ii] <- if (dependent)
      sqrt(.np_bias_score_variance(rows$factor, blocklen, law)) else rows$se
    if (stochastic) {
      stream <- .np_bias_stream(start,
        if (is.null(cached)) .np_bias_process_draws(rows$factor, B, law, blocklen)
        else .np_udist_cache_draws(rows$factor, cached, B, law, blocklen))
      if (!is.null(end) && !identical(end, stream$end))
        stop("absolute CDF tile RNG mismatch", call. = FALSE)
      end <- stream$end
      limits <- .np_bias_limits_update(limits, abs(stream$value), ii)
      if (signed.calibration)
        tails <- .np_band_tail_update(tails, stream$value,
          map[ii, , drop = FALSE], seq_len(B))
    }
  }
  if (any(!is.finite(se)) || !any(se > 0))
    stop("absolute CDF process unavailable: no positive observed variance", call. = FALSE)
  raw <- if (stochastic) unname(quantile(limits$maximum, 1-alpha, type = 1L)) else NULL
  if (law == "gaussian") {
    point <- outer(se, qnorm(alpha/2, lower.tail = FALSE))
    bonf <- outer(se, qnorm(alpha/(2*m), lower.tail = FALSE))
    joint <- if (stochastic) pmin(apply(bonf, 2L, max),
      pmax(apply(point, 2L, max), raw)) else apply(point, 2L, max)
  }
  if (law != "gaussian") {
    point <- limits$point; bonf <- limits$bonf; joint <- raw
  }
  joint <- matrix(rep(joint, each = m), m, length(alpha))
  raw.bonf <- bonf; bonf <- pmax(bonf, joint)
  widths <- list(pointwise = point, simultaneous = joint, bonferroni = bonf)
  intervals <- lapply(widths[bands], function(w) {
    out <- array(0, c(m, 2L, length(alpha)),
      dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    out[, 1L, ] <- center-w; out[, 2L, ] <- center+w
    out
  })
  if (stochastic) assign(".Random.seed", end, .GlobalEnv) else end <- entry
  answer <- list(center = center, original.center = original, se = se,
    intervals = intervals, critical = widths[bands], raw.max = raw,
    raw.bonf = raw.bonf, envelope = bonf > raw.bonf,
    process.maximum = if (stochastic) drop(limits$maximum) else NULL,
    zero.marginals = which(se == 0), degeneracy = integer(m), cancellation = integer(m),
    settings = list(construction = "cdf-absolute-process-v1", critical.units = "CDF",
      bwtype = state$type, h = state$h, law = law, n = state$n, alpha = alpha,
      B = B, corrected = corrected, budget = budget, selection = selection,
      representation = "bounded observation factors", metadata = state$metadata,
      target.ids = rownames(grid), scope = "declared finite family"),
    rng = list(kind = kind, entry = entry, start = start, end = end))
  if (nearest) {
    answer$settings$h <- NULL
    answer$settings$counts <- state$requested
    answer$settings$radius <- "sample-specific query kth distance"
    answer$settings$calibration <- if (gnn.bootstrap)
      "literal exact GNN center displacement" else "complete GNN influence process"
    if (gnn.bootstrap) answer$settings$scale <- "selected-law empirical SD of exact center displacements"
  }
  if (isTRUE(corrected) && identical(state$type, "fixed") &&
      state$q %in% 1:2 && state$s == 0L && law %in% c("gaussian", "inid"))
    answer <- .np_udist_marginal_envelope(answer, widths, state$n, state$q, alpha)
  if (signed.calibration || score.calibration) {
    tail.result <- if (signed.calibration) .np_band_tail_finish(tails) else NULL
    score <- isTRUE(corrected) && identical(state$type, "fixed") &&
      state$q %in% 1:2 && state$s == 0L && law %in% c("gaussian", "inid")
    answer <- .np_udist_asymmetric_finish(answer, widths, tail.result, calibration,
      state$n, state$q, alpha, score)
  } else if (!is.null(calibration)) {
    answer$settings$calibration.requested <- calibration
    answer$settings$calibration.method <- calibration
  }
  empirical.bands <- if (law == "gaussian") {
    if (m > 1L) intersect(bands, "simultaneous") else character()
  } else bands
  if (length(empirical.bands))
    answer$resolution <- .np_band_resolution(B, rep("cdf", m), alpha,
      se > 0, empirical.bands)
  if (signed.calibration) {
    answer$absolute.resolution <- answer$resolution
    answer$resolution <- .np_band_resolution(B, rep("cdf", m), alpha, se > 0,
      bands, signed = TRUE, simultaneous.index = matrix(1 + floor(B*alpha/2), 1L))
    answer$resolution$calibration <- "separate-tails"
  }
  success <- TRUE
  answer
}
