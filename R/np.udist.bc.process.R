# Thin CDF coordinator; existing generic calibration and literal inid transport.
.np_udist_bc_inference <- function(state, grid, law = "gaussian", alpha = c(.05, .1),
                                    B = 9999L, tile = 25L,
                                    bands = c("pointwise", "simultaneous", "bonferroni"),
                                    corrected = TRUE, blocklen = NULL) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  entry <- if (existed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind(); success <- FALSE
  on.exit({if (!success) {
    do.call(RNGkind, as.list(kind))
    if (existed) assign(".Random.seed", entry, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }}, add = TRUE)
  law <- match.arg(law, c("gaussian", "inid", "fixed", "geom"))
  if (law == "inid" && !is.null(blocklen) &&
      (!is.numeric(blocklen) || length(blocklen) != 1L || !is.finite(blocklen) || blocklen != 1))
    stop("inid CDF inference requires block length one", call. = FALSE)
  block.request <- if (law %in% c("fixed", "geom")) .np_udist_bc_block_request(state$x, law, blocklen) else NULL
  bands <- match.arg(bands, c("pointwise", "simultaneous", "bonferroni"), several.ok = TRUE)
  if (!is.numeric(alpha) || !length(alpha) || any(!is.finite(alpha)) ||
      any(alpha <= 0 | alpha >= .5) || anyDuplicated(alpha)) stop("invalid CDF alpha", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "CDF draws", .Machine$integer.max)
  if (B < 2L) stop("CDF draws must be at least two", call. = FALSE)
  .np_udist_bc_grid(state, grid); m <- nrow(grid)
  budget <- .np_udist_bc_budget(state, m, B, tile)
  stochastic <- law != "gaussian" || "simultaneous" %in% bands
  if (stochastic && !existed) runif(1L)
  start <- if (stochastic) get(".Random.seed", .GlobalEnv, inherits = FALSE) else entry
  .np_bias_normal_replay_check(law,
    if (is.null(block.request)) NULL else block.request$length,
    stochastic && m > tile && !(law == "gaussian" && m <= 32L), seed = start)
  end <- start; covariance <- NULL; raw.max <- NULL
  degenerate <- cancellation <- integer(m)
  if (law == "gaussian" && m <= 32L) {
    z <- .np_udist_bc_covariance(state, grid, tile, corrected)
    center <- z$center; original <- z$original.center; se <- z$se; covariance <- z$covariance
    if (any(!is.finite(se)) || any(se <= 0)) stop("CDF inference requires positive observed scales", call. = FALSE)
    point <- matrix(rep(qnorm(alpha/2, lower.tail = FALSE), each = m), m)
    bonf <- matrix(rep(qnorm(alpha/(2*m), lower.tail = FALSE), each = m), m)
    joint <- point
    if (stochastic) {
      draws <- .np_plot_gaussian_calibrate_alphas(covariance,
        alpha = alpha, draws = B, max.bytes = state$max.bytes)
      end <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
      joint <- matrix(rep(vapply(draws, function(d) d$critical, numeric(1L)), each = m), m)
      raw.max <- vapply(draws, function(d) d$critical.raw, numeric(1L))
    }
    critical <- list(pointwise = point, simultaneous = joint, bonferroni = bonf)
    raw.bonf <- bonf
  } else {
    limits <- .np_bias_limits_prepare(rep("cdf", m), alpha, law == "gaussian", B)
    center <- original <- se <- numeric(m); end <- if (stochastic) NULL else start
    for (block in .np_bias_target_blocks(m, 1L, anchors = tile)) {
      ii <- block$grid
      z <- .np_udist_bc_rows(state, grid[ii, , drop = FALSE], corrected)
      scale <- if (is.null(block.request)) z$se else sqrt(.np_bias_score_variance(z$factor, block.request$length, law))
      if (any(!is.finite(scale)) || any(scale <= 0)) stop("CDF inference requires positive observed scales", call. = FALSE)
      if (law == "gaussian") {
        if (stochastic) {
          stream <- .np_bias_stream(start,
            .np_bias_process_draws(z$factor, B, "gaussian", NULL))
          if (!is.null(end) && !identical(end, stream$end))
            stop("CDF tile RNG mismatch", call. = FALSE)
          end <- stream$end
          limits <- .np_bias_limits_update(limits, abs(stream$value / scale), ii)
        }
        d <- list(degeneracy = integer(length(ii)), cancellation = integer(length(ii)))
      } else {
        stream <- .np_bias_stream(start, if (!identical(state$type, "fixed") && state$q > 0L)
          .np_udist_bc_nn_dispatch(state, grid[ii, , drop = FALSE], z, B, law,
            if (is.null(block.request)) 1L else block.request$length, corrected) else if (is.null(block.request))
          .np_udens_bc_inid_dispatch(z$factor, z$corrected, B, accurate = TRUE) else
          .np_udist_bc_block_dispatch(z$factor, z$corrected, B, law, block.request$length))
        if (!is.null(end) && !identical(end, stream$end)) stop("CDF tile RNG mismatch", call. = FALSE)
        end <- stream$end; d <- stream$value
        limits <- .np_udens_bc_limits_update(limits, d$absolute, ii, d$degeneracy)
      }
      center[ii] <- z$center; original[ii] <- z$original.center; se[ii] <- scale
      degenerate[ii] <- d$degeneracy; cancellation[ii] <- d$cancellation
    }
    finished <- .np_bias_limits_finish(limits)
    critical <- lapply(finished$limits, function(z) matrix(z[, 2L, ], m, length(alpha)))
    raw.bonf <- if (law == "gaussian") critical$bonferroni else limits$bonf
    raw.max <- if (stochastic) drop(finished$raw_max_critical) else NULL
  }
  intervals <- lapply(critical[bands], function(k) {
    out <- array(0, c(m, 2L, length(alpha)), dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    out[, 1L, ] <- center - se*k; out[, 2L, ] <- center + se*k
    out
  })
  if (stochastic) assign(".Random.seed", end, .GlobalEnv)
  success <- TRUE
  list(center = center, original.center = original, se = se, covariance = covariance,
    intervals = intervals, critical = critical[bands], raw.max = raw.max, raw.bonf = raw.bonf,
    envelope = critical$bonferroni > raw.bonf, degeneracy = degenerate, cancellation = cancellation,
    settings = list(construction = state$construction, bwtype = state$type, h = state$h,
      law = law, block = block.request, n = state$n, alpha = alpha, B = B, corrected = corrected, budget = budget,
      representation = if (law == "gaussian" && m <= 32L) "bounded covariance" else "bounded observation factors",
      metadata = state$metadata, target.ids = rownames(grid), scope = "declared finite family"),
    rng = list(kind = kind, entry = entry, start = start, end = end))
}
