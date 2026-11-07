# Joint scale envelope: preserve the F3 center and F4 observed scale, and
# separately account for uncertainty in that scale. beta is deliberately
# fixed across confidence levels; it is not a second coverage guarantee.
.np_udens_scale_rows <- function(state, grid) {
  n <- state$n
  G <- matrix(0, nrow(grid), n)
  scale <- numeric(nrow(grid))
  for (j in seq_len(nrow(grid))) {
    target <- grid[j, , drop = FALSE]
    a <- .np_udens_b2_pilot_moment(state, target)
    b <- .np_udens_b2_pilot_moment(state, target, target)
    aw <- .np_udens_b4_pilot_moment(state, target, tau = 2)
    bw <- .np_udens_b4_pilot_moment(state, target, target, tau = 2)
    v <- (mean(b) - mean(a)^2) / n
    w <- (mean(bw) - mean(aw)^2) / n
    if (any(!is.finite(c(v, w))) || v <= 0 || w <= 0 || n < 2L)
      stop("density scale envelope requires positive finite pilot variances", call. = FALSE)
    iv <- (b - mean(b) - 2*mean(a)*(a - mean(a))) / n
    iw <- (bw - mean(bw) - 2*mean(aw)*(aw - mean(aw))) / n
    influence <- iv/v - iw/(2*w)
    G[j, ] <- (influence - mean(influence)) / sqrt(n*(n-1))
    scale[j] <- v/sqrt(w)
  }
  if (any(!is.finite(G)))
    stop("nonfinite density log-scale influence", call. = FALSE)
  list(factor = G, t = sqrt(rowSums(G^2)), scale = scale)
}

# One-sided reduction of the canonical bounded Gaussian symmetric-root
# construction, generated once from the currently active normal stream.
.np_udens_scale_gaussian <- function(G, B, max.bytes) {
  C <- tcrossprod(G)
  # Preserve the initialized-stream precondition, including for a zero family.
  get(".Random.seed", .GlobalEnv)
  process <- .np_plot_gaussian_maxima(C, alpha = .05, draws = B,
    max.bytes = max.bytes, reduction = "negative")
  end <- get(".Random.seed", .GlobalEnv)
  if (!process$metadata$positive.variance)
    return(list(maxima = rep(0, B), end = end))
  list(maxima = process$maxima, end = end, rank = process$rank)
}

.np_udens_scale_inference <- function(state, grid, law = c("gaussian", "inid", "fixed", "geom"),
    alpha = .05, B = 9999L, tile = 25L,
    bands = c("pointwise", "simultaneous", "bonferroni"), calibration = NULL,
    blocklen = NULL) {
  law <- match.arg(law)
  if (!identical(state$metadata$type, "fixed"))
    stop("joint-scale density bounds require fixed bandwidths", call. = FALSE)
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (existed) get(".Random.seed", .GlobalEnv) else NULL
  success <- FALSE
  on.exit(if (!success) {
    if (existed) assign(".Random.seed", previous, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  m <- nrow(grid); small <- law == "gaussian" && m <= 512L
  budget <- .np_udens_b2_budget(state, m, B, tile, small)
  extra <- 8*(if (small) 4*as.double(m)*state$n else 0) + 8*32*m
  if (!is.null(calibration) && calibration %in% c("separate-tails", "tail-envelope", "percentile-tails"))
    extra <- extra + .np_udens_tail_bytes(B, m, alpha, tile)
  if (identical(calibration, "percentile-tails")) extra <- extra + 48*as.double(m)*length(alpha)
  budget$per.process <- budget$per.process + extra
  budget$bytes <- budget$per.process * budget$processes
  if (budget$bytes > state$max.bytes)
    stop("density scale envelope live-byte forecast exceeds the requested budget", call. = FALSE)
  dependent <- law %in% c("fixed", "geom")
  cache <- NULL
  if (dependent) {
    request <- .np_udens_block_request(state$x, law, blocklen)
    blocklen <- request$blocklen
    if (!existed) runif(1L)
    cache <- .np_udens_block_cache(state$n, B, m, tile, law, blocklen,
                                  budget, state$max.bytes)
    if (!is.null(cache)) {
      budget$per.process <- budget$per.process + cache$plan$extra.bytes
      budget$bytes <- cache$plan$bytes
      budget$innovation.cache <- cache$plan
    }
  }
  # The location and scale passes replay even for a single target tile.
  .np_bias_replay_check(law, blocklen, !small,
    retained = !is.null(cache))
  base <- .np_udens_reference_inference(state, grid, law, alpha, B, tile, bands, calibration,
                                        blocklen = blocklen, process.cache = cache)
  if (dependent) base$rng$entry <- previous
  start <- if (small) base$rng$end else base$rng$start
  point <- tau <- numeric(m); maximum <- rep(-Inf, B); end <- NULL
  for (first in seq.int(1L, m, by = if (small) m else tile)) {
    ii <- seq.int(first, min(m, first + if (small) m-1L else tile-1L))
    z <- .np_udens_scale_rows(state, grid[ii, , drop = FALSE])
    if (any(abs(z$scale - base$se[ii]) > 2e-12*pmax(1, base$se[ii])))
      stop("density scale envelope changed the F4 scale", call. = FALSE)
    tau[ii] <- z$t
    if (small) {
      # The reference producer left this stream live at base$rng$end.
      draw <- .np_udens_scale_gaussian(z$factor, B, state$max.bytes)
      maximum <- draw$maxima; end <- draw$end
      point[ii] <- ifelse(z$t > 0, qnorm(.95), 0)
    } else {
      draw <- .np_bias_stream(start,
        if (dependent) .np_udens_block_draws(z$factor, B, law, blocklen, cache)
        else .np_bias_process_draws(z$factor, B, law, NULL))
      if (!is.null(end) && !identical(end, draw$end))
        stop("density scale tiles consumed inconsistent streams", call. = FALSE)
      end <- draw$end
      divisor <- z$t; zero <- which(divisor == 0)
      if (length(zero) && any(draw$value[zero, , drop = FALSE] != 0))
        stop("zero log-scale variance requires an exact zero process", call. = FALSE)
      divisor[zero] <- 1 # exact deterministic zero, never a variance floor
      signed <- -draw$value/divisor
      if (any(!is.finite(signed))) stop("nonfinite density log-scale process", call. = FALSE)
      maximum <- pmax(maximum, apply(signed, 2L, max))
      point[ii] <- if (law == "gaussian") ifelse(z$t > 0, qnorm(.95), 0) else
        apply(signed, 1L, quantile, probs = .95, type = 1L, names = FALSE)
    }
    .np_plot_progress_notify(stage = "Calibrating density scale uncertainty",
      done = max(ii), total = m, unit = "targets")
  }
  raw <- unname(quantile(maximum, .95, type = 1L))
  joint <- if (law == "gaussian") min(qnorm(1-.05/m), max(qnorm(.95), raw)) else max(raw, point)
  if (m == 1L) joint <- point[1L]
  inflation <- list(pointwise = exp(pmax(0, point)*tau),
                    simultaneous = exp(max(0, joint)*tau))
  if (any(!is.finite(unlist(inflation))))
    stop("density joint scale envelope is unbounded", call. = FALSE)
  base$base.critical <- base$critical
  base$base.raw.critical <- base$raw.critical
  base$scale.envelope <- list(beta = .05, t = tau, point.critical = point,
    joint.critical = joint, joint.raw = raw, inflation = inflation,
    rng = list(start = start, end = end), reference.rng = base$rng)
  base$critical$pointwise <- base$critical$pointwise * inflation$pointwise
  base$critical$simultaneous <- pmax(base$critical$simultaneous * inflation$simultaneous,
                                    base$critical$pointwise)
  if (!is.null(base$critical$bonferroni))
    base$critical$bonferroni <- pmax(base$critical$bonferroni * inflation$simultaneous,
                                   base$critical$simultaneous)
  base$bands <- lapply(base$critical, function(k) {
    out <- array(0, c(m, 2L, length(alpha)),
      dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    out[, 1L, ] <- base$center - base$se*k
    out[, 2L, ] <- base$center + base$se*k
    if (any(!is.finite(out))) stop("nonfinite density scale-envelope endpoints", call. = FALSE)
    out
  })
  base$intervals <- base$absolute.intervals <- base$bands
  if (!is.null(base[["tail", exact = TRUE]])) {
    base$location.raw.intervals <- base$raw.intervals
    endpoint.roots <- if (identical(calibration, "percentile-tails"))
      .np_band_percentile_roots(base$tail$limits) else base$tail$limits
    base$raw.intervals <- .np_band_scale_union(base$center, base$se, endpoint.roots,
      list(pointwise = inflation$pointwise, simultaneous = inflation$simultaneous,
           bonferroni = inflation$simultaneous))
    base$raw.intervals$simultaneous <- .np_band_tail_protect(
      list(band = base$raw.intervals$simultaneous),
      list(band = base$raw.intervals$pointwise))$limits$band
    base$raw.intervals$bonferroni <- .np_band_tail_protect(
      list(band = base$raw.intervals$bonferroni),
      list(band = base$raw.intervals$simultaneous))$limits$band
    base$bands <- if (calibration == "tail-envelope")
      .np_band_tail_protect(base$raw.intervals, base$absolute.intervals)$limits else base$raw.intervals
    base$intervals <- base$bands
  }
  base$settings$construction <- "density-joint-scale-envelope-S-v1"
  base$settings$budget <- budget
  base$settings$calibration.status <- "bounded full-range qualification; conservative and local limitations documented"
  if (law != "gaussian")
    base$resolution <- rbind(base$resolution, .np_udens_resolution_scale(B, m, tau))
  if (!is.null(base[["tail", exact = TRUE]]))
    base$settings$calibration.status <- "private asymmetric candidate; coverage unqualified"
  if (dependent) {
    base$settings$block.selection <- request
    base$settings$calibration.status <- "private block-scale transport; coverage unqualified"
  }
  base$rng$end <- end
  assign(".Random.seed", end, .GlobalEnv)
  success <- TRUE
  base
}
