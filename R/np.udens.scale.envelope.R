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
# construction. The shared two-sided calibrator is intentionally unchanged.
.np_udens_scale_gaussian <- function(G, B, max.bytes) {
  C <- tcrossprod(G)
  # Validate the complete covariance using the existing bounded calibrator,
  # replaying its stream for the signed reduction below.
  start <- get(".Random.seed", .GlobalEnv)
  canonical <- .np_bias_stream(start,
    .np_plot_gaussian_calibrate(C, alpha = .05, draws = B, max.bytes = max.bytes))
  active <- which(diag(C) > 0)
  if (!length(active)) return(list(maxima = rep(0, B), end = canonical$end))
  sdv <- sqrt(diag(C)[active])
  R <- C[active, active, drop = FALSE] / outer(sdv, sdv)
  R <- (R + t(R))/2
  eig <- eigen(R, symmetric = TRUE)
  root <- tcrossprod(sweep(eig$vectors, 2L, sqrt(pmax(eig$values, 0)), "*"), eig$vectors)
  m <- nrow(C)
  block <- min(4096L, B, floor((max.bytes - 64*as.double(m)^2 - 24*B)/(64*m)))
  if (block < 1L) stop("density scale covariance exceeds memory budget", call. = FALSE)
  maximum <- numeric(B)
  for (first in seq.int(1L, B, by = block)) {
    last <- min(B, first + block - 1L)
    Z <- root %*% matrix(rnorm(length(active)*(last-first+1L)), length(active))
    maximum[first:last] <- apply(-Z, 2L, max)
    if (length(active) < m) maximum[first:last] <- pmax(0, maximum[first:last])
  }
  end <- get(".Random.seed", .GlobalEnv)
  if (!identical(end, canonical$end))
    stop("density signed Gaussian replay consumed inconsistent stream", call. = FALSE)
  list(maxima = maximum, end = end, rank = canonical$value$rank)
}

.np_udens_scale_inference <- function(state, grid, law = c("gaussian", "inid"),
    alpha = .05, B = 9999L, tile = 25L,
    bands = c("pointwise", "simultaneous", "bonferroni")) {
  law <- match.arg(law)
  if (!state$q %in% 1:2 || state$r != 0L || state$s != 0L || !identical(state$metadata$type, "fixed"))
    stop("joint-scale density bounds require fixed bandwidths and one or two numeric variables", call. = FALSE)
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
  budget$per.process <- budget$per.process + extra
  budget$bytes <- budget$per.process * budget$processes
  if (budget$bytes > state$max.bytes)
    stop("density scale envelope live-byte forecast exceeds the requested budget", call. = FALSE)
  base <- .np_udens_reference_inference(state, grid, law, alpha, B, tile, bands)
  start <- if (small) base$rng$end else base$rng$start
  point <- tau <- numeric(m); maximum <- rep(-Inf, B); end <- NULL
  for (first in seq.int(1L, m, by = if (small) m else tile)) {
    ii <- seq.int(first, min(m, first + if (small) m-1L else tile-1L))
    z <- .np_udens_scale_rows(state, grid[ii, , drop = FALSE])
    if (any(abs(z$scale - base$se[ii]) > 2e-12*pmax(1, base$se[ii])))
      stop("density scale envelope changed the F4 scale", call. = FALSE)
    tau[ii] <- z$t
    if (small) {
      draw <- .np_bias_stream(start, .np_udens_scale_gaussian(z$factor, B, state$max.bytes))
      maximum <- draw$value$maxima; end <- draw$end
      point[ii] <- ifelse(z$t > 0, qnorm(.95), 0)
    } else {
      draw <- .np_bias_stream(start, .np_bias_process_draws(z$factor, B, law, NULL))
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
  base$settings$construction <- "density-joint-scale-envelope-S-v1"
  base$settings$budget <- budget
  base$settings$calibration.status <- "bounded full-range qualification; conservative and local limitations documented"
  base$rng$end <- end
  assign(".Random.seed", end, .GlobalEnv)
  success <- TRUE
  base
}
