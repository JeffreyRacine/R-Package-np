# Side-specific outward composition in CDF units. Gaussian process symmetry
# is unchanged; any asymmetry comes from signed empirical tails or the score.
.np_udist_side_envelope <- function(intervals, center, n, q, alpha, score) {
  required <- c("pointwise", "simultaneous", "bonferroni")
  m <- length(center)
  if (!identical(names(intervals), required) || !is.numeric(center) ||
      any(!is.finite(center)) || !is.logical(score) || length(score) != 1L || is.na(score))
    stop("invalid side-specific CDF envelope", call. = FALSE)
  for (z in intervals) {
    if (!identical(dim(z), c(m, 2L, length(alpha))) || any(!is.finite(z)) ||
        any(z[, 1L, ] > z[, 2L, ]))
      stop("invalid CDF envelope endpoints", call. = FALSE)
  }
  out <- intervals
  expansion <- matrix(0, length(alpha), 2L,
    dimnames = list(as.character(alpha), c("lower", "upper")))
  for (a in seq_along(alpha)) {
    point <- matrix(intervals$pointwise[, , a], m, 2L)
    if (score) {
      bound <- .np_udist_marginal_score(center, n, q, alpha[a])
      point[, 1L] <- pmin(point[, 1L], bound[, 1L])
      point[, 2L] <- pmax(point[, 2L], bound[, 2L])
    }
    lower <- max(center-intervals$simultaneous[, 1L, a], center-point[, 1L])
    upper <- max(intervals$simultaneous[, 2L, a]-center, point[, 2L]-center)
    simultaneous <- cbind(center-lower, center+upper)
    bonferroni <- cbind(pmin(intervals$bonferroni[, 1L, a], simultaneous[, 1L]),
                        pmax(intervals$bonferroni[, 2L, a], simultaneous[, 2L]))
    if (m == 1L) simultaneous <- bonferroni <- point
    expansion[a, ] <- c(max(intervals$simultaneous[, 1L, a]-simultaneous[, 1L]),
                         max(simultaneous[, 2L]-intervals$simultaneous[, 2L, a]))
    out$pointwise[, , a] <- point
    out$simultaneous[, , a] <- simultaneous
    out$bonferroni[, , a] <- bonferroni
  }
  list(intervals = out, expansion = expansion)
}

.np_udist_asymmetric_finish <- function(answer, widths, tails, calibration,
                                         n, q, alpha, score) {
  m <- length(answer$center)
  absolute <- answer$intervals
  available <- names(widths)
  process <- lapply(widths, function(w) {
    z <- array(0, c(m, 2L, length(alpha)),
      dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    z[, 1L, ] <- answer$center-w
    z[, 2L, ] <- answer$center+w
    z
  })
  if (!is.null(tails)) {
    for (band in available) dimnames(tails$limits[[band]]) <- dimnames(process[[band]])
    process <- .np_band_rank_invert(answer$center, rep(1, m), tails$limits)
    answer$tail <- tails
    answer$tail$units <- "CDF"
  }
  composed <- .np_udist_side_envelope(process, answer$center, n, q, alpha, score)
  raw <- composed$intervals[names(absolute)]
  answer$absolute.intervals <- absolute
  answer$raw.process.intervals <- process[names(absolute)]
  answer$raw.intervals <- raw
  answer$intervals <- if (identical(calibration, "tail-envelope"))
    .np_band_tail_protect(raw, absolute)$limits else raw
  answer$critical <- lapply(answer$intervals, function(z)
    matrix(pmax(answer$center-z[, 1L, ], z[, 2L, ]-answer$center), m, length(alpha)))
  answer$settings$construction <- paste0("cdf-", calibration, "-v1")
  answer$settings$calibration.requested <- calibration
  answer$settings$calibration.method <- calibration
  answer$settings$side.nesting.expansion <- composed$expansion
  if (score) {
    answer$settings$marginal.calibration <- "side-specific signed-contribution score envelope"
    answer$settings$nesting.expansion <- apply(composed$expansion, 1L, max)
  }
  answer
}

# Reject only the new explicit CDF score choice at unsupported public consumers.
.np_udist_score_scope <- function(control, supported) {
  if (!supported && inherits(control, "np_udist_bias_control") &&
      identical(control[["calibration", exact = TRUE]], "score-envelope"))
    stop("CDF score-envelope requires supported distribution inference", call. = FALSE)
  invisible(NULL)
}
