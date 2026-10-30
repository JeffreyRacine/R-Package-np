# Union of root-inverted intervals over a positive scale range [s, s*I].
# Max/min root displacements preserve one-sided roots without zero clamping.
.np_band_scale_union <- function(center, scale, limits, inflation) {
  if (!identical(names(limits), names(inflation)) || !length(limits) ||
      !is.numeric(scale) || any(!is.finite(scale)) || any(scale <= 0))
    stop("invalid scale-union bands or scales", call. = FALSE)
  protected <- limits
  for (band in names(limits)) {
    x <- limits[[band]]; I <- inflation[[band]]
    if (!is.numeric(I) || length(I) != length(scale) || any(!is.finite(I)) ||
        any(I < 1) || length(dim(x)) != 3L || dim(x)[1L] != length(scale) ||
        dim(x)[2L] != 2L || any(!is.finite(x)))
      stop("invalid scale-union inflation or roots", call. = FALSE)
    for (a in seq_len(dim(x)[3L])) {
      if (any(x[, 1L, a] > x[, 2L, a]))
        stop("reversed scale-union roots", call. = FALSE)
      protected[[band]][, 1L, a] <- pmin(x[, 1L, a], I*x[, 1L, a])
      protected[[band]][, 2L, a] <- pmax(x[, 2L, a], I*x[, 2L, a])
    }
  }
  out <- .np_band_rank_invert(center, scale, protected)
  if (any(!is.finite(unlist(out)))) stop("nonfinite scale-union endpoints", call. = FALSE)
  out
}

# Re-express percentile endpoints through the canonical inverse/union owners.
# This changes the confidence construction, not the signed process orientation.
.np_band_percentile_roots <- function(limits) {
  lapply(limits, function(x) {
    if (length(dim(x)) != 3L || dim(x)[2L] != 2L ||
        any(!is.finite(x)) || any(x[, 1L, ] > x[, 2L, ]))
      stop("invalid percentile root endpoints", call. = FALSE)
    out <- x
    out[, 1L, ] <- -x[, 2L, ]
    out[, 2L, ] <- -x[, 1L, ]
    out
  })
}
