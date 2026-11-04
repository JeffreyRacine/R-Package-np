# Default selection is density-local. Explicit controls retain their versions.
.np_udens_d8_default <- function(bws, control, center, method, alpha = .05) {
  is.null(control) && identical(center, "bias-corrected") &&
    identical(method, "asymptotic") && inherits(bws, "bandwidth") &&
    identical(bws[["type", exact = TRUE]], "fixed") &&
    identical(bws[["ckertype", exact = TRUE]], "gaussian") &&
    isTRUE(bws[["ckerorder", exact = TRUE]] == 2L) &&
    identical(bws[["ckerbound", exact = TRUE]], "none") &&
    isTRUE(bws[["ndim", exact = TRUE]] %in% 1:3) &&
    isTRUE(bws[["ncon", exact = TRUE]] == bws[["ndim", exact = TRUE]]) &&
    is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
    alpha > 1e-4 && alpha < .5
}

.np_udens_d8_active <- function(state) {
  identical(state[["construction", exact = TRUE]], "density-d8-positive-pilot-v2")
}

# A multi-panel plot is one inferential family: calibrate the supplied union
# once, then select rows without recalibration or independent panel innovations.
.np_udens_d8_plot_subset <- function(raw, rows) {
  raw$boot.err <- raw$boot.err[rows, , drop = FALSE]
  if (!is.null(raw$boot.all.err))
    raw$boot.all.err <- lapply(raw$boot.all.err, function(z) z[rows, , drop = FALSE])
  z <- raw$inference
  for (key in c("center", "original.center", "se")) z[[key]] <- z[[key]][rows]
  for (key in c("bands", "intervals", "absolute.intervals"))
    z[[key]] <- lapply(z[[key]], function(a) a[rows, , , drop = FALSE])
  z$half.width <- lapply(z$half.width, function(a) a[rows, , drop = FALSE])
  z$settings$family.rows <- rows
  z$settings$family <- "union of all supplied plot slices"
  raw$inference <- z
  raw
}
