# Requested-only paired empirical-mass covariance for the realized conditional
# LP fit. Endpoint construction and published point results remain with callers.
.np_conditional_lp_pair_se <- function(bws, txdat, tydat, upper, lower,
                                       upper.y, lower.y = upper.y, cdf = FALSE,
                                       quantile = FALSE) {
  m <- nrow(upper)
  if (nrow(lower) != m || nrow(upper.y) != m || nrow(lower.y) != m)
    stop("internal conditional LP covariance endpoint mismatch", call. = FALSE)
  out <- matrix(NA_real_, m, 6L)
  if (!m) return(out)
  xbw <- .npcdhat_make_xbw(bws, txdat)
  ybw <- .npcdhat_make_ybw(bws, tydat)
  spec <- npConditionalRegEngineSpec(xbw, where = "conditional LP covariance",
                                    ncon.field = "ncon")
  z <- as.matrix(W.lp(txdat[, xbw$icon, drop = FALSE],
    degree = spec$degree.engine, basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.engine))
  storage.mode(z) <- "double"
  # At most eight endpoint pairs per tile; additionally cap the kernel rows
  # near one MiB each. No observation-by-all-evaluation result is retained.
  tile <- max(1L, min(8L, floor(131072 / nrow(txdat))))
  for (first in seq.int(1L, m, by = tile)) {
    rows <- seq.int(first, min(m, first + tile - 1L))
    ex <- rbind(upper[rows, , drop = FALSE], lower[rows, , drop = FALSE])
    ey <- rbind(upper.y[rows, , drop = FALSE], lower.y[rows, , drop = FALSE])
    d <- as.matrix(W.lp(txdat[, xbw$icon, drop = FALSE],
      exdat = ex[, xbw$icon, drop = FALSE], degree = spec$degree.engine,
      basis = spec$basis.engine, bernstein.basis = spec$bernstein.engine))
    storage.mode(d) <- "double"
    w <- .np_kernel_weights_direct(xbw, txdat, exdat = ex,
      bandwidth.divide = TRUE,
      int.do.tree = .npreg_fit_tree_code(xbw, xbw$ncon, xbw$nuno + xbw$nord))
    r <- t(.npcdhat_make_kernel_matrix(ybw, tydat, ey,
      operator = if (cdf) "integral" else "normal"))
    slope <- if (quantile) t(.npcdhat_make_kernel_matrix(
      ybw, tydat, ey, operator = "normal")) else NULL
    out[rows, ] <- .Call("C_np_conditional_lp_pair_se", z, d,
      as.matrix(w), r, slope, PACKAGE = "np")
  }
  out
}
