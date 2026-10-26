# One bounded contribution tile; no response, residual, pilot or random draw.
.np_udens_bc_contribution <- function(state, grid, corrected = TRUE) {
  .np_udens_bc_grid(state, grid)
  if (!is.logical(corrected) || length(corrected) != 1L || is.na(corrected))
    stop("corrected must be TRUE or FALSE", call. = FALSE)
  if (nrow(grid) > 25L)
    stop("density contribution requests must use at most 25 target rows", call. = FALSE)
  .np_udens_bc_budget(state, nrow(grid), 1L, nrow(grid), FALSE)
  K <- .np_plot_density_matrix(state$kbw, state$x, grid, rep("normal", state$q))
  radius <- matrix(0, nrow(grid), state$n)
  for (j in seq_len(state$q))
    radius <- radius + (outer(grid[[j]], state$x[[j]], "-") / state$h[j])^2
  L <- if (corrected) K * (1 + state$q/2 - radius/2) else K
  if (any(!is.finite(K)) || any(!is.finite(L)))
    stop("nonfinite density correction contributions", call. = FALSE)
  # Anchor before averaging: identical observations certify exact zero
  # uncertainty, rather than a small positive scale from summation roundoff.
  center <- L[, 1L] + rowMeans(L - L[, 1L])
  F <- (L - center) / state$n
  list(original = K, corrected = L, original.center = rowMeans(K),
       center = center, factor = F, se = sqrt(rowSums(F^2)))
}

# Complete cross-tile covariance, bounded by the existing fixed grid cap.
.np_udens_bc_covariance <- function(state, grid, tile = 25L, corrected = TRUE) {
  .np_udens_bc_grid(state, grid)
  m <- nrow(grid)
  .np_udens_bc_budget(state, m, 1L, tile, TRUE)
  blocks <- .np_bias_target_blocks(m, 1L, anchors = tile)
  C <- matrix(0, m, m); center <- original <- se <- numeric(m)
  for (a in seq_along(blocks)) {
    ii <- blocks[[a]]$grid
    left <- .np_udens_bc_contribution(state, grid[ii, , drop = FALSE], corrected)
    center[ii] <- left$center; original[ii] <- left$original.center; se[ii] <- left$se
    C[ii, ii] <- tcrossprod(left$factor)
    if (a > 1L) for (b in seq_len(a - 1L)) {
      jj <- blocks[[b]]$grid
      right <- .np_udens_bc_contribution(state, grid[jj, , drop = FALSE], corrected)
      cross <- tcrossprod(left$factor, right$factor)
      C[ii, jj] <- cross; C[jj, ii] <- t(cross)
    }
  }
  if (any(!is.finite(se)) || any(se <= 0))
    stop("density inference requires positive finite observed scales", call. = FALSE)
  list(covariance = C, center = center, original.center = original, se = se)
}
