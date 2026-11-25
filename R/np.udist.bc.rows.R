# Sum-of-products Gaussian primitive; zero-aware recurrence never divides by Phi.
.np_udist_bc_continuous <- function(x, grid, h, corrected = TRUE) {
  P <- matrix(1, nrow(grid), nrow(x)); Q <- P * 0
  for (j in seq_len(ncol(x))) {
    u <- outer(grid[[j]], x[[j]], "-") / h[j]
    phi <- stats::pnorm(u)
    moment <- u * stats::dnorm(u)
    moment[is.infinite(u)] <- 0
    Q <- Q * phi + P * moment
    P <- P * phi
  }
  list(original = P, corrected = if (corrected) P + .5 * Q else P)
}

# Fourth-order Epanechnikov primitive on its standardized support [-1, 1].
.np_udist_epanechnikov_cdf4 <- function(t, t2 = t*t) {
  .5 + t*(45 + t2*(21*t2 - 50))/32
}

.np_udist_uniform_cdf4 <- function(t) {
  .5 + t*(9 - 5*t*t)/8
}

.np_udist_bc_rows <- function(state, grid, corrected = TRUE) {
  if (!identical(state$type, "fixed")) return(.np_udist_bc_gnn_rows(state, grid, corrected))
  .np_udist_bc_grid(state, grid)
  if (!is.logical(corrected) || length(corrected) != 1L || is.na(corrected))
    stop("corrected must be TRUE or FALSE", call. = FALSE)
  if (nrow(grid) > 25L) stop("CDF row provider requires bounded target tiles", call. = FALSE)
  .np_udist_bc_budget(state, nrow(grid), 1L, nrow(grid))
  halfline <- identical(state$metadata$kernel, "gaussian2-halfline")
  if (halfline) {
    z <- .np_udist_halfline_rows(state$x[[1L]], grid[[1L]], state$h[1L], state$metadata$support)
    if (!corrected) z$corrected <- z$original
  } else if (identical(state$metadata$kernel, "beta2")) {
    component <- .np_udist_beta_rows(grid[[1L]], state$beta, state$beta$left, state$beta$right)
    z <- list(original = component$ordinary,
              corrected = if (corrected) 2*component$ordinary-component$wide else component$ordinary)
  } else if (identical(state$metadata$kernel, "uniform2")) {
    t <- outer(grid[[1L]], state$x[[1L]], "-") / state$h[1L]
    t[] <- pmax(-1, pmin(1, t))
    original <- .5 + .5*t
    z <- list(original = original, corrected = if (corrected) .np_udist_uniform_cdf4(t) else original)
  } else if (identical(state$metadata$kernel, "epanechnikov2")) {
    t <- outer(grid[[1L]], state$x[[1L]], "-") / state$h[1L] / sqrt(5)
    t[] <- pmax(-1, pmin(1, t))
    t2 <- t*t
    original <- .5 + .25*t*(3-t2)
    # The primitive of the moment-corrected kernel has signed range outside [0, 1].
    z <- list(original = original, corrected = if (corrected) .np_udist_epanechnikov_cdf4(t, t2) else original)
  } else {
    z <- .np_udist_bc_continuous(state$x[state$numeric], grid[state$numeric], state$h[state$numeric], corrected)
  }
  for (category in state$category) {
    r <- outer(as.double(toMatrix(grid[category$index])), category$training, "-")
    z$original <- z$original * .np_udist_bc_ordered(r, category$lambda, FALSE)
    z$corrected <- z$corrected * .np_udist_bc_ordered(r, category$lambda, corrected)
  }
  if (any(!is.finite(z$corrected))) stop("nonfinite CDF contributions", call. = FALSE)
  shifted <- z$corrected - z$corrected[, 1L]
  mean.shift <- rowMeans(shifted)
  z$center <- z$corrected[, 1L] + mean.shift
  divisor <- if (halfline) sqrt(as.double(state$n)*(state$n-1)) else state$n
  z$factor <- (shifted - mean.shift) / divisor
  z$original.center <- z$original[, 1L] + rowMeans(z$original - z$original[, 1L])
  z$se <- sqrt(rowSums(z$factor^2))
  z
}

.np_udist_bc_covariance <- function(state, grid, tile = 25L, corrected = TRUE) {
  m <- nrow(grid); .np_udist_bc_budget(state, m, 1L, tile)
  if (m > 32L) stop("CDF dense covariance is bounded at 32 targets", call. = FALSE)
  blocks <- .np_bias_target_blocks(m, 1L, anchors = tile)
  C <- matrix(0, m, m); center <- original <- numeric(m)
  for (a in seq_along(blocks)) {
    ii <- blocks[[a]]$grid
    left <- .np_udist_bc_rows(state, grid[ii, , drop = FALSE], corrected)
    center[ii] <- left$center; original[ii] <- left$original.center
    C[ii, ii] <- tcrossprod(left$factor)
    if (a > 1L) for (b in seq_len(a - 1L)) {
      jj <- blocks[[b]]$grid
      right <- .np_udist_bc_rows(state, grid[jj, , drop = FALSE], corrected)
      cross <- tcrossprod(left$factor, right$factor)
      C[ii, jj] <- cross; C[jj, ii] <- t(cross)
    }
  }
  list(covariance = C, center = center, original.center = original,
       se = sqrt(diag(C)))
}
