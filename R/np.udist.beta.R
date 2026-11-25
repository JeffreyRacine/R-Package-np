# Fixed-support beta CDF contributions. Donors determine the shape; queries
# determine the integration limit. Compute each physical tail before scaling.
.np_udist_beta_support <- function(bws, x, h) {
  lower <- bws[["ckerlb", exact = TRUE]]
  upper <- bws[["ckerub", exact = TRUE]]
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 1L ||
      length(upper) != 1L || !is.finite(upper-lower) || upper <= lower ||
      any(x < lower | x > upper))
    stop("beta CDF requires finite known bounds containing the training observations", call. = FALSE)
  if (any(x == lower | x == upper))
    stop("beta corrected CDF inference requires observations strictly inside the known support; boundary mass or rounding is not covered", call. = FALSE)
  width <- unname(upper-lower)
  precision <- (width / h)^2
  if (!is.finite(precision) || precision <= 0)
    stop("beta CDF concentration is not representable", call. = FALSE)
  list(lower = unname(lower), upper = unname(upper), width = width,
       precision = unname(precision), left = (x-lower)/width, right = (upper-x)/width)
}

.np_udist_beta_rows <- function(grid, support, left, right) {
  precision <- support$precision
  evaluate <- function(t) {
    upper <- t > support$lower + .5 * support$width
    if (upper) {
      coordinate <- (support$upper-t)/support$width
      a <- right; b <- left
    } else {
      coordinate <- (t-support$lower)/support$width
      a <- left; b <- right
    }
    ordinary <- pbeta(coordinate, 1+a*precision, 1+b*precision, lower.tail = !upper)
    wide <- pbeta(coordinate, 1+a*(precision/2), 1+b*(precision/2), lower.tail = !upper)
    cbind(ordinary, wide)
  }
  values <- array(vapply(grid, evaluate, matrix(0, length(left), 2L)),
                  c(length(left), 2L, length(grid)))
  get.scale <- function(k) t(matrix(values[, k, ], length(left), length(grid)))
  list(ordinary = get.scale(1L), wide = get.scale(2L))
}

# Each primitive is decreasing in the donor coordinate at a fixed query
# (monotone likelihood ratio). Bound their signed difference cell by cell.
.np_udist_beta_ranges <- function(grid, support) {
  edges <- seq(0, 1, length.out = 65L)
  cells <- function(t) {
    v <- .np_udist_beta_rows(t, support, edges, 1-edges)
    lower <- min(2*v$ordinary[-1L] - v$wide[-65L])
    upper <- max(2*v$ordinary[-65L] - v$wide[-1L])
    # Round outward and intersect with the exact probability enclosure.
    c(max(-1, lower-2e-13), min(2, upper+2e-13))
  }
  enclosure <- vapply(grid, cells, numeric(2L))
  list(lower = enclosure[1L, ], upper = enclosure[2L, ])
}
