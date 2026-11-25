# Map every pair of coordinates, retaining factor levels and formula order.
.np_conditional_profile_rows <- function(ex) {
  first <- match(ex[[1L]], unique(ex[[1L]]))
  second <- match(ex[[2L]], unique(ex[[2L]]))
  key <- first + as.double(nrow(ex)) * (second - 1)
  list(index = match(key, unique(key)), first = !duplicated(key))
}

.np_conditional_profile_weights <- function(bws, x, ex) {
  h <- .npcdhat_physical_bandwidth(bws, "x")
  g <- untangle(x)
  wide <- h * ifelse(g$icon, sqrt(2), 2)
  upper <- rep(Inf, length(h))
  upper[g$iuno] <- vapply(x[g$iuno], function(v) 1 - 1 / nlevels(v), numeric(1L))
  upper[g$iord] <- 1
  if (any(!is.finite(wide)) || any(wide > upper))
    stop("widened categorical bandwidth is outside the fitted kernel's domain", call. = FALSE)
  weights <- function(width) {
    kbw <- .npcdhat_retained_kbandwidth(bw = width, ckertype = bws$cxkertype,
      ckerorder = bws$cxkerorder, bwtype = bws$type, bwscaling = FALSE,
      ukertype = bws$uxkertype, okertype = .np_density_okertype(bws$oxkertype),
      nobs = nrow(x), xdati = g, xnames = names(x))
    mat <- .np_conditional_literal_density_matrix(kbw, x, ex, "normal")
    total <- rowSums(mat)
    if (any(!is.finite(total)) || any(total <= 0) || any(mat < 0))
      stop("undefined conditioning mass on the requested panel", call. = FALSE)
    t(mat / total)
  }
  list(targets = seq_len(nrow(ex)), main = weights(h), wide = weights(wide))
}

.np_conditional_profile_inference <- function(bws, x, y, ex, ey, B, alpha, pointwise) {
  coordinates <- .np_conditional_profile_rows(ex)
  grid <- data.frame(x = coordinates$index, y = ey[[1L]])
  hx <- .npcdhat_physical_bandwidth(bws, "x")
  hy <- .npcdhat_physical_bandwidth(bws, "y")
  continuous <- which(vapply(x, is.numeric, logical(1L)))
  h <- c(hx[continuous[1L]], hy)
  training <- x[[continuous[1L]]]
  allocation <- .np_conditional_workspace(training, y[[1L]], grid, h, B, TRUE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= .5)
    stop("alpha must lie in (0, 0.5)", call. = FALSE)
  if (npLogicalOption("np.largeh", TRUE)) {
    tol <- getOption("np.largeh.rel.tol", 1e-3)
    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= .1)
      stop("invalid np.largeh.rel.tol", call. = FALSE)
    cutoff <- sqrt(2) * sqrt(-2 * log(1 - tol))
    near <- vapply(continuous, function(j) diff(range(x[[j]], ex[[j]])) <= cutoff * hx[j], logical(1L))
    if (any(near) || diff(range(y[[1L]], ey[[1L]])) <= cutoff * hy)
      stop("experimental conditional CDF inference does not admit the near-constant kernel substitution", call. = FALSE)
  }
  .np_conditional_runtime_admit(nrow(x), allocation$B)
  prepare <- function() .np_conditional_profile_weights(bws, x, ex[coordinates$first, , drop = FALSE])
  answer <- .np_conditional_pairs_score(training, y[[1L]], grid, h, allocation$B, alpha,
    original = NULL, bytes = allocation$estimated.bytes, tile.size = allocation[["tile.size"]],
    pointwise = pointwise, profile.weights = prepare)
  answer$physical.bandwidth <- c(hx, hy)
  answer$grid <- cbind(ex, ey)
  answer
}
