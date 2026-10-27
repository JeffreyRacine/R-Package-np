# Internal fixed-model state for the opt-in density construction.
.np_udens_bc_prepare <- function(bws, x, max.bytes = 256 * 1024^2) {
  if (!inherits(bws, "bandwidth") ||
      !identical(bws[["type", exact = TRUE]], "fixed") ||
      !identical(bws[["ckertype", exact = TRUE]], "gaussian") ||
      !isTRUE(bws[["ckerorder", exact = TRUE]] == 2L) ||
      !identical(bws[["ckerbound", exact = TRUE]], "none"))
    stop("density correction requires fixed unbounded second-order Gaussian bandwidths", call. = FALSE)
  if (!is.data.frame(x) || nrow(x) < 2L || ncol(x) < 1L ||
      anyNA(x) || anyDuplicated(names(x)) || any(!nzchar(names(x))))
    stop("invalid density correction training frame", call. = FALSE)
  if (!is.numeric(max.bytes) || length(max.bytes) != 1L ||
      !is.finite(max.bytes) || max.bytes <= 0)
    stop("density correction memory budget must be positive and finite", call. = FALSE)
  if (!all(vapply(x, is.numeric, logical(1L))))
    return(.np_udens_bc_mixed_prepare(bws, x, max.bytes))
  if (any(!is.finite(as.matrix(x))))
    stop("density correction requires finite numeric coordinates", call. = FALSE)
  h <- .np_physical_bandwidth(bws)
  if (!is.numeric(h) || length(h) != ncol(x) || any(!is.finite(h)) || any(h <= 0))
    stop("density correction requires positive finite physical bandwidths", call. = FALSE)
  kbw <- .np_make_kbandwidth_unconditional(bws, x)
  if (!all(kbw[["icon", exact = TRUE]]))
    stop("density correction type metadata disagrees with training frame", call. = FALSE)
  list(x = x, kbw = kbw, n = nrow(x), q = ncol(x), r = 0L, s = 0L,
       h = h, max.bytes = max.bytes, construction = "density-bc-v1",
       metadata = list(type = "fixed", continuous = "gaussian2", ratio = 1,
                       support = "Lebesgue", normalization = "n-squared"))
}

.np_udens_bc_grid <- function(state, grid) {
  if (state$r + state$s > 0L) return(.np_udens_bc_mixed_grid(state, grid))
  if (!is.data.frame(grid) || nrow(grid) < 1L ||
      !identical(names(grid), names(state$x)) || anyNA(grid) ||
      !all(vapply(grid, is.numeric, logical(1L))) || any(!is.finite(as.matrix(grid))))
    stop("invalid density correction evaluation frame", call. = FALSE)
  # The incumbent fit may substitute a constant Gaussian factor near infinite
  # smoothing. Raw matrix helpers intentionally do not. Exclude that new
  # inference class rather than changing the incumbent option or its owner.
  if (npLogicalOption("np.largeh", TRUE)) {
    tolerance <- getOption("np.largeh.rel.tol", 1e-3)
    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
        !is.finite(tolerance) || tolerance <= 0 || tolerance >= .1)
      stop("invalid np.largeh.rel.tol", call. = FALSE)
    threshold <- sqrt(-2 * log(1 - tolerance))
    span <- vapply(seq_len(state$q), function(j)
      diff(range(state$x[[j]], grid[[j]])), numeric(1L))
    if (any(span <= state$h * threshold))
      stop("near-constant Gaussian owner is outside this density correction admission", call. = FALSE)
  }
  invisible(TRUE)
}

.np_udens_bc_budget <- function(state, m, B, tile, covariance, processes = 1L) {
  tile <- .np_plot_gaussian_integer(tile, "density target tile", 25L)
  B <- .np_plot_gaussian_integer(B, "density draws", .Machine$integer.max)
  if (covariance && m > 512L)
    stop("density covariance representation exceeds the 512-point cap", call. = FALSE)
  n <- as.double(state$n); p <- ncol(state$x); t <- min(tile, m)
  k <- if (covariance) m else 0
  processes <- .np_plot_gaussian_integer(processes, "density resident processes", .Machine$integer.max)
  levels <- sum(vapply(state$x, function(z) length(levels(z)), integer(1L)))
  support.bytes <- sum(vapply(state$x, function(z)
    if (is.factor(z)) as.numeric(utils::object.size(levels(z))) else 0, numeric(1L)))
  per.process <- 8 * (4*n*p + 12*n*t + 4*n*32 + 4*t*B + 8*B + 12*k^2 + 4*m*p + 4*levels) + 4*support.bytes
  # Conservatively charge the full coordinator forecast to every resident
  # process, including replicated inputs, serialized tasks and R copies.
  bytes <- processes * per.process
  if (!is.finite(bytes) || bytes > state$max.bytes)
    stop("density correction live-byte forecast exceeds the requested budget", call. = FALSE)
  list(bytes = bytes, per.process = per.process, processes = processes,
       retained.levels = levels, tile = tile, draw.chunk = 32L, covariance.dimension = k)
}
