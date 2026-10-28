# Private CDF prototype: immutable retained-data context, no public dispatch.
.np_udist_bc_prepare <- function(bws, x, max.bytes = 256 * 1024^2) {
  if (is.list(bws) && bws[["type", exact = TRUE]] %in% c("generalized_nn", "adaptive_nn"))
    return(.np_udist_bc_nn_prepare(bws, x, max.bytes))
  if (!is.list(bws) || !identical(bws[["type", exact = TRUE]], "fixed") ||
      !identical(bws[["ckertype", exact = TRUE]], "gaussian") ||
      !isTRUE(bws[["ckerorder", exact = TRUE]] == 2L) ||
      !identical(bws[["ckerbound", exact = TRUE]], "none"))
    stop("CDF correction requires fixed unbounded second-order Gaussian bandwidths", call. = FALSE)
  if (!is.data.frame(x) || nrow(x) < 2L || !ncol(x) || anyNA(x) ||
      anyDuplicated(names(x)) || any(!nzchar(names(x))) ||
      !identical(names(x), bws[["xnames", exact = TRUE]]) ||
      !isTRUE(bws[["nobs", exact = TRUE]] == nrow(x)))
    stop("invalid CDF correction training frame", call. = FALSE)
  numeric <- vapply(x, is.numeric, logical(1L))
  ordered <- vapply(x, is.ordered, logical(1L))
  if (any(!(numeric | ordered)) || any(!is.finite(as.matrix(x[numeric]))) ||
      !identical(unname(bws[["icon", exact = TRUE]]), unname(numeric)) ||
      !identical(unname(bws[["iord", exact = TRUE]]), unname(ordered)))
    stop("CDF correction accepts numeric and ordered variables only", call. = FALSE)
  if (any(ordered) && !identical(bws[["okertype", exact = TRUE]], "liracine"))
    stop("CDF correction requires the normalized ordered LR kernel", call. = FALSE)
  h <- .np_physical_bandwidth(bws)
  if (length(h) != ncol(x) || any(!is.finite(h)) || any(h[numeric] <= 0) ||
      any(h[ordered] < 0 | h[ordered] >= 1))
    stop("invalid CDF physical bandwidths", call. = FALSE)
  if (length(max.bytes) != 1L || !is.numeric(max.bytes) ||
      !is.finite(max.bytes) || max.bytes <= 0)
    stop("invalid CDF memory budget", call. = FALSE)
  category <- lapply(which(ordered), function(j) {
    support <- dlev(x[[j]])
    .np_ordered_levels_contract(support, "nliracine")
    if (any(!is.finite(support) | support != round(support)))
      stop("ordered CDF correction requires integer lattice distances", call. = FALSE)
    probe <- data.frame(z = ordered(levels(x[[j]])[c(1L, nlevels(x[[j]]))], levels = levels(x[[j]])))
    kb <- kbandwidth.numeric(bw = h[j], bwscaling = FALSE, bwtype = "fixed",
      okertype = "nliracine", nobs = 2L, xdati = untangle(probe), xnames = "z")
    r <- outer(as.double(toMatrix(probe)), as.double(toMatrix(probe)), "-")
    raw <- .np_udist_bc_ordered(r, h[j], FALSE)
    effective <- .np_plot_density_factor_weights(raw, kb, probe, probe, "integral")
    if (!identical(raw, effective)) stop("effective ordered CDF owner is outside correction admission", call. = FALSE)
    list(index = j, lambda = h[j], levels = levels(x[[j]]), support = support,
      training = as.double(toMatrix(x[j])), raw.probe = raw, effective.probe = effective)
  })
  list(x = x, n = nrow(x), q = sum(numeric), s = sum(ordered), numeric = numeric,
    h = h, bws = bws, type = "fixed", max.bytes = max.bytes,
    row.ids = rownames(x), category = category, construction = "cdf-bc-v1",
    metadata = list(normalization = "n-squared", support = "Lebesgue x integer lattice",
      coefficient = .5, ordered.operator = "R(2S-S*S)", kernel = "gaussian2"))
}

.np_udist_bc_grid <- function(state, grid) {
  if (!is.data.frame(grid) || !nrow(grid) || anyNA(grid) ||
      !identical(names(grid), names(state$x)) ||
      !all(vapply(grid[state$numeric], is.numeric, logical(1L))) ||
      any(!is.finite(as.matrix(grid[state$numeric]))))
    stop("invalid CDF evaluation frame", call. = FALSE)
  for (z in state$category) if (!is.ordered(grid[[z$index]]) ||
      !identical(levels(grid[[z$index]]), z$levels))
    stop("CDF evaluation must retain fitted ordered levels", call. = FALSE)
  if (identical(state$type, "fixed") && state$q && npLogicalOption("np.largeh", TRUE)) {
    threshold <- sqrt(-2*log(1-getOption("np.largeh.rel.tol", 1e-3)))
    span <- vapply(which(state$numeric), function(j) diff(range(state$x[[j]], grid[[j]])), numeric(1L))
    if (any(span <= state$h[state$numeric]*threshold))
      stop("large-h/global-limit CDF correction is outside this investigation", call. = FALSE)
  }
  invisible(TRUE)
}

# Cumulative normalized lattice kernel and its full convolution correction.
.np_udist_bc_ordered <- function(r, lambda, corrected = TRUE) {
  if (lambda == 0) return(1 * (r >= 0))
  tail <- function(k) if (corrected)
    lambda^k * (1 + 3*lambda - k*(1-lambda^2))/(1+lambda)^3 else
    lambda^k/(1+lambda)
  ifelse(r < 0, tail(abs(r)), 1 - tail(pmax(1, r+1)))
}

.np_udist_bc_budget <- function(state, m, B, tile, processes = max(1L, mpi.comm.size(1L))) {
  tile <- .np_plot_gaussian_integer(tile, "CDF target tile", 25L)
  t <- min(tile, m); n <- as.double(state$n); p <- ncol(state$x)
  per.process <- 8 * (8*n*p + 20*n*t + 8*n*32 + 6*t*B + 12*B + 16*min(m, 32L)^2 + 12*max(0, m-32L))
  bytes <- processes * per.process
  if (!is.finite(bytes) || bytes > state$max.bytes)
    stop("CDF resident memory forecast exceeds budget", call. = FALSE)
  list(bytes = bytes, per.process = per.process, processes = processes,
       target.tile = tile, draw.chunk = 32L, covariance.cap = 32L)
}
