# Fixed-support density operators. Never enumerate joint profiles or form a
# support-square matrix. The two-point probe retains all marginal levels.
.np_udens_bc_mixed_prepare <- function(bws, x, max.bytes, admit.effective.constant = FALSE,
                                      admit.ordered.literal = FALSE) {
  if (!is.logical(admit.effective.constant) || length(admit.effective.constant) != 1L ||
      is.na(admit.effective.constant))
    stop("invalid effective categorical operator admission", call. = FALSE)
  if (!is.logical(admit.ordered.literal) || length(admit.ordered.literal) != 1L ||
      is.na(admit.ordered.literal))
    stop("invalid literal ordered operator admission", call. = FALSE)
  numeric <- vapply(x, is.numeric, logical(1L))
  ordered <- vapply(x, is.ordered, logical(1L))
  unordered <- vapply(x, is.factor, logical(1L)) & !ordered
  if (any(!(numeric | ordered | unordered)) ||
      any(!is.finite(as.matrix(x[numeric]))))
    stop("density correction requires numeric, factor or ordered coordinates", call. = FALSE)
  if ((any(unordered) && !identical(bws[["ukertype", exact = TRUE]], "aitchisonaitken")) ||
      (any(ordered) && !identical(bws[["okertype", exact = TRUE]], "liracine")))
    stop("mixed density correction requires default AA and ordered LR kernels", call. = FALSE)
  if (!identical(unname(bws[["icon", exact = TRUE]]), unname(numeric)) ||
      !identical(unname(bws[["iuno", exact = TRUE]]), unname(unordered)) ||
      !identical(unname(bws[["iord", exact = TRUE]]), unname(ordered)))
    stop("density correction type metadata disagrees with training frame", call. = FALSE)
  h <- .np_physical_bandwidth(bws)
  if (length(h) != ncol(x) || any(!is.finite(h)) || any(h[numeric] <= 0))
    stop("invalid mixed density physical bandwidths", call. = FALSE)
  .np_udens_bc_budget(list(x = x, n = nrow(x), max.bytes = max.bytes), 1L, 1L, 1L, FALSE)
  con <- NULL
  if (any(numeric)) {
    xx <- x[numeric]
    kb <- kbandwidth.numeric(bw = h[numeric], bwscaling = FALSE, bwtype = "fixed",
      ckertype = "gaussian", ckerorder = 2L, ckerbound = "none",
      nobs = nrow(x), xdati = untangle(xx), xnames = names(xx))
    con <- list(x = xx, kbw = kb, n = nrow(x), q = ncol(xx), r = 0L, s = 0L,
                h = h[numeric], max.bytes = max.bytes)
  }
  category <- vector("list", sum(!numeric)); indices <- which(!numeric)
  for (a in seq_along(indices)) {
    j <- unname(indices[a]); k <- nlevels(x[[j]]); lambda <- h[j]
    if (k < 1L || lambda < 0 || (unordered[j] && lambda > (k-1)/k) ||
        (ordered[j] && lambda >= 1))
      stop("categorical density bandwidth is outside the admitted range", call. = FALSE)
    support <- dlev(x[[j]])
    if (ordered[j]) .np_ordered_levels_contract(support, "nliracine")
    if (ordered[j] && any(!is.finite(support) | support != round(support)))
      stop("ordered density correction requires declared integer-lattice support", call. = FALSE)
    probe <- data.frame(z = factor(levels(x[[j]])[c(1L, k)],
                                  levels = levels(x[[j]]), ordered = ordered[j]))
    kb <- kbandwidth.numeric(bw = lambda, bwscaling = FALSE, bwtype = "fixed",
      ckertype = "gaussian", ckerorder = 2L, ckerbound = "none",
      ukertype = "aitchisonaitken", okertype = "nliracine",
      nobs = 2L, xdati = untangle(probe), xnames = names(probe))
    codes <- as.double(toMatrix(probe))
    distance <- abs(outer(codes, codes, "-"))
    raw <- if (unordered[j]) {
      if (k == 1L) matrix(1, 2L, 2L) else
        ifelse(distance == 0, 1-lambda, lambda/(k-1))
    } else (1-lambda)/(1+lambda) * lambda^distance
    effective <- .np_plot_density_factor_weights(raw, kb, probe, probe, "normal")
    pooling <- unordered[j] && k > 1L && lambda == (k-1)/k
    constant <- if (pooling) effective[1L, 1L] else NULL
    owner <- if (pooling) "exact AA pooling" else "raw normalized density kernel"
    certificate <- lattice <- NULL
    if (!identical(raw, effective) && !pooling) {
      if (admit.effective.constant && unordered[j])
        certificate <- .np_udens_bc_effective_aa_constant(lambda, k, effective)
      if (admit.ordered.literal && ordered[j])
        lattice <- .np_udens_bc_literal_lr_certificate(lambda, support, raw, effective)
      if (is.null(certificate) && is.null(lattice))
        stop("near-upper effective categorical operator is outside density correction admission", call. = FALSE)
      if (!is.null(certificate)) {
        constant <- certificate$constant
        owner <- "native effective AA constant"
      } else owner <- "literal normalized LR lattice"
    }
    category[[a]] <- list(index = j, ordered = unname(ordered[j]), levels = levels(x[[j]]),
      support = support, k = k, lambda = lambda,
      training = as.double(toMatrix(x[j])), pooling = pooling,
      constant = constant, owner = owner,
      raw.probe = raw, effective.probe = effective)
    if (!is.null(certificate)) category[[a]]$effective.constant <- certificate
    if (!is.null(lattice)) category[[a]]$literal.lattice <- lattice
  }
  list(x = x, n = nrow(x), q = sum(numeric), r = sum(unordered), s = sum(ordered),
    h = h, max.bytes = max.bytes, continuous = con, numeric = numeric,
    category = category, construction = "density-bc-v1",
    owner.options = list(largelambda = npLogicalOption("np.largelambda", TRUE),
                         tolerance = npDiscUpperRelTol()),
    metadata = list(type = "fixed", continuous = "gaussian2", ratio = 1,
      support = "Lebesgue x declared AA support x integer lattice",
      normalization = "n-squared", categorical.operator = "2S-S-squared",
      categorical = lapply(category, function(z) z[c("index", "ordered", "levels", "support",
        "k", "lambda", "pooling", "constant", "owner", "raw.probe", "effective.probe")])))
}

.np_udens_bc_mixed_grid <- function(state, grid) {
  if (!is.data.frame(grid) || nrow(grid) < 1L ||
      !identical(names(grid), names(state$x)) || anyNA(grid))
    stop("invalid mixed density correction evaluation frame", call. = FALSE)
  if (!identical(state$owner.options, list(largelambda = npLogicalOption("np.largelambda", TRUE),
                                           tolerance = npDiscUpperRelTol())))
    stop("categorical owner options changed within the fixed density fit", call. = FALSE)
  if (!is.null(state$continuous))
    .np_udens_bc_grid(state$continuous, grid[state$numeric])
  for (z in state$category) {
    value <- grid[[z$index]]
    if (!is.factor(value) || !identical(is.ordered(value), z$ordered) ||
        !identical(levels(value), z$levels))
      stop("density evaluation categories must retain the fitted support and type", call. = FALSE)
  }
  invisible(TRUE)
}

.np_udens_bc_categorical_rows <- function(category, evaluation) {
  distance <- abs(outer(as.double(toMatrix(data.frame(evaluation))), category$training, "-"))
  lambda <- category$lambda; k <- category$k
  if (category$ordered) {
    c <- (1-lambda)/(1+lambda)
    S <- c * lambda^distance
    Q <- S * (2 - c*(distance+1) - 2*lambda^2/(1+lambda)^2)
  } else if (k == 1L) {
    S <- Q <- matrix(1, nrow(distance), ncol(distance))
  } else if (!is.null(category$constant)) {
    a <- category$constant
    S <- matrix(a, nrow(distance), ncol(distance))
    Q <- matrix(2*a-k*a^2, nrow(distance), ncol(distance))
  } else {
    effective <- k*lambda^2/(k-1)
    S <- ifelse(distance == 0, 1-lambda, lambda/(k-1))
    Q <- ifelse(distance == 0, 1-effective, effective/(k-1))
  }
  list(original = S, corrected = Q)
}

.np_udens_bc_mixed_contribution <- function(state, grid, corrected) {
  .np_udens_bc_mixed_grid(state, grid)
  .np_udens_bc_budget(state, nrow(grid), 1L, nrow(grid), FALSE)
  if (is.null(state$continuous)) {
    K <- L <- matrix(1, nrow(grid), state$n)
  } else {
    continuous <- .np_udens_bc_contribution(state$continuous, grid[state$numeric], corrected)
    K <- continuous$original; L <- continuous$corrected
  }
  for (category in state$category) {
    value <- .np_udens_bc_categorical_rows(category, grid[[category$index]])
    K <- K * value$original
    L <- L * if (corrected) value$corrected else value$original
  }
  if (any(!is.finite(K)) || any(!is.finite(L)))
    stop("nonfinite mixed density correction contributions", call. = FALSE)
  center <- L[, 1L] + rowMeans(L - L[, 1L])
  F <- (L-center)/state$n
  list(original = K, corrected = L, original.center = rowMeans(K),
       center = center, factor = F, se = sqrt(rowSums(F^2)))
}

# Opt-in certificate for finite AA support only. Query the unchanged native
# owner; never infer a constant operator by rounding a nonconstant witness.
.np_udens_bc_effective_aa_constant <- function(lambda, k, effective) {
  if (k < 2L) return(NULL)
  probe <- matrix(c(1, 2), 2L, 1L)
  empty <- matrix(numeric(), 2L, 0L)
  native <- .np_categorical_profile_kernel_tile(probe, empty, probe, empty,
    0L, integer(), 0L, lambda, k, list(NULL))
  c <- native[1L, 1L]
  if (any(!is.finite(native)) || c <= 0 || !all(native == c)) return(NULL)
  # This check concerns only the existing R multiply/divide round trip.
  # Constancy itself has already been established exactly by the native tile.
  tol <- 64 * .Machine$double.eps * max(1, abs(c))
  if (any(!is.finite(effective)) || any(abs(effective-c) > tol)) return(NULL)
  corrected <- 2*c-k*c^2
  if (!is.finite(corrected))
    stop("effective categorical correction arithmetic is undefined", call. = FALSE)
  list(constant = c, corrected.constant = corrected,
       effective.mass = k*c, corrected.mass = k*corrected,
       sampling = "unchanged normalized raw AA law", native.probe = native)
}

# Private corrected curves use the normalized infinite-lattice operator.  The
# native product fit's finite-geometry constant is retained only as a diagnostic;
# it cannot be convolved as a constant over all integers.  No native option changes.
.np_udens_bc_literal_lr_certificate <- function(lambda, support, raw, effective) {
  probe <- matrix(support[c(1L, length(support))], 2L, 1L)
  empty <- matrix(numeric(), 2L, 0L)
  native <- .np_categorical_profile_kernel_tile(empty, probe, empty, probe,
    integer(), 2L, 0L, lambda, length(support), list(as.double(support)))
  tolerance <- 64 * .Machine$double.eps * pmax(abs(raw), .Machine$double.xmin)
  if (any(!is.finite(raw)) || any(!is.finite(native)) ||
      any(abs(native-raw) > tolerance)) return(NULL)
  list(kernel = "normalized LR", measure = "integer lattice", mass = 1,
       correction = "2K-K*K", native.scalar.probe = native,
       fit.factor.probe = effective, fit.probe.finite = all(is.finite(effective)),
       sampling = "unchanged difference of geometricals")
}
