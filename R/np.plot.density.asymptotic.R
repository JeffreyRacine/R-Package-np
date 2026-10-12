# Joint variability for density-family plots only. The incumbent fit and marginal SEs
# are unchanged; these empirical influence rows supply cross-grid dependence.  Kernels,
# basis preparation and accepted solves retain their canonical owners.
.np_plot_density_joint_supported <- function(bws, conditional) {
  second <- function(kernel, order) kernel %in% c("gaussian", "epanechnikov") && identical(as.integer(order),
    2L)
  if (!conditional) {
    return(bws[["ncon", exact = TRUE]] == 0L || second(bws[["ckertype", exact = TRUE]],
      bws[["ckerorder", exact = TRUE]]) || (identical(bws[["type", exact = TRUE]],
      "adaptive_nn") && identical(bws[["ckertype", exact = TRUE]], "uniform")))
  }
  xok <- bws[["xncon", exact = TRUE]] == 0L || second(bws[["cxkertype", exact = TRUE]],
    bws[["cxkerorder", exact = TRUE]])
  yok <- bws[["yncon", exact = TRUE]] == 0L || second(bws[["cykertype", exact = TRUE]],
    bws[["cykerorder", exact = TRUE]])
  if (xok && yok)
    return(TRUE)
  spec <- npConditionalRegEngineSpec(bws, where = "Gaussian density plot")
  lc <- identical(spec$reg.engine, "lc") || all(spec$degree.engine == 0L)
  # Bounded beta derivatives have separate normalized kernel semantics.  Admission is
  # restricted to the qualified continuous fixed LC contract.
  identical(bws[["type", exact = TRUE]], "fixed") && lc && bws[["xncon", exact = TRUE]] ==
    bws[["xndim", exact = TRUE]] && bws[["yncon", exact = TRUE]] == bws[["yndim", exact = TRUE]] &&
    (xok || identical(bws[["cxkertype", exact = TRUE]], "beta")) && (yok || identical(bws[["cykertype",
    exact = TRUE]], "beta"))
}

.np_plot_density_matrix <- function(kbw, x, ex, operator) {
  .npcdhat_make_kernel_matrix(kbw, x, ex, operator = operator)
}
# Query/donor NN normalization differs. Recover GNN physical scaling using a positive
# native contraction, never a possibly cancelling derivative sum.
.np_plot_density_derivative_matrix <- function(kbw, x, ex, index) {
  op <- rep("normal", ncol(x))
  op[index] <- "derivative"
  n <- nrow(x)
  m <- nrow(ex)
  raw <- .np_local_operator_ksum(kbw, x, ex, op, matrix(1, n, 1L), TRUE, "conditional Gaussian derivative")$kernel.weights
  if (identical(kbw[["type", exact = TRUE]], "adaptive_nn"))
    return(t(raw))
  if (identical(kbw[["type", exact = TRUE]], "fixed") && !identical(kbw[["ckertype", exact = TRUE]],
    "beta")) {
    h <- kbw[["bw", exact = TRUE]][kbw[["icon", exact = TRUE]]]
    cop <- op[kbw[["icon", exact = TRUE]]]
    return(t(raw)/(prod(h[cop == "normal"]) * prod(h[cop == "derivative"]^2)))
  }
  D <- matrix(0, m, n)
  for (first in seq.int(1L, m, by = 16L)) {
    ix <- seq.int(first, min(m, first + 15L))
    w <- raw[, ix, drop = FALSE]
    mass <- colSums(abs(w))
    keep <- which(mass > 0)
    if (length(keep)) {
      q <- ix[keep]
      z <- .np_local_operator_ksum(kbw, x, ex[q, , drop = FALSE], op, sign(w[, keep,
        drop = FALSE]), FALSE, "conditional Gaussian derivative scale")$ksum
      D[q, ] <- t(w[, keep, drop = FALSE]) * (diag(z)/mass[keep])
    }
  }
  D
}
.np_plot_density_center <- function(A) A - rowMeans(A)
.np_plot_density_rescale <- function(A, se) {
  if (!is.matrix(A) || ncol(A) < 2L || length(se) != nrow(A) || any(!is.finite(A)) || any(!is.finite(se)) ||
    any(se < 0))
    stop("Invalid influence or unavailable marginal uncertainty")
  # Row-wise normalization bounds temporary memory even on a surface.
  norm <- numeric(nrow(A))
  F <- A
  for (i in seq_len(nrow(A))) {
    scale <- max(abs(A[i, ]))
    safe <- if (scale > 0)
      scale else 1
    z <- A[i, ]/safe
    sz <- sqrt(sum(z * z))
    norm[i] <- sz * scale
    if (sz == 0 && se[i] > 0)
      stop("Unidentified correlation with positive marginal uncertainty")
    F[i, ] <- if (sz > 0)
      z * (se[i]/sz) else 0
  }
  list(factor = F, se = se, empirical.se = norm)
}
# Obtain the accepted inverse and its response-intercept correction through the
# existing canonical ranked native solve. No new ridge/rank policy is selected.
.np_plot_density_lp_system <- function(W, w, n) {
  p <- ncol(W)
  M <- crossprod(W, w * W)
  packed <- unlist(lapply(seq_len(p), function(a) M[a, a:p]))
  identity <- .np_inid_lp_design_identity(W)
  rank <- .np_inid_lp_rank_bounds(w, p, identity = identity)
  Z <- lapply(seq_len(p), function(a) matrix(diag(p)[a, ], 1L, p))
  V <- matrix(NA_real_, p, p)
  ridge <- steps <- NULL
  for (a in seq_len(p)) {
    out <- .np_inid_lp_batch_project(matrix(packed, 1L), Z, diag(p)[a, ], n, diagnostics = TRUE,
      rank.bounds = rank)
    V[a, ] <- out$values[1L, ]
    if (a == 1L) {
      ridge <- out$ridge_total[1L]
      steps <- out$ridge_steps[1L]
    } else stopifnot(identical(out$ridge_total[1L], ridge), identical(out$ridge_steps[1L],
      steps))
  }
  anchor <- M[1L, 1L]
  c0 <- if (ridge == 0)
    1 else 1 + ridge/anchor
  inverse <- V
  inverse[, 1L] <- inverse[, 1L]/c0
  dr <- numeric(n)
  if (ridge > 0) {
    d <- abs(diag(M))
    mask <- diag(p) == 1
    mx <- max(d)
    if (mx == 0) {
      d <- abs(M)
      mask <- matrix(TRUE, p, p)
      mx <- max(d)
    }
    ix <- which(abs(M) == mx & mask, arr.ind = TRUE)
    if (!nrow(ix))
      stop("Unidentified ridge scale")
    da <- sign(M[ix[1, 1], ix[1, 2]]) * w * W[, ix[1, 1]] * W[, ix[1, 2]]
    for (k in seq_len(nrow(ix))) if (any(da != sign(M[ix[k, 1], ix[k, 2]]) * w * W[,
      ix[k, 1]] * W[, ix[k, 2]]))
      stop("Nonsmooth accepted ridge scale; influence not defined")
    dr <- steps/n * da
  }
  list(M = M, response.inverse = V, inverse = inverse, ridge = ridge, steps = steps, dr = dr,
    anchor = anchor)
}
.np_plot_density_lp_score <- function(sys, W, w, r, e) {
  b <- drop(crossprod(W, w * r))
  theta <- drop(sys$response.inverse %*% b)
  v <- drop(t(sys$inverse) %*% e)
  vz <- drop(W %*% v)
  zt <- drop(W %*% theta)
  a <- w * vz * (r - zt)
  if (sys$ridge > 0) {
    rho <- sys$ridge/sys$anchor
    da <- w * W[, 1L]^2
    a <- a + w * v[1L] * rho * W[, 1L] * r + v[1L] * b[1L] * (sys$dr/sys$anchor - rho *
      da/sys$anchor) - sys$dr * sum(v * theta)
  }
  list(estimate = sum(e * theta), score = a, theta = theta, ridge = sys$ridge)
}
.np_plot_density_conditional <- function(bws, x, y, ex, ey, cdf = FALSE, gradient = NULL,
  order = 1L) {
  n <- nrow(x)
  m <- nrow(ex)
  if (nrow(y) != n || nrow(ey) != m || n < 2L || !m)
    stop("Unaligned conditional data")
  if (!is.null(gradient)) {
    if (length(gradient) != 1L || !is.finite(gradient) || gradient != as.integer(gradient) ||
      gradient < 1L || gradient > ncol(x) || length(order) != 1L || order < 1L || order !=
      as.integer(order))
      stop("Invalid gradient selection")
    if (!bws[["ixcon", exact = TRUE]][gradient]) {
      if (order != 1L)
        stop("Categorical gradients require first differences")
      frames <- npCategoricalFirstDifferenceFrames(ex, index = gradient, where = "Gaussian density plot")
      u <- .np_plot_density_conditional(bws, x, y, frames$upper, ey, cdf)
      l <- .np_plot_density_conditional(bws, x, y, frames$lower, ey, cdf)
      return(list(estimate = u$estimate - l$estimate, A = u$A - l$A, target = "paired endpoint empirical-mass contrast; realized geometry",
        ridge = c(u$ridge, l$ridge)))
    }
  }
  xb <- .npcdhat_make_xkbw(bws, x)
  yb <- .npcdhat_make_ybw(bws, y)
  K <- .np_plot_density_matrix(xb, x, ex, rep("normal", ncol(x)))
  Y <- .np_plot_density_matrix(yb, y, ey, rep(if (cdf) "integral" else "normal", ncol(y)))
  spec <- npConditionalRegEngineSpec(bws, where = "Gaussian density plot")
  lc <- identical(spec$reg.engine, "lc") || all(spec$degree.engine == 0)
  if (lc) {
    den <- rowSums(K)
    mu <- rowSums(K * Y)/den
    A <- K * (Y - mu)/den
    if (!is.null(gradient)) {
      if (order != 1L)
        stop("LC supports first continuous derivatives only")
      D <- .np_plot_density_derivative_matrix(xb, x, ex, gradient)
      H <- D/den - (K/den) * (rowSums(D)/den)
      deriv <- rowSums(H * Y)
      A <- H * (Y - mu) - (K/den) * deriv
      mu <- deriv
    }
    return(list(estimate = mu, A = .np_plot_density_center(A), target = "conditional ratio empirical mass; realized geometry",
      ridge = 0))
  }
  if (!is.null(gradient) && identical(bws[["type", exact = TRUE]], "generalized_nn") &&
    order > 1L)
    stop("Incumbent generalized-NN higher derivative operator is unavailable")
  cont <- which(bws[["ixcon", exact = TRUE]])
  grad <- NULL
  if (!is.null(gradient) && order > spec$degree.engine[match(gradient, cont)])
    stop("Requested derivative exceeds the selected LP degree")
  if (!is.null(gradient)) {
    grad <- integer(length(cont))
    grad[match(gradient, cont)] <- order
  }
  basisargs <- list(xdat = x[, cont, drop = FALSE], degree = spec$degree.engine, basis = spec$basis.engine,
    bernstein.basis = spec$bernstein.engine)
  W <- as.matrix(do.call(W.lp, basisargs))
  E <- as.matrix(do.call(W.lp, c(basisargs, list(exdat = ex[, cont, drop = FALSE], gradient.vec = grad))))
  groups <- .np_plot_exact_row_groups(ex)$groups
  A <- matrix(NA_real_, m, n)
  mu <- ridge <- numeric(m)
  for (ix in groups) {
    w <- K[ix[1L], ]
    sys <- .np_plot_density_lp_system(W, w, n)
    for (q in ix) {
      v <- .np_plot_density_lp_score(sys, W, w, Y[q, ], E[q, ])
      mu[q] <- v$estimate
      A[q, ] <- v$score
      ridge[q] <- v$ridge
    }
  }
  list(estimate = mu, A = .np_plot_density_center(A), target = "LP empirical mass including smooth accepted ridge; realized geometry",
    ridge = ridge)
}
# Same-sample random-radius term. Geometry, kernels and boundary-face pilots come from
# their current owners; vectorized interval accumulation avoids n^2.
.np_plot_density_ann <- function(bws, x, ex, pdf) {
  n <- nrow(x)
  m <- nrow(ex)
  reason <- .np_ann_se_domain(bws, n, pdf)
  if (!is.null(reason))
    .np_plot_joint_unavailable(reason)
  st <- .np_ann_se_geometry(bws, x)
  if (!st$valid)
    .np_plot_joint_unavailable("nonregular ANN spacing geometry")
  kb <- .np_make_kbandwidth_unconditional(bws, x)
  p <- ncol(st$x)
  e <- as.matrix(ex[, bws[["icon", exact = TRUE]], drop = FALSE])
  # The ANN native kernel owner returns physically normalized weights. Use its
  # object-fed local entry, including for zero derivative row sums; avoid
  # reconstructing normalization from a sum that may legitimately be zero.
  local.matrix <- function(op) .np_local_operator_ksum(kb, x, ex, op, matrix(1, n, 1L),
    return.kernel.weights = TRUE, where = "ANN plot influence")$kernel.weights
  op <- rep(if (pdf) "normal" else "integral", ncol(x))
  K <- local.matrix(op)
  D <- array(NA_real_, c(n, m, p))
  for (j in seq_len(p)) {
    dop <- op
    dop[which(bws[["icon", exact = TRUE]])[j]] <- if (pdf)
      "derivative" else "normal"
    D[, , j] <- local.matrix(dop)
  }
  A <- matrix(NA_real_, m, n)
  faces <- NULL
  if (pdf && identical(bws[["ckertype", exact = TRUE]], "uniform")) {
    st$uniform <- .np_ann_uniform_prepare(kb, x)
    # Reuse the unchanged support-face formula with a local kernel dependency.  Only
    # this copied closure is adapted; the package namespace stays intact.
    face.owner <- .np_ann_uniform_faces
    deps <- new.env(parent = environment(face.owner))
    deps$npksum <- function(bws, txdat, exdat, operator, bandwidth.divide, return.kernel.weights,
      .np.internal.bandwidth.divide.weights) {
      stopifnot(identical(bws[["type", exact = TRUE]], "adaptive_nn"), isTRUE(bandwidth.divide),
        isTRUE(return.kernel.weights), isTRUE(.np.internal.bandwidth.divide.weights))
      value <- .np_local_operator_ksum(bws, txdat, exdat, operator, matrix(1, nrow(txdat),
        1L), TRUE, "ANN plot support faces")
      list(kw = value$kernel.weights)
    }
    environment(face.owner) <- deps
    faces <- face.owner(st, x, ex, e)
    if (any(!is.finite(c(faces$cut, faces$coefficient))))
      .np_plot_joint_unavailable("nonregular ANN support face")
  }
  for (q in seq_len(m)) {
    phi <- K[, q]
    for (j in seq_len(p)) {
      g <- st$geometry[[j]]
      o <- st$order[[j]]
      w <- (e[q, j] - st$x[o, j]) * D[o, q, j]
      if (pdf)
        w <- w + K[o, q]
      w <- w/(g$radius * g$density.sum)
      d <- numeric(n + 1L)
      agg <- rowsum(c(w, -w), c(g$left, g$right) + 1L, reorder = FALSE)
      d[as.integer(rownames(agg))] <- agg[, 1L]
      phi[o] <- phi[o] + cumsum(d)[seq_len(n)]/n
      if (!is.null(faces))
        for (h in 1:2) {
          k <- 2L * (j - 1L) + h
          phi <- phi + faces$coefficient[q, k] * ((st$x[, j] <= e[q, j]) - (st$x[,
          j] <= faces$cut[q, k]))
        }
    }
    A[q, ] <- (phi - mean(phi))/n * sqrt(n/(n - 1))
  }
  list(estimate = colMeans(K), A = A, target = "same-sample ANN random-radius influence",
    ridge = 0)
}
.np_plot_density_unconditional <- function(bws, x, ex, cdf = FALSE) {
  if (identical(bws[["type", exact = TRUE]], "adaptive_nn") && any(bws[["icon", exact = TRUE]]))
    return(.np_plot_density_ann(bws, x, ex, !cdf))
  kb <- .np_make_kbandwidth_unconditional(bws, x)
  K <- .np_plot_density_matrix(kb, x, ex, rep(if (cdf) "integral" else "normal", ncol(x)))/nrow(x)
  list(estimate = rowSums(K), A = .np_plot_density_center(K), target = "kernel empirical mass; fixed or frozen-query geometry",
    ridge = 0)
}
.np_plot_density_joint <- function(bws, x, ex, se, y = NULL, ey = NULL, cdf = FALSE, gradient = NULL,
  order = 1L, alpha = 0.05, gaussian.draws = 9999L, max.bytes = 128 * 1024^2, estimate = NULL) {
  n <- nrow(x)
  m <- nrow(ex)
  incumbent <- estimate
  if (!.np_plot_density_joint_supported(bws, !is.null(y)))
    .np_plot_joint_unavailable("this density kernel specification is not supported")
  if (length(se) != m || any(!is.finite(se)) || any(se < 0))
    .np_plot_joint_unavailable("marginal uncertainty is unavailable on this grid")
  # Admission budget includes kernel, influence, derivative and calibration arrays with
  # conservative simultaneous lifetime. No unbounded allocation.
  p <- sum(vapply(x, is.numeric, logical(1)))
  bytes <- 8 * as.double(n) * m
  tile <- min(m, max(1L, floor((4 * 1024^2)/(8 * as.double(n) * (8 + p)))))
  # Preparation retains A plus a bounded kernel tile; rescaling has at most three
  # factor-sized arrays. Calibration separately budgets its factor.
  lp.bytes <- 0
  if (!is.null(y) && bws[["xncon", exact = TRUE]] > 0L) {
    spec <- npConditionalRegEngineSpec(bws, where = "Gaussian density plot memory")
    if (!all(spec$degree.engine == 0L)) {
      width <- as.double(npLpBasisNcol(spec$basis.engine, spec$degree.engine))
      # Basis copies, packed systems, accepted inverses and native solve workspace.
      lp.bytes <- 64 * (n + tile) * width + 256 * width^2
    }
  }
  forecast <- max(3 * bytes + 8 * n * tile * (8 + p) + lp.bytes,
                  2 * bytes + 64 * 1024^2)
  if (length(max.bytes) != 1L || !is.finite(max.bytes) || max.bytes <= 0 || forecast >
    max.bytes || !is.finite(forecast))
    .np_plot_joint_unavailable("density influence exceeds the plot memory budget")
  if (!is.null(gradient) && is.null(y))
    stop("Unconditional public gradients are not exposed")
  A <- matrix(NA_real_, m, n)
  estimate <- numeric(m)
  target <- NULL
  progress <- .np_plot_bootstrap_progress_begin(m, "Computing density plot influence",
    unit = "row")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (first in seq.int(1L, m, by = tile)) {
    ix <- seq.int(first, min(m, first + tile - 1L))
    value <- if (is.null(y))
      .np_plot_density_unconditional(bws, x, ex[ix, , drop = FALSE], cdf) else .np_plot_density_conditional(bws, x, y, ex[ix, , drop = FALSE], ey[ix, , drop = FALSE],
      cdf, gradient, order)
    A[ix, ] <- value$A
    estimate[ix] <- value$estimate
    target <- value$target
    progress <- .np_plot_progress_tick(progress, max(ix))
  }
  value <- NULL
  if (!is.null(incumbent)) {
    tolerance <- 0.01 * se + 128 * .Machine$double.eps * pmax(1, abs(incumbent))
    if (length(incumbent) != m || any(!is.finite(estimate)) || any(abs(estimate - incumbent) >
      tolerance))
      stop("density influence does not reproduce the plotted target within 1% of its SE",
        call. = FALSE)
  }
  P <- .np_plot_density_rescale(A, se)
  if (m <= 512L)
    cal <- .np_plot_gaussian_calibrate(tcrossprod(P$factor), alpha = alpha, draws = gaussian.draws) else cal <- .np_plot_gaussian_factor_calibrate(P$factor, alpha = alpha, draws = gaussian.draws)
  cal$target <- target
  cal$marginal.scale <- "incumbent asymptotic standard errors"
  cal$dependence <- "empirical influence correlation"
  cal$forecast.bytes <- forecast
  cal
}

# Called only after the existing interval owner has constructed its payload.  Explicit
# statistical/resource unavailability retains the previous NA joint family; malformed
# inputs or failures of admitted arithmetic remain errors.
.np_plot_density_interval <- function(interval, bws, x, ex, se, estimate, y = NULL, ey = NULL,
  cdf = FALSE, gradient = NULL, gradient.order = 1L, alpha = 0.05, band.type = "all", gaussian.draws = 9999L) {
  order <- 1L
  if (!is.null(gradient) && bws[["ixcon", exact = TRUE]][gradient]) {
    order <- rep_len(gradient.order, sum(bws[["ixcon", exact = TRUE]]))[match(gradient,
      which(bws[["ixcon", exact = TRUE]]))]
  }
  calibration <- .np_plot_joint_try(.np_plot_density_joint(bws, x, ex, se, y, ey, cdf,
    gradient, order, alpha, gaussian.draws, estimate = estimate))
  if (!is.null(calibration)) {
    simultaneous <- cbind(se, se, deparse.level = 0L) * calibration$critical
    if (identical(band.type, "all"))
      interval$all.err$simultaneous <- simultaneous else interval$err[, 1:2] <- simultaneous
    interval$calibration <- calibration
  }
  interval
}

.np_plot_density_band_fields <- function(object, interval, gradients = FALSE, component = 1L) {
  .np_plot_joint_band_fields(object, interval$calibration, interval$all.err,
                            if (gradients) component else NULL)
}
