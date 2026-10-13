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
# Some density fitting owners use native near-upper categorical shortcuts;
# generic kernel matrices deliberately retain raw kernels. Match that owner only
# where requested by the density plot policy, without changing shared kernels.
.np_plot_density_factor_weights <- function(K, kbw, x, ex, operator, ordered = TRUE,
  constant.response = FALSE) {
  constant <- if (constant.response && all(kbw[["iuno", exact = TRUE]] &
    operator == "normal")) 1 else NA_real_
  for (j in which(kbw[["iuno", exact = TRUE]] & operator == "normal")) {
    nc <- length(kbw[["xdati", exact = TRUE]]$all.dlev[[j]])
    if (nc < 2L)
      next
    lambda <- kbw[["bw", exact = TRUE]][j]
    ku <- if (identical(kbw[["ukertype", exact = TRUE]], "aitchisonaitken")) 0L else 1L
    probe <- matrix(c(1, 2), 2L, 1L)
    empty <- matrix(numeric(), 2L, 0L)
    effective <- .np_categorical_profile_kernel_tile(probe, empty, probe, empty,
      ku, integer(), 0L, lambda, nc, list(NULL))[1L, ]
    if (is.finite(constant))
      constant <- if (effective[1L] == effective[2L]) constant * effective[1L] else NA_real_
    # Literal normalized kernels, used solely to remove their factor from K.
    # The native tile owns the shortcut predicate, options and effective values.
    raw <- if (ku == 0L) c(1 - lambda, lambda/(nc - 1L)) else
      c(1, lambda)/(1 + (nc - 1L) * lambda)
    if (identical(effective, raw))
      next
    if (any(raw == 0 & effective != 0))
      stop("Unavailable categorical kernel correction")
    ratio <- rep(1, 2L)
    active <- raw != 0
    ratio[active] <- effective[active]/raw[active]
    different <- outer(as.integer(ex[[j]]), match(x[[j]], levels(ex[[j]])), "!=")
    K <- K * matrix(ratio[1L + different], nrow(K), ncol(K))
  }
  # Reproduce a certified donor-independent response exactly, avoiding the
  # roundoff introduced by removing raw factors and multiplying them back.
  if (is.finite(constant))
    K[,] <- constant
  ordered <- if (ordered) which(kbw[["iord", exact = TRUE]] &
    operator %in% c("normal", "integral")) else integer()
  if (length(ordered) && !identical(kbw[["okertype", exact = TRUE]], "wangvanryzin") &&
      npLogicalOption("np.largelambda", TRUE)) {
    tol <- npDiscUpperRelTol()
    for (j in ordered) {
      lambda <- kbw[["bw", exact = TRUE]][j]
      # Match the existing product-owner certificate in jksum.c:
      # np_disc_ordered_near_upper / np_disc_near_const_kernel.
      # The scalar ordered tile exposes raw values, not that product shortcut.
      if (abs(lambda - 1) > max(tol, 16 * .Machine$double.eps))
        next
      cats <- as.double(kbw[["xdati", exact = TRUE]]$all.dlev[[j]])
      nc <- length(cats)
      ko <- if (identical(kbw[["okertype", exact = TRUE]], "racineliyan")) 3L else 2L
      op <- if (operator[j] == "integral") 3L else 0L
      tile <- if (ko == 3L && op == 0L) function(train, eval)
        t(.Call(C_np_ordered_rly_matrix, as.double(train), as.double(eval),
          as.double(lambda), cats)) else function(train, eval)
        .np_categorical_profile_kernel_tile(
          matrix(numeric(), length(train), 0L), matrix(train, ncol = 1L),
          matrix(numeric(), length(eval), 0L), matrix(eval, ncol = 1L),
          integer(), ko, op, lambda, nc, list(cats))
      if (ko == 3L && op == 3L) {
        # With >=2 categories the retained-support CDF range is at least 1/2:
        # the last-donor mass at the first category is <=1/nc, while F(last)=1.
        # The canonical tolerance is strictly <1/2, so no shortcut is possible.
        # With one category the raw kernel already equals its constant value.
        next
      } else if (ko == 3L) {
        # The native range certificate uses each donor's farthest endpoint
        # and diagonal. Reuse the cached-normalizer bridge, in bounded tiles,
        # avoiding a cubic scan of all donor/evaluation/support triples.
        bounds <- c(min(tile(cats, cats[c(1L, nc)])), -Inf)
        for (first in seq.int(1L, nc, by = 512L)) {
          ix <- seq.int(first, min(nc, first + 511L))
          bounds[2L] <- max(bounds[2L], diag(tile(cats[ix], cats[ix])))
        }
      } else bounds <- as.vector(tile(cats[1L], cats[c(1L, nc)]))
      if (abs(diff(bounds)) > tol * max(1, abs(bounds)) || bounds[1L] == bounds[2L])
        next
      tx <- as.double(toMatrix(x[j]))
      eq <- as.double(toMatrix(ex[j]))
      raw <- t(tile(tx, eq))
      K <- K * (mean(bounds)/raw)
    }
  }
  K
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
# Batch direct response and adjoint solves through the unchanged canonical
# ranked owner; do not form an inverse and multiply it by response moments.
.np_plot_density_lp_score <- function(W, w, R, e, n, identity) {
  p <- ncol(W); R <- as.matrix(R); nr <- ncol(R)
  M <- crossprod(W, w * W)
  packed <- matrix(unlist(lapply(seq_len(p), function(a) M[a, a:p])), 1L)
  rank <- .np_inid_lp_rank_bounds(w, p, identity = identity)
  # Intercept equivariance holds also for the accepted ridge correction:
  # C0 W' diag(w) (r+c) = C0 W' diag(w) r + c (M+ridge I) e1.
  # Centering preserves that estimator and avoids cancellation at CDF endpoints.
  center <- R[which.max(w), ]
  R <- R - rep(center, each = n)
  B <- crossprod(W, w * R)
  U <- diag(p)
  Z <- cbind(B, e, U[, 1L])
  Z <- lapply(seq_len(p), function(a) matrix(Z[a, ], 1L))
  C <- matrix(NA_real_, p, nr + 2L)
  ridge <- steps <- NULL
  for (a in seq_len(p)) {
    ans <- .np_inid_lp_batch_project(packed, Z, U[a, ], n, TRUE, rank)
    C[a, ] <- ans$values[1L, ]
    if (a == 1L) {
      ridge <- ans$ridge_total[1L]; steps <- ans$ridge_steps[1L]
    } else stopifnot(identical(ans$ridge_total[1L], ridge),
                     identical(ans$ridge_steps[1L], steps))
  }
  anchor <- M[1L, 1L]
  c0 <- if (ridge == 0) 1 else 1 + ridge / anchor
  theta <- C[, seq_len(nr), drop = FALSE]
  v <- C[, nr + 1L] - (1 - 1 / c0) * e[1L] * C[, nr + 2L]
  # The existing plot-only derivative of the accepted ridge scale, unchanged.
  dr <- numeric(n)
  if (ridge > 0) {
    d <- abs(diag(M)); mask <- diag(p) == 1; mx <- max(d)
    if (mx == 0) { d <- abs(M); mask <- matrix(TRUE, p, p); mx <- max(d) }
    ix <- which(abs(M) == mx & mask, arr.ind = TRUE)
    if (!nrow(ix)) stop("Unidentified ridge scale")
    da <- sign(M[ix[1, 1], ix[1, 2]]) * w * W[, ix[1, 1]] * W[, ix[1, 2]]
    for (k in seq_len(nrow(ix)))
      if (any(da != sign(M[ix[k, 1], ix[k, 2]]) * w * W[, ix[k, 1]] * W[, ix[k, 2]]))
        stop("Nonsmooth accepted ridge scale; influence not defined")
    dr <- steps / n * da
  }
  vz <- drop(W %*% v); zt <- W %*% theta
  A <- (w * vz) * (R - zt)
  if (ridge > 0) {
    rho <- ridge / anchor; da <- w * W[, 1L]^2
    A <- A + (w * v[1L] * rho * W[, 1L]) * R +
      outer(v[1L] * (dr / anchor - rho * da / anchor), B[1L, ]) -
      outer(dr, drop(crossprod(v, theta)))
  }
  list(estimate = drop(crossprod(e, theta)) + center * e[1L], score = t(A), ridge = rep(ridge, nr))
}
.np_plot_density_conditional <- function(bws, x, y, ex, ey, cdf = FALSE, gradient = NULL,
  order = 1L, factor.owner = NULL) {
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
      spec <- npConditionalRegEngineSpec(bws, where = "Gaussian density contrast")
      # Match npcdhat's current owner: LC contrasts gate unordered weights;
      # LP hats use raw weights except GNN degree > 1, which uses the core.
      # These need not equal differences of the ordinary gated LP level fits.
      gated <- if (identical(spec$reg.engine, "lc") || all(spec$degree.engine == 0L))
        "unordered" else if (identical(bws[["type", exact = TRUE]], "generalized_nn") &&
          any(spec$degree.engine > 1L)) "all" else "none"
      u <- .np_plot_density_conditional(bws, x, y, frames$upper, ey, cdf, factor.owner = gated)
      l <- .np_plot_density_conditional(bws, x, y, frames$lower, ey, cdf, factor.owner = gated)
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
  # LC continuous gradients use the native permutation owner. Its unordered
  # explanatory factors are gated in both numerator and denominator (including
  # derivative products); ordered factors retain their literal kernel values.
  # The tile bridge observes that owner's effective native state, including
  # its cached options; it does not create a second selection policy.
  lc.gradient <- lc && !is.null(gradient)
  if (is.null(factor.owner))
    factor.owner <- if (lc.gradient) "unordered" else if (lc) "none" else "all"
  if (!identical(factor.owner, "none"))
    K <- .np_plot_density_factor_weights(K, xb, x, ex, rep("normal", ncol(x)),
      ordered = identical(factor.owner, "all"))
  if (lc.gradient) {
    yop <- rep(if (cdf) "integral" else "normal", ncol(y))
    # Unordered CDF factors use the normal categorical kernel, as in the
    # native OP_UFUN_OFFSETS[OP_INTEGRAL] == 0 contract.
    yop[yb[["iuno", exact = TRUE]]] <- "normal"
    Y <- .np_plot_density_factor_weights(Y, yb, y, ey, yop,
      ordered = FALSE, constant.response = TRUE)
  }
  if (lc) {
    den <- rowSums(K)
    # Translation cancels from a ratio derivative. For categorical responses
    # it also preserves exact zero influence when the effective Y is constant.
    if (lc.gradient && any(yb[["iuno", exact = TRUE]]))
      Y <- Y - Y[, 1L]
    mu <- rowSums(K * Y)/den
    A <- K * (Y - mu)/den
    if (!is.null(gradient)) {
      if (order != 1L)
        stop("LC supports first continuous derivatives only")
      D <- .np_plot_density_derivative_matrix(xb, x, ex, gradient)
      dop <- rep("normal", ncol(x)); dop[gradient] <- "derivative"
      D <- .np_plot_density_factor_weights(D, xb, x, ex, dop, ordered = FALSE)
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
  identity <- .np_inid_lp_design_identity(W)
  E <- as.matrix(do.call(W.lp, c(basisargs, list(exdat = ex[, cont, drop = FALSE], gradient.vec = grad))))
  groups <- .np_plot_exact_row_groups(ex)$groups
  A <- matrix(NA_real_, m, n)
  mu <- ridge <- numeric(m)
  for (ix in groups) {
    w <- K[ix[1L], ]
    v <- .np_plot_density_lp_score(W, w, t(Y[ix, , drop = FALSE]), E[ix[1L], ], n, identity)
    mu[ix] <- v$estimate
    A[ix, ] <- v$score
    ridge[ix] <- v$ridge
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
  local.matrix <- function(op) {
    value <- .np_local_operator_ksum(kb, x, ex, op, matrix(1, n, 1L),
      return.kernel.weights = TRUE, where = "ANN plot influence")$kernel.weights
    if (any(kb[["iuno", exact = TRUE]] & op == "normal") ||
        any(kb[["iord", exact = TRUE]] & op %in% c("normal", "integral")))
      value <- t(.np_plot_density_factor_weights(t(value), kb, x, ex, op))
    value
  }
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
  op <- rep(if (cdf) "integral" else "normal", ncol(x))
  K <- .np_plot_density_matrix(kb, x, ex, op)
  K <- .np_plot_density_factor_weights(K, kb, x, ex, op)
  K <- K/nrow(x)
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
  lp <- FALSE
  if (!is.null(y) && bws[["xncon", exact = TRUE]] > 0L) {
    spec <- npConditionalRegEngineSpec(bws, where = "Gaussian density plot memory")
    if (!all(spec$degree.engine == 0L)) {
      lp <- TRUE
      width <- as.double(npLpBasisNcol(spec$basis.engine, spec$degree.engine))
      # Basis copies, packed systems, batched response/adjoint and native solve workspace.
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
  # Surfaces revisit the same X at many Y values. Group before bounded tiling
  # so each solve can handle several responses, then scatter to the original grid.
  schedule <- if (lp) unlist(.np_plot_exact_row_groups(ex)$groups, use.names = FALSE) else seq_len(m)
  for (first in seq.int(1L, m, by = tile)) {
    last <- min(m, first + tile - 1L)
    ix <- schedule[seq.int(first, last)]
    value <- if (is.null(y))
      .np_plot_density_unconditional(bws, x, ex[ix, , drop = FALSE], cdf) else .np_plot_density_conditional(bws, x, y, ex[ix, , drop = FALSE], ey[ix, , drop = FALSE],
      cdf, gradient, order)
    A[ix, ] <- value$A
    estimate[ix] <- value$estimate
    target <- value$target
    progress <- .np_plot_progress_tick(progress, last)
  }
  value <- NULL
  if (!is.null(incumbent)) {
    if (length(incumbent) != m || any(!is.finite(incumbent)) || any(!is.finite(estimate)))
      stop("Invalid density plot target", call. = FALSE)
    # Accumulated rounding is independent of sampling uncertainty. In particular,
    # CDF endpoint SEs can be zero. This allowance changes only the comparison,
    # never the incumbent curve, SEs, influence rows or calibrated covariance.
    roundoff <- max(128, 8 * as.double(n)) * .Machine$double.eps *
      pmax(1, abs(incumbent), abs(estimate))
    tolerance <- 0.01 * se + roundoff
    if (any(abs(estimate - incumbent) > tolerance))
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
