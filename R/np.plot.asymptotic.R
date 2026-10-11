# Finite-grid Gaussian calibration. These helpers have no estimator callers.
.np_plot_gaussian_integer <- function(value, name, maximum) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != floor(value) || value > maximum)
    stop(sprintf("%s must be an integer in [1, %d]", name, maximum),
         call. = FALSE)
  as.integer(value)
}

.np_plot_gaussian_calibrate <- function(covariance, alpha = 0.05,
                                        draws = 50000L, block = 4096L,
                                        max.bytes = 64 * 1024^2) {
  draws <- .np_plot_gaussian_integer(draws, "gaussian.draws", 1000000L)
  block <- .np_plot_gaussian_integer(block, "Gaussian draw block", 1000000L)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5)
    stop("alpha must lie in (0, 0.5)", call. = FALSE)
  if (!is.numeric(max.bytes) || length(max.bytes) != 1L ||
      !is.finite(max.bytes) || max.bytes <= 0)
    stop("Gaussian memory budget must be positive and finite", call. = FALSE)
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      nrow(covariance) < 1L || nrow(covariance) != ncol(covariance))
    stop("Gaussian covariance must be a nonempty square numeric matrix",
         call. = FALSE)
  g <- nrow(covariance)
  if (g > 512L || 64 * as.double(g)^2 + 8 * draws + 64 * g > max.bytes)
    stop("finite-grid Gaussian covariance exceeds the plot memory/grid budget",
         call. = FALSE)
  if (any(!is.finite(covariance)))
    stop("Gaussian covariance must be finite", call. = FALSE)
  tolerance <- 64 * .Machine$double.eps * g * max(abs(covariance))
  if (max(abs(covariance - t(covariance))) > tolerance)
    stop("Gaussian covariance is materially asymmetric", call. = FALSE)
  variance <- diag(covariance)
  if (any(variance < 0))
    stop("Gaussian covariance has a negative variance", call. = FALSE)
  active <- which(variance > 0)
  zero <- which(variance == 0)
  if (length(zero) && (any(covariance[zero, , drop = FALSE] != 0) ||
                       any(covariance[, zero, drop = FALSE] != 0)))
    stop("zero Gaussian variance requires a zero covariance row and column",
         call. = FALSE)
  metadata <- list(method = "asymptotic Gaussian", scope = "finite grid within panel",
                   draws = draws, critical = 0, grid.size = g,
                   positive.variance = length(active), rank = 0L)
  if (!length(active)) return(metadata)
  s <- sqrt(variance[active])
  correlation <- covariance[active, active, drop = FALSE] / outer(s, s)
  correlation <- (correlation + t(correlation)) / 2
  eig <- eigen(correlation, symmetric = TRUE)
  tolerance <- 64 * .Machine$double.eps * length(active) * max(abs(correlation))
  if (min(eig$values) < -tolerance)
    stop("Gaussian covariance is materially indefinite", call. = FALSE)
  values <- pmax(eig$values, 0)
  # The symmetric square root avoids eigenvector sign conventions. Every
  # column of the normal matrix is one draw in canonical coordinate order.
  root <- tcrossprod(sweep(eig$vectors, 2L, sqrt(values), "*"), eig$vectors)
  block <- min(block, draws, floor((max.bytes - 64 * as.double(g)^2 -
                                     8 * draws) / (64 * g)))
  maxima <- numeric(draws)
  for (start in seq.int(1L, draws, by = block)) {
    last <- min(draws, start + block - 1L)
    normal <- matrix(rnorm(length(active) * (last - start + 1L)),
                     nrow = length(active))
    simulated <- root %*% normal
    maxima[start:last] <- apply(abs(simulated), 2L, max)
  }
  metadata$critical <- unname(quantile(maxima, 1 - alpha, type = 1L))
  metadata$rank <- sum(values > tolerance)
  metadata
}

.np_plot_regression_joint_covariance <- function(bws, xdat, ydat, exdat,
                                                residual = NULL,
                                                block.bytes = 4 * 1024^2,
                                                max.bytes = 64 * 1024^2) {
  n <- nrow(xdat)
  g <- nrow(exdat)
  ydat <- .np_hat_response(ydat, bws[["ydati", exact = TRUE]])
  if (!is.numeric(block.bytes) || length(block.bytes) != 1L ||
      !is.finite(block.bytes) || block.bytes <= 0 ||
      !is.numeric(max.bytes) || length(max.bytes) != 1L ||
      !is.finite(max.bytes) || max.bytes <= 0)
    stop("joint covariance memory budgets must be positive and finite",
         call. = FALSE)
  spec <- npValidatedConditionalRegSpec(bws, where = "plot joint covariance",
                                         ncon.field = "ncon")
  width <- if (identical(spec$regtype.engine, "lc")) 1 else
    npLpBasisNcol(spec$basis.engine, spec$degree.engine)
  # Conservative workspace allowance includes the unchanged operator's basis,
  # solve, coordinate copies and tree, as well as the plot-owned arrays.
  setup.bytes <- 64 * n * (width + ncol(xdat) + 1) + 128 * width^2
  available <- min(block.bytes, max.bytes - 64 * as.double(g)^2 -
                     16 * n - setup.bytes)
  rows <- floor(available / (64 * as.double(n)))
  if (n < 2L || g < 1L || g > 512L || rows < 1 ||
      length(ydat) != n || any(!is.finite(ydat)) ||
      anyNA(xdat) || anyNA(exdat))
    stop("joint covariance has invalid data or exceeds the plot memory/grid budget",
         call. = FALSE)
  rows <- min(rows, n, g)
  operator <- function(evaluation) {
    npreghat(bws = bws, txdat = xdat, exdat = evaluation,
             output = "matrix", .np.require.finite = TRUE)
  }
  z <- residual
  if (is.null(z)) {
  z <- numeric(n)
  for (start in seq.int(1L, n, by = rows)) {
    index <- seq.int(start, min(n, start + rows - 1L))
    hat <- operator(xdat[index, , drop = FALSE])
    ridge <- attr(hat, "ridge.used", exact = TRUE)
    if (is.null(ridge) || length(ridge) != length(index))
      stop("joint covariance requires the influence owner's ridge certificate",
           call. = FALSE)
    for (j in seq_along(index)) {
      i <- index[j]
      if (ridge[j] == 0) {
        weights <- hat[j, -i]
        scale <- max(abs(weights))
        if (scale == 0)
          stop("joint covariance has unidentified training residual information",
               call. = FALSE)
        weights <- weights / scale
        z[i] <- sum(weights * (ydat[i] - ydat[-i])) /
          sqrt(sum(weights)^2 + sum(weights^2))
      } else {
        residual.row <- -hat[j, ]
        residual.row[i] <- 1 + residual.row[i]
        scale <- max(abs(residual.row))
        if (scale == 0)
          stop("joint covariance has unidentified training residual information",
               call. = FALSE)
        residual.row <- residual.row / scale
        z[i] <- sum(residual.row * ydat) / sqrt(sum(residual.row^2))
      }
    }
  }
  }
  if (length(z) != n || any(!is.finite(z)))
    stop("joint covariance normalized residuals are nonfinite", call. = FALSE)
  covariance <- matrix(0, g, g)
  starts <- seq.int(1L, g, by = rows)
  for (a in seq_along(starts)) {
    ia <- seq.int(starts[a], min(g, starts[a] + rows - 1L))
    left <- sweep(operator(exdat[ia, , drop = FALSE]), 2L, z, "*")
    for (b in seq_len(a)) {
      ib <- seq.int(starts[b], min(g, starts[b] + rows - 1L))
      right <- if (a == b) left else
        sweep(operator(exdat[ib, , drop = FALSE]), 2L, z, "*")
      cell <- tcrossprod(left, right)
      covariance[ia, ib] <- cell
      covariance[ib, ia] <- t(cell)
    }
  }
  list(covariance = covariance, residual = z, block.rows = rows)
}

.np_plot_regression_joint_supported <- function(bws, gradients = FALSE) {
  !isTRUE(gradients) && identical(bws[["type", exact = TRUE]], "fixed") &&
    bws[["ckertype", exact = TRUE]] %in% c("gaussian", "epanechnikov") &&
    identical(as.integer(bws[["ckerorder", exact = TRUE]]), 2L) &&
    isTRUE(bws[["ncon", exact = TRUE]] > 0L)
}

.np_plot_scoef_joint_covariance <- function(bws, xdat, ydat, zdat,
                                           exdat, ezdat, residual,
                                           block.bytes = 4 * 1024^2,
                                           max.bytes = 64 * 1024^2) {
  n <- nrow(xdat)
  g <- nrow(exdat)
  ydat <- .np_hat_response(ydat, bws[["ydati", exact = TRUE]])
  separate.z <- !is.null(zdat)
  if (!separate.z) zdat <- xdat
  if (is.null(ezdat)) ezdat <- exdat
  spec <- .npscoef_canonical_spec(bws, zdat, where = "plot joint covariance")
  width <- (ncol(xdat) + 1) * if (identical(spec$regtype.engine, "lc")) 1 else
    npLpBasisNcol(spec$basis.engine, spec$degree.engine)
  if (!is.numeric(block.bytes) || length(block.bytes) != 1L ||
      !is.finite(block.bytes) || block.bytes <= 0 ||
      !is.numeric(max.bytes) || length(max.bytes) != 1L ||
      !is.finite(max.bytes) || max.bytes <= 0)
    stop("joint covariance memory budgets must be positive and finite",
         call. = FALSE)
  setup.bytes <- 64 * n * (width + ncol(zdat) + 1) + 128 * width^2
  available <- min(block.bytes, max.bytes - 64 * as.double(g)^2 -
                     16 * n - setup.bytes)
  rows <- floor(available / (64 * as.double(n)))
  if (n < 2L || g < 1L || g > 512L || rows < 1 ||
      nrow(zdat) != n || nrow(ezdat) != g ||
      length(ydat) != n || any(!is.finite(ydat)) ||
      anyNA(xdat) || anyNA(exdat) || anyNA(zdat) || anyNA(ezdat))
    stop("joint covariance has invalid data or exceeds the plot memory/grid budget",
         call. = FALSE)
  rows <- min(rows, n, g)
  operator <- function(x, z) {
    args <- list(bws = bws, txdat = xdat, exdat = x,
                 output = "matrix", iterate = FALSE)
    if (separate.z) {
      args$tzdat <- zdat
      args$ezdat <- z
    }
    hat <- do.call(npscoefhat, args)
    if (any(!is.finite(hat)))
      stop("joint covariance requires a finite smooth-coefficient influence map",
           call. = FALSE)
    hat
  }
  z <- residual
  if (length(z) != n || any(!is.finite(z)))
    stop("joint covariance smooth-coefficient residuals are nonfinite", call. = FALSE)
  covariance <- matrix(0, g, g)
  mean <- numeric(g)
  starts <- seq.int(1L, g, by = rows)
  for (a in seq_along(starts)) {
    ia <- seq.int(starts[a], min(g, starts[a] + rows - 1L))
    hat <- operator(exdat[ia, , drop = FALSE], ezdat[ia, , drop = FALSE])
    mean[ia] <- drop(hat %*% ydat)
    left <- sweep(hat, 2L, z, "*")
    for (b in seq_len(a)) {
      ib <- seq.int(starts[b], min(g, starts[b] + rows - 1L))
      right <- if (a == b) left else
        sweep(operator(exdat[ib, , drop = FALSE], ezdat[ib, , drop = FALSE]),
              2L, z, "*")
      cell <- tcrossprod(left, right)
      covariance[ia, ib] <- cell
      covariance[ib, ia] <- t(cell)
    }
  }
  list(covariance = covariance, residual = z, block.rows = rows, mean = mean)
}
