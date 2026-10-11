# A missing joint interval must retain one row per retained plotted center.
# Other interval routes preserve their established omission behavior.
.np_plot_joint_interval_rows <- function(errors, center, joint, type) {
  if (joint && identical(type, "simultaneous"))
    errors[!is.na(center), , drop = FALSE]
  else
    na.omit(errors)
}

.np_plot_joint_panel_range <- function(estimate, err, plotOnEstimate = TRUE) {
  limits <- .np_plot_panel_error_range(estimate, err,
    plot.errors.type = "simultaneous", plotOnEstimate = plotOnEstimate)
  if (all(is.finite(limits))) limits else c(-1, 1)
}

.np_plot_joint_unavailable <- function(message, call. = FALSE) {
  stop(structure(list(message = message, call = NULL),
                 class = c("np_plot_joint_unavailable", "error", "condition")))
}

# Only declared joint-band limitations are recoverable here. Invalid arguments,
# programming errors, invalid covariance objects and interrupts still propagate.
.np_plot_joint_try <- function(expr, unavailable = NULL) {
  expression <- substitute(expr)
  context <- parent.frame()
  tryCatch(eval(expression, context), np_plot_joint_unavailable = function(e) {
    .np_warning(paste0("asymptotic simultaneous variability bands are unavailable: ",
      conditionMessage(e), "; returning NA simultaneous limits and retaining any ",
      "pointwise/Bonferroni intervals."), call. = FALSE)
    unavailable
  })
}

# Finite-grid Gaussian calibration. These helpers have no estimator callers.
.np_plot_gaussian_integer <- function(value, name, maximum = 2^53 - 1) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != floor(value) || value > maximum)
    stop(sprintf("%s must be an integer in [1, %.0f]", name, maximum),
         call. = FALSE)
  if (value <= .Machine$integer.max) as.integer(value) else as.double(value)
}

.np_plot_gaussian_draws <- function(draws, alpha, m) {
  # The simulated statistic is one maximum: its upper tail has mass alpha,
  # irrespective of grid size. This is a probability-scale Monte Carlo target,
  # not an endpoint-accuracy or statistical-coverage guarantee.
  minimum <- ceiling(qnorm(.975)^2 * (1 - alpha) / (.1^2 * alpha))
  draws <- .np_plot_gaussian_integer(draws, "gaussian.draws")
  if (draws < minimum)
    .np_warning(sprintf(paste0("gaussian.draws=%.0f is below the recommended Monte Carlo ",
      "precision count %.0f (alpha=%g, m=%.0f evaluation points). This approximate ",
      "95%% target limits relative simulation error in the maximum's upper-tail ",
      "probability to 10%%; it does not guarantee endpoint accuracy or statistical ",
      "coverage. Increase gaussian.draws for greater precision; retaining the requested count."),
      draws, minimum, alpha, m), call. = FALSE)
  draws
}

# Surface grids use the observation-level factor directly: A %*% e has
# covariance A %*% t(A), without allocating or diagonalizing a g-by-g matrix.
.np_plot_gaussian_factor_calibrate <- function(influence, alpha = 0.05,
                                               draws = 9999L, block = 4096L,
                                               max.bytes = 64 * 1024^2) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5)
    stop("alpha must lie in (0, 0.5)", call. = FALSE)
  if (!is.matrix(influence) || !is.numeric(influence) ||
      nrow(influence) < 1L || ncol(influence) < 1L || any(!is.finite(influence)))
    stop("Gaussian influence factor must be a finite nonempty numeric matrix",
         call. = FALSE)
  if (!is.numeric(max.bytes) || length(max.bytes) != 1L ||
      !is.finite(max.bytes) || max.bytes <= 0)
    stop("Gaussian memory budget must be positive and finite", call. = FALSE)
  g <- nrow(influence)
  n <- ncol(influence)
  draws <- .np_plot_gaussian_draws(draws, alpha, g)
  block <- .np_plot_gaussian_integer(block, "Gaussian draw block", 1000000L)
  available <- max.bytes - 8 * as.double(g) * n - 24 * draws - 64 * g
  block <- min(block, draws, floor(available / (64 * (as.double(g) + n))))
  if (block < 1)
    .np_plot_joint_unavailable("finite-grid Gaussian influence factor exceeds the plot memory budget",
         call. = FALSE)
  # Row sums avoid a second full factor-sized squared array.
  se <- vapply(seq_len(g), function(i) sqrt(sum(influence[i, ]^2)), numeric(1L))
  if (any(!is.finite(se)))
    stop("Gaussian influence factor has nonfinite variance", call. = FALSE)
  active <- se > 0
  metadata <- list(method = "asymptotic Gaussian", scope = "finite grid within panel",
                   draws = draws, critical = 0, grid.size = g,
                   positive.variance = sum(active), rank = NA_integer_,
                   representation = "observation influence factor")
  if (!any(active)) return(metadata)
  divisor <- ifelse(active, se, 1)
  maxima <- numeric(draws)
  progress <- .np_plot_bootstrap_progress_begin(draws, "Gaussian plot calibration",
                                                 unit = "draw")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (start in seq.int(1L, draws, by = block)) {
    last <- min(draws, start + block - 1L)
    normal <- matrix(rnorm(n * (last - start + 1L)), nrow = n)
    simulated <- (influence %*% normal) / divisor
    maxima[start:last] <- apply(abs(simulated), 2L, max)
    progress <- .np_plot_progress_tick(progress, last)
  }
  metadata$critical <- unname(quantile(maxima, 1 - alpha, type = 1L))
  metadata
}

.np_plot_gaussian_calibrate <- function(covariance, alpha = 0.05,
                                        draws = 9999L, block = 4096L,
                                        max.bytes = 64 * 1024^2) {
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
  draws <- .np_plot_gaussian_draws(draws, alpha, g)
  if (g > 512L || 64 * as.double(g)^2 + 24 * draws + 64 * g > max.bytes)
    .np_plot_joint_unavailable("finite-grid Gaussian covariance exceeds the plot memory/grid budget",
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
                                     24 * draws) / (64 * g)))
  maxima <- numeric(draws)
  progress <- .np_plot_bootstrap_progress_begin(draws, "Gaussian plot calibration",
                                                 unit = "draw")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (start in seq.int(1L, draws, by = block)) {
    last <- min(draws, start + block - 1L)
    normal <- matrix(rnorm(length(active) * (last - start + 1L)),
                     nrow = length(active))
    simulated <- root %*% normal
    maxima[start:last] <- apply(abs(simulated), 2L, max)
    progress <- .np_plot_progress_tick(progress, last)
  }
  metadata$critical <- unname(quantile(maxima, 1 - alpha, type = 1L))
  metadata$rank <- sum(values > tolerance)
  metadata
}

.np_plot_regression_joint_covariance <- function(bws, xdat, ydat, exdat,
                                                residual = NULL, deriv = NULL,
                                                categorical.index = NULL,
                                                factor = FALSE,
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
  joint.bytes <- if (factor) 16 * as.double(g) * n else 64 * as.double(g)^2
  available <- min(block.bytes, max.bytes - joint.bytes -
                     16 * n - setup.bytes)
  rows <- floor(available / (64 * as.double(n)))
  if (n < 2L || g < 1L ||
      length(ydat) != n || any(!is.finite(ydat)) ||
      anyNA(xdat) || anyNA(exdat))
    stop("joint covariance has invalid data", call. = FALSE)
  if ((!factor && g > 512L) || rows < 1)
    .np_plot_joint_unavailable("joint covariance exceeds the plot memory/grid budget")
  rows <- min(rows, n, g)
  operator <- function(evaluation, derivative = NULL) {
    npreghat(bws = bws, txdat = xdat, exdat = evaluation, deriv = derivative,
             output = "matrix", .np.require.finite = TRUE)
  }
  z <- residual
  if (is.null(z) && bws[["type", exact = TRUE]] %in% c("generalized_nn", "adaptive_nn"))
    z <- .np_plot_response_residuals(bws, xdat, ydat, rows)$residual
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
          .np_plot_joint_unavailable("joint covariance has unidentified training residual information",
               call. = FALSE)
        weights <- weights / scale
        z[i] <- sum(weights * (ydat[i] - ydat[-i])) /
          sqrt(sum(weights)^2 + sum(weights^2))
      } else {
        residual.row <- -hat[j, ]
        residual.row[i] <- 1 + residual.row[i]
        scale <- max(abs(residual.row))
        if (scale == 0)
          .np_plot_joint_unavailable("joint covariance has unidentified training residual information",
               call. = FALSE)
        residual.row <- residual.row / scale
        z[i] <- sum(residual.row * ydat) / sqrt(sum(residual.row^2))
      }
    }
  }
  }
  if (length(z) != n || any(!is.finite(z)))
    .np_plot_joint_unavailable("joint covariance normalized residuals are nonfinite", call. = FALSE)
  evaluation.operator <- operator
  if (!is.null(categorical.index)) {
    if (!is.null(deriv))
      stop("joint covariance cannot combine continuous and categorical derivatives",
           call. = FALSE)
    evaluation.operator <- function(evaluation, derivative = NULL) {
      frames <- npCategoricalFirstDifferenceFrames(evaluation, categorical.index,
                                                   "plot joint covariance")
      operator(frames$upper) - operator(frames$lower)
    }
  }
  covariance <- if (factor) NULL else matrix(0, g, g)
  influence <- if (factor) matrix(0, g, n) else NULL
  estimate <- if (!is.null(deriv) || !is.null(categorical.index)) numeric(g) else NULL
  starts <- seq.int(1L, g, by = rows)
  for (a in seq_along(starts)) {
    ia <- seq.int(starts[a], min(g, starts[a] + rows - 1L))
    hat <- evaluation.operator(exdat[ia, , drop = FALSE], deriv)
    if (!is.null(estimate)) estimate[ia] <- drop(hat %*% ydat)
    left <- sweep(hat, 2L, z, "*")
    if (factor) {
      influence[ia, ] <- left
      next
    }
    for (b in seq_len(a)) {
      ib <- seq.int(starts[b], min(g, starts[b] + rows - 1L))
      right <- if (a == b) left else
        sweep(evaluation.operator(exdat[ib, , drop = FALSE], deriv), 2L, z, "*")
      cell <- tcrossprod(left, right)
      covariance[ia, ib] <- cell
      covariance[ib, ia] <- t(cell)
    }
  }
  list(covariance = covariance, residual = z, block.rows = rows,
       estimate = estimate, influence = influence)
}

.np_plot_regression_joint_supported <- function(bws) {
  bws[["type", exact = TRUE]] %in% c("fixed", "generalized_nn", "adaptive_nn") &&
    (isTRUE(bws[["ncon", exact = TRUE]] == 0L) ||
     (bws[["ckertype", exact = TRUE]] %in% c("gaussian", "epanechnikov") &&
      identical(as.integer(bws[["ckerorder", exact = TRUE]]), 2L)))
}

.np_plot_scoef_joint_covariance <- function(bws, xdat, ydat, zdat,
                                           exdat, ezdat, residual,
                                           gradient.index = NULL, factor = FALSE,
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
  joint.bytes <- if (factor) 16 * as.double(g) * n else 64 * as.double(g)^2
  available <- min(block.bytes, max.bytes - joint.bytes -
                     16 * n - setup.bytes)
  rows <- floor(available / (64 * as.double(n)))
  if (n < 2L || g < 1L ||
      nrow(zdat) != n || nrow(ezdat) != g ||
      length(ydat) != n || any(!is.finite(ydat)) ||
      anyNA(xdat) || anyNA(exdat) || anyNA(zdat) || anyNA(ezdat))
    stop("joint covariance has invalid data", call. = FALSE)
  if ((!factor && g > 512L) || rows < 1)
    .np_plot_joint_unavailable("joint covariance exceeds the plot memory/grid budget")
  if (!is.null(gradient.index))
    rows <- min(rows, floor(available / (64 * (as.double(n) + width^2))))
  if (rows < 1)
    .np_plot_joint_unavailable("joint coefficient covariance exceeds the plot memory budget", call. = FALSE)
  rows <- min(rows, n, g)
  operator <- function(x, z) {
    args <- list(bws = bws, txdat = xdat, exdat = x,
                 output = "matrix", iterate = FALSE)
    if (separate.z) {
      args$tzdat <- zdat
      args$ezdat <- z
    }
    if (is.null(gradient.index)) {
      hat <- do.call(npscoefhat, args)
    } else {
      args$output <- NULL
      args$tydat <- ydat
      args$se <- FALSE
      args$.np_coefficient_projection <- gradient.index
      state <- do.call(.np_scoef_fit_internal, args)
      hat <- .npscoef_coefficient_projection_block(state, 1L, nrow(x))
    }
    if (any(!is.finite(hat)))
      .np_plot_joint_unavailable("joint covariance requires a finite smooth-coefficient influence map",
           call. = FALSE)
    hat
  }
  z <- residual
  if (length(z) != n || any(!is.finite(z)))
    .np_plot_joint_unavailable("joint covariance smooth-coefficient residuals are nonfinite", call. = FALSE)
  covariance <- if (factor) NULL else matrix(0, g, g)
  influence <- if (factor) matrix(0, g, n) else NULL
  mean <- numeric(g)
  starts <- seq.int(1L, g, by = rows)
  for (a in seq_along(starts)) {
    ia <- seq.int(starts[a], min(g, starts[a] + rows - 1L))
    hat <- operator(exdat[ia, , drop = FALSE], ezdat[ia, , drop = FALSE])
    mean[ia] <- drop(hat %*% ydat)
    left <- sweep(hat, 2L, z, "*")
    if (factor) {
      influence[ia, ] <- left
      next
    }
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
  list(covariance = covariance, residual = z, block.rows = rows,
       mean = mean, influence = influence)
}

.np_plot_response_residuals <- function(bws, xdat, ydat, rows,
                                          resx = NULL, projection = NULL) {
  n <- nrow(xdat)
  partial <- !is.null(projection)
  K <- if (partial) matrix(0, nrow(projection), n) else NULL
  score <- square <- scale <- numeric(n)
  if (bws[["type", exact = TRUE]] %in% c("fixed", "adaptive_nn")) {
    # These geometries use the same weights at training coordinates whether
    # evaluation is implicit or explicit. Generalized NN does not.
    operator <- function(index) .npreghat_complete(bws, txdat = xdat,
      exdat = xdat[index, , drop = FALSE], output = "matrix")
    progress <- .np_plot_bootstrap_progress_begin(n * (1L + as.integer(partial)),
      "Computing plot training influence", unit = "row")
    on.exit(.np_plot_progress_end(progress), add = TRUE)
    if (partial) {
      K <- projection
      for (start in seq.int(1L, n, by = rows)) {
        index <- seq.int(start, min(n, start + rows - 1L))
        K <- K - projection[, index, drop = FALSE] %*% operator(index)
        progress <- .np_plot_progress_tick(progress, max(index))
      }
    }
    z <- numeric(n)
    for (start in seq.int(1L, n, by = rows)) {
      index <- seq.int(start, min(n, start + rows - 1L))
      E <- -operator(index)
      if (partial) E <- E - resx[index, , drop = FALSE] %*% K
      E[cbind(seq_along(index), index)] <- 1 + E[cbind(seq_along(index), index)]
      scale <- apply(abs(E), 1L, max)
      if (any(!is.finite(scale)) || any(scale == 0))
        .np_plot_joint_unavailable("plot covariance has unidentified training residual information", call. = FALSE)
      E <- E / scale
      z[index] <- drop(E %*% ydat) / sqrt(rowSums(E^2))
      progress <- .np_plot_progress_tick(progress, n * as.integer(partial) + max(index))
    }
    return(list(residual = z, coefficient.map = K))
  }
  spec <- npValidatedConditionalRegSpec(bws, where = "plot response influence", ncon.field = "ncon")
  progress <- .np_plot_bootstrap_progress_begin(n, "Computing plot training influence", unit = "column")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  # Response-column blocks retain implicit training-evaluation radii. External
  # evaluation at the same coordinates is different for generalized NN.
  for (start in seq.int(1L, n, by = rows)) {
    index <- seq.int(start, min(n, start + rows - 1L))
    response <- matrix(0, n, length(index))
    response[cbind(index, seq_along(index))] <- 1
    # Use the unchanged batched owner's declared capability. Other training
    # maps use canonical scalar response fits, keeping workspace bounded even
    # where the generic hat API would materialize a full training matrix.
    if (.npreghat_native_lp_mean_capability(bws, spec$regtype.engine,
        spec$degree.engine, spec$basis.engine, spec$bernstein.basis.engine,
        integer(bws[["ncon", exact = TRUE]]), FALSE)) {
      Sy <- .npreghat_exact_lp_apply_from_regression_core(bws, txdat = xdat,
        y = response, basis = spec$basis.engine, degree = spec$degree.engine,
        bernstein.basis = spec$bernstein.basis.engine)
    } else {
      Sy <- vapply(seq_along(index), function(j)
        .np_regression_direct(bws, txdat = xdat, tydat = response[, j])$mean,
        numeric(n))
    }
    Sy <- matrix(Sy, n, length(index))
    E <- -Sy
    if (partial) {
      K[, index] <- projection[, index, drop = FALSE] - projection %*% Sy
      E <- E - resx %*% K[, index, drop = FALSE]
    }
    E[cbind(index, seq_along(index))] <- 1 + E[cbind(index, seq_along(index))]
    next.scale <- pmax(scale, apply(abs(E), 1L, max))
    divisor <- ifelse(next.scale > 0, next.scale, 1)
    E <- E / divisor
    score <- score * (scale / divisor) + drop(E %*% ydat[index])
    square <- square * (scale / divisor)^2 + rowSums(E^2)
    scale <- next.scale
    progress <- .np_plot_progress_tick(progress, max(index))
  }
  if (any(!is.finite(square)) || any(square == 0))
    .np_plot_joint_unavailable("plot covariance has unidentified training residual information", call. = FALSE)
  list(residual = score / sqrt(square), coefficient.map = K)
}

.np_plot_plreg_influence <- function(bws, xdat, ydat, zdat, exdat, ezdat,
                                     factor = FALSE, block.bytes = 4 * 1024^2,
                                     max.bytes = 64 * 1024^2) {
  n <- nrow(xdat); p <- ncol(xdat); g <- nrow(exdat)
  bandwidths <- bws[["bw", exact = TRUE]]
  ydat <- .np_hat_response(ydat, bws[["ydati", exact = TRUE]])
  if (any(!is.finite(c(block.bytes, max.bytes))) ||
      length(block.bytes) != 1L || length(max.bytes) != 1L ||
      block.bytes <= 0 || max.bytes <= 0)
    stop("partially linear plot memory budgets must be positive and finite", call. = FALSE)
  widths <- vapply(bandwidths, function(bw) {
    spec <- npValidatedConditionalRegSpec(bw, where = "partially linear plot influence", ncon.field = "ncon")
    if (identical(spec$regtype.engine, "lc")) 1 else npLpBasisNcol(spec$basis.engine, spec$degree.engine)
  }, numeric(1L))
  width <- max(widths)
  setup <- 64 * n * (width + p + ncol(zdat) + 1) + 128 * (width^2 + p^2) + 64 * g * (p + 1)
  retained <- if (factor) 16 * as.double(g) * n else 0
  joint.unavailable <- NULL
  if (factor && max.bytes - setup - retained < 64 * as.double(n)) {
    factor <- FALSE
    retained <- 0
    joint.unavailable <- "partially linear joint influence exceeds the plot memory budget"
  }
  rows <- min(n, g, floor(min(block.bytes, max.bytes - setup - retained) / (64 * as.double(n))))
  if (n < 2L || p < 1L || g < 1L || rows < 1 ||
      length(ydat) != n || any(!is.finite(ydat)))
    stop("partially linear plot influence exceeds its data or memory budget", call. = FALSE)
  xnum <- function(dat, j, evaluation = FALSE) {
    if (!is.factor(dat[[j]])) return(as.double(dat[[j]]))
    yi <- bandwidths[[j + 1L]][["ydati", exact = TRUE]]
    adjusted <- .np_regression_response_levels(dat[, j, drop = FALSE], yi, allowNewCells = evaluation)
    yi$all.dlev[[1L]][as.integer(adjusted[[1L]])]
  }
  rx <- matrix(0, n, p); re <- matrix(0, g, p); formation <- numeric(p)
  for (j in seq_len(p)) {
    x <- xnum(xdat, j)
    sx <- .np_regression_direct(bandwidths[[j + 1L]], txdat = zdat, tydat = x)$mean
    lx <- .np_regression_direct(bandwidths[[j + 1L]], txdat = zdat, exdat = ezdat, tydat = x)$mean
    rx[, j] <- x - sx; re[, j] <- xnum(exdat, j, TRUE) - lx
    formation[j] <- .np_plreg_residual_formation_error(x, sx)
  }
  qrR <- qr(rx, tol = .Machine$double.eps)
  .np_plreg_check_residualized_rank(qrR, p, "partially linear plot influence", formation)
  # R^+ without an n-by-n identity or normal-equation inversion.
  B <- matrix(0, p, n)
  B[qrR$pivot, ] <- backsolve(qr.R(qrR), t(qr.Q(qrR)))
  training <- .np_plot_response_residuals(bandwidths[["yzbw", exact = TRUE]],
                                         zdat, ydat, rows, rx, B)
  K <- training$coefficient.map
  z <- training$residual
  operator <- function(evaluation) .npreghat_complete(
    bandwidths[["yzbw", exact = TRUE]], txdat = zdat, exdat = evaluation, output = "matrix")
  influence <- if (factor) matrix(0, g, n) else NULL
  se <- mean <- numeric(g)
  for (start in seq.int(1L, g, by = rows)) {
    index <- seq.int(start, min(g, start + rows - 1L))
    H <- operator(ezdat[index, , drop = FALSE]) + re[index, , drop = FALSE] %*% K
    mean[index] <- drop(H %*% ydat)
    A <- sweep(H, 2L, z, "*")
    se[index] <- sqrt(rowSums(A^2))
    if (factor) influence[index, ] <- A
  }
  list(mean = mean, se = se, influence = influence, residual = z, block.rows = rows,
       joint.unavailable = joint.unavailable)
}
