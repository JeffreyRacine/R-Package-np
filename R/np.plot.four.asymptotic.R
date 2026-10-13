# Four-family Gaussian plot adapters. Existing estimators, marginal uncertainty,
# regression influence owners and Gaussian calibrators remain unchanged.
.np_plot_four_calibrate <- function(A, se, alpha, draws, target,
                                      exact.diagonal = FALSE) {
  if (!is.matrix(A) || any(!is.finite(A)))
    stop("invalid four-family influence factor", call. = FALSE)
  if (length(se) != nrow(A) || any(!is.finite(se)) || any(se < 0))
    .np_plot_joint_unavailable("marginal uncertainty is unavailable on the declared grid")
  norm <- numeric(nrow(A))
  zero <- logical(nrow(A))
  correlation.factor <- A
  for (i in seq_len(nrow(A))) {
    # A squared norm can underflow even when the influence is nonzero. Only
    # exact zero entries certify a degenerate row. Standardize before squaring.
    zero[i] <- all(A[i, ] == 0)
    if (zero[i]) next
    scale <- max(abs(A[i, ]))
    row <- A[i, ] / scale
    size <- sqrt(sum(row^2))
    norm[i] <- scale * size
    if (!is.finite(norm[i]) || norm[i] <= 0)
      stop("influence variance is not representable", call. = FALSE)
    correlation.factor[i, ] <- row / size
  }
  if (any(zero != (se == 0)))
    stop("influence and marginal uncertainty disagree about zero variance", call. = FALSE)
  if (exact.diagonal && any(abs(norm[!zero] / se[!zero] - 1) >
      256 * .Machine$double.eps))
    stop("joint influence does not reproduce the marginal variance", call. = FALSE)
  # This unit-row factor represents correlation only. Displayed widths retain
  # the separately justified marginal SE approximation. The unscaled factor
  # and its independent covariance oracle are never changed to fit a diagonal.
  cal <- if (nrow(A) <= 512L)
    .np_plot_gaussian_calibrate(tcrossprod(correlation.factor), alpha = alpha, draws = draws) else
    .np_plot_gaussian_factor_calibrate(correlation.factor, alpha = alpha, draws = draws)
  cal$target <- target
  cal$marginal.scale <- "incumbent asymptotic standard errors"
  cal$dependence <- "observation influence correlation"
  cal
}

.np_plot_four_interval <- function(interval, calibration, se, band.type) {
  if (is.null(calibration)) return(interval)
  simultaneous <- cbind(se, se, deparse.level = 0L) * calibration$critical
  if (identical(band.type, "all"))
    interval$all.err$simultaneous <- simultaneous else
    interval$err[, 1:2] <- simultaneous
  interval$calibration <- calibration
  interval
}

.np_plot_lsq_joint <- function(bws, x, ex, tau, se, gradient = NULL,
                                gradient.order = 1L, alpha = .05,
                                gaussian.draws = 9999L) {
  index <- .nplsqreg_match_plot_tau(tau, bws[["tau", exact = TRUE]])
  one <- if (length(bws[["tau", exact = TRUE]]) == 1L) bws else
    bws[["tau.bws", exact = TRUE]][[index]]
  if (is.null(one)) stop("missing per-tau location-scale state", call. = FALSE)
  reg <- one[["reg.bws", exact = TRUE]]
  if (!.np_plot_regression_joint_supported(reg))
    .np_plot_joint_unavailable("this location-scale regression kernel is not qualified")
  deriv <- categorical <- NULL
  if (!is.null(gradient)) {
    continuous <- which(reg[["icon", exact = TRUE]])
    if (gradient %in% continuous) {
      deriv <- integer(length(continuous))
      deriv[match(gradient, continuous)] <- rep_len(gradient.order,
        length(continuous))[match(gradient, continuous)]
    } else categorical <- gradient
  }
  joint <- .np_plot_regression_joint_covariance(reg, x,
    one[["qdat", exact = TRUE]], ex, deriv = deriv,
    categorical.index = categorical, factor = TRUE)
  .np_plot_four_calibrate(joint$influence, se, alpha, gaussian.draws,
    "location-scale transformed-response regression; frozen pilot, delta and bandwidths")
}

.np_plot_quantile_interval <- function(interval, bws, x, y, ex, tau,
                                       estimate, se, gradient = NULL,
                                       gradient.order = 1L, alpha = .05,
                                       band.type = "all", gaussian.draws = 9999L) {
  calibration <- .np_plot_joint_try(if (inherits(bws, "lsqregressionbandwidth"))
    .np_plot_lsq_joint(bws, x, ex, tau, se, gradient, gradient.order,
      alpha, gaussian.draws) else
    .np_plot_qreg_joint(bws, x, y, ex, tau, estimate, se, gradient,
      alpha, gaussian.draws))
  .np_plot_four_interval(interval, calibration, se, band.type)
}

.np_plot_quantile_band_fields <- function(object, interval, gradient = NULL,
                                          tau.index = NULL) {
  if (is.null(tau.index))
    return(.np_plot_joint_band_fields(object, interval$calibration,
      interval$all.err, gradient))
  key <- if (is.null(gradient)) "quantile" else paste0("gradient", gradient)
  if (!is.null(interval$calibration))
    object$band.calibration[[key]][[tau.index]] <- interval$calibration
  if (!is.null(interval$all.err))
    object$band.errors[[key]][[tau.index]] <- interval$all.err
  object
}

.np_plot_four_budget <- function(n, m, p, width = 1L,
                                  max.bytes = 128 * 1024^2) {
  if (n < 2L || m < 1L || !is.finite(width) || width < 1L)
    stop("invalid four-family influence dimensions", call. = FALSE)
  tile <- min(m, max(1L, floor(4 * 1024^2 / (8 * as.double(n) * (12 + p)))))
  forecast <- 24 * as.double(n) * m + 64 * 1024^2 +
    8 * as.double(n) * tile * (12 + p) +
    64 * (as.double(n) + tile) * width + 256 * width^2
  if (!is.finite(forecast) || forecast > max.bytes)
    .np_plot_joint_unavailable("four-family influence exceeds the plot memory budget")
  list(tile = tile, forecast.bytes = forecast)
}

.np_plot_four_conditional_admit <- function(bws, family) {
  if (!identical(bws[["type", exact = TRUE]], "fixed"))
    .np_plot_joint_unavailable(paste(family, "joint bands currently require fixed bandwidths"))
  if (!.np_plot_density_joint_supported(bws, TRUE) ||
      identical(bws[["cxkertype", exact = TRUE]], "beta") ||
      identical(bws[["cykertype", exact = TRUE]], "beta"))
    .np_plot_joint_unavailable(paste(family, "joint bands require qualified second-order Gaussian/Epanechnikov kernels"))
  spec <- npConditionalRegEngineSpec(bws, where = "four-family joint plot")
  if (all(spec$degree.engine == 0L)) 1L else
    npLpBasisNcol(spec$basis.engine, spec$degree.engine)
}

.np_plot_qreg_influence <- function(bws, x, y, ex, tau, quantile, se) {
  width <- .np_plot_four_conditional_admit(bws, "quantile")
  n <- nrow(x); m <- nrow(ex)
  budget <- .np_plot_four_budget(n, m, ncol(x), width)
  if (length(quantile) != m || any(!is.finite(quantile)) ||
      length(se) != m || any(!is.finite(se)) || any(se < 0))
    .np_plot_joint_unavailable("quantiles or marginal uncertainty are unavailable on the declared grid")
  if (ncol(y) != 1L || !is.numeric(y[[1L]]))
    .np_plot_joint_unavailable("regular quantile joint bands require one continuous response")
  # The inverter clamps to observed support when the requested probability
  # lies outside its bracket. A small residual alone does not make that root
  # regular. Metadata can be lost during plot row selection, so retain the
  # conservative endpoint certificate as well; never alter the fitted root.
  support <- range(y[[1L]], finite = TRUE)
  if (any(.npqreg_quantile_clamp(quantile) != "none") ||
      any(quantile <= support[1L] | quantile >= support[2L]))
    .np_plot_joint_unavailable("quantile joint bands require interior, unclamped roots")
  A <- matrix(NA_real_, m, n)
  density <- root <- numeric(m)
  progress <- .np_plot_bootstrap_progress_begin(m, "Computing quantile plot influence", unit = "row")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (first in seq.int(1L, m, by = budget$tile)) {
    rows <- seq.int(first, min(m, first + budget$tile - 1L))
    ey <- setNames(data.frame(quantile[rows]), names(y))
    value <- .np_plot_density_conditional(bws, x, y, ex[rows, , drop = FALSE], ey, cdf = TRUE)
    dens <- .np_plot_conditional_eval(bws = bws, xdat = x, ydat = y,
      exdat = ex[rows, , drop = FALSE], eydat = ey, cdf = FALSE,
      gradients = FALSE, se = FALSE, allow.external = TRUE)[["condens", exact = TRUE]]
    if (any(!is.finite(dens)) || any(dens <= 0))
      .np_plot_joint_unavailable("quantile joint bands require finite positive density slopes")
    residual <- value$estimate - tau
    # Inversion error must be negligible on the displayed uncertainty scale.
    # Do not refit, floor a slope, or accept a clamped tail root as regular.
    tolerance <- .01 * se[rows] * dens +
      max(128, 8 * as.double(n)) * .Machine$double.eps
    if (any(abs(residual) > tolerance))
      .np_plot_joint_unavailable("quantile inversion residual exceeds 1% of its asymptotic uncertainty")
    A[rows, ] <- -value$A / dens
    density[rows] <- dens; root[rows] <- residual
    progress <- .np_plot_progress_tick(progress, max(rows))
  }
  list(A = A, density = density, root.residual = root,
    forecast.bytes = budget$forecast.bytes)
}

.np_plot_qreg_joint <- function(bws, x, y, ex, tau, estimate, se,
                                 gradient = NULL, alpha = .05,
                                 gaussian.draws = 9999L) {
  if (!is.null(gradient))
    .np_plot_joint_unavailable("quantile derivative joint bands require a separately qualified joint delta approximation")
  value <- .np_plot_qreg_influence(bws, x, y, ex, tau, estimate, se)
  cal <- .np_plot_four_calibrate(value$A, se, alpha, gaussian.draws,
    "inverse conditional CDF quantiles; fixed bandwidths, within one tau curve")
  cal$maximum.root.residual <- max(abs(value$root.residual))
  cal$forecast.bytes <- value$forecast.bytes
  cal
}

.np_plot_conmode_influence <- function(object, x, y, ex, level, levels,
                                        estimate, se) {
  bws <- object[["bws", exact = TRUE]]
  width <- .np_plot_four_conditional_admit(bws, "class-probability")
  n <- nrow(x); m <- nrow(ex)
  budget <- .np_plot_four_budget(n, m, ncol(x), width)
  y.frame <- toFrame(y)
  if (length(bws[["ynames", exact = TRUE]]) == ncol(y.frame))
    names(y.frame) <- bws[["ynames", exact = TRUE]]
  if (length(se) != m || any(!is.finite(se)) || any(se < 0))
    .np_plot_joint_unavailable("class-probability uncertainty is unavailable on the declared grid")
  chosen <- match(level, levels)
  if (is.na(chosen)) stop("unknown class-probability influence level", call. = FALSE)
  A <- matrix(NA_real_, m, n); fitted <- numeric(m)
  proper <- isTRUE(object[["proper.requested", exact = TRUE]])
  tol <- object[["proper.info", exact = TRUE]][["tol", exact = TRUE]]
  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)
  progress <- .np_plot_bootstrap_progress_begin(m, "Computing class-probability plot influence", unit = "row")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (first in seq.int(1L, m, by = budget$tile)) {
    rows <- seq.int(first, min(m, first + budget$tile - 1L))
    total <- numeric(length(rows)); total.A <- matrix(0, length(rows), n)
    selected <- NULL; minimum <- rep(Inf, length(rows))
    # All levels use the same observations; stream their factors. No
    # level-by-grid-by-observation tensor or independent class draws.
    for (k in seq_along(levels)) {
      ey <- .np_plot_conmode_level_factor(y, levels[k], length(rows))
      value <- .np_plot_density_conditional(bws, x, y.frame, ex[rows, , drop = FALSE], toFrame(ey))
      total <- total + value$estimate; total.A <- total.A + value$A
      minimum <- pmin(minimum, value$estimate)
      if (k == chosen) selected <- value
    }
    if (proper) {
      if (any(!is.finite(total)) || any(total <= 0) || any(minimum < 0) ||
          any(abs(total - 1) > tol))
        .np_plot_joint_unavailable("projected or clipped class probabilities do not have a qualified smooth joint approximation")
      selected$estimate <- selected$estimate / total
      selected$A <- (selected$A - selected$estimate * total.A) / total
    }
    A[rows, ] <- selected$A; fitted[rows] <- selected$estimate
    progress <- .np_plot_progress_tick(progress, max(rows))
  }
  roundoff <- max(128, 8 * as.double(n)) * .Machine$double.eps *
    pmax(1, abs(fitted), abs(estimate))
  if (length(estimate) != m || any(!is.finite(fitted)) ||
      any(abs(fitted - estimate) > .01 * se + roundoff))
    stop("class-probability influence does not reproduce the plotted target", call. = FALSE)
  list(A = A, estimate = fitted, forecast.bytes = budget$forecast.bytes)
}

.np_plot_conmode_joint_columns <- function(dat, object, x, y, ex, level,
                                            levels, se, alpha, band,
                                            gaussian.draws) {
  cal <- .np_plot_joint_try({
    value <- .np_plot_conmode_influence(object, x, y, ex, level, levels,
      dat$probability, se)
    result <- .np_plot_four_calibrate(value$A, se, alpha, gaussian.draws,
      "named class probability; fixed bandwidths and regular proper normalization")
    result$forecast.bytes <- value$forecast.bytes
    result
  })
  if (is.null(cal)) return(dat)
  lower <- dat$center - cal$critical * se
  upper <- dat$center + cal$critical * se
  if (identical(band, "all")) {
    dat$simultaneous.lower <- lower; dat$simultaneous.upper <- upper
  } else {
    dat$lower <- lower; dat$upper <- upper
  }
  attr(dat, "band.calibration") <- cal
  dat
}

.np_plot_copula_influence <- function(object, data, ex, se) {
  bws <- object[["bws", exact = TRUE]]
  if (!identical(bws[["type", exact = TRUE]], "fixed"))
    .np_plot_joint_unavailable("copula joint bands currently require fixed bandwidths")
  n <- nrow(data); m <- nrow(ex); p <- ncol(data)
  budget <- .np_plot_four_budget(n, m, p)
  density <- identical(object[["target", exact = TRUE]], "density")
  if (!density && !.np_plot_density_joint_supported(bws, FALSE))
    .np_plot_joint_unavailable("this copula distribution kernel is not qualified for joint bands")
  if (length(se) != m || any(!is.finite(se)) || any(se < 0))
    .np_plot_joint_unavailable("copula uncertainty is unavailable on the declared grid")
  A <- matrix(NA_real_, m, n); fitted <- numeric(m)
  if (density) {
    marginal <- lapply(seq_len(p), function(j)
      .npcopula_marginal_bw(bws, data, j, target = "density", kbandwidth = TRUE))
    train <- lapply(seq_len(p), function(j) .npcopula_marginal_data(bws, data, j))
    evaluate <- lapply(seq_len(p), function(j) .npcopula_marginal_eval_data(bws, ex, j))
  }
  progress <- .np_plot_bootstrap_progress_begin(m, "Computing copula plot influence", unit = "row")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (first in seq.int(1L, m, by = budget$tile)) {
    rows <- seq.int(first, min(m, first + budget$tile - 1L))
    if (density) {
      product <- matrix(1, length(rows), n); total <- matrix(0, length(rows), n)
      for (j in seq_len(p)) {
        W <- .np_plot_density_matrix(marginal[[j]], train[[j]],
          evaluate[[j]][rows, , drop = FALSE], "normal")
        mass <- rowMeans(W)
        if (any(!is.finite(mass)) || any(mass == 0))
          .np_plot_joint_unavailable("copula density has a zero or nonfinite marginal denominator")
        W <- W / mass
        product <- product * W; total <- total + W
      }
      fitted[rows] <- rowMeans(product)
      phi <- product - fitted[rows] * (total - p + 1)
      A[rows, ] <- (phi - rowMeans(phi)) / n
    } else {
      value <- .np_plot_density_unconditional(bws, data, ex[rows, , drop = FALSE], cdf = TRUE)
      A[rows, ] <- value$A; fitted[rows] <- value$estimate
    }
    progress <- .np_plot_progress_tick(progress, max(rows))
  }
  estimate <- object[["copula", exact = TRUE]]
  roundoff <- max(128, 8 * as.double(n)) * .Machine$double.eps *
    pmax(1, abs(fitted), abs(estimate))
  if (length(estimate) != m || any(!is.finite(fitted)) ||
      any(abs(fitted - estimate) > .01 * se + roundoff))
    stop("copula influence does not reproduce the plotted target", call. = FALSE)
  list(A = A, estimate = fitted, forecast.bytes = budget$forecast.bytes,
    exact.diagonal = density)
}

.np_plot_copula_interval <- function(interval, object, data, ex, se,
                                      alpha, band.type, gaussian.draws) {
  cal <- .np_plot_joint_try({
    value <- .np_plot_copula_influence(object, data, ex, se)
    result <- .np_plot_four_calibrate(value$A, se, alpha, gaussian.draws,
      if (value$exact.diagonal) "copula density ratio; frozen marginal inverse coordinates and bandwidths" else
        "copula distribution; frozen marginal inverse coordinates and bandwidths",
      exact.diagonal = value$exact.diagonal)
    result$forecast.bytes <- value$forecast.bytes
    result
  })
  .np_plot_four_interval(interval, cal, se, band.type)
}
