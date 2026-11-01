# Public adapter: reuse the already captured fit, never refit selected geometry.
.np_reg_percentile_admit <- function(bws, control, alpha, B, bootstrap, blocklen) {
  spec <- npValidatedConditionalRegSpec(bws, where = 'regression percentile', ncon.field = 'ncon')
  if (!identical(bootstrap, 'residual') || !is.null(blocklen))
    stop('regression smoothed-percentile uses bootstrap="residual"; omit bootstrap or select residual explicitly, without blocklen', call. = FALSE)
  if (!identical(bws[['type', exact = TRUE]], 'fixed') ||
      !all(bws[['icon', exact = TRUE]]) || !length(bws$icon) %in% 1:2 ||
      bws$ckertype != 'gaussian' || bws$ckerorder != 2L || bws$ckerbound != 'none' ||
      !bws$regtype %in% c('lc', 'll', 'lp') ||
      (bws$regtype == 'lp' && (any(spec$degree.engine > 3L) || !isTRUE(bws$bernstein.basis))))
    stop('regression smoothed-percentile currently requires fixed Gaussian2 numeric q=1 or q=2, LC/LL or Bernstein LP degree at most 3', call. = FALSE)
  if (!identical(control$adaptation, 'none') || control$iterations != 16L ||
      !isTRUE(alpha == .05))
    stop('regression smoothed-percentile currently requires alpha=0.05, 16 iterations and bias.control=np_bias_control(adaptation="none")', call. = FALSE)
  .np_plot_gaussian_integer(B, 'regression percentile draws', 9999L)
  invisible(NULL)
}
.np_reg_percentile_prepare <- function(state, fit, grid, B = NULL,
                                        targets = c('mean', paste0('d', seq_len(state$d)))) {
  valid <- c('mean', paste0('d', seq_len(state$d)))
  if (!is.character(targets) || !length(targets) || anyNA(targets) ||
      anyDuplicated(targets) || any(!targets %in% valid))
    stop('invalid regression percentile target family', call.=FALSE)
  full.map <- expand.grid(row=seq_len(nrow(grid)), target=valid, stringsAsFactors=FALSE)
  map <- full.map[unlist(lapply(targets, function(t) which(full.map$target == t)),
                        use.names=FALSE), , drop=FALSE]
  n <- length(state$y); m <- nrow(map)
  if (!is.null(B)) .np_reg_percentile_forecast(n, m, B)
  if (8 * (6 * as.double(m) * n + 10 * n * 39 + 20 * n * 256) > 256 * 1024^2)
    stop('regression percentile preparation exceeds its 256 MiB memory budget', call. = FALSE)
  # Preserve the qualified bounded arithmetic tile shape. Reducing its RHS
  # width changes BLAS rounding, even for rows retained by both preparations.
  capture <- .np_bias_capture_prepare(state, fit, grid, valid)
  A <- matrix(0, m, n); center <- se <- numeric(m)
  key <- paste(map$target, map$row)
  for (batch in .np_bias_target_blocks(nrow(grid), length(valid))) {
    jj <- batch$grid; lm <- full.map[batch$rows, , drop=FALSE]
    at <- match(paste(lm$target, lm$row), key)
    lm$row <- match(lm$row, jj); v <- capture$rows(jj, lm)
    keep <- which(!is.na(at)); rr <- at[keep]
    A[rr, ] <- v$W[keep,,drop=FALSE]
    center[rr] <- v$center[keep]; se[rr] <- v$se[keep]
  }
  sp <- .np_bias_student_prepare(state); pilot <- numeric(m); eg <- grid
  for (k in seq_len(state$d)) eg[[k]] <- (eg[[k]] - state$offset[k]) / state$scale[k]
  for (j in seq_len(nrow(grid))) {
    v <- .np_bias_pilot_rows(state, eg[j, , drop = FALSE], derivatives = TRUE)
    for (t in seq_along(targets)) {
      k <- if (targets[t] == 'mean') 0L else as.integer(substring(targets[t], 2L))
      pilot[j + (t - 1L) * nrow(grid)] <- if (k == 0L)
        drop(v$mean %*% state$y[v$indices]) else
        drop(v$derivatives[[k]] %*% state$y[v$indices]) / state$scale[k]
    }
  }
  list(state = state, A = A, center = center, se = se, y = state$y,
    map = map, mu = state$y - state$residual_raw,
    sigma = .np_bias_student_error_scale(sp, 'mammen'), pilot.eval = pilot,
    residual = state$residual_raw, row.variance = state$residual_row_variance)
}
.np_reg_percentile_plot <- function(state, fit, grid, target, B, bands, resolution.context = NULL) {
  names(grid) <- names(state$x)
  p <- .np_reg_percentile_prepare(state, fit, grid, B = B, targets = target)
  start <- get('.Random.seed', .GlobalEnv, inherits = FALSE)
  seeds <- .np_bias_stream(start, sample.int(2147483646L, 79L, replace = FALSE))$value
  z <- .np_reg_percentile_infer(p, B, seeds[1L], seeds[2:40], seeds[41:79])
  ii <- which(p$map$target == target)
  if (!length(ii)) stop('regression percentile supports means and first derivatives only', call. = FALSE)
  z$bands <- lapply(z$bands[bands], function(x) x[ii, , , drop = FALSE])
  z$center <- z$center[ii]; z$se <- z$se[ii]; z$expectation <- z$expectation[ii]
  z$variance <- z$variance[ii]; z$map <- z$map[ii, , drop = FALSE]
  z$ranks <- z$ranks[target]; z$calibration <- z$calibration[target]
  # Retain the original full-family diagnostic row identities after reducing
  # persistent storage to the requested family. Arithmetic and RNG are unchanged.
  families <- sort(c('mean', paste0('d', seq_len(state$d))))
  rownames(z$resolution) <- (match(target, families) - 1L) * nrow(z$resolution) +
    seq_len(nrow(z$resolution))
  z$resolution <- z$resolution[z$resolution$family == target & z$resolution$band %in% bands, , drop = FALSE]
  z$settings <- list(construction = 'regression-calibrated-residual-percentile-E1-E6',
    calibration.requested = 'smoothed-percentile', law = 'standardized-empirical-residual',
    residual.smoothing = FALSE, B = B, alpha = .05, outer = 39L, rank = 38L,
    iterations = 16L, geometry = 'fixed selected bandwidths and degrees',
    endpoint.construction = 'unshifted corrected curves; E1 mean / E6 derivative companion ranks',
    coverage.status = 'experimental; q1 evidence and conservative q2 cells; sensitivity and control-loss failures retained')
  z$warnings <- character()
  if (isTRUE(z$calibration[[target]]$support.used))
    z$warnings <- c(z$warnings, 'Full bootstrap-support endpoint used; bounds may be very wide')
  if (!isTRUE(z$calibration[[target]]$attainable))
    z$warnings <- c(z$warnings, 'Calibration remains unattainable; all bands retained as diagnostic only')
  z$calibration.seeds <- list(original = seeds[1L], outer.data = seeds[2:40], outer.inner = seeds[41:79])
  if (is.null(resolution.context)) {
    messages <- z$warnings
    withCallingHandlers(.np_band_resolution_warn(z$resolution), warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart('muffleWarning')
    })
    if (length(messages)) warning(paste(unique(messages), collapse = '; '), call. = FALSE)
  } else resolution.context$percentile.messages <- unique(c(resolution.context$percentile.messages, z$warnings))
  z
}
