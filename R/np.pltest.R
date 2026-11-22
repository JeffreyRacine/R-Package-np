# Statistic: Sun, Chen, Liang and Ruppert (2021), doi:10.5705/ss.202019.0238.
.pltest_geometry <- function(predictors, standardize) {
  columns <- lapply(predictors, function(x) {
    if (is.factor(x) && !is.ordered(x)) {
      code <- as.integer(droplevels(x))
      return(outer(code, seq_len(max(code)), "==") / sqrt(2))
    }
    v <- as.numeric(x)
    if (standardize) {
      magnitude <- max(abs(v))
      if (magnitude > 0) v <- v/magnitude
      v <- v-mean(v)
      scale <- sd(v)
      if (scale > 0) v <- v/scale
    }
    v
  })
  V <- do.call(cbind, unname(columns))
  if (any(!is.finite(V))) stop("Non-finite projection coordinates; rescale the predictors.", call. = FALSE)
  storage.mode(V) <- "double"
  .np_with_compiled_fit_progress("Computing projection weights", nrow(V),
    expr = .Call(C_np_pltest_angles, V))
}

.pltest_sun_bandwidth <- function(X, y, Z, continuous, cells) {
  .pltest_local({
    q <- sum(continuous)
    if (q >= 4L)
      stop("Sun local-linear undersmoothing requires fewer than four continuous nuisance coordinates. Use bw.method = 'cv' or 'object'.", call. = FALSE)
    h <- numeric(ncol(Z))
    if (q) {
      C <- as.matrix(Z[, continuous, drop = FALSE])
      scales <- apply(C, 2L, function(v) min(sd(v), 0.75 * IQR(v)))
      if (any(!is.finite(scales) | scales <= 0)) stop("Each continuous nuisance coordinate needs a finite positive Sun scale (SD and IQR).", call. = FALSE)
      h[continuous] <- 1.06 * scales * length(y)^(-3 / (q + 8))
      if (any(!is.finite(h[continuous]) | h[continuous] <= 0))
        stop("The Sun bandwidth is not representable in the nuisance units. Rescale Z.", call. = FALSE)
      for (idx in split(seq_along(y), cells)) {
        # The intercept makes each local rank invariant to the origin of Z.
        cell <- sweep(C[idx, , drop = FALSE], 2L, C[idx[1L], ], "-")
        cell <- sweep(cell, 2L, h[continuous], "/")
        if (any(!is.finite(cell)))
          stop("The standardized local differences are not finite. Rescale Z.", call. = FALSE)
        design <- cbind(1, cell)
        if (length(idx) < q + 1L || qr(design)$rank < q + 1L)
          stop("The continuous local-linear design is rank deficient in factor cell ", cells[idx[1L]], ".", call. = FALSE)
        for (j in seq_along(idx)) {
          if (j %% 64L == 1L) .np_fit_progress_step()
          delta <- sweep(cell, 2L, cell[j, ], "-")
          weighted <- cbind(1, delta) * exp(-rowSums(delta^2) / 4)
          if (qr(weighted)$rank < q + 1L)
            stop("The Gaussian-weighted local-linear design is singular at row ", idx[j], ".", call. = FALSE)
        }
      }
    }
  })
  .pltest_with_local_fit(npplregbw(xdat = X, ydat = y, zdat = Z, bwtype = "fixed",
                bws = matrix(h, nrow = ncol(X) + 1L, ncol = ncol(Z), byrow = TRUE),
                regtype = if (q) "ll" else "lc", bandwidth.compute = FALSE))
}

# Restricted pilot for the bootstrap only. Linear factors use np's slope coding.
.pltest_pilot <- function(X, y, Z, continuous, cells, bw) {
  linear <- do.call(cbind, lapply(seq_along(X), function(j) {
    .np_plreg_numeric_response(X[[j]], bw[["bw"]][[j + 1L]][["ydati"]])
  }))
  groups <- outer(as.integer(cells), seq_len(nlevels(cells)), "==") * 1
  base <- cbind(groups, linear)
  n <- length(y)
  candidate <- function(design, degree, shape) {
    if (ncol(design) >= n) return(NULL)
    fit <- qr(design)
    if (fit$rank != ncol(design)) return(NULL)
    leverage <- rowSums(qr.Q(fit)^2)
    if (any(1 - leverage <= sqrt(.Machine$double.eps))) return(NULL)
    residual <- qr.resid(fit, y)
    list(center = y - residual, amplitude = residual / sqrt(1 - leverage),
         degree = degree, shape = shape, rank = fit$rank, leverage = leverage,
         cv = mean((residual / (1 - leverage))^2))
  }
  best <- candidate(base, 0L, "cell intercepts")
  C <- as.matrix(Z[, continuous, drop = FALSE])
  q <- ncol(C)
  if (q > 0L) {
    # Affine scaling preserves each polynomial space and makes QR more stable.
    C <- scale(C, scale = FALSE)
    span <- apply(abs(C), 2L, max)
    C <- sweep(C, 2L, ifelse(span > 0, span, 1), "/")
    powers <- list()
    for (degree in 1:3) {
      if (ncol(base) + choose(q + degree, degree) - 1 >= n) break
      # Nondecreasing coordinate indices enumerate each monomial exactly once.
      terms <- if (degree == 1L) as.list(seq_len(q)) else
        do.call(c, lapply(terms, function(t) lapply(seq.int(t[length(t)], q),
                                                    function(j) c(t, j))))
      powers <- c(powers, lapply(terms, function(t) Reduce(`*`, lapply(t, function(j) C[, j]))))
      block <- do.call(cbind, powers)
      for (shape in if (ncol(groups) > 1L) c("common", "cell-specific") else "common") {
        if (shape == "cell-specific" && ncol(base) + ncol(groups)*ncol(block) >= n) next
        extra <- if (shape == "common") block else
          do.call(cbind, lapply(seq_len(ncol(groups)), function(j) block * groups[, j]))
        fit <- candidate(cbind(base, extra), degree, shape)
        if (!is.null(fit) && (is.null(best) || fit$cv < best$cv)) best <- fit
      }
    }
  }
  if (is.null(best)) stop("No identifiable restricted pilot has leave-one-out residuals at every row.", call. = FALSE)
  best
}

.pltest_with_local_fit <- function(expr) force(expr)

# Reuse the regression fit's accepted categorical gates when constructing rows.
# The ordinary hat and the fit can use different weights near a factor endpoint.
.pltest_mean <- function(bw, Z, y) {
  categorical <- .np_regression_cat_profile_mean(bw, Z, y)
  if (!is.null(categorical)) return(as.numeric(categorical))
  as.numeric(.np_regression_direct(bw, txdat = Z, tydat = y, local.mode = TRUE)$mean)
}

.pltest_response_operator <- function(bw, Z, y, shared.sun) {
  n <- length(y)
  if (shared.sun && all(vapply(Z, is.numeric, logical(1))))
    return(.npreghat_complete(bws = bw, txdat = Z, output = "matrix"))

  if (!bw$ncon && npUseCategoricalCompress(ncon = 0L, ncat = ncol(Z)) &&
      identical(bw$type, "fixed")) {
    z <- adjustLevels(Z, bw$xdati)
    codes <- .np_cat_profile_code_matrix(z)
    keys <- .np_cat_profile_keys(codes)
    first <- which(!duplicated(keys))
    id <- match(keys, keys[first])
    profiles <- codes[first, , drop = FALSE]
    K <- .np_regression_cat_profile_kernel_matrix(profiles, profiles,
                                                    z[first, , drop = FALSE], bw)
    den <- drop(K %*% tabulate(id, nbins = length(first)))
    if (all(is.finite(den)) && all(abs(den) > .Machine$double.xmin))
      return((K / den)[id, id, drop = FALSE])
  }

  # For non-Gaussian or NN fits, retain the training solve's projection and
  # replay complete-design kernel gates. Slicing the training sample would
  # change both near-constant decisions and implicit NN occurrence radii.
  rows <- .np_plot_lsq_gate_level_operator(bw, Z, y, max.bytes = .Machine$double.xmax)
  S <- matrix(NA_real_, n, n)
  width <- .np_plot_plreg_apply_block_rows(n, n)
  for (first in seq.int(1L, n, by = width)) {
    ii <- seq.int(first, min(n, first + width - 1L))
    S[ii, ] <- rows(ii)
    .np_fit_progress_step()
  }
  S
}

.pltest_training_solve <- function(children, X, y, Z, shared = NULL) {
  n <- length(y)
  linear <- vapply(seq_along(X), function(j)
    .np_plreg_numeric_response(X[[j]], children[[j + 1L]][["ydati"]]), numeric(n))
  smooth <- if (!is.null(shared)) shared %*% linear else
    vapply(seq_along(X), function(j)
      .pltest_mean(children[[j + 1L]], Z, X[[j]]), numeric(n))
  formation <- vapply(seq_along(X), function(j)
    .np_plreg_residual_formation_error(linear[, j], smooth[, j]), numeric(1))
  yhat <- .pltest_mean(children[[1L]], Z, y)
  .np_plreg_linear_solve(y - yhat, linear - smooth, y, yhat, ncol(Z),
                         where = "nppltest", se = FALSE, formation.error = formation)
}

.pltest_operator <- function(bw, X, y, Z, shared.sun = FALSE) {
  .pltest_local(.pltest_with_local_fit({
    children <- bw[["bw"]]
    S <- .pltest_response_operator(children[[1L]], Z, y, shared.sun)
    solved <- .pltest_training_solve(children, X, y, Z, if (shared.sun) S)
    Q <- qr.Q(solved$qr)
    H <- S + Q %*% (t(Q) - crossprod(Q, S))
    mu <- solved$train.fit
    replay <- drop(H %*% y)
    if (any(!is.finite(H)) || any(!is.finite(mu)) || any(!is.finite(replay)) ||
        max(abs(replay - mu)) > 1e-7 * max(abs(y), .Machine$double.xmin))
      stop("The selected bandwidth object did not reproduce its public fixed-bandwidth fit.",
           call. = FALSE)
    list(H = H, mean = mu)
  }))
}

.pltest_residual_norms <- function(H) {
  n <- nrow(H)
  width <- .np_plot_plreg_apply_block_rows(n, n)
  norms <- numeric(n)
  for (first in seq.int(1L, n, by = width)) {
    rows <- seq.int(first, min(n, first + width - 1L))
    residual <- -H[rows, , drop = FALSE]
    diagonal <- cbind(seq_along(rows), rows)
    residual[diagonal] <- 1 + residual[diagonal]
    norms[rows] <- rowSums(residual^2)
  }
  norms
}

.pltest_cv_pilot <- function(X, y, Z, adjust) {
  bw <- npplregbw(xdat = X, ydat = y, zdat = Z, bwmethod = "cv.ls",
                  regtype = if (any(!vapply(Z, is.factor, logical(1)))) "ll" else "lc",
                  bwtype = "fixed", nmulti = 1L)
  if (!adjust) {
    center <- .pltest_with_local_fit(.pltest_training_solve(bw[["bw"]], X, y, Z)$train.fit)
    if (length(center) != length(y) || any(!is.finite(center)))
      stop("The prediction-CV pilot returned an invalid training fit.", call. = FALSE)
    return(list(center = center, amplitude = y - center, bws = bw))
  }
  pilot <- .pltest_operator(bw, X, y, Z)
  .pltest_local({
    variance.factor <- .pltest_residual_norms(pilot$H)
    if (any(!is.finite(variance.factor) | variance.factor <= .Machine$double.eps))
      stop("CV-pilot HC2 is undefined because a residual row has zero or non-finite squared norm.",
           call. = FALSE)
  })
  list(center = pilot$mean, amplitude = (y - pilot$mean)/sqrt(variance.factor),
       leverage = diag(pilot$H), variance.factor = variance.factor, bws = bw)
}

# Draw order is observation within replication, irrespective of block ownership.
.pltest_draw_block <- function(size, H, A, center, amplitude, boot.method) {
  u <- matrix(stats::runif(as.double(length(center))*size), nrow = length(center))
  .pltest_block_statistics(u, H, A, center, amplitude, boot.method)
}

.pltest_block_statistics <- function(u, H, A, center, amplitude, boot.method) {
  xi <- if (boot.method == "wild-rademacher") ifelse(u < 0.5, -1, 1) else
    ifelse(u < (sqrt(5)+1)/(2*sqrt(5)), (1-sqrt(5))/2, (1+sqrt(5))/2)
  ystar <- center + amplitude*xi
  e <- ystar - H %*% ystar
  colSums(e*(A %*% e))
}

nppltest <- function(object, bw.method = c("sun", "cv", "object"),
                     boot.method = c("wild-rademacher", "wild"), B = 399L,
                     random.seed = 42L, projection.scale = c("original", "standardized"), ...,
                     boot.correction = c("none", "cv-pilot", "pilot-hc2", "cv-raw")) {
  started <- proc.time()[["elapsed"]]
  call <- match.call()
  dot.names <- names(match.call(expand.dots = FALSE)$...)
  no.object <- missing(object)
  .pltest_local({
    bw.method <- match.arg(bw.method, c("sun", "cv", "object"))
    boot.method <- match.arg(boot.method, c("wild-rademacher", "wild"))
    boot.correction <- match.arg(boot.correction, c("none", "cv-pilot", "pilot-hc2", "cv-raw"))
    projection.scale <- match.arg(projection.scale, c("original", "standardized"))
    scalar.integer <- function(x) {
      is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) == 1L &&
        is.finite(x) && abs(x) <= .Machine$integer.max && x == trunc(x)
    }
    if (!scalar.integer(B) || B < 9L)
      stop("B must be a single integer from 9 through .Machine$integer.max.", call. = FALSE)
    B <- as.integer(B)
    if (!is.null(random.seed) && !scalar.integer(random.seed))
      stop("random.seed must be NULL or a single representable integer.", call. = FALSE)
    if (no.object || !inherits(object, c("plregression", "plbandwidth")))
      stop("object must be an npplreg() fit or an npplregbw() bandwidth object.", call. = FALSE)
    controls <- c("bwtype", "bwscaling", "ckertype", "ckerorder", "ckerbound", "ckerlb",
                  "ckerub", "ukertype", "okertype", "nmulti", "itmax", "ftol", "tol",
                  "small", "powell.remin", "scale.factor.search.lower")
    if (...length() && (bw.method != "cv" || length(dot.names) != ...length() ||
        anyDuplicated(dot.names) || !all(dot.names %in% controls)))
      stop("Supply only named, distinct documented CV controls in ... with bw.method = 'cv'.",
           call. = FALSE)
    dots <- list(...)
    .pltest_preflight(B)
  })
  .pltest_same(list(bw.method, boot.method, boot.correction, projection.scale, B, random.seed, dots))
  if (is.null(random.seed)) {
    .pltest_current_stream()
    on.exit(.pltest_current_stream(), add = TRUE)
  } else {
    state <- .pltest_rng_capture()
    on.exit(.np_rng_scope_restore(state), add = TRUE)
    .np_rng_scope_seed(state, random.seed)
  }

  .pltest_local({
    # Use the installed package's own training-sample recovery, including fit overrides.
    if (inherits(object, "plregression")) {
      input.bw <- .np_plreg_bws(object, where = "nppltest")
      data <- .np_plreg_predict_se_data(input.bw, object)
    } else {
      input.bw <- object
      data <- .np_plreg_predict_train_data(object)
    }

    X <- as.data.frame(data$txdat)
    Z <- as.data.frame(data$tzdat)
    y <- data$tydat
    if (is.matrix(y) && ncol(y) == 1L) y <- drop(y)
    if (!is.numeric(y) || !is.null(dim(y)) || nrow(X) != length(y) || nrow(Z) != length(y))
      stop("A scalar numeric training response aligned with X and Z is required.", call. = FALSE)
    if (!ncol(X) || !ncol(Z)) stop("Both X and Z must contain predictors.", call. = FALSE)
    keep <- complete.cases(X, Z, y)
    if (sum(keep) < 2L) stop("At least two complete training observations are required.", call. = FALSE)
    X <- X[keep, , drop = FALSE]
    Z <- Z[keep, , drop = FALSE]
    y <- as.numeric(y[keep])
    predictors <- c(X, Z)
    if (any(!is.finite(y))) stop("The training response contains infinite values.", call. = FALSE)
    valid.predictor <- function(v) {
      is.null(dim(v)) && (is.factor(v) || (is.numeric(v) && all(is.finite(v))))
    }
    if (!all(vapply(predictors, valid.predictor, logical(1))))
      stop("Each training predictor must be a finite numeric vector or factor.", call. = FALSE)
    if (length(y) <= ncol(X) + 1)
      stop("Too few complete observations for the linear slopes and nuisance intercept.", call. = FALSE)
    continuous <- !vapply(Z, is.factor, logical(1))
    categories <- Z[, !continuous, drop = FALSE]
    cells <- if (ncol(categories)) {
      do.call(interaction, c(unname(lapply(categories, as.integer)), list(drop = TRUE)))
    } else factor(rep(1L, length(y)))
  })
  .pltest_same(list(lapply(X, identity), lapply(Z, identity), y))
  bw <- switch(bw.method,
               sun = .np_progress_activity_run("Preparing Sun bandwidths",
                 .pltest_sun_bandwidth(X, y, Z, continuous, cells)),
               cv = do.call(npplregbw, c(list(xdat = X, ydat = y, zdat = Z,
                                      regtype = if (any(continuous)) "ll" else "lc",
                                      bwmethod = "cv.ls"), dots)),
               object = input.bw)
  if (bw.method == "object") {
    # Check shared kernel, basis and bandwidth settings before collective fits.
    bw <- .pltest_selected_bw(bw)
  }
  fitted.null <- .np_progress_activity_run("Fitting the partially linear null",
    .pltest_operator(bw, X, y, Z, shared.sun = bw.method == "sun"))
  .pltest_local({
    H <- fitted.null$H
    mu <- fitted.null$mean
    e <- y - mu
    n <- length(y)
    if (sum(.pltest_residual_norms(H)) <= n * .Machine$double.eps ||
        max(abs(e)) <= 64 * .Machine$double.eps * max(abs(y)))
      stop("The fitted null leaves no numerical residual variation to test.", call. = FALSE)
  })
  A <- .np_progress_activity_run("Computing projection weights",
    .pltest_angles(predictors, projection.scale == "standardized"))
  .pltest_local({
    statistic <- function(E) colSums(E * (A %*% E))
    observed <- unname(statistic(as.matrix(e)))
    if (any(!is.finite(A)) || !is.finite(observed))
      stop("Non-finite projection statistic. Check the scales of the response and predictors.", call. = FALSE)
    if (observed == 0 && max(abs(e)) < sqrt(.Machine$double.xmin))
      stop("The observed quadratic form underflows in these response units. Rescale Y.", call. = FALSE)
  })
  pilot <- NULL
  if (boot.correction != "none") {
    pilot <- .np_progress_activity_run("Fitting the bootstrap pilot",
      if (boot.correction %in% c("cv-pilot","cv-raw"))
      .pltest_local(.pltest_cv_pilot(X, y, Z, adjust = boot.correction == "cv-pilot")) else
       .pltest_local(.pltest_pilot(X, y, Z, continuous, cells, bw)))
  }
  .pltest_local({
    center <- if (is.null(pilot)) mu else pilot$center
    amplitude <- if (is.null(pilot)) e else pilot$amplitude
    if (any(!is.finite(amplitude))) stop("Non-finite bootstrap residual amplitude.", call. = FALSE)
  })
  draws <- .pltest_bootstrap(H, A, center, amplitude, B, boot.method)
  boot <- draws$statistics
  if (length(boot) != B || any(!is.finite(boot)))
    stop("Expected B finite bootstrap statistics. Check response and predictor scales.", call. = FALSE)
  exceed <- .np_bootstrap_upper_tail_count(boot, observed)
  bandwidths <- do.call(rbind, lapply(bw[["bw"]], function(child) as.numeric(child[["bw"]])))
  dimnames(bandwidths) <- list(c("Y", names(X)), names(Z))
  structure(list(call = call, model.call = object[["call"]], statistic = observed, P = exceed / B,
                 exceedances = exceed, critical.values = quantile(boot, c(.90,.95,.99), type = 1),
                 boot.correction = boot.correction, pilot = pilot,
                 bootstrap = boot, rank.counts = draws$rank.counts, boot.method = boot.method, boot.num = B, random.seed = random.seed,
                 bws = bw, bandwidths = bandwidths, bw.method = bw.method,
                 projection.scale = projection.scale, n = n, x.names = names(X), z.names = names(Z),
                 continuous.z = sum(continuous), cell.counts = table(cells),
                 factor.levels = lapply(categories, levels),
                 total.time = proc.time()[["elapsed"]] - started), class = "pltest")
}

print.pltest <- function(x, ...) {
  cat("Projection test of partial linearity\n")
  .np_spec_test_print_result(x)
  invisible(x)
}

summary.pltest <- function(object, ...) {
  print(object)
  cat("Linear predictors:", paste(object$x.names, collapse = ", "), "\n")
  cat("Nuisance predictors:", paste(object$z.names, collapse = ", "), "\n")
  cat("Bandwidth method:", object$bw.method, "   Projection scale:", object$projection.scale, "\n")
  cat("Wild bootstrap:", object$boot.method, "  Draws:", object$boot.num, "\n")
  cat("Residual correction:", object$boot.correction, "\n")
  if (object$boot.correction %in% c("cv-pilot", "cv-raw")) {
    if (identical(object$boot.correction, "cv-pilot"))
      cat(sprintf("CV npplreg pilot, mean leverage %.3f, minimum residual row norm %.3g\n",
                mean(object$pilot$leverage), sqrt(min(object$pilot$variance.factor))))
    widths <- do.call(rbind, lapply(object$pilot$bws[["bw"]], function(b) as.numeric(b[["bw"]])))
    dimnames(widths) <- dimnames(object$bandwidths)
    cat("Pilot", object$pilot$bws[["pscaling"]], "(", object$pilot$bws[["type"]], ")\n")
    print(widths)
  }
  if (identical(object$boot.correction, "pilot-hc2"))
    cat(sprintf("Pilot-centred HC2 refit, %s, degree %d, rank %d\n",
                object$pilot$shape, object$pilot$degree, object$pilot$rank))
  cat("Test", object$bws[["pscaling"]], "(", object$bws[["type"]], ")\n")
  print(object$bandwidths)
  if (length(object$factor.levels)) {
    cat("Observed nuisance cells (indices in declared level order):\n")
    print(object$cell.counts)
  }
  .np_test_print_time(object)
  invisible(object)
}

.pltest_rng_capture <- function() .np_rng_scope_capture()

.pltest_preflight <- function(B, where = "nppltest") invisible(NULL)

.pltest_bootstrap <- function(H, A, center, amplitude, B, boot.method) {
  progress <- .np_progress_activity_begin(B, "Bootstrapping partial linearity", detail = NULL)
  complete <- FALSE
  on.exit(.np_progress_activity_end(progress, complete), add = TRUE)
  values <- numeric(B)
  for (first in seq.int(1, B, by = 64)) {
    index <- seq.int(first, min(B, as.double(first) + 63))
    values[index] <- .pltest_draw_block(length(index), H, A, center, amplitude, boot.method)
    .np_progress_activity_step(progress, done = index[length(index)])
  }
  complete <- TRUE
  list(statistics = values, rank.counts = NULL)
}

.pltest_local <- function(expr, where = "nppltest") force(expr)
.pltest_same <- function(value) invisible(NULL)
.pltest_current_stream <- function() invisible(NULL)
.pltest_angles <- function(predictors, standardize) .pltest_geometry(predictors, standardize)
.pltest_selected_bw <- function(bw) bw

# Use the same strict cutoffs as npsigtest, including attainable boundary p-values.
.np_spec_test_print_result <- function(x) {
  stars <- isTRUE(getOption("show.signif.stars", TRUE))
  symbol <- if (stars)
    c(" ***", " **", " *", " .", "")[1L + findInterval(x$P, c(.001, .01, .05, .1))]
  else ""
  cat(sprintf("T: %.6g    P Value: %s%s    n: %d\n",
              x$statistic, .np_test_format_pvalue(x$P, x$boot.num), symbol, x$n))
  if (stars)
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(NULL)
}
