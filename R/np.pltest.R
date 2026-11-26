# Statistic: Sun, Chen, Liang and Ruppert (2021), doi:10.5705/ss.202019.0238.
.pltest_coordinates <- function(predictors, standardize, compact = FALSE) {
  columns <- lapply(predictors, function(x) {
    if (is.factor(x) && !is.ordered(x)) {
      code <- as.integer(droplevels(x))
      if (compact) {
        levels <- max(code)
        if (levels == 1L) return(NULL)
        if (levels == 2L) return(code - 1.5)
        if (levels > 1112L)
          stop("Too many unordered levels for Sobol integration; use integration = 'exact'.", call. = FALSE)
        j <- seq_len(levels - 1L)
        basis <- sweep(stats::contr.helmert(levels), 2L, sqrt(2 * j * (j + 1)), "/")
        return(basis[code, , drop = FALSE])
      }
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
  V
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

.pltest_with_local_fit <- function(expr) .npRmpi_with_local_bootstrap(.np_plot_with_local_compiled_eval(expr))

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
  e <- .pltest_residual_block(u, H, center, amplitude, boot.method)
  colSums(e*(A %*% e))
}

# A fixed Cranley-Patterson shift preserves the bootstrap stream and makes
# direction budgets nested. Generating numbers and shift provenance are in
# data-raw/make-pltest-sobol.R.
.pltest_sobol_directions <- function(m, d) {
  table <- .pltest_sobol_table
  if (d > ncol(table$words))
    stop("Sobol integration supports up to 1111 projection coordinates; use integration = 'exact'.", call. = FALSE)
  i <- seq_len(m)
  gray <- bitwXor(i, bitwShiftR(i, 1L))
  u <- matrix(0L, m, d)
  for (k in seq_len(min(31L, floor(log2(m)) + 1L))) {
    rows <- which(bitwAnd(gray, bitwShiftL(1L, k - 1L)) != 0L)
    for (j in seq_len(d))
      u[rows, j] <- bitwXor(u[rows, j], table$words[k, j])
  }
  u <- sweep(u / 2^31, 2L, table$shift[seq_len(d)], "+") %% 1
  # A shifted binary fraction can land exactly on the left endpoint.
  u[u == 0] <- .Machine$double.eps / 2
  z <- stats::qnorm(u)
  z / sqrt(rowSums(z^2))
}

# One rank owns all directions for each residual column it receives. Order
# the observations once per direction and reuse that order across its draws.
.pltest_direction_scores <- function(E, V, directions) {
  n <- nrow(E)
  magnitude <- matrixStats::colMaxs(abs(E), useNames = FALSE)
  if (any(!is.finite(magnitude)))
    stop("Non-finite bootstrap or observed residuals. Rescale the response.", call. = FALSE)
  extreme <- magnitude > 1e100 | (magnitude > 0 & magnitude < 1e-100)
  if (any(extreme)) {
    scale <- 2^pmin(1023, floor(log2(magnitude[extreme])))
    E[, extreme] <- sweep(E[, extreme, drop = FALSE], 2L, scale, "/")
  }
  ans <- numeric(ncol(E))
  for (j in seq_len(nrow(directions))) {
    if (j %% 16L == 1L) .np_fit_progress_step()
    z <- drop(V %*% directions[j, ])
    if (any(!is.finite(z))) stop("Non-finite projected coordinates; rescale the predictors.", call. = FALSE)
    ord <- order(z, method = "radix")
    sorted <- z[ord]
    ends <- if (anyDuplicated(sorted)) c(which(sorted[-1L] != sorted[-n]), n) else integer()
    ans <- ans + .Call(C_np_pltest_prefix_scores, E, ord, ends)
  }
  ans <- ans / (as.double(n)^2 * nrow(directions))
  if (any(extreme)) ans[extreme] <- (ans[extreme] * scale) * scale
  ans
}

.pltest_residual_block <- function(u, H, center, amplitude, boot.method) {
  xi <- if (boot.method == "wild-rademacher") ifelse(u < 0.5, -1, 1) else
    ifelse(u < (sqrt(5)+1)/(2*sqrt(5)), (1-sqrt(5))/2, (1+sqrt(5))/2)
  ystar <- center + amplitude*xi
  ystar - H %*% ystar
}

.pltest_sobol_bootstrap <- function(V, directions, H, e, center, amplitude, B, boot.method) {
  n <- length(e)
  # Bound the residual bank independently of B. The dense null operator
  # stays with its owner and is never a projection-scoring payload.
  width <- max(1L, min(512L, floor(4e6 / n)))
  values <- numeric(B)
  owners <- integer(B)
  for (first in seq.int(1, B, by = width)) {
    ids <- seq.int(first, min(B, as.double(first) + width - 1))
    include.observed <- first == 1
    E <- .np_progress_activity_run("Preparing bootstrap residuals", .pltest_local(if (.pltest_owns_residuals()) {
      E <- matrix(0, n, length(ids) + include.observed)
      if (include.observed) E[, 1L] <- e
      for (start in seq.int(1L, length(ids), by = 64L)) {
        cols <- seq.int(start, min(start + 63L, length(ids)))
        u <- matrix(stats::runif(as.double(n) * length(cols)), n)
        E[, cols + include.observed] <- .pltest_residual_block(u, H, center, amplitude, boot.method)
        .np_fit_progress_step()
      }
      E
    }))
    scored <- .pltest_score_bank(E, V, directions, length(ids) + include.observed)
    if (include.observed) observed <- scored[1L, 1L]
    rows <- seq_along(ids) + include.observed
    values[ids] <- scored[rows, 1L]
    owners[ids] <- scored[rows, 2L]
  }
  ranks <- table(owners)
  list(observed = observed, statistics = values,
       rank.counts = if (all(owners == -1L)) NULL else setNames(as.integer(ranks), names(ranks)))
}

.pltest_owns_residuals <- function()
  !.npRmpi_cms_collective_context() || mpi.comm.rank(1L) == 0L

.pltest_score_task <- function(task, V, directions) {
  cbind(.pltest_direction_scores(task$E, V, directions), mpi.comm.rank(1L))
}

.pltest_score_bank <- function(E, V, directions, count) {
  size <- mpi.comm.size(1L)
  tasks <- .pltest_local(if (.pltest_owns_residuals()) {
    parts <- .npRmpi_bootstrap_chunk_tasks(count, ceiling(count / size), with.seeds = FALSE)
    lapply(parts, function(task) {
      ids <- seq.int(task$start, length.out = task$bsz)
      task$E <- E[, ids, drop = FALSE]
      task
    })
  })
  if (!.npRmpi_cms_collective_context())
    return(.np_spec_test_fanout(count, tasks, .pltest_score_task, ncol.out = 2L,
      what = "Sobol projection scores", progress.label = "Scoring Sobol projections",
      progress.total = NULL,
      profile.where = "mpi.applyLB:nppltest.sobol", V = V, directions = directions))

  # In a manual collective, scatter only each rank's marked residual columns.
  # Agree any serialization error before starting the scatter.
  prepared <- .pltest_local(if (mpi.comm.rank(1L) == 0L) {
    length(tasks) <- size
    .npRmpi_scatter_prepare(tasks, comm = 1L)
  })
  task <- .npRmpi_scatter_prepared(prepared, root = 0L, comm = 1L)
  local <- .pltest_local({
    out <- matrix(0, count, 2L)
    if (!is.null(task)) {
      ids <- seq.int(task$start, length.out = task$bsz)
      out[ids, ] <- if (mpi.comm.rank(1L) == 0L)
        .np_progress_activity_run("Scoring Sobol projections", .pltest_score_task(task, V, directions))
      else .pltest_score_task(task, V, directions)
    }
    out
  })
  matrix(mpi.allreduce(as.double(local), type = 2L, op = "sum", comm = 1L), count, 2L)
}

nppltest <- function(object, bw.method = c("sun", "cv", "object"),
                     boot.method = c("wild-rademacher", "wild"), B = 399L,
                     random.seed = 42L, projection.scale = c("original", "standardized"), ...,
                     boot.correction = c("none", "cv-pilot", "pilot-hc2", "cv-raw"),
                     integration = c("sobol", "exact"), directions = 2048L) {
  started <- proc.time()[["elapsed"]]
  call <- match.call()
  dot.names <- names(match.call(expand.dots = FALSE)$...)
  no.object <- missing(object)
  supplied.count <- !missing(directions)
  .pltest_local({
    bw.method <- match.arg(bw.method, c("sun", "cv", "object"))
    boot.method <- match.arg(boot.method, c("wild-rademacher", "wild"))
    boot.correction <- match.arg(boot.correction, c("none", "cv-pilot", "pilot-hc2", "cv-raw"))
    projection.scale <- match.arg(projection.scale, c("original", "standardized"))
    scalar.integer <- function(x) {
      is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) == 1L &&
        is.finite(x) && abs(x) <= .Machine$integer.max && x == trunc(x)
    }
    integration <- match.arg(integration, c("sobol", "exact"))
    if (integration == "exact" && supplied.count)
      stop("directions applies only to integration = 'sobol'.", call. = FALSE)
    if (!scalar.integer(directions) || directions < 1L)
      stop("directions must be a single positive representable integer.", call. = FALSE)
    directions <- as.integer(directions)
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
  .pltest_same(list(bw.method, boot.method, boot.correction, projection.scale, integration, directions, B, random.seed, dots))
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
  if (integration == "exact") {
    A <- .pltest_angles(predictors, projection.scale == "standardized")
    .pltest_local({
      statistic <- function(E) colSums(E * (A %*% E))
      observed <- unname(statistic(as.matrix(e)))
      if (any(!is.finite(A)) || !is.finite(observed))
        stop("Non-finite projection statistic. Check the scales of the response and predictors.", call. = FALSE)
      if (observed == 0 && max(abs(e)) < sqrt(.Machine$double.xmin))
        stop("The observed quadratic form underflows in these response units. Rescale Y.", call. = FALSE)
    })
  }
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
  if (integration == "sobol") {
    V <- .pltest_local(.pltest_coordinates(predictors, projection.scale == "standardized", compact = TRUE))
    net <- .pltest_local(.pltest_sobol_directions(directions, ncol(V)))
    draws <- .pltest_sobol_bootstrap(V, net, H, e, center, amplitude, B, boot.method)
    observed <- draws$observed
    if (!is.finite(observed))
      stop("Non-finite projection statistic. Check the scales of the response and predictors.", call. = FALSE)
    if (observed == 0 && max(abs(e)) < sqrt(.Machine$double.xmin))
      stop("The observed quadratic form underflows in these response units. Rescale Y.", call. = FALSE)
  } else {
    draws <- .pltest_bootstrap(H, A, center, amplitude, B, boot.method)
  }
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
                 integration = integration, directions = if (integration == "sobol") directions else NULL,
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
  if (identical(object$integration, "sobol"))
    cat("Projection integration: shifted Sobol (", object$directions, " directions)\n", sep = "")
  else cat("Projection integration: exact\n")
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

.pltest_rng_capture <- function() {
  .np_rng_scope_capture(collective = .npRmpi_cms_collective_context(), same.kind = TRUE)
}

# Agree local validation before any participant enters the next MPI phase.
.pltest_local <- function(expr, where = "nppltest") {
  if (!.npRmpi_cms_collective_context()) return(force(expr))
  result <- tryCatch(expr, error = identity, interrupt = identity)
  failed <- inherits(result, "condition")
  rank <- mpi.comm.rank(1L)
  size <- mpi.comm.size(1L)
  owner <- mpi.allreduce(if (failed) rank else size, type = 1L, op = "min", comm = 1L)
  if (owner < size) {
    cause <- if (owner == rank) conditionMessage(result) else NULL
    shared <- mpi.bcast.Robj(cause, rank = owner, comm = 1L)
    if (owner != rank) cause <- shared
    stop(where, " failed on MPI rank ", owner, ": ", cause, call. = FALSE)
  }
  result
}

.pltest_root_value <- function(value) {
  answer <- mpi.bcast.Robj(value, rank = 0L, comm = 1L)
  if (mpi.comm.rank(1L) == 0L) value else answer
}

.pltest_same <- function(value) {
  if (.npRmpi_cms_collective_context()) {
    root <- .pltest_root_value(value)
    .pltest_local(if (!identical(value, root))
      stop("all ranks must supply the same test controls, training data and bandwidth specification.", call. = FALSE))
  }
  invisible(NULL)
}

# Under NULL seeding the master owns the current stream. A state transfer
# requires a fully represented RNG on every rank; certify before assignment.
.pltest_current_stream <- function() {
  if (.npRmpi_cms_collective_context()) {
    scope <- .np_rng_scope_capture(collective = TRUE, same.kind = TRUE)
    seed <- .pltest_root_value(scope$seed)
    .pltest_local(.np_rng_scope_state(scope, seed))
  }
  invisible(NULL)
}

.pltest_angle_columns <- function(V, ids) {
  if (mpi.comm.rank(1L) == 0L)
    return(.np_with_compiled_fit_progress("Computing projection weights", nrow(V),
      expr = .Call(C_np_pltest_angle_columns, V, ids)))
  .Call(C_np_pltest_angle_columns, V, ids)
}

.pltest_angles <- function(predictors, standardize) {
  V <- .pltest_local(.pltest_coordinates(predictors, standardize))
  n <- nrow(V)
  collective <- .npRmpi_cms_collective_context()
  ranks <- mpi.comm.size(1L)
  # A second dispatch and matrix transfer dominate short angular calculations.
  if (n < 512L || ranks <= 1L) {
    A <- .pltest_local(if (!collective || mpi.comm.rank(1L) == 0L)
      .np_progress_activity_run("Computing projection weights",
        .np_with_compiled_fit_progress("Computing projection weights", n,
          expr = .Call(C_np_pltest_angles, V))))
    return(if (collective) .pltest_root_value(A) else A)
  }

  if (collective) {
    # Strided columns balance the triangular work. Every rank already has V
    # and needs the completed matrix for its subsequent bootstrap draws.
    own <- .np_progress_activity_run("Computing projection weights",
      .npRmpi_bootstrap_collective_apply(n, n, function(ids, pos)
        .pltest_angle_columns(V, ids)))[[1L]]
    rank <- mpi.comm.rank(1L)
    A <- .pltest_local(matrix(0, n, n))
    for (owner in seq_len(min(n, ranks)) - 1L) {
      ids <- seq.int(owner + 1L, n, by = ranks)
      block <- mpi.bcast.Robj(if (rank == owner) own, rank = owner, comm = 1L)
      .pltest_local(A[, ids] <- if (rank == owner) own else block)
    }
  } else {
    # Pair long and short columns so each task has comparable work and the
    # master participates under the same scheduler as the bootstrap.
    half <- seq_len(ceiling(n/2))
    columns <- as.integer(rbind(half, n - half + 1L))[seq_len(n)]
    width <- 2L * ceiling(n/(2 * ranks))
    tasks <- .npRmpi_bootstrap_chunk_tasks(n, width, with.seeds = FALSE)
    worker <- function(task, V, columns) {
      ids <- columns[seq.int(task$start, length.out = task$bsz)]
      t(.pltest_angle_columns(V, ids))
    }
    rows <- .np_spec_test_fanout(n, tasks, worker, ncol.out = n,
      what = "projection weights", progress.label = "Computing projection weights",
      profile.where = "mpi.applyLB:nppltest.geometry", progress.unit = "column",
      V = V, columns = columns)
    A <- t(rows[order(columns), , drop = FALSE])
  }
  # Mirror only after all owners have finished accumulating the lower triangle.
  for (j in seq_len(n - 1L)) {
    below <- seq.int(j + 1L, n)
    A[j, below] <- A[below, j]
  }
  A
}

.pltest_preflight <- function(B, where = "nppltest") {
  if (isTRUE(.npRmpi_session_fact("npRmpi.profile.active", FALSE)) &&
      !.npRmpi_cms_collective_context())
    stop("Under the MPI startup profile, use mpi.bcast.cmd(..., caller.execute = TRUE).",
         call. = FALSE)
  .npRmpi_require_active_slave_pool(where = paste0(where, "()"))
  if (!.npRmpi_cms_collective_context())
    .npRmpi_bootstrap_fanout_enabled(B = B, what = where)
}

.pltest_bootstrap <- function(H, A, center, amplitude, B, boot.method) {
  if (.npRmpi_cms_collective_context())
    return(.pltest_collective_bootstrap(H, A, center, amplitude, B, boot.method))
  # Generate the same stream as serial R, then replay disjoint draw blocks.
  workers <- .npRmpi_bootstrap_worker_count()
  size <- min(64L, max(1L, ceiling(B/(workers + 1))))
  tasks <- .npRmpi_bootstrap_chunk_tasks(B, size, with.seeds = FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks, function(task)
    invisible(stats::runif(as.double(length(center))*task$bsz)))
  worker <- function(task, H, A, center, amplitude, boot.method) {
    .npRmpi_bootstrap_use_task_rng(task)
    v <- .pltest_draw_block(task$bsz, H, A, center, amplitude, boot.method)
    cbind(statistic = v, rank = mpi.comm.rank(1L))
  }
  out <- .np_spec_test_fanout(B, tasks, worker, ncol.out = 2L,
    what = "partial linearity", progress.label = "Bootstrapping partial linearity",
    profile.where = "mpi.applyLB:nppltest",
    required.bindings = list(.pltest_draw_block = .pltest_draw_block),
    H = H, A = A, center = center, amplitude = amplitude, boot.method = boot.method)
  counts <- table(out[,2])
  list(statistics = as.numeric(out[,1]), rank.counts = setNames(as.integer(counts), names(counts)))
}

# Manual profile broadcasts run on every rank. Rank zero supplies the serial
# draw stream in bounded blocks. The collective batch helper agrees errors
# before any rank enters the result reduction, including ranks with no draws.
.pltest_collective_bootstrap <- function(H, A, center, amplitude, B, boot.method) {
  rank <- mpi.comm.rank(1L)
  size <- mpi.comm.size(1L)
  values <- numeric(B)
  counts <- integer(size)
  progress <- .pltest_local(if (rank == 0L)
    .np_progress_activity_begin(B, "Bootstrapping partial linearity", detail = NULL))
  complete <- FALSE
  on.exit(if (rank == 0L) .np_progress_activity_end(progress, complete), add = TRUE)
  for (first in seq.int(1, B, by = 64)) {
    cols <- seq.int(first, min(B, as.double(first) + 63))
    k <- length(cols)
    if (rank == 0L) {
      u <- tryCatch(matrix(stats::runif(as.double(length(center))*k), nrow = length(center)),
                    error = identity, interrupt = identity)
      mpi.bcast.Robj(u, rank = 0L, comm = 1L)
    } else {
      u <- mpi.bcast.Robj(rank = 0L, comm = 1L)
    }
    if (inherits(u, "condition")) stop(u)
    local <- numeric(k)
    pieces <- .npRmpi_bootstrap_collective_apply(k, k, function(ids, pos) {
      list(ids = ids, statistics = .pltest_block_statistics(
        u[, ids, drop = FALSE], H, A, center, amplitude, boot.method))
    })
    for (piece in pieces) {
      local[piece$ids] <- piece$statistics
      counts[rank+1L] <- counts[rank+1L] + length(piece$ids)
    }
    values[cols] <- mpi.allreduce(local, type = 2L, op = "sum", comm = 1L)
    .pltest_local(if (rank == 0L) .np_progress_activity_step(progress, done = tail(cols, 1L)))
  }

  complete <- TRUE
  counts <- mpi.allreduce(counts, type = 1L, op = "sum", comm = 1L)
  names(counts) <- seq_len(size)-1L
  list(statistics = values, rank.counts = counts[counts > 0L])
}

.pltest_selected_bw <- function(bw) {
  if (.npRmpi_cms_collective_context()) {
    # Calls and model-frame environments are local to each R process.
    fields <- c("bw", "type", "scaling", "bandwidth", "sfactor", "sdev", "nconfac", "ncatfac",
                "regtype", "basis", "degree", "bernstein.basis", "ckertype", "ckerorder",
                "ckerbound", "ckerlb", "ckerub", "ukertype", "okertype", "xdati", "ydati")
    signature <- .pltest_local(lapply(bw[["bw"]], function(b) unclass(b)[fields]))
    .pltest_same(signature)
  }
  bw
}

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

.np_spec_test_fanout <- function(B, ..., progress.label, progress.total = B) {
  context <- .np_progress_activity_begin(progress.total, progress.label, detail = NULL)
  complete <- FALSE
  on.exit(.np_progress_activity_end(context, complete), add = TRUE)
  if (!is.null(context)) {
    context$use.bootstrap.done <- !is.null(progress.total)
    # The master can still be fitting a local draw between MPI completions.
    .np_progress_runtime$fit_forward <- function()
      .np_progress_activity_step(context, if (context$use.bootstrap.done) context$state$last_done)
  }
  out <- .npRmpi_bootstrap_run_fanout(..., progress.label = progress.label,
                                     progress.context = context)
  complete <- TRUE
  out
}
