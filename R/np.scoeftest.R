.sct_setup <- function(object, method, jlwy.bandwidth, h, trim) {
  fit <- inherits(object, "smoothcoefficient")
  if (!fit && !inherits(object, "scbandwidth"))
    stop("Supply an npscoef training fit or npscoefbw bandwidth object.", call. = FALSE)
  if (fit && !isTRUE(object[["trainiseval"]]))
    stop("Supply a training-sample fit without external evaluation data.", call. = FALSE)
  if (fit && isTRUE(object[["fit.controls"]][["leave.one.out"]]))
    stop("Use an in-sample fit with leave.one.out = FALSE.", call. = FALSE)
  b <- if (fit) object[["bws"]] else object
  if (!identical(b[["type"]], "fixed") || !identical(b[["regtype"]], "lc") ||
      !identical(b[["ckertype"]], "gaussian") || !isTRUE(b[["ckerorder"]] == 2) ||
      !identical(b[["ckerbound"]], "none") || !is.null(b[["bw.fitted"]]))
    stop("Use a non-iterated, unbounded second-order Gaussian local-constant model with fixed bandwidths.", call. = FALSE)
  tr <- b[[".np.native.training"]]
  if (is.null(tr)) {
    frame <- b[[".np.formula.training"]][["frame"]]
    roles <- b[["chromoly"]]
    if (is.null(frame) || length(roles) != 3L)
      stop("The bandwidth object must retain training data with separate X and Z roles.", call. = FALSE)
    tr <- list(xdat = frame[roles[[2L]]], zdat = frame[roles[[3L]]], ydat = stats::model.response(frame))
  }
  x <- as.data.frame(tr[["xdat"]])
  z <- as.data.frame(tr[["zdat"]])
  y <- tr[["ydat"]]
  n <- length(y)
  plain <- function(v) is.numeric(v) && !is.object(v) && is.null(dim(v)) && all(is.finite(v))
  if (!plain(y) || nrow(x) != n || nrow(z) != n || !ncol(x) || !ncol(z) ||
      n <= ncol(x) + 1L || !all(vapply(c(x, z), plain, TRUE)))
    stop("Supply finite numeric Y and separate continuous X and Z columns with enough training observations.", call. = FALSE)
  w <- cbind(x, z)
  scales <- vapply(w, stats::sd, 0.0)
  if (any(!is.finite(scales) | scales <= 0))
    stop("Remove constant predictors or rescale nonfinite predictor scales.", call. = FALSE)
  student <- method == "jlwy"
  input.a <- unlist(b[["bandwidth"]], use.names = FALSE)
  a <- input.a
  if (student) a <- vapply(z, stats::sd, 0.0) * n^(-1/6)
  if (length(a) != ncol(z) || !plain(a) || any(a <= 0))
    stop("The null bandwidths must be finite and positive.", call. = FALSE)
  if (is.null(h)) {
    scale <- if (!student) sqrt(2)*n^(-1/5) else if (jlwy.bandwidth == "application") n^(-1/6) else n^(-1/5)
    h <- scales * scale
  } else if (!is.null(names(h))) {
    if (anyDuplicated(names(w)) || anyDuplicated(names(h)) || !setequal(names(h), names(w)))
      stop("Named h must uniquely match all X and Z names; use an unnamed vector for positional matching.", call. = FALSE)
    h <- h[names(w)]
  }
  if (!plain(h) || length(h) != ncol(w) || any(h <= 0))
    stop("Supply one positive finite pair bandwidth per coordinate, ordered as X then Z.", call. = FALSE)
  if (!is.null(trim) && (!student || !is.logical(trim) || !is.null(dim(trim)) ||
      length(trim) != n || anyNA(trim) || sum(trim) < 2L))
    stop("For method = 'jlwy', trim must retain at least two rows in a logical vector aligned to the training data.", call. = FALSE)
  bw <- .npscoefbw_build_scbandwidth(x, y, z, a, FALSE,
             list(regtype = "lc", bwmethod = "manual", bandwidth = a))
  list(x = x, z = z, y = y, bw = bw, pooled = .npscoefbw_fast_eligible(bw, z),
       w = w, a = a, h = h, student = student, trim = trim, input.a = input.a,
       input.fit = if (fit && !student) as.numeric(fitted(object)))
}

.sct_multipliers <- function(n, B) {
  phi <- sqrt(5)
  u <- stats::runif(as.double(n)*B)
  weights <- c((1-phi)/2, (1+phi)/2)
  weights[1L + (u >= (1+phi)/(2*phi))]
}

.sct_wild <- function(n, B) {
  draws <- .sct_multipliers(n, B)
  dim(draws) <- c(n, B)
  draws
}

npscoeftest <- function(object, method = c("wide", "jlwy"), B = 399L,
                        random.seed = 42L, h = NULL,
                        jlwy.bandwidth = c("application", "simulation"), trim = NULL, ...) {
  call <- match.call()
  jlwy.control <- !missing(jlwy.bandwidth)
  started <- proc.time()[["elapsed"]]
  ints <- function(v) is.numeric(v) && !is.object(v) && is.null(dim(v)) && length(v) == 1L &&
    is.finite(v) && abs(v) <= .Machine$integer.max && v == trunc(v)
  .pltest_local({
    if (...length()) stop("Unused arguments to npscoeftest().", call. = FALSE)
    method <- match.arg(method, c("wide", "jlwy"))
    jlwy.bandwidth <- match.arg(jlwy.bandwidth, c("application", "simulation"))
    if (jlwy.control && method != "jlwy")
      stop("jlwy.bandwidth applies only to method = 'jlwy'.", call. = FALSE)
    if (!ints(B) || B < 9L || B >= .Machine$integer.max)
      stop("B must be an integer from 9 through 2147483646.", call. = FALSE)
    if (!is.null(random.seed) && !ints(random.seed))
      stop("random.seed must be NULL or a representable integer.", call. = FALSE)
    .pltest_preflight(B, "npscoeftest")
    s <- .sct_setup(object, method, jlwy.bandwidth, h, trim)
  }, "npscoeftest")
  .sct_agree(list(method, jlwy.bandwidth, as.integer(B), random.seed, s$x, s$z, s$y,
                  s$a, s$h, s$pooled, s$trim, s$input.fit, s$input.a))
  if (is.null(random.seed)) {
    .pltest_current_stream()
    on.exit(.pltest_current_stream(), add = TRUE)
  } else {
    state <- .pltest_rng_capture()
    on.exit(.np_rng_scope_restore(state), add = TRUE)
    .np_rng_scope_seed(state, random.seed)
  }
  n <- length(s$y)
  # Scale before forming quadratic products; apply the common units only once.
  # The null includes an intercept. Remove the response location before
  # scaling so a change of origin does not look like zero residual variation.
  centre <- min(s$y)/2 + max(s$y)/2
  y <- s$y - centre
  size <- max(abs(y))
  if (!size) stop("The null has zero residual variation.", call. = FALSE)
  y <- y/size
  if (!is.null(s$input.fit)) {
    s$fit.tol <- sqrt(.Machine$double.eps) + 64*.Machine$double.eps*max(abs(s$y))/size
    s$input.fit <- (s$input.fit-centre)/size
  }
  scored <- .sct_compute(s, y, as.integer(B))
  t <- scored$values
  K <- .np_bootstrap_upper_tail_count(t[-1L], t[1L])
  reported <- t
  log.multiplier <- 0
  if (!s$student) {
    # The raw ordered-pair sum can be negative. Preserve signs even when
    # conversion back to response units would overflow or underflow.
    log.factor <- 2*log(size)
    physical <- sign(t)*exp(log(abs(t))+log.factor)
    if (all(is.finite(physical)) && all(t == 0 | abs(physical) >= .Machine$double.xmin) &&
        .np_bootstrap_upper_tail_count(physical[-1L], physical[1L]) == K)
      reported <- physical
    else log.multiplier <- log.factor
  }
  structure(list(call = call, statistic = reported[1L], bootstrap = reported[-1L],
    log.multiplier = log.multiplier, rng.kind = RNGkind(),
    P = K/B, exceedances = K, boot.num = as.integer(B), random.seed = random.seed,
    critical.values = stats::quantile(reported[-1L], c(.90,.95,.99), type = 1),
    n = n, method = method, jlwy.bandwidth = if (s$student) jlwy.bandwidth,
    n.retained = if (is.null(trim)) n else sum(trim),
    input.bandwidth = s$input.a, null.bandwidth = s$a,
    effective.df = scored$effective.df, h = setNames(as.numeric(s$h), names(s$w)),
    rank.counts = scored$counts, total.time = proc.time()[["elapsed"]] - started), class = "scoeftest")
}

print.scoeftest <- function(x, ...) {
  cat(if (x$method == "wide") "Unstudentized smooth-coefficient specification test\n"
      else "Jansen-Li-Wang-Yang smooth-coefficient specification test\n")
  if (x$log.multiplier != 0)
    cat("Reported scores are scaled. Log multiplier for original units:", format(x$log.multiplier), "\n")
  .np_spec_test_print_result(x)
  invisible(x)
}

summary.scoeftest <- function(object, ...) {
  print(object)
  cat("Null bandwidths:", format(object$null.bandwidth), "\n")
  cat("Pair bandwidths:", paste(names(object$h), format(object$h), collapse = ", "), "\n")
  if (object$method == "jlwy")
    cat("JLWY convention:", object$jlwy.bandwidth, "   Retained observations:", object$n.retained, "\n")
  cat("Mammen wild bootstrap draws:", object$boot.num, "\n")
  .np_test_print_time(object)
  invisible(object)
}

# Split the residual bank once. Align donor blocks to target panels so a
# target belongs to exactly one block. Earlier donors can go straight to BLAS.
.sct_score_bank <- function(residuals, width, student) {
  n <- nrow(residuals)
  block <- width * max(1L, floor(1024L/width))
  starts <- seq.int(1L, n, by = block)
  pieces <- if (length(starts) == 1L) list(residuals) else
    lapply(starts, function(first)
      residuals[seq.int(first, min(n, first+block-1L)), , drop = FALSE])
  list(pieces = pieces, squared = if (student) lapply(pieces, function(x) x*x), starts = starts, block = block, n = n, columns = ncol(residuals))
}

.sct_panel_score <- function(s, bank, panel) {
  first <- (panel - 1) * s$width + 1
  last <- min(bank$n, panel * s$width)
  target <- seq.int(first, last)
  final <- (first - 1L) %/% bank$block + 1L
  local.target <- target - bank$starts[final] + 1L
  target.residuals <- bank$pieces[[final]][local.target, , drop = FALSE]
  numerator <- numeric(bank$columns)
  if (s$student) denominator <- numeric(bank$columns)
  for (part in seq_len(final)) {
    residuals <- bank$pieces[[part]]
    start <- bank$starts[part]
    used <- min(nrow(residuals), last - start + 1L)
    donors <- seq.int(start, length.out = used)
    if (used < nrow(residuals)) residuals <- residuals[seq_len(used), , drop = FALSE]
    distance <- matrix(0, used, length(target))
    for (j in seq_along(s$w)) {
      v <- s$w[[j]]
      delta <- outer(v[donors], v[target], "-") / s$h[j]
      distance <- distance + delta*delta
    }
    kernel <- exp(-0.5*distance)
    if (part == final) kernel[cbind(local.target, seq_along(target))] <- 0
    # Earlier donor blocks occur once in this triangular traversal. The
    # target-by-target square already includes both orders of its pairs.
    if (s$student) kernel2 <- kernel*kernel
    if (part == final) {
      kernel[local.target, ] <- 0.5*kernel[local.target, , drop = FALSE]
      if (s$student) kernel2[local.target, ] <- 0.5*kernel2[local.target, , drop = FALSE]
    }
    numerator <- numerator + colSums(target.residuals * crossprod(kernel, residuals))
    if (s$student) {
      squared <- bank$squared[[part]]
      if (used < nrow(squared)) squared <- squared[seq_len(used), , drop = FALSE]
      denominator <- denominator + colSums(target.residuals^2 * crossprod(kernel2, squared))
    }
  }
  if (s$student) c(2*numerator, 2*denominator) else 2*numerator
}

.sct_check_hat_rows <- function(s, ids, weights, adjoint, moments) {
  # A second calculation is only needed for difficult designs. Compare the
  # entire row, so the check also covers responses yet to be bootstrapped.
  diagonal <- moments[diag(s$gram.map), , drop = FALSE]
  diagonal <- (diagonal / s$x.scale) / s$x.scale
  targets <- t(s$design[ids, , drop = FALSE])
  directional <- colSums(diagonal) * sqrt(colSums((adjoint*s$x.scale)^2) /
                                           colSums((targets/s$x.scale)^2))
  risky <- s$weak.design | !is.finite(directional) | directional > .Machine$double.eps^(-1/4)
  for (j in which(risky)) {
    root <- sqrt(weights[,j])
    local <- qr(s$scaled.design * root, LAPACK = TRUE)
    rhs <- s$scaled.design[ids[j], local$pivot]
    coordinate <- forwardsolve(t(qr.R(local)), rhs)
    reference <- root * drop(qr.Q(local) %*% coordinate)
    actual <- weights[,j] * drop(s$design %*% adjoint[,j])
    relative <- sum(abs(actual-reference)) / max(1, sum(abs(reference)))
    if (!is.finite(relative) || relative > sqrt(.Machine$double.eps))
      stop("The local coefficient design is numerically ill-conditioned. Reparameterize X or increase the null bandwidth.", call. = FALSE)
  }
  invisible(NULL)
}

.sct_task <- function(task, phase, s, bank, geometry) {
  ids <- seq.int(task$start, length.out = task$bsz)
  if (phase == "score") {
    scores <- t(vapply(ids, function(i) .sct_panel_score(s, bank, i), numeric(bank$columns * (1L + s$student))))
    return(cbind(scores, rank = 0))
  }
  targets <- if (s$pooled) rep.int(1L, length(ids)) else ids
  weights <- .np_kernel_weights_direct(s$kbw, txdat = s$z,
                 exdat = s$z[targets, , drop = FALSE], bandwidth.divide = TRUE)
  if (phase == "fit") {
    p <- ncol(s$design)
    moments <- crossprod(geometry, weights)
    gram <- array(moments[as.vector(s$gram.map), , drop = FALSE], c(p, p, length(ids)))
    # Solve M_i v_i = D_i with npscoef's own LAPACK rank check. Then
    # H_ij = v_i' D_j K_ji. The same adjoint applies to every wild draw.
    adjoint <- .npscoef_batch_zero_solve(t(s$design[ids, , drop = FALSE]), gram)
    if (is.null(adjoint))
      stop("A local coefficient design requires ridge regularization. The prototype supports nonsingular local fits; increase the null bandwidth or simplify X.", call. = FALSE)
    .sct_check_hat_rows(s, ids, weights, adjoint, moments)
    fitted <- colSums(adjoint * moments[s$gram.terms + seq_len(p), , drop = FALSE])
    leverage <- weights[cbind(ids, seq_along(ids))] *
                colSums(t(s$design[ids, , drop = FALSE]) * adjoint)
    return(cbind(fitted, t(adjoint), leverage, rank = 0))
  }
  hat <- weights * tcrossprod(s$design, geometry[ids, , drop = FALSE])
  cbind(bank[ids, , drop = FALSE] - crossprod(hat, bank), rank = 0)
}

.sct_engine <- function(s, y, B) {
  n <- length(y)
  .pltest_local({
    design <- cbind(1, toMatrix(s$x))
    p <- ncol(design)
    s$design <- design
    s$x.scale <- apply(abs(design), 2L, max)
    s$scaled.design <- sweep(design, 2L, s$x.scale, "/")
    # Normal equations square the condition number. When their potential
    # loss exceeds sqrt(eps), verify rather than changing the estimator.
    s$weak.design <- rcond(qr.R(qr(s$scaled.design, LAPACK = TRUE))) < .Machine$double.eps^(1/4)
    s$kbw <- kbandwidth(s$bw)
    # Up to 128 targets fit an efficient BLAS panel. Independent caps keep a
    # single geometry panel under 32 MiB and a response block under 128 MiB
    # whenever one column itself fits those limits.
    s$width <- .sct_width(n, max(1L, min(128L, floor(2^22/n))))
    terms <- seq_len(p)
    # Each symmetric Gram product is accumulated once. Expand only the
    # small target-panel Grams when passing them to the package solver.
    index <- matrix(0L, p, p)
    lower <- which(lower.tri(index, diag = TRUE), arr.ind = TRUE)
    s$gram.terms <- nrow(lower)
    index[lower] <- seq_len(s$gram.terms)
    index[upper.tri(index)] <- t(index)[upper.tri(index)]
    s$gram.map <- index
    features <- cbind(design[, lower[,1L], drop = FALSE] *
                      design[, lower[,2L], drop = FALSE], design*y)
  }, "npscoeftest")
  fit <- .sct_map("fit", s, NULL, features, n, s$width, p+3L, "Fitting the smooth-coefficient null")
  .pltest_local({
    mu <- fit[, 1L]
    if (!is.null(s$input.fit) && (length(s$input.fit) != n ||
        !all(is.finite(s$input.fit)) || max(abs(mu-s$input.fit)) > s$fit.tol))
      stop("The test cannot reproduce the supplied training fit within numerical tolerance. Check the training inputs or reparameterize ill-conditioned X.", call. = FALSE)
    e <- y - mu
    if (max(abs(e)) <= 64*.Machine$double.eps)
      stop("The null has numerically zero residual variation.", call. = FALSE)
    adjoint <- fit[, 1L + terms, drop = FALSE]
    batch <- max(1L, min(B, floor(2^24/n)))
    scores <- numeric(as.double(B) + 1)
    owners <- integer()
  }, "npscoeftest")
  for (first in seq.int(1L, B, by = batch)) {
    cols <- seq.int(first, min(as.double(B), first + batch - 1))
    wild <- .sct_weights(n, length(cols))
    response <- .pltest_local(mu + e*wild, "npscoeftest")
    rm(wild)
    res <- .sct_map("refit", s, response, adjoint, n, s$width, length(cols)+1L,
                    "Refitting bootstrap responses")
    .pltest_local({
      rm(response)
      residuals <- res[, seq_along(cols), drop = FALSE]
      if (first == 1L) residuals <- cbind(e, residuals)
      rm(res)
      columns <- ncol(residuals)
      if (!is.null(s$trim)) residuals[!s$trim, ] <- 0
      bank <- .sct_score_bank(residuals, s$width, s$student)
      rm(residuals)
    }, "npscoeftest")
    components <- columns * (1L + s$student)
    reduced <- .sct_map("score", s, bank, NULL, ceiling(n/s$width), 1L,
                        components+1L, "Computing off-diagonal kernel scores")
    rm(bank)
    .pltest_local({
      value <- colSums(reduced[, seq_len(columns), drop = FALSE])
      if (s$student) {
        variance <- 2*colSums(reduced[, columns + seq_len(columns), drop = FALSE])
        if (any(!is.finite(variance) | variance <= 0))
          stop("The JLWY variance is unresolved. Check pair bandwidths and retained observations.", call. = FALSE)
        value <- value/sqrt(variance)
      }
      if (any(!is.finite(value)))
        stop("Nonfinite off-diagonal score. Check predictor and response scales.", call. = FALSE)
      target <- if (first == 1L) c(1L, cols+1) else cols+1
      scores[target] <- value
      if (first == 1L) owners <- table(reduced[, components+1L])
    }, "npscoeftest")
  }
  if (all(scores == 0)) stop("No resolved off-diagonal scores. Check bandwidths and design.", call. = FALSE)
  list(values = scores, counts = owners, effective.df = sum(fit[, p+2L]))
}

.sct_agree <- function(value) invisible(NULL)
.sct_width <- function(n, limit) limit
.sct_weights <- .sct_wild

.sct_map <- function(phase, s, bank, geometry, total, chunk, columns, label) {
  result <- matrix(NA_real_, total, columns)
  progress <- .np_progress_activity_begin(total, label, detail = NULL)
  complete <- FALSE
  on.exit(.np_progress_activity_end(progress, complete), add = TRUE)
  for (first in seq.int(1L, total, by = chunk)) {
    count <- min(chunk, total - first + 1L)
    task <- list(start = first, bsz = count)
    rows <- seq.int(first, length.out = count)
    result[rows, ] <- .sct_task(task, phase, s, bank, geometry)
    .np_progress_activity_step(progress, done = first + count - 1L)
  }
  complete <- TRUE
  result
}

.sct_compute <- .sct_engine
