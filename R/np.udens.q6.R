# Plot-only compact correction under an explicit Student6-smoothed empirical law.
# This owner never fits a public estimator, searches h, or changes an incumbent.
.np_udens_q6_check <- function(a, b, stage) {
  if (!identical(dim(a), dim(b)) || length(a) != length(b) ||
      any(!is.finite(a)) || any(!is.finite(b)))
    stop(paste("Q6 nonfinite or inconsistent", stage), call. = FALSE)
  limit <- 1e-12 + 1e-9 * pmax(abs(a), abs(b))
  error <- max(abs(a - b) / limit)
  if (!is.finite(error) || error > 1)
    stop(sprintf("Q6 %s refinement failed: normalized error %.17g", stage, error),
         call. = FALSE)
  invisible(error)
}

.np_udens_q6_rule <- function(order) {
  j <- seq_len(order - 1L)
  M <- matrix(0, order, order)
  M[cbind(j, j + 1L)] <- j / sqrt(4 * j * j - 1)
  e <- eigen(M + t(M), symmetric = TRUE)
  ix <- order(e$values)
  list(x = e$values[ix], w = 2 * e$vectors[1L, ix]^2)
}

.np_udens_q6_kernel <- function(r, kernel) {
  if (identical(kernel, "uniform"))
    return(ifelse(abs(r) < 1, (225 - 1050 * r*r + 945 * r^4) / 128, 0))
  if (!identical(kernel, "epanechnikov"))
    stop("Q6 requires a uniform or Epanechnikov kernel", call. = FALSE)
  span <- sqrt(5)
  ifelse(abs(r) < span, 3 / (4 * span) * (1 - r*r / 5) *
           (2.734375 - 3.28125 * r*r + .721875 * r^4), 0)
}

.np_udens_q6_budget <- function(n, m, B = 0, block = 32L, order = 96L) {
  sizes <- as.double(c(n, m, B, block, order))
  if (any(!is.finite(sizes)) || any(sizes != floor(sizes)) ||
      n < 10 || m < 1 || B < 0 || block < 1 || order < 2 ||
      any(sizes > .Machine$integer.max))
    stop("invalid Q6 workspace dimensions", call. = FALSE)
  # Worst-case support partition; include simultaneous refinement/factor copies.
  nodes <- as.double(order) * (2 * as.double(m) - 1)
  bytes <- 8 * (12 * as.double(m) * nodes + 20 * as.double(m)^2 +
    8 * as.double(n) * block + 512 * n + 10 * nodes + 4 * n + 4 * B +
    4 * order^2 + 2500 * ceiling(B / block))
  if (!is.finite(bytes) || bytes > 192 * 1024^2)
    stop("Q6 inference exceeds its 192 MiB workspace bound", call. = FALSE)
  bytes
}

.np_udens_q6_geometry <- function(x, h, grid, kernel) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) < 10L ||
      any(!is.finite(x)) || !is.numeric(h) || length(h) != 1L ||
      !is.finite(h) || h <= 0 || !is.numeric(grid) || !length(grid) ||
      any(!is.finite(grid)) || !is.character(kernel) || length(kernel) != 1L ||
      is.na(kernel) || !kernel %in% c("uniform", "epanechnikov"))
    stop("invalid Q6 data, bandwidth or target geometry", call. = FALSE)
  mu <- mean(x)
  v <- mean((x - mu)^2)
  if (!is.finite(v) || v <= 0 || !is.finite(v + h*h))
    stop("Q6 requires positive finite guide variance", call. = FALSE)
  .np_udens_q6_budget(length(x), length(unique(grid)))
  invisible(NULL)
}

.np_udens_q6_density <- function(x, h, u) {
  scale <- h * sqrt(2/3)
  out <- numeric(length(u))
  for (first in seq.int(1L, length(u), by = 128L)) {
    ix <- first:min(length(u), first + 127L)
    out[ix] <- rowMeans(dt(outer(u[ix], x, "-") / scale, df = 6)) / scale
  }
  out
}

.np_udens_q6_convolution <- function(t, mu, sd, h, kernel, order, rule = NULL) {
  if (is.null(rule)) rule <- .np_udens_q6_rule(order)
  span <- if (kernel == "uniform") 1 else sqrt(5)
  r <- span * rule$x
  w <- span * rule$w * .np_udens_q6_kernel(r, kernel)
  out <- numeric(length(t))
  for (j in seq_along(r)) out <- out + w[j] * dnorm(t + h*r[j], mu, sd)
  out
}

# Columns are complete resamples; every column refits its own normal guide.
.np_udens_q6_curve <- function(Y, target, h, kernel, order, guide = NULL, rule = NULL) {
  mu <- if (is.null(guide)) colMeans(Y) else guide$mu
  sd <- if (is.null(guide)) sqrt(colMeans(sweep(Y, 2L, mu)^2)) else guide$sd
  if (any(!is.finite(sd)) || any(sd <= 0))
    stop("Q6 conditional guide variance is not positive finite", call. = FALSE)
  target <- rep(target, length.out = ncol(Y))
  raw <- colMeans(.np_udens_q6_kernel(sweep(Y, 2L, target, "-") / h, kernel)) / h
  out <- raw - .np_udens_q6_convolution(target, mu, sd, h, kernel, order, rule) +
    dnorm(target, mu, sd)
  if (any(!is.finite(out))) stop("nonfinite Q6 corrected curve", call. = FALSE)
  out
}

# T_h(Q6), not q6 itself and not T_h(P_n).
.np_udens_q6_reference_center <- function(x, h, grid, kernel, order) {
  span <- if (kernel == "uniform") 1 else sqrt(5)
  rule <- .np_udens_q6_rule(order)
  r <- span * rule$x
  w <- span * rule$w * .np_udens_q6_kernel(r, kernel)
  u <- outer(grid, h*r, "+")
  mu <- mean(x)
  v <- mean((x - mu)^2) + h*h
  drop(matrix(.np_udens_q6_density(x, h, as.vector(u)) -
                dnorm(as.vector(u), mu, sqrt(v)), nrow = length(grid)) %*% w) +
    dnorm(grid, mu, sqrt(v))
}

.np_udens_q6_covariance <- function(x, h, grid, kernel, order) {
  .np_udens_q6_geometry(x, h, grid, kernel)
  .np_udens_q6_budget(length(x), length(grid), order = order)
  n <- length(x); mu <- mean(x); v <- mean((x - mu)^2); vq <- v + h*h
  span <- if (kernel == "uniform") 1 else sqrt(5)
  cuts <- sort(unique(c(grid - span*h, grid + span*h)))
  if (length(cuts) < 2L || any(!is.finite(cuts)) || any(diff(cuts) <= 0))
    stop("unrepresentable Q6 support partition", call. = FALSE)
  rule <- .np_udens_q6_rule(order)
  u <- w <- numeric()
  for (i in seq_len(length(cuts) - 1L)) {
    l <- cuts[i]; r <- cuts[i + 1L]; mid <- (l + r) / 2
    if (!any(abs(mid - grid) < span*h)) next
    u <- c(u, mid + (r-l) / 2 * rule$x)
    w <- c(w, (r-l) / 2 * rule$w)
  }
  A <- .np_udens_q6_kernel(outer(grid, u, function(t, z) (z-t)/h), kernel) / h
  q <- .np_udens_q6_density(x, h, u)
  ma <- drop(A %*% (w*q))
  CAA <- tcrossprod(sweep(A, 2L, sqrt(w*q), "*")) - tcrossprod(ma)
  IF <- cbind(u-mu, ((u-mu)^2/vq - 1)/2)
  X <- A %*% (IF * (w*q))
  pq <- dnorm(u, mu, sqrt(vq)); gq <- dnorm(grid, mu, sqrt(vq))
  D <- cbind(gq*(grid-mu)/vq, gq*((grid-mu)^2/vq - 1)) -
    A %*% (cbind((u-mu)/vq, (u-mu)^2/vq - 1) * (w*pq))
  m3 <- mean((x-mu)^3)
  m4 <- mean((x-mu)^4) + 6*h*h*v + 6*h^4
  V <- matrix(c(vq, m3/(2*vq), m3/(2*vq), (m4/vq^2 - 1)/4), 2L)
  C <- (CAA + X %*% t(D) + D %*% t(X) + D %*% V %*% t(D)) / n
  list(C = C, mu = mu, v = v, vq = vq, ma = ma, D = D, V = V, X = X,
       m3 = m3, m4 = m4, nodes = length(u), cuts = cuts)
}

.np_udens_q6_reference <- function(x, h, grid, kernel, progress = NULL) {
  .np_udens_q6_geometry(x, h, grid, kernel)
  targets <- sort(unique(grid)); index <- match(grid, targets)
  low <- .np_udens_q6_covariance(x, h, targets, kernel, 64L)
  if (!is.null(progress)) progress("reference covariance", 1L, 2L)
  high <- .np_udens_q6_covariance(x, h, targets, kernel, 96L)
  refinement <- .np_udens_q6_check(low$C, high$C, "covariance")
  .np_udens_q6_check(high$C, t(high$C), "symmetry")
  e <- eigen(high$C, symmetric = TRUE)
  if (any(!is.finite(e$values)) || any(e$values <= 0))
    stop("Q6 covariance is not strictly positive definite", call. = FALSE)
  factor <- sweep(e$vectors, 2L, sqrt(e$values), "*")
  reconstruction <- .np_udens_q6_check(tcrossprod(factor), high$C, "factor reconstruction")
  lo <- .np_udens_q6_reference_center(x, h, targets, kernel, 48L)
  center.q <- .np_udens_q6_reference_center(x, h, targets, kernel, 64L)
  centering <- .np_udens_q6_check(lo, center.q, "reference center")
  center <- vapply(targets, function(t)
    .np_udens_q6_curve(matrix(x, length(x), 1L), t, h, kernel, 64L), 0.0)
  se <- sqrt(diag(high$C))
  if (any(!is.finite(se)) || any(se <= 0))
    stop("Q6 original scale is not positive finite", call. = FALSE)
  if (!is.null(progress)) progress("reference covariance", 2L, 2L)
  list(x = x, h = h, kernel = kernel, grid = targets, index = index,
       center = center, reference.center = center.q, se = se, covariance = high$C,
       diagnostics = c(covariance = refinement, center = centering,
                       reconstruction = reconstruction),
       construction = "density-compact-q6-symmetric-v1")
}

.np_udens_q6_maximum <- function(Y, reference) {
  maximum <- numeric(ncol(Y)); error <- 0
  # The guide and quadrature rules depend on this replicate block, not its targets.
  mu <- colMeans(Y)
  guide <- list(mu = mu, sd = sqrt(colMeans(sweep(Y, 2L, mu)^2)))
  low.rule <- .np_udens_q6_rule(48L); high.rule <- .np_udens_q6_rule(64L)
  for (j in seq_along(reference$grid)) {
    low <- .np_udens_q6_curve(Y, reference$grid[j], reference$h, reference$kernel, 48L, guide, low.rule)
    high <- .np_udens_q6_curve(Y, reference$grid[j], reference$h, reference$kernel, 64L, guide, high.rule)
    error <- max(error, .np_udens_q6_check(low, high, "conditional curve"))
    root <- (high - reference$reference.center[j]) / reference$se[j]
    maximum <- pmax(maximum, abs(root))
  }
  list(maximum = maximum, refinement = error)
}

.np_udens_q6_intervals <- function(reference, maximum, alpha) {
  if (!is.numeric(maximum) || !length(maximum) || any(!is.finite(maximum)) ||
      any(maximum < 0) || !is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= .5)
    stop("invalid Q6 maximum calibration", call. = FALSE)
  critical <- unname(quantile(maximum, 1-alpha, type = 7))
  index <- reference$index
  center <- reference$center[index]; se <- reference$se[index]
  list(center = center, se = se, critical = critical,
       lower = center - critical*se, upper = center + critical*se)
}

# Admission is explicit and precedes data preparation, RNG and any dispatch.
.np_udens_q6_options <- function(dots, allow = FALSE) {
  value <- function(public, internal) {
    z <- dots[[public, exact = TRUE]]
    if (is.null(z)) dots[[internal, exact = TRUE]] else z
  }
  law <- value("bootstrap", "plot.errors.boot.method")
  if (!identical(law, "smoothed-q6")) return(invisible(FALSE))
  if (!isTRUE(allow)) stop("smoothed-q6 is available only for supported unconditional density plots", call. = FALSE)
  ctrl <- dots[["boot.control", exact = TRUE]]
  calibration <- if (!is.null(ctrl)) .np_plot_boot_calibration(ctrl) else
    dots[["plot.errors.boot.calibration", exact = TRUE]]
  if (!identical(value("errors", "plot.errors.method"), "bootstrap") ||
      !identical(value("center", "plot.errors.center"), "bias-corrected") ||
      !isTRUE(value("band", "plot.errors.type") %in% c("pointwise", "simultaneous", "bonferroni", "all")) ||
      !isTRUE(calibration %in% c("symmetric-sd", "asymmetric-sd")))
    stop("smoothed-q6 requires explicit bootstrap errors, bias-corrected center, pointwise, simultaneous, bonferroni or all band and symmetric-sd or asymmetric-sd calibration", call. = FALSE)
  if (!is.null(dots[["bias.control", exact = TRUE]]) ||
      isTRUE(dots[["gradients", exact = TRUE]]) || isTRUE(dots[["gradient", exact = TRUE]]) ||
      "gaussian.draws" %in% names(dots) ||
      !is.null(dots[["plot.errors.boot.blocklen", exact = TRUE]]) ||
      (!is.null(ctrl) && (!is.null(ctrl$blocklen) ||
        !identical(ctrl$nonfixed, "exact") || !identical(ctrl$wild, "mammen"))) ||
      (!is.null(dots[["plot.errors.boot.nonfixed", exact = TRUE]]) &&
        !identical(dots[["plot.errors.boot.nonfixed", exact = TRUE]], "exact")) ||
      (!is.null(dots[["plot.errors.boot.wild", exact = TRUE]]) &&
        !identical(dots[["plot.errors.boot.wild", exact = TRUE]], "mammen")))
    stop("smoothed-q6 does not accept custom bias, gradient, Gaussian-draw or dependent/nonfixed/wild controls", call. = FALSE)
  B <- value("B", "plot.errors.boot.num")
  if (!is.null(B) && (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B) || B > .Machine$integer.max))
    stop("Q6 B must be a positive representable integer", call. = FALSE)
  invisible(TRUE)
}

.np_udens_q6_spec <- function(bws) {
  if (!inherits(bws, "bandwidth") ||
      !identical(bws[["type", exact = TRUE]], "fixed") ||
      !identical(bws[["ckerbound", exact = TRUE]], "none") ||
      !isTRUE(bws[["ckertype", exact = TRUE]] %in% c("uniform", "epanechnikov")) ||
      !isTRUE(bws[["ckerorder", exact = TRUE]] == 2L) ||
      !isTRUE(bws[["ndim", exact = TRUE]] == 1L) ||
      !isTRUE(bws[["ncon", exact = TRUE]] == 1L))
    stop("smoothed-q6 requires q=1 continuous, whole-line, fixed-bandwidth order-two uniform or Epanechnikov density", call. = FALSE)
  invisible(NULL)
}

.np_udens_q6_prepare <- function(bws, x, seed) {
  .np_udens_q6_spec(bws)
  if (!is.data.frame(x) || ncol(x) != 1L || !nzchar(names(x)[1L]) ||
      is.object(x[[1L]]) || !is.numeric(x[[1L]]))
    stop("Q6 requires one ordinary numeric continuous coordinate", call. = FALSE)
  h <- as.double(.np_physical_bandwidth(bws))
  .np_udens_q6_geometry(x[[1L]], h, x[[1L]][1L], bws$ckertype)
  npValidateNonNegativeInteger(seed, "random.seed")
  list(x=x, h=h, n=nrow(x), q=1L, kernel=bws$ckertype, plot.random.seed=seed,
    plot.law="smoothed-q6", plot.calibration="symmetric-sd",
    construction="density-compact-q6-symmetric-v1")
}

.np_udens_q6_active <- function(state)
  identical(state[["construction", exact=TRUE]], "density-compact-q6-symmetric-v1")

.np_udens_q6_public_errors <- function(state, grid, alpha, B, band) {
  if (isTRUE(band %in% c("pointwise","bonferroni","all")))
    return(.np_udens_q6_pointwise_public(state, grid, alpha, B, band))
  if (!.np_udens_q6_active(state) || !identical(band,"simultaneous") ||
      !is.data.frame(grid) || ncol(grid)!=1L || !nrow(grid) ||
      !identical(names(grid), names(state$x)))
    stop("invalid Q6 simultaneous-band request", call.=FALSE)
  B <- .np_plot_gaussian_integer(B, "Q6 B", .Machine$integer.max)
  .np_udens_q6_geometry(state$x[[1L]], state$h, grid[[1L]], state$kernel)
  .np_udens_q6_budget(state$n, length(unique(grid[[1L]])), B)
  progress <- function(stage,done,total) .np_plot_progress_notify(
    stage=paste("Q6 density:",stage),done=done,total=total,unit="work")
  inference <- .np_plot_corrected_activity_run("Constructing Q6 density simultaneous intervals",
    .np_with_seed(state$plot.random.seed, {
      reference <- .np_udens_q6_reference(state$x[[1L]],state$h,grid[[1L]],state$kernel,progress)
      tasks <- .np_udens_q6_plan(get(".Random.seed",.GlobalEnv,inherits=FALSE),state$n,B)
      draws <- .np_udens_q6_dispatch(tasks,reference)
      if (!is.matrix(draws) || !identical(dim(draws),c(B,2L)) || any(!is.finite(draws)))
        stop("incomplete Q6 replicate collection",call.=FALSE)
      z <- .np_udens_q6_intervals(reference,draws[,1L],alpha)
      limits <- array(c(z$lower,z$upper),c(nrow(grid),2L,1L),
        dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
      bands <- list(simultaneous=limits)
      resolution <- .np_band_resolution(B,rep("density",nrow(grid)),alpha,
        rep(TRUE,nrow(grid)),"simultaneous")
      resolution$calibration <- "smoothed-q6 symmetric-sd"
      resolution$quantile.type <- 7L
      resolution$interpolated.index <- 1+(B-1)*(1-alpha)
      list(center=z$center,se=z$se,critical=z$critical,maximum=draws[,1L],
        bands=bands,intervals=bands,absolute.intervals=bands,covariance=NULL,
        resolution=resolution,
        checks=c(reference$diagnostics,conditional=max(draws[,2L])),
        settings=list(construction=state$construction,law="smoothed-q6",h=state$h,n=state$n,q=1L,
          alpha=alpha,B=B,grid.size=length(reference$grid),family.rows=reference$index,
          family="unique supplied finite grid",available.bands="simultaneous",primary="absolute",
          quantile.type=7L,calibration="symmetric-sd",
          maximum.statistic="complete corrected bootstrap root divided by fixed original scale",
          reference="empirical mixture plus h*sqrt(2/3)*Student6 noise",
          root.center="full corrected functional T_h(Q6)",
          covariance="complete corrected-curve influence under Q6",
          correction="order-six compact operator with refitted normal guide",
          calibration.status="explicit simultaneous-only Q6 density plot inference"))
    }))
  .np_band_resolution_warn(inference$resolution)
  width <- inference$critical*inference$se
  list(boot.err=cbind(width,width,inference$center),boot.all.err=NULL,inference=inference,bxp=list())
}
