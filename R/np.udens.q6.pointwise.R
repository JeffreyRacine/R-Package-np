# Marginal Q6 calibration shares the conditional law with simultaneous inference.
# Target blocks bound return/quantile storage; primitive streams are replayed exactly.
.np_udens_q6_pointwise_budget <- function(n, m, B, targets = 8L) {
  base <- .np_udens_q6_budget(n, m, B)
  bytes <- base + 8 * (6 * as.double(B) * (min(m, targets) + 1) + 4 * m)
  if (!is.finite(bytes) || bytes > 192 * 1024^2)
    stop("Q6 pointwise inference exceeds its 192 MiB workspace bound", call. = FALSE)
  bytes
}

.np_udens_q6_pointwise_chunk <- function(task, reference) {
  Y <- .np_udens_q6_samples(task, reference$x, reference$h)
  mu <- colMeans(Y)
  guide <- list(mu = mu, sd = sqrt(colMeans(sweep(Y, 2L, mu)^2)))
  low.rule <- .np_udens_q6_rule(48L); high.rule <- .np_udens_q6_rule(64L)
  out <- matrix(0, task$bsz, length(reference$grid) + 1L)
  error <- 0
  for (j in seq_along(reference$grid)) {
    low <- .np_udens_q6_curve(Y, reference$grid[j], reference$h, reference$kernel, 48L, guide, low.rule)
    high <- .np_udens_q6_curve(Y, reference$grid[j], reference$h, reference$kernel, 64L, guide, high.rule)
    error <- max(error, .np_udens_q6_check(low, high, "conditional curve"))
    out[, j] <- abs((high - reference$reference.center[j]) / reference$se[j])
  }
  out[, ncol(out)] <- error
  out
}

.np_udens_q6_pointwise_reduce <- function(reference, tasks, alpha) {
  B <- sum(vapply(tasks, function(z) z$bsz, 0L))
  .np_udens_q6_pointwise_budget(length(reference$x), length(reference$grid), B)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= .5)
    stop("invalid Q6 pointwise alpha", call. = FALSE)
  critical <- numeric(length(reference$grid)); error <- 0
  for (first in seq.int(1L, length(reference$grid), by = 8L)) {
    ix <- first:min(first + 7L, length(reference$grid))
    part <- reference
    for (field in c("grid", "center", "reference.center", "se")) part[[field]] <- reference[[field]][ix]
    part$covariance <- NULL
    draws <- .np_udens_q6_pointwise_dispatch(tasks, part)
    if (!is.matrix(draws) || !identical(dim(draws), c(B, length(ix) + 1L)) ||
        any(!is.finite(draws)) || any(draws < 0))
      stop("incomplete Q6 pointwise replicate collection", call. = FALSE)
    critical[ix] <- vapply(seq_along(ix), function(j)
      unname(quantile(draws[, j], 1-alpha, type = 7)), 0.0)
    error <- max(error, draws[, ncol(draws)])
    .np_plot_progress_notify(stage = "Calibrating Q6 pointwise bands",
      done = max(ix), total = length(reference$grid), unit = "targets")
  }
  index <- reference$index
  list(center = reference$center[index], se = reference$se[index],
       critical = critical[index], refinement = error)
}

.np_udens_q6_pointwise_public <- function(state, grid, alpha, B) {
  if (!.np_udens_q6_active(state) || !is.data.frame(grid) || ncol(grid) != 1L ||
      !nrow(grid) || !identical(names(grid), names(state$x)))
    stop("invalid Q6 pointwise request", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "Q6 B", .Machine$integer.max)
  .np_udens_q6_geometry(state$x[[1L]], state$h, grid[[1L]], state$kernel)
  .np_udens_q6_pointwise_budget(state$n, length(unique(grid[[1L]])), B)
  progress <- function(stage, done, total) .np_plot_progress_notify(
    stage = paste("Q6 density:", stage), done = done, total = total, unit = "work")
  inference <- .np_plot_corrected_activity_run("Constructing Q6 density pointwise intervals",
    .np_with_seed(state$plot.random.seed, {
      reference <- .np_udens_q6_reference(state$x[[1L]], state$h, grid[[1L]], state$kernel, progress)
      tasks <- .np_udens_q6_plan(get(".Random.seed", .GlobalEnv, inherits = FALSE), state$n, B)
      z <- .np_udens_q6_pointwise_reduce(reference, tasks, alpha)
      width <- z$critical * z$se
      limits <- array(c(z$center-width, z$center+width), c(nrow(grid), 2L, 1L),
        dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
      bands <- list(pointwise = limits)
      resolution <- .np_band_resolution(B, rep("density", nrow(grid)), alpha,
        rep(TRUE, nrow(grid)), "pointwise")
      resolution$calibration <- "smoothed-q6 symmetric-sd"
      resolution$quantile.type <- 7L
      resolution$interpolated.index <- 1+(B-1)*(1-alpha)
      list(center=z$center, se=z$se, critical=z$critical, bands=bands,
        intervals=bands, absolute.intervals=bands, covariance=NULL, resolution=resolution,
        checks=c(reference$diagnostics, conditional=z$refinement),
        settings=list(construction=state$construction, law="smoothed-q6", h=state$h,
          n=state$n, q=1L, alpha=alpha, B=B, grid.size=length(reference$grid),
          family.rows=reference$index, family="marginal intervals at supplied fixed targets",
          available.bands="pointwise", primary="absolute", quantile.type=7L,
          calibration="symmetric-sd",
          marginal.statistic="absolute complete corrected bootstrap root divided by fixed original scale",
          reference="empirical mixture plus h*sqrt(2/3)*Student6 noise",
          root.center="full corrected functional T_h(Q6)",
          covariance="complete corrected-curve influence under Q6",
          correction="order-six compact operator with refitted normal guide",
          calibration.status="explicit pointwise Q6 density plot inference"))
    }))
  .np_band_resolution_warn(inference$resolution)
  width <- inference$critical * inference$se
  list(boot.err=cbind(width, width, inference$center), boot.all.err=NULL, inference=inference, bxp=list())
}
