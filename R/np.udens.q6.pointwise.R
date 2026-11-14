# Marginal Q6 calibration shares the conditional law with simultaneous inference.
# Target blocks bound return/quantile storage; primitive streams are replayed exactly.
.np_udens_q6_pointwise_budget <- function(n, m, B, targets = 8L) {
  base <- .np_udens_q6_budget(n, m, B)
  bytes <- base + 8 * (6 * as.double(B) * (min(m, targets) + 1) + 4 * m)
  if (!is.finite(bytes) || bytes > 192 * 1024^2)
    stop("Q6 pointwise inference exceeds its 192 MiB workspace bound", call. = FALSE)
  bytes
}

.np_udens_q6_pointwise_block <- function(n, m, B) {
  # Preserve the original admission boundary, then use the remaining envelope.
  .np_udens_q6_pointwise_budget(n, m, B)
  if (B < 1L) stop("invalid Q6 pointwise replicate count", call. = FALSE)
  base <- .np_udens_q6_budget(n, m, B)
  targets <- min(m, floor((192 * 1024^2 - base - 32 * m) / (48 * B)) - 1)
  if (targets < 1L)
    stop("Q6 pointwise inference exceeds its 192 MiB workspace bound", call. = FALSE)
  .np_udens_q6_pointwise_budget(n, m, B, targets)
  as.integer(targets)
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

.np_udens_q6_pointwise_reduce <- function(reference, tasks, alpha, companions=FALSE) {
  B <- sum(vapply(tasks, function(z) z$bsz, 0L))
  block <- .np_udens_q6_pointwise_block(length(reference$x), length(reference$grid), B)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= .5)
    stop("invalid Q6 pointwise alpha", call. = FALSE)
  critical <- numeric(length(reference$grid)); error <- 0
  if(companions) { bonferroni <- critical; maximum <- numeric(B) }
  for (first in seq.int(1L, length(reference$grid), by = block)) {
    ix <- first:min(first + block - 1L, length(reference$grid))
    part <- reference
    for (field in c("grid", "center", "reference.center", "se")) part[[field]] <- reference[[field]][ix]
    part$covariance <- NULL
    draws <- .np_udens_q6_pointwise_dispatch(tasks, part)
    if (!is.matrix(draws) || !identical(dim(draws), c(B, length(ix) + 1L)) ||
        any(!is.finite(draws)) || any(draws < 0))
      stop("incomplete Q6 pointwise replicate collection", call. = FALSE)
    critical[ix] <- vapply(seq_along(ix), function(j)
      unname(quantile(draws[, j], 1-alpha, type = 7)), 0.0)
    if(companions) {
      bonferroni[ix] <- vapply(seq_along(ix),function(j)
        unname(quantile(draws[,j],1-alpha/length(reference$grid),type=7)),0.0)
      maximum <- pmax(maximum,apply(draws[,seq_along(ix),drop=FALSE],1L,max))
    }
    error <- max(error, draws[, ncol(draws)])
    .np_plot_progress_notify(stage = "Calibrating Q6 pointwise bands",
      done = max(ix), total = length(reference$grid), unit = "targets")
  }
  index <- reference$index
  out <- list(center = reference$center[index], se = reference$se[index],
       critical = critical[index], refinement = error)
  if(companions) {
    out$maximum <- maximum
    out$simultaneous <- unname(quantile(maximum,1-alpha,type=7))
    out$bonferroni.raw <- bonferroni[index]
    out$bonferroni <- pmax(out$bonferroni.raw,out$simultaneous)
  }
  out
}

.np_udens_q6_pointwise_public <- function(state, grid, alpha, B, band="pointwise") {
  if (!.np_udens_q6_active(state) || !isTRUE(band %in% c("pointwise","bonferroni","all")) || !is.data.frame(grid) || ncol(grid) != 1L ||
      !nrow(grid) || !identical(names(grid), names(state$x)))
    stop("invalid Q6 pointwise request", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "Q6 B", .Machine$integer.max)
  .np_udens_q6_geometry(state$x[[1L]], state$h, grid[[1L]], state$kernel)
  .np_udens_q6_pointwise_budget(state$n, length(unique(grid[[1L]])), B)
  progress <- function(stage, done, total) .np_plot_progress_notify(
    stage = paste("Q6 density:", stage), done = done, total = total, unit = "work")
  companions <- band %in% c("bonferroni","all")
  inference <- .np_plot_corrected_activity_run("Constructing Q6 density marginal intervals",
    .np_with_seed(state$plot.random.seed, {
      reference <- .np_udens_q6_reference(state$x[[1L]], state$h, grid[[1L]], state$kernel, progress)
      tasks <- .np_udens_q6_plan(get(".Random.seed", .GlobalEnv, inherits = FALSE), state$n, B)
      z <- .np_udens_q6_pointwise_reduce(reference, tasks, alpha, companions)
      make <- function(critical) {
        width <- critical*z$se
        array(c(z$center-width,z$center+width),c(nrow(grid),2L,1L),
          dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
      }
      bands <- list(pointwise=make(z$critical))
      if(companions) {
        bands$simultaneous <- make(z$simultaneous)
        bands$bonferroni <- make(z$bonferroni)
      }
      requested <- if(band=="all")names(bands) else band
      bands <- bands[requested]
      family.size <- if(companions)length(reference$grid) else nrow(grid)
      resolution <- .np_band_resolution(B,rep("density",family.size),alpha,
        rep(TRUE,family.size),requested)
      resolution$calibration <- "smoothed-q6 symmetric-sd"
      resolution$quantile.type <- 7L
      resolution$interpolated.index <- 1+(B-1)*(1-alpha/ifelse(resolution$band=="bonferroni",family.size,1L))
      list(center=z$center, se=z$se, critical=z$critical, bands=bands,
        intervals=bands, absolute.intervals=bands, covariance=NULL, resolution=resolution,
        maximum=if(companions)z$maximum else NULL,
        raw.bonferroni=if(companions)make(z$bonferroni.raw) else NULL,
        companion.critical=if(companions)list(pointwise=z$critical,simultaneous=z$simultaneous,
          bonferroni=z$bonferroni,bonferroni.raw=z$bonferroni.raw) else NULL,
        checks=c(reference$diagnostics, conditional=z$refinement),
        settings=list(construction=state$construction, law="smoothed-q6", h=state$h,
          n=state$n, q=1L, alpha=alpha, B=B, grid.size=length(reference$grid),
          family.rows=reference$index, family=if(companions)"unique supplied finite grid" else "marginal intervals at supplied fixed targets",
          available.bands=requested, primary="absolute", quantile.type=7L,
          bonferroni.construction=if(companions)"type-7 marginal absolute quantile at alpha/m with simultaneous hull; raw endpoints retained" else NULL,
          calibration="symmetric-sd",
          marginal.statistic="absolute complete corrected bootstrap root divided by fixed original scale",
          reference="empirical mixture plus h*sqrt(2/3)*Student6 noise",
          root.center="full corrected functional T_h(Q6)",
          covariance="complete corrected-curve influence under Q6",
          correction="order-six compact operator with refitted normal guide",
          calibration.status="explicit symmetric Q6 density plot inference"))
    }))
  .np_band_resolution_warn(inference$resolution)
  widths <- if(companions) lapply(inference$companion.critical[names(inference$bands)],function(k){width<-k*inference$se;cbind(width,width)}) else NULL
  if(!companions) {
    # Preserve the incumbent directly multiplied offsets exactly.
    width <- inference$critical*inference$se
    return(list(boot.err=cbind(width,width,inference$center),boot.all.err=NULL,inference=inference,bxp=list()))
  }
  list(boot.err=cbind(widths[[if(band=="all")"pointwise" else band]],inference$center),
       boot.all.err=if(band=="all")widths else NULL,inference=inference,bxp=list())
}
