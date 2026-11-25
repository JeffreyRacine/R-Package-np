# Explicit experimental conditional inference; incumbent NULL controls are unchanged.
np_conditional_bias_control <- function(construction = c("analytic", "jackknife", "jackknife-mp1", "pairs-score"),
                                        pointwise = NULL) {
  construction <- match.arg(construction)
  if (is.null(pointwise))
    return(structure(list(construction = construction, version = 1L),
                     class = "np_conditional_bias_control"))
  if (!identical(construction, "pairs-score") || !identical(pointwise, "local-wilson"))
    stop("pointwise = 'local-wilson' requires construction = 'pairs-score'", call. = FALSE)
  structure(list(construction = construction, pointwise = pointwise, version = 2L),
            class = "np_conditional_bias_control")
}

.np_conditional_control <- function(control, bws, cdf, method, center,
                                    calibration, law, gradients, quantreg, proper) {
  if (is.null(control)) {
    if (identical(law, "conditional-pilot"))
      stop('bootstrap="conditional-pilot" requires np_conditional_bias_control()', call. = FALSE)
    return(FALSE)
  }
  pointwise <- if (is.list(control)) control[["pointwise", exact = TRUE]]
  if (!inherits(control, "np_conditional_bias_control") ||
      !is.list(control) || !identical(control[["construction", exact = TRUE]] %in%
        c("analytic", "jackknife", "jackknife-mp1", "pairs-score"), TRUE) ||
      (!is.null(pointwise) && !(identical(pointwise, "local-wilson") &&
        identical(control[["construction", exact = TRUE]], "pairs-score"))) ||
      !identical(control, np_conditional_bias_control(control[["construction", exact = TRUE]], pointwise)))
    stop("invalid np_conditional_bias_control object", call. = FALSE)
  if (!identical(center, "bias-corrected") ||
      !method %in% c("bootstrap", "asymptotic"))
    stop('conditional bias.control requires center="bias-corrected" and errors="asymptotic" or "bootstrap"', call. = FALSE)
  compact.response <- identical(control$construction, "pairs-score") && bws$cykertype %in% c("epanechnikov", "uniform")
  beta.response <- identical(control$construction, "pairs-score") && identical(bws$cykertype, "beta") &&
    identical(bws$cykerbound, "fixed") && length(bws$cykerlb) == 1L && length(bws$cykerub) == 1L &&
    is.finite(bws$cykerlb) && is.finite(bws$cykerub) && bws$cykerlb < bws$cykerub
  spec <- npConditionalRegEngineSpec(bws, where = "conditional bias inference")
  if (isTRUE(gradients) || isTRUE(quantreg) || isTRUE(proper) ||
      !identical(spec$reg.engine, "lc") || bws$xndim != 1L || bws$yndim != 1L ||
      bws$xncon != 1L || bws$yncon != 1L || !identical(bws$type, "fixed") ||
      !identical(bws$cxkertype, "gaussian") || (!identical(bws$cykertype, "gaussian") && !compact.response && !beta.response) ||
      bws$cxkerorder != 2L || bws$cykerorder != 2L ||
      any(is.finite(c(bws$cxkerlb, bws$cxkerub))) ||
      (!beta.response && any(is.finite(c(bws$cykerlb, bws$cykerub)))))
    stop(if (identical(control$construction, "pairs-score"))
      "conditional CDF pairs-score requires LC, one continuous X and Y, fixed unbounded Gaussian2 X and Gaussian2, Epanechnikov2, uniform or known-support beta2 Y, no gradients/quantiles/proper projection"
    else "conditional bias inference requires LC, one continuous X and Y, fixed unbounded Gaussian2 kernels, no gradients/quantiles/proper projection", call. = FALSE)
  if (identical(control$construction, "pairs-score")) {
    if (!cdf || method != "bootstrap" || !identical(law, "inid") ||
        !identical(calibration, "default"))
      stop('pairs-score requires conditional CDF, bootstrap="inid" and percentile intervals', call. = FALSE)
    return(TRUE)
  }
  if (identical(control$construction,"jackknife-mp1")) {
    if (cdf || method!="bootstrap" || !identical(law,"conditional-pilot") ||
        !identical(calibration,"asymmetric-sd"))
      stop('jackknife-mp1 requires density bootstrap, bootstrap="conditional-pilot" and interval="asymmetric-sd"',call.=FALSE)
    return(TRUE)
  }
  construction <- if (!cdf && method == "asymptotic") "jackknife" else "analytic"
  if (!identical(control$construction, construction))
    stop(sprintf('this conditional inference route requires construction="%s"', construction), call. = FALSE)
  if (method == "bootstrap" && (!identical(law, "conditional-pilot") ||
      !calibration %in% c("symmetric-sd", "asymmetric-sd")))
    stop('conditional bootstrap requires explicit bootstrap="conditional-pilot" and interval="symmetric-sd" or interval="asymmetric-sd"', call. = FALSE)
  TRUE
}

# Analytic moments grow quadratically in the number of response targets.
# Pairs draws use only donor weights, numerators and the returned root matrix.
.np_conditional_grid_admit <- function(m, ny, pairs.score) {
  cap <- if (pairs.score) c(2500L, 50L) else c(256L, 25L)
  if (m < 1L || m > cap[1L] || ny > cap[2L])
    stop(sprintf("experimental conditional inference supports at most %d target points and %d distinct response targets; reduce neval explicitly",
                 cap[1L], cap[2L]), call. = FALSE)
}

.np_conditional_workspace <- function(x, y, grid, h, B, pairs.score = FALSE) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y) || length(x) < 3L ||
      !is.data.frame(grid) || !identical(names(grid), c("x", "y")) ||
      !is.numeric(grid$x) || !is.numeric(grid$y) || length(h) != 2L ||
      any(!is.finite(c(x, y, unlist(grid), h))) || any(h <= 0))
    stop("invalid fixed conditional inference inputs", call. = FALSE)
  B <- .np_plot_gaussian_integer(B, "conditional inference draws")
  if (B < 2L) stop("conditional inference requires at least 2 draws", call.=FALSE)
  n <- length(x); m <- nrow(grid); ny <- length(unique(grid$y)); nx <- length(unique(grid$x))
  .np_conditional_grid_admit(m, ny, pairs.score)
  limit <- 256 * 1024^2
  within.old.grid <- m <= 256L && ny <= 25L
  bytes <- if (within.old.grid)
    8 * (4*n*(3*ny+1)^2 + 12*n*(nx+ny) +
      2*n*B + 8*m*B + 12*m*m + 12*n*128 + 256*n)
  else Inf
  if (pairs.score && bytes > limit) {
    # The pairs provider builds donor and draw matrices, not analytic moments.
    # Preserve the reported budget on old admissions; use the same pairs bound
    # for any other candidate so shrinking a successful grid cannot refuse it.
    bytes <- 8 * (24 * n * m + 3 * n * B + 12 * m * B + 12 * n * 128 + 256 * n)
  }
  if (pairs.score && bytes > limit) {
    if (B > 9999L)
      stop("conditional pairs numerical workspace exceeds the 256 MiB budget and tiling supports at most 9,999 draws; reduce the grid or draws explicitly", call. = FALSE)
    allocation <- .np_conditional_pairs_allocation(n, m, B)
    if (!is.null(allocation)) return(allocation)
  }
  if (!is.finite(bytes) || bytes > limit)
    stop("conditional inference numerical workspace exceeds the 256 MiB budget; reduce the requested grid or draws explicitly", call. = FALSE)
  list(B = B, estimated.bytes = bytes)
}

# Same donor order and logical blocks as the reference; no worker RNG needed.
.np_conditional_tasks <- function(z, B, gaussian = FALSE) {
  starts <- seq.int(1L, B, by = 128L)
  tasks <- vector("list", length(starts))
  for (a in seq_along(starts)) {
    first <- starts[a]; size <- min(128L, B-first+1L)
    if (gaussian) value <- matrix(stats::rnorm(z$m*size), z$m, size) else {
      donor <- matrix(0L, size, z$n)
      for (i in seq.int(1L, z$n, by = 32L)) {
        ii <- i:min(z$n, i+31L); p <- .np_conditional_den_probabilities(z, ii)
        for (j in seq_along(ii)) donor[, ii[j]] <- sample.int(z$n, size, TRUE, p[j, ])
      }
      value <- matrix(z$y[t(donor)], z$n, size) + sqrt(2)*z$h[2]*matrix(stats::rnorm(z$n*size), z$n, size)
    }
    tasks[[a]] <- list(start = first, bsz = size, value = value)
  }
  tasks
}

.np_conditional_root_block <- function(task, plan, family) {
  if (family == "jackknife-mp1") return(.np_conditional_mp1_block(task,plan))
  if (family == "pairs-score") return(.np_conditional_pairs_block(task, plan))
  if (family == "gaussian") return(t(plan$factor %*% task$value))
  if (family == "cdf") return(t(.np_conditional_cp_apply(plan, task$value)))
  z <- plan$z; Y <- task$value; out <- matrix(0, z$m, task$bsz)
  for (j in seq_len(z$ny)) {
    v <- (z$gy[j]-Y)/z$h[2]; K <- stats::dnorm(v)/z$h[2]
    R <- crossprod(z$A, K-plan$mean[,j]) + crossprod(z$D, v*v*K-plan$mean[,z$ny+j])
    ii <- which(z$iy == j); out[ii, ] <- R[z$ix[ii], , drop=FALSE]
  }
  t(out)
}

.np_conditional_dispatch <- function(tasks, plan, family, m) {
  out <- matrix(0, sum(vapply(tasks, function(x) x$bsz, 0L)), m)
  label <- plan[["tile.label"]]
  if (is.null(label)) label <- "Conditional inference"
  progress <- .np_plot_bootstrap_progress_begin(nrow(out), label=label, unit="rep")
  on.exit(.np_plot_progress_end(progress), add=TRUE)
  for (task in tasks) {
    ii <- seq.int(task$start, length.out=task$bsz)
    out[ii, ] <- .np_conditional_root_block(task, plan, family)
    progress <- .np_plot_progress_tick(progress, max(ii))
  }
  t(out)
}

.np_conditional_gaussian <- function(center, C, group, B, alpha, z, cdf) {
  s <- sqrt(diag(C))
  if (any(!is.finite(s)) || any(s <= 0)) stop("nonpositive conditional pilot variance", call.=FALSE)
  scale <- unname(vapply(split(s, group), max, 0.0)[as.character(group)])
  e <- eigen(C, symmetric=TRUE)
  if (min(e$values) < -max(e$values)*1e-10) stop("conditional pilot covariance is not positive semidefinite", call.=FALSE)
  F <- sweep(e$vectors, 2L, sqrt(pmax(e$values, 0)), `*`)
  tasks <- .np_conditional_tasks(z, B, gaussian=TRUE)
  roots <- .np_conditional_dispatch(tasks, list(factor=F), "gaussian", z$m)
  critical <- max(stats::qnorm(1-alpha/2), as.numeric(stats::quantile(apply(abs(roots)/scale, 2L, max), 1-alpha, type=1L)))
  make <- function(w) cbind(lower=center-w, upper=center+w)
  bands <- list(pointwise=make(stats::qnorm(1-alpha/2)*s), simultaneous=make(critical*scale),
    bonferroni=make(if(cdf) max(critical,stats::qnorm(1-alpha/(2*z$m)))*scale else pmax(stats::qnorm(1-alpha/(2*z$m))*s,critical*scale)))
  list(center=center, se=s, scale=scale, critical=critical, bands=bands)
}

.np_conditional_inference <- function(x, y, grid, h, cdf, method, calibration,
                                      B=9999L, alpha=.05, cdf.centers=NULL, construction=NULL, response.kernel="gaussian", pointwise=NULL, response.bounds=NULL) {
  allocation <- .np_conditional_workspace(x, y, grid, h, B, identical(construction, "pairs-score"))
  B <- allocation$B
  if (!is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=.5)
    stop("alpha must lie in (0, 0.5)", call.=FALSE)
  if (identical(construction, "pairs-score")) {
    if (!isTRUE(cdf) || method != "bootstrap" || !identical(calibration, "default") ||
        !is.matrix(cdf.centers) || !identical(dim(cdf.centers), c(nrow(grid), if (response.kernel %in% c("epanechnikov", "uniform", "beta")) 1L else 2L)) ||
        any(!is.finite(cdf.centers[, 1L])))
      stop("invalid conditional CDF pairs-score inputs", call. = FALSE)
    .np_conditional_runtime_admit(length(x), B)
    return(.np_conditional_pairs_score(x, y, grid, h, B, alpha, cdf.centers[, 1L],
                                       allocation$estimated.bytes, allocation[["tile.size"]], response.kernel, pointwise, response.bounds))
  }
  if (!is.logical(cdf) || length(cdf)!=1L || is.na(cdf) ||
      !method %in% c("asymptotic", "bootstrap") ||
      (method=="bootstrap" && !calibration %in% c("symmetric-sd", "asymmetric-sd")))
    stop("invalid conditional inference method", call.=FALSE)
  if (identical(construction,"jackknife-mp1")) {
    if (cdf || method!="bootstrap" || !identical(calibration,"asymmetric-sd"))
      stop("invalid MP1 inference route",call.=FALSE)
    .np_conditional_mp1_workspace(x,y,grid,h,B)
    .np_conditional_runtime_admit(length(x),B)
    return(.np_conditional_mp1_inference(x,y,grid,h,B,alpha))
  }
  .np_conditional_runtime_admit(length(x),B)
  if (cdf) {
    z <- .np_conditional_context(x,y,grid,h)
    if (!is.matrix(cdf.centers) || !identical(dim(cdf.centers), c(nrow(grid),2L)) ||
        any(!is.finite(cdf.centers))) stop("CDF inference requires native-compatible centers", call.=FALSE)
    center <- cdf.centers[,2L]; original <- cdf.centers[,1L]
    bounds <- .np_conditional_weighted_range(x,grid,h[1L])
    if (method=="asymptotic") {
      mm <- .np_conditional_moments(z,sqrt(2),96L)
      covariance <- .np_conditional_conditional_cov(z,mm,"analytic")
      base <- .np_conditional_gaussian(center,covariance,grid$x,B,alpha,z,TRUE)
      out <- .np_conditional_range_envelope(base,grid$x,bounds,"sharp-score",alpha)
    } else {
      plan <- .np_conditional_cp_plan(z,B,256*1024^2)
      roots <- .np_conditional_dispatch(.np_conditional_tasks(z,B),plan,"cdf",z$m)
      base <- .np_conditional_bootstrap_slice(center,roots,grid$x,calibration,alpha)
      out <- .np_conditional_bootstrap_score(base,grid$x,bounds,calibration,alpha)
      covariance <- NULL
    }
  } else {
    z <- .np_conditional_den_context(x,y,grid,h)
    ctr <- .np_conditional_den_center(z); original <- ctr$original
    if (method=="asymptotic") {
      z <- .np_conditional_jack_context(z); center <- .np_conditional_jack_center(z)
      p <- .np_conditional_jack_pilot(z); covariance <- p$covariance
      base <- .np_conditional_gaussian(center,covariance,seq_len(z$m),B,alpha,z,FALSE)
      out <- .np_conditional_marginal_guard(base,.np_conditional_jack_range(z),
        list(range_interval=.np_conditional_range_interval),alpha)
    } else {
      center <- ctr$corrected; p <- .np_conditional_den_pilot(z)
      roots <- .np_conditional_dispatch(.np_conditional_tasks(z,B),list(z=z,mean=p$mean),"density",z$m)
      base <- .np_conditional_bootstrap_slice(center,roots,grid$x,calibration,alpha)
      out <- .np_conditional_bootstrap_score(base,grid$x,.np_conditional_den_range(z),calibration,alpha)
      covariance <- NULL
    }
  }
  out$original <- original; out$se <- base$se; out$scale <- base$scale
  out$covariance <- covariance
  out$covariance.provenance <- if(method=="asymptotic") "complete conditional-pilot moments, including correction cross terms" else "complete corrected conditional-pilot roots; empirical marginal scales and family extremes"
  out$alpha <- alpha; out$B <- B
  out$law <- "conditional-pilot"; out$correction <- if(!cdf && method=="asymptotic") "jackknife" else "analytic"
  out$calibration <- if(method=="asymptotic") "gaussian-score" else calibration
  out$physical.bandwidth <- h; out$grid <- grid
  out$scope <- "experimental LC; one continuous X and Y; fixed Gaussian2; requested finite target family"
  out$estimated.numerical.bytes <- allocation$estimated.bytes
  out$resolution.limited <- method=="bootstrap" && B*alpha/(2*z$m)<1
  out
}

# Native-compatible original/analytic CDF centers. The fixed-X covariance and
# response law retain their frozen Gaussian moment owners. In particular this
# must not replace the native integral by pnorm: the native approximation depends
# on its computational route, and the retained coverage cases used that center.
.np_conditional_cdf_centers <- function(bws, x, y, ex, ey) {
  hx <- .npcdhat_physical_bandwidth(bws,"x"); hy <- .npcdhat_physical_bandwidth(bws,"y")
  if (npLogicalOption("np.largeh",TRUE)) {
    tol <- getOption("np.largeh.rel.tol",1e-3)
    if (!is.numeric(tol) || length(tol)!=1L || !is.finite(tol) || tol<=0 || tol>=.1)
      stop("invalid np.largeh.rel.tol",call.=FALSE)
    t <- sqrt(-2*log(1-tol))
    ty <- switch(bws$cykertype, uniform = 1, epanechnikov = sqrt(tol), t)
    # The new response-kernel route uses both ordinary and sqrt(2)-wide ratios.
    enlargement <- if (bws$cykertype %in% c("epanechnikov", "uniform", "beta")) sqrt(2) else 1
    if(diff(range(x[[1L]],ex[[1L]]))<=enlargement*hx*t ||
       (!identical(bws$cykertype, "beta") && diff(range(y[[1L]],ey[[1L]]))<=enlargement*hy*ty))
      stop("experimental conditional CDF inference does not admit the near-constant kernel substitution",call.=FALSE)
  }
  K <- .np_conditional_literal_density_matrix(.npcdhat_make_xkbw(bws,x),x,ex,"normal")
  G <- .np_conditional_literal_density_matrix(.npcdhat_make_ybw(bws,y),y,ey,"integral")
  if (bws$cykertype %in% c("epanechnikov", "uniform", "beta")) {
    mass <- rowSums(K)
    if (any(!is.finite(mass)) || any(mass <= 0))
      stop("conditional CDF denominator is undefined", call. = FALSE)
    return(cbind(original = rowSums(K * G) / mass))
  }
  u <- outer(ex[[1L]],x[[1L]],"-")/hx; v <- outer(ey[[1L]],y[[1L]],"-")/hy
  a <- (u^2-1)/2; R <- -v*stats::dnorm(v)/2
  g0 <- G[,1L]; r0 <- R[,1L]; g <- G-g0; rr <- R-r0; ds <- rowSums(K)
  if(any(!is.finite(ds)) || any(ds<=0)) stop("conditional CDF denominator is undefined",call.=FALSE)
  theta <- rowSums(K*g)/ds; ar <- rowSums(K*(a*g+rr))/ds; cr <- rowSums(K*a)/ds
  cbind(original=g0+theta,analytic=g0-r0+theta-ar+theta*cr)
}

.np_conditional_plot_interval <- function(bws, x, y, ex, ey, cdf, method,
                                          calibration, B, alpha, band, construction=NULL, pointwise=NULL) {
  .np_plot_corrected_activity_run("Preparing conditional inference", {
  ex <- toFrame(ex); ey <- toFrame(ey)
  grid <- data.frame(x=ex[[1L]], y=ey[[1L]])
  h <- c(.npcdhat_physical_bandwidth(bws,"x"), .npcdhat_physical_bandwidth(bws,"y"))
  .np_conditional_workspace(x[[1L]], y[[1L]], grid, h, B, identical(construction, "pairs-score"))
  cc <- if(cdf) .np_conditional_cdf_centers(bws,x,y,ex,ey) else NULL
  a <- .np_conditional_inference(x[[1L]],y[[1L]],grid,h,cdf,method,calibration,B,alpha,cc,construction,response.kernel=bws$cykertype,pointwise=pointwise,
                                 response.bounds=if (identical(bws$cykertype, "beta")) c(bws$cykerlb, bws$cykerub))
  if(isTRUE(a$resolution.limited)) .np_warning(paste(a$law, "Bonferroni tails have limited finite-draw resolution; increase B for better tail resolution"))
  widths <- lapply(a$bands,function(v)cbind(lower=a$center-v[,1L],upper=v[,2L]-a$center,center=a$center))
  chosen <- if(band=="all") "simultaneous" else band
  list(err=widths[[chosen]], all.err=widths, inference=a)
  })
}

.np_conditional_plot_fields <- function(object, interval, band) {
  a <- interval$inference; object$inference <- a
  object$bias.corrected <- a$center; object$bias <- a$original-a$center
  chosen <- if(band=="all") "simultaneous" else band
  object$conderr <- sweep(a$bands[[chosen]],1L,a$center,`-`)
  object
}

.np_conditional_runtime_admit <- function(n,B) invisible(NULL)
