# Normal companions share the complete covariance owner of the corrected curve.
# Bonferroni uses the unique supplied finite grid, without a simultaneous hull.
.np_udens_whole_normal_assemble <- function(inference, alpha, band) {
  m <- inference$settings$grid.size
  make <- function(k) {
    width <- k*inference$se
    array(c(inference$center-width,inference$center+width),
      c(length(inference$center),2L,1L),
      dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
  }
  normal <- list(pointwise=make(qnorm(1-alpha/2)),
    bonferroni=make(qnorm(1-alpha/(2*m))))
  requested <- if (identical(band,"all")) c("pointwise","simultaneous","bonferroni") else band
  bands <- c(normal,inference$bands)[requested]
  if (any(vapply(bands,is.null,logical(1L))))
    stop("unavailable compact normal component",call.=FALSE)
  inference$bands <- inference$intervals <- inference$absolute.intervals <- bands
  inference$settings$available.bands <- requested
  inference$settings$normal.calibration <- "two-sided normal quantile"
  inference$settings$bonferroni.construction <- "raw normal quantile at alpha/m; no simultaneous hull"
  inference$settings$calibration.status <- "normal companions of the complete corrected-curve covariance"
  widths <- lapply(bands,function(z)cbind(inference$center-z[,1L,1L],
    z[,2L,1L]-inference$center))
  list(boot.err=cbind(widths[[if (band=="all") "pointwise" else band]],inference$center),
    boot.all.err=if (band=="all") widths else NULL,inference=inference,bxp=list())
}

.np_udens_m4_normal_public <- function(state, grid, alpha, B, band) {
  if (!.np_udens_m4_active(state) || !isTRUE(band %in% c("pointwise","bonferroni","all")) ||
      !is.data.frame(grid) || ncol(grid)!=1L || !nrow(grid) ||
      !identical(names(grid),names(state$x)) || !is.numeric(grid[[1L]]) ||
      any(!is.finite(grid[[1L]])) || !is.numeric(alpha) || length(alpha)!=1L ||
      !is.finite(alpha) || alpha<=0 || alpha>=1)
    stop("invalid M4 normal-band request",call.=FALSE)
  progress <- function(stage,done,total) .np_plot_progress_notify(
    stage=paste("Compact density M4:",stage),done=done,total=total,unit="work")
  inference <- .np_plot_corrected_activity_run("Constructing compact density normal intervals", {
    if (band=="all") {
      # Reuse the incumbent maximum and its exact stream/endpoints once.
      .np_with_seed(state$plot.random.seed,
        .np_udens_m4_inference(state,grid,alpha,B,progress))
    } else {
      targets <- sort(unique(grid[[1L]])); index <- match(grid[[1L]],targets)
      low <- .np_udens_m4_moments(state$x[[1L]],state$h,targets,24L,progress=progress)
      high <- .np_udens_m4_moments(state$x[[1L]],state$h,targets,48L,progress=progress)
      error <- max(abs(low$C-high$C)); limit <- 1e-12+1e-9*max(abs(high$C))
      center.error <- max(abs(low$center-high$center))
      if (error>limit || center.error>1e-9)
        stop("M4 quadrature refinement failed",call.=FALSE)
      eig <- eigen(high$C,symmetric=TRUE,only.values=TRUE)$values
      if (max(eig)<=0 || min(eig)< -1e-10*max(eig))
        stop("M4 covariance is not positive semidefinite",call.=FALSE)
      list(center=high$center[index],se=sqrt(diag(high$C))[index],bands=list(),
        covariance=NULL,checks=list(covariance.error=error,covariance.limit=limit,
          center.error=center.error,minimum.eigenvalue=min(eig)),
        settings=list(construction=state$construction,law="gaussian",h=state$h,n=state$n,q=1L,
          alpha=alpha,B=0L,grid.size=length(targets),family="unique supplied finite grid",
          family.rows=index,primary="absolute",
          covariance="complete corrected-curve influence under equal-weight Gaussian mixture pilot",
          correction="order-six compact operator with fitted normal guide"))
    }
  })
  .np_udens_whole_normal_assemble(inference,alpha,band)
}


.np_udens_uniform_normal_request <- function(bws,control,center,method,band,calibration=NULL) {
  is.null(control) && is.null(calibration) && identical(center,"bias-corrected") &&
    identical(method,"asymptotic") && isTRUE(band %in% c("pointwise","bonferroni","simultaneous","all")) &&
    inherits(bws,"bandwidth") && identical(bws[["ckertype",exact=TRUE]],"uniform") &&
    identical(bws[["type",exact=TRUE]],"fixed") &&
    identical(bws[["ckerbound",exact=TRUE]],"none") &&
    isTRUE(bws[["ckerorder",exact=TRUE]]==2L) &&
    isTRUE(bws[["ndim",exact=TRUE]]==1L) && isTRUE(bws[["ncon",exact=TRUE]]==1L)
}

.np_udens_uniform_normal_active <- function(state)
  identical(state[["construction",exact=TRUE]],"density-uniform-q6-normal-v1")

.np_udens_uniform_normal_prepare <- function(bws,x,seed) {
  if (!.np_udens_uniform_normal_request(bws,NULL,"bias-corrected","asymptotic","pointwise"))
    stop("invalid uniform normal specification",call.=FALSE)
  state <- .np_udens_q6_prepare(bws,x,seed)
  state$construction <- "density-uniform-q6-normal-v1"
  state$plot.law <- "gaussian"
  state$plot.calibration <- NULL
  state
}

# Observation influence of the signed compact curve, including the fitted normal
# guide. QR accumulates centered rows without retaining an n-by-grid matrix.
.np_udens_uniform_empirical_factor <- function(x, h, grid) {
  n <- length(x)
  location <- mean(x)
  variance <- mean((x - location)^2)
  sigma <- sqrt(variance)
  rule <- .np_udens_q6_rule(64L)
  weights <- rule$w * .np_udens_q6_kernel(rule$x, "uniform")
  nodes <- outer(grid, h * rule$x, "+")
  score <- function(t, scale) {
    z <- (t - location) / sigma
    dnorm(t, location, sigma) * if (scale) z^2 - 1 else z / sigma
  }
  guide <- cbind(score(grid, FALSE) - drop(score(nodes, FALSE) %*% weights),
                 score(grid, TRUE) - drop(score(nodes, TRUE) %*% weights))
  targets <- length(grid)
  count <- 0L
  average <- numeric(targets)
  triangular <- matrix(numeric(), 0L, targets)
  for (first in seq.int(1L, n, by = 64L)) {
    y <- x[first:min(n, first + 63L)]
    size <- length(y)
    u <- (y - location) / sigma
    rows <- .np_udens_q6_kernel(outer(y, grid, "-") / h, "uniform") / h +
      cbind(y - location, (u^2 - 1) / 2) %*% t(guide)
    block.mean <- colMeans(rows)
    delta <- block.mean - average
    joined <- rbind(triangular, sweep(rows, 2L, block.mean, "-"),
                    delta * sqrt(as.double(count) * size / (count + size)))
    if (any(!is.finite(joined)))
      stop("nonfinite empirical uniform density influence", call. = FALSE)
    decomposition <- qr(joined, LAPACK = TRUE)
    triangular <- qr.R(decomposition)[, order(decomposition$pivot), drop = FALSE]
    average <- average + delta * (size / (count + size))
    count <- count + size
  }
  # The symmetric square root fixes the draw orientation even at repeated or
  # zero singular values, without flooring variance or inverting the covariance.
  s <- svd(triangular / n, nu = 0L)
  sweep(s$v, 2L, s$d, "*") %*% t(s$v)
}

.np_udens_gaussian_absolute_maximum <- function(factor, B, progress, empirical = NULL) {
  m <- nrow(factor)
  maxima <- numeric(B)
  if (!is.null(empirical)) empirical.maxima <- numeric(B)
  for (first in seq.int(1L, B, by = 128L)) {
    ids <- first:min(B, first + 127L)
    z <- matrix(rnorm(m * length(ids)), m)
    values <- factor %*% z
    if (any(!is.finite(values))) stop("nonfinite compact density Gaussian process", call. = FALSE)
    maxima[ids] <- apply(abs(values), 2L, max)
    if (!is.null(empirical)) {
      values <- empirical %*% z
      if (any(!is.finite(values)))
        stop("nonfinite empirical uniform density Gaussian process", call. = FALSE)
      empirical.maxima[ids] <- apply(abs(values), 2L, max)
    }
    progress("Gaussian maximum", max(ids), B)
  }
  if (is.null(empirical)) maxima else cbind(Q6 = maxima, empirical = empirical.maxima)
}

.np_udens_uniform_normal_public <- function(state,grid,alpha,band,B = NULL) {
  if (!.np_udens_uniform_normal_active(state) || !isTRUE(band %in% c("pointwise","bonferroni","simultaneous","all")) ||
      !is.data.frame(grid) || ncol(grid)!=1L || !nrow(grid) ||
      !identical(names(grid),names(state$x)) || !is.numeric(grid[[1L]]) ||
      any(!is.finite(grid[[1L]])) || !is.numeric(alpha) || length(alpha)!=1L ||
      !is.finite(alpha) || alpha<=0 || alpha>=1)
    stop("invalid uniform normal-band request",call.=FALSE)
  progress <- function(stage,done,total) .np_plot_progress_notify(
    stage=paste("Compact density Q6 normal:",stage),done=done,total=total,unit="work")
  joint <- band %in% c("simultaneous", "all")
  if (joint) {
    B <- .np_plot_gaussian_integer(B, "density draws", .Machine$integer.max)
    m <- length(unique(grid[[1L]]))
    .np_udens_q6_budget(state$n, m, B = B, block = 128L)
    .np_plot_gaussian_draws(B, alpha, m)
  }
  ref <- .np_plot_corrected_activity_run("Constructing uniform density normal intervals",
    .np_udens_q6_reference(state$x[[1L]],state$h,grid[[1L]],"uniform",progress,keep.factor=joint))
  inference <- list(center=ref$center[ref$index],se=ref$se[ref$index],bands=list(),
    covariance=NULL,checks=ref$diagnostics,
    settings=list(construction=state$construction,law="gaussian",h=state$h,n=state$n,q=1L,
      alpha=alpha,B=0L,grid.size=length(ref$grid),family="unique supplied finite grid",
      family.rows=ref$index,primary="absolute",
      reference="empirical mixture plus h*sqrt(2/3)*Student6 noise",
      covariance="complete corrected-curve influence under Q6",
      correction="order-six compact operator with fitted normal guide"))
  if (joint) {
    empirical <- .np_plot_corrected_activity_run("Constructing empirical uniform density covariance",
      .np_udens_uniform_empirical_factor(state$x[[1L]], state$h, ref$grid))
    maximum <- .np_plot_corrected_activity_run("Calibrating the uniform density band",
      .np_with_seed(state$plot.random.seed, .np_udens_gaussian_absolute_maximum(ref$factor, B, progress, empirical)))
    component <- apply(maximum, 2L, quantile, probs = 1-alpha, type = 7L, names = FALSE)
    critical <- max(component)
    inference$settings$covariance <- "Q6 normal marginals; joint Q6 and complete empirical corrected-curve influences"
    inference$settings$critical.components <- component
    center <- inference$center
    inference$bands$simultaneous <- array(c(center-critical, center+critical),
      c(length(center), 2L, 1L), dimnames = list(NULL, c("lower", "upper"), as.character(alpha)))
    inference$critical <- critical
    inference$settings$B <- B
    inference$settings$quantile.type <- 7L
    inference$settings$critical.units <- "density"
    inference$settings$maximum.statistic <- "unstudentized absolute Gaussian maximum"
    inference$settings$calibration <- "Gaussian maximum and normal quantiles"
  }
  result <- .np_udens_whole_normal_assemble(inference,alpha,band)
  if (joint) result$inference$settings$calibration.status <- "unstudentized finite-grid Gaussian maximum in density units using the larger Q6 and empirical critical values"
  result
}
