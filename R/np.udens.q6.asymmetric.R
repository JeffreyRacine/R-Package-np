# Selected-coordinate Q6 is a separate, explicitly requested signed process.
.np_udens_q6_asymmetric_active <- function(state)
  identical(state[["construction", exact=TRUE]], "density-compact-q6-asymmetric-selected-v1")

.np_udens_q6_asymmetric_design <- function(design) {
  vary <- if (is.list(design)) design[["vary", exact=TRUE]] else NULL
  if (!is.numeric(vary) || length(vary)!=1L || !is.finite(vary) || vary!=1)
    stop("asymmetric Q6 requires an explicit q=1 empirical-range plot family", call.=FALSE)
  trim <- design[["xtrim", exact=TRUE]]
  full <- if (is.list(trim)) length(trim)==1L && is.numeric(trim[[1L]]) &&
    identical(as.double(trim[[1L]]), c(0,1)) else
    is.numeric(trim) && length(trim)==1L && is.finite(trim) && trim==0
  if (!isTRUE(full))
    stop("asymmetric Q6 currently requires the full empirical-range plot family (xtrim=0)", call.=FALSE)
  invisible(NULL)
}

.np_udens_q6_asymmetric_prepare <- function(bws, x, seed, design) {
  .np_udens_q6_asymmetric_design(design)
  state <- .np_udens_q6_prepare(bws, x, seed)
  state$plot.calibration <- "asymmetric-sd"
  state$construction <- "density-compact-q6-asymmetric-selected-v1"
  state
}

.np_udens_q6_asymmetric_family <- function(state, grid, design) {
  .np_udens_q6_asymmetric_design(design)
  if (!.np_udens_q6_asymmetric_active(state) || !is.data.frame(grid) || ncol(grid)!=1L ||
      !nrow(grid) || !identical(names(grid), names(state$x)) ||
      !is.numeric(grid[[1L]]) || any(!is.finite(grid[[1L]])))
    stop("invalid asymmetric Q6 empirical-range grid", call.=FALSE)
  targets <- sort(unique(as.double(grid[[1L]])))
  expected <- seq(min(state$x[[1L]]),max(state$x[[1L]]),length.out=length(targets))
  if (!identical(targets, expected))
    stop("asymmetric Q6 grid does not match its declared empirical-range family", call.=FALSE)
  list(targets=targets, fractions=seq(0,1,length.out=length(targets)))
}

.np_udens_q6_asymmetric_maximum <- function(Y, reference) {
  low.end <- apply(Y, 2L, min); high.end <- apply(Y, 2L, max)
  plus <- minus <- rep(-Inf, ncol(Y)); error <- 0
  mu <- colMeans(Y)
  guide <- list(mu=mu, sd=sqrt(colMeans(sweep(Y,2L,mu)^2)))
  low.rule <- .np_udens_q6_rule(48L); high.rule <- .np_udens_q6_rule(64L)
  for (j in seq_along(reference$grid)) {
    target <- low.end + reference$fractions[j]*(high.end-low.end)
    low <- .np_udens_q6_curve(Y,target,reference$h,reference$kernel,48L,guide,low.rule)
    high <- .np_udens_q6_curve(Y,target,reference$h,reference$kernel,64L,guide,high.rule)
    error <- max(error,.np_udens_q6_check(low,high,"selected conditional curve"))
    # Density at each selected coordinate, not the corrected fixed-target T_h(Q6).
    truth <- .np_udens_q6_density(reference$x,reference$h,target)
    root <- (high-truth)/reference$se[j]
    if (any(!is.finite(root))) stop("nonfinite selected Q6 root",call.=FALSE)
    plus <- pmax(plus,root); minus <- pmax(minus,-root)
  }
  list(plus=plus,minus=minus,refinement=error)
}

.np_udens_q6_asymmetric_intervals <- function(reference, plus, minus, alpha) {
  if (!is.numeric(plus) || !length(plus) || !is.numeric(minus) ||
      length(plus)!=length(minus) || any(!is.finite(c(plus,minus))) ||
      !is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=.5)
    stop("invalid asymmetric Q6 maximum calibration",call.=FALSE)
  critical <- c(lower=unname(quantile(plus,1-alpha/2,type=7L)),
                upper=unname(quantile(minus,1-alpha/2,type=7L)))
  index <- reference$index;center <- reference$center[index];scale <- reference$se[index]
  lower <- center-critical[["lower"]]*scale;upper <- center+critical[["upper"]]*scale
  if (any(!is.finite(c(lower,upper))) || any(lower>upper))
    stop("invalid asymmetric Q6 interval ordering",call.=FALSE)
  list(center=center,se=scale,critical=critical,lower=lower,upper=upper)
}

.np_udens_q6_asymmetric_budget <- function(n, m, B) {
  # Retain the incumbent bound and reserve signed-output/reduction copies.
  bytes <- .np_udens_q6_budget(n,m,B) + 8*(16*as.double(B)+512+4*as.double(m))
  if (!is.finite(bytes) || bytes > 192*1024^2)
    stop("asymmetric Q6 inference exceeds its 192 MiB workspace bound",call.=FALSE)
  bytes
}

.np_udens_q6_asymmetric_public <- function(state, grid, alpha, B, band, design) {
  if (!identical(band,"simultaneous"))
    stop("asymmetric Q6 currently supports simultaneous bands only",call.=FALSE)
  family <- .np_udens_q6_asymmetric_family(state,grid,design)
  B <- .np_plot_gaussian_integer(B,"Q6 B",.Machine$integer.max)
  .np_udens_q6_geometry(state$x[[1L]],state$h,grid[[1L]],state$kernel)
  .np_udens_q6_asymmetric_budget(state$n,length(family$targets),B)
  progress <- function(stage,done,total) .np_plot_progress_notify(
    stage=paste("Asymmetric Q6 density:",stage),done=done,total=total,unit="work")
  inference <- .np_plot_corrected_activity_run("Constructing asymmetric Q6 density intervals",
    .np_with_seed(state$plot.random.seed, {
      reference <- .np_udens_q6_reference(state$x[[1L]],state$h,grid[[1L]],state$kernel,progress)
      reference$fractions <- family$fractions
      tasks <- .np_udens_q6_plan(get(".Random.seed",.GlobalEnv,inherits=FALSE),state$n,B)
      draws <- .np_udens_q6_asymmetric_dispatch(tasks,reference)
      if (!is.matrix(draws) || !identical(dim(draws),c(B,3L)) || any(!is.finite(draws)))
        stop("incomplete asymmetric Q6 replicate collection",call.=FALSE)
      z <- .np_udens_q6_asymmetric_intervals(reference,draws[,1L],draws[,2L],alpha)
      limits <- array(c(z$lower,z$upper),c(nrow(grid),2L,1L),
        dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
      bands <- list(simultaneous=limits)
      resolution <- .np_udens_resolution(B,length(reference$grid),alpha,
        "simultaneous","smoothed-q6 selected-coordinate asymmetric-sd")
      resolution$quantile.type <- 7L
      resolution$interpolated.index <- 1+(B-1)*(1-alpha/2)
      list(center=z$center,se=z$se,critical=z$critical,maximum=draws[,1:2,drop=FALSE],
        bands=bands,intervals=bands,absolute.intervals=bands,covariance=NULL,resolution=resolution,
        checks=c(reference$diagnostics,conditional=max(draws[,3L])),
        settings=list(construction=state$construction,law="smoothed-q6",h=state$h,n=state$n,q=1L,
          alpha=alpha,B=B,grid.size=length(reference$grid),family.rows=reference$index,
          family="selected full empirical-range coordinates in each resample",family.design=design,
          available.bands="simultaneous",primary="absolute",quantile.type=7L,
          calibration="asymmetric-sd",root.center="Q6 density at the resample-selected coordinate",
          scale="original complete fixed-coordinate influence scale used as a fixed weighting function",
          maximum.statistic="signed complete corrected selected-coordinate root divided by original weighting scale",
          reference="empirical mixture plus h*sqrt(2/3)*Student6 noise",
          correction="order-six compact operator with refitted normal guide",
          calibration.status="explicit q1 selected-coordinate simultaneous Q6; pointwise/trimmed/q2 unqualified"))
    }))
  .np_band_resolution_warn(inference$resolution)
  list(boot.err=cbind(inference$critical[["lower"]]*inference$se,
    inference$critical[["upper"]]*inference$se,inference$center),
    boot.all.err=NULL,inference=inference,bxp=list())
}
