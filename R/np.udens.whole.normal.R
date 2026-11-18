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
    identical(method,"asymptotic") && isTRUE(band %in% c("pointwise","bonferroni")) &&
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

.np_udens_uniform_normal_public <- function(state,grid,alpha,band) {
  if (!.np_udens_uniform_normal_active(state) || !isTRUE(band %in% c("pointwise","bonferroni")) ||
      !is.data.frame(grid) || ncol(grid)!=1L || !nrow(grid) ||
      !identical(names(grid),names(state$x)) || !is.numeric(grid[[1L]]) ||
      any(!is.finite(grid[[1L]])) || !is.numeric(alpha) || length(alpha)!=1L ||
      !is.finite(alpha) || alpha<=0 || alpha>=1)
    stop("invalid uniform normal-band request",call.=FALSE)
  progress <- function(stage,done,total) .np_plot_progress_notify(
    stage=paste("Compact density Q6 normal:",stage),done=done,total=total,unit="work")
  ref <- .np_plot_corrected_activity_run("Constructing uniform density normal intervals",
    .np_udens_q6_reference(state$x[[1L]],state$h,grid[[1L]],"uniform",progress))
  inference <- list(center=ref$center[ref$index],se=ref$se[ref$index],bands=list(),
    covariance=NULL,checks=ref$diagnostics,
    settings=list(construction=state$construction,law="gaussian",h=state$h,n=state$n,q=1L,
      alpha=alpha,B=0L,grid.size=length(ref$grid),family="unique supplied finite grid",
      family.rows=ref$index,primary="absolute",
      reference="empirical mixture plus h*sqrt(2/3)*Student6 noise",
      covariance="complete corrected-curve influence under Q6",
      correction="order-six compact operator with fitted normal guide"))
  .np_udens_whole_normal_assemble(inference,alpha,band)
}
