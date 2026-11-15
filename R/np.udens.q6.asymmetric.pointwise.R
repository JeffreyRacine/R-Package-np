# Signed marginal roots on the same selected-coordinate Q6 family.
.np_udens_q6_asymmetric_pointwise_block <- function(n, m, B) {
  .np_udens_q6_pointwise_budget(n,m,B)
  base <- .np_udens_q6_budget(n,m,B)
  targets <- min(m,floor((192*1024^2-base-64*m-16*B)/(48*B))-1)
  if (!is.finite(targets) || targets<1)
    stop("asymmetric Q6 inference exceeds its 192 MiB workspace bound",call.=FALSE)
  as.integer(targets)
}

.np_udens_q6_asymmetric_pointwise_chunk <- function(task, reference) {
  Y <- .np_udens_q6_samples(task,reference$x,reference$h)
  low.end <- apply(Y,2L,min);high.end <- apply(Y,2L,max)
  mu <- colMeans(Y);guide <- list(mu=mu,sd=sqrt(colMeans(sweep(Y,2L,mu)^2)))
  low.rule <- .np_udens_q6_rule(48L);high.rule <- .np_udens_q6_rule(64L)
  out <- matrix(0,task$bsz,length(reference$grid)+1L);error <- 0
  for(j in seq_along(reference$grid)) {
    target <- low.end+reference$fractions[j]*(high.end-low.end)
    low <- .np_udens_q6_curve(Y,target,reference$h,reference$kernel,48L,guide,low.rule)
    high <- .np_udens_q6_curve(Y,target,reference$h,reference$kernel,64L,guide,high.rule)
    error <- max(error,.np_udens_q6_check(low,high,"selected conditional curve"))
    out[,j] <- high-.np_udens_q6_density(reference$x,reference$h,target)
  }
  out[,ncol(out)] <- error
  out
}

.np_udens_q6_asymmetric_pointwise_reduce <- function(reference, tasks, alpha, companions=FALSE) {
  B <- sum(vapply(tasks,function(z)z$bsz,0L));m <- length(reference$grid)
  block <- .np_udens_q6_asymmetric_pointwise_block(length(reference$x),m,B)
  if (!is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=.5)
    stop("invalid asymmetric Q6 marginal alpha",call.=FALSE)
  cuts <- matrix(0,m,2L);error <- 0
  if(companions){raw <- cuts;plus <- minus <- rep(-Inf,B)}
  for(first in seq.int(1L,m,by=block)) {
    ix <- first:min(first+block-1L,m);part <- reference
    for(field in c("grid","fractions","center","reference.center","se"))part[[field]]<-reference[[field]][ix]
    part$covariance <- NULL
    draws <- .np_udens_q6_asymmetric_pointwise_dispatch(tasks,part)
    if(!is.matrix(draws) || !identical(dim(draws),c(B,length(ix)+1L)) || any(!is.finite(draws)))
      stop("incomplete asymmetric Q6 marginal replicate collection",call.=FALSE)
    for(j in seq_along(ix)) {
      cuts[ix[j],] <- quantile(draws[,j],c(alpha/2,1-alpha/2),type=7L,names=FALSE)
      if(companions) {
        raw[ix[j],] <- quantile(draws[,j],c(alpha/(2*m),1-alpha/(2*m)),type=7L,names=FALSE)
        standardized <- draws[,j]/reference$se[ix[j]]
        plus <- pmax(plus,standardized);minus <- pmax(minus,-standardized)
      }
    }
    error <- max(error,draws[,ncol(draws)])
    .np_plot_progress_notify(stage="Calibrating asymmetric Q6 marginal bands",done=max(ix),total=m,unit="targets")
  }
  make <- function(q)list(lower=reference$center-q[,2L],upper=reference$center-q[,1L])
  result <- list(pointwise=make(cuts),refinement=error)
  if(companions) {
    result$raw.bonferroni <- make(raw)
    critical <- c(lower=unname(quantile(plus,1-alpha/2,type=7L)),upper=unname(quantile(minus,1-alpha/2,type=7L)))
    result$simultaneous <- list(lower=reference$center-critical[["lower"]]*reference$se,
      upper=reference$center+critical[["upper"]]*reference$se)
    result$bonferroni <- list(lower=pmin(result$raw.bonferroni$lower,result$simultaneous$lower),
      upper=pmax(result$raw.bonferroni$upper,result$simultaneous$upper))
    result$maximum <- cbind(plus=plus,minus=minus);result$critical <- critical
  }
  result
}

.np_udens_q6_asymmetric_pointwise_public <- function(state,grid,alpha,B,band,design) {
  if(!isTRUE(band %in% c("pointwise","bonferroni","all")))
    stop("invalid asymmetric Q6 marginal request",call.=FALSE)
  family <- .np_udens_q6_asymmetric_family(state,grid,design)
  B <- .np_plot_gaussian_integer(B,"Q6 B",.Machine$integer.max)
  .np_udens_q6_geometry(state$x[[1L]],state$h,grid[[1L]],state$kernel)
  .np_udens_q6_asymmetric_pointwise_block(state$n,length(family$targets),B)
  progress <- function(stage,done,total).np_plot_progress_notify(stage=paste("Asymmetric Q6 density:",stage),done=done,total=total,unit="work")
  companions <- band %in% c("bonferroni","all")
  inference <- .np_plot_corrected_activity_run("Constructing asymmetric Q6 density marginal intervals",
    .np_with_seed(state$plot.random.seed,{
      reference <- .np_udens_q6_reference(state$x[[1L]],state$h,grid[[1L]],state$kernel,progress)
      reference$fractions <- family$fractions
      tasks <- .np_udens_q6_plan(get(".Random.seed",.GlobalEnv,inherits=FALSE),state$n,B)
      z <- .np_udens_q6_asymmetric_pointwise_reduce(reference,tasks,alpha,companions)
      make <- function(interval) {
        lo <- interval$lower[reference$index];hi <- interval$upper[reference$index]
        if(any(!is.finite(c(lo,hi))) || any(lo>hi))stop("invalid asymmetric Q6 interval ordering",call.=FALSE)
        array(c(lo,hi),c(nrow(grid),2L,1L),dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
      }
      requested <- if(band=="all")c("pointwise","simultaneous","bonferroni")else band
      bands <- lapply(z[requested],make)
      resolution <- .np_udens_resolution(B,length(reference$grid),alpha,requested,
        "smoothed-q6 selected-coordinate asymmetric-sd")
      resolution$quantile.type <- 7L
      resolution$interpolated.index <- 1+(B-1)*(1-alpha/ifelse(resolution$band=="bonferroni",2*length(reference$grid),2L))
      list(center=reference$center[reference$index],se=reference$se[reference$index],bands=bands,
        intervals=bands,absolute.intervals=bands,covariance=NULL,resolution=resolution,
        maximum=z$maximum,critical=z$critical,raw.bonferroni=if(companions)make(z$raw.bonferroni)else NULL,
        checks=c(reference$diagnostics,conditional=z$refinement),
        settings=list(construction=state$construction,law="smoothed-q6",h=state$h,n=state$n,q=1L,
          alpha=alpha,B=B,grid.size=length(reference$grid),family.rows=reference$index,
          family="selected full empirical-range coordinates in each resample",family.design=design,
          available.bands=requested,primary="absolute",quantile.type=7L,calibration="asymmetric-sd",
          root.center="Q6 density at the resample-selected coordinate",
          marginal.statistic="signed complete corrected selected-coordinate root",
          scale="original complete fixed-coordinate influence scale used as a fixed weighting function",
          bonferroni.construction=if(companions)"type-7 signed equal-tail marginal quantiles at alpha/(2m) with simultaneous hull; raw endpoints retained"else NULL,
          reference="empirical mixture plus h*sqrt(2/3)*Student6 noise",
          correction="order-six compact operator with refitted normal guide",
          calibration.status="explicit q1 selected-coordinate Q6; local coverage and finite-B limitations retained"))
    }))
  .np_band_resolution_warn(inference$resolution)
  widths <- lapply(inference$bands,function(v)cbind(inference$center-v[,1L,1L],v[,2L,1L]-inference$center))
  list(boot.err=cbind(widths[[if(band=="all")"pointwise"else band]],inference$center),
    boot.all.err=if(band=="all")widths else NULL,inference=inference,bxp=list())
}
