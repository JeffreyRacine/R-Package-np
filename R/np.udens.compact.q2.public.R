# Exact new-owner admission; q1 and q2 normal retain their incumbent predicates.
.np_udens_compact_q2_request <- function(bws,control,center,method,band,calibration=NULL,sd.explicit=FALSE) {
  inherits(bws,'bandwidth') && is.null(control) &&
    isTRUE(bws[['ckertype',exact=TRUE]]%in%c('uniform','epanechnikov')) &&
    isTRUE(bws[['ndim',exact=TRUE]]==2L) && isTRUE(bws[['ncon',exact=TRUE]]==2L) &&
    identical(bws[['type',exact=TRUE]],'fixed') && isTRUE(bws[['ckerorder',exact=TRUE]]==2L) &&
    isTRUE(bws[['ckerbound',exact=TRUE]]%in%c('fixed','range')) &&
    identical(center,'bias-corrected') && isTRUE(band%in%c('pointwise','simultaneous','bonferroni','all')) &&
    ((identical(method,'asymptotic')&&is.null(calibration)&&isTRUE(band%in%c('pointwise','bonferroni'))) ||
     (identical(method,'bootstrap')&&isTRUE(sd.explicit)&&isTRUE(calibration%in%c('symmetric-sd','asymmetric-sd'))))
}
.np_udens_compact_q2_active <- function(state)
  identical(state[['construction',exact=TRUE]],'density-fixed-compact-q2-calibrated-v1')

.np_udens_compact_q2_public_prepare <- function(bws,x,method,bootstrap,band,calibration,alpha,seed) {
  if(!identical(bws[['ckerbound',exact=TRUE]],'fixed'))
    stop("Compact q2 corrected bands require known ckerbound='fixed' bounds; range bounds are not supported",call.=FALSE)
  if(!.np_udens_compact_q2_request(bws,NULL,'bias-corrected',method,band,calibration,identical(bootstrap,'inid')))
    stop('invalid compact q2 calibrated request',call.=FALSE)
  if(!is.data.frame(x)||nrow(x)<10L||ncol(x)!=2L||anyDuplicated(names(x))||any(!nzchar(names(x)))||
     !all(vapply(x,is.numeric,logical(1L)))||any(!is.finite(as.matrix(x))))
    stop('Compact q2 corrected bands require at least ten finite numeric observations',call.=FALSE)
  lower<-bws[['ckerlb',exact=TRUE]];upper<-bws[['ckerub',exact=TRUE]]
  if(!is.numeric(lower)||!is.numeric(upper)||length(lower)!=2L||length(upper)!=2L||
     any(!is.finite(c(lower,upper)))||any(upper<=lower))
    stop('Compact q2 corrected bands require finite ordered coordinatewise bounds',call.=FALSE)
  width<-upper-lower;jacobian<-prod(width)
  if(any(!is.finite(width))||!is.finite(jacobian)||jacobian<=0||!is.finite(1/jacobian))
    stop('Compact q2 corrected-band support units are not representable',call.=FALSE)
  normalized<-sweep(sweep(as.matrix(x),2,lower),2,width,'/')
  physical.h<-as.double(.np_physical_bandwidth(bws));h<-physical.h/width
  .np_udens_compact_q2_validate(normalized,h,normalized[1L,,drop=FALSE],bws[['ckertype',exact=TRUE]])
  if(!is.numeric(alpha)||length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1)
    stop('Compact q2 corrected bands require alpha in (0,1)',call.=FALSE)
  npValidateNonNegativeInteger(seed,'random.seed')
  list(x=x,normalized=normalized,h=h,physical.h=physical.h,lower=lower,width=width,jacobian=jacobian,kernel=bws[['ckertype',exact=TRUE]],
    q=2L,n=nrow(x),variant=if(method=='asymptotic')'asymptotic' else calibration,
    plot.law=if(method=='asymptotic')'gaussian' else 'inid',plot.random.seed=seed,
    construction='density-fixed-compact-q2-calibrated-v1')
}

.np_udens_compact_q2_hidden <- function(x,h,grid,kernel,alpha,band,B=9999L,seed=42L,
                                       progress=function(...)NULL,variant='asymptotic',physical=NULL) {
  .np_udens_compact_q2_validate(x,h,grid,kernel)
  if(!is.numeric(alpha)||length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1||
     length(variant)!=1L||!variant%in%c('asymptotic','symmetric-sd','asymmetric-sd')||
     length(band)!=1L||!band%in%c('pointwise','simultaneous','bonferroni','all'))
    stop('invalid compact q2 band request',call.=FALSE)
  if(is.null(physical)) {
    original<-grid;grid<-unique(grid);grid<-grid[do.call(order,as.data.frame(grid)),,drop=FALSE]
    key<-function(x)apply(x,1,function(z)paste(format(z,digits=17,scientific=TRUE),collapse=':'))
    restore<-match(key(original),key(grid))
  } else {
    canonical<-.np_udens_compact_physical_grid(grid,physical)
    grid<-canonical$grid;physical<-canonical$physical;restore<-canonical$restore
  }
  n<-nrow(x);m<-nrow(grid)
  sd<-variant!='asymptotic';draws<-sd||band%in%c('simultaneous','all')
  if(draws).np_udens_beta_random_request(if(sd)'multinomial' else 'normal',if(sd)n else m,B,seed,RNGkind())
  disk<-if(sd)4*as.double(n)*B+8*as.double(m)*B+262144 else if(draws)8*as.double(m)*B+8*as.double(m)*m+262144 else 0
  if(!is.finite(disk)||disk>512*1024^2||(sd&&8*as.double(min(64L,m))*B>32*1024^2))
    stop('compact q2 request exceeds owned-storage budget',call.=FALSE)
  store<-if(draws).np_udens_beta_store() else NULL
  if(draws)on.exit(.np_udens_beta_store_close(store),add=TRUE)
  prepared<-.np_udens_compact_q2_prepare(x,h,grid,kernel,!sd&&draws,store,progress,physical)
  region<-.np_udens_beta_regions(prepared$theta,prepared$Veta,grid,alpha,progress)
  rhat<-prepared$center-prepared$fitted;scales<-list(sqrt(prepared$pilot_variance),sqrt(prepared$empirical_variance))
  if(draws)stream<-.np_udens_beta_random_stream(store,if(sd)'multinomial' else 'normal',if(sd)n else m,B,seed)
  if(sd) {
    bootstrap<-.np_udens_compact_q2_sd_roots(prepared,stream,store,progress)
    if(band%in%c('bonferroni','all'))result<-.np_udens_bounded_sd(bootstrap$roots,rhat,scales,region,prepared$center,alpha,band,variant,progress)
    else {
      parts<-lapply(scales,function(s).np_udens_beta_stored_reduce(bootstrap$roots,rhat,s,region,prepared$center,alpha,band,variant))
      result<-parts[[1L]];result$lower<-pmin(parts[[1L]]$lower,parts[[2L]]$lower);result$upper<-pmax(parts[[1L]]$upper,parts[[2L]]$upper)
    }
    prepared$checks$bootstrap<-bootstrap$checks
  } else {
    companions<-.np_udens_bounded_analytic(rhat,region,scales,alpha,m)
    if(draws) {
      covariance<-.np_udens_compact_q2_covariance_read(prepared)
      critical.pilot<-.np_udens_beta_gaussian_critical(covariance,B,alpha,stream,progress);rm(covariance);gc(FALSE)
      T<-.np_udens_beta_basis(x);T<-sweep(T,2,colMeans(T));K<-prepared$B%*%solve(prepared$fit$state$I)
      residual<-prepared$observed-rowMeans(prepared$observed)-K%*%t(T)
      covariance<-tcrossprod(residual)/(n*(n-1));rm(residual,K,T);gc(FALSE)
      critical.empirical<-.np_udens_beta_gaussian_critical(covariance,B,alpha,stream,progress);rm(covariance);gc(FALSE)
      critical<-c(critical.pilot,critical.empirical)
      lo<-lapply(1:2,function(j)rhat+region$lower-critical[j]*scales[[j]])
      hi<-lapply(1:2,function(j)rhat+region$upper+critical[j]*scales[[j]])
      simultaneous<-list(lower=pmin(lo[[1L]],lo[[2L]]),upper=pmax(hi[[1L]],hi[[2L]]),critical=critical)
    }
    result<-switch(band,pointwise=companions$pointwise,bonferroni=companions$bonferroni,simultaneous=simultaneous,all=simultaneous)
    if(band=='all')result$bands<-list(pointwise=companions$pointwise,simultaneous=simultaneous,bonferroni=companions$bonferroni)
  }
  result$center<-prepared$center;result$pilot_se<-scales[[1L]];result$empirical_se<-scales[[2L]]
  for(k in c('lower','upper','center','pilot_se','empirical_se'))result[[k]]<-result[[k]][restore]
  if(band=='all')result$bands<-lapply(result$bands,function(z){z$lower<-z$lower[restore];z$upper<-z$upper[restore];z})
  result$B<-if(draws)B else 0L;result$family_size<-m;result$theta<-prepared$theta;result$checks<-prepared$checks
  result$construction<-switch(variant,asymptotic='compact-q2-pilot-empirical-envelope-plus-guide-region',
    'symmetric-sd'='compact-q2-residual-absolute-containing-hull-envelope','asymmetric-sd'='compact-q2-signed-residual-inversion-envelope')
  low<-if(sd&&band%in%c('bonferroni','all'))result$sd_companion$resolution.limited else draws&&(if(sd&&variant=='asymmetric-sd')floor(B*alpha/4)else floor(B*alpha/2))<20L
  if(low)warning(if(sd&&band%in%c('bonferroni','all'))
    'compact q2 Bonferroni calibration has low finite-draw tail resolution; returning finite empirical limits with the simultaneous hull' else
    'compact q2 calibration has low finite-draw tail resolution; returning the finite band',call.=FALSE)
  result
}
