.np_udens_fgd_request <- function(bws,control,center,method,calibration=NULL) {
  inherits(bws,'bandwidth') && is.null(control) &&
    !identical(calibration,'legacy') &&
    identical(bws[['ckertype',exact=TRUE]],'gaussian') &&
    isTRUE(bws[['ckerbound',exact=TRUE]]%in%c('fixed','range')) &&
    identical(center,'bias-corrected') && isTRUE(method%in%c('asymptotic','bootstrap'))
}

.np_udens_fgd_active <- function(state) {
  identical(state[['construction',exact=TRUE]],'density-fixed-gaussian-fgd-v1')
}

.np_udens_fgd_public_prepare <- function(bws,x,method,bootstrap,band,calibration,alpha,seed) {
  if(!identical(bws[['ckerbound',exact=TRUE]],'fixed'))
    stop("Gaussian corrected bands require known ckerbound='fixed' bounds; range bounds are not supported",call.=FALSE)
  if(!identical(bws[['type',exact=TRUE]],'fixed') || !isTRUE(bws[['ckerorder',exact=TRUE]]==2L) ||
     !isTRUE(bws[['ndim',exact=TRUE]]%in%1:2) || !isTRUE(bws[['ncon',exact=TRUE]]==bws[['ndim',exact=TRUE]]))
    stop('Gaussian corrected bands require fixed order-2 bandwidths and one or two continuous coordinates',call.=FALSE)
  normal <- .np_udens_bounded_normal_request(bws,method,band,calibration,alpha)
  if(!normal && !band%in%c('pointwise','simultaneous','bonferroni','all'))
    stop("Gaussian corrected bands require explicit band='pointwise', 'simultaneous', 'bonferroni' or 'all'",call.=FALSE)
  variant<-if(identical(method,'asymptotic'))'asymptotic' else calibration
  if(identical(method,'asymptotic')&&!is.null(calibration))
    stop('Gaussian asymptotic bands do not accept bootstrap calibration',call.=FALSE)
  if(identical(method,'bootstrap')&&(!identical(bootstrap,'inid')||is.null(variant)||!variant%in%c('symmetric-sd','asymmetric-sd')))
    stop("Gaussian bootstrap bands require bootstrap='inid' and calibration='symmetric-sd' or 'asymmetric-sd'",call.=FALSE)
  if(!is.data.frame(x)||nrow(x)<10L||ncol(x)!=bws[['ndim',exact=TRUE]]||
     anyDuplicated(names(x))||any(!nzchar(names(x)))||!all(vapply(x,is.numeric,logical(1L)))||any(!is.finite(as.matrix(x))))
    stop('Gaussian corrected bands require at least ten finite numeric observations',call.=FALSE)
  lower<-bws[['ckerlb',exact=TRUE]];upper<-bws[['ckerub',exact=TRUE]]
  if(!is.numeric(lower)||!is.numeric(upper)||length(lower)!=ncol(x)||length(upper)!=ncol(x)||
     any(!is.finite(c(lower,upper)))||any(upper<=lower))
    stop('Gaussian corrected bands require finite ordered coordinatewise bounds',call.=FALSE)
  width<-upper-lower;jacobian<-prod(width)
  if(any(!is.finite(width))||!is.finite(jacobian)||jacobian<=0||!is.finite(1/jacobian))
    stop('Gaussian corrected-band support units are not representable',call.=FALSE)
  normalized<-sweep(sweep(as.matrix(x),2,lower),2,width,'/')
  h<-as.double(.np_physical_bandwidth(bws))/width
  if(length(h)!=ncol(x)||any(!is.finite(h))||any(h<=0)||any(normalized<=0|normalized>=1))
    stop('Gaussian corrected bands require positive bandwidths and strictly interior observations',call.=FALSE)
  if(!is.numeric(alpha)||length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1)
    stop('Gaussian corrected bands require alpha in (0,1)',call.=FALSE)
  npValidateNonNegativeInteger(seed,'random.seed')
  state<-list(x=x,normalized=normalized,h=h,lower=lower,width=width,jacobian=jacobian,
    q=ncol(x),n=nrow(x),variant=variant,plot.law=if(method=='asymptotic')'gaussian' else 'inid',
    plot.random.seed=seed,construction='density-fixed-gaussian-fgd-v1')
  if(normal)state$normal.pmzsd<-TRUE
  state
}

.np_udens_fgd_public_errors <- function(state,grid,alpha,B,band) {
  if(identical(state[['normal.pmzsd',exact=TRUE]],TRUE))
    return(.np_udens_bounded_normal_public(state,grid,alpha,band,'gaussian'))
  if(!is.data.frame(grid)||nrow(grid)<1L||!identical(names(grid),names(state$x))||
     !all(vapply(grid,is.numeric,logical(1L)))||any(!is.finite(as.matrix(grid))))
    stop('invalid Gaussian corrected-band evaluation frame',call.=FALSE)
  normalized<-sweep(sweep(as.matrix(grid),2,state$lower),2,state$width,'/')
  if(any(normalized<=0|normalized>=1))
    stop('Gaussian corrected bands require strictly interior evaluation points; literal support endpoints are excluded',call.=FALSE)
  progress<-function(stage,done,total).np_plot_progress_notify(stage=paste('Gaussian density:',stage),done=done,total=total,unit='work')
  result<-.np_plot_corrected_activity_run('Constructing fixed-bound Gaussian density intervals',
    .np_udens_fgd_hidden(state$normalized,state$h,normalized,alpha,state$variant,band,B,state$plot.random.seed,progress))
  for(key in c('center','lower','upper','pilot_se','empirical_se'))result[[key]]<-result[[key]]/state$jacobian
  if(any(!is.finite(unlist(result[c('center','lower','upper','pilot_se','empirical_se')]))))
    stop('Gaussian corrected-band output units are not representable',call.=FALSE)
  interval<-array(c(result$lower,result$upper),c(nrow(grid),2L,1L),dimnames=list(NULL,c('lower','upper'),as.character(alpha)))
  intervals<-setNames(list(interval),band)
  inference<-list(center=result$center,bands=intervals,intervals=intervals,absolute.intervals=intervals,
    pilot.residual.scale=result$pilot_se,empirical.residual.scale=result$empirical_se,
    resolution=list(draws=result$B,signed.lower.index=result$k_signed,
      signed.upper.index=if(!is.null(result$k_signed))result$B+1L-result$k_signed else NULL,
      absolute.index=result$k_absolute,guide.risk=alpha/2,residual.risk=alpha/2,
      signed.tail.probability=alpha/(4*if(band=='bonferroni'&&state$variant!='asymptotic')result$family_size else 1),
      absolute.tail.probability=alpha/(2*if(band=='bonferroni'&&state$variant!='asymptotic')result$family_size else 1),
      gaussian.quantile.type=if(state$variant=='asymptotic'&&band=='simultaneous')7L else NULL),
    settings=list(construction=state$construction,calibration=state$variant,endpoint.construction=result$construction,law=state$plot.law,
      family='one supplied finite interior grid',grid.size=nrow(grid),q=state$q,n=state$n,
      alpha=alpha,B=result$B,resolution.limited=result$B>0L&&(if(state$variant=='asymmetric-sd')floor(result$B*alpha/4) else floor(result$B*alpha/2))<20L,
      critical=list(pilot=result$pilot_critical,empirical=result$empirical_critical),available.bands=band,
      bounds=cbind(state$lower,state$lower+state$width),bandwidth=state$h*state$width,
      guide.risk=alpha/2,residual.risk=alpha/2,
      calibration.status='fixed-bound Gaussian complete pilot/empirical residual envelope plus guide region'),
    numerical.checks=result$checks)
  out<-list(boot.err=cbind(result$center-result$lower,result$upper-result$center,result$center),
    boot.all.err=NULL,inference=inference,bxp=list())
  if(band%in%c('bonferroni','all'))out<-.np_udens_bounded_public_companions(out,result,band,state$jacobian,alpha)
  out
}
