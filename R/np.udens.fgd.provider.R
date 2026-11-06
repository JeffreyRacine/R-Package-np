# Combine complete constituent intervals. In particular signed roots cannot be
# combined by replacing the two scales with their pointwise maximum.
.np_udens_fgd_envelope <- function(pilot,empirical) {
  if(length(pilot$lower)!=length(empirical$lower)||
     length(pilot$upper)!=length(empirical$upper)||
     any(!is.finite(c(pilot$lower,pilot$upper,empirical$lower,empirical$upper)))||
     any(pilot$lower>pilot$upper)||any(empirical$lower>empirical$upper))
    stop('invalid Gaussian plot constituent intervals',call.=FALSE)
  list(lower=pmin(pilot$lower,empirical$lower),upper=pmax(pilot$upper,empirical$upper),
    pilot_critical=pilot$critical,empirical_critical=empirical$critical)
}

.np_udens_fgd_hidden <- function(x,h,grid,alpha,variant,band,B=9999L,seed=42L,
                                 progress=function(...)NULL) {
  if(!is.matrix(x)||!is.matrix(grid)||!ncol(x)%in%1:2||ncol(x)!=ncol(grid)||nrow(x)<10L||!nrow(grid)||
     any(!is.finite(c(x,grid,h)))||any(x<=0|x>=1)||any(grid<=0|grid>=1)||
     length(h)!=ncol(x)||any(h<=0)||length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1||
     length(variant)!=1L||!variant%in%c('asymptotic','symmetric-sd','asymmetric-sd')||
     length(band)!=1L||!band%in%c('pointwise','simultaneous'))stop('invalid Gaussian plot hidden provider request',call.=FALSE)
  original<-grid;unique.grid<-unique(grid);ordering<-do.call(order,as.data.frame(unique.grid))
  grid<-unique.grid[ordering,,drop=FALSE]
  keys<-function(x)apply(x,1,function(z)paste(format(z,digits=17,scientific=TRUE),collapse=':'))
  restore<-match(keys(original),keys(grid))
  sd<-variant!='asymptotic';draws<-sd||band=='simultaneous';m<-nrow(grid)
  if(draws).np_udens_beta_random_request(if(sd)'multinomial' else 'normal',if(sd)nrow(x) else m,B,seed,RNGkind())
  .np_udens_fgd_disk(nrow(x),m,ncol(x),B,sd,draws)
  store<-if(draws).np_udens_beta_store() else NULL
  if(draws)on.exit(.np_udens_beta_store_close(store),add=TRUE)
  prepared<-.np_udens_fgd_prepare(x,h,grid,!sd&&band=='simultaneous',store,progress,B,0L)
  region<-.np_udens_beta_regions(prepared$theta,prepared$Veta,grid,alpha,progress)
  scales<-list(sqrt(prepared$pilot_variance),sqrt(prepared$empirical_variance))
  rhat<-prepared$center-prepared$fitted
  if(any(!is.finite(unlist(scales)))||any(unlist(scales)<=0))stop('Gaussian plot residual scales must be positive',call.=FALSE)
  if(draws)stream<-.np_udens_beta_random_stream(store,if(sd)'multinomial' else 'normal',
    if(sd)nrow(x) else m,B,seed)
  constituents<-vector('list',2L)
  if(sd) {
    bootstrap<-.np_udens_beta_sd_roots(x,h,grid,prepared,stream,store,progress,
      .np_udens_fgd_score,.np_udens_fgd_axis_cache,c(24L,48L))
    for(j in 1:2)constituents[[j]]<-.np_udens_beta_stored_reduce(bootstrap$roots,rhat,scales[[j]],region,prepared$center,alpha,band,variant)
  } else for(j in 1:2) {
    if(band=='pointwise')critical<-qnorm(1-alpha/4) else {
      path<-prepared$covariance_files[[c('covariance','empirical')[j]]]
      covariance<-.np_udens_fgd_read(path,m,m)
      critical<-.np_udens_beta_gaussian_critical(covariance,B,alpha,stream,progress)
      rm(covariance);gc(FALSE)
    }
    constituents[[j]]<-list(lower=rhat+region$lower-critical*scales[[j]],
      upper=rhat+region$upper+critical*scales[[j]],critical=critical)
  }
  result<-.np_udens_fgd_envelope(constituents[[1]],constituents[[2]])
  if(sd) {
    result$k_signed<-constituents[[1]]$k_signed
    result$k_absolute<-constituents[[1]]$k_absolute
  }
  result$center<-prepared$center;result$pilot_se<-scales[[1]];result$empirical_se<-scales[[2]]
  for(name in c('lower','upper','center','pilot_se','empirical_se'))result[[name]]<-result[[name]][restore]
  result$theta<-prepared$theta;result$variant<-variant;result$band<-band;result$alpha<-alpha
  result$B<-if(draws)B else 0L
  result$construction<-'complete pilot/empirical residual interval envelope with guide region'
  result$checks<-prepared$checks
  if(draws&&(if(sd&&variant=='asymmetric-sd')floor(B*alpha/4) else floor(B*alpha/2))<20)
    warning('Gaussian plot calibration has low finite-draw tail resolution; returning the finite band',call.=FALSE)
  result
}
