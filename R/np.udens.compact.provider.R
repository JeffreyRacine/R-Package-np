# A single finite target family, with common innovations for both covariance laws.
.np_udens_compact_hidden <- function(x,h,grid,kernel,alpha,band,B=9999L,seed=42L,
                                    progress=function(...)NULL,variant='asymptotic') {
  if(!identical(variant,'asymptotic'))
    return(.np_udens_compact_sd_hidden(x,h,grid,kernel,alpha,variant,band,B,seed,progress))
  if(!is.numeric(x)||length(x)<10L||!is.numeric(grid)||!length(grid)||
     any(!is.finite(c(x,grid,h)))||any(x<=0|x>=1)||any(grid<=0|grid>=1)||
     length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1||
     length(band)!=1L||!band%in%c('pointwise','bonferroni','simultaneous','all'))
    stop('invalid compact corrected-band request',call.=FALSE)
  original<-grid;grid<-sort(unique(grid));restore<-match(original,grid)
  draws<-band%in%c('simultaneous','all');m<-length(grid)
  if(draws).np_udens_beta_random_request('normal',m,B,seed,RNGkind())
  prepared<-.np_udens_compact_prepare(x,h,grid,kernel,draws,progress)
  region<-.np_udens_beta_regions(prepared$theta,prepared$Veta,matrix(grid,ncol=1L),alpha,progress)
  rhat<-prepared$center-prepared$fitted
  scales<-list(sqrt(prepared$pilot_variance),sqrt(prepared$empirical_variance))
  companions<-.np_udens_bounded_analytic(rhat,region,scales,alpha,m)
  if(draws) {
    store<-.np_udens_beta_store();on.exit(.np_udens_beta_store_close(store),add=TRUE)
    stream<-.np_udens_beta_random_stream(store,'normal',m,B,seed)
    critical<-vapply(list(prepared$covariance,prepared$empirical),
      function(v).np_udens_beta_gaussian_critical(v,B,alpha,stream,progress),numeric(1L))
    lo<-lapply(1:2,function(j)rhat+region$lower-critical[j]*scales[[j]])
    hi<-lapply(1:2,function(j)rhat+region$upper+critical[j]*scales[[j]])
    simultaneous<-list(lower=pmin(lo[[1]],lo[[2]]),upper=pmax(hi[[1]],hi[[2]]),critical=critical)
  }
  result<-switch(band,pointwise=companions$pointwise,bonferroni=companions$bonferroni,
    simultaneous=simultaneous,all=simultaneous)
  if(band=='all')result$bands<-list(pointwise=companions$pointwise,
    simultaneous=simultaneous,bonferroni=companions$bonferroni)
  result$center<-prepared$center;result$pilot_se<-scales[[1]];result$empirical_se<-scales[[2]]
  for(k in c('lower','upper','center','pilot_se','empirical_se'))result[[k]]<-result[[k]][restore]
  if(band=='all')result$bands<-lapply(result$bands,function(v) {
    v$lower<-v$lower[restore];v$upper<-v$upper[restore];v
  })
  result$B<-if(draws)B else 0L;result$family_size<-m
  result$theta<-prepared$theta;result$checks<-prepared$checks
  result$construction<-'compact-cubic-pilot-empirical-envelope-plus-guide-region'
  if(draws&&floor(B*alpha/2)<20L)
    warning('compact plot calibration has low finite-draw tail resolution; returning the finite band',call.=FALSE)
  result
}
