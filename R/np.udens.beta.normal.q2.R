# Complete fixed-bound beta normal variance; the original guide owns the center.
.np_udens_beta_normal_q2_request <- function(bws) {
  identical(bws[['ckertype',exact=TRUE]],'beta') &&
    identical(bws[['type',exact=TRUE]],'fixed') &&
    identical(bws[['ckerbound',exact=TRUE]],'fixed') &&
    isTRUE(bws[['ckerorder',exact=TRUE]]==2L) &&
    isTRUE(bws[['ndim',exact=TRUE]]==2L) && isTRUE(bws[['ncon',exact=TRUE]]==2L)
}

.np_udens_beta_normal_q2_fit <- function(target,theta) {
  for(it in 0:49) {
    st<-.np_udens_beta_joint_rule(theta,16L);grad<-st$mean-target
    ev<-eigen(st$I,symmetric=TRUE,only.values=TRUE)$values
    if(any(!is.finite(ev))||min(ev)<=0||max(ev)/min(ev)>1e13)
      stop('beta normal q2 reference information conditioning',call.=FALSE)
    d<-solve(st$I,grad);dec<-sum(grad*d)
    if(dec<1e-14&&max(abs(grad))<2e-8&&
       max(abs(grad))<=1e-12+2e-8*max(abs(target)))break
    step<-1
    for(j in 1:4)if(d[j]>0)step<-min(step,.99*(theta[j]+1)/d[j])
    accepted<-FALSE
    for(ls in 1:35) {
      trial<-theta-step*d
      delta<-.np_udens_beta_joint_objective_delta(theta,trial-theta,target)
      if(delta<=-1e-4*step*dec+1e-14){theta<-trial;accepted<-TRUE;break}
      step<-step*.5
    }
    if(!accepted||it==49L)
      stop('beta normal q2 reference score did not converge',call.=FALSE)
  }
  s<-.np_udens_beta_joint_rule(theta,32L)
  .np_udens_beta_check(s$mean,target,'q2 reference stationarity',2e-8)
  list(theta=theta,state=s,iterations=it,score_error=max(abs(s$mean-target)))
}

.np_udens_beta_normal_q2_budget <- function(n,m,theta,outer) {
  nodes<-vapply(list(theta[1:2],theta[3:4]),function(a)
    .np_udens_beta_marginal_rule(a,outer)$nodes,integer(1L))
  # Both guide preparations, row caches, ratio arrays and native work coexist.
  bytes<-8*(12*as.double(n)*max(nodes)+4*prod(as.double(nodes))+
    96*as.double(m)*sum(nodes)+40*as.double(m)^2+200*m+60*n)
  if(max(nodes)>1600L||prod(as.double(nodes))>2560000L||
     !is.finite(bytes)||bytes>192*1024^2)
    stop('beta normal q2 plot workspace exceeds 192 MiB',call.=FALSE)
  bytes
}

.np_udens_beta_normal_q2_moments <- function(x,h,grid,progress=function(...)NULL) {
  if(!is.matrix(x)||!is.numeric(x)||ncol(x)!=2L||nrow(x)<10L||
     !is.matrix(grid)||!is.numeric(grid)||ncol(grid)!=2L||!nrow(grid)||
     !is.numeric(h)||length(h)!=2L||any(!is.finite(c(x,h,grid)))||
     any(h<=0)||any(x<=0|x>=1)||any(grid<=0|grid>=1))
    stop('invalid beta normal q2 moment request',call.=FALSE)
  n<-nrow(x);m<-nrow(grid)
  # Cheap input-size refusal precedes guide or kernel preparation.
  if(8*(4*as.double(m)*n+40*as.double(m)^2)>192*1024^2)
    stop('beta normal q2 plot workspace exceeds 192 MiB',call.=FALSE)
  progress('original joint guide',0L,1L)
  fit<-.np_udens_beta_q2_fit(x);theta<-fit$theta
  progress('original joint guide',1L,1L);levels<-vector('list',2L);workspace<-0
  raw<-function(p).np_udens_beta_joint_raw_native(p$cache,p$axes,p$theta[5L],p$logz,p$ratio,FALSE)
  blocks<-.np_udens_beta_blocks(m,16L)
  for(level in 1:2) {
    inner<-c(64L,128L)[level];outer<-c(32L,64L)[level]
    width<-min(16L,m)
    workspace<-max(workspace,.np_udens_beta_normal_q2_budget(n,width,theta,outer))
    progress(paste('joint reference pilot',level),0L,1L)
    first<-blocks[[1L]]
    p<-.np_udens_beta_q2_moment_prepare(x,h,grid[first,,drop=FALSE],theta,inner,outer,fit$state)
    initial<-raw(p)[[2L]];target<-initial$mean[length(first)+seq_len(5L)]
    fq<-.np_udens_beta_normal_q2_fit(target,theta)
    workspace<-max(workspace,.np_udens_beta_normal_q2_budget(n,width,fq$theta,outer))
    pq<-.np_udens_beta_q2_moment_prepare(x,h,grid[first,,drop=FALSE],fq$theta,inner,outer,fq$state)
    z<-list(center=numeric(m),se=numeric(m),theta=theta,thetaQ=fq$theta,
      meanQ=numeric(m),D=matrix(0,m,5L),variance.checks=list())
    for(rows in blocks) {
      if(!identical(rows,first)) {
        p$grid<-pq$grid<-grid[rows,,drop=FALSE]
        p$cache<-.np_udens_beta_axis_cache(p$grid,p$axes,h,n^(-13/54),inner)
        pq$cache<-.np_udens_beta_axis_cache(pq$grid,pq$axes,h,n^(-13/54),inner)
      }
      ii<-seq_along(rows);tt<-length(rows)+seq_len(5L)
      qr<-if(identical(rows,first))initial else raw(p)[[2L]]
      .np_udens_beta_check(qr$mean[tt],target,'q2 blocked reference score',1e-10)
      gq<-raw(pq)[[1L]]
      Tg<-.np_udens_beta_basis(grid[rows,,drop=FALSE]);gf<-exp(drop(Tg%*%fq$theta)-fq$state$A)
      G<-sweep(sweep(Tg,2L,fq$state$mean),1L,gf,'*')
      D<-G-gq$cross[ii,,drop=FALSE];J<-D%*%solve(fq$state$I)
      cross<-qr$cross[ii,,drop=FALSE];V<-qr$cross[tt,,drop=FALSE]
      components<-cbind(kernel=qr$variance[ii],cross=2*rowSums(J*cross),
        guide=rowSums((J%*%V)*J))
      variance<-(components[,1L]+components[,2L]+components[,3L])/n
      if(any(!is.finite(variance))||any(variance<0))
        stop('beta normal q2 complete variance is not nonnegative',call.=FALSE)
      base<-.np_udens_beta_q2_moment_block(p,as.integer(ii),FALSE)
      z$center[rows]<-base$center;z$se[rows]<-sqrt(variance)
      z$meanQ[rows]<-qr$mean[ii];z$D[rows,]<-D
      if(is.null(z$reference))z$reference<-list(normal=exp(qr$logmass),V=V,I=fq$state$I,
        score_error=fq$score_error,iterations=fq$iterations)
      else .np_udens_beta_check(V,z$reference$V,'q2 blocked reference covariance')
      z$variance.checks[[length(z$variance.checks)+1L]]<-list(components=components,identities=base$identities)
      progress(paste('complete joint variance',level),max(rows),m)
    }
    levels[[level]]<-z
    rm(p,pq,initial,qr,gq,base)
  }
  fields<-c('center','se','thetaQ','meanQ','D')
  levels[[1L]]$checks<-setNames(lapply(fields,function(k).np_udens_beta_check(
    levels[[1L]][[k]],levels[[2L]][[k]],paste('N3 refinement',k))),fields)
  levels[[1L]]$workspace<-workspace
  levels[[1L]]
}

.np_udens_beta_normal_q2_public <- function(state,grid,alpha,band) {
  if(!identical(state$q,2L)||!identical(state$variant,'asymptotic')||!identical(band,'pmzsd'))
    stop('invalid beta normal q2 plot state',call.=FALSE)
  if(!is.data.frame(grid)||!nrow(grid)||!identical(names(grid),names(state$x))||
     !all(vapply(grid,is.numeric,logical(1L)))||any(!is.finite(as.matrix(grid))))
    stop('invalid beta normal q2 evaluation frame',call.=FALSE)
  normalized<-sweep(sweep(as.matrix(grid),2L,state$lower),2L,state$width,'/')
  if(any(normalized<=0|normalized>=1))
    stop('beta normal q2 bands require strictly interior evaluation points',call.=FALSE)
  progress<-function(stage,done,total).np_plot_progress_notify(
    stage=paste('Beta density normal:',stage),done=done,total=total,unit='work')
  z<-.np_plot_corrected_activity_run('Constructing fixed-bound beta normal intervals',
    .np_udens_beta_normal_q2_moments(state$normalized,state$h,normalized,progress))
  inference<-list(center=z$center/state$jacobian,se=z$se/state$jacobian,bands=list(),
    settings=list(construction='density-fixed-q2-beta-ratio-pilot-normal-v1',
      calibration='normal-quantile',law='asymptotic',kernel='beta',q=2L,n=state$n,
      alpha=alpha,B=0L,available.bands='pmzsd',grid.size=nrow(grid),
      family='marginal intervals at supplied interior targets',
      bounds=cbind(state$lower,state$lower+state$width),bandwidth=state$h*state$width,
      calibration.status='complete joint ratio-pilot variance; unchanged corrected center'),
    reference=list(original.guide=z$theta,variance.guide=z$thetaQ,
      normalization=z$reference$normal,score.error=z$reference$score_error,
      score.iterations=z$reference$iterations,score.covariance=z$reference$V,
      guide.information=z$reference$I),
    numerical.checks=list(refinement=z$checks,variance=z$variance.checks,workspace=z$workspace))
  inference<-.np_plot_bias_normal(inference,alpha)
  inference$intervals<-inference$absolute.intervals<-inference$bands
  endpoints<-inference$bands$pmzsd
  list(boot.err=cbind(inference$center-endpoints[,1L,1L],endpoints[,2L,1L]-inference$center,inference$center),
    boot.all.err=NULL,inference=inference,bxp=list())
}
