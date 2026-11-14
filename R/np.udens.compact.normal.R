# Compact fixed-bound normal intervals: complete influence under the ratio pilot.
# The observation operator and the original corrected center are unchanged.
.np_udens_compact_normal_moments <- function(x,h,grid,kernel,
                                            progress=function(...) NULL,physical=NULL) {
  if (!is.numeric(x) || length(x)<10L || !is.numeric(grid) || !length(grid) ||
      any(!is.finite(c(x,grid,h))) || any(x<=0|x>=1) || any(grid<=0|grid>=1) ||
      length(h)!=1L || h<=0 || length(kernel)!=1L ||
      !kernel%in%c('uniform','epanechnikov'))
    stop('invalid compact normal q1 moment request',call.=FALSE)
  n<-length(x);m<-length(grid)
  if(8*(4*as.double(m)*n+60*m)>192*1024^2)
    stop('compact normal plot workspace exceeds 192 MiB',call.=FALSE)
  plan<-.np_udens_compact_operator_plan(grid,h,kernel)
  progress('original guide',0L,1L)
  theta<-.np_udens_beta_q1_fit(x)$theta
  progress('original guide',1L,1L)
  Tg<-cbind(log(grid),log1p(-grid));result<-vector('list',2L)
  for(level in 1:2) {
    order<-c(128L,256L)[level]
    progress(paste('reference pilot',level),0L,1L)
    p<-.np_udens_bounded_normal_pilot(x,theta,order,m)
    st<-.np_udens_compact_rule(theta,order,plan)
    sq<-.np_udens_compact_rule(p$fit$theta,order,plan,keep.log=TRUE)
    for(key in c('A','mean','I')) {
      .np_udens_beta_check(st[[key]],p$st[[key]],paste('compact normal original',key),
        if(key=='mean')1e-10 else 1e-9)
      .np_udens_beta_check(sq[[key]],p$sq[[key]],paste('compact normal reference',key),
        if(key=='mean')1e-10 else 1e-9)
    }
    # Both guide rules integrate the same normalized Q, not two pilot laws.
    w<-.np_udens_q1_ratio_weights(x,theta,st$A,st)
    cross.rule<-sq
    cross.rule$prob<-exp(sq$log.prob+drop(sq$T%*%(theta-p$fit$theta))+sq$A-st$A)
    .np_udens_beta_check(sum(cross.rule$prob),1,'compact normal cross-guide mass',1e-10)
    wq<-.np_udens_q1_ratio_weights(x,theta,st$A,cross.rule)
    .np_udens_beta_check(w$normal,p$normal,'compact normal normalization')
    .np_udens_beta_check(wq$normal,p$normal,'compact normal cross-guide normalization')
    bytes<-8*(12*32*max(n,st$nodes,sq$nodes)+48*(st$nodes+sq$nodes)+
      4*as.double(m)*n+60*m+sq$nodes)
    if(!is.finite(bytes)||bytes>192*1024^2)
      stop('compact normal plot workspace exceeds 192 MiB',call.=FALSE)
    # Retain the incumbent compact center's block/reduction order. Only its
    # center moments are needed; residual envelopes and random stores are not.
    tile<-max(1L,min(512L,floor(4*1024^2/(8*max(n,m)))))
    mu<-numeric(m)
    for(ids in .np_udens_beta_blocks(st$nodes,tile)) {
      A<-.np_udens_compact_operator(plan,exp(st$lv[ids]))-1
      mu<-mu+drop(A%*%st$prob[ids])
      progress(paste('original center',level),max(ids),st$nodes)
    }
    observed<-matrix(0,m,n)
    for(ids in .np_udens_beta_blocks(n,tile))
      observed[,ids]<-.np_udens_compact_operator(plan,x[ids],
        .np_udens_compact_geometry_slice(physical,observations=ids))-1
    fitted<-exp(drop(Tg%*%theta)-st$A)
    fittedQ<-exp(drop(Tg%*%p$fit$theta)-sq$A)
    G<-sweep(sweep(Tg,2,sq$mean),1,fittedQ,'*')
    z<-list(center=rowMeans(observed)-mu+fitted,se=numeric(m),mu=mu,
      fitted=fitted,meanQ=numeric(m),theta=theta,thetaQ=p$fit$theta,
      reference=list(normal=w$normal,V=p$V,I=sq$I,score_error=p$fit$score_error,
        iterations=p$fit$iterations),workspace=max(bytes,p$workspace),variance.checks=list())
    rm(observed,A)
    for(first in seq.int(1L,m,by=32L)) {
      ix<-seq.int(first,min(m,first+31L));pl<-plan
      for(key in c('target','mass','condition'))pl[[key]]<-plan[[key]][ix]
      pl$coefficients<-plan$coefficients[ix,,drop=FALSE]
      A<-.np_udens_compact_operator(pl,exp(st$lv))
      Aq<-.np_udens_compact_operator(pl,exp(sq$lv))
      B<-Aq%*%(sweep(sq$T,2,sq$mean)*sq$prob)
      D<-G[ix,,drop=FALSE]-B
      v<-.np_udens_bounded_normal_variance(A,st$T,w$prob,D,sq$I,n)
      meanQ<-drop(Aq%*%wq$prob)
      .np_udens_beta_check(v$meanA,meanQ,'compact normal cross-rule mean')
      z$se[ix]<-sqrt(v$value);z$meanQ[ix]<-meanQ
      z$variance.checks[[length(z$variance.checks)+1L]]<-v$check
      progress(paste('complete variance',level),max(ix),m)
    }
    result[[level]]<-z
  }
  fields<-c('center','se','mu','fitted','meanQ','thetaQ')
  result[[1L]]$checks<-setNames(lapply(fields,function(k).np_udens_beta_check(
    result[[1L]][[k]],result[[2L]][[k]],paste('compact normal refinement',k),
    if(k=='mu')1e-10 else 1e-9)),fields)
  result[[1L]]
}
