# Dedicated plot-only beta numerics. Public admission remains unwired.
.np_udens_beta_rule <- function(shape, order) {
  .Call(C_np_plot_udens_beta_rule, as.double(shape), as.integer(order))
}

.np_udens_beta_axis <- function(targets, lv, l1v, rh, rb, order) {
  .Call(C_np_plot_udens_beta_axis, as.double(targets), as.double(lv),
        as.double(l1v), c(as.double(rh), as.double(rb)), as.integer(order))
}

.np_udens_beta_check <- function(a, b, name, relative = 1e-9) {
  error <- max(abs(a-b)); limit <- 1e-12 + relative * max(abs(b))
  if (!is.finite(error) || error > limit)
    stop(sprintf("beta plot %s refinement failed: %.17g > %.17g", name, error, limit), call.=FALSE)
  c(error=error,limit=limit)
}

.np_udens_beta_logsum <- function(x) {
  a<-max(x); if(!is.finite(a)) stop("invalid beta plot log mass",call.=FALSE)
  a+log(sum(exp(x-a)))
}

.np_udens_beta_marginal_rule <- function(theta, order) {
  a<-theta[1]+1;b<-theta[2]+1
  if(length(theta)!=2L || any(!is.finite(theta)) || min(a,b)<=0)
    stop("invalid beta plot guide shapes",call.=FALSE)
  budget<-70+abs(lbeta(a,b));TL<-max(40,budget/a);TR<-max(40,budget/b)
  if(max(TL,TR)>1e7)stop("unresolved beta plot guide boundary",call.=FALSE)
  fixed<-c(0,.5,1,2,3,5,10,20,40,80,160,320,640,1280,2560,5120,10240,
    20480,40960,81920,163840,327680,655360,1310720,2621440,5242880)
  edges<-sort(unique(c(-TL,-fixed[fixed<TL],fixed[fixed<TR],TR)))
  gl<-.np_udens_beta_rule(c(1,1),order);z<-2*gl[,1]-1;w<-2*gl[,2]
  t<-unlist(lapply(seq_len(length(edges)-1L),function(i)(edges[i]+edges[i+1])/2+diff(edges[i+0:1])*z/2))
  qw<-unlist(lapply(seq_len(length(edges)-1L),function(i)diff(edges[i+0:1])*w/2))
  if(length(t)>8192L)stop("beta plot marginal quadrature node budget",call.=FALSE)
  # Stable log-logistic coordinates, including underflowing ordinary values.
  softplus<-function(t)pmax(t,0)+log1p(exp(-abs(t)))
  lv<--softplus(-t);l1v<--softplus(t);T<-cbind(lv,l1v)
  lw<-drop(T%*%theta)+lv+l1v+log(qw);A<-.np_udens_beta_logsum(lw)
  p<-exp(lw-A);p<-p/sum(p);mean<-drop(crossprod(p,T));D<-sweep(T,2,mean)
  information<-crossprod(D,D*p)
  tail<-sum(vapply(list(c(a,TL),c(b,TR)),function(v)
    exp(-A-v[1]*v[2])*((v[2]+1)^2/v[1]+2*(v[2]+1)/v[1]^2+2/v[1]^3),numeric(1)))
  list(A=A,mean=mean,I=information,lv=lv,l1v=l1v,T=T,prob=p,tail=tail,nodes=length(t))
}

.np_udens_beta_q1_fit <- function(x, weights = rep(1,length(x))) {
  n<-length(x)
  if(n<10L || any(!is.finite(x)) || any(x<=0 | x>=1) ||
      length(weights)!=n || any(!is.finite(weights)) || any(weights<0) || sum(weights)<=0)
    stop("invalid beta plot guide sample/weights",call.=FALSE)
  w<-weights/sum(weights);target<-c(sum(w*log(x)),sum(w*log1p(-x)))
  mean<-sum(w*x);variance<-sum(w*(x-mean)^2);s<-mean*(1-mean)/variance-1
  if(!is.finite(s)||s<=0)stop("invalid beta plot guide moment start",call.=FALSE)
  theta<-c(mean*s,(1-mean)*s)-1;history<-list()
  for(it in 0:79) {
    st<-.np_udens_beta_marginal_rule(theta,64L);grad<-st$mean-target
    ev<-eigen(st$I,symmetric=TRUE,only.values=TRUE)$values
    if(min(ev)<=0 || max(ev)/min(ev)>1e13)stop("unqualified beta plot information conditioning",call.=FALSE)
    direction<-solve(st$I,grad);decrement<-sum(grad*direction);obj<-st$A-sum(theta*target)
    history[[it+1L]]<-c(objective=obj,decrement=decrement,gradient=max(abs(grad)))
    if(decrement<1e-14 && max(abs(grad))<2e-8)break
    step<-1
    for(j in 1:2)if(direction[j]>0)step<-min(step,.99*(theta[j]+1)/direction[j])
    accepted<-FALSE
    for(ls in 1:40) {
      trial<-theta-step*direction;next.st<-.np_udens_beta_marginal_rule(trial,64L)
      if(next.st$A-sum(trial*target)<=obj-1e-4*step*decrement+1e-14) {theta<-trial;accepted<-TRUE;break}
      step<-step*.5
    }
    if(!accepted)stop("beta plot guide likelihood line-search failed",call.=FALSE)
    if(it==79L)stop("beta plot guide likelihood stationarity not reached",call.=FALSE)
  }
  refined<-.np_udens_beta_marginal_rule(theta,128L)
  .np_udens_beta_check(refined$mean,st$mean,"guide mean",1e-10)
  .np_udens_beta_check(refined$I,st$I,"guide information")
  if(refined$tail>1e-12)stop("beta plot guide tail certificate failed",call.=FALSE)
  list(theta=theta,state=refined,history=history)
}

.np_udens_beta_q1_score <- function(grid, lv, l1v, rh, rb, order) {
  if(8*length(grid)*length(lv)>32*1024^2)stop("beta plot score tile budget",call.=FALSE)
  coef<-c(4,-6,4,-1);coefficient<-c(coef,coef,-as.vector(t(outer(coef,coef))))
  ans<-matrix(0,length(grid),length(lv))
  tile<-max(1L,min(8L,floor(32*1024^2/(24*8*length(lv)))))
  for(start in seq.int(1L,length(grid),by=tile)) {
    ix<-seq.int(start,min(length(grid),start+tile-1L))
    parts<-.np_udens_beta_axis(grid[ix],lv,l1v,rh,rb,order)
    for(j in seq_len(24L))ans[ix,]<-ans[ix,,drop=FALSE]+coefficient[j]*matrix(parts[,,j],length(ix))
  }
  ans
}

.np_udens_beta_q1_moments <- function(x,h,grid,theta,inner,outer,full=TRUE,
                                      score=.np_udens_beta_q1_score,variance.plan=NULL) {
  n<-length(x);rb<-n^(-7/17);st<-.np_udens_beta_marginal_rule(theta,outer)
  if(is.null(variance.plan))variance.plan<-.np_udens_q1_operator_plan(theta,h,grid)
  vr<-.np_udens_q1_operator_rule(theta,variance.plan,outer)
  for(key in c('A','mean','I')) .np_udens_beta_check(vr[[key]],st[[key]],
    paste('q1 operator guide',key),if(key=='mean')1e-10 else 1e-9)
  coordinates<-list(vr,list(lv=log(x),l1v=log1p(-x)))
  columns<-.np_udens_q1_score_columns(grid,coordinates,h,rb,inner,score)
  P<-Pvar<-columns[[1L]];obs<-columns[[2L]];rm(columns,coordinates)
  mu<-drop(P%*%vr$prob);S<-sweep(vr$T,2,st$mean);B<-P%*%(S*vr$prob)
  Td<-cbind(log(x),log1p(-x));U<-t(solve(st$I,t(sweep(Td,2,colMeans(Td)))))
  Rdata<-obs-rowMeans(obs)-B%*%t(U);Ueta<-sweep(U,2,theta+1,"/")
  Veta<-crossprod(Ueta)/(n*(n-1));Ugrid<-solve(st$I,t(sweep(vr$T,2,st$mean)))
  R<-Pvar-mu-B%*%Ugrid;lgx<-drop(Td%*%theta)-st$A
  rho<-n^(-2/5);lp<-numeric(vr$nodes)
  for(start in seq.int(1L,vr$nodes,by=32L)) {
    ix<-seq.int(start,min(vr$nodes,start+31L));v<-exp(vr$lv[ix]);a<-v/rho;b<-(1-v)/rho
    logs<-outer(a,log(x))+outer(b,log1p(-x))-lbeta(a+1,b+1)
    logs<-sweep(logs,2,lgx)
    lp[ix]<-apply(logs,1,.np_udens_beta_logsum)-log(n)
  }
  lw<-log(vr$prob)+lp;normal<-.np_udens_beta_logsum(lw);pp<-exp(lw-normal);pp<-pp/sum(pp)
  Rc<-R-drop(R%*%pp);weighted<-sweep(Rc,2,sqrt(pp),"*")
  covariance<-if(full)tcrossprod(weighted)/n else NULL
  empirical<-if(full)tcrossprod(Rdata)/(n*(n-1)) else NULL
  guide.mean<-drop(R%*%vr$prob)
  if(max(abs(guide.mean))>1e-9)stop('beta plot guide residual is not centered',call.=FALSE)
  guide.rows<-sweep(R-guide.mean,2,sqrt(vr$prob),'*')
  guide.covariance<-if(full)tcrossprod(guide.rows)/n else NULL
  guide.variance<-rowSums(guide.rows^2)/n
  pilot.variance<-rowSums(weighted^2)/n;empirical.variance<-rowSums(Rdata^2)/(n*(n-1))
  if(any(!is.finite(c(pilot.variance,empirical.variance)))||any(pilot.variance<=0)||any(empirical.variance<=0))
    stop("beta plot residual scale is not positive",call.=FALSE)
  fitted<-exp(drop(cbind(log(grid),log1p(-grid))%*%theta)-st$A)
  G<-sweep(sweep(cbind(log(grid),log1p(-grid)),2,st$mean),1,fitted,'*')
  identity.error<-0
  for(first in seq.int(1L,length(grid),by=32L)) {
    ix<-seq.int(first,min(length(grid),first+31L))
    corrected<-R[ix,,drop=FALSE]+G[ix,,drop=FALSE]%*%Ugrid
    direct<-Pvar[ix,,drop=FALSE]-mu[ix]+(G[ix,,drop=FALSE]-B[ix,,drop=FALSE])%*%Ugrid
    identity.error<-max(identity.error,.np_udens_beta_check(corrected,direct,'complete q1 influence')[1L])
  }
  if(full)for(cv in list(covariance,empirical,guide.covariance))
    if(min(eigen(cv,symmetric=TRUE,only.values=TRUE)$values)< -1e-10)
      stop('beta plot q1 residual covariance PSD',call.=FALSE)
  list(center=rowMeans(obs)-mu+fitted,fitted=fitted,mu=mu,B=B,observed=obs,
    covariance=covariance,empirical=empirical,Veta=Veta,ratio_normalization=exp(normal),
    pilot_variance=pilot.variance,empirical_variance=empirical.variance,
    guide_covariance=guide.covariance,guide_variance=guide.variance,identity_error=identity.error)
}
