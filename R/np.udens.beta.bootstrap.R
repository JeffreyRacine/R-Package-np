.np_udens_beta_q1_refits <- function(x,counts) {
  .Call(C_np_plot_udens_beta_refits,as.double(x),matrix(as.double(counts),nrow(counts)))
}

.np_udens_beta_q1_common <- function(x,h,grid,thetas,inner,outer,
                                     score=.np_udens_beta_q1_score,moment.plan=NULL) {
  base<-apply(thetas,2,min)
  if(is.null(moment.plan))moment.plan<-.np_udens_q1_operator_plan(base,h,grid)
  st<-.np_udens_q1_operator_rule(base,moment.plan,outer)
  columns<-.np_udens_q1_score_columns(grid,list(st,list(lv=log(x),l1v=log1p(-x))),h,length(x)^(-7/17),inner,score)
  P<-columns[[1L]];obs<-columns[[2L]];rm(columns)
  mu<-fitted<-matrix(0,nrow(thetas),length(grid));Tg<-cbind(log(grid),log1p(-grid))
  mass.error<-score.error<-0
  for(start in seq.int(1L,nrow(thetas),by=32L)) {
    ix<-seq.int(start,min(start+31L,nrow(thetas)));th<-thetas[ix,,drop=FALSE]
    a<-th[,1]+1;b<-th[,2]+1;A<-lbeta(a,b)
    lw<-log(st$prob)+st$T%*%t(sweep(th,2,base))
    lw<-sweep(lw,2,st$A-A,'+');zm<-apply(lw,2,.np_udens_beta_logsum)
    p<-exp(sweep(lw,2,zm));mass<-max(abs(expm1(zm)))
    exact<-cbind(digamma(a)-digamma(a+b),digamma(b)-digamma(a+b))
    err<-max(abs(crossprod(p,st$T)-exact))
    if(mass>1e-10 || err>1e-12+1e-10*max(abs(exact)))
      stop("beta plot common quadrature mass/score",call.=FALSE)
    mass.error<-max(mass.error,mass);score.error<-max(score.error,err)
    mu[ix,]<-t(P%*%p);fitted[ix,]<-exp(sweep(th%*%t(Tg),1,A))
  }
  list(mu=mu,fitted=fitted,observed=obs,diagnostics=c(mass=mass.error,score=score.error,nodes=st$nodes))
}


# P1 count-injection bridge. Public storage/streaming is a later gated provider.
