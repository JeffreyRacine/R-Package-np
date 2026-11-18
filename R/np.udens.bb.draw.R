# Dual-log transport is restricted to the explicit beta-binomial plot law.
.np_udens_bb_draw <- function(a,b) {
  if(!is.numeric(a) || !is.numeric(b) || length(a)!=length(b))
    stop("invalid beta-binomial shape vectors",call.=FALSE)
  .Call(C_np_plot_udens_bb_draw,as.double(a),as.double(b))
}

.np_udens_bb_fit <- function(x, lx, lq, weights = rep(1,length(x))) {
  n<-length(x)
  if(n<10L || any(!is.finite(x)) || any(x<0 | x>1) ||
      length(weights)!=n || any(!is.finite(weights)) || any(weights<0) || sum(weights)<=0)
    stop("invalid beta plot guide sample/weights",call.=FALSE)
  if(!is.numeric(lx) || !is.numeric(lq) || length(lx)!=n || length(lq)!=n ||
      any(!is.finite(c(lx,lq))) || any(lx>=0|lq>=0) ||
      max(abs(exp(lx)+exp(lq)-1))>1e-12 || max(abs(exp(lx)-x))>1e-14)
    stop("invalid beta-binomial paired log coordinates",call.=FALSE)
  w<-weights/sum(weights);target<-c(sum(w*lx),sum(w*lq))
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
