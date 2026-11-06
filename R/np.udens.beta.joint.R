# Joint guide and full residual covariance. Hidden until numerical/resource proof.
.np_udens_beta_basis <- function(x) {
  cbind(log(x[,1]),log1p(-x[,1]),log(x[,2]),log1p(-x[,2]),(2*x[,1]-1)*(2*x[,2]-1))
}

.np_udens_beta_joint_rule <- function(theta,order,moments=TRUE) {
  axes<-list(.np_udens_beta_marginal_rule(theta[1:2],order),
             .np_udens_beta_marginal_rule(theta[3:4],order))
  a<-axes[[1]];b<-axes[[2]];n1<-a$nodes;n2<-b$nodes
  if(n1*n2>600000L)stop("beta plot joint guide node resource gate",call.=FALSE)
  # The same normalization is shared by information and mass-only callers.
  # Region certificates do not need the target-by-five basis/information arrays.
  interaction<-as.vector(outer(2*exp(a$lv)-1,2*exp(b$lv)-1))
  lw<-rep(log(a$prob),n2)+rep(log(b$prob),each=n1)+theta[5]*interaction
  z<-.np_udens_beta_logsum(lw)
  if(!moments)return(list(axes=axes,A=a$A+b$A+z,nodes=n1*n2))
  T<-cbind(rep(a$lv,n2),rep(a$l1v,n2),rep(b$lv,each=n1),rep(b$l1v,each=n1),unname(interaction))
  p<-exp(lw-z);p<-p/sum(p)
  mean<-drop(crossprod(p,T));D<-sweep(T,2,mean);I<-crossprod(D,D*p)
  list(axes=axes,mean=mean,I=I,A=a$A+b$A+z,nodes=n1*n2)
}

# Stable exponential-family likelihood increment. Same objective, shape bounds,
# Armijo coefficient and stopping criteria; avoids subtracting two near-equal
# normalizers. This is the same identity used by the frozen SD guide refits.
.np_udens_beta_joint_objective_delta <- function(theta,delta,target,order=16L) {
  a<-.np_udens_beta_marginal_rule(theta[1:2],order)
  b<-.np_udens_beta_marginal_rule(theta[3:4],order)
  n1<-a$nodes;n2<-b$nodes
  if(n1*n2>600000L)stop('beta plot joint guide node resource gate',call.=FALSE)
  interaction<-as.vector(outer(2*exp(a$lv)-1,2*exp(b$lv)-1))
  lw<-rep(log(a$prob),n2)+rep(log(b$prob),each=n1)+theta[5]*interaction
  p<-exp(lw-.np_udens_beta_logsum(lw));p<-p/sum(p)
  z<-delta[1]*rep(a$lv,n2)+delta[2]*rep(a$l1v,n2)+
     delta[3]*rep(b$lv,each=n1)+delta[4]*rep(b$l1v,each=n1)+
     delta[5]*interaction-sum(target*delta)
  if(max(abs(z))<.1)log1p(sum(p*expm1(z))) else
    .np_udens_beta_logsum(log(p)+z)
}

.np_udens_beta_q2_fit <- function(x) {
  theta<-c(.np_udens_beta_q1_fit(x[,1])$theta,.np_udens_beta_q1_fit(x[,2])$theta,0)
  target<-colMeans(.np_udens_beta_basis(x));history<-list()
  for(it in 0:49) {
    st<-.np_udens_beta_joint_rule(theta,16L);grad<-st$mean-target
    direction<-solve(st$I,grad);decrement<-sum(grad*direction);obj<-st$A-sum(theta*target)
    history[[it+1L]]<-c(objective=obj,decrement=decrement)
    if(decrement<1e-14 && max(abs(grad))<2e-8)break
    step<-1
    for(j in 1:4)if(direction[j]>0)step<-min(step,.99*(theta[j]+1)/direction[j])
    accepted<-FALSE
    for(ls in 1:35) {
      trial<-theta-step*direction;next.st<-.np_udens_beta_joint_rule(trial,16L)
      if(.np_udens_beta_joint_objective_delta(theta,trial-theta,target)<= -1e-4*step*decrement+1e-14) {theta<-trial;accepted<-TRUE;break}
      step<-step*.5
    }
    if(!accepted)stop("beta plot joint guide likelihood line-search",call.=FALSE)
    if(it==49L)stop("beta plot joint guide likelihood stationarity",call.=FALSE)
  }
  refined<-.np_udens_beta_joint_rule(theta,32L)
  .np_udens_beta_check(refined$mean,st$mean,"joint guide mean",1e-10)
  .np_udens_beta_check(refined$I,st$I,"joint guide information")
  list(theta=theta,state=refined,history=history)
}

.np_udens_beta_axis_cache <- function(grid,axes,h,rb,inner) {
  parts<-index<-vector('list',2)
  for(j in 1:2) {
    values<-sort(unique(grid[,j]));index[[j]]<-match(grid[,j],values)
    parts[[j]]<-.np_udens_beta_axis(values,axes[[j]]$lv,axes[[j]]$l1v,h[j]^2,rb,inner)
  }
  list(parts=parts,index=index)
}



