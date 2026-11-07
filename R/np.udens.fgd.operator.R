# Explicit Gaussian operator supplied to kernel-independent guide/moment leaves.
# The value is A-1: constants cancel exactly from the guide residual.
.np_udens_fgd_logmass <- function(x,h) {
  if(h<1)return(log(sqrt(2*pi)*h*(pnorm((1-x)/h)-pnorm(-x/h))))
  gl<-.np_udens_beta_rule(c(1,1),16L)
  log1p(drop(expm1(-outer(x,gl[,1],'-')^2/(2*h*h))%*%gl[,2]))
}

.np_udens_fgd_operator <- function(x,u,h) {
  if(!is.numeric(x)||!length(x)||!is.numeric(u)||!length(u)||
     any(!is.finite(c(x,u,h)))||length(h)!=1L||h<=0||
     any(x<0|x>1)||any(u<0|u>1))stop('invalid Gaussian plot operator geometry',call.=FALSE)
  coefficients<-matrix(0,length(x),4L)
  # The frozen local-polynomial operator evaluates log mass per target.
  # A matrix BLAS contraction can round differently near the uniform limit.
  mass<-if(h<1).np_udens_fgd_logmass(x,h) else
    vapply(x,function(xx).np_udens_fgd_logmass(xx,h),numeric(1))
  if(h>=1)gl<-.np_udens_beta_rule(c(1,1),32L)
  for(i in seq_along(x)) {
    if(h<1) {
      lo<--x[i]/h;hi<-(1-x[i])/h;den<-pnorm(hi)-pnorm(lo)
      mom<-numeric(7L);mom[1]<-1;mom[2]<-(dnorm(lo)-dnorm(hi))/den
      edge<-function(z,k){p<-dnorm(z);if(p==0)0 else z^k*p}
      for(k in 2:6)mom[k+1]<-(edge(lo,k-1)-edge(hi,k-1))/den+(k-1)*mom[k-1]
    } else {
      z<-gl[,1]-x[i];w<-gl[,2]*(1+expm1(-z^2/(2*h*h)-mass[i]))
      mom<-vapply(0:6,function(k)sum(w*z^k),numeric(1))
    }
    M<-outer(0:3,0:3,function(a,b)mom[a+b+1L])
    if(rcond(M)<1e-12)stop('Gaussian local-polynomial moment conditioning',call.=FALSE)
    coefficients[i,]<-solve(M,c(1,0,0,0))
  }
  .Call(C_np_plot_fgd_evaluate,as.double(x),as.double(u),as.double(h),coefficients,as.double(mass))
}

.np_udens_fgd_score <- function(grid,lv,l1v,rh,rb,order) {
  a<-.np_udens_fgd_operator(grid,exp(lv),sqrt(rh))
  # Retain the frozen reference's arithmetic ordering as well as its operator.
  # Although algebraically A-1, this combination affects tiny factor eigenmodes.
  ((4*a-6*a)+4*a)-a
}

.np_udens_fgd_axis_cache <- function(grid,axes,h,rb=NULL,inner=24L) {
  parts<-index<-vector('list',2L)
  for(j in 1:2) {
    target<-sort(unique(grid[,j]));index[[j]]<-as.integer(match(grid[,j],target))
    u<-exp(axes[[j]]$lv);a<-.np_udens_fgd_operator(target,u,h[j])
    z<-array(0,c(length(target),length(u),24L))
    for(s in 1:4) {
      z[,,s]<-a;z[,,s+4L]<-if(j==1L)a else 1
      z[,,9L+4L*(s-1L)]<-if(j==1L)1 else -.25*a
    }
    parts[[j]]<-z
  }
  list(parts=parts,index=index)
}
