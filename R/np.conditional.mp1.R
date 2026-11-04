# Asymmetric MP1 only. Frozen LCD12 arithmetic; no RNG or resampling in worker leaves.
# Private source-loaded exact conditional-mixture variance. No RNG or fitting.
.np_conditional_mp1_prepare <- function(z,tile=32L) {
  A<-2*z$w;D<- -z$jw
  Q1<-Q2<-Q3<-matrix(0,z$nx,z$n)
  for(first in seq.int(1L,z$n,by=tile)) {
    ii<-first:min(z$n,first+tile-1L);P<-.np_conditional_den_probabilities(z,ii)
    Q1<-Q1+crossprod(A[ii,,drop=FALSE]^2,P)
    Q2<-Q2+crossprod(D[ii,,drop=FALSE]^2,P)
    Q3<-Q3+crossprod(2*A[ii,,drop=FALSE]*D[ii,,drop=FALSE],P)
  }
  list(A=A,D=D,Q1=Q1,Q2=Q2,Q3=Q3,tile=tile)
}
.np_conditional_mp1_variance <- function(z,Y,prep=.np_conditional_mp1_prepare(z)) {
  stopifnot(is.matrix(Y),nrow(Y)==z$n,all(is.finite(Y)))
  nb<-ncol(Y);ny<-z$ny;h<-z$h[2];H<-z$H[2];s2<-2*h*h
  # Packed columns: all replicates for y1, then all replicates for y2, etc.
  F1<-F2<-G1<-G2<-G3<-matrix(0,z$n,ny*nb)
  for(j in seq_len(ny)) {
    cc<-(j-1L)*nb+seq_len(nb);v<-z$gy[j]-Y
    F1[,cc]<-dnorm(v,sd=sqrt(s2+h*h))
    F2[,cc]<-dnorm(v,sd=sqrt(s2+H*H))
    G1[,cc]<-dnorm(0,sd=sqrt(2)*h)*dnorm(v,sd=sqrt(s2+h*h/2))
    G2[,cc]<-dnorm(0,sd=sqrt(2)*H)*dnorm(v,sd=sqrt(s2+H*H/2))
    G3[,cc]<-dnorm(0,sd=sqrt(h*h+H*H))*dnorm(v,sd=sqrt(s2+h*h*H*H/(h*h+H*H)))
  }
  out<-prep$Q1%*%G1+prep$Q2%*%G2+prep$Q3%*%G3
  for(first in seq.int(1L,z$n,by=prep$tile)) {
    ii<-first:min(z$n,first+prep$tile-1L);P<-.np_conditional_den_probabilities(z,ii)
    m1<-P%*%F1;m2<-P%*%F2
    out<-out-crossprod(prep$A[ii,,drop=FALSE]^2,m1*m1)-
      crossprod(prep$D[ii,,drop=FALSE]^2,m2*m2)-
      crossprod(2*prep$A[ii,,drop=FALSE]*prep$D[ii,,drop=FALSE],m1*m2)
  }
  V<-matrix(0,z$m,nb)
  for(j in seq_len(ny)) {
    ii<-which(z$iy==j);cc<-(j-1L)*nb+seq_len(nb)
    V[ii,]<-out[z$ix[ii],cc,drop=FALSE]
  }
  if(any(!is.finite(V))||any(V<=0))stop('nonpositive/nonfinite replicate variance; no floor or fallback')
  V
}

# MP1: independent per-observation equal pilot mixtures; no fits or RNG here.
.np_conditional_mp1_between <- function(z,Y,p0,tile=32L) {
 nb<-ncol(Y);ny<-z$ny;h<-z$h[2];H<-z$H[2]
 F1<-F2<-matrix(0,z$n,ny*nb)
 for(j in seq_len(ny)){cc<-(j-1L)*nb+seq_len(nb);v<-z$gy[j]-Y
  F1[,cc]<-dnorm(v,sd=sqrt(3)*h);F2[,cc]<-dnorm(v,sd=2*h)}
 out<-matrix(0,z$nx,ny*nb);which_y<-rep(seq_len(ny),each=nb)
 for(first in seq.int(1L,z$n,by=tile)){
  ii<-first:min(z$n,first+tile-1L);P<-.np_conditional_den_probabilities(z,ii)
  d1<-P%*%F1-p0$mean[ii,which_y,drop=FALSE]
  d2<-P%*%F2-p0$mean[ii,ny+which_y,drop=FALSE]
  for(x in seq_len(z$nx)){
   dg<-2*z$w[ii,x]*d1-z$jw[ii,x]*d2
   out[x,]<-out[x,]+colSums(dg*dg)
  }
 }
 ans<-matrix(0,z$m,nb)
 for(j in seq_len(ny)){ii<-which(z$iy==j);cc<-(j-1L)*nb+seq_len(nb);ans[ii,]<-out[z$ix[ii],cc,drop=FALSE]}
 stopifnot(all(is.finite(ans)),all(ans>=0));ans
}

.np_conditional_mp1_reduce <- function(center,s,T,alpha) {
  if (!is.matrix(T) || any(!is.finite(T)) || length(s)!=nrow(T) ||
      any(!is.finite(s)) || any(s<=0) || length(center)!=nrow(T))
    stop("invalid MP1 corrected-root reduction",call.=FALSE)
  m<-nrow(T);B<-ncol(T)
    k<-1L+floor(B*alpha/2);kb<-1L+floor(B*alpha/(2*m))
    pair<-function(x,k)sort(x,partial=c(k,B+1L-k))[c(k,B+1L-k)]
    point<-t(apply(T,1,pair,k=k));bonf<-t(apply(T,1,pair,k=kb))
    critical<-c(sort(apply(T,2,min),partial=k)[k],sort(apply(T,2,max),partial=B+1L-k)[B+1L-k])
    invert<-function(x)cbind(lower=center-s*x[,2],upper=center-s*x[,1])
    bands<-list(pointwise=invert(point),simultaneous=cbind(lower=center-s*critical[2],upper=center-s*critical[1]),bonferroni=invert(bonf))
  list(center=center,se=s,scale=s,critical=critical,bands=bands)
}

# Deterministic worker leaf. The coordinator owns every response draw.
.np_conditional_mp1_block <- function(task,plan) {
  z<-plan$z;Y<-task$value;b<-task$bsz;p0<-plan$pilot
  J<-matrix(0,z$m,b)
  for(j in seq_len(z$ny)) {
    K<-stats::dnorm((z$gy[j]-Y)/z$h[2])/z$h[2]
    KH<-stats::dnorm((z$gy[j]-Y)/z$H[2])/z$H[2]
    RJ<-2*crossprod(z$w,K-p0$mean[,j])-crossprod(z$jw,KH-p0$mean[,z$ny+j])
    ii<-which(z$iy==j);J[ii,]<-RJ[z$ix[ii],,drop=FALSE]
  }
  # sqrt(Vb)^2 is retained deliberately for exact frozen-reference arithmetic.
  replicate.se<-sqrt(.np_conditional_mp1_variance(z,Y,plan$prep))
  between<-.np_conditional_mp1_between(z,Y,p0)
  baseV<-diag(p0$covariance)
  V<-.5*(baseV+replicate.se^2)+.25*between
  if(any(!is.finite(V)) || any(V<=0) || any(V<.5*baseV))
    stop("invalid MP1 mixture variance; no floor or fallback",call.=FALSE)
  se<-sqrt(V);T<-J/se;ratios<-se/sqrt(baseV)
  cbind(t(T),minimum.scale.ratio=apply(ratios,2L,min),maximum.scale.ratio=apply(ratios,2L,max))
}

.np_conditional_mp1_workspace <- function(x,y,grid,h,B) {
  a<-.np_conditional_workspace(x,y,grid,h,B)
  n<-length(x);m<-nrow(grid);ny<-length(unique(grid$y));nx<-length(unique(grid$x))
  # Explicitly budget retained responses, gather/reduction copies, full pilot
  # moments, plan matrices and the fixed 128-draw worker temporaries.
  bytes<-8*(2*n*a$B+4*(m+2)*a$B+4*n*(2*ny)^2+
    10*n*ny*128+16*n*32+10*n*nx+12*m*m)
  if(!is.finite(bytes) || bytes>256*1024^2)
    stop("MP1 numerical workspace exceeds the 256 MiB budget; reduce the requested grid or draws explicitly",call.=FALSE)
  a$estimated.bytes<-max(a$estimated.bytes,bytes);a
}

.np_conditional_mp1_inference <- function(x,y,grid,h,B,alpha) {
  allocation<-.np_conditional_mp1_workspace(x,y,grid,h,B)
  z<-.np_conditional_jack_context(.np_conditional_den_context(x,y,grid,h))
  p<-.np_conditional_jack_pilot(z);s0<-sqrt(diag(p$covariance))
  if(any(!is.finite(s0)) || any(s0<=0))
    stop("nonpositive original MP1 pilot variance",call.=FALSE)
  plan<-list(z=z,pilot=p,prep=.np_conditional_mp1_prepare(z))
  draws<-.np_conditional_dispatch(.np_conditional_tasks(z,B),plan,"jackknife-mp1",z$m+2L)
  T<-draws[seq_len(z$m),,drop=FALSE]
  base<-.np_conditional_mp1_reduce(.np_conditional_jack_center(z),s0,T,alpha)
  out<-.np_conditional_marginal_guard(base,.np_conditional_jack_range(z),
    list(range_interval=.np_conditional_range_interval),alpha)
  out$gaussian_critical<-NULL
  out$bootstrap.critical<-base$critical;out$reference.bands<-base$bands
  out$original<-.np_conditional_den_center(z)$original
  out$se<-s0;out$scale<-s0;out$covariance<-p$covariance
  out$covariance.provenance<-"complete original conditional-pilot jackknife covariance; resample MP1 mixture studentizers retain correction cross terms"
  out$studentization<-"MP1 independent per-observation equal pilot mixtures"
  out$mixture.weight<-.5
  out$studentizer.ratio.range<-c(min(draws[z$m+1L,]),max(draws[z$m+2L,]))
  out$alpha<-alpha;out$B<-B;out$law<-"conditional-pilot";out$correction<-"jackknife"
  out$construction<-"jackknife-mp1";out$calibration<-"asymmetric-sd"
  out$physical.bandwidth<-h;out$grid<-grid
  out$scope<-"experimental LC MP1; one continuous X and Y; fixed Gaussian2; requested finite target family"
  out$estimated.numerical.bytes<-allocation$estimated.bytes
  out$resolution.limited<-B*alpha/(2*z$m)<1
  out
}
