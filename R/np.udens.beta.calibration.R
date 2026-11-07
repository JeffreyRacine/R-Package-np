# Beta-local calibration: no generic calibration alias or density contract change.
.np_udens_beta_hull <- function(center, lower, upper) {
  .Call(C_np_plot_udens_beta_hull,as.double(center),as.double(lower),as.double(upper))
}

.np_udens_beta_q1_region <- function(theta, Veta, grid, alpha, census) {
  eta<-log(theta+1);L<-t(chol(Veta));radius<-sqrt(qchisq(1-alpha/2,2))
  angle<-seq(0,2*pi,length.out=census+1L)
  transform<-radius*L;shape<-exp(eta+transform%*%rbind(cos(angle),sin(angle)))
  tangent<-transform%*%rbind(-sin(angle),cos(angle))
  a<-shape[1,];b<-shape[2,];psi<-digamma(a+b)
  logx<-log(grid);log1x<-log1p(-grid)
  derivative<-sweep(outer(logx,-digamma(a)+psi,"+"),2,a*tangent[1,],"*")+
    sweep(outer(log1x,-digamma(b)+psi,"+"),2,b*tangent[2,],"*")
  at<-function(t,j,derivative=FALSE) {
    ab<-exp(eta+drop(transform%*%c(cos(t),sin(t))))
    if(!derivative)return((ab[1]-1)*logx[j]+(ab[2]-1)*log1x[j]-lbeta(ab[1],ab[2]))
    du<-drop(transform%*%c(-sin(t),cos(t)))
    sum((c(logx[j],log1x[j])-digamma(ab)+digamma(sum(ab)))*ab*du)
  }
  lower<-upper<-numeric(length(grid));counts<-integer(length(grid))
  for(j in seq_along(grid)) {
    ix<-which(derivative[j,-ncol(derivative)]*derivative[j,-1L]<=0)
    roots<-vapply(ix,function(k)uniroot(function(t)at(t,j,TRUE),angle[k+0:1],
      tol=1e-13,maxiter=1000L,check.conv=TRUE)$root,numeric(1))
    if(length(roots)<2L)stop("beta plot ellipse extrema root census",call.=FALSE)
    vals<-vapply(roots,at,numeric(1),j=j)
    lower[j]<-exp(min(vals));upper[j]<-exp(max(vals));counts[j]<-length(roots)
  }
  if(any(!is.finite(c(lower,upper)))||any(lower>upper))
    stop("beta plot ellipse endpoints not finite/ordered",call.=FALSE)
  list(lower=lower,upper=upper,counts=counts)
}

.np_udens_beta_region_q1 <- function(theta,Veta,grid,alpha) {
  main<-.np_udens_beta_q1_region(theta,Veta,grid,alpha,512L)
  refined<-.np_udens_beta_q1_region(theta,Veta,grid,alpha,1024L)
  if(!identical(main$counts,refined$counts))stop("beta plot ellipse census refinement",call.=FALSE)
  .np_udens_beta_check(c(main$lower,main$upper),c(refined$lower,refined$upper),
    "ellipse extrema",1e-10)
  main
}

.np_udens_beta_sd_reduce <- function(roots,rhat,scale,region,center,alpha,band,calibration) {
  if(!is.matrix(roots)||!nrow(roots)||!ncol(roots)||any(!is.finite(roots))||
     length(scale)!=nrow(roots)||length(rhat)!=nrow(roots)||length(center)!=nrow(roots)||
     length(region$lower)!=nrow(roots)||length(region$upper)!=nrow(roots)||
     any(!is.finite(c(scale,rhat,center,region$lower,region$upper)))||any(scale<=0)||
     length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1||
     !band%in%c("pointwise","simultaneous")||
     !calibration%in%c("symmetric-sd","asymmetric-sd"))stop("invalid beta plot SD reduction",call.=FALSE)
  B<-ncol(roots);k<-floor(B*alpha/4)+1L;upper.rank<-B+1L-k;ka<-ceiling(B*(1-alpha/2))
  if(band=="pointwise") {
    if(calibration=="symmetric-sd") c<-apply(abs(roots),1,function(z)sort(z,partial=ka)[ka]) else {
      l<-apply(roots,1,function(z)sort(z,partial=k)[k])
      u<-apply(roots,1,function(z)sort(z,partial=upper.rank)[upper.rank])
    }
  } else {
    if(calibration=="symmetric-sd")c<-sort(apply(abs(roots),2,max),partial=ka)[ka] else {
      l<-sort(apply(roots,2,min),partial=k)[k]
      u<-sort(apply(roots,2,max),partial=upper.rank)[upper.rank]
    }
  }
  if(calibration=="symmetric-sd") {dl<--c;du<-c} else {dl<--u;du<--l}
  lower<-rhat+region$lower+scale*dl;upper<-rhat+region$upper+scale*du
  if(any(!is.finite(c(lower,upper)))||any(lower>upper))stop("invalid beta plot SD endpoints",call.=FALSE)
  if(calibration=="symmetric-sd") {h<-.np_udens_beta_hull(center,lower,upper);lower<-h[,1];upper<-h[,2]}
  list(lower=lower,upper=upper,k_signed=k,k_absolute=ka,B=B,
       construction=if(calibration=="symmetric-sd")"residual-absolute-containing-hull" else "signed-residual-inversion")
}
