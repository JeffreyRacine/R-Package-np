# R77 conditional-CDF Gaussian8 deleted-fold reference

The existing n72 seed2026072702 fixture, raw GLP(1,1), response k25 and explanatory k29/27, has CVLS3.1542637393732895 from independent R order statistics, scalar dnorm/pnorm, deleted weighted least squares in the space (1,x1,x2), and empirical loss across j != i. Every radius is recomputed on donors excluding only i. No package radius, basis, solve or objective is used in the reference.

The established Gaussian8 kernel uses decimal coefficient0.02083333333 for z^6. Replacing it with exact1/48 gives3.1542650432424995 on this signed ill-conditioned fixture; this pre-existing kernel-definition difference is retained, not silently confused with the CV repair. The independent reference uses the package's existing decimal kernel definition. Candidate native3.1542637389592714 differs by4.14e-10, below the unchanged5e-8 test tolerance. Remaining density/fixed/ANN transcripts stay unchanged.

Evidence and independently implemented equations: Development/tmp/np_r77_conditional_cdf_gnn_cv_20261114_r1/scripts/signed_oracle.R and evidence/signed_fixture.csv. The initial exact-rational comparison failure is retained in logs/signed_independent_oracle.log.

Self-contained reproduction (base R only):

```r
set.seed(2026072702L);n<-72L
x<-data.frame(x1=runif(n,.05,.95),x2=runif(n,.05,.95));y<-.4+.3*sin(2*pi*x$x1)*cos(pi*x$x2)+rnorm(n,sd=.09)
f<-function(coef){mean(vapply(1:n,function(i){
 d<-setdiff(1:n,i);hx<-c(sort(abs(x$x1[d]-x$x1[i]))[29],sort(abs(x$x2[d]-x$x2[i]))[27]);z<-sweep(sweep(as.matrix(x[d,]),2,as.numeric(x[i,]),'-'),2,hx,'/')
 K<-function(z)dnorm(z)*(35/16-35/16*z^2+7/16*z^4-coef*z^6);w<-K(z[,1])*K(z[,2]);A<-cbind(1,z);v<-as.vector(A%*%solve(crossprod(A,w*A),c(1,0,0)))*w
 mean(vapply(d,function(j){hy<-sort(abs(y[d]-y[j]))[25];p<-sum(v*pnorm((y[j]-y[d])/hy));(as.numeric(y[i]<=y[j])-p)^2},0))
},0))}
f(.02083333333) # independent reference for the existing kernel
```
