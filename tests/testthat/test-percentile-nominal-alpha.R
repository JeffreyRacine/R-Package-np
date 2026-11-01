test_that("nominal percentile ranks preserve the supplied alpha", {
new <- .np_udens_mx1
B<-1:9999
for(m in c(1,2,3,50)) {
 stopifnot(identical(1+floor(B*.05/2),1+floor(B*(1-.95)/2)),
  identical(1+floor(B*.05/(2*m)),1+floor(B*(1-.95)/(2*m))))
}
rows<-list();j<-0L
for(B in c(19L,200L,999L))for(m in c(1L,2L,3L,50L))for(alpha in c(.05,.1,.2,.1-.Machine$double.eps,.1+.Machine$double.eps)){
 v<-seq(-1,1,length.out=B);T<-matrix(rep(v,m),nrow=m,byrow=TRUE)
 p<-list(m=m,B=B,tile=17L,bound=c(-2,2),get=function(ii)T[ii,,drop=FALSE])
 nominal<-1-alpha;z<-new$scan(p,nominal=nominal,nominal.alpha=alpha)
 at<-match(nominal,z$q);kp<-1+floor(B*alpha/2);kb<-1+floor(B*alpha/(2*m))
 stopifnot(z$raw.k$pointwise[at]==kp,z$raw.k$bonferroni[at]==kb)
 if(alpha==.05)stopifnot(identical(z,new$scan(p)))
 # Recover complete reduced bands without using the implementation's rank helpers.
 depth<-pmin(seq_len(B),rev(seq_len(B)))
 for(q in c(.501,nominal,1)){
  records<-rep(list(list(row=list(q.required=q),truth=rep(0,m),support=c(-2,2))),39)
  ans<-new$reduce(p,records,warn=FALSE,nominal=nominal,nominal.alpha=alpha)
  qtail<-if(identical(q,nominal))alpha else 1-q
  ks<-min(kp,1+floor(B*qtail/2),sort(depth)[B-ceiling(B*q)+1])
  k<-c(pointwise=if(m==1)ks else kp,simultaneous=ks,bonferroni=min(kb,ks))
  for(b in names(k)){
   expected<-matrix(rep(v[c(k[b],B+1-k[b])],each=m),m,2)
   if(q==1&&(b!='pointwise'||m==1))expected<-matrix(rep(c(-2,2),each=m),m,2)
   stopifnot(identical(unname(ans$bands[[b]]),expected))
  }
 }
 j<-j+1L;rows[[j]]<-data.frame(B,m,alpha,pointwise=kp,bonferroni=kb,pass=TRUE)
}
expect_equal(j,60L)
})
