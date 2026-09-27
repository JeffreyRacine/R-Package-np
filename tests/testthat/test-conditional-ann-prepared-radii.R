ann_lp_conditional_literal <- function(x,y,ex,ey,k,degree,kernel,cdf) {
 x<-as.matrix(x);ex<-as.matrix(ex);n<-nrow(x);p<-ncol(x)
 terms<-as.matrix(expand.grid(rep(list(0:degree),p)))
 terms<-terms[rowSums(terms)<=degree,,drop=FALSE]
 design<-function(z)vapply(seq_len(nrow(terms)),function(a)
   apply(sweep(z,2,terms[a,],'^'),1,prod),numeric(nrow(z)))
 K<-function(z)switch(kernel,gaussian=dnorm(z),uniform=.5*(abs(z)<1),
   epanechnikov=3/(4*sqrt(5))*pmax(0,1-z*z/5))
 P<-function(z) {
   if(kernel=="gaussian")return(pnorm(z))
   q<-pmax(-1,pmin(1,z/if(kernel=="uniform")1 else sqrt(5)))
   if(kernel=="uniform")(q+1)/2 else .5+.75*q-.25*q^3
 }
 radius<-function(z)vapply(seq_len(n),function(i)sort(abs(z[-i]-z[i]))[k],0.0)
 hx<-vapply(seq_len(p),function(l)radius(x[,l]),numeric(n))
 hy<-radius(y);B<-design(x);E<-design(ex)
 vapply(seq_len(nrow(ex)),function(j) {
   w<-rep(1,n)
   for(l in seq_len(p))w<-w*K((ex[j,l]-x[,l])/hx[,l])/hx[,l]
   v<-if(cdf)P((ey[j]-y)/hy)else K((ey[j]-y)/hy)/hy
   sum(E[j,]*solve(crossprod(B,w*B),crossprod(B,w*v)))
 },0.0)
}
test_that("adaptive conditional prepared radii preserve donor geometry and WLS", {
 old<-options(np.messages=FALSE);on.exit(options(old))
 set.seed(679);n<-128L
 x0<-data.frame(x1=runif(n,-1,1),x2=runif(n,-1,1));x0[2,]<-x0[1,]
 y<-sin(x0$x1)+.4*x0$x2+rnorm(n,sd=.4)
 for(cdf in c(FALSE,TRUE))for(kernel in c("epanechnikov","uniform","gaussian"))
 for(q in list(c(1,0),c(1,1),c(1,3),c(2,2))) {
  p<-q[1];degree<-q[2];x<-x0[,seq_len(p),drop=FALSE];k<-96L
  ix<-c(18,3,3,85,100);ex<-x[ix,,drop=FALSE];ey<-y[ix]+.001
  ctor<-if(cdf)npcdistbw else npcdensbw;fit<-if(cdf)npcdist else npcdens
  b<-ctor(xdat=x,ydat=y,bws=rep(k,p+1),bwtype="adaptive_nn",bwmethod="cv.ls",
    regtype="lp",degree=rep(degree,p),bernstein.basis=TRUE,
    cxkertype=kernel,cykertype=kernel,bandwidth.compute=FALSE)
  for(external in c(FALSE,TRUE)) {
   args<-list(bws=b,txdat=x,tydat=y,se=TRUE,gradients=degree>0)
   e<-if(external)ex else x;v<-if(external)ey else y
   if(external){args$exdat<-e;args$eydat<-v}
   options(np.tree=FALSE);off<-do.call(fit,args)
   options(np.tree=TRUE);on<-do.call(fit,args)
   ref<-ann_lp_conditional_literal(x,y,e,v,k,degree,kernel,cdf)
   info<-paste(cdf,kernel,p,degree,external)
   expect_true(max(abs(fitted(off)-ref))<2e-9,info=info)
   expect_true(max(abs(fitted(on)-ref))<2e-9,info=info)
   for(field in c(if(cdf)"condist"else"condens","conderr","congrad","congerr"))
     expect_equal(on[[field]],off[[field]],tolerance=2e-9,info=info)
  }
 }
})
test_that("adaptive prepared radii survive zero-radius refusal and donor permutation", {
 old<-options(np.messages=FALSE,np.tree=TRUE);on.exit(options(old))
 set.seed(680);x<-data.frame(x=c(rep(0,6),runif(26,.1,1)));y<-rnorm(32)
 for(cdf in c(FALSE,TRUE)) {
  ctor<-if(cdf)npcdistbw else npcdensbw;fit<-if(cdf)npcdist else npcdens
  make<-function(k)ctor(xdat=x,ydat=y,bws=c(24,k),bwtype="adaptive_nn",
   regtype="lp",degree=2L,bernstein.basis=TRUE,cxkertype="epanechnikov",
   cykertype="epanechnikov",bandwidth.compute=FALSE)
  bad<-make(3L)
  expect_error(fit(bws=bad,txdat=x,tydat=y,se=FALSE),"zero literal radius")
  b<-make(24L);ix<-c(9,9,14)
  a<-fit(bws=b,txdat=x,tydat=y,exdat=x[ix,,drop=FALSE],eydat=y[ix],se=TRUE)
  perm<-sample(seq_len(nrow(x)))
  z<-fit(bws=b,txdat=x[perm,,drop=FALSE],tydat=y[perm],
    exdat=x[ix,,drop=FALSE],eydat=y[ix],se=TRUE)
  expect_true(max(abs(fitted(a)-fitted(z)))<2e-9)
  expect_true(max(abs(se(a)-se(z)))<2e-9)
 }
})
