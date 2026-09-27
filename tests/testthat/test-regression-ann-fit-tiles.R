# Independent fitting oracle: full-sample, donor-centred ANN radii.
ann_literal<-function(x,y,e,k,degree,kernel){
 x<-as.matrix(x);e<-as.matrix(e);n<-nrow(x);p<-ncol(x)
 powers<-as.matrix(expand.grid(rep(list(0:degree),p)))
 powers<-powers[rowSums(powers)<=degree,,drop=FALSE]
 design<-function(z)vapply(seq_len(nrow(powers)),function(a)
  apply(sweep(z,2,powers[a,],'^'),1,prod),numeric(nrow(z)))
 B<-design(x);E<-design(e)
 radii<-vapply(seq_len(p),function(d)vapply(seq_len(n),function(i)
  sort(abs(x[-i,d]-x[i,d]))[k],0),numeric(n))
 K<-if(kernel=='uniform')function(u).5*(abs(u)<1)else
  function(u)3/(4*sqrt(5))*pmax(0,1-u*u/5)
 vapply(seq_len(nrow(e)),function(j){
  w<-rep(1,n)
  for(d in seq_len(p))w<-w*K((e[j,d]-x[,d])/radii[,d])/radii[,d]
  sum(E[j,]*solve(crossprod(B,w*B),crossprod(B,w*y)))
 },0)
}
test_that('ANN LP fitting retains literal donor-radius WLS and uncertainty',{
 old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old))
 set.seed(417);n<-144L;x<-data.frame(x=runif(n,-1,1),z=runif(n,-1,1))
 y<-sin(x$x)+x$z^2+rnorm(n,sd=.1)
 for(p in 1:2)for(degree in 1:2)for(kernel in c('epanechnikov','uniform'))
  for(bernstein in c(FALSE,TRUE))for(external in c(FALSE,TRUE)){
   xx<-x[,seq_len(p),drop=FALSE];ee<-if(external)xx[c(7,7,seq(9,141,by=2)),,drop=FALSE]else xx
   b<-npregbw(xdat=xx,ydat=y,bws=rep(112,p),bwtype='adaptive_nn',regtype='lp',
    degree=rep(degree,p),bernstein.basis=bernstein,ckertype=kernel,bandwidth.compute=FALSE)
   args<-list(bws=b,txdat=xx,tydat=y,se=TRUE,gradients=TRUE)
   if(external)args$exdat<-ee
   options(np.tree=FALSE);off<-do.call(npreg,args)
   options(np.tree=TRUE);on<-do.call(npreg,args)
   ref<-ann_literal(xx,y,ee,112,degree,kernel)
   check<-function(u,v)expect_true(max(abs(u-v))<2e-9,
    info=paste(p,degree,kernel,bernstein,external,max(abs(u-v))))
   check(fitted(off),ref);check(fitted(on),ref)
   check(se(on),se(off));check(gradients(on),gradients(off))
   check(gradients(on,se=TRUE),gradients(off,se=TRUE))
  }
})
test_that('ANN LP mixed external consumers preserve query and donor identity',{
 old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old))
 set.seed(418);n<-96L
 x<-data.frame(x=runif(n,-1,1),u=factor(rep(letters[1:3],32)),o=ordered(rep(1:3,32)))
 y<-sin(x$x)+as.integer(x$u)/8+rnorm(n,sd=.1)
 e<-x[c(5,5,29),]
 for(degree in 1:2)for(permute in c(FALSE,TRUE)){
  ix<-if(permute)sample.int(n)else seq_len(n);xx<-x[ix,];yy<-y[ix]
  b<-npregbw(xdat=xx,ydat=yy,bws=c(72,.2,.3),bwtype='adaptive_nn',regtype='lp',
   degree=degree,bernstein.basis=TRUE,ckertype='epanechnikov',bandwidth.compute=FALSE)
  out<-lapply(c(FALSE,TRUE),function(tree){
   options(np.tree=tree);fit<-npreg(bws=b,txdat=xx,tydat=yy,exdat=e,se=TRUE,gradients=TRUE)
   H<-npreghat(bws=b,txdat=xx,exdat=e,output='matrix')
   Y<-cbind(yy,cos(xx$x));saved<-Y
   A<-npreghat(bws=b,txdat=xx,exdat=e,y=Y,output='apply')
   L<-npreghat(bws=b,txdat=xx,y=Y,output='apply',leave.one.out=TRUE)
   expect_true(max(abs(A-H%*%Y))<2e-9);expect_identical(Y,saved)
   expect_true(max(abs(A[,1]-fitted(fit)))<2e-9)
   c(fitted(fit),se(fit),gradients(fit),gradients(fit,se=TRUE),H,A,L)
  })
  expect_true(max(abs(out[[1]]-out[[2]]))<2e-9)
 }
})
