test_that("adaptive LP fitting tiles preserve mixed frames and query identity", {
 old<-options(np.messages=FALSE);on.exit(options(old))
 set.seed(910);n<-192L
 x<-data.frame(x=runif(n,-1,1),u=factor(rep(letters[1:3],length.out=n)),
               o=ordered(rep(c(.5,1,2),length.out=n)))
 y<-sin(x$x)+rnorm(n,sd=.4);ix<-c(7,7,seq_len(67))
 for(cdf in c(FALSE,TRUE))for(kernel in c("epanechnikov","uniform"))
 for(basis in c(FALSE,TRUE)) {
  ctor<-if(cdf)npcdistbw else npcdensbw;fit<-if(cdf)npcdist else npcdens
  b<-ctor(xdat=x,ydat=y,bws=c(144,144,.2,.3),bwtype="adaptive_nn",
    regtype="lp",degree=2L,bernstein.basis=basis,cxkertype=kernel,
    cykertype=kernel,oxkertype="racineliyan",bandwidth.compute=FALSE)
  args<-list(bws=b,txdat=x,tydat=y,exdat=x[ix,],eydat=y[ix]+.01,
             se=TRUE,gradients=TRUE)
  options(np.tree=FALSE);off<-do.call(fit,args)
  options(np.tree=TRUE);on<-do.call(fit,args)
  for(field in c(if(cdf)"condist"else"condens","conderr","congrad","congerr"))
    expect_equal(on[[field]],off[[field]],tolerance=2e-9)
  expect_equal(fitted(on)[1],fitted(on)[2],tolerance=0)
  perm<-sample.int(n);args$txdat<-x[perm,];args$tydat<-y[perm]
  permuted<-do.call(fit,args)
  expect_equal(fitted(on),fitted(permuted),tolerance=2e-9)
  expect_equal(se(on),se(permuted),tolerance=2e-9)
 }
})
