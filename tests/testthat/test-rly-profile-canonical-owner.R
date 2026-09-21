test_that("RLY profile matrices use the same native kernel owner", {
  old<-options(np.messages=FALSE);on.exit(options(old))
  s<-c(1,3,7); x<-data.frame(o=ordered(rep(s,each=6),levels=s))
  y<-sin(seq_len(nrow(x)))
  owner<-getFromNamespace(".np_cat_profile_kernel_matrix","npRmpi")
  for(lambda in c(0,.42,1)) {
    b<-npregbw(xdat=x,ydat=y,bws=lambda,bandwidth.compute=FALSE,
                okertype="racineliyan")
    expected<-lambda^abs(outer(s,s,"-"));expected<-expected/rowSums(expected)
    got<-owner(matrix(s,ncol=1),matrix(s,ncol=1),x,b)
    expect_equal(got,t(expected),tolerance=2e-14)
    expect_equal(owner(matrix(numeric(),0,1),matrix(s,ncol=1),x,b),
                 matrix(numeric(),0,3))
    expect_error(owner(matrix(2,ncol=1),matrix(s,ncol=1),x,b),"outside retained support")
  }
  bridge<-function(t,e,h,s).Call("C_np_ordered_rly_matrix",t,e,h,s,PACKAGE="npRmpi")
  expect_error(bridge(1,1,NA_real_,1),"bandwidth")
  expect_error(bridge(1,1,.4,c(1,1.5)),"integer distances")
  expect_error(bridge(NA_real_,1,.4,1),"outside retained support")
  expect_error(bridge(1,1,.4,c(1,Inf)),"finite")
  oldgc<-gctorture(TRUE);on.exit(gctorture(oldgc),add=TRUE)
  value<-bridge(c(1,3),c(1,3),.4,c(1,3))
  gctorture(oldgc)
  expect_equal(value,
    matrix(c(1,.16,.16,1),2,2)/1.16,tolerance=2e-14)
})

test_that("RLY categorical bootstrap profiles preserve fixed-draw plot bands", {
  old<-options(np.messages=FALSE,np.tree=FALSE,np.categorical.compress=FALSE)
  on.exit(options(old))
  set.seed(648);d<-data.frame(o=ordered(rep(c(1,3,7),each=12),levels=c(1,3,7)))
  d$y<-sin(as.numeric(as.character(d$o)))+rnorm(nrow(d),sd=.2)
  fit<-npreg(y~o,data=d,bws=.42,okertype="racineliyan")
  for(grad in c(FALSE,TRUE)) {
    results<-list(); seeds<-list()
    for(compress in c(FALSE,TRUE)) {
      options(np.categorical.compress=compress);set.seed(831)
      expect_warning(results[[as.character(compress)]]<-plot(fit,plot.behavior="data",
        gradients=grad,errors="bootstrap",plot.errors.boot.method="inid",
        plot.errors.boot.num=9L,band="all"),"B=9 is too small")
      seeds[[as.character(compress)]]<-.Random.seed
    }
    expect_equal(results[["TRUE"]][[1]]$mean,results[["FALSE"]][[1]]$mean,tolerance=2e-12)
    expect_equal(results[["TRUE"]][[1]]$merr,results[["FALSE"]][[1]]$merr,tolerance=2e-12)
    expect_identical(seeds[["TRUE"]],seeds[["FALSE"]])
  }
})
